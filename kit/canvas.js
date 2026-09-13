// Drawing 2D arrays to a canvas, plus the scalebar convention the lab's
// matplotlib figures use.
//
// The Python renders through matplotlib's `imshow` at a fixed pixel width
// (600-680 px, dpi=72, transparent figure patch). Canvas gives us the same
// fixed-pixel contract with far less machinery.

import {applyColormap, complexToRGB} from "./color.js";

/**
 * Wrap an RGBA byte array as ImageData, with the first axis **vertical**.
 *
 * Our arrays are row-major `[ix * ny + iy]`, and here `ix` is the first
 * (vertical) axis, matching numpy — so the image is `ny` wide and `nx` tall,
 * exactly like `imshow`. This is the right choice when the array came from an
 * image in the first place.
 *
 * When the first axis is a physical **x** instead — a scan whose first index
 * runs along the field's width, a diffraction pattern indexed by `kx` — use
 * {@link transposeToImageData}, which puts it horizontal.
 */
export function toImageData(rgba, nx, ny) {
  return new ImageData(rgba, ny, nx);
}

/**
 * Wrap an RGBA byte array as ImageData, with the first axis **horizontal**.
 *
 * The counterpart to {@link toImageData}: this one transposes, so `nx` becomes
 * the image's width. Reach for it whenever the first array index is a physical
 * x — otherwise the picture comes out on its side, which on a square array is a
 * silent 90° rotation and on a rectangular one a stretch. Nothing throws either
 * way, so the mistake survives until someone notices the image looks wrong.
 *
 * @param {Uint8ClampedArray} rgba 4 bytes per sample, `[ix * ny + iy]`
 * @param {number} nx samples along x — the returned image's width
 * @param {number} ny samples along y — its height
 */
export function transposeToImageData(rgba, nx, ny) {
  const out = new Uint8ClampedArray(nx * ny * 4);
  for (let ix = 0; ix < nx; ix++) {
    for (let iy = 0; iy < ny; iy++) {
      const from = (ix * ny + iy) * 4;
      const to = (iy * nx + ix) * 4;
      out[to] = rgba[from];
      out[to + 1] = rgba[from + 1];
      out[to + 2] = rgba[from + 2];
      out[to + 3] = rgba[from + 3];
    }
  }
  return new ImageData(out, nx, ny);
}

/**
 * Draw a real-valued array to a canvas through a colormap.
 *
 * @param {HTMLCanvasElement} canvas
 * @param {ArrayLike<number>} array row-major, nx*ny
 * @param {number} nx rows
 * @param {number} ny columns
 * @param {object} [options] forwarded to {@link applyColormap}, plus:
 * @param {boolean} [options.smooth=false] whether to interpolate when scaling up
 */
export function drawArray(canvas, array, nx, ny, options = {}) {
  return blit(canvas, applyColormap(array, options), nx, ny, options.smooth ?? false);
}

/**
 * Draw a complex array to a canvas as a domain-coloured image.
 *
 * @param {HTMLCanvasElement} canvas
 * @param {{re: Float64Array, im: Float64Array}} array
 * @param {object} [options] forwarded to {@link complexToRGB}
 */
export function drawComplex(canvas, array, nx, ny, options = {}) {
  return blit(canvas, complexToRGB(array, options), nx, ny, options.smooth ?? false);
}

// How many backing-store pixels one CSS pixel is worth, per canvas. Overlays
// (scalebars, quivers) need this to size themselves in *display* units rather
// than array units — otherwise text drawn onto a 128x128 backing store gets
// scaled up with the image and comes out huge.
const canvasScale = new WeakMap();

/**
 * Size a canvas's backing store and record its scale.
 *
 * `size` matters because canvases are usually still detached when we draw into
 * them (they are being assembled inside an `htl.html` template), so `clientWidth`
 * reads 0 and we have nothing else to go on.
 */
function sizeCanvas(canvas, nx, ny, size) {
  const dpr = globalThis.devicePixelRatio ?? 1;
  const cssWidth = canvas.clientWidth || size || ny;
  const cssHeight = canvas.clientHeight || size || nx;
  canvas.width = Math.round(cssWidth * dpr);
  canvas.height = Math.round(cssHeight * dpr);
  canvasScale.set(canvas, canvas.width / cssWidth);
  return canvas.getContext("2d");
}

/**
 * Paint RGBA data onto a canvas, scaling it to the canvas's display size and
 * honouring devicePixelRatio so images stay crisp on retina displays.
 *
 * @param {HTMLCanvasElement} canvas
 * @param {Uint8ClampedArray} rgba
 * @param {number} nx rows
 * @param {number} ny columns
 * @param {object|boolean} [options] `{smooth, size}` — pass `size` (CSS px) when
 *   the canvas is not in the document yet. A bare boolean is read as `smooth`.
 */
export function blit(canvas, rgba, nx, ny, options = {}) {
  const {smooth = false, size} = typeof options === "boolean" ? {smooth: options} : options;

  // Render at native array resolution first, then scale onto the visible canvas.
  // Smoothing off by default reproduces matplotlib's nearest-neighbour `imshow`.
  const source = new OffscreenCanvas(ny, nx);
  source.getContext("2d").putImageData(toImageData(rgba, nx, ny), 0, 0);

  const ctx = sizeCanvas(canvas, nx, ny, size);
  ctx.imageSmoothingEnabled = smooth;
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.drawImage(source, 0, 0, canvas.width, canvas.height);
  return ctx;
}

/**
 * Resolve the element's inherited text colour to something canvas can use.
 *
 * The 2D context does NOT understand `currentColor` (or `color-mix`): assigning
 * one silently leaves the previous value, which is why hand-drawn overlays come
 * out black in dark mode. Read the computed colour instead.
 */
export function currentColor(element, fallback = "#000") {
  return globalThis.getComputedStyle?.(element).color || fallback;
}

/**
 * Draw a 2D displacement/vector field as arrows on a transparent canvas.
 *
 * Arrows are scaled so the largest one spans roughly one sampling step, which
 * keeps the field readable whatever the magnitude.
 *
 * @param {HTMLCanvasElement} canvas
 * @param {{uRow: ArrayLike<number>, uCol: ArrayLike<number>}} field row/column displacement
 * @param {number} nx
 * @param {number} ny
 * @param {object} [options]
 * @param {number} [options.step=16] sample every `step` points
 * @param {string} [options.color] defaults to the canvas's inherited text colour
 * @param {number} [options.backgroundOpacity=0.06] faint ground so an all-zero
 *   field still reads as a panel
 * @param {number} [options.size] CSS size in px, for when the canvas is still
 *   detached and `clientWidth` would read 0 — which is the normal case when
 *   building a figure inside an `htl.html` template
 */
export function drawQuiver(canvas, {uRow, uCol}, nx, ny, options = {}) {
  const {step = 16, color = currentColor(canvas), backgroundOpacity = 0.06, size} = options;

  const ctx = sizeCanvas(canvas, nx, ny, size);
  const dpr = canvasScale.get(canvas) ?? 1;
  ctx.clearRect(0, 0, canvas.width, canvas.height);

  if (backgroundOpacity > 0) {
    ctx.save();
    ctx.globalAlpha = backgroundOpacity;
    ctx.fillStyle = color;
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.restore();
  }

  let maxU = 0;
  for (let i = 0; i < uRow.length; i++) maxU = Math.max(maxU, Math.hypot(uRow[i], uCol[i]));
  if (!(maxU > 0) || !isFinite(maxU)) return ctx; // nothing to draw

  const pxPerSample = canvas.width / ny;
  const scale = (step * pxPerSample) / maxU;

  ctx.strokeStyle = color;
  ctx.fillStyle = color;
  ctx.lineWidth = Math.max(1, dpr);
  ctx.globalAlpha = 0.8;

  const half = step >> 1;
  for (let ix = half; ix < nx; ix += step) {
    for (let iy = half; iy < ny; iy += step) {
      const i = ix * ny + iy;
      const dx = uCol[i] * scale;
      const dy = uRow[i] * scale;
      if (!isFinite(dx) || !isFinite(dy)) continue;

      const x0 = (iy + 0.5) * pxPerSample;
      const y0 = (ix + 0.5) * pxPerSample;
      ctx.beginPath();
      ctx.moveTo(x0, y0);
      ctx.lineTo(x0 + dx, y0 + dy);
      ctx.stroke();

      const length = Math.hypot(dx, dy);
      if (length > 3) {
        const ux = dx / length;
        const uy = dy / length;
        const head = Math.min(5 * dpr, length * 0.4);
        ctx.beginPath();
        ctx.moveTo(x0 + dx, y0 + dy);
        ctx.lineTo(x0 + dx - head * (ux + uy * 0.5), y0 + dy - head * (uy - ux * 0.5));
        ctx.lineTo(x0 + dx - head * (ux - uy * 0.5), y0 + dy - head * (uy + ux * 0.5));
        ctx.closePath();
        ctx.fill();
      }
    }
  }
  return ctx;
}

/**
 * Lay an Observable `Inputs.form` out as a responsive grid.
 *
 * `Inputs.form` stacks its children one per row, which gets unwieldy past three
 * or four sliders. This flows them into `columns` columns instead, dropping to
 * fewer when the container gets too narrow for each one to stay legible.
 *
 * @param {HTMLElement} form the element returned by `Inputs.form`
 * @param {object} [options]
 * @param {number} [options.columns=2] columns to use when there is room
 * @param {number} [options.minWidth=240] below this, use fewer columns
 * @param {number} [options.gap=18] column gap in px
 */
export function gridForm(form, {columns = 2, minWidth = 240, gap = 18} = {}) {
  // The track width is whatever an even split would give, floored at `minWidth`
  // so `auto-fit` drops a column rather than squeezing the labels.
  const track = `max(${minWidth}px, calc((100% - ${(columns - 1) * gap}px) / ${columns}))`;
  form.style.display = "grid";
  form.style.gridTemplateColumns = `repeat(auto-fit, minmax(${track}, 1fr))`;
  form.style.columnGap = `${gap}px`;
  form.style.rowGap = "2px";
  form.style.alignItems = "center";
  return form;
}

/**
 * Choose a "nice" scalebar length: the largest 1/2/5 x 10^n that fits in about a
 * quarter of the field of view.
 *
 * Ported from `estimate_scalebar_length` in `13.upsampled-example.ipynb`.
 *
 * @param {number} extent full width of the image in physical units
 * @param {number} [fraction=0.25] target fraction of the width
 * @returns {number} bar length in the same physical units
 */
export function niceScalebarLength(extent, fraction = 0.25) {
  const target = extent * fraction;
  const decade = Math.pow(10, Math.floor(Math.log10(target)));
  const normalized = target / decade;
  const step = normalized >= 5 ? 5 : normalized >= 2 ? 2 : 1;
  return step * decade;
}

/**
 * Draw a scalebar in the lower-right corner, matching the look of
 * `ctf/visualize.py:add_scalebar`.
 *
 * @param {CanvasRenderingContext2D} ctx
 * @param {object} options
 * @param {number} options.extent physical width of the whole image
 * @param {string} [options.units="Å"]
 * @param {number} [options.length] bar length in physical units; defaults to a nice value
 * @param {string} [options.color="white"]
 */
export function drawScalebar(ctx, {extent, units = "Å", length, color = "white", pad = 12}) {
  const canvas = ctx.canvas;
  const barUnits = length ?? niceScalebarLength(extent);
  const pxPerUnit = canvas.width / extent;
  const barPx = barUnits * pxPerUnit;
  // Draw the bar and its label in display pixels, not array pixels — `blit`
  // recorded the ratio when it sized the backing store.
  const scale = canvasScale.get(canvas) ?? 1;

  const height = Math.max(2, Math.round(3 * scale));
  const padPx = pad * scale;
  const x = canvas.width - padPx - barPx;
  const y = canvas.height - padPx - height;

  ctx.save();
  ctx.fillStyle = color;
  ctx.fillRect(x, y, barPx, height);

  const fontPx = Math.round(12 * scale);
  ctx.font = `${fontPx}px ui-sans-serif, system-ui, sans-serif`;
  ctx.textAlign = "center";
  ctx.textBaseline = "bottom";
  ctx.fillText(`${formatLength(barUnits)} ${units}`, x + barPx / 2, y - 2 * scale);
  ctx.restore();
  return barUnits;
}

function formatLength(v) {
  if (v >= 100) return v.toFixed(0);
  if (v >= 10) return v.toFixed(0);
  if (v >= 1) return v.toFixed(1).replace(/\.0$/, "");
  return v.toPrecision(2);
}

/**
 * A horizontal colorbar strip for a named colormap, as a standalone canvas.
 *
 * @param {string} name colormap name
 * @param {number} [width=200]
 * @param {number} [height=12]
 * @returns {HTMLCanvasElement}
 */
export function colorbar(name, width = 200, height = 12) {
  const canvas = document.createElement("canvas");
  canvas.width = width;
  canvas.height = height;
  const ramp = new Float64Array(width);
  for (let i = 0; i < width; i++) ramp[i] = i / (width - 1);
  const rgba = applyColormap(ramp, {colormap: name});

  const strip = new Uint8ClampedArray(width * height * 4);
  for (let y = 0; y < height; y++) strip.set(rgba, y * width * 4);
  canvas.getContext("2d").putImageData(new ImageData(strip, width, height), 0, 0);
  return canvas;
}

/**
 * Convert a pointer event to array indices for a canvas showing an nx-by-ny array.
 *
 * Handles the CSS-vs-backing-store scaling, so widgets can wire up hover/drag
 * interaction (the `motion_notify_event` pattern used by most of the lab's
 * Python widgets) in one line.
 *
 * @returns {{ix: number, iy: number, inside: boolean}}
 */
export function pointerToIndices(event, canvas, nx, ny) {
  const rect = canvas.getBoundingClientRect();
  const fx = (event.clientX - rect.left) / rect.width;
  const fy = (event.clientY - rect.top) / rect.height;
  const iy = Math.floor(fx * ny);
  const ix = Math.floor(fy * nx);
  return {ix, iy, inside: ix >= 0 && ix < nx && iy >= 0 && iy < ny};
}
