// A small 3D scene renderer for canvas 2D.
//
// Everything here is a projection plus painter's algorithm — no WebGL, no
// three.js, no dependencies. That is also how the widgets this borrows from work
// (ophusgroup/landing's `stem-experiment.js`, `ptycho-ms.js`, `diffraction-sim.js`);
// only their decorative logo widgets pull three.js.
//
// CONVENTION: `z` is the optic axis and maps almost straight to screen-y, so a
// scene reads as a column — aperture at the top, specimen in the middle,
// detector at the bottom. Positive z is up-beam (toward the source).
//
// Draw order matters. To make a specimen occlude the beam, draw the cone below
// the sample, then the specimen, then the cone above it:
//
//   drawProbeCone(ctx, view, {...opts, part: "lower"});
//   drawSpheres(ctx, view, atoms);
//   drawProbeCone(ctx, view, {...opts, part: "upper"});

// ---------------------------------------------------------------------------
// Camera
// ---------------------------------------------------------------------------

/**
 * Build a camera: rotate about the optic axis by `azDeg`, tilt by `elDeg`, then
 * apply a weak perspective divide.
 *
 * Sensible values: azimuth 0 (the reference widgets never orbit), elevation
 * −14° to −30°, and a *small* focal length (45–140). The `+200` bias in
 * {@link project} keeps the perspective mild — 10–20% size falloff front to
 * back, which is what stops the scene looking flatly isometric. `fl = 0`
 * disables perspective entirely (orthographic).
 *
 * @param {number} [azDeg=0] azimuth [degrees]
 * @param {number} [elDeg=-26] elevation [degrees]; negative looks down on the scene
 * @param {number} [zoom=20] screen pixels per scene unit
 * @param {number} [cx=0] screen x of the origin
 * @param {number} [cy=0] screen y of the origin
 * @param {number} [fl=120] focal length; 0 for orthographic
 */
export function makeView(azDeg = 0, elDeg = -26, zoom = 20, cx = 0, cy = 0, fl = 120) {
  const az = (azDeg * Math.PI) / 180;
  const el = (elDeg * Math.PI) / 180;
  return {
    ca: Math.cos(az), sa: Math.sin(az),
    ce: Math.cos(el), se: Math.sin(el),
    zoom, cx, cy, fl
  };
}

/** Perspective bias — see {@link makeView}. */
const PERSPECTIVE_BIAS = 200;

/**
 * Project a scene point to screen coordinates.
 *
 * @returns {{sx: number, sy: number, depth: number}} `depth` **increases toward
 *   the camera** — a larger depth is nearer, and projects bigger. {@link depthSort}
 *   sorts ascending, which is therefore back-to-front, and the depth cueing in
 *   {@link drawSpheres} brightens and enlarges with increasing depth.
 */
export function project(x, y, z, v) {
  const rx = x * v.ca - y * v.sa;
  const ry = x * v.sa + y * v.ca;
  const depth = ry * v.ce - z * v.se;
  const scale = v.fl > 0 ? v.fl / (v.fl - depth + PERSPECTIVE_BIAS) : 1;
  return {
    sx: rx * v.zoom * scale + v.cx,
    sy: (-ry * v.se - z * v.ce) * v.zoom * scale + v.cy,
    depth
  };
}

/**
 * Invert {@link project} **on the z = 0 plane** — enough to turn a pointer
 * position into a specimen coordinate, which is the only inverse a scene needs.
 *
 * Solved exactly rather than approximated. The perspective scale depends on the
 * point's own depth, so it cannot simply be divided out; but on z = 0 the depth
 * is `ry·cos(el)`, which makes the screen-y equation linear in `ry`:
 *
 *     sy - cy = -ry·se·zoom·fl / (B - ry·ce)   with B = fl + bias
 *  => ry = (sy - cy)·B / (-se·zoom·fl + (sy - cy)·ce)
 *
 * (The reference widgets approximate this and accept ~0.2 scene units of error;
 * exact costs nothing and makes pointer picking pixel-accurate.)
 *
 * Requires a non-zero elevation: at `el = 0` the plane is edge-on and the
 * inverse genuinely does not exist.
 */
export function unproject(sx, sy, v) {
  const S = sy - v.cy;
  let rx, ry;
  if (Math.abs(v.se) < 1e-12) {
    ry = 0; // edge-on: the plane projects to a line, so y is unrecoverable
    rx = (sx - v.cx) / v.zoom;
  } else if (v.fl > 0) {
    const B = v.fl + PERSPECTIVE_BIAS;
    const denom = -v.se * v.zoom * v.fl + S * v.ce;
    ry = Math.abs(denom) < 1e-12 ? 0 : (S * B) / denom;
    const scale = v.fl / (B - ry * v.ce);
    rx = (sx - v.cx) / (v.zoom * scale);
  } else {
    ry = -S / (v.se * v.zoom);
    rx = (sx - v.cx) / v.zoom;
  }
  return {x: v.ca * rx + v.sa * ry, y: -v.sa * rx + v.ca * ry};
}

/** Sort drawables back to front. Items need a numeric `depth`. */
export function depthSort(items) {
  return items.sort((a, b) => a.depth - b.depth);
}

/** Normalised depth 0..1 across an array of projected points, for depth cueing. */
export function depthRange(items) {
  let min = Infinity, max = -Infinity;
  for (const it of items) {
    if (it.depth < min) min = it.depth;
    if (it.depth > max) max = it.depth;
  }
  const span = max - min || 1;
  return {min, max, t: (d) => (d - min) / span};
}

// ---------------------------------------------------------------------------
// Quaternion trackball, for tilting a specimen or crystal
// ---------------------------------------------------------------------------

/** Quaternion from an axis and angle [rad]. */
export function quatFromAxisAngle([x, y, z], angle) {
  const n = Math.hypot(x, y, z) || 1;
  const s = Math.sin(angle / 2) / n;
  return [Math.cos(angle / 2), x * s, y * s, z * s];
}

/** Hamilton product. */
export function quatMultiply(a, b) {
  const [aw, ax, ay, az] = a;
  const [bw, bx, by, bz] = b;
  return [
    aw * bw - ax * bx - ay * by - az * bz,
    aw * bx + ax * bw + ay * bz - az * by,
    aw * by - ax * bz + ay * bw + az * bx,
    aw * bz + ax * by - ay * bx + az * bw
  ];
}

/** Unit-normalise a quaternion, so repeated composition cannot drift. */
export function quatNormalize(q) {
  const n = Math.hypot(q[0], q[1], q[2], q[3]) || 1;
  return [q[0] / n, q[1] / n, q[2] / n, q[3] / n];
}

/** Quaternion to a row-major 3x3 rotation matrix. */
export function quatToMatrix([w, x, y, z]) {
  return [
    1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w),
    2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w),
    2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)
  ];
}

/** Apply a 3x3 row-major matrix to a vector. */
export function applyMatrix(m, [x, y, z]) {
  return [
    m[0] * x + m[1] * y + m[2] * z,
    m[3] * x + m[4] * y + m[5] * z,
    m[6] * x + m[7] * y + m[8] * z
  ];
}

/**
 * Trackball drag: compose `quat` with the rotation implied by a pointer motion
 * of (dx, dy) pixels across a control of `size` pixels.
 *
 * Ported from `diffraction-sim.js`, where dragging tilts the crystal so the near
 * face follows the pointer.
 */
export function trackball(quat, dx, dy, size) {
  const angle = (Math.hypot(dx, dy) * Math.PI) / size;
  if (!(angle > 0)) return quat;
  return quatNormalize(quatMultiply(quatFromAxisAngle([dy, dx, 0], angle), quat));
}

/** In-plane twist about the optic axis, for shift-drag / two-finger rotation. */
export function twist(quat, angle) {
  return quatNormalize(quatMultiply(quatFromAxisAngle([0, 0, 1], angle), quat));
}

// ---------------------------------------------------------------------------
// Primitives
// ---------------------------------------------------------------------------

/** Points of a circle of radius `r` at height `z`, already projected. */
export function ellipsePoints(view, x, y, z, r, segments = 48, projectFn = project) {
  const pts = [];
  for (let i = 0; i < segments; i++) {
    const a = (i / segments) * 2 * Math.PI;
    pts.push(projectFn(x + Math.cos(a) * r, y + Math.sin(a) * r, z, view));
  }
  return pts;
}

/**
 * Stroke a circle at height `z` — it projects to an ellipse for free.
 *
 * @param {object} [options]
 * @param {"near"|"far"|null} [options.half=null] draw only the half of the ring
 *   nearer to, or further from, the camera. Drawing the far half before an object
 *   and the near half after it puts the ring *around* that object rather than
 *   flatly in front of or behind it — which is how a beam can be shown passing
 *   through a specimen that was rendered in one piece.
 */
export function ring(ctx, view, x, y, z, r, {segments = 48, half = null, projectFn = project} = {}) {
  const pts = ellipsePoints(view, x, y, z, r, segments, projectFn);
  if (!half) {
    ctx.beginPath();
    pts.forEach((p, i) => (i ? ctx.lineTo(p.sx, p.sy) : ctx.moveTo(p.sx, p.sy)));
    ctx.closePath();
    ctx.stroke();
    return;
  }

  // The half is contiguous around the ellipse, so walk from the point that
  // starts it and stop when the depth test flips.
  const axisDepth = projectFn(x, y, z, view).depth;
  const wanted = (p) => (half === "near" ? p.depth >= axisDepth : p.depth <= axisDepth);
  let start = -1;
  for (let i = 0; i < pts.length; i++) {
    if (wanted(pts[i]) && !wanted(pts[(i - 1 + pts.length) % pts.length])) { start = i; break; }
  }
  if (start < 0) {
    // Every point is on the wanted side (the ring is edge-on or degenerate).
    if (!wanted(pts[0])) return;
    start = 0;
  }
  ctx.beginPath();
  for (let n = 0; n < pts.length; n++) {
    const p = pts[(start + n) % pts.length];
    if (!wanted(p) && n > 0) break;
    if (n === 0) ctx.moveTo(p.sx, p.sy);
    else ctx.lineTo(p.sx, p.sy);
  }
  ctx.stroke();
}

/** Twice the signed area of a polygon; its sign gives the winding direction. */
export function signedArea(pts) {
  let s = 0;
  for (let i = 0; i < pts.length; i++) {
    const a = pts[i];
    const b = pts[(i + 1) % pts.length];
    s += a.sx * b.sy - b.sx * a.sy;
  }
  return s;
}

/**
 * Append a subpath, flipping it if needed so its winding matches `refSign`.
 *
 * This is what lets a shape be built as a **union** of overlapping subpaths:
 * with a common winding, `fill("nonzero")` adds them together instead of
 * cancelling the overlaps into holes.
 */
export function emitSubpath(ctx, pts, refSign) {
  const seq = signedArea(pts) * refSign >= 0 ? pts : pts.slice().reverse();
  seq.forEach((p, i) => (i ? ctx.lineTo(p.sx, p.sy) : ctx.moveTo(p.sx, p.sy)));
}

// ---------------------------------------------------------------------------
// Specimens
// ---------------------------------------------------------------------------

/**
 * Draw atoms as glossy depth-cued spheres.
 *
 * The "3D ball" illusion is entirely a radial gradient with an off-centre
 * highlight plus a dark rim — no lighting model. Near atoms are drawn larger and
 * brighter, which is what reads as depth.
 *
 * @param {Array<{x,y,z, color?: [number,number,number], radius?: number}>} atoms
 * @param {object} [options]
 * @param {number} [options.radius=0.5] radius in **scene units** — so it can be
 *   set from a real atomic radius and stays consistent with the rest of the scene
 * @param {[number,number,number]} [options.color] default colour
 * @param {number} [options.fadeFar=0.55] brightness of the furthest atoms
 * @param {number} [options.depthScale=0.25] how much nearer atoms are enlarged
 * @param {number} [options.rim=0.45] alpha of the dark outline
 */
export function drawSpheres(ctx, view, atoms, options = {}) {
  const {
    radius: baseRadius = 0.5,
    color: baseColor = [150, 170, 200],
    fadeFar = 0.55,
    depthScale = 0.25,
    rim = 0.45,
    projectFn = project
  } = options;

  const items = atoms.map((a) => {
    const p = projectFn(a.x, a.y, a.z, view);
    return {...p, atom: a};
  });
  depthSort(items);
  const {t} = depthRange(items);

  for (const it of items) {
    const depthT = t(it.depth);
    const fade = fadeFar + (1 - fadeFar) * depthT;
    const [r0, g0, b0] = it.atom.color ?? baseColor;
    const cr = r0 * fade, cg = g0 * fade, cb = b0 * fade;
    // Radius is in scene units, so `view.zoom` converts it to pixels; nearer
    // atoms are then drawn slightly larger as a depth cue.
    const rad = baseRadius * (it.atom.radius ?? 1) * view.zoom *
      (1 - depthScale + 2 * depthScale * depthT);

    const grad = ctx.createRadialGradient(
      it.sx - rad * 0.38, it.sy - rad * 0.42, rad * 0.1,
      it.sx, it.sy, rad
    );
    grad.addColorStop(0, `rgb(${Math.min(255, cr * 1.55 + 40) | 0},${Math.min(255, cg * 1.55 + 40) | 0},${Math.min(255, cb * 1.55 + 40) | 0})`);
    grad.addColorStop(0.5, `rgb(${cr | 0},${cg | 0},${cb | 0})`);
    grad.addColorStop(1, `rgb(${(cr * 0.4) | 0},${(cg * 0.4) | 0},${(cb * 0.4) | 0})`);

    ctx.fillStyle = grad;
    ctx.strokeStyle = `rgba(0,0,0,${rim})`;
    ctx.lineWidth = 0.9;
    ctx.beginPath();
    ctx.arc(it.sx, it.sy, rad, 0, Math.PI * 2);
    ctx.fill();
    ctx.stroke();
  }
  return items;
}

/**
 * Draw bonds as round-capped lines. Call this *before* {@link drawSpheres} so
 * the spheres cover the bond ends.
 *
 * @param {Array<[number, number]>} bonds index pairs into `atoms`
 */
export function drawBonds(ctx, view, atoms, bonds, options = {}) {
  const {color = "rgba(116,130,160,0.62)", width = 2.1, projectFn = project} = options;
  ctx.save();
  ctx.lineCap = "round";
  ctx.strokeStyle = color;
  ctx.lineWidth = width;
  ctx.beginPath();
  for (const [i, j] of bonds) {
    const a = projectFn(atoms[i].x, atoms[i].y, atoms[i].z, view);
    const b = projectFn(atoms[j].x, atoms[j].y, atoms[j].z, view);
    ctx.moveTo(a.sx, a.sy);
    ctx.lineTo(b.sx, b.sy);
  }
  ctx.stroke();
  ctx.restore();
}

// ---------------------------------------------------------------------------
// Textured planes (detector images, potential slices)
// ---------------------------------------------------------------------------

/**
 * Draw an image onto the plane at height `z`, as an affine parallelogram.
 *
 * Fast and good enough whenever the tilt is mild, which it is for these scenes.
 * An affine map can only produce a parallelogram, so the far edge is not
 * foreshortened — use {@link drawPlaneBanded} when that shows.
 *
 * @param {CanvasImageSource} image
 * @param {number} half half-width of the plane in scene units
 * @returns {{sx: number, sy: number}} the drawn centre, which is generally NOT
 *   `project(x, y, z)` — pass it to {@link drawProbeCone} as `anchor` so the
 *   beam lands on the image rather than on the true projected point.
 */
export function drawPlane(ctx, view, image, x, y, z, half, options = {}) {
  const {alpha = 1, border, projectFn = project} = options;
  const c0 = projectFn(x - half, y - half, z, view);
  const c1 = projectFn(x + half, y - half, z, view);
  const c3 = projectFn(x - half, y + half, z, view);
  const c2 = {sx: c1.sx + c3.sx - c0.sx, sy: c1.sy + c3.sy - c0.sy};

  const w = image.width, h = image.height;
  ctx.save();
  ctx.globalAlpha = alpha;
  ctx.imageSmoothingEnabled = options.smooth ?? false;
  // `transform`, not `setTransform`: the context already carries a device-pixel
  // scale from `makeCanvas`, and replacing the matrix outright would draw the
  // image at 1/dpr in the wrong corner.
  ctx.save();
  ctx.transform(
    (c1.sx - c0.sx) / w, (c1.sy - c0.sy) / w,
    (c3.sx - c0.sx) / h, (c3.sy - c0.sy) / h,
    c0.sx, c0.sy
  );
  ctx.drawImage(image, 0, 0);
  ctx.restore();

  if (border) {
    ctx.globalAlpha = 1;
    ctx.strokeStyle = border;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(c0.sx, c0.sy);
    ctx.lineTo(c1.sx, c1.sy);
    ctx.lineTo(c2.sx, c2.sy);
    ctx.lineTo(c3.sx, c3.sy);
    ctx.closePath();
    ctx.stroke();
  }
  ctx.restore();

  return {sx: (c0.sx + c2.sx) / 2, sy: (c0.sy + c2.sy) / 2};
}

/**
 * Perspective-correct version of {@link drawPlane}: split the quad into
 * horizontal bands, each with its own affine map.
 *
 * This is the stacked-slices look in `ptycho-ms.js`. The residual trapezoid
 * error within a band is halved by averaging the top and bottom row vectors, and
 * bands overlap slightly in the source so no hairline seams appear.
 */
export function drawPlaneBanded(ctx, view, image, x, y, z, half, options = {}) {
  const {bands = 32, alpha = 1, border, projectFn = project} = options;
  const w = image.width, h = image.height;

  const left = [], right = [];
  for (let b = 0; b <= bands; b++) {
    const yy = y - half + (2 * half * b) / bands;
    left.push(projectFn(x - half, yy, z, view));
    right.push(projectFn(x + half, yy, z, view));
  }

  const srcBand = h / bands;
  ctx.save();
  ctx.globalAlpha = alpha;
  ctx.imageSmoothingEnabled = options.smooth ?? true;
  for (let b = 0; b < bands; b++) {
    const l0 = left[b], r0 = right[b], l1 = left[b + 1], r1 = right[b + 1];
    const a = ((r0.sx - l0.sx) + (r1.sx - l1.sx)) / 2 / w;
    const bb = ((r0.sy - l0.sy) + (r1.sy - l1.sy)) / 2 / w;
    const c = (l1.sx - l0.sx) / srcBand;
    const d = (l1.sy - l0.sy) / srcBand;
    const ox = l0.sx + ((r0.sx - l0.sx) - a * w) / 2;
    const oy = l0.sy + ((r0.sy - l0.sy) - bb * w) / 2;
    // See `drawPlane`: compose, never replace — the DPR scale lives here too.
    ctx.save();
    ctx.transform(a, bb, c, d, ox, oy);
    const overlap = b < bands - 1 ? 0.6 : 0;
    ctx.drawImage(image, 0, b * srcBand, w, srcBand + overlap, 0, 0, w, srcBand + overlap);
    ctx.restore();
  }

  if (border) {
    ctx.globalAlpha = 1;
    ctx.strokeStyle = border;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(left[0].sx, left[0].sy);
    ctx.lineTo(right[0].sx, right[0].sy);
    for (let b = 1; b <= bands; b++) ctx.lineTo(right[b].sx, right[b].sy);
    for (let b = bands; b >= 0; b--) ctx.lineTo(left[b].sx, left[b].sy);
    ctx.closePath();
    ctx.stroke();
  }
  ctx.restore();

  const mid = bands >> 1;
  return {sx: (left[mid].sx + right[mid].sx) / 2, sy: (left[mid].sy + right[mid].sy) / 2};
}

// ---------------------------------------------------------------------------
// The probe cone
// ---------------------------------------------------------------------------

/** Green tokens, matching the reference widgets. */
export const BEAM = {
  fill: "rgba(0,200,110,0.3)",
  stroke: "#00b863",
  bright: "#00ff88",
  accent: "#00cc66"
};

/**
 * Draw the converging/diverging electron beam.
 *
 * Call it **twice** around the specimen so the specimen occludes the beam — that
 * ordering is what creates the depth illusion:
 *
 *   drawProbeCone(ctx, view, {...opts, part: "lower"});   // sample -> detector
 *   drawSpheres(ctx, view, atoms);
 *   drawProbeCone(ctx, view, {...opts, part: "upper"});   // aperture -> sample
 *
 * The cone's radius is `|z - crossoverZ| * slope`, with the slope fixed by
 * `detectorRadius` at `zDetector`. So moving `crossoverZ` (which should track
 * defocus) moves the waist, and the beam footprint on the specimen follows.
 *
 * @param {object} opts
 * @param {number} [opts.x=0] probe position
 * @param {number} [opts.y=0]
 * @param {number} opts.zAperture top of the beam
 * @param {number} [opts.zSample=0] the specimen plane
 * @param {number} opts.zDetector bottom of the beam
 * @param {number} opts.detectorRadius beam radius at the detector, scene units
 * @param {number} [opts.crossoverZ=0] height of the waist; 0 means focused on the sample
 * @param {"lower"|"upper"|"both"} [opts.part="both"]
 * @param {(z: number) => {x: number, y: number}} [opts.axis] where the beam's
 *   axis sits at each height. The default is a vertical column through (x, y).
 *   Pass a function to model **descan**: below the specimen the post-specimen
 *   optics bring the pattern back onto the optic axis, so the lower half of the
 *   beam leans from the scan position back to the centre of the detector. That
 *   is also what makes an *angular* detector plane consistent — the pattern is
 *   centred on zero scattering angle wherever the probe happens to sit.
 * @param {(x,y,z,view) => {sx,sy}} [opts.projectFn] override, e.g. to align the
 *   beam with an image drawn by {@link drawPlane} (pass its returned centre).
 */
export function drawProbeCone(ctx, view, opts) {
  const {
    x = 0, y = 0,
    zAperture, zSample = 0, zDetector,
    detectorRadius,
    crossoverZ = 0,
    part = "both",
    fill = BEAM.fill,
    stroke = BEAM.stroke,
    lineWidth = 1.6,
    segments = 48,
    showCrossover = true,
    disks = 0,
    diskStroke = BEAM.stroke,
    axis = null,
    projectFn = project
  } = opts;

  const axisAt = axis ?? (() => ({x, y}));

  // The cone is pinned by its radius at the detector, so the beam always lands
  // on its own bright-field disc. Its half-angle then drifts as the waist moves,
  // which is why the widgets using this keep their defocus range short.
  const slope = Math.abs(zDetector - crossoverZ) > 1e-9
    ? detectorRadius / Math.abs(zDetector - crossoverZ)
    : 0;
  const radiusAt = (z) => Math.abs(z - crossoverZ) * slope;
  const P = (px, py, pz) => projectFn(px, py, pz, view);
  const edge = (z, sign) => {
    const a = axisAt(z);
    return P(a.x + sign * radiusAt(z), a.y, z);
  };

  /** Fill the frustum between two heights, splitting at the waist if it lies between. */
  function fillBetween(zTop, zBottom) {
    const stops = [zTop];
    if ((crossoverZ - zTop) * (crossoverZ - zBottom) < 0) stops.push(crossoverZ);
    stops.push(zBottom);

    const caps = [];
    for (const z of [zTop, zBottom]) {
      const r = radiusAt(z);
      const a = axisAt(z);
      if (r > 0) caps.push(ellipsePoints(view, a.x, a.y, z, r, segments, projectFn));
    }
    const reference = caps[0] ??
      ellipsePoints(view, axisAt(zTop).x, axisAt(zTop).y, zTop, 1, segments, projectFn);
    const refSign = signedArea(reference) >= 0 ? 1 : -1;

    ctx.fillStyle = fill;
    ctx.globalAlpha = 1;
    ctx.beginPath();
    for (const cap of caps) emitSubpath(ctx, cap, refSign);
    for (let i = 0; i < stops.length - 1; i++) {
      const z1 = stops[i], z2 = stops[i + 1];
      // A segment that touches the waist degenerates to a triangle, which is
      // exactly right — the two edges meet at a point.
      emitSubpath(ctx, [edge(z1, -1), edge(z1, +1), edge(z2, +1), edge(z2, -1)], refSign);
    }
    ctx.fill("nonzero");

    ctx.globalAlpha = 0.85;
    ctx.strokeStyle = stroke;
    ctx.lineWidth = lineWidth;
    ctx.beginPath();
    for (const sign of [-1, +1]) {
      const a = edge(zTop, sign);
      ctx.moveTo(a.sx, a.sy);
      for (let i = 1; i < stops.length; i++) {
        const p = edge(stops[i], sign);
        ctx.lineTo(p.sx, p.sy);
      }
    }
    ctx.stroke();

    ctx.strokeStyle = stroke;
    for (const z of [zTop, zBottom]) {
      const r = radiusAt(z);
      const a = axisAt(z);
      if (r > 0) ring(ctx, view, a.x, a.y, z, r, {segments, projectFn});
    }

    // Rings along the beam. The cone alone is ambiguous about which way it is
    // going; a ladder of rings reads as depth, and because the caller draws the
    // half below the specimen before the atoms and the half above after, the
    // rings pass behind and in front of the particle as the waist moves through
    // it — which is the whole point of showing defocus in 3D.
    if (disks > 0) {
      ctx.strokeStyle = diskStroke;
      ctx.globalAlpha = 0.5;
      ctx.lineWidth = 1;
      for (let d = 1; d <= disks; d++) {
        const z = zTop + ((zBottom - zTop) * d) / (disks + 1);
        const r = radiusAt(z);
        const a = axisAt(z);
        if (r > 0.02) ring(ctx, view, a.x, a.y, z, r, {segments, projectFn});
      }
      ctx.lineWidth = lineWidth;
    }
    ctx.globalAlpha = 1;
  }

  ctx.save();
  ctx.setLineDash([]);
  if (part === "lower" || part === "both") fillBetween(zSample, zDetector);
  if (part === "upper" || part === "both") {
    fillBetween(zAperture, zSample);
    if (showCrossover && crossoverZ > zSample && crossoverZ < zAperture) {
      const a = axisAt(crossoverZ);
      const c = P(a.x, a.y, crossoverZ);
      ctx.globalAlpha = 0.9;
      ctx.fillStyle = BEAM.bright;
      ctx.beginPath();
      ctx.arc(c.sx, c.sy, 2.5, 0, Math.PI * 2);
      ctx.fill();
    }
  }
  ctx.restore();
}

/** Beam radius at the specimen for a given crossover height — handy for captions. */
export function beamRadiusAt(z, {crossoverZ, zDetector, detectorRadius}) {
  const slope = Math.abs(zDetector - crossoverZ) > 1e-9
    ? detectorRadius / Math.abs(zDetector - crossoverZ)
    : 0;
  return Math.abs(z - crossoverZ) * slope;
}
