// Widget chrome: the layout and controls that make a notebook look like a
// designed widget rather than a stack of matplotlib panels.
//
// Two layouts:
//
//   scene()   ONE canvas holding the whole thing — specimen, beam, detector.
//             This is what lets a beam visually connect the specimen to the
//             detector, and it is the main reason the reference widgets read as
//             instruments rather than as subplots. Use it for anything with a
//             column or a beam.
//
//   panels()  small labelled canvases, in a row or a grid. Right for analysis
//             views — a CTF next to its PSF — where the panels are separate.
//
// `row()` and `column()` compose the two into a figure.
//
// Everything here is plain DOM; no htl, no framework, so the kit stays usable
// outside Observable too.

// ---------------------------------------------------------------------------
// Palette
// ---------------------------------------------------------------------------

/**
 * The lab's colours, from the Curious Beams mark.
 *
 * Chrome uses these, not `BEAM` in `scene3d.js`. The beam is green because that
 * is how an electron beam is drawn; while the sliders took their accent from it,
 * restyling the site meant editing the beam.
 */
export const UI = {
  accent: "#04a1cd",
  bright: "#00ccff",
  deep: "#0f2a44"
};

// ---------------------------------------------------------------------------
// Stylesheet
// ---------------------------------------------------------------------------

const styled = new WeakSet();
const rootStyled = new WeakSet();

const STYLES = `
.em-scene, .em-panels { margin: 0; }
.em-title {
  font-size: 15px; font-weight: 600; text-align: center;
  margin: 0 0 6px; color: var(--theme-foreground-alt, inherit);
}
.em-canvas {
  display: block; width: 100%; height: auto;
  border-radius: 6px; touch-action: none;
  border: 1px solid var(--theme-foreground-faintest, rgba(128,128,128,0.25));
}
.em-caption {
  font-size: 12px; text-align: center; margin-top: 5px;
  color: var(--theme-foreground-muted, inherit); opacity: 0.85;
}
.em-hint {
  font-size: 11px; text-align: center; margin-top: 4px;
  color: var(--theme-foreground-faint, inherit); opacity: 0.7;
}
.em-panels-row {
  display: flex; gap: 10px; flex-wrap: wrap; justify-content: center;
  align-items: flex-start;
}
.em-panels-grid { display: grid; gap: 10px; justify-content: center; }
.em-panel { margin: 0; position: relative; }
/* A label laid over the image instead of under it, for compact layouts. White
   with a dark halo, so it holds up over a light panel and a dark one alike. */
.em-panel-label {
  position: absolute; top: 5px; left: 7px; pointer-events: none;
  font-size: 10px; font-weight: 600; letter-spacing: 0.03em; color: #fff;
  text-shadow: 0 0 3px rgba(0,0,0,0.95), 0 0 7px rgba(0,0,0,0.7);
}
.em-row {
  display: flex; gap: 12px; flex-wrap: wrap; align-items: flex-start;
  justify-content: center;
}
.em-column { display: flex; flex-direction: column; gap: 10px; align-items: center; }

/* Controls. Observable's Inputs are native <input type=range>, so styling is
   mostly one accent-color. */
.em-controls {
  display: grid; column-gap: 18px; row-gap: 2px; align-items: center;
  margin: 8px 0 10px;
}
.em-controls input[type=range] { accent-color: ${UI.accent}; }
.em-controls form[class^="inputs-"] { font: 13px var(--monospace, ui-monospace, Menlo, monospace); }

/* Compact: the label above its input rather than beside it.
   Observable's Inputs reserve 120px for the label, which is most of a narrow
   column and leaves the slider a stub. Stacking gives the slider the full
   width, which is the part the reader actually drags. */
.em-controls--compact form[class^="inputs-"] {
  display: block; font-size: 11px; line-height: 1.35; margin: 0 0 5px;
}
.em-controls--compact form[class^="inputs-"] > label {
  display: block; width: auto; padding: 0; margin: 0 0 1px;
  color: var(--theme-foreground-muted, inherit);
}
.em-controls--compact form[class^="inputs-"] > div {
  display: flex; gap: 5px; align-items: center; width: auto;
}
.em-controls--compact form[class^="inputs-"] input[type=range] {
  flex: 1 1 auto; min-width: 0; width: auto;
}
.em-controls--compact form[class^="inputs-"] input[type=number] {
  flex: 0 0 auto; width: 54px;
}
.em-controls--compact form[class^="inputs-"] select { width: 100%; }

/* Inline: label beside its input, for a single control sharing a line with
   buttons. Compact stacks the label above the input, which is right for a
   slider and wrong for a lone checkbox. */
/* max-content on both: Observable's own form styling gives it a width, and a
   flex item inherits that rather than shrinking to its label and checkbox —
   a 110px control then reserves 360px and pushes its neighbours onto a new
   line. */
.em-controls--inline { display: inline-block; margin: 0; width: max-content; }
.em-controls--inline form[class^="inputs-"] {
  display: flex; align-items: center; gap: 7px; margin: 0; width: max-content;
  font-size: 11px; white-space: nowrap;
}
.em-controls--inline form[class^="inputs-"] > label {
  display: inline; width: auto; padding: 0; margin: 0;
  color: var(--theme-foreground-muted, inherit);
}
.em-controls--inline form[class^="inputs-"] > div { display: flex; width: auto; }

.em-toggle {
  padding: 5px 13px; border-radius: 6px; cursor: pointer;
  font-size: 12px; font-weight: 600; font-family: inherit;
  background: transparent; color: inherit;
  border: 1px solid var(--theme-foreground-faintest, rgba(128,128,128,0.35));
  transition: background 120ms, border-color 120ms, color 120ms;
}
.em-toggle:hover { border-color: ${UI.accent}; }
.em-toggle[aria-pressed="true"] {
  background: ${UI.deep}; border-color: ${UI.accent}; color: ${UI.bright};
}
.em-collapsible { margin: 4px 0 10px; }
.em-collapsible > summary {
  cursor: pointer; font-size: 12px; font-weight: 600; padding: 3px 0;
  color: var(--theme-foreground-muted, inherit); list-style-position: outside;
}
.em-collapsible > summary:hover { color: var(--theme-foreground, inherit); }
.em-collapsible[open] > summary { margin-bottom: 4px; }
.em-badge {
  font-size: 11px; font-variant-numeric: tabular-nums;
  color: var(--theme-foreground-muted, inherit); opacity: 0.75;
}
@media (prefers-reduced-motion: reduce) {
  .em-toggle { transition: none; }
}
`;

/**
 * Attach the kit stylesheet to a component, once per component.
 *
 * The style element goes *inside* the node and is hoisted to the shadow root on
 * the next frame. It cannot start there: a cell builds its DOM before Observable
 * attaches it, so `getRootNode()` returns the detached node and the styles would
 * land in `document.head`, which a shadow root cannot see. Hoisting afterwards
 * keeps it out of a button's `textContent` and leaves one stylesheet per widget.
 */
export function ensureStyles(node) {
  if (!node || styled.has(node)) return;
  styled.add(node);
  const style = document.createElement("style");
  style.textContent = STYLES;
  node.appendChild(style);

  if (typeof requestAnimationFrame !== "function") return;
  requestAnimationFrame(() => {
    const root = node.getRootNode();
    if (!(typeof ShadowRoot === "function" && root instanceof ShadowRoot)) return;
    if (rootStyled.has(root)) style.remove();
    else { rootStyled.add(root); root.appendChild(style); }
  });
}

// ---------------------------------------------------------------------------
// Canvas helpers
// ---------------------------------------------------------------------------

/**
 * Create a canvas whose backing store is sized for the display.
 *
 * Returns the canvas plus the context already scaled, so all drawing can be done
 * in CSS pixels and stays crisp on retina.
 */
export function makeCanvas(cssWidth, cssHeight, {dpr = globalThis.devicePixelRatio ?? 1, maxDpr = 2} = {}) {
  const ratio = Math.min(dpr, maxDpr);
  const canvas = document.createElement("canvas");
  canvas.className = "em-canvas";
  canvas.width = Math.round(cssWidth * ratio);
  canvas.height = Math.round(cssHeight * ratio);
  canvas.style.width = `${cssWidth}px`;
  canvas.style.height = `${cssHeight}px`;
  const ctx = canvas.getContext("2d");
  ctx.scale(ratio, ratio);
  return {canvas, ctx, width: cssWidth, height: cssHeight, ratio};
}

// ---------------------------------------------------------------------------
// Layouts
// ---------------------------------------------------------------------------

/**
 * One canvas for the whole scene, with optional title, caption and hint.
 *
 * @param {object} options
 * @param {number} options.width CSS width
 * @param {number} options.height CSS height
 * @param {string} [options.title]
 * @param {string} [options.caption]
 * @param {string} [options.hint] a quieter line for interaction instructions
 * @param {string} [options.cursor="crosshair"]
 * @returns {{node: HTMLElement, canvas: HTMLCanvasElement, ctx: CanvasRenderingContext2D,
 *           width: number, height: number, setCaption(text: string): void}}
 */
export function scene({width, height, title, caption, hint, cursor = "crosshair"}) {
  const {canvas, ctx} = makeCanvas(width, height);
  canvas.style.cursor = cursor;

  const node = document.createElement("figure");
  node.className = "em-scene";

  if (title) {
    const t = document.createElement("div");
    t.className = "em-title";
    t.textContent = title;
    node.appendChild(t);
  }
  node.appendChild(canvas);

  const cap = document.createElement("figcaption");
  cap.className = "em-caption";
  if (caption) cap.textContent = caption;
  else cap.style.display = "none";
  node.appendChild(cap);

  if (hint) {
    const h = document.createElement("div");
    h.className = "em-hint";
    h.textContent = hint;
    node.appendChild(h);
  }
  ensureStyles(node);

  return {
    node, canvas, ctx, width, height,
    setCaption(text) {
      cap.textContent = text;
      cap.style.display = text ? "" : "none";
    }
  };
}

/**
 * A row of labelled canvases, for analysis views.
 *
 * @param {Array<{label?: string, width: number, height: number}>} specs
 * @returns {{node: HTMLElement, panels: Array<{canvas, ctx, node, setLabel}>}}
 */
export function panels(specs, {columns = 0, inset = false} = {}) {
  const row = document.createElement("div");
  if (columns > 0) {
    row.className = "em-panels-grid";
    row.style.gridTemplateColumns = `repeat(${columns}, max-content)`;
  } else {
    row.className = "em-panels-row";
  }

  const made = specs.map((spec) => {
    const {canvas, ctx} = makeCanvas(spec.width, spec.height);
    if (spec.cursor) canvas.style.cursor = spec.cursor;
    const fig = document.createElement("figure");
    fig.className = "em-panel";
    fig.style.width = `${spec.width}px`;
    fig.appendChild(canvas);

    const cap = document.createElement(inset ? "div" : "figcaption");
    cap.className = inset ? "em-panel-label" : "em-caption";
    cap.textContent = spec.label ?? "";
    if (!spec.label) cap.style.display = "none";
    fig.appendChild(cap);

    row.appendChild(fig);
    return {
      canvas, ctx, node: fig, width: spec.width, height: spec.height,
      setLabel(text) {
        cap.textContent = text;
        cap.style.display = text ? "" : "none";
      }
    };
  });

  ensureStyles(row);
  return {node: row, panels: made};
}

// ---------------------------------------------------------------------------
// Controls
// ---------------------------------------------------------------------------

/**
 * Lay an Observable `Inputs.form` out as a responsive grid and apply the kit's
 * control styling.
 *
 * Supersedes `canvas.gridForm`, which stays for the existing notebooks.
 *
 * @param {HTMLElement} form the element returned by `Inputs.form`
 * @param {object} [options]
 * @param {number} [options.columns=2] columns when there is room
 * @param {number} [options.minWidth=240] below this, use fewer columns
 */
export function controls(form, {columns = 2, minWidth = 240, gap = 18, compact = false} = {}) {
  ensureStyles(form);
  form.classList.add("em-controls");
  if (compact) form.classList.add("em-controls--compact");
  // `minmax(0, ...)` at the floor when no minimum is asked for, so a control
  // grid inside a narrow column can shrink instead of overflowing it.
  const track = minWidth > 0
    ? `max(${minWidth}px, calc((100% - ${(columns - 1) * gap}px) / ${columns}))`
    : `calc((100% - ${(columns - 1) * gap}px) / ${columns})`;
  form.style.gridTemplateColumns = `repeat(auto-fit, minmax(${track}, 1fr))`;
  if (minWidth === 0) form.style.minWidth = "0";
  form.style.columnGap = `${gap}px`;
  return form;
}

/**
 * A toggle button that swaps its own label, e.g. `▶ Scan` ⇄ `■ Stop`.
 *
 * Behaves like an Observable input: it has a `.value` and emits `input` events,
 * so `view(toggleButton(...))` works.
 *
 * @param {object} options
 * @param {string} options.on label shown while the value is true
 * @param {string} options.off label shown while false
 * @param {boolean} [options.value=false] initial state
 */
export function toggleButton({on, off, value = false}) {
  const button = document.createElement("button");
  button.type = "button";
  button.className = "em-toggle";
  ensureStyles(button);

  // `value` must be a real boolean. HTMLButtonElement.value is a reflected
  // string attribute, so a plain assignment stores "false" — which is truthy,
  // and the toggle silently never turns anything off. Shadow it with an own
  // property so `view(toggleButton(...))` yields a boolean.
  let state = Boolean(value);
  Object.defineProperty(button, "value", {
    get: () => state,
    set: (v) => { state = Boolean(v); paint(); },
    configurable: true
  });

  // The label lives in its own span: `ensureStyles` puts a <style> element
  // inside the button, and setting `textContent` on the button would throw it
  // away on the first repaint.
  const label = document.createElement("span");
  button.appendChild(label);

  function paint() {
    label.textContent = state ? on : off;
    button.setAttribute("aria-pressed", String(state));
  }
  paint();

  button.addEventListener("click", () => {
    button.value = !button.value;
    button.dispatchEvent(new CustomEvent("input", {bubbles: true}));
  });
  return button;
}

/**
 * Wrap a node in a `<details>`, for controls that matter but should not crowd
 * the ones a reader reaches for first.
 *
 * @param {string} summary the always-visible label
 * @param {HTMLElement} node
 * @param {boolean} [open=false]
 */
export function collapsible(summary, node, open = false) {
  const details = document.createElement("details");
  details.className = "em-collapsible";
  details.open = open;
  const label = document.createElement("summary");
  label.textContent = summary;
  details.append(label, node);
  ensureStyles(details);

  // Forward the wrapped node's value, so `view(collapsible("specimen", form))`
  // behaves like `view(form)`.
  //
  // Not cosmetic: Observable's `Generators.input` only emits an initial value
  // when `element.value !== undefined`. A wrapper without one therefore never
  // yields at all — the cell hangs, silently and forever, taking every cell
  // downstream of it with it. `input` events already bubble out of the form, so
  // the getter is the whole of what is missing.
  if (node && "value" in node) {
    Object.defineProperty(details, "value", {
      get: () => node.value,
      set: (v) => { node.value = v; },
      configurable: true
    });
  }
  return details;
}

/**
 * Put nodes side by side, or one above another.
 *
 * Widgets that grew past a single canvas need somewhere to say "these belong
 * together" without each notebook hand-rolling a flexbox.
 */
export function row(nodes, {gap = 12, align = "flex-start"} = {}) {
  const el = document.createElement("div");
  el.className = "em-row";
  el.style.gap = `${gap}px`;
  el.style.alignItems = align;
  for (const n of nodes) if (n) el.appendChild(n);
  ensureStyles(el);
  return el;
}

/** {@link row}, stacked vertically. */
export function column(nodes, {gap = 10, align = "center"} = {}) {
  const el = document.createElement("div");
  el.className = "em-column";
  el.style.gap = `${gap}px`;
  el.style.alignItems = align;
  for (const n of nodes) if (n) el.appendChild(n);
  ensureStyles(el);
  return el;
}

/** A small right-aligned readout, for frame timings and derived quantities. */
export function badge(text = "") {
  const el = document.createElement("div");
  el.className = "em-badge";
  ensureStyles(el);
  el.textContent = text;
  return el;
}

/**
 * Run `cheap` on every `input` event and `expensive` only on `change`, i.e. when
 * the reader releases the slider.
 *
 * The reference widgets use this to keep dragging smooth: the geometry redraws
 * live, the simulation is rebuilt once at the end.
 */
export function liveAndSettled(input, {cheap, expensive}) {
  if (cheap) input.addEventListener("input", () => cheap(input.value));
  if (expensive) input.addEventListener("change", () => expensive(input.value));
  return input;
}
