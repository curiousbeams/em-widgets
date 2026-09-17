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

/**
 * Ordered categorical colours, for when a widget needs to tell several things
 * apart.
 *
 * These are the colours already in use across the widgets — the accent blue,
 * the orange that marks a specimen bundle, the beam green, a violet — collected
 * so that the third series in one widget is the same hue as the third in
 * another. All four carry enough contrast against both a white page and a dark
 * one, which rules out yellows and pale greens however well they separate.
 */
export const SERIES = ["#04a1cd", "#e8833a", "#00b863", "#7b52d3"];

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
  /* Observable's Inputs carry an inline max-width of 640px. That never bites in
     a multi-column grid, where no cell is that wide, but a control given a row
     of its own then wraps at 640 with empty space beside it.
     (No backticks in this block: it lives inside a template literal.) */
  width: auto; max-width: none;
}
.em-controls--compact form[class^="inputs-"] > label {
  display: block; width: auto; padding: 0; margin: 0 0 1px;
  color: var(--theme-foreground-muted, inherit);
}
/* Wrapping matters for radios. A slider's row is one range plus one number and
   always fits, but a radio's row is one label per option and will happily run
   past the end of a 240px control column and over whatever is beside it —
   which for a canvas that handles its own pointer events means the options
   underneath stop being clickable. */
.em-controls--compact form[class^="inputs-"] > div {
  display: flex; gap: 5px; align-items: center; width: auto;
  flex-wrap: wrap; min-width: 0; max-width: 100%;
}
.em-controls--compact form[class^="inputs-"] input[type=range] {
  flex: 1 1 auto; min-width: 0; width: auto;
}
.em-controls--compact form[class^="inputs-"] input[type=number] {
  /* Wide enough for a four-decimal value plus the spinner. At 54px a 0.875
     came out as "0.8…". */
  flex: 0 0 auto; width: 66px;
}
.em-controls--compact form[class^="inputs-"] select { width: 100%; }

/* Rows: label beside its input, on one line, in a fixed-width column.
   Compact stacks the label above the input, which doubles the height of every
   control — fine in a panel, but not when each slider has to sit at the height
   of the component it drives and two components are 150 mm apart. */
.em-controls--rows form[class^="inputs-"] {
  display: flex; align-items: center; gap: 6px; margin: 0;
  font-size: 10px; line-height: 1.2; width: auto;
}
.em-controls--rows form[class^="inputs-"] > label {
  display: block; flex: 0 0 112px; width: 112px; padding: 0; margin: 0;
  white-space: nowrap; overflow: hidden; text-overflow: ellipsis;
  color: var(--theme-foreground-muted, inherit);
}
.em-controls--rows form[class^="inputs-"] > div {
  display: flex; gap: 5px; align-items: center; flex: 1 1 auto; min-width: 0;
}
.em-controls--rows form[class^="inputs-"] input[type=range] {
  flex: 1 1 auto; min-width: 0; width: auto;
}
.em-controls--rows form[class^="inputs-"] input[type=number] {
  flex: 0 0 auto; width: 58px; font-size: 10px; padding: 1px 3px;
}

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

/* A diagram of panels joined by labelled arrows. */
.em-flow { position: relative; display: grid; justify-content: center; }
.em-flow > svg {
  position: absolute; inset: 0; width: 100%; height: 100%;
  pointer-events: none; overflow: visible;
}
.em-flow-node { margin: 0; position: relative; }
.em-flow-node > figcaption {
  font-size: 10px; text-align: center; margin-top: 4px; line-height: 1.25;
  color: var(--theme-foreground-muted, inherit);
}
.em-flow-link-label {
  font-size: 10px; font-family: var(--monospace, ui-monospace, Menlo, monospace);
}

/* A two-dimensional control: magnitude and direction in one gesture. */
.em-pad { margin: 0 0 6px; display: block; width: max-content; }
.em-pad > .em-pad-label {
  display: block; font-size: 11px; line-height: 1.35; margin: 0 0 2px;
  color: var(--theme-foreground-muted, inherit);
  font-family: var(--monospace, ui-monospace, Menlo, monospace);
}
.em-pad > svg { display: block; touch-action: none; cursor: crosshair; }
.em-pad > .em-pad-readout {
  display: block; font-size: 10px; margin-top: 2px; text-align: center;
  font-variant-numeric: tabular-nums;
  font-family: var(--monospace, ui-monospace, Menlo, monospace);
  color: var(--theme-foreground-muted, inherit);
}

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
export function controls(form, {columns = 2, minWidth = 240, gap = 18, compact = false, rows = false} = {}) {
  ensureStyles(form);
  form.classList.add("em-controls");
  if (compact) form.classList.add("em-controls--compact");
  if (rows) form.classList.add("em-controls--rows");
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

let flowCount = 0;

/**
 * Labelled canvases on a grid, joined by arrows drawn over them.
 *
 * `panels` lays out a row of images and leaves the reader to supply the arrows.
 * An algorithm is a sequence of operations on those images, and the operations
 * are the part worth naming — so this places each panel on a grid and draws the
 * connection between them, with the operation written on the arrow.
 *
 * Arrow geometry is measured from the laid-out elements rather than computed
 * from the grid, so it survives reflow, and a `ResizeObserver` redraws it.
 *
 * A node is a canvas by default. Give it `content: element` instead of a width
 * and height and that element is placed in the cell as-is, which is how a plot
 * or a legend joins the diagram without being drawn on a canvas.
 *
 * @param {Array<{key: string, label?: string, width?: number, height?: number,
 *   content?: HTMLElement, row?: number, column?: number, rowSpan?: number,
 *   colSpan?: number, cursor?: string}>} nodes
 * @param {Array<{from: string, to: string, label?: string,
 *   route?: "auto"|"h"|"v", bias?: number, anchor?: "centre"|"overlap"}>} links
 *   `from` and `to` are node keys. `bias` moves the corner of an elbow along the
 *   run, 0 being hard against the start and 1 against the end; the default 0.5
 *   is the midpoint. `anchor: "overlap"` meets each panel at the middle of the
 *   span the two share rather than at its own centre, which is what a panel
 *   spanning several rows needs when more than one link touches it.
 * @param {object} [options]
 * @param {number} [options.gap=30] space between grid cells, which is where the
 *   arrows and their labels go
 * @param {{before: string, label?: string}} [options.divider] a dashed rule down
 *   the left edge of a node's column, for marking a change of domain
 * @returns {{node: HTMLElement, panels: Map<string, object>, redraw(): void,
 *   setLinkLabel(from: string, to: string, text: string): void}}
 */
export function flowDiagram(nodes, links, {gap = 30, divider = null} = {}) {
  const svgNS = "http://www.w3.org/2000/svg";
  const make = (name, attrs = {}) => {
    const node = document.createElementNS(svgNS, name);
    for (const [key, value] of Object.entries(attrs)) node.setAttribute(key, value);
    return node;
  };

  const root = document.createElement("div");
  root.className = "em-flow";
  root.style.gap = `${gap}px`;
  // Somewhere for a divider's label to sit. Without it the label is drawn over
  // the top edge of whichever panel the rule runs past.
  if (divider?.label) root.style.paddingTop = "15px";
  ensureStyles(root);

  const id = `em-flow-${++flowCount}`;
  const made = new Map();

  for (const spec of nodes) {
    const figure = document.createElement("figure");
    figure.className = "em-flow-node";
    if (spec.width != null) figure.style.width = `${spec.width}px`;
    if (spec.row != null) figure.style.gridRow = `${spec.row + 1} / span ${spec.rowSpan ?? 1}`;
    if (spec.column != null) figure.style.gridColumn = `${spec.column + 1} / span ${spec.colSpan ?? 1}`;

    let canvas = null;
    let ctx = null;
    if (spec.content) {
      figure.appendChild(spec.content);
    } else {
      ({canvas, ctx} = makeCanvas(spec.width, spec.height));
      canvas.style.borderRadius = "4px";
      if (spec.cursor) canvas.style.cursor = spec.cursor;
      figure.appendChild(canvas);
    }

    const caption = document.createElement("figcaption");
    caption.textContent = spec.label ?? "";
    if (!spec.label) caption.style.display = "none";
    figure.appendChild(caption);

    root.appendChild(figure);
    made.set(spec.key, {
      canvas, ctx, node: figure, content: spec.content ?? null,
      width: spec.width, height: spec.height,
      setLabel(text) {
        caption.textContent = text;
        caption.style.display = text ? "" : "none";
      }
    });
  }

  // The arrow layer goes last so it paints over the panels.
  const svg = make("svg");
  const defs = make("defs");
  const marker = make("marker", {
    id: `${id}-head`, viewBox: "0 0 8 8", refX: 7, refY: 4,
    markerWidth: 6, markerHeight: 6, orient: "auto-start-reverse"
  });
  marker.appendChild(make("path", {d: "M0,0 L8,4 L0,8 z", fill: "currentColor"}));
  defs.appendChild(marker);
  svg.appendChild(defs);
  root.appendChild(svg);

  const labels = new Map();

  const redraw = () => {
    const frame = root.getBoundingClientRect();
    if (frame.width === 0) return;
    svg.setAttribute("viewBox", `0 0 ${frame.width} ${frame.height}`);
    while (svg.lastChild !== defs) svg.removeChild(svg.lastChild);

    // Horizontal links anchor on the canvas, so a caption under one panel and
    // not another does not pull the arrow off centre. Vertical links anchor on
    // the figure, because a caption sits between the canvas and the arrow and
    // the line would otherwise be drawn straight through the words.
    const box = (key, whole = false) => {
      const panel = made.get(key);
      const r = ((whole ? panel.node : panel.canvas) ?? panel.node).getBoundingClientRect();
      return {
        left: r.left - frame.left, right: r.right - frame.left,
        top: r.top - frame.top, bottom: r.bottom - frame.top,
        cx: (r.left + r.right) / 2 - frame.left,
        cy: (r.top + r.bottom) / 2 - frame.top
      };
    };

    if (divider) {
      const at = box(divider.before).left - gap / 2;
      svg.appendChild(make("line", {
        x1: at, y1: 0, x2: at, y2: frame.height,
        stroke: "currentColor", "stroke-opacity": 0.28,
        "stroke-width": 1, "stroke-dasharray": "5 4"
      }));
      if (divider.label) {
        const text = make("text", {
          x: at + 5, y: 9, "font-size": 10, fill: "currentColor",
          "fill-opacity": 0.55, class: "em-flow-link-label"
        });
        text.textContent = divider.label;
        svg.appendChild(text);
      }
    }

    for (const link of links) {
      const route = link.route ?? "auto";
      const horizontal = route === "h"
        || (route === "auto" && box(link.to).left >= box(link.from).right - 1);
      const a = box(link.from, !horizontal);
      const b = box(link.to, !horizontal);

      const bias = link.bias ?? 0.5;
      // Where a link meets a panel that spans more rows or columns than the one
      // at the other end. Centre to centre is the default and is right when the
      // two are the same size. It is wrong when a tall panel has two links, one
      // from a panel above and one from a panel below: both would arrive at the
      // same point on its edge, and the two paths would run along each other to
      // get there. `overlap` aims at the middle of the range the two panels
      // share instead, which makes each link a straight line to the panel it
      // actually comes from.
      const shared = (lo, hi) => (link.anchor === "overlap" && lo < hi ? (lo + hi) / 2 : null);
      let start;
      let end;
      let mid;
      if (horizontal) {
        const forward = b.cx >= a.cx;
        const y = shared(Math.max(a.top, b.top), Math.min(a.bottom, b.bottom));
        start = {x: forward ? a.right : a.left, y: y ?? a.cy};
        end = {x: forward ? b.left : b.right, y: y ?? b.cy};
      } else {
        const down = b.cy >= a.cy;
        const x = shared(Math.max(a.left, b.left), Math.min(a.right, b.right));
        start = {x: x ?? a.cx, y: down ? a.bottom : a.top};
        end = {x: x ?? b.cx, y: down ? b.top : b.bottom};
      }
      mid = {
        x: start.x + bias * (end.x - start.x),
        y: start.y + bias * (end.y - start.y)
      };

      // A straight line when the two are aligned, an elbow when they are not,
      // so an arrow never cuts diagonally across a panel between them.
      const aligned = horizontal
        ? Math.abs(start.y - end.y) < 2
        : Math.abs(start.x - end.x) < 2;
      const d = aligned
        ? `M ${start.x} ${start.y} L ${end.x} ${end.y}`
        : horizontal
          ? `M ${start.x} ${start.y} L ${mid.x} ${start.y} L ${mid.x} ${end.y} L ${end.x} ${end.y}`
          : `M ${start.x} ${start.y} L ${start.x} ${mid.y} L ${end.x} ${mid.y} L ${end.x} ${end.y}`;

      svg.appendChild(make("path", {
        d, fill: "none", stroke: "currentColor", "stroke-opacity": 0.5,
        "stroke-width": 1.2, "marker-end": `url(#${id}-head)`
      }));

      const text = labels.get(`${link.from}>${link.to}`) ?? link.label;
      if (text) {
        // Beside a straight vertical line, above a horizontal one or an elbow.
        // A vertical run has no side to sit above: centring the label on it puts
        // the line through the middle of the words.
        const beside = aligned && !horizontal;
        const node = make("text", {
          x: beside ? mid.x + 6 : mid.x,
          y: beside ? mid.y : mid.y - 5,
          "text-anchor": beside ? "start" : "middle",
          "dominant-baseline": beside ? "middle" : "auto",
          fill: "currentColor", "fill-opacity": 0.75,
          "paint-order": "stroke",
          stroke: "var(--theme-background, #fff)", "stroke-width": 3,
          "stroke-linejoin": "round", class: "em-flow-link-label"
        });
        node.textContent = text;
        svg.appendChild(node);
      }
    }
  };

  // Lay out first, then measure. A cell builds its DOM before Observable
  // attaches it, so every rect is zero until the next frame.
  requestAnimationFrame(redraw);
  if (typeof ResizeObserver === "function") new ResizeObserver(redraw).observe(root);

  return {
    node: root,
    panels: made,
    redraw,
    setLinkLabel(from, to, text) {
      labels.set(`${from}>${to}`, text);
      redraw();
    }
  };
}

/**
 * A two-dimensional input: drag inside a disc to set a magnitude and a
 * direction at once.
 *
 * Aberrations like astigmatism and coma are one physical quantity with a size
 * and an axis, and splitting them across two linear sliders makes the reader
 * reassemble them. Worse, a slider labelled "angle" leaves the mapping between
 * its number and the direction on screen entirely unstated.
 *
 * So the pad is drawn as a **miniature of the panel it acts on**: the first
 * array axis runs down it and the second across it, exactly as an image of a
 * corner-centered grid is drawn. Where the handle sits is where the aberration
 * points, and `angle` is `atan2(second, first)` — the `phi` convention of
 * `grid.polarCoordinates` — so it can be handed straight to `chi` as a
 * `phi_nm`.
 *
 * Behaves like an Observable input: it has a `.value` and emits `input`, so
 * `view(vectorPad(...))` and `Inputs.form({astigmatism: vectorPad(...)})` both
 * work.
 *
 * @param {object} [options]
 * @param {string} [options.label] shown above the disc
 * @param {number} [options.size=104] the disc's box, in CSS pixels
 * @param {number} [options.max=1] magnitude at the rim
 * @param {number} [options.magnitude=0] initial magnitude
 * @param {number} [options.angle=0] initial angle [degrees]
 * @param {number} [options.fold=1] rotational order of the aberration. Above
 *   one, the directions that give the same aberration are drawn as ghosts, so
 *   the degeneracy is visible rather than a surprise.
 * @param {string} [options.unit=""] appended to the magnitude in the readout
 * @param {number} [options.digits=1] decimals in the readout
 * @returns {HTMLElement} with `.value = {magnitude, angle}`
 */
export function vectorPad({
  label = "",
  size = 104,
  max = 1,
  magnitude = 0,
  angle = 0,
  fold = 1,
  unit = "",
  digits = 1
} = {}) {
  const svgNS = "http://www.w3.org/2000/svg";
  const make = (name, attrs) => {
    const node = document.createElementNS(svgNS, name);
    for (const [key, value] of Object.entries(attrs)) node.setAttribute(key, value);
    return node;
  };

  const centre = size / 2;
  const rim = centre - 7;

  const root = document.createElement("div");
  root.className = "em-pad";
  ensureStyles(root);

  if (label) {
    const text = document.createElement("span");
    text.className = "em-pad-label";
    text.textContent = label;
    root.appendChild(text);
  }

  const svg = make("svg", {width: size, height: size, viewBox: `0 0 ${size} ${size}`});
  const faint = "var(--theme-foreground-faintest, rgba(128,128,128,0.3))";

  // Two guide rings and a pair of axes, so a magnitude can be judged by eye.
  for (const fraction of [1, 0.5]) {
    svg.appendChild(make("circle", {
      cx: centre, cy: centre, r: rim * fraction,
      fill: "none", stroke: faint, "stroke-width": 1,
      "stroke-dasharray": fraction === 1 ? "none" : "2 3"
    }));
  }
  svg.appendChild(make("line", {
    x1: centre - rim, y1: centre, x2: centre + rim, y2: centre,
    stroke: faint, "stroke-width": 1
  }));
  svg.appendChild(make("line", {
    x1: centre, y1: centre - rim, x2: centre, y2: centre + rim,
    stroke: faint, "stroke-width": 1
  }));

  // Ghost handles first, so the live one draws over them.
  const ghosts = [];
  for (let i = 1; i < fold; i++) {
    const ghost = make("circle", {
      cx: centre, cy: centre, r: 4, fill: "none",
      stroke: UI.accent, "stroke-width": 1.5, "stroke-opacity": 0.4
    });
    ghosts.push(ghost);
    svg.appendChild(ghost);
  }

  const stem = make("line", {
    x1: centre, y1: centre, x2: centre, y2: centre,
    stroke: UI.accent, "stroke-width": 1.5, "stroke-opacity": 0.6
  });
  const handle = make("circle", {
    cx: centre, cy: centre, r: 5.5, fill: UI.accent, "fill-opacity": 0.85,
    stroke: "var(--theme-background, #fff)", "stroke-width": 1.5
  });
  svg.appendChild(stem);
  svg.appendChild(handle);
  root.appendChild(svg);

  const readout = document.createElement("span");
  readout.className = "em-pad-readout";
  root.appendChild(readout);

  // `angle` is measured in the phi convention: from the first axis, which runs
  // down the pad, towards the second, which runs across it.
  const place = (node, m, degrees) => {
    const radians = (degrees * Math.PI) / 180;
    const radius = (Math.min(m, max) / max) * rim;
    node.setAttribute("cx", centre + radius * Math.sin(radians));
    node.setAttribute("cy", centre + radius * Math.cos(radians));
  };

  const render = () => {
    const {magnitude: m, angle: a} = root.value;
    place(handle, m, a);
    stem.setAttribute("x2", handle.getAttribute("cx"));
    stem.setAttribute("y2", handle.getAttribute("cy"));
    ghosts.forEach((ghost, i) => place(ghost, m, a + ((i + 1) * 360) / fold));
    readout.textContent = `${m.toFixed(digits)}${unit} · ${Math.round(a)}°`;
  };

  root.value = {magnitude, angle};
  render();

  const pick = (event) => {
    const box = svg.getBoundingClientRect();
    const across = event.clientX - box.left - centre;
    const down = event.clientY - box.top - centre;
    const radius = Math.hypot(across, down);
    root.value = {
      magnitude: Math.min(1, radius / rim) * max,
      // Wrapped into [0, 360) so a value read back never carries a sign that
      // depends on which way the pointer came in.
      angle: ((Math.atan2(across, down) * 180) / Math.PI + 360) % 360
    };
    render();
    root.dispatchEvent(new CustomEvent("input", {bubbles: true}));
  };

  svg.addEventListener("pointerdown", (event) => {
    svg.setPointerCapture(event.pointerId);
    pick(event);
  });
  svg.addEventListener("pointermove", (event) => {
    if (event.buttons > 0) pick(event);
  });

  return root;
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
  //
  // Both labels are always laid out, stacked in one grid cell, with the inactive
  // one merely invisible. That keeps the button exactly one size: swapping the
  // text otherwise resizes it, and `play` and `pause` do not even share a line
  // height because their glyphs come from different fonts — enough that toggling
  // one visibly nudged the whole widget below it down the page.
  const label = document.createElement("span");
  label.style.display = "grid";
  label.style.alignItems = "center";
  label.style.justifyItems = "center";

  const faces = [off, on].map((text) => {
    const face = document.createElement("span");
    face.style.gridArea = "1 / 1";
    face.textContent = text;
    label.appendChild(face);
    return face;
  });
  button.appendChild(label);

  function paint() {
    faces[0].style.visibility = state ? "hidden" : "visible";
    faces[1].style.visibility = state ? "visible" : "hidden";
    button.setAttribute("aria-pressed", String(state));
    // For assistive technology and for anything reading the button's text.
    button.setAttribute("aria-label", state ? on : off);
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
