// observable-notebook.mjs
//
// A MyST {anywidget} that renders an Observable Notebook 2.0 (.html) source file
// entirely client-side: fetch -> deserialize -> transpile -> run.
//
// Usage from markdown:
//
//   :::{anywidget} https://curiousbeams.github.io/em-widgets/observable-notebook.mjs
//   {
//     notebook: "https://curiousbeams.github.io/em-widgets/notebooks/sem-ray-diagram.html",
//     cells: ["rayDiagram"],      // optional: cell ids or declared names to display
//     params: {semiangle: 0.15},  // optional: redefine notebook cells from markdown
//     hideCode: true              // optional: suppress echo of `pinned` cells (default true)
//   }
//   :::
//
// The body is JSON5, so comments, trailing commas and unquoted keys are fine.
//
// Notes for future editors:
//  - MyST's theme ALWAYS renders anywidgets inside an open shadow root, so every
//    stylesheet has to be injected into `el`. Anything appended to document.head
//    is invisible here.
//  - The MyST model shim only implements get/set/on. Calling model.off(), .send()
//    or .save_changes() throws, so cleanup must not touch them.
//  - Always use jsDelivr `/+esm` URLs. The raw dist/ files contain bare
//    `import "./inputs.css"` statements that a browser cannot load.

// Pinned deliberately: Notebooks 2.0 is still a technology preview, so we upgrade
// on purpose rather than tracking `latest`.
import {
  deserialize,
  transpile,
  resolveImportDefault
} from "https://cdn.jsdelivr.net/npm/@observablehq/notebook-kit@2.5.6/+esm";
import {
  NotebookRuntime,
  library,
  registerFile
} from "https://cdn.jsdelivr.net/npm/@observablehq/notebook-kit@2.5.6/dist/src/runtime/index.js/+esm";

// The generators namespace is not a named export of the runtime module — it is
// only reachable through the builtin thunk on `library`.
const Generators = library.Generators();

const INPUTS_CSS = "https://cdn.jsdelivr.net/npm/@observablehq/inputs/dist/index.css";

// Observable's whole palette derives from --theme-foreground and --theme-background-a
// via color-mix, so we only pin those two and inherit the rest of the chain. We
// deliberately do NOT load notebook-kit's global.css: it targets :root/html/body/main
// (which match nothing inside a shadow root) and would override the site's typography.
const THEME_CSS = `
@import url("${INPUTS_CSS}");

:host, .em-notebook {
  --theme-foreground: #1b1e23;
  --theme-background-a: #ffffff;
  --theme-foreground-focus: #3b5fc0;
  --theme-error: #e7040f;
  color-scheme: light;
}
.em-notebook[data-theme="dark"] {
  --theme-foreground: #f5f7fa;
  --theme-background-a: #1b1e23;
  --theme-foreground-focus: #8ba7ff;
  color-scheme: dark;
}
.em-notebook {
  --theme-background-b: color-mix(in srgb, var(--theme-foreground) 4%, var(--theme-background-a));
  --theme-background: var(--theme-background-a);
  --theme-background-alt: var(--theme-background-b);
  --theme-foreground-alt: color-mix(in srgb, var(--theme-foreground) 90%, var(--theme-background-a));
  --theme-foreground-muted: color-mix(in srgb, var(--theme-foreground) 60%, var(--theme-background-a));
  --theme-foreground-faint: color-mix(in srgb, var(--theme-foreground) 50%, var(--theme-background-a));
  --theme-foreground-fainter: color-mix(in srgb, var(--theme-foreground) 30%, var(--theme-background-a));
  --theme-foreground-faintest: color-mix(in srgb, var(--theme-foreground) 14%, var(--theme-background-a));
  --monospace: ui-monospace, Menlo, Consolas, monospace;

  color: var(--theme-foreground);
  font: inherit;              /* inherit the host site's typography */
  line-height: 1.5;
}
.em-notebook--cell:empty { min-height: 0; }
.em-notebook--cell { max-width: 100%; }
.em-notebook svg, .em-notebook canvas, .em-notebook img { max-width: 100%; }

/* Observable Inputs read these; keep them legible against the host background. */
.em-notebook form[class^="inputs-"] { font: 13px var(--monospace); color: var(--theme-foreground); }
.em-notebook figure { margin: 0; }

.em-notebook--error {
  color: var(--theme-error);
  font: 13px var(--monospace);
  white-space: pre-wrap;
  padding: 8px 10px;
  border: 1px solid currentColor;
  border-radius: 4px;
}
.em-notebook--code {
  font: 12px var(--monospace);
  overflow-x: auto;
  padding: 8px 10px;
  background: var(--theme-background-alt);
  border-radius: 4px;
}
`;

// Globals notebook-kit does not recognise but that widgets reasonably reach for.
// Providing them makes the site strictly more permissive than Observable Desktop,
// so anything that runs in Desktop also runs here (never the other way round).
const EXTRA_GLOBALS = {
  Event: () => Event,
  ResizeObserver: () => ResizeObserver,
  DOMParser: () => DOMParser,
  IntersectionObserver: () => IntersectionObserver,
  MutationObserver: () => MutationObserver,
  structuredClone: () => structuredClone,
  OffscreenCanvas: () => (typeof OffscreenCanvas === "undefined" ? undefined : OffscreenCanvas)
};

/** True when the host page is in dark mode. */
function isDark() {
  const root = document.documentElement;
  if (root.classList.contains("dark")) return true;               // MyST book-theme / Tailwind
  if (root.dataset.theme === "dark") return true;
  if (root.dataset.theme === "light") return false;
  return matchMedia?.("(prefers-color-scheme: dark)").matches ?? false;
}

/**
 * Keep `el` in sync with the host page's light/dark mode, and expose that as a
 * generator so notebooks can react to it.
 *
 * A single observer does both jobs, so the CSS custom properties are always
 * updated *before* cells re-run. That matters for anything hand-drawn on a
 * canvas, which has to resolve the inherited colour itself (canvas 2D does not
 * understand `currentColor`) — see `canvas.currentColor` in the kit.
 *
 * Notebooks opt in just by referencing `dark`:
 *
 *     drawQuiver(canvas, field, nx, ny, {color: dark ? "#eee" : "#222"})
 *
 * This replaces notebook-kit's built-in `dark`, which only watches
 * `prefers-color-scheme` and so would miss MyST's own theme toggle.
 */
function trackTheme(el, notify) {
  const apply = () => {
    const dark = isDark();
    el.dataset.theme = dark ? "dark" : "light";
    notify?.(dark);
  };
  apply();
  const mq = matchMedia?.("(prefers-color-scheme: dark)");
  mq?.addEventListener("change", apply);
  const mo = new MutationObserver(apply);
  mo.observe(document.documentElement, {attributes: true, attributeFilter: ["class", "data-theme"]});
  return () => {
    mq?.removeEventListener("change", apply);
    mo.disconnect();
  };
}

/**
 * Resolve a cell's import specifier.
 *
 * Relative specifiers resolve against the NOTEBOOK's url, not document.baseURI,
 * so a notebook that imports "./helper.js" keeps working when embedded in a MyST
 * page on a different domain.
 *
 * Note that the lab's own notebooks import the kit by ABSOLUTE url instead:
 * Observable Desktop sandboxes a notebook to its own folder, so a relative path
 * that climbs out of it ("../kit/index.js") cannot be resolved while editing.
 * This hook still matters for files that genuinely sit beside the notebook.
 */
function makeResolveImport(notebookUrl) {
  return (specifier) =>
    /^\.{0,2}\//.test(specifier)
      ? new URL(specifier, notebookUrl).href
      : resolveImportDefault(specifier);
}

/** Should this cell's output be displayed? `show` may hold cell ids or declared names. */
function isShown(show, cell, transpiled) {
  if (!show) return true;
  return show.some(
    (s) => String(s) === String(cell.id) || transpiled.outputs?.includes(String(s))
  );
}

/**
 * Turn a transpiled cell into a runnable definition.
 *
 * `transpile()` returns `body` as SOURCE TEXT — notebook-kit's Vite plugin emits
 * it into a generated module, where it becomes a function at load time. Doing the
 * transpile at runtime instead means we have to compile it ourselves.
 *
 * The body is a self-contained function expression whose every free reference is
 * passed in as a parameter, so a fresh `Function` scope is exactly right.
 *
 * Caveat: this needs `unsafe-eval`. MyST's themes do not set a restrictive CSP,
 * but a site that does would have to precompile notebooks at build time instead.
 */
function compile(transpiled, cellId) {
  try {
    return {...transpiled, body: new Function(`return (${transpiled.body});`)()};
  } catch (error) {
    throw new Error(`Cell ${cellId} failed to compile: ${error.message}`, {cause: error});
  }
}

function showError(el, error) {
  const pre = document.createElement("pre");
  pre.className = "em-notebook--error";
  pre.textContent = `Could not render notebook.\n\n${error?.stack ?? error}`;
  el.appendChild(pre);
}

export default {
  async render({model, el}) {
    const style = document.createElement("style");
    style.textContent = THEME_CSS;
    el.appendChild(style);

    const root = document.createElement("div");
    root.className = "em-notebook";
    el.appendChild(root);

    // One observer drives both the CSS variables and the `dark` builtin.
    let notifyDark;
    const untrackTheme = trackTheme(root, (value) => notifyDark?.(value));
    const darkGenerator = Generators.observe((notify) => {
      notifyDark = notify;
      notify(root.dataset.theme === "dark");
      return () => (notifyDark = undefined);
    });

    let runtime;

    try {
      const src = model.get("notebook");
      if (!src) throw new Error("No `notebook` given. Pass the URL of an Observable Notebook 2.0 .html file.");
      const notebookUrl = new URL(src, document.baseURI).href;

      const response = await fetch(notebookUrl);
      if (!response.ok) throw new Error(`Fetching ${notebookUrl} failed with HTTP ${response.status}.`);
      const notebook = deserialize(await response.text());

      const show = model.get("cells") ?? null;
      const hideCode = model.get("hideCode") ?? true;
      const resolveImport = makeResolveImport(notebookUrl);

      // `width` must observe our own container: notebook-kit's default watches
      // document.querySelector("main"), which is outside the shadow root.
      runtime = new NotebookRuntime({
        ...library,
        ...EXTRA_GLOBALS,
        width: () => Generators.width(root),
        dark: () => darkGenerator
      });

      for (const cell of notebook.cells) {
        const transpiled = transpile(cell, {resolveImport});

        // Resolve FileAttachment("x") against the notebook, not the host page.
        for (const name of transpiled.files ?? []) {
          registerFile(name, {path: name}, notebookUrl);
        }

        const cellRoot = document.createElement("div");
        cellRoot.className = "em-notebook--cell";

        // Every cell is DEFINED so dependencies resolve; only selected cells are
        // attached. The Observable runtime is lazy, so unobserved cells that
        // nothing depends on never run.
        if (isShown(show, cell, transpiled)) {
          root.appendChild(cellRoot);
          if (cell.pinned && !hideCode) {
            const pre = document.createElement("pre");
            pre.className = "em-notebook--code";
            pre.textContent = cell.value;
            root.appendChild(pre);
          }
        }

        runtime.define(
          {root: cellRoot, expanded: [], variables: []},
          {
            id: cell.id,
            ...compile(transpiled, cell.id),
            display: cell.mode === "js" || cell.mode === "ts" || cell.mode === "sql"
          }
        );
      }

      applyParams(runtime, model.get("params"));
      // MyST's model has no `off`, so this listener lives as long as the page.
      model.on("change:params", () => applyParams(runtime, model.get("params")));
    } catch (error) {
      console.error("observable-notebook:", error);
      showError(root, error);
    }

    return () => {
      untrackTheme();
      runtime?.runtime.dispose();
    };
  }
};

/** Override notebook cells with values supplied from the markdown call site. */
function applyParams(runtime, params) {
  if (!params || typeof params !== "object") return;
  for (const [name, value] of Object.entries(params)) {
    try {
      runtime.main.redefine(name, value);
    } catch (error) {
      console.warn(`observable-notebook: could not redefine "${name}":`, error);
    }
  }
}
