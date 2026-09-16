# em-widgets

Interactive electron-microscopy teaching widgets for the
[Curious Beams Lab](https://curiousbeams.github.io), built as
[Observable Notebooks 2.0](https://observablehq.com/notebook-kit/) and embedded in
MyST sites with a single `{anywidget}` wrapper.

No build step, no kernel, no server. Notebooks are plain HTML files in git, edited
live in Observable Desktop, and rendered client-side wherever they are embedded.

```
kit/                      shared primitives  ->  npm @curiousbeams/em-kit
notebooks/*.html          the widgets        ->  edit in Observable Desktop
observable-notebook.mjs   the MyST wrapper   ->  one file, used by every widget
test/                     standalone harness + numeric tests against the Python
data/*.zarr.zip           large arrays, chunked and read lazily
tools/                    fixture generation + npy -> zarr.zip conversion
demo/                     a throwaway MyST project for checking the real render
```

Eleven widgets so far: `sem-ray-diagram`, `geometric-aberrations`,
`aperture-autocorrelation`, `probe-aberrations`, `stem-measurements`,
`stem-experiment` — a scanning 4D-STEM instrument whose multislice is checked
against abtem to ~1e-6 — `paraxial-rays`, which integrates a real lens field to
find its focal length, `non-paraxial-rays`, which traces the coupled 3D ray
equations to show Larmor rotation and spherical aberration, `electron-column`,
which drives an SEM or a S/TEM column from one table of components,
`reciprocity`, which builds a bright-field STEM column by reversing a TEM one,
and `projection-sets`, which shows why alternating projections is not enough for
phase retrieval.

## Adding a widget to a page

```markdown
:::{anywidget} https://curiousbeams.github.io/em-widgets/observable-notebook.mjs
{
  // the body is JSON5 — comments, trailing commas and unquoted keys are fine
  notebook: "https://curiousbeams.github.io/em-widgets/notebooks/sem-ray-diagram.html",
  cells: ["rayDiagram"],       // optional: cell ids or names to show; omit for all
  params: {semiangle: 0.15},   // optional: override notebook cells from markdown
  hideCode: true               // optional: suppress echo of `pinned` cells
}
:::
```

`{anywidget}` is a built-in MyST directive, so this renders under a local
`myst start` as well as on the deployed site. Prefer it over Curvenote's
`{any:bundle}`, which only renders once deployed.

`params` is what makes one notebook reusable: the same widget can be embedded on
several pages with different defaults, via `main.redefine(name, value)`.

## Writing a widget

1. Open (or create) a notebook in `notebooks/` with Observable Desktop —
   or run `npm run preview` if you would rather use the browser.
2. Import what you need from the kit, **by absolute URL**:

   ```js
   import {chi, complexProbe, applyColormap} from "https://curiousbeams.github.io/em-widgets/kit/index.js"
   ```

   Not a relative path. Observable Desktop sandboxes a notebook to its own
   folder — its docs say you can open "a local file in the same folder" — so
   `../kit/index.js` cannot be reached and the notebook fails to run there.
   (Reproduce it with `npx notebooks preview --root notebooks`: `/kit/index.js`
   returns 404, because Vite will not serve above its root.)

   The consequence to be aware of: a notebook always runs against the
   **published** kit. Edit `kit/`, run `npm test`, push, then reload the notebook
   — local kit edits are not picked up by Desktop.
3. Give the cell you want to embed a name — `const rayDiagram = display(...)` —
   so a page can select it with `cells: ["rayDiagram"]`.

   **Prefer one display cell per piece** rather than one that assembles
   everything. A notebook that exposes `controlsView`, `sceneView`,
   `convergenceView` and `readoutView` separately can be embedded whole, or as
   just the diagram, or as the diagram with its main controls but without the
   parameter sliders — the page decides. Cells that are not displayed still run,
   so nothing about the computation depends on the choice. `projection-sets` is
   the example to copy. The harness takes the same option:
   `#nb=…&cells=controlsView,sceneView`.
4. Check it in the harness: `npm run serve`, then open
   <http://localhost:8080/test/>. The harness mounts the widget in a real open
   shadow root, exactly as MyST does, so styling behaves identically.

### Things worth knowing

- **Everything renders inside a shadow root.** Styles injected into
  `document.head` will not reach the widget. Put styles in `el`.
- The MyST model shim implements only `get`, `set` and `on`. `off`, `send` and
  `save_changes` throw.
- `npm:` imports inside cells resolve to jsDelivr automatically, so
  `import {mean} from "npm:d3-array"` just works.
- Observable's `md` cells do not render `$...$`. Use `${tex`\chi`}` instead.
- Imports and data URLs must be absolute (see above) — Desktop cannot see
  anything outside the notebook's own folder.
- **A cell that builds or re-parents DOM must not depend on a live input
  value.** Observable re-runs a cell whenever any input it reads changes, so if
  the cell that assembles the layout also reads a slider, every pixel of slider
  travel tears the widget down and rebuilds it. The slider being dragged is
  removed from the document and reinserted, which drops the pointer capture, and
  the symptom is that **sliders move one step per gesture and cannot be
  dragged**. Build the structure in one cell and mutate it (`style.top`,
  `replaceChildren`) from another that depends on the values. This has caused
  the same bug twice, in `non-paraxial-rays` (a legend that read the cone
  controls) and in `electron-column` (a slider parked at a back focal plane that
  moves). `_build/pw/drag.mjs <notebook>` drags every slider across its range
  and reports any that do not travel.
- Controls that only change the camera or re-sample an existing result are
  better kept out of the reactive graph entirely — a mutable state object plus a
  listener that repaints. See the tilt and defocus sliders in
  `non-paraxial-rays`.
- **Do not name a cell `view`.** It is an Observable builtin, and a cell
  declaring it does not shadow the builtin for other cells. They keep seeing the
  standard library's function, so `view.tilt` reads `undefined` and the failure
  surfaces much later as a `NaN` coordinate. `viewState` is a safe name.
- `Event` is not a recognised Observable global (`CustomEvent` is). The wrapper
  provides `Event`, `ResizeObserver`, `DOMParser`, `MutationObserver`,
  `IntersectionObserver` and `structuredClone` as builtins, so the site is
  strictly more permissive than Observable Desktop — anything that runs in
  Desktop runs when embedded.

## The kit

`kit/` is a dependency-free ES module library of the primitives that were
duplicated across the lab's Python notebooks.

| module | what's in it |
| --- | --- |
| `units.js` | `electronWavelength`, `interactionParameter` |
| `grid.js` | `fftfreq`, `spatialFrequencies`, `polarCoordinates`, `angularSpatialFrequencies`, `checkerboard` |
| `fft.js` | `fft2`/`ifft2` (radix-2 + Bluestein, any size), `fftshift`, `fourierShift`, complex helpers |
| `optics.js` | `chi` to 6th order, `rayDisplacement`, `softAperture`, `complexProbe`, `fresnelPropagator`, `multislice`, `antialiasAperture`, `scherzerDefocus` |
| `specimen.js` | `fccCluster`, `amorphousSupport`, `projectedPotential` — the Lobato/Fourier superposition abtem uses |
| `paraxial.js` | `schiskeField`, `glaserField`, `integrateRay`, `principalRays`, `imagingProperties`, `thinLensPower` — **SI units**, unlike the rest of the kit |
| `nonparaxial.js` | `laplaceExpansion`, `equationsOfMotion`, `traceRay3D`, `coneOfRays`, `larmorAngle`, `axialCrossing`, `spotDiagram` — also **SI units** |
| `detectors.js` | `annularMask`, `segmentedDetector`, `centreOfMassFromSegments`, `segmentIndex` |
| `image.js` | `histogramScaling`, `radialAverage`, `integrateGradient`, `gradient2d`, `warpNearest`, `poisson`, `cropCenter` |
| `color.js` | `complexToRGB` (inverse CIECAM02), `phaseWheel`, `applyColormap` — magma, gray, twilight, RdBu, PuOr, PiYG, eclipse |
| `scene3d.js` | `makeView`, `project`/`unproject`, `drawSpheres`, `drawPlane`, `drawProbeCone` — canvas 2D, no WebGL |
| `ui.js` | `scene`, `panels`, `row`/`column`, `controls`, `collapsible`, `toggleButton`, `UI` and `SERIES` (the lab's palettes) |
| `anim.js` | `frames`, `whenVisible`, `qualityBudget` — generator cells the Observable scheduler drives |
| `canvas.js` | `blit`, `drawQuiver`, `drawScalebar`, `colorbar`, `currentColor`, `pointerToIndices` |
| `projections.js` | `generalisedProjection`, `productProjection`, `namedParameters` (AP/DM/RRR/RAAR), `line`, `polarCurve`, `iteration`, `residual` |
| `raytrace.js` | `transferMatrix`, `traceRays`, `traceParallel`, `traceFrom`, `COLUMNS`/`buildColumn`, `reverseColumn`, `findCrossovers`/`findPlanes`, `pairedSeparationAt` |
| `data.js` | `openZarrZip`, `readAll`, `readSlice` — **not** re-exported from `index.js`, so widgets that use no data never load zarrita |

### Conventions

- Lengths in Angstroms, energies in eV, angles in radians (unless a name says mrad).
- 2D arrays are flat and row-major, `[ix * ny + iy]` — numpy's `indexing="ij"`.
  Whether `ix` is vertical (`imshow`) or horizontal (a scan whose first index runs
  along the field) is the caller's choice: `toImageData` assumes the first,
  `transposeToImageData` the second. Picking the wrong one is a silent 90° turn.
- Complex arrays are `{re, im}` pairs of `Float64Array`s.
- Reciprocal-space grids are **corner-centered** (fftfreq order). Only `fftshift`
  when you are about to draw.

### Large datasets

Arrays bigger than a megabyte or so are stored as zarr inside a ZipStore: one
file to host and version, chunked inside, read lazily over HTTP range requests.

```sh
python tools/npy-to-zarr-zip.py path/to/array.npy \
    --shape 7 244 242 --dtype float32 \
    --out data/fcc-slab-potential.zarr.zip
```

```js
import {openZarrZip, readAll, readSlice} from "../kit/data.js"
const potential = await openZarrZip("https://curiousbeams.github.io/em-widgets/data/fcc-slab-potential.zarr.zip");
const oneSlice = await readSlice(potential, 0);   // fetches one chunk, not the file
```

Measured on the 72×72×64×64 uint8 streaming dataset: **20.25 MB → 0.83 MB**
with zstd (it is sparse), a single chunk costs ~12 KB, and reading the whole
array costs 899 KB across 148 range requests. GitHub Pages, jsDelivr and MyST's
own dev server all support `HEAD` + `Range` with CORS, so this works wherever
you host it.

Two rules, both enforced by the tool:

- the zip must be **STORED**, not deflated, or entries cannot be range-read;
- pick a codec the browser can decode — `zstd` and `gzip` both work through
  zarrita; zstd compresses better on sparse data, gzip on dense.

### Tests

`npm test` checks the numerics against reference values generated from the lab's
own Python (numpy, py4DSTEM, colorspacious):

- `fft2` agrees with `numpy.fft.fft2` to ~1e-14, including non-power-of-two sizes.
- `jchToSrgb` agrees with `colorspacious.cspace_convert(JCh, "JCh", "sRGB1")` to ~1e-14.
- `chi`, `histogramScaling`, `radialAverage`, `electronWavelength` and the
  apertures agree with `aberration_utils.py` / `ctf/` to float32 precision (the
  precision of the Python reference).
- `complexProbe` and its inverse transform agree with **abtem's `CTF`** to 1.5e-7
  — the end-to-end convention test. It pins the fftfreq ordering, the sign of
  chi (`exp(-i chi)`, and the wrong sign is checked to *fail*), the `ifft2`
  normalisation and the `fftshift`. Note abtem's `defocus` is `-C10`.

`electronWavelength` differs from abtem by 1.3e-7 relative, because the kit uses
the lab's Python constants and abtem uses CODATA. The kit deliberately tracks the
notebooks it was ported from.

Fixtures live in `test/fixtures/`. Regenerate them with `tools/generate-fixtures.py`,
run from a repo whose venv has the lab's Python stack — the script's docstring has
the exact invocation.

## Publishing

- **Widgets and notebooks** are served from this repo's GitHub Pages, at
  `https://curiousbeams.github.io/em-widgets/`. Pushing to `main` publishes them.
- **The kit** is additionally published to npm as `@curiousbeams/em-kit`, so
  notebooks elsewhere can `import ... from "npm:@curiousbeams/em-kit"`.

## Checking the real MyST render

```sh
cd demo && myst start
```

`demo/myst.yml` uses `static_files` to serve `notebooks/` and `kit/` from the MyST
site itself, which is also a viable hosting arrangement if you ever prefer it to
this repo's Pages origin.
