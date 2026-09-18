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

Seventeen widgets so far: `sem-ray-diagram`, `geometric-aberrations`,
`aperture-autocorrelation`, `probe-aberrations`, `stem-measurements`,
`stem-experiment` — a scanning 4D-STEM instrument whose multislice is checked
against abtem to ~1e-6 — `paraxial-rays`, which integrates a real lens field to
find its focal length, `non-paraxial-rays`, which traces the coupled 3D ray
equations to show Larmor rotation and spherical aberration, `electron-column`,
which drives an SEM or a S/TEM column from one table of components,
`reciprocity`, which builds a bright-field STEM column by reversing a TEM one,
`projection-sets`, which shows why alternating projections is not enough for
phase retrieval, `aperture-overlap`, which draws both slices of the aperture
overlap function every direct phase-retrieval method is built on, and
`direct-ptychography`, which runs SSB, OBF, parallax and iCOM through one
pipeline where only the kernel changes, `iterative-ptychography`, which
draws ePIE as the closed loop it is and lets you take one scan position at a
time, `surface-diffraction`, which scans a focused beam across an Si(111)
surface and maps where it has reconstructed from the fractional spots in the
pattern, `ptycho-tomography`, which measures a three-dimensional specimen with a
ptychographic scan at every tilt and reconstructs the series one projection at a
time or all at once, and `research-overview`, a row of three square thumbnails that run
themselves and hand over to the pointer on hover — the landing-page form of
three of the others.

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
- **A quantity with a size and a direction wants one control, not two.** Two
  linear sliders make the reader reassemble a vector, and a slider labelled
  "angle" never says how its number relates to the picture. `ui.vectorPad` is
  drawn as a miniature of the panel it acts on — first array axis down, second
  across — so the handle points where the aberration points, and its `angle` is
  already in the `phi` convention `chi` expects. `aperture-overlap` uses it for
  astigmatism and coma and keeps a plain slider for defocus, which has no axis.
  `_build/pw/drag.mjs` only exercises range inputs, so a pad needs its own
  check.
- **A cell that reads `form.value.x` depends on the whole form.** Observable
  tracks the reference, not the field, so any other control in the same
  `Inputs.form` re-runs it — which in `direct-ptychography` made choosing a
  different kernel re-simulate a second of data the kernel has no effect on.
  Where a cell is expensive, hand out one value per field (a `Generators.observe`
  per key that stays quiet when its own value has not changed) so each consumer
  depends only on what it reads.
- **Observable's Inputs carry an inline `max-width: 640px`.** It never bites in
  a multi-column control grid, where no cell is that wide, so it surfaces only
  when a control is given a row of its own — a five-option radio then wraps to
  two lines with 200px of empty space beside it. `ui.controls`' compact mode
  clears it.
- **A control that swaps its own label must not resize.** `ui.toggleButton`
  lays both labels out in one grid cell and hides the inactive one, so the
  button is always the size of the larger. Letting the text change instead
  resizes it — `▶` and `⏸` do not even share a line height, since they come
  from different fonts — and a play button that grows by two pixels nudges
  every panel below it down the page.
- **A four-dimensional dataset has to be laid out for whichever axis the
  algorithm walks.** The direct methods read one detector pixel's virtual image
  across the whole scan, so `[k][R]` puts what they need in one contiguous run;
  ePIE reads one scan position's whole pattern, so it wants `[R][k]`. The
  numbers are the same either way and the wrong one turns every read into a
  cache miss, so `forwardModel` takes a `layout` rather than picking one.
- **Match the stretch to what the array holds.** `image.histogramScaling` is the
  right default for a one-signed picture with content in most of it — a
  projected potential, say, where clipping the tails is what gives the structure
  the rest of the range. A signed field wants limits symmetric about zero
  instead, from a high percentile of the magnitude, so that zero stays the middle
  of the diverging map rather than landing wherever the data happened to be
  centred. And an array holding nothing yet wants neither: a percentile stretch
  of an all-zero array is a stretch of rounding error, so
  `iterative-ptychography` paints a flat grey and says in the readout that it is
  waiting. All three appear in that one widget.
- **A dose is per unit area of specimen, so the same number is a different
  camera.** Electrons per diffraction pattern are the dose times the area one
  scan position covers, and the three specimens in the ptychography widgets
  differ by sixty in that area — one dose figure would be four electrons per
  detector pixel for strontium titanate and thirty thousand for apoferritin.
  Give each specimen its own range, and report the per-pattern count next to the
  dose so the reader can see which of the two they are actually setting.
- **A shifted wave's spectrum carries a phase ramp, and domain colouring shows
  it.** An exit wave built on the whole grid with the probe rolled to the scan
  position transforms to something multiplied by `exp(-2 pi i k . R)` — correct,
  and what brings the correction home to the right place, but also a rainbow that
  turns over once per pixel of displacement. A ptychography code never sees it,
  because it crops the object around the probe. `epie.centreSpectrum` multiplies
  the ramp back out on the way to the canvas, which is the same picture that crop
  would give; the sign is the only thing that can go wrong, and the test pins it
  by rolling a probe and asserting the centred transform matches the unrolled
  one to 1e-16.
- **`frames({fps})` throttles the value it yields, not the loop body.** The
  runtime calls a generator cell once per animation frame either way, so a cell
  that does its work unconditionally runs at the display's refresh rate whatever
  `fps` says — it only matters to a cell that keys its animation off the yielded
  index. The consequence is that per-frame work counts do not set a speed:
  halving the count per frame just lets the frame rate rise to meet it (measured
  in `direct-ptychography`: 96 pixels a frame ran at 27 fps, 24 a frame ran at
  53). When a rate is what you want, pace against the clock — spend a stated
  number of seconds on the sweep and take `elapsed * rate` each frame, capped by
  the quality budget so a slow machine degrades instead of lurching. Leave direct
  manipulation at full speed: a click should land at once.
- **A step size of zero is an inspect mode for free.** `epieStep` at `beta: 0`
  runs the whole forward and backward calculation, fills in every intermediate
  stage, and leaves the reconstruction bit-for-bit where it was — so a widget can
  offer "show me what this position sees" without a second code path and without
  the reader changing the answer by asking. `iterative-ptychography` puts it on a
  toggle beside play, and relabels the arrow back into the specimen to say the
  loop is open.
- **A hover that does something must act on the state changing, not on the
  event.** `pointermove` keeps firing while the pointer sits still, so a handler
  that takes a step per event — or per frame — hammers whatever is under the
  cursor. In `iterative-ptychography` that meant the reconstruction overfitting
  one scan position until its gradient collapsed, which looks like a bug in the
  gradient. Record where the pointer is, convert that to a position in the frame
  loop, and act only when it differs from the last one.
- **A cell that attaches an event listener must not depend on anything that
  changes.** Observable re-runs the cell and the old listener stays on the
  element, so they accumulate one per parameter change, each writing into a state
  object nothing reads any more. Keep the handler's dependencies to the element
  it is bound to and a constant channel — `iterative-ptychography` has a
  dependency-free `pending` cell that both the hover and the reset button write
  into, and the frame loop is the only thing that acts on it.
- **A thrown exception inside an animation generator is silent.** A cell that
  yields frames swallows the error: no error block, no `pageerror`, nothing in
  the console — the panels simply stop being drawn while everything around them
  looks healthy. `_build/pw/run.mjs` catches it anyway, because a widget whose
  loop never ran has an empty `.em-badge` and canvases still at their CSS size.
- **Reshuffle the visiting order, then the per-position error curve is worth
  plotting.** How well one position fits says as much about which position it is
  — over the specimen or over vacuum — as about how far the reconstruction has
  got. With a fixed order that spread is replayed identically on every pass, so
  the curve is a sawtooth whose period is exactly one pass and which says nothing
  about convergence. Reshuffled, the same spread becomes the width of a band and
  the band's fall is the convergence, which is more informative than the per-pass
  mean that hides it. Plot it thin, drop the dots past a couple of hundred
  points, and halve the stored series when it gets long — otherwise a widget left
  playing through a lecture ends up rebuilding a path with a hundred thousand
  points on it. Halving also decouples the point count from the step number, so
  the axis domain has to come from the last point's own x.
- **An arrow into a panel that spans several rows needs to know which row it
  came from.** `ui.flowDiagram` joins panels centre to centre, which is right
  until a tall panel has more than one link: they then meet it at the same point
  and the two paths run along each other to get there. `anchor: "overlap"` on a
  link meets each panel at the middle of the span the two share, which is what
  closes the loop in `iterative-ptychography`.
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
| `ui.js` | `scene`, `panels`, `row`/`column`, `controls`, `collapsible`, `toggleButton`, `vectorPad` (magnitude and axis in one gesture), `flowDiagram` (panels joined by labelled arrows, canvases or arbitrary content), `UI` and `SERIES` (the lab's palettes) |
| `anim.js` | `frames`, `whenVisible`, `hoverTracker`, `sweep`/`pingPong`, `qualityBudget` — generator cells the Observable scheduler drives |
| `canvas.js` | `blit`, `drawQuiver`, `drawScalebar`, `colorbar`, `currentColor`, `pointerToIndices` |
| `projections.js` | `generalisedProjection`, `productProjection`, `namedParameters` (AP/DM/RRR/RAAR), `line`, `polarCurve`, `iteration`, `residual` |
| `ptycho.js` | `overlapFunction` and `overlapFunctionAtQ` (the two slices of Γ), `overlapRegions`, `overlapSums`, `parallaxShifts`, `tileSpectrum`, `directAccumulator`/`accumulatePixel`/`directImage`, `directCTF`, `directSSNR` |
| `ptycho-sim.js` | `measurePosition`, `forwardModel` (resumable; a detector subset and either data layout), `scanSpectra` |
| `epie.js` | `epieState` (complex or potential object, with positivity), `epieStep` (one position, object and probe), `epieReset`, `reconstructedPhase`, `scanOrder`, `probeDiameter`, `centreSpectrum`, `epieLine` (the same step in one dimension) |
| `tomo.js` | `projectVolume`, `backProject` (its exact transpose, which the tests pin), `tiltSeries` |
| `raytrace.js` | `transferMatrix`, `traceRays`, `traceParallel`, `traceFrom`, `COLUMNS`/`buildColumn`, `reverseColumn`, `findCrossovers`/`findPlanes`, `pairedSeparationAt` |
| `surface.js` | `reciprocalVectors`, `latticeOrigins`, `splatAtoms` (depth-weighted, so it is a surface measurement), `beamWindow`, `diffractionPattern`, `superstructureSpots`, `spotSum` (with the local background subtracted) |
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

Smaller arrays do not need any of this. the three `data/*-potential.f32`
specimens are 36 kB each of raw little-endian float32, read with one `fetch` and
wrapped in a `Float32Array`; as JSON each would have been half a megabyte, and
base64 inside JSON would have needed a decode step in the notebook. Either way, resolve the URL
through `resolveAsset` so a local checkout does not silently read production
data.

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

The ptychography modules check each other rather than a fixture, because the
closed forms are exact:

- `Γ(q, k=0)` is the HRTEM transfer function `-2i ψ(0) A(q) sin χ(q)` to 4e-19,
  which pins the sign of `chi`, the probe conjugation, the `k±q` index order and
  the fftfreq wraparound at once.
- the simulated data satisfies `G(q,k) = i Γ(q,k) φ̂(q)` to 1.5e-5 at an object
  amplitude of 1e-4, and the residual falls linearly with the object strength —
  so it is the weak-phase truncation and nothing else. This is the only test
  that ties the forward model to the inversion.
- ePIE run on that same forward model recovers a known weak phase object to
  0.01% of its own rms, its residual falls on every pass, and refining a probe
  started at the wrong defocus fits the data fifty times better than holding it
  there.
- solving for the potential rather than the complex object gives the same answer,
  and adding positivity brings the free additive constant down from 94% of the
  object's own rms to 20% — the measurement fixes `exp(i phi)` and so leaves
  `phi` floating, and the constraint is what pins it.
- a known probe comes out of 2560 steps bit for bit unchanged, a step size of
  zero fills every stage while leaving object and probe bit for bit unchanged,
  and a reset returns a run that has been refining both to exactly the state it
  started in.

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
