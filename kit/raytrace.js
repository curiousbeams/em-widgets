// Paraxial ray tracing by transfer matrices.
//
// Lifted verbatim from the pinned helper cells of
// `observable-notebooks/sem-ray-diagram.html` so that notebook (and the AP3402
// w4.* ray-diagram pages, which currently ship static Mathematica output) can
// share one implementation.
//
// A ray is `[r, theta]`: off-axis distance and slope. Optical components are
// plain objects with an `axial_position` and a `type` of "lens", "aperture" or
// "sample"; everything is in consistent length units (the SEM notebook uses metres).

/** 2x2 paraxial transfer matrices. */
export const transferMatrix = {
  /** Free-space propagation over a distance d. */
  propagation: (d) => [
    [1, d],
    [0, 1]
  ],
  /** A thin lens of focal length f. */
  lens: (f) => [
    [1, 0],
    [-1 / f, 1]
  ]
};

/** Apply a 2x2 matrix to a ray. */
export const applyTransferMatrix = (M, ray) => [
  M[0][0] * ray[0] + M[0][1] * ray[1],
  M[1][0] * ray[0] + M[1][1] * ray[1]
];

/**
 * Sort components by axial position and interleave free-space propagation steps,
 * so the result can be walked in a single pass.
 *
 * Components above `zStart` are dropped rather than applied where they sit,
 * which is what lets a bundle be launched from part-way down the column — a ray
 * leaving the specimen has not passed through the condensers.
 */
export function postprocessOpticalSequence(components, zStart, zEnd) {
  const sorted = components
    .filter((c) => c.axial_position >= zStart)
    .sort((a, b) => a.axial_position - b.axial_position);
  const sequence = [];
  let zPrev = zStart;
  for (const component of sorted) {
    if (component.axial_position > zPrev) {
      sequence.push({type: "propagation", dz: component.axial_position - zPrev});
    }
    sequence.push(component);
    zPrev = component.axial_position;
  }
  if (zPrev < zEnd) sequence.push({type: "propagation", dz: zEnd - zPrev});
  return sequence;
}

/**
 * Trace a single ray through a prepared sequence.
 *
 * Returns the trajectory as `{z, r}` samples. A ray blocked by an aperture stops
 * there, which is what makes the aperture sliders visually meaningful.
 */
export function traceRaySequence(ray0, sequence, zStart = 0) {
  let ray = [...ray0];
  let z = zStart;
  // Each sample carries the slope as well as the height, so a bundle can be
  // picked up at a plane and continued — which is what a specimen does to it.
  const trajectory = [{z, r: ray[0], theta: ray[1]}];

  for (const element of sequence) {
    switch (element.type) {
      case "propagation":
        z += element.dz;
        ray = applyTransferMatrix(transferMatrix.propagation(element.dz), ray);
        trajectory.push({z, r: ray[0], theta: ray[1]});
        break;
      case "lens":
        ray = applyTransferMatrix(transferMatrix.lens(element.focal_length), ray);
        trajectory.push({z, r: ray[0], theta: ray[1]});
        break;
      case "deflector":
        // A pure angular kick, with no change of position. Two of them in
        // series shift the beam without tilting it, which is how a real column
        // scans a probe across the specimen.
        ray = [ray[0], ray[1] + element.deflection];
        trajectory.push({z, r: ray[0], theta: ray[1]});
        break;
      case "aperture":
        if (Math.abs(ray[0]) > element.aperture_radius) return trajectory; // blocked
        break;
    }
  }
  return trajectory;
}

/**
 * Continue a bundle from a plane, given each ray's state there.
 *
 * What a specimen does to a beam: every ray arrives with its own height and
 * slope, and leaves from the same height with its slope changed by the
 * scattering angle. Re-launching the whole illumination this way gives the
 * diffracted beam, in any illumination condition — parallel, where it becomes a
 * second spot in the back focal plane, or converged, where it becomes a second
 * overlapping disc.
 *
 * @param {object[]} components
 * @param {number} zStart the plane the states were taken at
 * @param {number} zEnd
 * @param {Array<[number, number]>} states one `[r, theta]` per ray
 * @returns {Array<Array<{z: number, r: number, theta: number}>>}
 */
export function traceFrom(components, zStart, zEnd, states) {
  const sequence = postprocessOpticalSequence(components, zStart, zEnd);
  return states.map((state) => traceRaySequence(state, sequence, zStart));
}

/**
 * Each ray's `[r, theta]` where it crosses a plane.
 *
 * Rays that stopped at an aperture before reaching it are dropped, so the
 * bundle handed to {@link traceFrom} is the one that actually got there.
 */
export function statesAt(trajectories, z, tolerance = 1e-9) {
  const out = [];
  for (const trajectory of trajectories) {
    const point = trajectory.find((p) => Math.abs(p.z - z) < tolerance);
    if (point) out.push([point.r, point.theta]);
  }
  return out;
}

/**
 * Trace a fan of rays leaving the origin over +/- `semiangle`.
 *
 * @param {object[]} components
 * @param {number} zStart
 * @param {number} zEnd
 * @param {number} semiangle half-angle of the fan [rad]
 * @param {number} [numRays=128]
 * @returns {Array<Array<{z: number, r: number}>>}
 */
export function traceRays(components, zStart, zEnd, semiangle, numRays = 128, options = {}) {
  const {height = 0, tilt = 0} = options;
  const sequence = postprocessOpticalSequence(components, zStart, zEnd);
  const angles = Array.from(
    {length: numRays + 1},
    (_, i) => -semiangle + (2 * semiangle * i) / numRays
  );
  return angles.map((alpha) =>
    traceRaySequence([height, Math.tan(alpha + tilt)], sequence, zStart)
  );
}

/**
 * Trace a collimated bundle: parallel rays spread over `halfWidth`, all at the
 * same slope.
 *
 * The counterpart of {@link traceRays}. A point source fans out in angle at one
 * position; parallel illumination spreads in position at one angle. Those are
 * the two halves of a TEM ray diagram, and which of them the projector lands on
 * the detector is the whole difference between imaging and diffraction mode.
 *
 * @param {object[]} components
 * @param {number} zStart
 * @param {number} zEnd
 * @param {number} halfWidth half-width of the bundle at `zStart`
 * @param {number} [numRays=128]
 * @param {object} [options]
 * @param {number} [options.tilt=0] common slope of every ray [rad]
 * @returns {Array<Array<{z: number, r: number}>>}
 */
export function traceParallel(components, zStart, zEnd, halfWidth, numRays = 128, options = {}) {
  const {tilt = 0} = options;
  const sequence = postprocessOpticalSequence(components, zStart, zEnd);
  const slope = Math.tan(tilt);
  return Array.from({length: numRays + 1}, (_, i) => {
    const r = -halfWidth + (2 * halfWidth * i) / numRays;
    return traceRaySequence([r, slope], sequence, zStart);
  });
}

// ---------------------------------------------------------------------------
// Reading the bundle
// ---------------------------------------------------------------------------

/**
 * The bundle's rms radius at one axial position.
 *
 * Trajectories are piecewise straight, so a segment is interpolated exactly
 * rather than approximated. Rays that stopped at an aperture are skipped past
 * their end, which is the point: a bundle that has lost its outer rays really
 * is narrower downstream.
 *
 * @param {Array<Array<{z: number, r: number}>>} trajectories
 * @param {number} z
 * @returns {number} rms radius, or 0 if no ray reaches `z`
 */
export function bundleRadiusAt(trajectories, z) {
  const heights = bundleAt(trajectories, z);
  if (heights.length === 0) return 0;
  let sum = 0;
  for (const r of heights) sum += r * r;
  return Math.sqrt(sum / heights.length);
}

/**
 * Every ray's height at one axial position.
 *
 * @param {Array<Array<{z: number, r: number}>>} trajectories
 * @param {number} z
 * @returns {number[]} one entry per ray that reaches `z`
 */
export function bundleAt(trajectories, z) {
  const out = [];
  for (const trajectory of trajectories) {
    const r = rayHeightAt(trajectory, z);
    if (r !== null) out.push(r);
  }
  return out;
}

/**
 * One ray's height at an axial position, or `null` if it never gets there.
 *
 * Trajectories are piecewise straight, so the segment containing `z` is
 * interpolated exactly. A ray stopped by an aperture simply ends, and returning
 * `null` rather than its last height is what lets two bundles be compared ray by
 * ray without counting rays that only one of them still has.
 */
export function rayHeightAt(trajectory, z) {
  for (let i = 0; i < trajectory.length - 1; i++) {
    const a = trajectory[i];
    const b = trajectory[i + 1];
    if (z < a.z || z > b.z) continue;
    const span = b.z - a.z;
    return span === 0 ? b.r : a.r + ((z - a.z) / span) * (b.r - a.r);
  }
  return null;
}

/**
 * How far apart two bundles are, ray by ray, at one axial position.
 *
 * The bundles must be index-aligned — the same rays, sent two different ways —
 * and only indices that survive in both are counted.
 *
 * Comparing the two bundles' *mean* heights instead looks equivalent and is
 * not, as soon as an aperture clips one of them. The aperture takes the rays
 * furthest off axis, and a tilted beam is not clipped symmetrically, so its
 * surviving mean shifts. The image plane then never quite registers, because
 * the two means never quite coincide. Per-ray, the shift cannot happen: at an
 * image plane every ray lands where its partner lands, however few are left.
 *
 * @returns {{separation: number, spread: number, count: number}} mean absolute
 *   per-ray difference, the spread of the first bundle over the same rays, and
 *   how many rays were common to both
 */
export function pairedSeparationAt(a, b, z) {
  let gap = 0;
  let sum = 0;
  let sumSquares = 0;
  let count = 0;
  const n = Math.min(a.length, b.length);
  for (let i = 0; i < n; i++) {
    const ra = rayHeightAt(a[i], z);
    if (ra === null) continue;
    const rb = rayHeightAt(b[i], z);
    if (rb === null) continue;
    gap += Math.abs(ra - rb);
    sum += ra;
    sumSquares += ra * ra;
    count++;
  }
  if (count === 0) return {separation: Infinity, spread: 0, count: 0};
  const mean = sum / count;
  return {
    separation: gap / count,
    spread: Math.sqrt(Math.max(0, sumSquares / count - mean * mean)),
    count
  };
}

/**
 * How tightly the bundle is gathered at one axial position, measured about its
 * own centre rather than about the optic axis.
 *
 * This is what "in focus" means for a bundle that left the specimen off axis.
 * Its image is off axis too — magnified and inverted — so the distance from the
 * axis there is large even though every ray has converged to the same point.
 * Measuring about the axis instead finds no image at all, which is exactly what
 * the first version of this did.
 *
 * For a symmetric fan about the axis the mean is zero and this equals
 * {@link bundleRadiusAt}.
 */
export function bundleSpreadAt(trajectories, z) {
  const heights = bundleAt(trajectories, z);
  if (heights.length === 0) return 0;
  let mean = 0;
  for (const r of heights) mean += r;
  mean /= heights.length;
  let sum = 0;
  for (const r of heights) sum += (r - mean) ** 2;
  return Math.sqrt(sum / heights.length);
}

/**
 * Where the bundle comes to a focus.
 *
 * Every local minimum of {@link bundleSpreadAt}, which for an illumination fan
 * is the source crossovers and for a bundle leaving the specimen is the image
 * planes. Marking them is what turns four presets from "different slider
 * values" into a picture of which plane the projector is throwing onto the
 * detector.
 *
 * A scan rather than anything cleverer, because the radius is only piecewise
 * smooth — it kinks at every lens — so a derivative-based search would stop at
 * the kinks.
 *
 * @param {Array<Array<{z: number, r: number}>>} trajectories
 * @param {object} options
 * @param {number} options.zStart
 * @param {number} options.zEnd
 * @param {number} [options.samples=1200]
 * @param {number} [options.tolerance=0.02] ignore minima wider than this
 *   fraction of the widest the bundle gets
 * @returns {Array<{z: number, radius: number}>}
 */
export function findCrossovers(trajectories, options) {
  return findPlanes((z) => bundleSpreadAt(trajectories, z), options);
}

/**
 * Every axial position where a measure of the beam reaches a local minimum.
 *
 * The general form of {@link findCrossovers}. Pass the bundle's spread to find
 * where it focuses; pass the distance between two bundles' centres to find
 * where they recombine, which is the image plane in a diagram drawn as a direct
 * beam plus a diffracted one.
 *
 * `edges` says which ends of the scan may themselves be reported. The far end
 * usually should be: a column is normally set up so that some plane lands on
 * the detector, and an interior-only scan silently failed to mark precisely the
 * plane the whole diagram is about. The near end usually should not, because a
 * scan started just past the specimen begins where the measure is trivially
 * small and would report that as a plane.
 *
 * @param {(z: number) => number} measure
 * @param {object} options
 * @param {number} options.zStart
 * @param {number} options.zEnd
 * @param {number} [options.samples=1200]
 * @param {number} [options.tolerance=0.02] ignore minima above this fraction of
 *   the measure's largest value
 * @param {"both"|"start"|"end"|"none"} [options.edges="end"]
 * @returns {Array<{z: number, radius: number}>}
 */
export function findPlanes(measure, {
  zStart, zEnd, samples = 1200, tolerance = 0.02, edges = "end"
}) {
  const allowStart = edges === "both" || edges === "start";
  const allowEnd = edges === "both" || edges === "end";
  const values = new Float64Array(samples + 1);
  let largest = 0;
  for (let i = 0; i <= samples; i++) {
    values[i] = measure(zStart + ((zEnd - zStart) * i) / samples);
    // A measure may report a plane as meaningless rather than as large — a
    // bundle that an aperture has stopped entirely has no position at all.
    // Such samples must not set the scale, or the cutoff becomes infinite and
    // every sample passes it.
    if (Number.isFinite(values[i]) && values[i] > largest) largest = values[i];
  }
  const cutoff = largest * tolerance;
  const step = (zEnd - zStart) / samples;

  const out = [];
  for (let i = 0; i <= samples; i++) {
    if (!Number.isFinite(values[i]) || values[i] > cutoff) continue;
    if (i === 0 && !allowStart) continue;
    if (i === samples && !allowEnd) continue;
    const before = i === 0 ? Infinity : values[i - 1];
    const after = i === samples ? Infinity : values[i + 1];
    if (values[i] > before || values[i] > after) continue;

    // Refine with a parabola through the three samples, so the plane is located
    // to better than the scan step. Not available at an endpoint, which has
    // only one neighbour.
    let shift = 0;
    if (i > 0 && i < samples) {
      const denom = before - 2 * values[i] + after;
      if (Math.abs(denom) > 1e-30) {
        shift = Math.max(-1, Math.min(1, (0.5 * (before - after)) / denom));
      }
    }
    const z = zStart + i * step + shift * step;
    // One entry per plane: a flat minimum spans several samples.
    if (out.length && z - out.at(-1).z < 4 * step) continue;
    out.push({z, radius: values[i]});
  }
  return out;
}

/** The bundle's mean height at one axial position. */
export function bundleCentreAt(trajectories, z) {
  const heights = bundleAt(trajectories, z);
  if (heights.length === 0) return 0;
  let sum = 0;
  for (const r of heights) sum += r;
  return sum / heights.length;
}

// ---------------------------------------------------------------------------
// Columns, as data
// ---------------------------------------------------------------------------

/**
 * Build a component list from a column description and a set of slider values.
 *
 * An element's `z` may be a function of the values rather than a number, which
 * is how the objective aperture follows the back focal plane as the objective
 * is refocused.
 *
 * @param {object} column an entry of {@link COLUMNS}
 * @param {Record<string, number>} values keyed by element `key`
 */
export function buildColumn(column, values) {
  return column.elements.map((element) => {
    const axial_position = typeof element.z === "function" ? element.z(values) : element.z;
    const component = {type: element.kind, name: element.name, axial_position};
    if (element.kind === "lens") component.focal_length = values[element.key];
    if (element.kind === "aperture") component.aperture_radius = values[element.key];
    if (element.kind === "deflector") component.deflection = values[element.key];
    return component;
  });
}

/** Default slider values for a column in a given mode. */
export function columnDefaults(column, mode) {
  const overrides = column.modes?.[mode] ?? {};
  const values = {};
  for (const element of column.elements) {
    if (element.key === undefined) continue;
    values[element.key] = overrides[element.key] ?? element.value;
  }
  return values;
}

/**
 * The lab's two teaching columns.
 *
 * Written as data rather than as a builder function per column, because the
 * only thing that separates the three S/TEM modes is a handful of default focal
 * lengths — the components are identical. The Observable original encoded that
 * as ternaries on the preset name inside each slider; a table says the same
 * thing once, and lets one widget generate the controls for either column by
 * mapping over whichever list is active.
 *
 * Positions are in metres and match the figures the two notebooks were built
 * from. The S/TEM column is the one from `@gvarnavi/stem-ray-diagrams`, itself
 * after the EPFL LSME TEM ray diagrams (CC-BY 4.0, Duncan T.L. Alexander).
 */
export const COLUMNS = {
  SEM: {
    label: "SEM",
    // Ends at the specimen. There is nothing below it in an SEM, and drawing
    // the fan continuing past it reads as the beam going through.
    zEnd: 1.5,
    semiangle: 0.15,
    modes: {},
    elements: [
      {kind: "lens", key: "c1", name: "1st condenser lens", label: "C1 focal length",
       z: 0.2, range: [0.025, 0.5], step: 0.005, value: 0.05},
      {kind: "aperture", key: "spray", name: "spray aperture", label: "spray aperture",
       z: 0.5, range: [0.025, 0.2], step: 0.005, value: 0.05},
      {kind: "lens", key: "c2", name: "2nd condenser lens", label: "C2 focal length",
       z: 0.7, range: [0.025, 0.5], step: 0.005, value: 0.15},
      {kind: "aperture", key: "objAperture", name: "objective aperture", label: "obj aperture",
       z: 1.2, range: [0.025, 0.2], step: 0.005, value: 0.05},
      // 0.13 rather than the 0.35 the standalone SEM notebook shipped with,
      // which leaves the probe 38 mm wide at the specimen. Starting focused
      // reads as a working instrument, and defocusing it is one slider away.
      {kind: "lens", key: "objective", name: "objective lens", label: "obj focal length",
       z: 1.3, range: [0.025, 0.5], step: 0.005, value: 0.13},
      {kind: "sample", name: "sample", z: 1.5}
    ]
  },

  "S/TEM": {
    label: "S/TEM",
    zEnd: 4.3,
    semiangle: 0.15,
    // The three modes differ only in these defaults. Everything else — the
    // components, their positions, the slider ranges — is shared.
    modes: {
      "TEM imaging": {c3: 0.07, objAperture: 0.01, intermediate: 0.16},
      "TEM diffraction": {c3: 0.07, objAperture: 0.01, intermediate: 0.26},
      "BF STEM": {c3: 1.0, objAperture: 0.15, intermediate: 0.295}
    },
    elements: [
      {kind: "lens", key: "c1", name: "1st condenser lens", label: "C1 focal length",
       z: 0.35, range: [0.025, 0.5], step: 0.005, value: 0.11},
      {kind: "lens", key: "c2", name: "2nd condenser lens", label: "C2 focal length",
       z: 0.85, range: [0.025, 0.5], step: 0.005, value: 0.21},
      {kind: "aperture", key: "c2Aperture", name: "C2 aperture", label: "C2 aperture",
       z: 1.0, range: [0.025, 0.2], step: 0.005, value: 0.05},
      {kind: "lens", key: "c3", name: "3rd condenser lens", label: "C3 focal length",
       z: 1.55, range: [0.025, 1.0], step: 0.005, value: 0.07},
      {kind: "lens", key: "objPre", name: "objective lens pre-field", label: "pre objective",
       z: 1.83, range: [0.025, 0.5], step: 0.005, value: 0.15},
      {kind: "sample", name: "sample", z: (1.83 + 2.28) / 2},
      {kind: "lens", key: "objPost", name: "objective lens post-field", label: "post objective",
       z: 2.28, range: [0.025, 0.5], step: 0.005, value: 0.15},
      // The back focal plane, which moves when the objective is refocused.
      {kind: "aperture", key: "objAperture", name: "objective aperture", label: "obj aperture",
       z: (v) => 2.28 + v.objPost, range: [0.005, 0.2], step: 0.005, value: 0.01},
      {kind: "lens", key: "intermediate", name: "intermediate lens", label: "int. focal length",
       z: 2.96, range: [0.025, 0.5], step: 0.005, value: 0.16},
      {kind: "lens", key: "projector", name: "projector lens", label: "proj. focal length",
       z: 3.7, range: [0.025, 0.5], step: 0.005, value: 0.165},
      {kind: "detector", name: "detector", z: 4.3}
    ]
  }
};

/**
 * The same column with its ends swapped: what was nearest the source is now
 * nearest the detector.
 *
 * This is the principle of reciprocity as an operation. A ray that runs from a
 * source at one end to a detector at the other follows the same path when the
 * two are exchanged, so a TEM — parallel illumination in, objective and
 * aperture after the specimen — and a bright-field STEM — aperture and
 * objective before the specimen, point detector after — are the same
 * instrument read in opposite directions.
 *
 * Every component keeps its strength, deflectors included. That looks wrong and
 * is not: the reversed ray reaches a deflector with slope `-(theta + delta)`
 * and has to leave it with `-theta`, so the kick it needs is `+delta`, the same
 * one. The sign flip the reversal needs is in the slope, applied once when the
 * ray is launched, not in the components.
 *
 * `test/raytrace.test.js` asserts the theorem the operation is named for.
 *
 * @param {object[]} components
 * @param {object} span `{zStart, zEnd}` — the ends being swapped
 */
export function reverseColumn(components, {zStart, zEnd}) {
  return components.map((component) => ({
    ...component,
    axial_position: zStart + zEnd - component.axial_position
  }));
}

/** Which column each selectable mode uses. */
export const COLUMN_MODES = {
  SEM: "SEM",
  "TEM imaging": "S/TEM",
  "TEM diffraction": "S/TEM",
  "BF STEM": "S/TEM"
};

/**
 * The two-condenser SEM column used by the sem-ray-diagram notebook.
 *
 * @param {[number, number, number]} focalLengths C1, C2 and objective focal lengths
 * @param {[number, number]} apertureRadii spray and objective aperture radii
 */
export const semElectroOpticalComponents = ([f1, f2, f3], [ap1, ap2]) => [
  {type: "lens", name: "1st condenser lens", axial_position: 0.2, focal_length: f1},
  {type: "aperture", name: "spray aperture", axial_position: 0.5, aperture_radius: ap1},
  {type: "lens", name: "2nd condenser lens", axial_position: 0.7, focal_length: f2},
  {type: "aperture", name: "objective aperture", axial_position: 1.2, aperture_radius: ap2},
  {type: "lens", name: "objective lens", axial_position: 1.3, focal_length: f3},
  {type: "sample", name: "sample", axial_position: 1.5}
];
