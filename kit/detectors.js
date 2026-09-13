// STEM detectors: which part of the diffraction pattern each one collects, and
// what signal that produces.
//
// Masks are built once per geometry and reused for every scan position, because
// the geometry changes when a slider moves but the scan runs thousands of times.
//
// All angles are in **mrad**; the grids come from `grid.angularSpatialFrequencies`,
// so they are corner-centered like the pattern straight out of `fft2`.

/**
 * A circular or annular detector.
 *
 * @param {{alpha: Float64Array}} angular from `angularSpatialFrequencies`
 * @param {number} innerMrad inner collection angle; 0 for a disc
 * @param {number} outerMrad outer collection angle
 * @returns {Float64Array} 1 inside, 0 outside
 */
export function annularMask({alpha}, innerMrad, outerMrad) {
  const inner = innerMrad * 1e-3;
  const outer = outerMrad * 1e-3;
  const mask = new Float64Array(alpha.length);
  for (let i = 0; i < alpha.length; i++) {
    mask[i] = alpha[i] >= inner && alpha[i] < outer ? 1 : 0;
  }
  return mask;
}

/**
 * A segmented detector: `rings` radial zones, each split into `sectors`
 * azimuthal segments.
 *
 * Each segment carries the mean scattering vector over its own area, which is
 * what lets a segmented detector estimate a centre of mass: the CoM is the
 * signal-weighted average of those segment centroids. That is a coarser
 * estimate than the first moment of a pixelated pattern — with 8 sectors the
 * azimuth is quantised to 45° — and the difference is worth showing rather
 * than hiding.
 *
 * The radii that matter are the ones straddling the **bright-field disc edge**,
 * around 0.5α to 1.05α: that is where a shifted disc moves intensity from one
 * sector to its opposite, and it is the geometry used in
 * `msa-em/stem-transfer-of-information` 08_segmented-icom. Segments far outside
 * the disc see dark-field scattering, whose azimuthal asymmetry is not the DPC
 * signal. Two rings and eight sectors is the Thermo Fisher Panther arrangement.
 *
 * @param {{alpha: Float64Array, phi: Float64Array, kx: Float64Array, ky: Float64Array}} angular
 * @param {object} options
 * @param {number} options.outer outer edge [mrad]
 * @param {number} [options.inner=0] inner edge [mrad]; a hole at the centre
 * @param {number} [options.rings=1] radial zones, evenly spaced in angle
 * @param {number} [options.sectors=8] azimuthal segments per zone
 * @param {number} [options.rotation] azimuthal offset of the segment boundaries,
 *   radians. Defaults to half a sector, which keeps the boundaries off the axes
 *   and the diagonals — where, on a square pixel grid, a great many pixel centres
 *   sit exactly on the dividing line and all get assigned to one side. Leaving
 *   them there makes sector areas differ by several percent, which biases the
 *   centre of mass. Real segmented detectors are rotated for the same reason.
 * @returns {{masks: Float64Array[], centroids: {kx: number, ky: number}[],
 *            sectors: number, rings: number, rotation: number, area: number[],
 *            edges: number[], bounds: {inner: number, outer: number}}}
 */
export function segmentedDetector(
  {alpha, phi, kx, ky},
  {outer, inner = 0, rings = 1, sectors = 8, rotation = Math.PI / 8}
) {
  const edges = [];
  for (let r = 0; r <= rings; r++) edges.push(inner + ((outer - inner) * r) / rings);

  const count = rings * sectors;
  const masks = Array.from({length: count}, () => new Float64Array(alpha.length));
  const sumKx = new Float64Array(count);
  const sumKy = new Float64Array(count);
  const area = new Float64Array(count);

  const twoPi = 2 * Math.PI;
  for (let i = 0; i < alpha.length; i++) {
    const mrad = alpha[i] * 1e3;
    let zone = -1;
    for (let r = 0; r < rings; r++) {
      if (mrad >= edges[r] && mrad < edges[r + 1]) { zone = r; break; }
    }
    if (zone < 0) continue;
    // atan2 gives (-π, π]; shift into [0, 2π) and rotate the boundaries off the
    // grid's symmetry directions.
    const angle = (((phi[i] - rotation) % twoPi) + twoPi) % twoPi;
    const sector = Math.min(sectors - 1, Math.floor((angle / twoPi) * sectors));
    const index = zone * sectors + sector;
    masks[index][i] = 1;
    sumKx[index] += kx[i];
    sumKy[index] += ky[i];
    area[index] += 1;
  }

  const centroids = [];
  for (let s = 0; s < count; s++) {
    const n = area[s] || 1;
    centroids.push({kx: sumKx[s] / n, ky: sumKy[s] / n});
  }

  return {
    masks, centroids, sectors, rings, rotation,
    area: Array.from(area), edges,
    bounds: {inner, outer}
  };
}

/**
 * Collapse a set of non-overlapping masks into one index map.
 *
 * A scanning widget integrates the same geometry over thousands of patterns, and
 * doing that mask by mask touches every pixel once per mask — 16 full passes for
 * a Panther-style detector. The segments do not overlap, so one pass suffices if
 * each pixel simply records which segment it belongs to.
 *
 * @param {Float64Array[]} masks non-overlapping, 1 inside and 0 outside
 * @param {number} length pixels per mask
 * @returns {Int16Array} the segment index per pixel, or -1 for none
 */
export function segmentIndex(masks, length) {
  const index = new Int16Array(length).fill(-1);
  for (let m = 0; m < masks.length; m++) {
    const mask = masks[m];
    for (let i = 0; i < length; i++) if (mask[i]) index[i] = m;
  }
  return index;
}

/**
 * {@link integrateSignals} in a single pass, using {@link segmentIndex}.
 *
 * @param {ArrayLike<number>} intensity
 * @param {Int16Array} index
 * @param {number} count number of segments
 * @param {Float64Array} [out] reused across scan positions
 */
export function integrateByIndex(intensity, index, count, out = new Float64Array(count)) {
  out.fill(0);
  for (let i = 0; i < index.length; i++) {
    const s = index[i];
    if (s >= 0) out[s] += intensity[i];
  }
  return out;
}

/**
 * Integrate an intensity pattern over a set of masks.
 *
 * @param {ArrayLike<number>} intensity
 * @param {Float64Array[]} masks
 * @returns {Float64Array} one sum per mask
 */
export function integrateSignals(intensity, masks) {
  const out = new Float64Array(masks.length);
  for (let m = 0; m < masks.length; m++) {
    const mask = masks[m];
    let sum = 0;
    for (let i = 0; i < mask.length; i++) sum += intensity[i] * mask[i];
    out[m] = sum;
  }
  return out;
}

/** Integrate over one mask — the common case, without allocating an array. */
export function integrateOne(intensity, mask) {
  let sum = 0;
  for (let i = 0; i < mask.length; i++) sum += intensity[i] * mask[i];
  return sum;
}

/**
 * The centre of mass a segmented detector would report: the signal-weighted
 * average of its segment centroids.
 *
 * @param {ArrayLike<number>} signals one per segment
 * @param {{kx: number, ky: number}[]} centroids
 * @param {number} [total] normalisation. Defaults to the sum over the segments,
 *   which makes the result a pure direction estimate. Pass the **whole
 *   pattern's** sum instead — as `08_segmented-icom` does — when the detector
 *   covers only part of the pattern: the signal then also carries how much
 *   intensity reached the detector at all, and stays dose-independent.
 * @returns {{kx: number, ky: number}} in the units of the angular grid (1/Å)
 */
export function centreOfMassFromSegments(signals, centroids, total = null) {
  let sum = 0, x = 0, y = 0;
  for (let s = 0; s < signals.length; s++) {
    sum += signals[s];
    x += signals[s] * centroids[s].kx;
    y += signals[s] * centroids[s].ky;
  }
  const norm = (total ?? sum) || 1;
  return {kx: x / norm, ky: y / norm};
}

/**
 * The first moment of the whole pattern — what a pixelated detector measures.
 *
 * Kept alongside the segmented version so the two can be compared; this is the
 * `py4DSTEM.process.utils.get_CoM` quantity, normalised by the total intensity
 * so it is dose-independent in expectation.
 */
export function centreOfMass(intensity, {kx, ky}, mask = null) {
  let total = 0, x = 0, y = 0;
  for (let i = 0; i < intensity.length; i++) {
    const w = mask ? intensity[i] * mask[i] : intensity[i];
    total += w;
    x += w * kx[i];
    y += w * ky[i];
  }
  const norm = total || 1;
  return {kx: x / norm, ky: y / norm};
}
