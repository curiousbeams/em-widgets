// Surface diffraction: a reconstructed surface, a beam focused on part of it,
// and the pattern that comes back.
//
// A clean surface reconstructs — the top layer rearranges into a larger unit
// cell than the bulk below it. Si(111) does it seven-fold, into the
// dimer-adatom-stacking-fault cell, and a larger cell in real space means
// extra spots in reciprocal space: six of them between every pair of substrate
// spots, at sevenths. Those fractional spots are the signature. They say that
// the surface has reconstructed, and nothing else in the pattern does.
//
// A beam wide enough to cover the whole surface sees the fractional spots of
// every reconstructed patch at once and cannot say where any of them is. Focus
// it down and scan it, and the pattern changes as the beam crosses from
// reconstructed to bulk-terminated ground — so the map of how much intensity
// sits in the fractional spots is an image of where the surface has
// reconstructed. That is the idea behind scanning surface diffraction, and this
// module is the arithmetic behind it.
//
// Conventions. Lengths are Angstroms and the field of view is square, `extent`
// across, sampled on an `n x n` grid; reciprocal-space quantities are in 1/A,
// so the DFT's own sample spacing is `1 / extent`. Patterns come back
// corner-centered like every other spectrum in this kit — `fftshift` only when
// you are about to draw. The surface is a finite patch rather than a periodic
// cell: nothing wraps, and the beam has to stay clear of the edges.

import {complex, fft2} from "./fft.js";

/**
 * The reciprocal vectors of a two-dimensional lattice: `a_i · b_j = δ_ij`.
 *
 * Crystallographer's convention, without the 2π — the same one the rest of this
 * kit uses for spatial frequencies, so `b` is directly a number of 1/Å.
 *
 * @param {[number, number]} a1
 * @param {[number, number]} a2
 * @returns {{b1: [number, number], b2: [number, number]}}
 */
export function reciprocalVectors(a1, a2) {
  const determinant = a1[0] * a2[1] - a1[1] * a2[0];
  if (Math.abs(determinant) < 1e-12) throw new Error("surface: the lattice vectors are parallel.");
  return {
    b1: [a2[1] / determinant, -a2[0] / determinant],
    b2: [-a1[1] / determinant, a1[0] / determinant]
  };
}

/**
 * Every lattice point whose cell can reach a square field of view.
 *
 * @param {[number, number]} a1
 * @param {[number, number]} a2
 * @param {number} extent field of view, Angstrom
 * @param {number} [margin=0] how far outside the field to keep going, so that
 *   cells hanging over the edge still contribute
 * @returns {{x: number, y: number}[]}
 */
export function latticeOrigins(a1, a2, extent, margin = 0) {
  const {b1, b2} = reciprocalVectors(a1, a2);
  // The corners of the field, in lattice coordinates: the range of whole cells
  // that can touch it is the box around them.
  const corners = [[-margin, -margin], [extent + margin, -margin],
                   [-margin, extent + margin], [extent + margin, extent + margin]];
  let iLo = Infinity, iHi = -Infinity, jLo = Infinity, jHi = -Infinity;
  for (const [x, y] of corners) {
    const i = x * b1[0] + y * b1[1];
    const j = x * b2[0] + y * b2[1];
    iLo = Math.min(iLo, i); iHi = Math.max(iHi, i);
    jLo = Math.min(jLo, j); jHi = Math.max(jHi, j);
  }
  const origins = [];
  for (let i = Math.floor(iLo) - 1; i <= Math.ceil(iHi) + 1; i++) {
    for (let j = Math.floor(jLo) - 1; j <= Math.ceil(jHi) + 1; j++) {
      origins.push({x: i * a1[0] + j * a2[0], y: i * a1[1] + j * a2[1], i, j});
    }
  }
  return origins;
}

/**
 * Atoms onto a grid, as Gaussians weighted by how deep they lie.
 *
 * The weight is what makes this a surface measurement rather than a bulk one:
 * an electron that scatters elastically and comes back out carries information
 * from the first few Angstroms, so atoms below that contribute less. Depth is
 * the third component of each atom, measured downwards from the topmost one.
 *
 * Atoms outside the field are dropped rather than wrapped. A surface patch is
 * not a periodic cell, and pretending it is puts a seam across the middle of
 * every pattern.
 *
 * @param {object} options
 * @param {number} options.n grid size
 * @param {number} options.extent field of view, Angstrom
 * @param {ArrayLike<[number, number, number]>} options.atoms `[x, y, depth]`, Angstrom
 * @param {number} [options.sigma=1.1] width of one atom, Angstrom
 * @param {number} [options.escapeDepth=3] Angstrom; `Infinity` weights every atom alike
 * @param {Float64Array} [out]
 * @returns {Float64Array} n*n, row-major `[ix * n + iy]`
 */
export function splatAtoms({n, extent, atoms, sigma = 1.1, escapeDepth = 3},
                           out = new Float64Array(n * n)) {
  out.fill(0);
  const sampling = extent / n;
  const spread = sigma / sampling;                  // in grid points
  const reach = Math.max(1, Math.ceil(3 * spread));
  const twoSpreadSq = 2 * spread * spread;
  for (const atom of atoms) {
    const cx = atom[0] / sampling;
    const cy = atom[1] / sampling;
    const weight = Number.isFinite(escapeDepth) ? Math.exp(-(atom[2] ?? 0) / escapeDepth) : 1;
    const x0 = Math.round(cx);
    const y0 = Math.round(cy);
    for (let dx = -reach; dx <= reach; dx++) {
      const ix = x0 + dx;
      if (ix < 0 || ix >= n) continue;
      const rx = ix - cx;
      for (let dy = -reach; dy <= reach; dy++) {
        const iy = y0 + dy;
        if (iy < 0 || iy >= n) continue;
        const ry = iy - cy;
        out[ix * n + iy] += weight * Math.exp(-(rx * rx + ry * ry) / twoSpreadSq);
      }
    }
  }
  return out;
}

/**
 * A focused beam, as the window it puts on the surface. Gaussian, normalised to
 * sum to one so the pattern's total intensity does not depend on where it sits.
 *
 * @param {object} options
 * @param {number} options.n grid size
 * @param {number} options.extent field of view, Angstrom
 * @param {number} options.x centre along the first axis, Angstrom
 * @param {number} options.y centre along the second axis, Angstrom
 * @param {number} [options.sigma=12] beam width, Angstrom. This is the
 *   resolution of the scan, and `1 / (2 pi sigma)` is how far each spot spreads
 *   — keep it well above the superstructure cell, or neighbouring fractional
 *   spots overlap and the reconstruction stops being readable however finely
 *   you scan.
 * @param {Float64Array} [out]
 * @returns {Float64Array} n*n, summing to 1
 */
export function beamWindow({n, extent, x, y, sigma = 12}, out = new Float64Array(n * n)) {
  out.fill(0);
  const sampling = extent / n;
  const spread = sigma / sampling;
  const reach = Math.ceil(4 * spread);
  const twoSpreadSq = 2 * spread * spread;
  const cx = x / sampling;
  const cy = y / sampling;
  const x0 = Math.round(cx);
  const y0 = Math.round(cy);
  let total = 0;
  for (let dx = -reach; dx <= reach; dx++) {
    const ix = x0 + dx;
    if (ix < 0 || ix >= n) continue;
    const rx = ix - cx;
    for (let dy = -reach; dy <= reach; dy++) {
      const iy = y0 + dy;
      if (iy < 0 || iy >= n) continue;
      const ry = iy - cy;
      const value = Math.exp(-(rx * rx + ry * ry) / twoSpreadSq);
      out[ix * n + iy] = value;
      total += value;
    }
  }
  if (total > 0) for (let i = 0; i < out.length; i++) out[i] /= total;
  return out;
}

/**
 * The far-field intensity from the illuminated patch — kinematical, which is the
 * honest approximation for a single scattering layer and no more than that.
 *
 * @param {Float64Array} density from {@link splatAtoms}
 * @param {Float64Array} window from {@link beamWindow}
 * @param {number} n
 * @param {{re: Float64Array, im: Float64Array}} [scratch] reused between positions
 * @param {Float64Array} [out]
 * @returns {Float64Array} n*n intensities, corner-centered
 */
export function diffractionPattern(density, window, n, scratch = complex(n * n),
                                   out = new Float64Array(n * n)) {
  const {re, im} = scratch;
  for (let i = 0; i < re.length; i++) {
    re[i] = density[i] * window[i];
    im[i] = 0;
  }
  fft2(scratch, n, n);
  for (let i = 0; i < out.length; i++) out[i] = re[i] * re[i] + im[i] * im[i];
  return out;
}

/**
 * The same pattern, from a patch around the beam rather than the whole field.
 *
 * A scan is many patterns from one surface, and transforming the entire field of
 * view for each of them is work thrown away: past four sigma the beam is down by
 * a factor of ten thousand, so everything further out contributes nothing to the
 * pattern and everything to the bill. Transforming a patch instead costs what the
 * patch costs — and `patch` is a power of two, which the FFT cares about a great
 * deal more than its size. (A 224-point transform is slower than a 256-point one,
 * because 224 is not a power of two and has to go the long way round.)
 *
 * Reciprocal sampling is coarser in proportion: the patch is the field of view
 * now, so one sample is `1 / (patch * sampling)`, which the returned `extent`
 * reports for {@link spotSum}. Keep the patch at least eight sigma across, or the
 * beam is being cut off rather than the empty surface around it.
 *
 * Near the edge of the field the patch is slid inwards rather than reading
 * outside, which leaves the beam off centre within it. That costs nothing: an
 * off-centre beam puts a phase ramp on the spectrum, and this returns intensity.
 *
 * The window is not normalised — the measure it feeds is a ratio.
 *
 * @param {object} options
 * @param {Float64Array} options.density from {@link splatAtoms}
 * @param {number} options.n grid size of the field
 * @param {number} options.extent field of view, Angstrom
 * @param {number} options.x beam centre along the first axis, Angstrom
 * @param {number} options.y beam centre along the second axis
 * @param {number} [options.sigma=12] beam width, Angstrom
 * @param {number} [options.patch=128] patch size in grid points; a power of two
 * @param {{re: Float64Array, im: Float64Array}} [scratch] reused between positions
 * @param {Float64Array} [out]
 * @returns {{pattern: Float64Array, n: number, extent: number}} corner-centered
 */
export function illuminatedPattern({density, n, extent, x, y, sigma = 12, patch = 128},
                                   scratch = complex(patch * patch),
                                   out = new Float64Array(patch * patch)) {
  const sampling = extent / n;
  const size = Math.min(patch, n);
  const half = size >> 1;
  const x0 = Math.max(0, Math.min(n - size, Math.round(x / sampling) - half));
  const y0 = Math.max(0, Math.min(n - size, Math.round(y / sampling) - half));
  const spread = sigma / sampling;
  const twoSpreadSq = 2 * spread * spread;

  // Separable, so this is 2n exponentials rather than n^2 of them.
  const wx = new Float64Array(size);
  const wy = new Float64Array(size);
  for (let i = 0; i < size; i++) {
    const rx = x0 + i - x / sampling;
    const ry = y0 + i - y / sampling;
    wx[i] = Math.exp(-(rx * rx) / twoSpreadSq);
    wy[i] = Math.exp(-(ry * ry) / twoSpreadSq);
  }

  const {re, im} = scratch;
  for (let ix = 0; ix < size; ix++) {
    const row = (x0 + ix) * n + y0;
    for (let iy = 0; iy < size; iy++) {
      const i = ix * size + iy;
      re[i] = wx[ix] * wy[iy] * density[row + iy];
      im[i] = 0;
    }
  }
  fft2(scratch, size, size);
  for (let i = 0; i < size * size; i++) out[i] = re[i] * re[i] + im[i] * im[i];
  return {pattern: out, n: size, extent: size * sampling};
}

/**
 * Which spots belong to the superstructure and which to the substrate under it.
 *
 * For an `m x m` reconstruction the substrate's own spots are the ones at whole
 * multiples of `m` in both indices; every other spot exists only because the
 * surface rearranged.
 *
 * @param {[number, number]} b1 reciprocal vectors of the SUPERSTRUCTURE cell, 1/Å
 * @param {[number, number]} b2
 * Spots too close to the direct beam are left out. However sharply a focused
 * beam is defined, its transform has a skirt, and within a step or two of the
 * origin that skirt is brighter than anything the surface itself does — so the
 * innermost fractional spots measure the beam rather than the reconstruction.
 * A real instrument has the same problem and answers it the same way.
 *
 * @param {object} options
 * @param {number} options.multiple the `m` of `m x m` (7 for Si(111))
 * @param {number} [options.rings=1] how many substrate orders out to go
 * @param {number} [options.exclude=1.5] drop spots within this many
 *   superstructure steps of the direct beam
 * @returns {{integer: [number, number][], fractional: [number, number][]}} in 1/Å
 */
export function superstructureSpots(b1, b2, {multiple, rings = 1, exclude = 1.5}) {
  const integer = [];
  const fractional = [];
  const reach = multiple * rings;
  const step = Math.min(Math.hypot(...b1), Math.hypot(...b2));
  const limit = rings * multiple * step * 1.01;
  for (let h = -reach; h <= reach; h++) {
    for (let k = -reach; k <= reach; k++) {
      if (h === 0 && k === 0) continue;   // the direct beam says nothing about the surface
      const q = [h * b1[0] + k * b2[0], h * b1[1] + k * b2[1]];
      const radius = Math.hypot(...q);
      if (radius > limit || radius < exclude * step) continue;
      if (h % multiple === 0 && k % multiple === 0) integer.push(q);
      else fractional.push(q);
    }
  }
  return {integer, fractional};
}

/**
 * How much of the pattern sits in a given set of spots, above the background
 * around them.
 *
 * The background matters more than it sounds. Anything that interrupts the
 * surface's periodicity — a boundary, a step, a defect — scatters diffusely, and
 * diffuse intensity lands everywhere, including on every spot position. Counted
 * raw, a beam sitting on a domain boundary reads as *more* ordered than one in
 * the middle of a domain, which is the opposite of the truth. Subtracting what
 * the annulus around each spot holds leaves the part that is genuinely a spot.
 *
 * @param {Float64Array} pattern corner-centered, from {@link diffractionPattern}
 * @param {number} n
 * @param {number} extent field of view, Angstrom — one DFT sample is `1 / extent`
 * @param {ArrayLike<[number, number]>} spots in 1/Å
 * @param {object} [options]
 * @param {number} [options.radius=2] how far around each spot to collect, in samples
 * @param {number} [options.background=4] outer radius of the annulus the local
 *   background is read from; 0 turns the subtraction off
 * @returns {number} never negative
 */
export function spotSum(pattern, n, extent, spots, {radius = 2, background = 4} = {}) {
  let sum = 0;
  const reach = Math.ceil(Math.max(radius, background));
  const inside = radius * radius;
  const outside = background * background;
  for (const spot of spots) {
    const px = Math.round(spot[0] * extent);
    const py = Math.round(spot[1] * extent);
    let peak = 0;
    let peakCount = 0;
    let ring = 0;
    let ringCount = 0;
    for (let dx = -reach; dx <= reach; dx++) {
      for (let dy = -reach; dy <= reach; dy++) {
        const r2 = dx * dx + dy * dy;
        const value = pattern[wrapIndex(px + dx, n) * n + wrapIndex(py + dy, n)];
        if (r2 <= inside) {
          peak += value;
          peakCount += 1;
        } else if (background > radius && r2 <= outside) {
          ring += value;
          ringCount += 1;
        }
      }
    }
    sum += Math.max(0, peak - (ringCount > 0 ? (ring / ringCount) * peakCount : 0));
  }
  return sum;
}

const wrapIndex = (i, n) => ((i % n) + n) % n;
