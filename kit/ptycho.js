// Direct (non-iterative) STEM phase retrieval: the aperture-overlap function and
// the family of estimators built on it.
//
// Under the weak phase object approximation the specimen transmission is
// `t(r) ~ 1 + i sigma V(r)`, and the diffraction intensities, Fourier
// transformed over the scan coordinate, separate into
//
//     G(q,k) = |psi(k)|^2 delta(q) + Gamma(q,k) phi(q)
//     Gamma(q,k) = psi*(k) psi(k-q) - psi(k) psi*(k+q)
//     psi(k) = A(k) exp(-i chi(k))
//
// so every direct method is an inversion of one linear operator. SSB, OBF,
// matched filter, parallax and iCOM differ only in the kernel `H(q,k)` they
// multiply `G` by before accumulating. That is the organising idea of this
// module, and of Algorithm 1 in the lab's MRS Bulletin review.
//
// Conventions follow the rest of the kit: corner-centered (fftfreq) grids,
// complex arrays as `{re, im}` pairs of Float64Array, row-major `[ix * ny + iy]`.
//
// One geometric assumption runs through everything here: the **scan field of
// view equals the object field of view**, so the scan frequency spacing equals
// the detector's own. The reconstruction grid is then a `q x q` sub-block of the
// probe's own `n x n` grid, `Gamma` is integer index arithmetic, and no
// interpolation is needed anywhere.

import {complex, fft2, ifft2, fftshift, ifftshift} from "./fft.js";
import {gradient2d} from "./image.js";

// ---------------------------------------------------------------------------
// The aperture overlap function
// ---------------------------------------------------------------------------

/**
 * `Gamma(q, k)` for one fixed detector pixel, over the whole `q` grid.
 *
 * Both shifted probes are read at wrapped indices. Wrapping is harmless — and
 * the reason this needs no interpolation — provided the aperture radius stays
 * below a quarter of the grid, so that a frequency pushed outside the grid
 * always lands back outside the aperture. {@link checkOverlapGeometry} asserts
 * it; the widgets size their apertures to satisfy it.
 *
 * @param {{re: Float64Array, im: Float64Array}} probe corner-centered, Fourier space
 * @param {number} n probe grid size (square)
 * @param {number} q output grid size, `q <= n`
 * @param {number} kIndex flat index of the detector pixel in the `n x n` grid
 * @param {{re: Float64Array, im: Float64Array}} [out]
 */
export function overlapFunction(probe, n, q, kIndex, out = complex(q * q)) {
  const ki = (kIndex / n) | 0;
  const kj = kIndex % n;
  // Signed frequency indices, so that adding and subtracting q is arithmetic on
  // frequencies rather than on array positions.
  const si = ki < n / 2 ? ki : ki - n;
  const sj = kj < n / 2 ? kj : kj - n;

  const pr = probe.re[kIndex];
  const pi = probe.im[kIndex];
  const half = q / 2;

  for (let a = 0; a < q; a++) {
    const fa = a < half ? a : a - q;
    const minusRow = (((si - fa) % n) + n) % n;
    const plusRow = (((si + fa) % n) + n) % n;
    const mRowBase = minusRow * n;
    const pRowBase = plusRow * n;
    const outBase = a * q;

    for (let b = 0; b < q; b++) {
      const fb = b < half ? b : b - q;
      const mIdx = mRowBase + ((((sj - fb) % n) + n) % n);
      const pIdx = pRowBase + ((((sj + fb) % n) + n) % n);

      // psi*(k) psi(k-q)
      const ar = probe.re[mIdx];
      const ai = probe.im[mIdx];
      const firstRe = pr * ar + pi * ai;
      const firstIm = pr * ai - pi * ar;

      // psi(k) psi*(k+q)
      const br = probe.re[pIdx];
      const bi = probe.im[pIdx];
      const secondRe = pr * br + pi * bi;
      const secondIm = pi * br - pr * bi;

      const o = outBase + b;
      out.re[o] = firstRe - secondRe;
      out.im[o] = firstIm - secondIm;
    }
  }
  return out;
}

/**
 * `Gamma(q, k)` for one fixed spatial frequency, over the whole detector plane.
 *
 * The complementary slice to {@link overlapFunction}, and the one the geometry
 * lives in: all three apertures — `A(k)` at the origin, `A(k-q)` and `A(k+q)`
 * displaced by `+q` and `-q` — are circles in this plane, so the double-overlap
 * regions and the triple overlap where they meet are visible directly.
 *
 * @param {{re, im}} probe corner-centered, Fourier space
 * @param {number} n grid size (square)
 * @param {number} qIndex flat index of the spatial frequency
 */
export function overlapFunctionAtQ(probe, n, qIndex, out = complex(n * n)) {
  const qi = (qIndex / n) | 0;
  const qj = qIndex % n;
  const dq0 = qi < n / 2 ? qi : qi - n;
  const dq1 = qj < n / 2 ? qj : qj - n;

  for (let a = 0; a < n; a++) {
    const minusRow = (((a - dq0) % n) + n) % n;
    const plusRow = (((a + dq0) % n) + n) % n;
    const rowBase = a * n;
    const mRowBase = minusRow * n;
    const pRowBase = plusRow * n;

    for (let b = 0; b < n; b++) {
      const i = rowBase + b;
      const mIdx = mRowBase + ((((b - dq1) % n) + n) % n);
      const pIdx = pRowBase + ((((b + dq1) % n) + n) % n);

      const pr = probe.re[i];
      const pi = probe.im[i];
      const ar = probe.re[mIdx];
      const ai = probe.im[mIdx];
      const br = probe.re[pIdx];
      const bi = probe.im[pIdx];

      out.re[i] = (pr * ar + pi * ai) - (pr * br + pi * bi);
      out.im[i] = (pr * ai - pi * ar) - (pi * br - pr * bi);
    }
  }
  return out;
}

/**
 * Where the shifted apertures overlap, as masks over the detector plane.
 *
 * `double` marks the two regions reached by `A(k)` and exactly one of its
 * displaced copies — the pair of lens shapes that carry the signal, and the
 * reason the literature calls them trotters. `triple` marks where all three
 * meet, which happens only below the aperture radius and which cancels exactly
 * in focus, carrying no information there.
 *
 * @param {Float64Array} aperture `A(k)`, corner-centered
 * @param {number} n
 * @param {number} qIndex
 * @param {number} [level=0.5] what counts as inside a soft-edged aperture
 * @returns {{double: Uint8Array, triple: Uint8Array}}
 */
export function overlapRegions(aperture, n, qIndex, level = 0.5) {
  const qi = (qIndex / n) | 0;
  const qj = qIndex % n;
  const dq0 = qi < n / 2 ? qi : qi - n;
  const dq1 = qj < n / 2 ? qj : qj - n;

  const double = new Uint8Array(n * n);
  const triple = new Uint8Array(n * n);

  for (let a = 0; a < n; a++) {
    const minusRow = (((a - dq0) % n) + n) % n;
    const plusRow = (((a + dq0) % n) + n) % n;
    for (let b = 0; b < n; b++) {
      const i = a * n + b;
      if (aperture[i] <= level) continue;
      const inMinus = aperture[minusRow * n + ((((b - dq1) % n) + n) % n)] > level;
      const inPlus = aperture[plusRow * n + ((((b + dq1) % n) + n) % n)] > level;
      if (inMinus && inPlus) triple[i] = 1;
      else if (inMinus || inPlus) double[i] = inMinus ? 1 : 2;
    }
  }
  return {double, triple};
}

/**
 * Throws unless the aperture is small enough that index wrapping in
 * {@link overlapFunction} cannot alias one aperture onto another.
 *
 * A frequency `|m| <= radius + q/2` is wrapped to `m - n`. For that to stay
 * outside the aperture we need `radius + q/2 < n - radius`, and since `q <= n`
 * the safe condition is `4 * radius < n`.
 *
 * @param {number} n grid size
 * @param {number} q reconstruction grid size
 * @param {number} radiusPixels aperture radius, in pixels of the `n` grid
 */
export function checkOverlapGeometry(n, q, radiusPixels) {
  if (q > n) throw new Error(`reconstruction grid q=${q} exceeds probe grid n=${n}`);
  if (n % q !== 0) throw new Error(`q=${q} must divide n=${n} so the grids share a spacing`);
  if (4 * radiusPixels >= n) {
    throw new Error(
      `aperture radius ${radiusPixels} px is too large for a ${n} px grid: ` +
        `the overlap function needs 4 * radius < n, so use at most ${Math.ceil(n / 4) - 1} px`
    );
  }
}

/**
 * Flat indices of the detector pixels inside the aperture.
 *
 * Ordered by scattering angle by default, so a reconstruction that consumes
 * them in order spirals outwards from the optic axis. That is both better to
 * watch and truer to what is happening: the pixels nearest the axis carry the
 * low spatial frequencies, so the image arrives coarse and sharpens, rather
 * than sweeping in as a raster of unrelated detail.
 *
 * @param {Float64Array} alpha scattering angle per detector pixel [rad]
 * @param {number} semiangleCutoff [mrad]
 * @param {"radial"|"raster"} [order="radial"]
 */
export function brightFieldIndices(alpha, semiangleCutoff, order = "radial") {
  const cutoff = semiangleCutoff * 1e-3;
  const out = [];
  for (let i = 0; i < alpha.length; i++) if (alpha[i] < cutoff) out.push(i);
  if (order === "radial") out.sort((i, j) => alpha[i] - alpha[j]);
  return Int32Array.from(out);
}

/**
 * `sum_k |Gamma(q,k)|` and `sum_k |Gamma(q,k)|^2` over the bright field.
 *
 * Both depend only on the probe — never on the object, the dose or the data —
 * so a caller can cache them on the aberrations, the aperture and the two grid
 * sizes. They are what the analytic contrast transfer functions are built from.
 *
 * `counts` is how many detector pixels contribute at each spatial frequency,
 * which is the double-overlap area counted twice plus the triple-overlap area.
 * It is the SSB variance: SSB adds the signal coherently and the noise
 * incoherently, so its noise grows as the square root of this.
 *
 * @param {number} [floor=1e-12] below this a `|Gamma|` counts as zero
 * @returns {{sum: Float64Array, energy: Float64Array, counts: Float64Array}}
 */
export function overlapSums(probe, n, q, bfIndices, {floor = 1e-12} = {}) {
  const sum = new Float64Array(q * q);
  const energy = new Float64Array(q * q);
  const counts = new Float64Array(q * q);
  const gamma = complex(q * q);

  for (let j = 0; j < bfIndices.length; j++) {
    overlapFunction(probe, n, q, bfIndices[j], gamma);
    for (let i = 0; i < sum.length; i++) {
      const magnitude2 = gamma.re[i] * gamma.re[i] + gamma.im[i] * gamma.im[i];
      const magnitude = Math.sqrt(magnitude2);
      sum[i] += magnitude;
      energy[i] += magnitude2;
      if (magnitude > floor) counts[i] += 1;
    }
  }
  return {sum, energy, counts};
}

// ---------------------------------------------------------------------------
// Fourier tiling
// ---------------------------------------------------------------------------

/**
 * Replicate an `m x m` scan spectrum `f` times along each axis, onto a
 * `q = m * f` grid, in fftfreq order.
 *
 * The scan samples the object on a coarser grid than the object itself, so its
 * spectrum is periodic with period `m`. Tiling states that periodicity
 * explicitly, which is what lets a direct method reconstruct onto a finer grid
 * than the scan step — the kernel differs across the tiles even though the data
 * does not.
 *
 * @param {{re: Float64Array, im: Float64Array}} g the `m x m` spectrum
 * @param {number} m
 * @param {number} f integer upsampling factor
 */
export function tileSpectrum(g, m, f, out = complex(m * f * m * f)) {
  const q = m * f;
  for (let a = 0; a < q; a++) {
    // fftfreq order tiles by taking the source frequency modulo m, which for a
    // corner-centered layout is exactly the array index modulo m.
    const srcRow = (a % m) * m;
    const dstRow = a * q;
    for (let b = 0; b < q; b++) {
      const src = srcRow + (b % m);
      out.re[dstRow + b] = g.re[src];
      out.im[dstRow + b] = g.im[src];
    }
  }
  return out;
}

/**
 * How far each bright-field pixel's virtual image is displaced, in pixels of
 * the reconstruction grid.
 *
 * This is the whole content of the parallax approximation. A virtual
 * bright-field image formed from an off-axis detector pixel, out of focus,
 * appears shifted sideways; putting them all back on top of each other and
 * adding them is the reconstruction. Formally the shift is the gradient of the
 * aberration surface,
 *
 *     dr = grad_k chi(k) / (2 pi)     [Angstrom]
 *
 * which follows from matching `exp[-i grad_k chi . q]` against the Fourier
 * shift `exp[-2 pi i q . dr]`. For pure defocus it reduces to
 * `dr = lambda C10 k`, the familiar statement that a beam tilted by `lambda k`
 * misses by the defocus times the tilt.
 *
 * Note this is not `optics.rayDisplacement`, which scales its gradient for
 * warping an image on the detector grid rather than for shifting one on the
 * reconstruction grid.
 *
 * @param {Float64Array} chiArray `chi` on the detector grid, corner-centered
 * @param {number} n detector grid size
 * @param {number} dk reciprocal sampling of the detector grid [1/Angstrom]
 * @param {number} qSampling real-space sampling of the reconstruction grid [Angstrom]
 * @returns {{shiftRow: Float64Array, shiftCol: Float64Array}} per detector pixel
 */
export function parallaxShifts(chiArray, n, dk, qSampling) {
  // Differentiate on the centred array, where neighbouring entries really are
  // neighbouring frequencies. On a corner-centered array the middle of each
  // axis is a wrap, and a central difference across it is meaningless.
  const centred = fftshift(chiArray, n, n);
  const {dRow, dCol} = gradient2d(centred, n, n, dk, dk);

  const scale = 1 / (2 * Math.PI * qSampling);
  for (let i = 0; i < dRow.length; i++) {
    dRow[i] *= scale;
    dCol[i] *= scale;
  }
  return {shiftRow: ifftshift(dRow, n, n), shiftCol: ifftshift(dCol, n, n)};
}

// ---------------------------------------------------------------------------
// The kernels
// ---------------------------------------------------------------------------

/**
 * The five direct-method kernel names, in the order the paper introduces them.
 */
export const PTYCHO_METHODS = ["ssb", "obf", "matchedFilter", "parallax", "icom"];

/** Display names, for radios and legends. */
export const PTYCHO_LABELS = {
  ssb: "SSB",
  obf: "OBF",
  matchedFilter: "matched filter",
  parallax: "parallax",
  icom: "iCOM"
};

/**
 * Whether a method's normalisation is a running sum over detector pixels
 * (`energy`) or a plain mean (`count`). See {@link directImage}.
 */
const NORMALISATION = {
  ssb: "count",
  obf: "energy-sqrt",
  matchedFilter: "energy",
  parallax: "count",
  icom: "count"
};

// ---------------------------------------------------------------------------
// Accumulation
// ---------------------------------------------------------------------------

/**
 * A Fourier-space accumulator: `num` collects `sum_k G H`, `den` collects
 * `sum_k |Gamma|^2`.
 *
 * Accumulating in Fourier space and inverse transforming once, at display time,
 * rather than transforming once per detector pixel, is about nine times faster
 * — one transform per frame instead of one per bright-field pixel. It also
 * gives a better partial result: because the normalisation is applied at
 * display time, the partial image is the exact estimate from the pixels seen so
 * far rather than a prefix of a sum that is not right until the last one lands.
 */
export function directAccumulator(q) {
  return {q, num: complex(q * q), den: new Float64Array(q * q), count: 0};
}

/** Zero an accumulator in place, keeping its buffers. */
export function resetAccumulator(acc) {
  acc.num.re.fill(0);
  acc.num.im.fill(0);
  acc.den.fill(0);
  acc.count = 0;
  return acc;
}

/**
 * Add one detector pixel's contribution. Deliberately does not transform.
 *
 * @param {object} acc from {@link directAccumulator}
 * @param {{re, im}} g the tiled scan spectrum `G(q', k)` for this pixel
 * @param {{re, im}} gamma `Gamma(q', k)` for this pixel
 * @param {string} method one of {@link PTYCHO_METHODS}
 * @param {object} ctx per-pixel geometry:
 *   `{shiftRow, shiftCol, phaseFlip, bandLimit, kRow, kCol, qRow, qCol, inverseQ2}`
 *   — `shiftRow`/`shiftCol` are `grad_k chi(k)` in pixels (parallax),
 *   `phaseFlip` is `sgn[sin chi(q')]` over the grid (parallax),
 *   `kRow`/`kCol` are the signed frequency indices of `k` (iCOM),
 *   `qRow`/`qCol`/`inverseQ2` describe the `q'` grid (iCOM).
 */
export function accumulatePixel(acc, g, gamma, method, ctx = {}) {
  const {num, den, q} = acc;
  const length = q * q;

  // The denominator is the same for every method, and OBF and MF read it back
  // at display time. Collecting it always costs one multiply-add and keeps the
  // accumulator method-agnostic, so a widget can change method mid-run.
  for (let i = 0; i < length; i++) {
    den[i] += gamma.re[i] * gamma.re[i] + gamma.im[i] * gamma.im[i];
  }

  switch (method) {
    case "ssb": {
      // H = -i Gamma* / |Gamma|
      for (let i = 0; i < length; i++) {
        const magnitude = Math.hypot(gamma.re[i], gamma.im[i]);
        if (magnitude === 0) continue;
        // -i * conj(gamma) = -i (gr - i gi) = -gi - i gr
        const hr = -gamma.im[i] / magnitude;
        const hi = -gamma.re[i] / magnitude;
        num.re[i] += g.re[i] * hr - g.im[i] * hi;
        num.im[i] += g.re[i] * hi + g.im[i] * hr;
      }
      break;
    }
    case "obf":
    case "matchedFilter": {
      // H = -i Gamma*, with the normalisation deferred to `directImage`.
      for (let i = 0; i < length; i++) {
        const hr = -gamma.im[i];
        const hi = -gamma.re[i];
        num.re[i] += g.re[i] * hr - g.im[i] * hi;
        num.im[i] += g.re[i] * hi + g.im[i] * hr;
      }
      break;
    }
    case "parallax": {
      // H = exp[-i grad_k chi(k) . q'], optionally times sgn[sin chi(q')].
      //
      // Unit modulus at every spatial frequency — a pure phase ramp, which is
      // the whole point of the approximation: it is the cheapest possible
      // kernel and needs nothing but the aberration gradient.
      //
      // A pure phase ramp: unit modulus at every spatial frequency the method
      // can transfer, which is the whole economy of it — nothing but the
      // aberration gradient is needed to build one.
      //
      // `bandLimit` is where "every frequency it can transfer" stops, in
      // detector-pixel units, and should be `2 * aperture radius`. The parallax
      // transfer function is `sin[chi]` times the aperture autocorrelation, and
      // that autocorrelation is identically zero past `2 alpha`, so nothing is
      // being thrown away. It suppresses noise rather than being load-bearing;
      // the reconstruction is sound without it. The overlap kernels get the
      // same cutoff for free, through `Gamma`.
      //
      // A constant-magnitude kernel does make this one sensitive to anything
      // spurious in the tiled spectrum, since it has no magnitude of its own to
      // suppress it with. See {@link scanSpectra}, which removes the one thing
      // that matters — the unscattered beam at zero frequency.
      const {shiftRow = 0, shiftCol = 0, phaseFlip = null, bandLimit = Infinity} = ctx;
      const twoPi = 2 * Math.PI;
      const half = q / 2;
      const band2 = bandLimit * bandLimit;
      for (let a = 0; a < q; a++) {
        const fa = a < half ? a : a - q;
        const rowPhase = (-twoPi * fa * shiftRow) / q;
        for (let b = 0; b < q; b++) {
          const i = a * q + b;
          const fb = b < half ? b : b - q;
          if (fa * fa + fb * fb > band2) continue;
          const phase = rowPhase - (twoPi * fb * shiftCol) / q;
          const flip = phaseFlip ? phaseFlip[i] : 1;
          const hr = Math.cos(phase) * flip;
          const hi = Math.sin(phase) * flip;
          num.re[i] += g.re[i] * hr - g.im[i] * hi;
          num.im[i] += g.re[i] * hi + g.im[i] * hr;
        }
      }
      break;
    }
    case "icom": {
      // H = -i (k . q') / |q'|^2
      // Band-limited for the same reason as parallax: the iCOM transfer
      // function is the probe autocorrelation, which also stops at `2 alpha`.
      const {kRow = 0, kCol = 0, qRow, qCol, inverseQ2, bandLimit = Infinity} = ctx;
      const limit2 = bandLimit * bandLimit;
      for (let i = 0; i < length; i++) {
        if (qRow[i] * qRow[i] + qCol[i] * qCol[i] > limit2) continue;
        // A real weight times -i, so the kernel is purely imaginary:
        // g * (-i w) = (g_im w) + i(-g_re w).
        const w = (kRow * qRow[i] + kCol * qCol[i]) * inverseQ2[i];
        num.re[i] += g.im[i] * w;
        num.im[i] += -g.re[i] * w;
      }
      break;
    }
    default:
      throw new Error(`unknown direct method "${method}"`);
  }

  acc.count += 1;
  return acc;
}

/**
 * Apply the method's normalisation, inverse transform, and take the real part.
 *
 * One transform, called once per animation frame on the running accumulator.
 *
 * @param {number} [epsilon=0] matched-filter regularisation, as a fraction of
 *   the largest accumulated `|Gamma|^2`
 * @returns {Float64Array} the reconstructed phase on the `q x q` grid
 */
export function directImage(acc, method, {epsilon = 0.01, out} = {}) {
  const {num, den, q, count} = acc;
  const length = q * q;
  const scratch = complex(length);

  const rule = NORMALISATION[method];
  let floor = 0;
  if (rule === "energy") {
    let peak = 0;
    for (let i = 0; i < length; i++) if (den[i] > peak) peak = den[i];
    floor = epsilon * peak;
  }

  for (let i = 0; i < length; i++) {
    let scale;
    if (rule === "count") scale = count > 0 ? 1 / count : 0;
    else if (rule === "energy-sqrt") scale = den[i] > 0 ? 1 / Math.sqrt(den[i]) : 0;
    else scale = den[i] + floor > 0 ? 1 / (den[i] + floor) : 0;
    scratch.re[i] = num.re[i] * scale;
    scratch.im[i] = num.im[i] * scale;
  }

  ifft2(scratch, q, q);
  const image = out ?? new Float64Array(length);
  for (let i = 0; i < length; i++) image[i] = scratch.re[i];
  return image;
}

// ---------------------------------------------------------------------------
// Analytic transfer
// ---------------------------------------------------------------------------

/** Real autocorrelation of a real array, via Wiener-Khinchin. */
function autocorrelation(array, n) {
  const f = {re: Float64Array.from(array), im: new Float64Array(array.length)};
  fft2(f, n, n);
  for (let i = 0; i < f.re.length; i++) {
    f.re[i] = f.re[i] * f.re[i] + f.im[i] * f.im[i];
    f.im[i] = 0;
  }
  ifft2(f, n, n);
  return f.re;
}

/**
 * The contrast transfer function of a direct method — analytic, no data.
 *
 * Ports `ctfs-figure-analytical.ipynb` cell 6, which is the reference
 * implementation for all five.
 *
 * **The probe is assumed L2-normalised**, which is what `optics.complexProbe`
 * returns. The Python reference normalises inside the CTF instead, by dividing
 * by `sum A^2`; doing both divides by the aperture weight twice and leaves
 * every CTF smaller than it should be by that factor. The aperture is still
 * needed, for the envelope parallax is built from.
 *
 * The factor of two follows the convention `psi_rec(q) = 2 phi(q) CTF(q)` used
 * throughout the review. In focus, `CTF_ssb` then peaks at 1.
 *
 * @param {string} method one of {@link PTYCHO_METHODS}
 * @param {object} context `{probe, aperture, chiArray, n, sums, phaseFlip}` — `sums`
 *   from {@link overlapSums} at `q = n`; `phaseFlip` applies only to parallax
 * @returns {Float64Array} `n * n`, corner-centered
 */
export function directCTF(method, {probe, aperture, chiArray, n, sums, epsilon = 0.01, phaseFlip = true}) {
  const length = n * n;
  const out = new Float64Array(length);
  const weight = apertureWeight(aperture);

  switch (method) {
    case "ssb":
      for (let i = 0; i < length; i++) out[i] = sums.sum[i] / 2;
      return out;
    case "obf":
      // The square root brings in one factor of the bright-field weight that the
      // coherent SSB sum does not have, so without it this lands a factor of
      // sqrt(N_BF) below every other curve and the shared axis is useless. The
      // Python reference applies the same normalisation, for the same reason.
      for (let i = 0; i < length; i++) out[i] = Math.sqrt(weight * sums.energy[i]) / 2;
      return out;
    case "matchedFilter": {
      // A ratio, so the probe normalisation cancels — but the regularisation
      // has to be measured against the same scale as the energy it guards.
      let peak = 0;
      for (let i = 0; i < length; i++) if (sums.energy[i] > peak) peak = sums.energy[i];
      const floor = epsilon * peak;
      for (let i = 0; i < length; i++) {
        out[i] = sums.energy[i] / (sums.energy[i] + floor);
      }
      return out;
    }
    case "parallax": {
      // The aperture autocorrelation envelope, from an L2-normalised aperture so
      // that it peaks at one, modulated by the axial transfer function.
      //
      // Without the phase flip this changes sign wherever `sin chi` does, which
      // is the contrast reversal parallax is known for. The flip replaces that
      // with its magnitude, leaving a transfer function of one sign that still
      // touches zero at every `sin chi` zero — the flip fixes the sign, not the
      // missing information.
      const envelope = autocorrelation(normalise(aperture), n);
      for (let i = 0; i < length; i++) {
        const transfer = -Math.sin(chiArray[i]) * envelope[i];
        out[i] = phaseFlip ? -Math.abs(Math.sin(chiArray[i])) * envelope[i] : transfer;
      }
      return out;
    }
    case "icom": {
      // The probe autocorrelation, which is the filter iCOM convolves the phase
      // with. The probe is complex and already normalised, so this is
      // |fft2(psi)|^2 transformed back.
      const f = {re: Float64Array.from(probe.re), im: Float64Array.from(probe.im)};
      fft2(f, n, n);
      for (let i = 0; i < length; i++) {
        f.re[i] = f.re[i] * f.re[i] + f.im[i] * f.im[i];
        f.im[i] = 0;
      }
      ifft2(f, n, n);
      return f.re;
    }
    default:
      throw new Error(`unknown direct method "${method}"`);
  }
}

/**
 * `sum A^2` over the aperture — the effective number of bright-field pixels.
 *
 * This is the scale the review's transfer functions are written in: its `A` is
 * a top-hat of height one, where `optics.complexProbe` divides by this weight's
 * square root. The kit works with the normalised probe throughout and puts the
 * weight back where a formula needs it, which is only OBF's transfer function
 * and SSB's noise.
 */
function apertureWeight(aperture) {
  let weight = 0;
  for (let i = 0; i < aperture.length; i++) weight += aperture[i] * aperture[i];
  return weight;
}

/** A copy of a real array scaled to unit L2 norm. */
function normalise(array) {
  let norm = 0;
  for (let i = 0; i < array.length; i++) norm += array[i] * array[i];
  norm = Math.sqrt(norm) || 1;
  return Float64Array.from(array, (v) => v / norm);
}

/**
 * The spectral signal-to-noise ratio of a direct method, under the unit-variance
 * additive noise model.
 *
 * The result worth knowing is that OBF and the matched filter come out equal:
 * the matched filter's higher CTF is exactly offset by its higher variance, so
 * the statistically reliable information is the same. SSB is close but pays for
 * treating every bright-field pixel alike.
 *
 * @param {Float64Array} [q] frequency magnitudes, required for iCOM
 */
export function directSSNR(method, context, q = null) {
  const {n, sums} = context;
  const length = n * n;
  const ctf = directCTF(method, context);
  const out = new Float64Array(length);

  switch (method) {
    case "ssb": {
      // SSB adds signal coherently and noise incoherently, so its variance is
      // simply how many detector pixels contribute at each frequency — the
      // double-overlap area twice over plus the triple overlap.
      //
      // By Cauchy-Schwarz `sum|Gamma| <= sqrt(counts * sum|Gamma|^2)`, so this
      // can never exceed the OBF result. That is the cost of weighting every
      // bright-field pixel alike.
      //
      // The aperture weight appears for the same reason it does in the OBF
      // transfer function — see {@link apertureWeight}.
      const weight = apertureWeight(context.aperture);
      for (let i = 0; i < length; i++) {
        const noise = Math.sqrt(sums.counts[i] / weight);
        out[i] = noise > 0 ? Math.abs(ctf[i]) / noise : 0;
      }
      return out;
    }
    case "obf":
    case "matchedFilter": {
      // Both normalise the variance to unity, so the SSNR is the OBF CTF.
      const obf = directCTF("obf", context);
      for (let i = 0; i < length; i++) out[i] = Math.abs(obf[i]);
      return out;
    }
    case "parallax":
      for (let i = 0; i < length; i++) out[i] = Math.abs(ctf[i]);
      return out;
    case "icom":
      if (!q) throw new Error("iCOM SSNR needs the frequency magnitudes");
      for (let i = 0; i < length; i++) out[i] = q[i] * Math.abs(ctf[i]);
      return out;
    default:
      throw new Error(`unknown direct method "${method}"`);
  }
}
