// Array statistics, contrast scaling and correlation helpers.
//
// Ported from `ctf/visualize.py:histogram_scaling`, `ctf/utils.py:radially_average_ctf`,
// and the `cross_correlate` / `sigmoid_taper` / `white_noise_object_2D` helpers
// that were duplicated across `08.calibration-disk-detection.ipynb`,
// `10.2.nanobeam-stitching.ipynb` and `12.upsampled-direct-ptychography-nb.ipynb`.

import {fft2, ifft2, complex, toComplex} from "./fft.js";

/**
 * Clip an array to its [vmin, vmax] quantiles and optionally renormalise to 0..1.
 *
 * This is the contrast scaling applied before essentially every `imshow` in the
 * lab's Python figures, so JS panels match them by default.
 *
 * Ported from `ctf/visualize.py:histogram_scaling`; also equivalent to
 * `py4DSTEM.visualize.return_scaled_histogram_ordering`.
 *
 * @param {ArrayLike<number>} array
 * @param {object} [options]
 * @param {number} [options.vmin=0.02] lower quantile
 * @param {number} [options.vmax=0.98] upper quantile
 * @param {boolean} [options.normalize=true] rescale the clipped result to 0..1
 * @returns {Float64Array}
 */
export function histogramScaling(array, {vmin = 0.02, vmax = 0.98, normalize = true} = {}) {
  const out = Float64Array.from(array);
  const n = out.length;
  if (n === 0) return out;

  let lo = Infinity, hi = -Infinity;
  for (let i = 0; i < n; i++) {
    const v = out[i];
    if (Number.isNaN(v)) continue;
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }

  if (lo === 0 && hi === 0) return out;
  if (Math.abs(hi - lo) <= 1e-8 * Math.max(Math.abs(hi), Math.abs(lo), 1e-300)) {
    if (normalize && hi !== 0) for (let i = 0; i < n; i++) out[i] /= hi;
    return out;
  }

  // numpy sorts only the non-NaN values, then indexes by rounded quantile
  const vals = Float64Array.from([...out].filter((v) => !Number.isNaN(v))).sort();
  const iMin = Math.max(0, Math.round((vals.length - 1) * vmin));
  const iMax = Math.min(vals.length - 1, Math.round((vals.length - 1) * vmax));
  const clipLo = vals[iMin];
  const clipHi = vals[iMax];

  for (let i = 0; i < n; i++) {
    if (out[i] < clipLo) out[i] = clipLo;
    else if (out[i] > clipHi) out[i] = clipHi;
  }

  if (normalize) {
    let mn = Infinity, mx = -Infinity;
    for (let i = 0; i < n; i++) {
      if (out[i] < mn) mn = out[i];
      if (out[i] > mx) mx = out[i];
    }
    for (let i = 0; i < n; i++) out[i] -= mn;
    mx -= mn;
    if (mx !== 0) for (let i = 0; i < n; i++) out[i] /= mx;
  }
  return out;
}

/**
 * Radial average of a **corner-centered** 2D array, using the fractional-weight
 * binning of `ctf/utils.py:radially_average_ctf`.
 *
 * Each sample contributes to its two neighbouring radial bins in proportion to
 * its distance from them, which is smoother than nearest-bin assignment.
 *
 * @param {ArrayLike<number>} array row-major, length nx*ny
 * @param {number} nx
 * @param {number} ny
 * @param {[number, number]} sampling real-space sampling [Angstrom]
 * @returns {{k: Float64Array, I: Float64Array}} bin centres and averages
 */
export function radialAverage(array, nx, ny, [sx, sy]) {
  const kxv = new Float64Array(nx);
  const kyv = new Float64Array(ny);
  for (let i = 0; i < nx; i++) kxv[i] = (i < (nx + 1) >> 1 ? i : i - nx) / (nx * sx);
  for (let i = 0; i < ny; i++) kyv[i] = (i < (ny + 1) >> 1 ? i : i - ny) / (ny * sy);

  const binSize = kxv.length > 1 ? kxv[1] - kxv[0] : 1;
  let kMax = 0;
  for (let ix = 0; ix < nx; ix++) {
    for (let iy = 0; iy < ny; iy++) {
      const k = Math.hypot(kxv[ix], kyv[iy]);
      if (k > kMax) kMax = k;
    }
  }

  // np.arange(0, k.max() + bin_size, bin_size)
  const nBins = Math.max(1, Math.ceil((kMax + binSize - 1e-12) / binSize));
  const weight = new Float64Array(nBins + 2);
  const total = new Float64Array(nBins + 2);

  for (let ix = 0; ix < nx; ix++) {
    for (let iy = 0; iy < ny; iy++) {
      const k = Math.hypot(kxv[ix], kyv[iy]);
      const ind = k / binSize;
      const f = Math.floor(ind);
      const d = ind - f;
      const v = array[ix * ny + iy];
      weight[f] += 1 - d;
      weight[f + 1] += d;
      total[f] += v * (1 - d);
      total[f + 1] += v * d;
    }
  }

  let kxAbsMax = 0;
  for (let i = 0; i < nx; i++) kxAbsMax = Math.max(kxAbsMax, Math.abs(kxv[i]));

  const kOut = [], iOut = [];
  for (let i = 0; i < nBins; i++) {
    const kBin = i * binSize;
    if (kBin > kxAbsMax) break;
    kOut.push(kBin);
    iOut.push(weight[i] === 0 ? NaN : total[i] / weight[i]);
  }
  return {k: Float64Array.from(kOut), I: Float64Array.from(iOut)};
}

/**
 * A radially symmetric sigmoid taper, used to soften a probe template before
 * correlation so its edges do not ring.
 *
 * Ported from `sigmoid_taper` in `08.calibration-disk-detection.ipynb`.
 */
export function sigmoidTaper(k, cutoff, width) {
  const out = new Float64Array(k.length);
  for (let i = 0; i < k.length; i++) {
    out[i] = 1 / (1 + Math.exp((k[i] - cutoff) / width));
  }
  return out;
}

/**
 * Cross- / phase-correlation of an image with a template.
 *
 * `power` interpolates between plain cross-correlation (0) and pure phase
 * correlation (1), matching the `correlation power` slider in the disk-detection
 * widgets (`08.calibration-disk-detection.ipynb`, `10.2.nanobeam-stitching.ipynb`).
 *
 * @param {ArrayLike<number>} image row-major, nx*ny
 * @param {ArrayLike<number>} template same shape
 * @returns {Float64Array} the real part of the correlogram, corner-centered
 */
export function crossCorrelate(image, template, nx, ny, power = 1) {
  const F = fft2(toComplex(image), nx, ny);
  const G = fft2(toComplex(template), nx, ny);
  const n = nx * ny;
  const P = complex(n);

  for (let i = 0; i < n; i++) {
    // F * conj(G)
    const re = F.re[i] * G.re[i] + F.im[i] * G.im[i];
    const im = F.im[i] * G.re[i] - F.re[i] * G.im[i];
    if (power === 0) {
      P.re[i] = re;
      P.im[i] = im;
    } else {
      const mag = Math.hypot(re, im);
      const scale = mag === 0 ? 0 : Math.pow(mag, -power);
      P.re[i] = re * scale;
      P.im[i] = im * scale;
    }
  }
  ifft2(P, nx, ny);
  return P.re;
}

/**
 * A random-phase, Hermitian-symmetric "white noise" object — a synthetic sample
 * whose power spectrum is flat, used to visualise a CTF without needing real data.
 *
 * Ported from `white_noise_object_2D`, duplicated in
 * `08.analytical-ctf-ingredients-nb.ipynb` and `12.upsampled-direct-ptychography-nb.ipynb`.
 *
 * @param {number} n grid size (square)
 * @param {number} [scale=1] standard deviation of the resulting phase
 * @param {() => number} [random=Math.random]
 */
export function whiteNoiseObject2D(n, scale = 1, random = Math.random) {
  const size = n * n;
  const spectrum = complex(size);
  for (let i = 0; i < size; i++) {
    const phase = 2 * Math.PI * random();
    spectrum.re[i] = Math.cos(phase);
    spectrum.im[i] = Math.sin(phase);
  }
  // enforce Hermitian symmetry so the inverse transform is real
  for (let ix = 0; ix < n; ix++) {
    for (let iy = 0; iy < n; iy++) {
      const i = ix * n + iy;
      const j = ((n - ix) % n) * n + ((n - iy) % n);
      if (j < i) {
        spectrum.re[i] = spectrum.re[j];
        spectrum.im[i] = -spectrum.im[j];
      }
    }
  }
  spectrum.im[0] = 0;
  ifft2(spectrum, n, n);

  let mean = 0;
  for (let i = 0; i < size; i++) mean += spectrum.re[i];
  mean /= size;
  let variance = 0;
  for (let i = 0; i < size; i++) variance += (spectrum.re[i] - mean) ** 2;
  const std = Math.sqrt(variance / size) || 1;

  const out = new Float64Array(size);
  for (let i = 0; i < size; i++) out[i] = ((spectrum.re[i] - mean) / std) * scale;
  return out;
}

/**
 * Sample from a Poisson distribution — Knuth's method for small means, and a
 * normal approximation above 30 where Knuth becomes slow and the approximation
 * is accurate to well under a count.
 *
 * Used for shot noise, which the Python gets from `abtem`'s `.poisson_noise(dose)`.
 */
export function poisson(lambda, random = Math.random) {
  if (lambda <= 0) return 0;
  if (lambda < 30) {
    const limit = Math.exp(-lambda);
    let k = 0;
    let p = 1;
    do {
      k++;
      p *= random();
    } while (p > limit);
    return k - 1;
  }
  // Box-Muller normal approximation, rounded and floored at zero
  const u = Math.max(random(), Number.MIN_VALUE);
  const v = random();
  const z = Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  return Math.max(0, Math.round(lambda + Math.sqrt(lambda) * z));
}

/** Apply Poisson shot noise to an intensity array at a given dose. */
export function applyShotNoise(intensity, dose, random = Math.random) {
  const out = new Float64Array(intensity.length);
  for (let i = 0; i < intensity.length; i++) {
    out[i] = poisson(Math.max(0, intensity[i]) * dose, random) / dose;
  }
  return out;
}

/**
 * Integrate a gradient field back to the scalar it came from, in Fourier space.
 *
 * This is the step that turns a centre-of-mass vector field into a phase — iCoM,
 * or integrated DPC. Solving `∇·g = ∇²T` gives
 *
 *     T̂(k) = (k₀·ĝ₀ + k₁·ĝ₁) / (2πi |k|²)
 *
 * Ported from abtem's `_integrate_gradient_2d` (`measurements.py:2239`),
 * including its DC guard (|k|² floored at 1e-12) and its final shift so the
 * minimum is zero. There is no other regularisation: the result is dominated by
 * low frequencies and will look noisy until a scan is reasonably complete.
 *
 * @param {ArrayLike<number>} g0 gradient along the first axis
 * @param {ArrayLike<number>} g1 gradient along the second axis
 * @param {number} nx
 * @param {number} ny
 * @param {[number, number]} sampling
 * @returns {Float64Array} the integrated scalar field, minimum zero
 */
export function integrateGradient(g0, g1, nx, ny, [sx, sy]) {
  const k0 = new Float64Array(nx);
  const k1 = new Float64Array(ny);
  for (let i = 0; i < nx; i++) k0[i] = (i < (nx + 1) >> 1 ? i : i - nx) / (nx * sx);
  for (let i = 0; i < ny; i++) k1[i] = (i < (ny + 1) >> 1 ? i : i - ny) / (ny * sy);

  const F0 = fft2(toComplex(g0), nx, ny);
  const F1 = fft2(toComplex(g1), nx, ny);

  const out = complex(nx * ny);
  for (let ix = 0; ix < nx; ix++) {
    for (let iy = 0; iy < ny; iy++) {
      const i = ix * ny + iy;
      const k2 = Math.max(k0[ix] * k0[ix] + k1[iy] * k1[iy], 1e-12);
      const nr = F0.re[i] * k0[ix] + F1.re[i] * k1[iy];
      const ni = F0.im[i] * k0[ix] + F1.im[i] * k1[iy];
      // Divide by 2πi|k|²: multiplying by -i swaps the parts and negates one.
      const d = 2 * Math.PI * k2;
      out.re[i] = ni / d;
      out.im[i] = -nr / d;
    }
  }
  ifft2(out, nx, ny);

  let min = Infinity;
  for (let i = 0; i < out.re.length; i++) if (out.re[i] < min) min = out.re[i];
  const result = new Float64Array(out.re.length);
  for (let i = 0; i < result.length; i++) result[i] = out.re[i] - min;
  return result;
}

/**
 * Bin a 2D array by an integer factor, averaging each block.
 *
 * The cheapest way to buy frame budget: binning by 2 quarters the pixel count
 * and cuts an FFT-bound render to roughly a quarter of its cost, which is what
 * makes a live multislice scene animate at all.
 *
 * @param {ArrayLike<number>} array row-major, nx*ny
 * @param {number} [factor=2]
 * @returns {{data: Float64Array, nx: number, ny: number}}
 */
export function binArray(array, nx, ny, factor = 2) {
  const outX = Math.floor(nx / factor);
  const outY = Math.floor(ny / factor);
  const data = new Float64Array(outX * outY);
  const norm = 1 / (factor * factor);
  for (let ix = 0; ix < outX; ix++) {
    for (let iy = 0; iy < outY; iy++) {
      let sum = 0;
      for (let dx = 0; dx < factor; dx++) {
        const row = (ix * factor + dx) * ny;
        for (let dy = 0; dy < factor; dy++) sum += array[row + iy * factor + dy];
      }
      data[ix * outY + iy] = sum * norm;
    }
  }
  return {data, nx: outX, ny: outY};
}

/**
 * Crop a `size` x `size` window whose ORIGIN is at (row, col), wrapping at the
 * edges.
 *
 * The origin, rather than the centre, is at the probe: a probe from `ifft2` sits
 * at index (0,0), so the two line up with no shifting.
 *
 * **This is an approximation, and not always a small one.** Propagating a probe
 * through a window instead of the full field is cheap — a 64x64 window costs
 * 13x less than a 244x242 field — and it is exact only while the probe's support
 * stays within `size / 2` of the origin. A converged probe looks compact, but a
 * hard aperture gives it Airy tails that fall off as 1/r, and on the torus the
 * window defines those tails wrap back onto the specimen.
 *
 * Measured on a 22 Å gold particle at 0.2 Å sampling, a 20 mrad probe and 80 kV:
 * a 64-wide window (12.8 Å) **inverted** the sign of the integrated centre of
 * mass, and 96 still did; only at 128 (25.6 Å) did the contrast come out the
 * right way round, and even then the image correlated poorly with the full-field
 * answer. Nothing about the result looks wrong — it is a plausible image of the
 * right object with the contrast reversed. If a widget can afford the full grid,
 * it should use the full grid.
 */
export function cropWrap(array, nx, ny, row, col, size) {
  const out = new Float64Array(size * size);
  for (let i = 0; i < size; i++) {
    const src = (((row + i) % nx) + nx) % nx;
    for (let j = 0; j < size; j++) {
      out[i * size + j] = array[src * ny + ((((col + j) % ny) + ny) % ny)];
    }
  }
  return out;
}

/** {@link cropWrap} for a complex `{re, im}` array. */
export function cropWrapComplex({re, im}, nx, ny, row, col, size) {
  return {
    re: cropWrap(re, nx, ny, row, col, size),
    im: cropWrap(im, nx, ny, row, col, size)
  };
}

/**
 * Crop the centre out of a 2D array.
 *
 * A focused probe occupies a handful of pixels in a 244 x 242 field, and a
 * diffraction pattern's useful content sits near the origin — the Python
 * notebooks crop both before display (e.g. `[74:-74, 73:-73]`).
 *
 * @param {ArrayLike<number>} array row-major, nx*ny
 * @returns {{data: Float64Array, nx: number, ny: number}}
 */
export function cropCenter(array, nx, ny, cropX, cropY = cropX) {
  const outX = Math.min(cropX, nx);
  const outY = Math.min(cropY, ny);
  const x0 = (nx - outX) >> 1;
  const y0 = (ny - outY) >> 1;
  const data = new Float64Array(outX * outY);
  for (let ix = 0; ix < outX; ix++) {
    for (let iy = 0; iy < outY; iy++) {
      data[ix * outY + iy] = array[(ix + x0) * ny + (iy + y0)];
    }
  }
  return {data, nx: outX, ny: outY};
}

/** {@link cropCenter} for a complex `{re, im}` array. */
export function cropCenterComplex({re, im}, nx, ny, cropX, cropY = cropX) {
  const r = cropCenter(re, nx, ny, cropX, cropY);
  const i = cropCenter(im, nx, ny, cropX, cropY);
  return {re: r.data, im: i.data, nx: r.nx, ny: r.ny};
}

/**
 * Central-difference gradient of a 2D field, matching
 * `numpy.gradient(field, dx, dy, edge_order=2)`.
 *
 * Interior points use the second-order central difference; boundaries use the
 * second-order one-sided formula, which is what `edge_order=2` selects and what
 * `aberration_utils.py` relies on.
 *
 * @returns {{dRow: Float64Array, dCol: Float64Array}} derivatives along axis 0 and 1
 */
export function gradient2d(field, nx, ny, dx = 1, dy = 1) {
  const dRow = new Float64Array(nx * ny);
  const dCol = new Float64Array(nx * ny);
  const at = (ix, iy) => field[ix * ny + iy];

  for (let iy = 0; iy < ny; iy++) {
    for (let ix = 0; ix < nx; ix++) {
      let v;
      if (nx < 2) v = 0;
      else if (ix === 0) v = (-3 * at(0, iy) + 4 * at(1, iy) - (nx > 2 ? at(2, iy) : at(1, iy))) / (2 * dx);
      else if (ix === nx - 1) {
        v = (3 * at(nx - 1, iy) - 4 * at(nx - 2, iy) + (nx > 2 ? at(nx - 3, iy) : at(nx - 2, iy))) / (2 * dx);
      } else v = (at(ix + 1, iy) - at(ix - 1, iy)) / (2 * dx);
      dRow[ix * ny + iy] = v;
    }
  }

  for (let ix = 0; ix < nx; ix++) {
    for (let iy = 0; iy < ny; iy++) {
      let v;
      if (ny < 2) v = 0;
      else if (iy === 0) v = (-3 * at(ix, 0) + 4 * at(ix, 1) - (ny > 2 ? at(ix, 2) : at(ix, 1))) / (2 * dy);
      else if (iy === ny - 1) {
        v = (3 * at(ix, ny - 1) - 4 * at(ix, ny - 2) + (ny > 2 ? at(ix, ny - 3) : at(ix, ny - 2))) / (2 * dy);
      } else v = (at(ix, iy + 1) - at(ix, iy - 1)) / (2 * dy);
      dCol[ix * ny + iy] = v;
    }
  }

  return {dRow, dCol};
}

/**
 * Nearest-neighbour pull-warp of an image by a displacement field, with wrapping.
 *
 * Ported from `aberration_utils.py:compute_warp`, which uses
 * `scipy.ndimage.map_coordinates(..., order=0, mode="wrap")`.
 *
 * @returns {{warped: Float64Array, inside: Float64Array}} the warped image and a
 *   mask that is 1 where the source sample came from within the original frame
 */
export function warpNearest(image, nx, ny, uRow, uCol) {
  const warped = new Float64Array(nx * ny);
  const inside = new Float64Array(nx * ny);
  for (let ix = 0; ix < nx; ix++) {
    for (let iy = 0; iy < ny; iy++) {
      const i = ix * ny + iy;
      const sr = Math.round(ix - uRow[i]);
      const sc = Math.round(iy - uCol[i]);
      inside[i] = sr >= 0 && sr < nx && sc >= 0 && sc < ny ? 1 : 0;
      const wr = ((sr % nx) + nx) % nx;
      const wc = ((sc % ny) + ny) % ny;
      warped[i] = image[wr * ny + wc];
    }
  }
  return {warped, inside};
}
