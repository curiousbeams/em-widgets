// A 4D-STEM forward model, sized and shaped for a browser widget.
//
// Simulating the experiment rather than downloading it is what makes a dose
// slider mean anything, and dose is what separates the direct estimators from
// each other. The model is single-slice — one transmission function,
// `O(r) = exp(i phi(r))` — so it shows the weak phase approximation breaking
// down as the object strengthens, but not multiple scattering. Widgets built on
// it should say so.
//
// Two decisions keep it affordable. The probe moves by an integer index shift
// rather than a Fourier phase ramp, which is exact because the scan divides the
// object grid and saves a transform per position. And only the bright-field
// pixels are retained: the full `I(R,k)` at 128 x 128 over a 32 x 32 scan is
// 134 MB, where the bright-field subset is 3.3 MB and is all the direct methods
// ever read.

import {complex, fft2} from "./fft.js";
import {poisson} from "./image.js";

/**
 * One scan position: `exit(r) = exp(i phi(r)) psi(r - R)`, then `|FFT2|^2`.
 *
 * @param {Float64Array} objectPhase `n * n`, radians
 * @param {{re, im}} probeReal real-space probe, centred at index (0,0)
 * @param {number} n grid size (square)
 * @param {number} row probe position, in pixels
 * @param {number} col
 * @param {{re, im}} scratch an `n * n` complex workspace, reused across calls
 * @param {Float64Array} out `n * n`, receives the diffraction intensity
 */
export function measurePosition(objectPhase, probeReal, n, row, col, scratch, out) {
  for (let ix = 0; ix < n; ix++) {
    // Reading the probe at a rolled index is the shift. Exact for whole pixels.
    const sx = (((ix - row) % n) + n) % n;
    const dst = ix * n;
    const src = sx * n;
    for (let iy = 0; iy < n; iy++) {
      const sy = (((iy - col) % n) + n) % n;
      const p = src + sy;
      const i = dst + iy;
      const phase = objectPhase[i];
      const cos = Math.cos(phase);
      const sin = Math.sin(phase);
      const pr = probeReal.re[p];
      const pi = probeReal.im[p];
      scratch.re[i] = cos * pr - sin * pi;
      scratch.im[i] = cos * pi + sin * pr;
    }
  }
  fft2(scratch, n, n);
  for (let i = 0; i < n * n; i++) {
    out[i] = scratch.re[i] * scratch.re[i] + scratch.im[i] * scratch.im[i];
  }
  return out;
}

/**
 * A resumable scan. `step(count)` measures the next `count` positions and
 * returns how many remain, so a notebook can spend a frame's worth of time on
 * it and yield.
 *
 * Shot noise is applied only to the retained pixels, with the electron count
 * taken from the whole pattern's total — noising the full pattern and then
 * discarding 95% of it costs twenty times as much for the same answer.
 *
 * @param {object} options
 * @param {Float64Array} options.objectPhase `n * n`, radians
 * @param {{re, im}} options.probeReal real-space probe, centred at (0,0)
 * @param {number} options.n object and detector grid size
 * @param {number} options.m scan grid size; must divide `n`
 * @param {Int32Array|null} [options.bfIndices] detector pixels to keep, or null
 *   for all of them (which is what ePIE needs)
 * @param {number} [options.dose=Infinity] electrons per square Angstrom
 * @param {[number, number]} [options.sampling=[1,1]] object sampling [Angstrom]
 * @param {boolean} [options.amplitude=false] store `sqrt(I)` rather than `I`
 * @param {() => number} [options.random=Math.random]
 * @returns {{step(count: number): number, data: Float32Array, total: number,
 *   done: number, nDetector: number, m: number}}
 *   `data` is laid out `[k][R]` — k-major, so each virtual image is contiguous
 */
export function forwardModel({
  objectPhase,
  probeReal,
  n,
  m,
  bfIndices = null,
  dose = Infinity,
  sampling = [1, 1],
  amplitude = false,
  random = Math.random
}) {
  if (n % m !== 0) {
    throw new Error(`scan grid m=${m} must divide the object grid n=${n}`);
  }
  const stride = n / m;
  const indices = bfIndices ?? Int32Array.from({length: n * n}, (_, i) => i);
  const nDetector = indices.length;
  const total = m * m;

  // Electrons per pattern, from the dose and the area one scan position covers.
  const stepArea = stride * sampling[0] * stride * sampling[1];
  const electrons = Number.isFinite(dose) ? dose * stepArea : Infinity;

  const data = new Float32Array(nDetector * total);
  const scratch = complex(n * n);
  const pattern = new Float64Array(n * n);
  let done = 0;

  const step = (count) => {
    const stop = Math.min(total, done + count);
    for (; done < stop; done++) {
      const scanRow = (done / m) | 0;
      const scanCol = done % m;
      measurePosition(objectPhase, probeReal, n, scanRow * stride, scanCol * stride, scratch, pattern);

      let patternTotal = 0;
      if (Number.isFinite(electrons)) {
        for (let i = 0; i < pattern.length; i++) patternTotal += pattern[i];
      }
      const scale = patternTotal > 0 ? electrons / patternTotal : 0;

      const scanIndex = scanRow * m + scanCol;
      for (let j = 0; j < nDetector; j++) {
        let value = pattern[indices[j]];
        if (Number.isFinite(electrons)) {
          value = poisson(value * scale, random) / scale;
        }
        data[j * total + scanIndex] = amplitude ? Math.sqrt(Math.max(0, value)) : value;
      }
    }
    return total - done;
  };

  return {step, data, total, done: 0, nDetector, m, get progress() { return done / total; }};
}

/**
 * Fourier transform every virtual bright-field image over the scan axes:
 * `G(q,k) = F_{R -> q}[ I(R,k) ]`.
 *
 * **The zero-frequency term is discarded.** In the forward model it is the
 * `|psi(k)|^2 delta(q)` term — the unscattered beam — which carries nothing
 * about the specimen and is larger than everything that does.
 *
 * Keeping it is harmless until the spectrum is tiled, and then it is not. Tiling
 * replicates the whole spectrum, so that one enormous spike reappears at `f^2`
 * positions across the reconstruction grid, all but one of them at a high
 * spatial frequency where it is pure fiction. The overlap kernels mostly
 * suppress the copies, because `Gamma` is small or zero there; parallax and iCOM
 * have constant magnitude and pass every one, which puts a periodic lattice of
 * spikes through the reconstruction. Zeroing it here removes the problem at
 * source for every method.
 *
 * @param {Float32Array} data from {@link forwardModel}, laid out `[k][R]`
 * @param {number} m scan grid size
 * @param {number} nDetector how many detector pixels were retained
 * @returns {{re: Float64Array, im: Float64Array}} laid out `[k][q]`
 */
export function scanSpectra(data, m, nDetector) {
  const area = m * m;
  const out = complex(nDetector * area);
  const scratch = complex(area);

  for (let j = 0; j < nDetector; j++) {
    const base = j * area;
    for (let i = 0; i < area; i++) {
      scratch.re[i] = data[base + i];
      scratch.im[i] = 0;
    }
    fft2(scratch, m, m);
    scratch.re[0] = 0;
    scratch.im[0] = 0;
    out.re.set(scratch.re, base);
    out.im.set(scratch.im, base);
  }
  return out;
}

/** One detector pixel's spectrum, as a view-shaped copy for {@link tileSpectrum}. */
export function spectrumAt(spectra, m, index, out = complex(m * m)) {
  const area = m * m;
  const base = index * area;
  out.re.set(spectra.re.subarray(base, base + area));
  out.im.set(spectra.im.subarray(base, base + area));
  return out;
}
