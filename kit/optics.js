// Electron-optical wave functions: the aberration function, apertures, probes,
// and Fresnel/multislice propagation.
//
// Ported from `aberration_utils.py` (chi, apertures), `ctf/utils.py`,
// `py4DSTEM.process.phase.utils.ComplexProbe.evaluate_chi`, and the
// `return_propagator_array` / `propagate_wavefunction` / `multislice_propagation`
// helpers in `02.stem-measurements-nb.ipynb`.
//
// Angles are in radians unless a name says mrad. Lengths are Angstroms.
// Grids are corner-centered (fftfreq order) unless stated otherwise.

import {electronWavelength} from "./units.js";
import {fft2, ifft2, complex, expi, multiply} from "./fft.js";
import {gradient2d} from "./image.js";
import {angularSpatialFrequencies, angularSampling, spatialFrequencies} from "./grid.js";

/**
 * The wave aberration function chi(alpha, phi), to sixth order.
 *
 * Coefficients follow the Krivanek notation used throughout the lab's Python:
 * `Cnm` in Angstroms and `phinm` in radians. Only the orders you supply are
 * evaluated, exactly as in `aberration_utils.py:aberration_surface`.
 *
 * @param {Float64Array} alpha scattering angle [rad] (= k * wavelength)
 * @param {Float64Array} phi azimuth [rad]
 * @param {number} wavelength [Angstrom]
 * @param {Record<string, number>} coefficients e.g. `{C10: -50, C30: 1e4, C12: 12, phi12: 0.3}`
 * @returns {Float64Array} chi [rad]
 */
export function chi(alpha, phi, wavelength, coefficients = {}) {
  const c = (name) => coefficients[name] ?? 0;
  const has = (...names) => names.some((n) => n in coefficients);
  const prefactor = (2 * Math.PI) / wavelength;
  const out = new Float64Array(alpha.length);

  const order2 = has("C10", "C12", "phi12");
  const order3 = has("C21", "phi21", "C23", "phi23");
  const order4 = has("C30", "C32", "phi32", "C34", "phi34");
  const order5 = has("C41", "phi41", "C43", "phi43", "C45", "phi45");
  const order6 = has("C50", "C52", "phi52", "C54", "phi54", "C56", "phi56");

  for (let i = 0; i < alpha.length; i++) {
    const a = alpha[i];
    const p = phi[i];
    const a2 = a * a;
    let v = 0;

    if (order2) v += 0.5 * a2 * (c("C10") + c("C12") * Math.cos(2 * (p - c("phi12"))));
    if (order3) {
      v += (1 / 3) * a2 * a *
        (c("C21") * Math.cos(p - c("phi21")) + c("C23") * Math.cos(3 * (p - c("phi23"))));
    }
    if (order4) {
      v += (1 / 4) * a2 * a2 *
        (c("C30") +
          c("C32") * Math.cos(2 * (p - c("phi32"))) +
          c("C34") * Math.cos(4 * (p - c("phi34"))));
    }
    if (order5) {
      v += (1 / 5) * a2 * a2 * a *
        (c("C41") * Math.cos(p - c("phi41")) +
          c("C43") * Math.cos(3 * (p - c("phi43"))) +
          c("C45") * Math.cos(5 * (p - c("phi45"))));
    }
    if (order6) {
      v += (1 / 6) * a2 * a2 * a2 *
        (c("C50") +
          c("C52") * Math.cos(2 * (p - c("phi52"))) +
          c("C54") * Math.cos(4 * (p - c("phi54"))) +
          c("C56") * Math.cos(6 * (p - c("phi56"))));
    }
    out[i] = v * prefactor;
  }
  return out;
}

/**
 * Geometric ray displacement from the aberration function, in **pixels**.
 *
 * The displacement of a ray is the gradient of chi with respect to spatial
 * frequency; dividing by the reciprocal sampling again converts it to pixels,
 * which is what a nearest-neighbour warp wants.
 *
 * Expects a **centred** (fftshifted) chi, matching `aberration_utils.py`.
 *
 * @param {Float64Array} chiCentered
 * @param {[number, number]} reciprocalSampling `1 / (gpts * sampling)` per axis
 * @returns {{uRow: Float64Array, uCol: Float64Array}} displacement in pixels
 */
export function rayDisplacement(chiCentered, nx, ny, [dkx, dky]) {
  const {dRow, dCol} = gradient2d(chiCentered, nx, ny, dkx, dky);
  const uRow = new Float64Array(nx * ny);
  const uCol = new Float64Array(nx * ny);
  for (let i = 0; i < uRow.length; i++) {
    uRow[i] = dRow[i] / dkx;
    uCol[i] = dCol[i] / dky;
  }
  return {uRow, uCol};
}

/** Scherzer defocus for a given spherical aberration: `-sign(Cs) sqrt(1.5 |Cs| lambda)`. */
export function scherzerDefocus(Cs, wavelength) {
  return -Math.sign(Cs) * Math.sqrt(1.5 * Math.abs(Cs) * wavelength);
}

/**
 * Hard-edged aperture: 1 inside the cutoff, 0 outside.
 * @param {number} semiangleCutoff [mrad]
 */
export function hardAperture(alpha, semiangleCutoff) {
  const cutoff = semiangleCutoff * 1e-3;
  const out = new Float64Array(alpha.length);
  for (let i = 0; i < alpha.length; i++) out[i] = alpha[i] < cutoff ? 1 : 0;
  return out;
}

/**
 * Antialiased aperture whose edge is one pixel wide.
 *
 * Ported from `aberration_utils.py:soft_aperture_function`.
 *
 * @param {number} semiangleCutoff [mrad]
 * @param {[number, number]} angularSampling [mrad] per pixel — see grid.angularSampling
 */
export function softAperture(alpha, phi, semiangleCutoff, [dax, day]) {
  const cutoff = semiangleCutoff * 1e-3;
  const out = new Float64Array(alpha.length);
  for (let i = 0; i < alpha.length; i++) {
    const denominator = Math.hypot(Math.cos(phi[i]) * dax * 1e-3, Math.sin(phi[i]) * day * 1e-3);
    const v = (cutoff - alpha[i]) / denominator + 0.5;
    out[i] = v < 0 ? 0 : v > 1 ? 1 : v;
  }
  return out;
}

/** {@link softAperture} or {@link hardAperture}, chosen by `softEdges`. */
export function aperture(alpha, phi, semiangleCutoff, angularSampling, softEdges = true) {
  return softEdges
    ? softAperture(alpha, phi, semiangleCutoff, angularSampling)
    : hardAperture(alpha, semiangleCutoff);
}

/**
 * A converged STEM probe in Fourier space: `A(k) exp(-i chi(k))`, L2-normalised.
 *
 * Equivalent to `py4DSTEM.process.phase.utils.ComplexProbe(...).build()`.
 *
 * @param {object} options
 * @param {[number, number]} options.gpts
 * @param {[number, number]} options.sampling [Angstrom]
 * @param {number} options.energy [eV]
 * @param {number} options.semiangleCutoff [mrad]
 * @param {Record<string, number>} [options.aberrations] e.g. `{C10: -80, C30: 1e4}`
 * @param {boolean} [options.softEdges=true]
 * @returns {{re: Float64Array, im: Float64Array}} corner-centered Fourier-space probe
 */
export function complexProbe({
  gpts,
  sampling,
  energy,
  semiangleCutoff,
  aberrations = {},
  softEdges = true
}) {
  const wavelength = electronWavelength(energy);
  const [nx, ny] = gpts;
  const {alpha, phi} = angularSpatialFrequencies(gpts, sampling, wavelength);
  const A = aperture(alpha, phi, semiangleCutoff, angularSampling(gpts, sampling, wavelength), softEdges);

  const chiArray = chi(alpha, phi, wavelength, aberrations);
  const out = complex(nx * ny);
  let norm = 0;
  for (let i = 0; i < out.re.length; i++) {
    const c = -chiArray[i];
    out.re[i] = A[i] * Math.cos(c);
    out.im[i] = A[i] * Math.sin(c);
    norm += out.re[i] * out.re[i] + out.im[i] * out.im[i];
  }
  norm = Math.sqrt(norm) || 1;
  for (let i = 0; i < out.re.length; i++) {
    out.re[i] /= norm;
    out.im[i] /= norm;
  }
  return out;
}

/**
 * Fresnel free-space propagator for one slice: `exp(i pi lambda dz k^2)`.
 *
 * Sign convention taken from `return_propagator_array` in
 * `02.stem-measurements-nb.ipynb`. Multiply into a wavefunction's Fourier
 * transform, then inverse transform.
 *
 * @param {[number, number]} gpts
 * @param {[number, number]} sampling [Angstrom]
 * @param {number} energy [eV]
 * @param {number} dz slice thickness [Angstrom]
 */
export function fresnelPropagator(gpts, sampling, energy, dz) {
  const wavelength = electronWavelength(energy);
  const {kx, ky} = spatialFrequencies(gpts, sampling);
  const prefactor = wavelength * Math.PI * dz;
  const phase = new Float64Array(kx.length);
  for (let i = 0; i < kx.length; i++) {
    phase[i] = (kx[i] * kx[i] + ky[i] * ky[i]) * prefactor;
  }
  return expi(phase);
}

/** Propagate a wavefunction one slice through free space. */
export function propagate(wave, propagator, nx, ny) {
  const F = fft2({re: Float64Array.from(wave.re), im: Float64Array.from(wave.im)}, nx, ny);
  multiply(F, propagator, F);
  return ifft2(F, nx, ny);
}

/**
 * Multislice propagation of a wavefunction through a stack of transmission
 * functions.
 *
 * Ported from `multislice_propagation()` in `02.stem-measurements-nb.ipynb`:
 * transmit through each slice, propagate between slices (but not after the last).
 *
 * @param {{re: Float64Array, im: Float64Array}} wave real-space entrance wave
 * @param {Array<{re: Float64Array, im: Float64Array}>} transmission per-slice `exp(i sigma V)`
 * @param {{re: Float64Array, im: Float64Array}} propagator from {@link fresnelPropagator}
 * @returns {{re: Float64Array, im: Float64Array}} exit wave
 */
export function multislice(wave, transmission, propagator, nx, ny) {
  // Copy once, then work in place: the propagation step is the hot loop, and a
  // 244 x 242 grid churns ~1 MB per copy.
  const psi = {re: Float64Array.from(wave.re), im: Float64Array.from(wave.im)};
  for (let s = 0; s < transmission.length; s++) {
    multiply(psi, transmission[s], psi);
    if (s + 1 < transmission.length) {
      fft2(psi, nx, ny);
      multiply(psi, propagator, psi);
      ifft2(psi, nx, ny);
    }
  }
  return psi;
}

/**
 * Transmission function `exp(i sigma V)` for a projected potential slice.
 *
 * @param {ArrayLike<number>} potential projected potential of the slice
 * @param {number} [sigma=1] interaction parameter; the lab's cached potentials
 *   are already in radians, so the default leaves them untouched
 */
export function transmissionFunction(potential, sigma = 1) {
  const phase = new Float64Array(potential.length);
  for (let i = 0; i < potential.length; i++) phase[i] = potential[i] * sigma;
  return expi(phase);
}
