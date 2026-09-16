import test from "node:test";
import assert from "node:assert/strict";

import {electronWavelength} from "../kit/units.js";
import {angularSpatialFrequencies, angularSampling} from "../kit/grid.js";
import {complexProbe, chi, softAperture, hardAperture} from "../kit/optics.js";
import {whiteNoiseObject2D} from "../kit/image.js";
import {mulberry32} from "../kit/specimen.js";
import {complex, fft2, ifft2} from "../kit/fft.js";
import {
  overlapFunction, overlapSums, brightFieldIndices, tileSpectrum,
  checkOverlapGeometry, directAccumulator, accumulatePixel, directImage,
  overlapFunctionAtQ, overlapRegions, parallaxShifts,
  directCTF, directSSNR, PTYCHO_METHODS
} from "../kit/ptycho.js";
import {forwardModel, scanSpectra, spectrumAt} from "../kit/ptycho-sim.js";

// A small grid throughout: every identity here is exact or has a stated
// truncation order, so size buys nothing but seconds.
const N = 64;
const SAMPLING = [0.25, 0.25];
const ENERGY = 300e3;
const WAVELENGTH = electronWavelength(ENERGY);
const DK = 1 / (N * SAMPLING[0]);

/** A semiangle giving exactly `pixels` aperture radii on the N grid. */
const semiangleFor = (pixels) => pixels * DK * WAVELENGTH * 1e3;

const RADIUS = 8;
const SEMIANGLE = semiangleFor(RADIUS);

const geometry = angularSpatialFrequencies([N, N], SAMPLING, WAVELENGTH);
const SAMPLING_MRAD = angularSampling([N, N], SAMPLING, WAVELENGTH);

const maxAbs = (a) => a.reduce((m, v) => Math.max(m, Math.abs(v)), 0);

/** Signed frequency index of array position `i` on a length-`n` fftfreq axis. */
const signed = (i, n) => (i < n / 2 ? i : i - n);

function l2Normalised(array) {
  let norm = 0;
  for (const v of array) norm += v * v;
  norm = Math.sqrt(norm) || 1;
  return Float64Array.from(array, (v) => v / norm);
}

// ---------------------------------------------------------------------------
// The aperture overlap function
// ---------------------------------------------------------------------------

test("axial illumination reproduces the HRTEM transfer function", (t) => {
  // Gamma(q, k = 0) = -2i psi(0) A(q) sin[chi(q)] — the reciprocity link between
  // a STEM detector pixel on axis and a TEM contrast transfer function. It pins
  // the sign of chi, the probe conjugation, and the k-q versus k+q index order
  // all at once, and it needs no data and no transform.
  //
  // Even aberrations only: chi must be even in q for psi(-q) = psi(q), and C21
  // is odd.
  const aberrations = {C10: -60, C30: 1e4, C12: 12};
  const probe = complexProbe({
    gpts: [N, N], sampling: SAMPLING, energy: ENERGY, semiangleCutoff: SEMIANGLE, aberrations
  });
  const chiArray = chi(geometry.alpha, geometry.phi, WAVELENGTH, aberrations);
  const aperture = l2Normalised(softAperture(geometry.alpha, geometry.phi, SEMIANGLE, SAMPLING_MRAD));

  const gamma = overlapFunction(probe, N, N, 0, complex(N * N));
  const psi0 = probe.re[0];

  let worst = 0;
  for (let i = 0; i < N * N; i++) {
    worst = Math.max(
      worst,
      Math.abs(gamma.re[i]),
      Math.abs(gamma.im[i] - -2 * psi0 * aperture[i] * Math.sin(chiArray[i]))
    );
  }
  t.diagnostic(`max error ${worst.toExponential(2)}`);
  assert.ok(worst < 1e-15, `axial Gamma is off by ${worst}`);
});

test("the overlap function is antisymmetric and vanishes at zero frequency", () => {
  const probe = complexProbe({
    gpts: [N, N], sampling: SAMPLING, energy: ENERGY, semiangleCutoff: SEMIANGLE,
    aberrations: {C10: -60, C21: 5e3, phi21: 0.4}
  });

  // Every k, not a sample: these are the tests that catch an fftfreq-ordering
  // slip, and a wraparound error need not show up at the k you happened to pick.
  for (const kIndex of [0, 1, N + 3, 5 * N + 11, N * N - 1]) {
    const gamma = overlapFunction(probe, N, N, kIndex, complex(N * N));
    assert.ok(Math.hypot(gamma.re[0], gamma.im[0]) < 1e-18,
      `Gamma(0, k=${kIndex}) does not vanish`);

    let worst = 0;
    for (let a = 0; a < N; a++) {
      for (let b = 0; b < N; b++) {
        const i = a * N + b;
        const mirrored = ((N - a) % N) * N + ((N - b) % N);
        worst = Math.max(
          worst,
          Math.abs(gamma.re[mirrored] + gamma.re[i]),
          Math.abs(gamma.im[mirrored] - gamma.im[i])
        );
      }
    }
    assert.ok(worst < 1e-18, `Gamma(-q,k) != -conj Gamma(q,k) at k=${kIndex}, off by ${worst}`);
  }
});

test("the overlap function is the same on a tiled sub-block as on the full grid", () => {
  // The reconstruction evaluates Gamma on a q x q sub-block of the probe's own
  // grid, which is where the index wraparound is easiest to get wrong and where
  // nothing else looks. Every sub-block entry must equal the corresponding
  // full-grid entry.
  const probe = complexProbe({
    gpts: [N, N], sampling: SAMPLING, energy: ENERGY, semiangleCutoff: SEMIANGLE,
    aberrations: {C10: -60}
  });
  const full = overlapFunction(probe, N, N, 3 * N + 2, complex(N * N));

  for (const q of [16, 32]) {
    const block = overlapFunction(probe, N, q, 3 * N + 2, complex(q * q));
    let worst = 0;
    for (let a = 0; a < q; a++) {
      for (let b = 0; b < q; b++) {
        const fa = signed(a, q);
        const fb = signed(b, q);
        const i = (((fa % N) + N) % N) * N + (((fb % N) + N) % N);
        const o = a * q + b;
        worst = Math.max(worst, Math.abs(block.re[o] - full.re[i]), Math.abs(block.im[o] - full.im[i]));
      }
    }
    assert.ok(worst < 1e-18, `q=${q} sub-block disagrees with the full grid by ${worst}`);
  }
});

test("the two slices through the overlap function agree where they cross", () => {
  // overlapFunction fixes k and varies q; overlapFunctionAtQ fixes q and varies
  // k. They are slices through the same function of two variables, so at every
  // (q, k) they must return the same number. Two independent index derivations,
  // checked against each other.
  const probe = complexProbe({
    gpts: [N, N], sampling: SAMPLING, energy: ENERGY, semiangleCutoff: SEMIANGLE,
    aberrations: {C10: -60, C12: 20, phi12: 0.7}
  });

  let worst = 0;
  for (const kIndex of [0, N + 3, 7 * N + 2, N * N - 5]) {
    const overQ = overlapFunction(probe, N, N, kIndex, complex(N * N));
    for (const qIndex of [0, 1, 2 * N + 5, 9 * N + 9, N * N - 1]) {
      const overK = overlapFunctionAtQ(probe, N, qIndex, complex(N * N));
      worst = Math.max(
        worst,
        Math.abs(overQ.re[qIndex] - overK.re[kIndex]),
        Math.abs(overQ.im[qIndex] - overK.im[kIndex])
      );
    }
  }
  assert.ok(worst < 1e-18, `the two slices disagree by ${worst}`);
});

test("the overlap regions are where the shifted apertures actually meet", () => {
  const aperture = hardAperture(geometry.alpha, SEMIANGLE);

  // Beyond twice the aperture radius nothing overlaps at all, so both regions
  // must be empty — this is the transfer limit of every direct method.
  const beyond = Math.round(2 * RADIUS + 4);
  const far = overlapRegions(aperture, N, beyond * N);
  assert.ok(far.double.every((v) => v === 0), "something overlaps past twice the aperture radius");
  assert.ok(far.triple.every((v) => v === 0), "a triple overlap past twice the aperture radius");

  // Between one and two radii the two lens regions exist and never meet.
  const single = Math.round(1.4 * RADIUS);
  const annulus = overlapRegions(aperture, N, single * N);
  const doubleArea = annulus.double.reduce((s, v) => s + (v > 0 ? 1 : 0), 0);
  assert.ok(doubleArea > 0, "no double overlap inside the single-overlap annulus");
  assert.ok(annulus.triple.every((v) => v === 0), "a triple overlap where there should be none");
  // Both lenses, and by symmetry the same size.
  const first = annulus.double.reduce((s, v) => s + (v === 1 ? 1 : 0), 0);
  const second = annulus.double.reduce((s, v) => s + (v === 2 ? 1 : 0), 0);
  assert.equal(first, second, "the two double-overlap regions came out different sizes");

  // Below the aperture radius all three discs meet.
  const low = overlapRegions(aperture, N, Math.round(0.4 * RADIUS) * N);
  assert.ok(low.triple.reduce((s, v) => s + v, 0) > 0, "no triple overlap at low frequency");
});

test("an aperture too large for the grid is refused rather than silently aliased", () => {
  // The wraparound in overlapFunction is only harmless while 4 * radius < n.
  assert.doesNotThrow(() => checkOverlapGeometry(64, 64, 15));
  assert.throws(() => checkOverlapGeometry(64, 64, 16), /4 \* radius < n/);
  assert.throws(() => checkOverlapGeometry(64, 128, 8), /exceeds probe grid/);
  assert.throws(() => checkOverlapGeometry(64, 24, 8), /must divide/);
});

test("in focus, the overlap sum is the aperture autocorrelation — in the annulus", (t) => {
  // With chi = 0 and a HARD aperture, A(k-q) and A(k+q) are indicators, so
  // |Gamma| is their exclusive-or and
  //
  //     sum_k |Gamma(q,k)| = 2 [ (A*A)(q) - T(q) ] / |A|
  //
  // with T the triple overlap. Only where the triple overlap vanishes — the
  // annulus R < |q| < 2R — does this reduce to the plain autocorrelation. The
  // bare claim "in-focus CTF is the aperture autocorrelation" is wrong below R,
  // by order unity, which is why this test states its domain.
  //
  // A soft aperture breaks the exclusive-or and the identity degrades to ~4e-2.
  const probe = complexProbe({
    gpts: [N, N], sampling: SAMPLING, energy: ENERGY,
    semiangleCutoff: SEMIANGLE, aberrations: {}, softEdges: false
  });
  const bf = brightFieldIndices(geometry.alpha, SEMIANGLE);
  const sums = overlapSums(probe, N, N, bf);

  const aperture = l2Normalised(hardAperture(geometry.alpha, SEMIANGLE));
  const f = {re: Float64Array.from(aperture), im: new Float64Array(N * N)};
  fft2(f, N, N);
  for (let i = 0; i < N * N; i++) {
    f.re[i] = f.re[i] * f.re[i] + f.im[i] * f.im[i];
    f.im[i] = 0;
  }
  ifft2(f, N, N);

  let annulus = 0;
  let everywhere = 0;
  for (let a = 0; a < N; a++) {
    for (let b = 0; b < N; b++) {
      const i = a * N + b;
      const q = Math.hypot(signed(a, N), signed(b, N));
      const error = Math.abs(sums.sum[i] - 2 * f.re[i]);
      everywhere = Math.max(everywhere, error);
      if (q > RADIUS + 1.5 && q < 2 * RADIUS - 1.5) annulus = Math.max(annulus, error);
    }
  }
  t.diagnostic(`annulus ${annulus.toExponential(2)}, whole grid ${everywhere.toExponential(2)}`);
  assert.ok(annulus < 1e-13, `single-overlap annulus is off by ${annulus}`);
  // And the triple-overlap correction really is order unity at low q, so a test
  // that asserted the bare identity everywhere would be asserting something false.
  assert.ok(everywhere > 1, "the triple-overlap correction has gone missing");
});

test("the parallax shift reduces to lambda C10 k for pure defocus", (t) => {
  // dr = grad chi / (2 pi), which for chi = pi lambda C10 k^2 is exactly
  // lambda C10 k. Getting the 2 pi wrong, or the reciprocal sampling, puts
  // every virtual bright-field image in the wrong place and parallax quietly
  // reconstructs blur.
  const C10 = -400;
  const aberrations = {C10};
  const chiArray = chi(geometry.alpha, geometry.phi, WAVELENGTH, aberrations);
  const qSampling = SAMPLING[0];

  const {shiftRow, shiftCol} = parallaxShifts(chiArray, N, DK, qSampling);

  let worst = 0;
  let scale = 0;
  for (let a = 0; a < N; a++) {
    for (let b = 0; b < N; b++) {
      const i = a * N + b;
      // Inside the bright field, which is the only place parallax uses it.
      if (Math.hypot(signed(a, N), signed(b, N)) > RADIUS) continue;
      const expectedRow = (WAVELENGTH * C10 * signed(a, N) * DK) / qSampling;
      const expectedCol = (WAVELENGTH * C10 * signed(b, N) * DK) / qSampling;
      worst = Math.max(worst, Math.abs(shiftRow[i] - expectedRow), Math.abs(shiftCol[i] - expectedCol));
      scale = Math.max(scale, Math.abs(expectedRow), Math.abs(expectedCol));
    }
  }
  t.diagnostic(`max error ${worst.toExponential(2)} px against a largest shift of ${scale.toFixed(2)} px`);
  assert.ok(worst < 1e-10, `parallax shifts are off by ${worst} pixels`);

  // And zero defocus must give no shift at all.
  const flat = parallaxShifts(chi(geometry.alpha, geometry.phi, WAVELENGTH, {}), N, DK, qSampling);
  assert.ok(maxAbs(flat.shiftRow) < 1e-12, "an unaberrated probe produced a parallax shift");
});

// ---------------------------------------------------------------------------
// Tiling
// ---------------------------------------------------------------------------

test("Fourier tiling replicates the spectrum exactly", () => {
  const m = 32;
  const g = complex(m * m);
  for (let i = 0; i < m * m; i++) {
    g.re[i] = Math.sin(i * 1.7);
    g.im[i] = Math.cos(i * 0.3);
  }

  for (const f of [1, 2, 4]) {
    const q = m * f;
    const tiled = tileSpectrum(g, m, f);

    assert.equal(tiled.re[0], g.re[0], "tiling moved the DC term");
    assert.equal(tiled.im[0], g.im[0], "tiling moved the DC term");

    let source = 0;
    let target = 0;
    for (let i = 0; i < m * m; i++) source += g.re[i] ** 2 + g.im[i] ** 2;
    for (let i = 0; i < q * q; i++) target += tiled.re[i] ** 2 + tiled.im[i] ** 2;
    assert.ok(Math.abs(target / source - f * f) < 1e-12,
      `f=${f} tiled energy ratio ${target / source}, want ${f * f}`);

    // Every entry is the source at the frequency taken modulo m.
    for (let a = 0; a < q; a++) {
      for (let b = 0; b < q; b++) {
        assert.equal(tiled.re[a * q + b], g.re[(a % m) * m + (b % m)]);
      }
    }
  }
});

// ---------------------------------------------------------------------------
// The forward model
// ---------------------------------------------------------------------------

test("the simulated dataset obeys the weak-phase forward model", (t) => {
  // G(q,k) = i Gamma(q,k) phi(q), the identity that ties the simulation to the
  // inversion. Nothing else here can catch a transposed scan axis or a
  // probe-shift sign.
  //
  // Measured: the fitted complex scale is i to five decimals at every object
  // strength, and the residual falls linearly with the object phase — 3.4e-3 at
  // max|phi| = 1e-2, 3.7e-4 at 1e-3 — which is the O(phi^2) truncation and
  // nothing else. Below about 1e-3 it stops improving and settles near 1e-3,
  // where float64 round-off on an intensity of order one takes over.
  const n = 32;
  const m = 32;
  const sampling = [0.25, 0.25];
  const wavelength = electronWavelength(ENERGY);
  const grid = angularSpatialFrequencies([n, n], sampling, wavelength);
  const semiangle = 4 * (1 / (n * sampling[0])) * wavelength * 1e3;

  const probe = complexProbe({
    gpts: [n, n], sampling, energy: ENERGY, semiangleCutoff: semiangle,
    aberrations: {C10: -40, C30: 1e4}
  });
  const probeReal = ifft2({re: Float64Array.from(probe.re), im: Float64Array.from(probe.im)}, n, n);

  const results = [];
  for (const amplitude of [1e-2, 1e-3]) {
    const phase = new Float64Array(n * n);
    let peak = 0;
    for (let a = 0; a < n; a++) {
      for (let b = 0; b < n; b++) {
        const v = Math.sin(0.7 * a + 1.3) * Math.cos(0.4 * b - 0.2) + 0.5 * Math.sin(0.21 * a * b);
        phase[a * n + b] = v;
        peak = Math.max(peak, Math.abs(v));
      }
    }
    for (let i = 0; i < phase.length; i++) phase[i] *= amplitude / peak;

    const phiHat = {re: Float64Array.from(phase), im: new Float64Array(n * n)};
    fft2(phiHat, n, n);

    const bf = brightFieldIndices(grid.alpha, semiangle);
    const model = forwardModel({objectPhase: phase, probeReal, n, m, bfIndices: bf, dose: Infinity});
    assert.equal(model.step(model.total), 0, "the scan did not finish");
    const G = scanSpectra(model.data, m, bf.length);

    const gamma = complex(m * m);
    let worstResidual = 0;
    let worstScale = 0;
    for (const j of [0, 3, 11, bf.length - 1]) {
      overlapFunction(probe, n, m, bf[j], gamma);
      const base = j * m * m;
      let numRe = 0, numIm = 0, denom = 0, total = 0;
      for (let i = 1; i < m * m; i++) {
        const mr = gamma.re[i] * phiHat.re[i] - gamma.im[i] * phiHat.im[i];
        const mi = gamma.re[i] * phiHat.im[i] + gamma.im[i] * phiHat.re[i];
        const gr = G.re[base + i];
        const gi = G.im[base + i];
        numRe += mr * gr + mi * gi;
        numIm += mr * gi - mi * gr;
        denom += mr * mr + mi * mi;
        total += gr * gr + gi * gi;
      }
      const cr = numRe / denom;
      const ci = numIm / denom;
      worstScale = Math.max(worstScale, Math.hypot(cr, ci - 1));

      let residual = 0;
      for (let i = 1; i < m * m; i++) {
        const mr = gamma.re[i] * phiHat.re[i] - gamma.im[i] * phiHat.im[i];
        const mi = gamma.re[i] * phiHat.im[i] + gamma.im[i] * phiHat.re[i];
        residual += (G.re[base + i] - (cr * mr - ci * mi)) ** 2
          + (G.im[base + i] - (cr * mi + ci * mr)) ** 2;
      }
      worstResidual = Math.max(worstResidual, Math.sqrt(residual / total));
    }

    results.push({amplitude, worstScale, worstResidual});
    t.diagnostic(
      `max|phi|=${amplitude.toExponential(0)}  |c - i|=${worstScale.toExponential(2)}` +
      `  residual=${worstResidual.toExponential(2)}`
    );
    assert.ok(worstScale < 1e-3, `the fitted scale is not i: off by ${worstScale}`);
    assert.ok(worstResidual < 5 * amplitude, `residual ${worstResidual} at phi=${amplitude}`);
  }

  // The residual is truncation, so weakening the object by ten must reduce it.
  const [strong, weak] = results;
  assert.ok(weak.worstResidual < strong.worstResidual / 3,
    `residual did not fall with the object strength: ${strong.worstResidual} -> ${weak.worstResidual}`);
});

test("shot noise is Poisson at the requested dose", () => {
  const n = 32;
  const m = 8;
  const sampling = [0.25, 0.25];
  const wavelength = electronWavelength(ENERGY);
  const semiangle = 4 * (1 / (n * sampling[0])) * wavelength * 1e3;
  const probe = complexProbe({
    gpts: [n, n], sampling, energy: ENERGY, semiangleCutoff: semiangle, aberrations: {}
  });
  const probeReal = ifft2({re: Float64Array.from(probe.re), im: Float64Array.from(probe.im)}, n, n);
  const phase = new Float64Array(n * n);

  // Deterministic stream, so a failure here is a real change and not luck.
  let seed = 12345;
  const random = () => {
    seed = (seed * 1664525 + 1013904223) >>> 0;
    return seed / 4294967296;
  };

  const clean = forwardModel({objectPhase: phase, probeReal, n, m, dose: Infinity});
  clean.step(clean.total);
  const noisy = forwardModel({
    objectPhase: phase, probeReal, n, m, dose: 1e3, sampling, random
  });
  noisy.step(noisy.total);

  let signal = 0;
  let noise = 0;
  for (let i = 0; i < clean.data.length; i++) {
    signal += clean.data[i];
    noise += (noisy.data[i] - clean.data[i]) ** 2;
  }
  assert.ok(signal > 0, "the noiseless scan is empty");
  assert.ok(noise > 0, "the dose slider changed nothing");

  // Raising the dose must reduce the noise, roughly as 1/dose in variance.
  seed = 12345;
  const brighter = forwardModel({
    objectPhase: phase, probeReal, n, m, dose: 1e5, sampling, random
  });
  brighter.step(brighter.total);
  let brighterNoise = 0;
  for (let i = 0; i < clean.data.length; i++) brighterNoise += (brighter.data[i] - clean.data[i]) ** 2;
  assert.ok(brighterNoise < noise / 10,
    `a hundredfold dose barely helped: ${noise} -> ${brighterNoise}`);
});

// ---------------------------------------------------------------------------
// The estimators
// ---------------------------------------------------------------------------

test("every direct method recovers a weak object from noiseless data", (t) => {
  const n = 32;
  const m = 32;
  const sampling = [0.25, 0.25];
  const wavelength = electronWavelength(ENERGY);
  const grid = angularSpatialFrequencies([n, n], sampling, wavelength);
  const semiangle = 4 * (1 / (n * sampling[0])) * wavelength * 1e3;
  const aberrations = {C10: -40};

  const probe = complexProbe({
    gpts: [n, n], sampling, energy: ENERGY, semiangleCutoff: semiangle, aberrations
  });
  const probeReal = ifft2({re: Float64Array.from(probe.re), im: Float64Array.from(probe.im)}, n, n);
  const chiArray = chi(grid.alpha, grid.phi, wavelength, aberrations);

  const phase = new Float64Array(n * n);
  for (let a = 0; a < n; a++) {
    for (let b = 0; b < n; b++) {
      phase[a * n + b] = 1e-3 * (Math.sin(0.55 * a) * Math.cos(0.42 * b) + 0.4 * Math.sin(0.3 * (a + b)));
    }
  }

  const bf = brightFieldIndices(grid.alpha, semiangle);
  const model = forwardModel({objectPhase: phase, probeReal, n, m, bfIndices: bf, dose: Infinity});
  model.step(model.total);
  const spectra = scanSpectra(model.data, m, bf.length);

  // iCOM and parallax need per-pixel geometry; build it once.
  const qRow = new Float64Array(m * m);
  const qCol = new Float64Array(m * m);
  const inverseQ2 = new Float64Array(m * m);
  const phaseFlip = new Float64Array(m * m);
  for (let a = 0; a < m; a++) {
    for (let b = 0; b < m; b++) {
      const i = a * m + b;
      qRow[i] = signed(a, m);
      qCol[i] = signed(b, m);
      const q2 = qRow[i] * qRow[i] + qCol[i] * qCol[i];
      inverseQ2[i] = q2 > 0 ? 1 / q2 : 0;
      phaseFlip[i] = Math.sign(Math.sin(chiArray[i])) || 1;
    }
  }

  const gamma = complex(m * m);
  const g = complex(m * m);

  for (const method of PTYCHO_METHODS) {
    const acc = directAccumulator(m);
    for (let j = 0; j < bf.length; j++) {
      const kIndex = bf[j];
      overlapFunction(probe, n, m, kIndex, gamma);
      spectrumAt(spectra, m, j, g);

      const kRow = signed((kIndex / n) | 0, n);
      const kCol = signed(kIndex % n, n);
      // Parallax shifts by the aberration gradient; for pure defocus that is
      // -lambda * C10 * k, in pixels of the reconstruction grid.
      const shift = -wavelength * aberrations.C10 * (1 / (n * sampling[0]));
      accumulatePixel(acc, g, gamma, method, {
        kRow, kCol, qRow, qCol, inverseQ2, phaseFlip,
        shiftRow: shift * kRow, shiftCol: shift * kCol
      });
    }
    const image = directImage(acc, method);

    // Correlate against the truth: these estimators recover the phase up to a
    // method-dependent scale and an additive constant, which is what a CTF is.
    let meanA = 0;
    let meanB = 0;
    for (let i = 0; i < m * m; i++) {
      meanA += phase[i];
      meanB += image[i];
    }
    meanA /= m * m;
    meanB /= m * m;
    let cov = 0;
    let varA = 0;
    let varB = 0;
    for (let i = 0; i < m * m; i++) {
      const da = phase[i] - meanA;
      const db = image[i] - meanB;
      cov += da * db;
      varA += da * da;
      varB += db * db;
    }
    const correlation = cov / Math.sqrt(varA * varB);
    t.diagnostic(`${method}: correlation with the truth ${correlation.toFixed(4)}`);

    // Parallax reconstructs the phase with reversed contrast: its transfer
    // function is `-sin[chi] (A*A)`, and the global phase flip makes that
    // `-|sin chi| (A*A)`, which is negative everywhere. So compare magnitudes,
    // and check the sign separately rather than letting it pass unnoticed.
    assert.ok(Math.abs(correlation) > 0.9,
      `${method} recovered a phase only ${correlation.toFixed(3)} correlated with the object`);
    const expectedSign = method === "parallax" ? -1 : 1;
    assert.equal(Math.sign(correlation), expectedSign,
      `${method} recovered the phase with the wrong sign`);
  }
});

test("every method survives Fourier-tiling upsampling", async (t) => {
  // A coarse scan reconstructed onto a finer grid. Tiling replicates the scan
  // spectrum past its own Nyquist frequency, and the kernel has to reject the
  // copies it did not measure.
  //
  // Tiling replicates the scan spectrum past its own Nyquist frequency, and
  // every kernel has to reject the copies it did not measure.
  //
  // Parallax is the method this can break, because its kernel has constant
  // magnitude and so nothing in it suppresses a copy on size alone. What used
  // to break it was not the copies of the signal but the copies of the
  // unscattered beam: leaving the zero-frequency term in the virtual
  // bright-field spectra put an enormous spike at every tile corner, and
  // parallax passed all of them. `scanSpectra` drops that term, and with it
  // gone parallax tracks the others.
  const n = 64;
  const m = 32;
  const sampling = [0.25, 0.25];
  const wavelength = electronWavelength(ENERGY);
  const dk = 1 / (n * sampling[0]);
  const radius = 8;
  const semiangle = radius * dk * wavelength * 1e3;
  const gridOf = angularSpatialFrequencies([n, n], sampling, wavelength);
  const bf = brightFieldIndices(gridOf.alpha, semiangle);
  const aberrations = {C10: -600};

  const probe = complexProbe({
    gpts: [n, n], sampling, energy: ENERGY, semiangleCutoff: semiangle, aberrations
  });
  const probeReal = ifft2({re: Float64Array.from(probe.re), im: Float64Array.from(probe.im)}, n, n);
  const chiDetector = chi(gridOf.alpha, gridOf.phi, wavelength, aberrations);

  const truth = whiteNoiseObject2D(n, 0.03, mulberry32(9));
  const model = forwardModel({
    objectPhase: truth, probeReal, n, m, bfIndices: bf, dose: Infinity
  });
  model.step(model.total);
  const spectra = scanSpectra(model.data, m, bf.length);

  for (const upsample of [1, 2]) {
    const Q = m * upsample;
    const qSampling = (n * sampling[0]) / Q;
    const qGrid = angularSpatialFrequencies([Q, Q], [qSampling, qSampling], wavelength);
    const chiQ = chi(qGrid.alpha, qGrid.phi, wavelength, aberrations);
    const qRow = new Float64Array(Q * Q);
    const qCol = new Float64Array(Q * Q);
    const inverseQ2 = new Float64Array(Q * Q);
    const phaseFlip = new Float64Array(Q * Q);
    for (let a = 0; a < Q; a++) {
      for (let b = 0; b < Q; b++) {
        const i = a * Q + b;
        qRow[i] = signed(a, Q);
        qCol[i] = signed(b, Q);
        const q2 = qRow[i] * qRow[i] + qCol[i] * qCol[i];
        inverseQ2[i] = q2 > 0 ? 1 / q2 : 0;
        phaseFlip[i] = Math.sign(Math.sin(chiQ[i])) || 1;
      }
    }
    const shifts = parallaxShifts(chiDetector, n, dk, qSampling);

    // The truth, low-passed at twice the semiangle and binned onto the
    // reconstruction grid.
    //
    // The band limit is what makes the number mean anything. A white-noise
    // object carries power at every frequency, most of it past the transfer
    // limit of every method here, so comparing against the raw object measures
    // mostly what none of them could ever recover. Against what is actually
    // reachable, a good reconstruction scores near one.
    const bandLimited = {re: Float64Array.from(truth), im: new Float64Array(n * n)};
    fft2(bandLimited, n, n);
    for (let a = 0; a < n; a++) {
      for (let b = 0; b < n; b++) {
        if (Math.hypot(signed(a, n), signed(b, n)) <= 2 * radius) continue;
        const i = a * n + b;
        bandLimited.re[i] = 0;
        bandLimited.im[i] = 0;
      }
    }
    ifft2(bandLimited, n, n);

    const factor = n / Q;
    const reference = new Float64Array(Q * Q);
    for (let a = 0; a < Q; a++) {
      for (let b = 0; b < Q; b++) {
        let total = 0;
        for (let i = 0; i < factor; i++) {
          for (let j = 0; j < factor; j++) {
            total += bandLimited.re[(a * factor + i) * n + (b * factor + j)];
          }
        }
        reference[a * Q + b] = total / (factor * factor);
      }
    }

    const gamma = complex(Q * Q);
    const spectrum = complex(m * m);
    const tiled = complex(Q * Q);
    for (const method of PTYCHO_METHODS) {
      const acc = directAccumulator(Q);
      for (let j = 0; j < bf.length; j++) {
        const kIndex = bf[j];
        overlapFunction(probe, n, Q, kIndex, gamma);
        spectrumAt(spectra, m, j, spectrum);
        tileSpectrum(spectrum, m, upsample, tiled);
        accumulatePixel(acc, tiled, gamma, method, {
          kRow: signed((kIndex / n) | 0, n), kCol: signed(kIndex % n, n),
          qRow, qCol, inverseQ2, phaseFlip, bandLimit: 2 * radius,
          shiftRow: shifts.shiftRow[kIndex], shiftCol: shifts.shiftCol[kIndex]
        });
      }
      const image = directImage(acc, method);

      let meanA = 0;
      let meanB = 0;
      for (let i = 0; i < Q * Q; i++) {
        meanA += reference[i];
        meanB += image[i];
      }
      meanA /= Q * Q;
      meanB /= Q * Q;
      let cov = 0;
      let varA = 0;
      let varB = 0;
      for (let i = 0; i < Q * Q; i++) {
        const da = reference[i] - meanA;
        const db = image[i] - meanB;
        cov += da * db;
        varA += da * da;
        varB += db * db;
      }
      const correlation = Math.abs(cov / Math.sqrt(varA * varB));
      t.diagnostic(`f=${upsample} ${method}: |r| = ${correlation.toFixed(3)}`);
      // Against a band-limited reference, so this is how much of the reachable
      // object each method actually recovered.
      //
      // iCOM gets its own bar, for a reason that is about the object rather
      // than the method. Its kernel carries a factor of `1 / |q|`, so against a
      // flat-spectrum object the result is dominated by the few lowest
      // frequencies and correlates poorly with the whole. On a specimen whose
      // own spectrum falls away — apoferritin, in the widget — the same code
      // reaches 0.93.
      const floor = method === "icom" ? 0.15 : 0.6;
      assert.ok(correlation > floor,
        `${method} at upsampling ${upsample} recovered |r| = ${correlation.toFixed(3)}`);
    }
  }
});

test("OBF and the matched filter carry the same information", (t) => {
  // The paper's non-obvious claim: the matched filter's higher transfer is
  // exactly offset by its higher variance, so the two have identical SSNR. A
  // wrong normalisation in either one breaks this.
  const probe = complexProbe({
    gpts: [N, N], sampling: SAMPLING, energy: ENERGY, semiangleCutoff: SEMIANGLE,
    aberrations: {C10: -60, C30: 1e4}
  });
  const aperture = softAperture(geometry.alpha, geometry.phi, SEMIANGLE, SAMPLING_MRAD);
  const chiArray = chi(geometry.alpha, geometry.phi, WAVELENGTH, {C10: -60, C30: 1e4});
  const bf = brightFieldIndices(geometry.alpha, SEMIANGLE);
  const sums = overlapSums(probe, N, N, bf);
  const context = {probe, aperture, chiArray, n: N, sums};

  const obf = directSSNR("obf", context);
  const mf = directSSNR("matchedFilter", context);
  let worst = 0;
  for (let i = 0; i < N * N; i++) worst = Math.max(worst, Math.abs(obf[i] - mf[i]));
  t.diagnostic(`max difference ${worst.toExponential(2)}, scale ${maxAbs(obf).toExponential(2)}`);
  assert.ok(worst < 1e-15, `OBF and matched-filter SSNR differ by ${worst}`);

  // And SSB must be close but never better: by Cauchy-Schwarz its coherent sum
  // over equally weighted pixels cannot beat the noise-flattened one. An SSB
  // curve that crept above OBF would mean a normalisation error in one of them.
  const ssb = directSSNR("ssb", context);
  let excess = 0;
  let difference = 0;
  for (let i = 0; i < N * N; i++) {
    excess = Math.max(excess, ssb[i] - obf[i]);
    difference = Math.max(difference, Math.abs(ssb[i] - obf[i]));
  }
  t.diagnostic(`SSB exceeds OBF by at most ${excess.toExponential(2)}`);
  assert.ok(excess < 1e-12, `SSB SSNR beats OBF by ${excess}, which cannot happen`);
  assert.ok(difference > 1e-6, "SSB and OBF came out identical, which they should not be");
});

test("the transfer functions share a usable scale", (t) => {
  // The five curves go on one axis, so they have to be within an order of
  // magnitude of each other. They are not naturally: the review writes each
  // estimator in its own normalisation, and OBF's differs from SSB's by the
  // square root of the bright-field pixel count — a factor of fourteen here,
  // which would flatten four of the five curves against the axis.
  const aberrations = {C10: -60, C30: 1e4};
  const probe = complexProbe({
    gpts: [N, N], sampling: SAMPLING, energy: ENERGY, semiangleCutoff: SEMIANGLE, aberrations
  });
  const context = {
    probe,
    aperture: softAperture(geometry.alpha, geometry.phi, SEMIANGLE, SAMPLING_MRAD),
    chiArray: chi(geometry.alpha, geometry.phi, WAVELENGTH, aberrations),
    n: N,
    sums: overlapSums(probe, N, N, brightFieldIndices(geometry.alpha, SEMIANGLE))
  };

  const peaks = PTYCHO_METHODS.map((method) => maxAbs(directCTF(method, context)));
  t.diagnostic(PTYCHO_METHODS.map((m, i) => `${m} ${peaks[i].toFixed(3)}`).join(", "));
  const smallest = Math.min(...peaks);
  const largest = Math.max(...peaks);
  assert.ok(largest / smallest < 10,
    `CTF peaks span ${(largest / smallest).toFixed(1)}x: ` +
    PTYCHO_METHODS.map((m, i) => `${m}=${peaks[i].toExponential(2)}`).join(" "));
});

test("the analytic CTFs match the estimators' definitions", (t) => {
  const aberrations = {C10: -60};
  const probe = complexProbe({
    gpts: [N, N], sampling: SAMPLING, energy: ENERGY, semiangleCutoff: SEMIANGLE, aberrations
  });
  const aperture = softAperture(geometry.alpha, geometry.phi, SEMIANGLE, SAMPLING_MRAD);
  const chiArray = chi(geometry.alpha, geometry.phi, WAVELENGTH, aberrations);
  const bf = brightFieldIndices(geometry.alpha, SEMIANGLE);
  const sums = overlapSums(probe, N, N, bf);
  const context = {probe, aperture, chiArray, n: N, sums};

  let weight = 0;
  for (const v of aperture) weight += v * v;

  const ssb = directCTF("ssb", context);
  const obf = directCTF("obf", context);
  let worst = 0;
  for (let i = 0; i < N * N; i++) {
    worst = Math.max(
      worst,
      Math.abs(ssb[i] - sums.sum[i] / 2),
      Math.abs(obf[i] - Math.sqrt(weight * sums.energy[i]) / 2)
    );
  }
  assert.ok(worst < 1e-15, `CTF definitions disagree by ${worst}`);

  // The scale of the transfer functions, checked where it is pinned: in focus,
  // in the single-overlap annulus, CTF_ssb is the autocorrelation of the
  // L2-normalised aperture. Normalising by the aperture weight a second time
  // would put this three orders of magnitude low, and nothing about the shape
  // of the curve would give it away.
  //
  // Note this does not say the CTF peaks at one. It does not: the triple
  // overlap cancels in focus, so CTF_ssb falls to zero at zero frequency and
  // peaks near 0.41 just outside the aperture radius.
  const inFocusProbe = complexProbe({
    gpts: [N, N], sampling: SAMPLING, energy: ENERGY,
    semiangleCutoff: SEMIANGLE, aberrations: {}, softEdges: false
  });
  const hard = hardAperture(geometry.alpha, SEMIANGLE);
  const inFocus = directCTF("ssb", {
    probe: inFocusProbe,
    aperture: hard,
    chiArray: new Float64Array(N * N),
    n: N,
    sums: overlapSums(inFocusProbe, N, N, bf)
  });

  const envelope = {re: l2Normalised(hard), im: new Float64Array(N * N)};
  fft2(envelope, N, N);
  for (let i = 0; i < N * N; i++) {
    envelope.re[i] = envelope.re[i] * envelope.re[i] + envelope.im[i] * envelope.im[i];
    envelope.im[i] = 0;
  }
  ifft2(envelope, N, N);

  let annulus = 0;
  for (let a = 0; a < N; a++) {
    for (let b = 0; b < N; b++) {
      const i = a * N + b;
      const q = Math.hypot(signed(a, N), signed(b, N));
      if (q <= RADIUS + 1.5 || q >= 2 * RADIUS - 1.5) continue;
      annulus = Math.max(annulus, Math.abs(inFocus[i] - envelope.re[i]));
    }
  }
  t.diagnostic(
    `in-focus SSB CTF: peak ${maxAbs(inFocus).toFixed(4)}, ` +
    `annulus error against the aperture autocorrelation ${annulus.toExponential(2)}`
  );
  assert.ok(annulus < 1e-13, `the SSB CTF is mis-scaled: off by ${annulus} in the annulus`);

  // Every method must produce something finite and non-trivial.
  for (const method of PTYCHO_METHODS) {
    const ctf = directCTF(method, context);
    assert.ok(ctf.every(Number.isFinite), `${method} CTF has a non-finite entry`);
    const peak = maxAbs(ctf);
    t.diagnostic(`${method}: peak |CTF| ${peak.toFixed(4)}`);
    assert.ok(peak > 1e-6, `${method} CTF is identically zero`);
  }
});

test("parallax approximates SSB for pure defocus, in the single-overlap annulus", (t) => {
  // Parallax keeps only the first-order term of Gamma, so
  // CTF_prlx = |sin chi| CTF_SSB where a single overlap contributes. Outside
  // that annulus the approximation is 20% or so wrong, which is the honest
  // statement and the reason the restriction is part of the test.
  const aberrations = {C10: -60};
  const probe = complexProbe({
    gpts: [N, N], sampling: SAMPLING, energy: ENERGY,
    semiangleCutoff: SEMIANGLE, aberrations, softEdges: false
  });
  const aperture = hardAperture(geometry.alpha, SEMIANGLE);
  const chiArray = chi(geometry.alpha, geometry.phi, WAVELENGTH, aberrations);
  const bf = brightFieldIndices(geometry.alpha, SEMIANGLE);
  const sums = overlapSums(probe, N, N, bf);
  // Without the phase flip, so the comparison is against `-sin chi` rather
  // than its magnitude.
  const context = {probe, aperture, chiArray, n: N, sums, phaseFlip: false};

  const ssb = directCTF("ssb", context);
  const parallax = directCTF("parallax", context);

  let annulus = 0;
  let scale = 0;
  for (let a = 0; a < N; a++) {
    for (let b = 0; b < N; b++) {
      const i = a * N + b;
      const q = Math.hypot(signed(a, N), signed(b, N));
      if (q <= RADIUS + 1.5 || q >= 2 * RADIUS - 1.5) continue;
      const expected = -Math.sign(Math.sin(chiArray[i])) * Math.abs(Math.sin(chiArray[i])) * ssb[i];
      annulus = Math.max(annulus, Math.abs(parallax[i] - expected));
      scale = Math.max(scale, Math.abs(expected));
    }
  }
  t.diagnostic(`annulus error ${annulus.toExponential(2)} against a scale of ${scale.toExponential(2)}`);
  assert.ok(annulus < 1e-13, `parallax and SSB disagree by ${annulus} inside the annulus`);
});
