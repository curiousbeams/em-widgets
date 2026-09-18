import test from "node:test";
import assert from "node:assert/strict";

import {electronWavelength} from "../kit/units.js";
import {angularSpatialFrequencies} from "../kit/grid.js";
import {complexProbe} from "../kit/optics.js";
import {mulberry32} from "../kit/specimen.js";
import {poisson} from "../kit/image.js";
import {complex, fft2, ifft2} from "../kit/fft.js";
import {forwardModel} from "../kit/ptycho-sim.js";
import {
  epieState, epieStep, epieReset, reconstructedPhase, scanOrder, probeDiameter,
  centreSpectrum, epieLine
} from "../kit/epie.js";

// Small enough that a few sweeps are milliseconds, big enough that the probe
// covers several pixels and the scan genuinely overlaps.
const N = 32;
const M = 16;                     // stride 2
const SAMPLING = [0.25, 0.25];
const ENERGY = 300e3;
const WAVELENGTH = electronWavelength(ENERGY);
const DK = 1 / (N * SAMPLING[0]);
const SEMIANGLE = 4 * DK * WAVELENGTH * 1e3;   // a four-pixel aperture radius

const geometry = angularSpatialFrequencies([N, N], SAMPLING, WAVELENGTH);

function realProbe(aberrations) {
  const probe = complexProbe({
    gpts: [N, N], sampling: SAMPLING, energy: ENERGY, semiangleCutoff: SEMIANGLE, aberrations
  });
  return ifft2({re: Float64Array.from(probe.re), im: Float64Array.from(probe.im)}, N, N);
}

/** A smooth weak phase object, band-limited well inside the detector. */
function testObject(amplitude = 0.15) {
  const phase = new Float64Array(N * N);
  let peak = 0;
  for (let a = 0; a < N; a++) {
    for (let b = 0; b < N; b++) {
      const x = (2 * Math.PI * a) / N;
      const y = (2 * Math.PI * b) / N;
      const v =
        Math.cos(2 * x + 0.3) * Math.sin(3 * y) +
        0.6 * Math.cos(5 * x - 1.1) +
        0.4 * Math.sin(4 * y + 0.7) * Math.cos(x);
      phase[a * N + b] = v;
      peak = Math.max(peak, Math.abs(v));
    }
  }
  let mean = 0;
  for (const v of phase) mean += v;
  mean /= phase.length;
  for (let i = 0; i < phase.length; i++) phase[i] = ((phase[i] - mean) * amplitude) / peak;
  return phase;
}

function amplitudesFor(phase, probeReal, dose = Infinity, seed = 7) {
  const model = forwardModel({
    objectPhase: phase, probeReal, n: N, m: M,
    amplitude: true, layout: "rk", dose, sampling: SAMPLING,
    random: mulberry32(seed)
  });
  assert.equal(model.step(model.total), 0, "the scan did not finish");
  return model.data;
}

/** Root-mean-square difference of two phase maps, each mean-subtracted. */
function phaseError(a, b) {
  const mean = (x) => x.reduce((s, v) => s + v, 0) / x.length;
  const ma = mean(a);
  const mb = mean(b);
  let sum = 0;
  for (let i = 0; i < a.length; i++) {
    const d = a[i] - ma - (b[i] - mb);
    sum += d * d;
  }
  return Math.sqrt(sum / a.length);
}

// ---------------------------------------------------------------------------
// The data layout
// ---------------------------------------------------------------------------

test("the two forward-model layouts hold the same numbers", () => {
  const phase = testObject();
  const probeReal = realProbe({C10: -30});
  const options = {objectPhase: phase, probeReal, n: N, m: M, amplitude: true, dose: Infinity};

  const kr = forwardModel({...options, layout: "kr"});
  const rk = forwardModel({...options, layout: "rk"});
  kr.step(kr.total);
  rk.step(rk.total);

  const size = N * N;
  const total = M * M;
  let worst = 0;
  for (let position = 0; position < total; position++) {
    for (let pixel = 0; pixel < size; pixel++) {
      worst = Math.max(worst, Math.abs(rk.data[position * size + pixel] - kr.data[pixel * total + position]));
    }
  }
  assert.equal(worst, 0, `the layouts disagree by ${worst}`);
});

// ---------------------------------------------------------------------------
// One step
// ---------------------------------------------------------------------------

test("a vacuum object is already the answer, so the first step corrects nothing", (t) => {
  // The simulator and the reconstruction have to agree on the probe's roll and
  // on the transform normalisation. If either differed, the modelled amplitude
  // would not match the measured one even when the guess is exactly right.
  const probeReal = realProbe({C10: -30});
  const vacuum = new Float64Array(N * N);
  const amplitudes = amplitudesFor(vacuum, probeReal);

  const state = epieState({n: N, m: M, probe: probeReal});
  let worst = 0;
  for (const position of [0, 5, M + 3, M * M - 1]) {
    const {error} = epieStep(state, amplitudes, position, {beta: 0.9});
    worst = Math.max(worst, error);
  }
  // Not exactly zero: the simulator stores amplitudes as Float32, so the
  // measurement the reconstruction is compared against has been rounded to
  // about seven digits. That rounding is the whole of what is left here.
  t.diagnostic(`residual ${worst.toExponential(2)}`);
  assert.ok(worst < 1e-6, `a correct guess still moved, by ${worst}`);
});

test("the correction is the amplitude residual, and it is zero only when matched", () => {
  const phase = testObject();
  const probeReal = realProbe({C10: -30});
  const amplitudes = amplitudesFor(phase, probeReal);

  const state = epieState({n: N, m: M, probe: probeReal});
  const {error} = epieStep(state, amplitudes, 0, {beta: 0});   // measure, do not move
  assert.ok(error > 1e-3, `a vacuum guess against a real object should be wrong, got ${error}`);

  // beta = 0 means the object is untouched, so the same step repeats exactly.
  const {error: again} = epieStep(state, amplitudes, 0, {beta: 0});
  assert.equal(again, error);
});

// ---------------------------------------------------------------------------
// Convergence
// ---------------------------------------------------------------------------

test("ePIE recovers a known weak phase object from noiseless data", (t) => {
  const phase = testObject();
  const probeReal = realProbe({C10: -30});
  const amplitudes = amplitudesFor(phase, probeReal);

  const state = epieState({n: N, m: M, probe: probeReal});
  const order = scanOrder(M * M, mulberry32(3));
  for (let sweep = 0; sweep < 12; sweep++) {
    for (const position of order) epieStep(state, amplitudes, position, {beta: 0.9});
  }

  const recovered = reconstructedPhase(state);
  const error = phaseError(recovered, phase);
  let rms = 0;
  for (const v of phase) rms += v * v;
  rms = Math.sqrt(rms / phase.length);

  t.diagnostic(`residual ${(error / rms * 100).toFixed(2)}% of the object's own rms`);
  assert.ok(error < 0.02 * rms, `recovered phase is off by ${(error / rms * 100).toFixed(1)}%`);
});

test("the amplitude error falls sweep on sweep", (t) => {
  const phase = testObject();
  const probeReal = realProbe({C10: -30});
  const amplitudes = amplitudesFor(phase, probeReal);

  const state = epieState({n: N, m: M, probe: probeReal});
  const order = scanOrder(M * M, mulberry32(11));
  const history = [];
  for (let sweep = 0; sweep < 8; sweep++) {
    let total = 0;
    for (const position of order) {
      total += epieStep(state, amplitudes, position, {beta: 0.5}).error;
    }
    history.push(total / order.length);
  }

  t.diagnostic(history.map((e) => e.toExponential(2)).join(" -> "));
  for (let i = 1; i < history.length; i++) {
    assert.ok(history[i] < history[i - 1], `sweep ${i} did not improve on sweep ${i - 1}`);
  }
});

test("shot noise leaves a floor the reconstruction cannot get under", (t) => {
  const phase = testObject();
  const probeReal = realProbe({C10: -30});

  const errors = [];
  for (const dose of [1e3, Infinity]) {
    const amplitudes = amplitudesFor(phase, probeReal, dose, 19);
    const state = epieState({n: N, m: M, probe: probeReal});
    const order = scanOrder(M * M, mulberry32(5));
    for (let sweep = 0; sweep < 12; sweep++) {
      for (const position of order) epieStep(state, amplitudes, position, {beta: 0.7});
    }
    errors.push(phaseError(reconstructedPhase(state), phase));
  }

  t.diagnostic(`1e3 e/A^2: ${errors[0].toExponential(2)} rad, noiseless: ${errors[1].toExponential(2)} rad`);
  assert.ok(errors[0] > 5 * errors[1], "noise made no difference, which cannot be right");
});

// ---------------------------------------------------------------------------
// The potential parameterisation
// ---------------------------------------------------------------------------

/** A non-negative phase object, as a projected potential is. */
function positiveObject(peak = 0.5) {
  const phase = testObject(1);
  let min = Infinity;
  let max = -Infinity;
  for (const v of phase) {
    if (v < min) min = v;
    if (v > max) max = v;
  }
  for (let i = 0; i < phase.length; i++) phase[i] = ((phase[i] - min) / (max - min)) * peak;
  return phase;
}

const rmsOf = (array) => Math.sqrt(array.reduce((s, v) => s + v * v, 0) / array.length);

/** Root-mean-square difference with the constant left in. */
function absoluteError(a, b) {
  let sum = 0;
  for (let i = 0; i < a.length; i++) sum += (a[i] - b[i]) ** 2;
  return Math.sqrt(sum / a.length);
}

function runPotential(phase, probeReal, amplitudes, options, sweeps = 40) {
  const state = epieState({n: N, m: M, probe: probeReal, ...options});
  const order = scanOrder(M * M, mulberry32(3));
  for (let sweep = 0; sweep < sweeps; sweep++) {
    for (const position of order) epieStep(state, amplitudes, position, {beta: 0.9});
  }
  return state;
}

test("solving for the potential recovers the same object as solving for O", (t) => {
  const phase = positiveObject();
  const probeReal = realProbe({C10: -30});
  const amplitudes = amplitudesFor(phase, probeReal);
  const rms = rmsOf(phase);

  const complexObject = phaseError(
    reconstructedPhase(runPotential(phase, probeReal, amplitudes, {}, 12)), phase);
  const potential = phaseError(
    reconstructedPhase(runPotential(phase, probeReal, amplitudes, {objectType: "potential"}, 12)),
    phase);

  t.diagnostic(
    `complex ${(complexObject / rms * 100).toFixed(2)}%, potential ${(potential / rms * 100).toFixed(2)}%`
  );
  assert.ok(potential < 0.02 * rms, `the potential is off by ${(potential / rms * 100).toFixed(1)}%`);
  assert.ok(complexObject < 0.02 * rms, "the complex object stopped agreeing with it");
});

test("positivity recovers the constant that ptychography cannot otherwise determine", (t) => {
  // The data fixes `exp(i phi)`, so it fixes `phi` only up to an additive
  // constant — every estimate of the shape is equally good, and an
  // unconstrained run settles on whichever the iteration drifts to. Saying the
  // specimen cannot advance the wave picks one out. It costs a little accuracy
  // in the shape and buys back the absolute scale, which is the trade worth
  // knowing about.
  const phase = positiveObject();
  const probeReal = realProbe({C10: -30});
  const amplitudes = amplitudesFor(phase, probeReal);
  const rms = rmsOf(phase);

  const measure = (state) => {
    const recovered = reconstructedPhase(state);
    let offset = 0;
    for (let i = 0; i < recovered.length; i++) offset += recovered[i] - phase[i];
    return {
      offset: offset / recovered.length,
      shape: phaseError(recovered, phase),
      absolute: absoluteError(recovered, phase)
    };
  };

  const free = measure(runPotential(phase, probeReal, amplitudes, {objectType: "potential"}));
  const held = measure(
    runPotential(phase, probeReal, amplitudes, {objectType: "potential", positivity: true}));

  for (const [name, r] of [["free", free], ["positive", held]]) {
    t.diagnostic(
      `${name.padEnd(8)} offset ${r.offset.toFixed(4)} rad, ` +
        `shape ${(r.shape / rms * 100).toFixed(2)}%, absolute ${(r.absolute / rms * 100).toFixed(1)}%`
    );
  }
  assert.ok(Math.abs(free.offset) > 0.5 * rms, "the unconstrained run happened not to drift");
  assert.ok(free.shape < 0.01 * rms, "the unconstrained run should have the shape exactly right");
  assert.ok(held.absolute < 0.3 * rms, `positivity left an absolute error of ${held.absolute}`);
  assert.ok(held.absolute < free.absolute / 3, "positivity did not recover the constant");
});

test("positivity holds the phase at or above zero", () => {
  const phase = positiveObject();
  const probeReal = realProbe({C10: -30});
  const amplitudes = amplitudesFor(phase, probeReal, 1e3, 31);   // noisy, so it pushes

  const state = epieState({n: N, m: M, probe: probeReal, objectType: "potential", positivity: true});
  const order = scanOrder(M * M, mulberry32(3));
  for (let sweep = 0; sweep < 6; sweep++) {
    for (const position of order) epieStep(state, amplitudes, position, {beta: 0.9});
  }

  let least = Infinity;
  for (const v of state.phase) if (v < least) least = v;
  assert.ok(least >= 0, `positivity let the phase reach ${least}`);

  // And the constraint has to be doing something: without it, noise takes the
  // phase negative somewhere.
  const free = epieState({n: N, m: M, probe: probeReal, objectType: "potential"});
  for (let sweep = 0; sweep < 6; sweep++) {
    for (const position of order) epieStep(free, amplitudes, position, {beta: 0.9});
  }
  let freeLeast = Infinity;
  for (const v of free.phase) if (v < freeLeast) freeLeast = v;
  assert.ok(freeLeast < 0, "unconstrained, the phase never went negative — the test proves nothing");
});

test("a potential object keeps its phase unwrapped", () => {
  // `phase` is the array being solved for and `object` is `exp(i phase)` kept in
  // step with it, so a phase past pi is held as itself. Reading it back off the
  // complex object — which is all a `"complex"` reconstruction can do — folds it
  // into (-pi, pi], and a panel of that shows a seam where the specimen is
  // thickest rather than a peak.
  const phase = positiveObject(0.5);
  const probeReal = realProbe({C10: -30});
  const amplitudes = amplitudesFor(phase, probeReal);

  const state = epieState({n: N, m: M, probe: probeReal, objectType: "potential"});
  epieStep(state, amplitudes, 0, {beta: 0.9});

  // The invariant, wherever the reconstruction has got to.
  let worst = 0;
  for (let i = 0; i < state.phase.length; i++) {
    worst = Math.max(
      worst,
      Math.abs(state.object.re[i] - Math.cos(state.phase[i])),
      Math.abs(state.object.im[i] - Math.sin(state.phase[i]))
    );
  }
  assert.ok(worst < 1e-15, `exp(i phi) drifted from phi by ${worst}`);

  // And it survives a value the complex object could not carry.
  state.phase[0] = 4;
  state.object.re[0] = Math.cos(4);
  state.object.im[0] = Math.sin(4);
  assert.equal(reconstructedPhase(state)[0], 4);
  assert.ok(Math.atan2(state.object.im[0], state.object.re[0]) < 0,
    "atan2 was supposed to wrap this one");
});

// ---------------------------------------------------------------------------
// The probe
// ---------------------------------------------------------------------------

test("a zero step size computes everything and changes nothing", () => {
  // What a widget's inspect mode is: run the whole forward and backward
  // calculation at a position, fill in the stages, and leave the reconstruction
  // exactly where it was.
  const phase = positiveObject();
  const probeReal = realProbe({C10: -30});
  const amplitudes = amplitudesFor(phase, probeReal);

  const state = epieState({n: N, m: M, probe: probeReal, objectType: "potential", positivity: true});
  const order = scanOrder(M * M, mulberry32(3));
  for (let sweep = 0; sweep < 3; sweep++) {
    for (const position of order) epieStep(state, amplitudes, position, {beta: 0.9, probeBeta: 0.9});
  }

  const before = {
    phase: Float64Array.from(state.phase),
    objectRe: Float64Array.from(state.object.re),
    objectIm: Float64Array.from(state.object.im),
    probeRe: Float64Array.from(state.probe.re)
  };
  const {error} = epieStep(state, amplitudes, 7, {beta: 0, probeBeta: 0});

  assert.deepEqual(state.phase, before.phase);
  assert.deepEqual(state.object.re, before.objectRe);
  assert.deepEqual(state.object.im, before.objectIm);
  assert.deepEqual(state.probe.re, before.probeRe);

  // And it is not a no-op: the stages hold this position's answer.
  assert.ok(error > 0, "a zero step size should still measure the residual");
  assert.ok(state.stages.gradient.some((v) => v !== 0), "the gradient was not computed");
  assert.ok(state.stages.exit.re.some((v) => v !== 0), "the exit wave was not computed");
});

test("a known probe is left exactly alone", () => {
  // The exit-wave panel of a widget changes as the object does, which can read
  // as the probe drifting. It does not: at `probeBeta` zero the probe is
  // untouched, bit for bit, however many positions have been visited.
  const phase = positiveObject();
  const probeReal = realProbe({C10: -30});
  const amplitudes = amplitudesFor(phase, probeReal, 1e3, 13);   // noisy, so it would drift

  const state = epieState({n: N, m: M, probe: probeReal, objectType: "potential", positivity: true});
  const before = {re: Float64Array.from(state.probe.re), im: Float64Array.from(state.probe.im)};
  const order = scanOrder(M * M, mulberry32(3));
  for (let sweep = 0; sweep < 10; sweep++) {
    for (const position of order) epieStep(state, amplitudes, position, {beta: 0.9});
  }
  assert.ok(state.visited > 2000, "the run was too short to prove anything");
  assert.deepEqual(state.probe.re, before.re);
  assert.deepEqual(state.probe.im, before.im);
  // And the exit wave really has moved, so the test is not passing vacuously.
  assert.ok(state.stages.exit.re.some((v, i) => v !== probeReal.re[i]));
});

test("refining the probe improves a reconstruction started from the wrong one", (t) => {
  // A wrong probe is not fatal, because the object can absorb some of the
  // difference — that ambiguity is why ePIE converges at all from a bad guess,
  // and why the improvement here is a factor rather than orders of magnitude.
  // Both the residual it minimises and the object it recovers get better.
  const phase = testObject();
  const truth = realProbe({C10: -30});
  const wrong = realProbe({C10: -90});
  const amplitudes = amplitudesFor(phase, truth);

  const run = (probeBeta) => {
    const state = epieState({n: N, m: M, probe: wrong});
    const order = scanOrder(M * M, mulberry32(23));
    let residual = 0;
    for (let sweep = 0; sweep < 25; sweep++) {
      residual = 0;
      for (const position of order) {
        residual += epieStep(state, amplitudes, position, {beta: 0.6, probeBeta}).error;
      }
      residual /= order.length;
    }
    return {residual, phase: phaseError(reconstructedPhase(state), phase)};
  };

  const held = run(0);
  const refined = run(0.6);
  t.diagnostic(
    `probe held: residual ${held.residual.toExponential(2)}, phase ${held.phase.toExponential(2)} rad`
  );
  t.diagnostic(
    `refined:    residual ${refined.residual.toExponential(2)}, phase ${refined.phase.toExponential(2)} rad`
  );
  assert.ok(refined.residual < held.residual / 4, "refining the probe did not fit the data better");
  assert.ok(refined.phase < 0.8 * held.phase, "refining the probe did not improve the object");
});

test("the probe diameter follows the defocus", (t) => {
  const diameters = [0, -40, -120].map((C10) => probeDiameter(realProbe({C10}), N));
  t.diagnostic(`C10 = 0, -40, -120 A -> ${diameters.join(", ")} px`);
  assert.ok(diameters[0] > 0 && diameters[0] < N, "a focused probe should be a few pixels across");
  assert.ok(diameters[1] > diameters[0], "defocus should spread the probe");
  assert.ok(diameters[2] > diameters[1], "more defocus should spread it further");
});

test("a reset returns a reconstruction to exactly where it started", () => {
  const phase = positiveObject();
  const truth = realProbe({C10: -30});
  const wrong = realProbe({C10: -90});
  const amplitudes = amplitudesFor(phase, truth);

  const state = epieState({n: N, m: M, probe: wrong, objectType: "potential", positivity: true});
  const fresh = epieState({n: N, m: M, probe: wrong, objectType: "potential", positivity: true});

  // Refining, so that the probe moves too and the reset has to put it back.
  const order = scanOrder(M * M, mulberry32(3));
  for (const position of order) epieStep(state, amplitudes, position, {beta: 0.9, probeBeta: 0.9});
  assert.ok(state.phase.some((v) => v !== 0), "the run did nothing to reset");

  epieReset(state, {probe: wrong});
  assert.deepEqual(state.phase, fresh.phase);
  assert.deepEqual(state.object.re, fresh.object.re);
  assert.deepEqual(state.object.im, fresh.object.im);
  assert.deepEqual(state.probe.re, fresh.probe.re);
  assert.deepEqual(state.stages.gradient, fresh.stages.gradient);
  assert.equal(state.visited, 0);

  // And it runs again from there, to the same place.
  for (const position of order) epieStep(state, amplitudes, position, {beta: 0.9, probeBeta: 0.9});
  for (const position of order) epieStep(fresh, amplitudes, position, {beta: 0.9, probeBeta: 0.9});
  assert.deepEqual(state.phase, fresh.phase);
});

test("centring a spectrum undoes the shift the probe's position put on it", (t) => {
  // A rolled probe and an unrolled one differ by `exp(-2 pi i k . R)` and by
  // nothing else, so removing the ramp from the first has to give the second
  // back exactly. That pins the sign, which is the whole of what can go wrong.
  const probeReal = realProbe({C10: -30, C30: 1e4});
  const reference = {re: Float64Array.from(probeReal.re), im: Float64Array.from(probeReal.im)};
  fft2(reference, N, N);

  let worst = 0;
  for (const [row, col] of [[0, 0], [2, 0], [0, 6], [4, 10], [N - 3, N - 1]]) {
    const rolled = complex(N * N);
    for (let ix = 0; ix < N; ix++) {
      const sx = (((ix - row) % N) + N) % N;
      for (let iy = 0; iy < N; iy++) {
        const sy = (((iy - col) % N) + N) % N;
        rolled.re[ix * N + iy] = probeReal.re[sx * N + sy];
        rolled.im[ix * N + iy] = probeReal.im[sx * N + sy];
      }
    }
    fft2(rolled, N, N);
    const centred = centreSpectrum(rolled, N, row, col);
    for (let i = 0; i < N * N; i++) {
      worst = Math.max(
        worst,
        Math.abs(centred.re[i] - reference.re[i]),
        Math.abs(centred.im[i] - reference.im[i])
      );
    }
  }
  t.diagnostic(`max error ${worst.toExponential(2)}`);
  assert.ok(worst < 1e-14, `centring left ${worst} behind — check the sign of the ramp`);
});

test("the scan order is a permutation", () => {
  const order = scanOrder(M * M, mulberry32(2));
  assert.equal(order.length, M * M);
  const seen = new Set(order);
  assert.equal(seen.size, M * M, "the shuffle lost or repeated a position");
});

// ---------------------------------------------------------------------------
// One dimension: the step a tilt series takes at every angle.
// ---------------------------------------------------------------------------

/** A line of measurements from a known phase, at a given dose. */
function lineExperiment(n, truth, {width = 4, stride = 4, dose = Infinity, seed = 1} = {}) {
  const probe = complex(n);
  for (let i = 0; i < n; i++) {
    const s = i > n / 2 ? i - n : i;
    probe.re[i] = Math.exp(-(s * s) / (2 * width * width));
  }
  const positions = Int32Array.from({length: n / stride}, (_, i) => i * stride);
  const amplitudes = new Float64Array(positions.length * n);
  const buffer = complex(n);
  const random = mulberry32(seed);
  for (let p = 0; p < positions.length; p++) {
    const shift = positions[p];
    for (let i = 0; i < n; i++) {
      const s = ((i - shift) % n + n) % n;
      buffer.re[i] = Math.cos(truth[i]) * probe.re[s];
      buffer.im[i] = Math.sin(truth[i]) * probe.re[s];
    }
    fft2(buffer, n, 1);
    for (let i = 0; i < n; i++) {
      let intensity = buffer.re[i] ** 2 + buffer.im[i] ** 2;
      if (Number.isFinite(dose)) intensity = poisson(intensity * dose, random) / dose;
      amplitudes[p * n + i] = Math.sqrt(intensity);
    }
  }
  return {probe, positions, amplitudes};
}

function smoothLine(n, amplitude = 0.8) {
  const truth = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const x = (2 * Math.PI * i) / n;
    truth[i] = amplitude * (Math.sin(x) + 0.4 * Math.cos(3 * x + 0.7) + 0.25 * Math.sin(5 * x));
  }
  return truth;
}

test("one-dimensional ePIE recovers a line, up to a constant", () => {
  const n = 64;
  const truth = smoothLine(n);
  const {probe, positions, amplitudes} = lineExperiment(n, truth);

  const phase = new Float64Array(n);
  let residual = 1;
  for (let sweep = 0; sweep < 200; sweep++) {
    residual = epieLine(phase, probe, amplitudes, positions, {beta: 0.9});
  }

  const mean = (a) => a.reduce((sum, v) => sum + v, 0) / a.length;
  const offsetTruth = mean(truth);
  const offsetPhase = mean(phase);
  let worst = 0;
  for (let i = 0; i < n; i++) {
    worst = Math.max(worst, Math.abs((phase[i] - offsetPhase) - (truth[i] - offsetTruth)));
  }
  assert.ok(residual < 1e-3, `residual ${residual}`);
  assert.ok(worst < 0.02, `worst phase error ${worst} rad`);
});

test("the constant is the one thing it cannot recover", () => {
  // Two runs from different starting offsets land on the same shape at different
  // heights, because the measurement is blind to a constant added to the phase.
  // A tilt series reconstructed one projection at a time therefore comes back
  // with a different offset on every line, which is the artefact a joint
  // reconstruction exists to avoid.
  const n = 64;
  const truth = smoothLine(n);
  const {probe, positions, amplitudes} = lineExperiment(n, truth);

  const runs = [0, 1.3].map((start) => {
    const phase = new Float64Array(n).fill(start);
    for (let sweep = 0; sweep < 200; sweep++) {
      epieLine(phase, probe, amplitudes, positions, {beta: 0.9});
    }
    return phase;
  });

  const differences = Float64Array.from(runs[0], (v, i) => v - runs[1][i]);
  const mean = differences.reduce((sum, v) => sum + v, 0) / n;
  let spread = 0;
  for (const d of differences) spread = Math.max(spread, Math.abs(d - mean));
  assert.ok(Math.abs(mean) > 0.5, `the two runs should differ by a constant, got ${mean}`);
  assert.ok(spread < 0.02, `and by nothing else, got ${spread}`);
});

test("shot noise puts a floor under the residual", () => {
  const n = 64;
  const truth = smoothLine(n);
  const quiet = lineExperiment(n, truth, {dose: 1e6, seed: 4});
  const noisy = lineExperiment(n, truth, {dose: 1e3, seed: 4});

  const settle = (experiment) => {
    const phase = new Float64Array(n);
    let residual = 1;
    for (let sweep = 0; sweep < 150; sweep++) {
      residual = epieLine(phase, experiment.probe, experiment.amplitudes,
        experiment.positions, {beta: 0.9});
    }
    return residual;
  };
  assert.ok(settle(noisy) > 3 * settle(quiet), "a noisy measurement cannot be fitted as closely");
});
