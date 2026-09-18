// Iterative ptychography: the extended ptychographical iterative engine.
//
// Direct methods (see `ptycho.js`) invert one linear operator once. ePIE gives
// that up and solves the nonlinear problem by repetition: guess an object, work
// out what the detector would have seen, replace the modelled amplitude with the
// measured one, and push the difference back into the guess. Then do it again at
// the next scan position.
//
// The step at one position is
//
//     psi(r)      = O(r) P(r - R)                     the exit wave
//     Psi(k)      = F[psi]                            what the detector models
//     Psi'(k)     = sqrt(I(k)) Psi(k) / |Psi(k)|      the measured modulus, modelled phase
//     dpsi(r)     = F^-1[Psi' - Psi]                  the correction
//     O(r)       += beta  conj(P(r-R)) dpsi(r) / max|P|^2
//     P(r-R)     += beta' conj(O(r))   dpsi(r) / max|O|^2
//
// The third line is a projection onto the set of waves with the measured
// modulus — the same `Pi_f` the `projections.js` module and the projection-sets
// widget explore in two dimensions, here in `n^2` of them.
//
// The object can be solved for in two parameterisations. `"complex"` treats
// `O(r)` as an unconstrained complex field, which is the general case and what
// the equations above are written in. `"potential"` writes it as
// `O(r) = exp[i phi(r)]` and solves for the real phase instead,
//
//     phi(r)     += beta Im{ conj(O) conj(P(r-R)) dpsi(r) } / max|P|^2
//
// which is the same gradient projected onto the manifold of pure phase objects.
// It is worth having for three reasons: a specimen thin enough for this model is
// a phase object, the update is real and so can be drawn without domain
// colouring, and the phase is no longer read back through `atan2`, so an object
// stronger than pi radians does not wrap. `positivity` then clamps it at zero,
// which says the specimen may retard the wave but not advance it — and, as a
// side effect, fixes the constant that ptychography otherwise cannot determine.
//
// Conventions follow the rest of the kit: corner-centered grids, complex arrays
// as `{re, im}` pairs of Float64Array, row-major `[ix * ny + iy]`. The probe is
// real-space and centred at index (0,0), as `ifft2(complexProbe(...))` returns
// it, and it is moved by an integer index roll — exact, because the scan grid
// divides the object grid.

import {complex, fft2, ifft2} from "./fft.js";

/**
 * Undo the phase ramp that a real-space shift puts on a spectrum.
 *
 * The exit wave here is built on the whole grid with the probe rolled to the
 * scan position, so its transform carries `exp(-2 pi i k . R)`. That is
 * physically right, and it is what brings the correction back to the right place
 * when it is transformed home — but domain-coloured it is a rainbow that turns
 * over once per pixel of displacement, and the structure underneath it is
 * unreadable. A ptychography code never shows it, because it crops the object
 * around the probe and so works with an exit wave already at the origin.
 *
 * This is that crop's view, for a widget that does not crop: multiply the
 * ramp back out. Nothing the reconstruction does changes — call it on a copy,
 * on the way to a panel.
 *
 * The shift is a whole number of pixels, so the ramp is a table lookup rather
 * than a transcendental, and the result is exact to rounding.
 *
 * @param {{re: Float64Array, im: Float64Array}} spectrum corner-centered
 * @param {number} n grid size (square)
 * @param {number} row shift along the first axis, in pixels
 * @param {number} col shift along the second axis
 * @param {{re, im}} [out] defaults to a fresh array; may not alias `spectrum`
 */
export function centreSpectrum(spectrum, n, row, col, out = complex(n * n)) {
  const rampRe = new Float64Array(n);
  const rampIm = new Float64Array(n);
  for (let j = 0; j < n; j++) {
    rampRe[j] = Math.cos((2 * Math.PI * j) / n);
    rampIm[j] = Math.sin((2 * Math.PI * j) / n);
  }

  const half = n / 2;
  // The ramp separates along the two axes, so one length-n pass sets up the
  // columns and the rows are hoisted out of the inner loop.
  const colRe = new Float64Array(n);
  const colIm = new Float64Array(n);
  for (let b = 0; b < n; b++) {
    const fb = b < half ? b : b - n;
    const j = (((fb * col) % n) + n) % n;
    colRe[b] = rampRe[j];
    colIm[b] = rampIm[j];
  }

  for (let a = 0; a < n; a++) {
    const fa = a < half ? a : a - n;
    const j = (((fa * row) % n) + n) % n;
    const ar = rampRe[j];
    const ai = rampIm[j];
    const base = a * n;
    for (let b = 0; b < n; b++) {
      const tr = ar * colRe[b] - ai * colIm[b];
      const ti = ar * colIm[b] + ai * colRe[b];
      const i = base + b;
      const sr = spectrum.re[i];
      const si = spectrum.im[i];
      out.re[i] = sr * tr - si * ti;
      out.im[i] = sr * ti + si * tr;
    }
  }
  return out;
}

/** The largest `|z|^2` in a complex array. */
function peakIntensity({re, im}) {
  let peak = 0;
  for (let i = 0; i < re.length; i++) {
    const value = re[i] * re[i] + im[i] * im[i];
    if (value > peak) peak = value;
  }
  return peak;
}

/**
 * A reconstruction in progress.
 *
 * The object starts at one everywhere — vacuum, zero phase — which is the honest
 * starting guess when nothing is known, and which makes the first exit wave the
 * probe itself.
 *
 * The `stages` are the intermediate quantities of one step, kept as buffers
 * rather than returned fresh so that a widget can draw them without allocating
 * five `n x n` arrays thirty times a second. `gradient` is real and is filled
 * only in the `"potential"` parameterisation, where the object update is a real
 * field; `delta` is the complex correction wave, and is what a `"complex"`
 * object is updated with.
 *
 * @param {object} options
 * @param {number} options.n object and detector grid size (square)
 * @param {number} options.m scan grid size; must divide `n`
 * @param {{re: Float64Array, im: Float64Array}} options.probe real-space probe,
 *   centred at index (0,0)
 * @param {"complex"|"potential"} [options.objectType="complex"] see the module
 *   comment. `"potential"` also makes `phase` the primary array.
 * @param {boolean} [options.positivity=false] clamp the phase at zero. Only
 *   meaningful for a `"potential"` object.
 */
export function epieState({n, m, probe, objectType = "complex", positivity = false}) {
  if (n % m !== 0) {
    throw new Error(`scan grid m=${m} must divide the object grid n=${n}`);
  }
  if (objectType !== "complex" && objectType !== "potential") {
    throw new Error(`unknown object type "${objectType}"`);
  }
  const size = n * n;
  const object = complex(size);
  object.re.fill(1);

  return {
    n,
    m,
    stride: n / m,
    objectType,
    positivity,
    object,
    // Kept alongside `object` rather than derived from it: `atan2` wraps, and a
    // reconstruction is free to walk past pi.
    phase: objectType === "potential" ? new Float64Array(size) : null,
    probe: {re: Float64Array.from(probe.re), im: Float64Array.from(probe.im)},
    stages: {
      exit: complex(size),
      modelled: complex(size),
      correction: complex(size),
      delta: complex(size),
      gradient: objectType === "potential" ? new Float64Array(size) : null
    },
    // Where the shifted probe is read from, rebuilt each step. A buffer rather
    // than a fresh array because a step runs a few hundred times a second.
    shifted: new Int32Array(size),
    visited: 0,
    error: 0
  };
}

/**
 * Put a reconstruction back to where it started, keeping its buffers.
 *
 * The object goes back to vacuum and the stages to zero. The probe goes back to
 * `probe` if one is given, which is what a run that has been refining it needs;
 * left out, the probe is kept as it now stands.
 *
 * @param {object} state from {@link epieState}
 * @param {{probe?: {re: Float64Array, im: Float64Array}}} [options]
 */
export function epieReset(state, {probe = null} = {}) {
  state.object.re.fill(1);
  state.object.im.fill(0);
  if (state.phase) state.phase.fill(0);
  for (const stage of Object.values(state.stages)) {
    if (!stage) continue;
    if (stage.re) {
      stage.re.fill(0);
      stage.im.fill(0);
    } else {
      stage.fill(0);
    }
  }
  if (probe) {
    state.probe.re.set(probe.re);
    state.probe.im.set(probe.im);
  }
  state.visited = 0;
  state.error = 0;
  return state;
}

/**
 * One position's update, in place.
 *
 * `amplitudes` is `sqrt(I)` for every scan position and detector pixel, laid out
 * `[R][k]` — that is, `forwardModel({amplitude: true, layout: "rk"})`, so that
 * this position's pattern is one contiguous run.
 *
 * The object is updated first and the probe second, from the object as it now
 * stands. Doing both from the pre-update object differs at order `beta^2` and
 * needs a second copy of the grid; every implementation this is checked against
 * makes the same choice.
 *
 * @param {object} state from {@link epieState}
 * @param {Float32Array} amplitudes
 * @param {number} scanIndex which scan position, `0 .. m*m - 1`
 * @param {object} [options]
 * @param {number} [options.beta=0.9] object step size
 * @param {number} [options.probeBeta=0] probe step size; zero holds the probe at
 *   whatever it was initialised with, which is what to do when it is known
 * @returns {{error: number, row: number, col: number}} `error` is the amplitude
 *   residual at this position, relative to the measured amplitude
 */
export function epieStep(state, amplitudes, scanIndex, {beta = 0.9, probeBeta = 0} = {}) {
  const {n, m, stride, object, phase, probe, stages, shifted, positivity} = state;
  const size = n * n;
  const {exit, modelled, correction, delta, gradient} = stages;
  const potential = state.objectType === "potential";

  const row = ((scanIndex / m) | 0) * stride;
  const col = (scanIndex % m) * stride;
  const base = scanIndex * size;

  // The probe's index map for this position. Reading `P` at a rolled index is
  // the shift, exact for whole pixels — the same trick `ptycho-sim.js` uses to
  // simulate the data, so the model and the reconstruction agree by
  // construction rather than to within an interpolation.
  for (let ix = 0; ix < n; ix++) {
    const sx = (((ix - row) % n) + n) % n;
    for (let iy = 0; iy < n; iy++) {
      const sy = (((iy - col) % n) + n) % n;
      shifted[ix * n + iy] = sx * n + sy;
    }
  }

  for (let i = 0; i < size; i++) {
    const p = shifted[i];
    const or = object.re[i];
    const oi = object.im[i];
    const pr = probe.re[p];
    const pi = probe.im[p];
    exit.re[i] = or * pr - oi * pi;
    exit.im[i] = or * pi + oi * pr;
  }

  modelled.re.set(exit.re);
  modelled.im.set(exit.im);
  fft2(modelled, n, n);

  // The Fourier modulus projection, kept as the difference it makes rather than
  // as the projected wave. The projected wave's modulus is the measured one by
  // construction, so a panel of it would be a panel of the data; the difference
  // is the part that carries new information, and it shrinks as the
  // reconstruction converges.
  let residual = 0;
  let measured = 0;
  for (let i = 0; i < size; i++) {
    const re = modelled.re[i];
    const im = modelled.im[i];
    const magnitude = Math.hypot(re, im);
    const target = amplitudes[base + i];
    const difference = magnitude - target;
    residual += difference * difference;
    measured += target * target;
    if (magnitude > 0) {
      const scale = target / magnitude - 1;
      correction.re[i] = re * scale;
      correction.im[i] = im * scale;
    } else {
      // Nothing modelled here, so there is no phase to keep. The measured
      // modulus goes in along the real axis; any choice is as good, and this one
      // at least does not depend on floating-point noise in an empty pixel.
      correction.re[i] = target;
      correction.im[i] = 0;
    }
  }

  delta.re.set(correction.re);
  delta.im.set(correction.im);
  ifft2(delta, n, n);

  const probePeak = peakIntensity(probe) || 1;
  const gain = beta / probePeak;
  if (potential) {
    for (let i = 0; i < size; i++) {
      const p = shifted[i];
      const or = object.re[i];
      const oi = object.im[i];
      const pr = probe.re[p];
      const pi = probe.im[p];
      const dr = delta.re[i];
      const di = delta.im[i];
      // Im{ conj(O) conj(P) dpsi }, with conj(O) conj(P) = (cr + i ci).
      const cr = or * pr - oi * pi;
      const ci = -(or * pi + oi * pr);
      const step = (cr * di + ci * dr) / probePeak;
      // Stored before the step size and before the clamp, so a widget showing
      // this panel is showing the gradient rather than what survived of it.
      gradient[i] = step;
      let value = phase[i] + beta * step;
      if (positivity && value < 0) value = 0;
      phase[i] = value;
      // `object` is kept in step so that the next exit wave, and the probe
      // update below, need no second pass over the grid.
      object.re[i] = Math.cos(value);
      object.im[i] = Math.sin(value);
    }
  } else {
    for (let i = 0; i < size; i++) {
      const p = shifted[i];
      // conj(P) * dpsi
      const pr = probe.re[p];
      const pi = probe.im[p];
      const dr = delta.re[i];
      const di = delta.im[i];
      object.re[i] += gain * (pr * dr + pi * di);
      object.im[i] += gain * (pr * di - pi * dr);
    }
  }

  if (probeBeta > 0) {
    const objectPeak = peakIntensity(object) || 1;
    const probeGain = probeBeta / objectPeak;
    for (let i = 0; i < size; i++) {
      const p = shifted[i];
      const or = object.re[i];
      const oi = object.im[i];
      const dr = delta.re[i];
      const di = delta.im[i];
      probe.re[p] += probeGain * (or * dr + oi * di);
      probe.im[p] += probeGain * (or * di - oi * dr);
    }
  }

  const error = measured > 0 ? Math.sqrt(residual / measured) : 0;
  state.visited += 1;
  state.error = error;
  return {error, row, col};
}

/**
 * The object's phase, which is what a thin specimen is.
 *
 * A `"potential"` object already holds it, unwrapped. A `"complex"` one has to
 * be read back through `atan2`, which is only the phase modulo 2 pi.
 */
export function reconstructedPhase(state, out = new Float64Array(state.n * state.n)) {
  if (state.phase) {
    out.set(state.phase);
    return out;
  }
  const {re, im} = state.object;
  for (let i = 0; i < out.length; i++) out[i] = Math.atan2(im[i], re[i]);
  return out;
}

/**
 * A visiting order for the scan positions.
 *
 * Shuffled, and not as a detail. Sweeping the scan in raster order lets each
 * update undo most of the last one at the same edge, and the reconstruction
 * develops a stripe pattern along the scan direction. Visiting at random breaks
 * that correlation, and is why ePIE is normally described with a random order.
 *
 * @param {number} total how many positions
 * @param {() => number} [random=Math.random]
 */
export function scanOrder(total, random = Math.random) {
  const order = Int32Array.from({length: total}, (_, i) => i);
  for (let i = total - 1; i > 0; i--) {
    const j = Math.floor(random() * (i + 1));
    const swap = order[i];
    order[i] = order[j];
    order[j] = swap;
  }
  return order;
}

/**
 * The diameter enclosing a given fraction of the probe's intensity, in pixels.
 *
 * The one number that says whether a ptychographic scan can work: neighbouring
 * probes have to overlap, and how much they overlap is this against the scan
 * step. Measured from the probe rather than from a formula for it, so it stays
 * right when the probe is defocused, aberrated or apodised.
 *
 * @param {{re, im}} probeReal centred at index (0,0)
 * @param {number} n
 * @param {number} [fraction=0.9]
 */
export function probeDiameter(probeReal, n, fraction = 0.9) {
  const half = n / 2;
  const bins = new Float64Array(half + 1);
  let total = 0;
  for (let ix = 0; ix < n; ix++) {
    const fx = ix < half ? ix : ix - n;
    for (let iy = 0; iy < n; iy++) {
      const fy = iy < half ? iy : iy - n;
      const i = ix * n + iy;
      const intensity = probeReal.re[i] * probeReal.re[i] + probeReal.im[i] * probeReal.im[i];
      const bin = Math.min(half, Math.round(Math.hypot(fx, fy)));
      bins[bin] += intensity;
      total += intensity;
    }
  }
  if (total === 0) return 0;
  let running = 0;
  for (let b = 0; b <= half; b++) {
    running += bins[b];
    if (running >= fraction * total) return 2 * b;
  }
  return n;
}

/**
 * One sweep of ePIE over a one-dimensional object.
 *
 * The same step as {@link epieStep} in its `"potential"` parameterisation, with
 * one axis instead of two. It is here because a projection is a line: a
 * tomographic tilt series measures one of these per angle, and solving them is
 * what a ptychographic-tomographic reconstruction does inside its outer loop.
 *
 * The phase it returns is fixed only up to an additive constant. A measurement
 * of `|F[e^{i phi} P]|` cannot see a constant added to `phi`, because that
 * multiplies the exit wave by a constant phase factor and the detector records
 * modulus. Every position in the sweep has the same blind spot, so sweeping does
 * not remove it — which is exactly why a tilt series reconstructed one
 * projection at a time comes back with a different offset on every line.
 *
 * @param {Float64Array} phase the estimate, `n` long, updated in place
 * @param {{re: Float64Array, im: Float64Array}} probe real-space, centred at index 0
 * @param {Float64Array} amplitudes measured `sqrt(I)`, laid out `[position][k]`
 * @param {ArrayLike<number>} positions integer offsets, one per measurement
 * @param {object} [options]
 * @param {number} [options.beta=0.9] step size
 * @param {{re: Float64Array, im: Float64Array}} [scratch]
 * @returns {number} the amplitude residual, as a fraction of the measured total
 */
export function epieLine(phase, probe, amplitudes, positions, {beta = 0.9} = {},
                         scratch = complex(phase.length)) {
  const n = phase.length;
  const {re, im} = scratch;
  let probePeak = 0;
  for (let i = 0; i < n; i++) {
    probePeak = Math.max(probePeak, probe.re[i] * probe.re[i] + probe.im[i] * probe.im[i]);
  }
  if (probePeak === 0) return 0;

  let residual = 0;
  let measured = 0;

  for (let p = 0; p < positions.length; p++) {
    const shift = positions[p];
    const base = p * n;

    // The exit wave: the object times the probe, moved to this position.
    for (let i = 0; i < n; i++) {
      const s = ((i - shift) % n + n) % n;
      const cos = Math.cos(phase[i]);
      const sin = Math.sin(phase[i]);
      re[i] = cos * probe.re[s] - sin * probe.im[s];
      im[i] = cos * probe.im[s] + sin * probe.re[s];
    }
    // `fft2` over an n x 1 grid is a one-dimensional transform.
    fft2(scratch, n, 1);

    // Replace the modelled modulus with the measured one.
    for (let i = 0; i < n; i++) {
      const modulus = Math.hypot(re[i], im[i]);
      const target = amplitudes[base + i];
      residual += Math.abs(modulus - target);
      measured += target;
      const scale = modulus > 1e-12 ? target / modulus : 0;
      re[i] = re[i] * scale - re[i];
      im[i] = im[i] * scale - im[i];
    }
    ifft2(scratch, n, 1);

    // And push the correction into the phase.
    for (let i = 0; i < n; i++) {
      const s = ((i - shift) % n + n) % n;
      const cos = Math.cos(phase[i]);
      const sin = Math.sin(phase[i]);
      // Im{ conj(O P) dpsi } — the same gradient the 2D potential branch takes,
      // written out.
      const cr = cos * probe.re[s] - sin * probe.im[s];
      const ci = -(cos * probe.im[s] + sin * probe.re[s]);
      phase[i] += (beta * (cr * im[i] + ci * re[i])) / probePeak;
    }
  }
  return measured > 0 ? residual / measured : 0;
}
