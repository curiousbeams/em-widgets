import test from "node:test";
import assert from "node:assert/strict";
import {readFileSync} from "node:fs";
import {fileURLToPath} from "node:url";

import {electronWavelength, interactionParameter} from "../kit/units.js";
import {
  complexProbe, fresnelPropagator, multislice, multisliceSpectrum,
  antialiasAperture, bandlimit
} from "../kit/optics.js";
import {fft2, ifft2, abs, expi, fftshift} from "../kit/fft.js";
import {spatialFrequencies} from "../kit/grid.js";

const fixture = JSON.parse(
  readFileSync(fileURLToPath(new URL("./fixtures/abtem-multislice.json", import.meta.url)), "utf8")
);

function decode(b64) {
  const buffer = Buffer.from(b64, "base64");
  return new Float32Array(buffer.buffer, buffer.byteOffset, buffer.length / 4);
}

const {gpts, sampling, energy, semiangle, slice_thickness: dz, n_slices: nSlices} = fixture;
const [nx, ny] = gpts;
const wavelength = electronWavelength(energy);

/**
 * Compare two complex fields up to a global complex scale.
 *
 * The kit normalises the probe in Fourier space and abtem in real space, so the
 * two differ by a constant factor that says nothing about the physics. Solving
 * for the best-fit scale removes it and leaves only the shape.
 */
function relativeError(actual, expected) {
  let nr = 0, ni = 0, dd = 0;
  for (let i = 0; i < expected.re.length; i++) {
    // <expected, actual> / <actual, actual>
    nr += expected.re[i] * actual.re[i] + expected.im[i] * actual.im[i];
    ni += expected.im[i] * actual.re[i] - expected.re[i] * actual.im[i];
    dd += actual.re[i] * actual.re[i] + actual.im[i] * actual.im[i];
  }
  const sr = nr / dd, si = ni / dd;
  let worst = 0, scale = 0;
  for (let i = 0; i < expected.re.length; i++) {
    const re = sr * actual.re[i] - si * actual.im[i];
    const im = sr * actual.im[i] + si * actual.re[i];
    worst = Math.max(worst, Math.hypot(re - expected.re[i], im - expected.im[i]));
    scale = Math.max(scale, Math.hypot(expected.re[i], expected.im[i]));
  }
  return worst / (scale || 1);
}

function realRelativeError(actual, expected) {
  let na = 0, ne = 0;
  for (let i = 0; i < expected.length; i++) { na += actual[i]; ne += expected[i]; }
  const s = ne / (na || 1);
  let worst = 0, scale = 0;
  for (let i = 0; i < expected.length; i++) {
    worst = Math.max(worst, Math.abs(s * actual[i] - expected[i]));
    scale = Math.max(scale, Math.abs(expected[i]));
  }
  return worst / (scale || 1);
}

/** The transmission functions, built from abtem's own potential slices. */
function transmissions() {
  const V = decode(fixture.V);
  const sigma = interactionParameter(energy);
  const per = nx * ny;
  return Array.from({length: nSlices}, (_, s) => {
    const phase = new Float64Array(per);
    for (let i = 0; i < per; i++) phase[i] = V[s * per + i] * sigma;
    return bandlimit(expi(phase), gpts, sampling);
  });
}

/** The kit's probe, shifted to `position` by a Fourier-space phase ramp. */
function entranceWave(defocus) {
  const P = complexProbe({
    gpts, sampling, energy, semiangleCutoff: semiangle,
    aberrations: {C10: -defocus}
  });
  const {kx, ky} = spatialFrequencies(gpts, sampling);
  const [x0, y0] = fixture.position;
  const shifted = {re: new Float64Array(nx * ny), im: new Float64Array(nx * ny)};
  for (let i = 0; i < kx.length; i++) {
    const t = -2 * Math.PI * (kx[i] * x0 + ky[i] * y0);
    const c = Math.cos(t), s = Math.sin(t);
    shifted.re[i] = P.re[i] * c - P.im[i] * s;
    shifted.im[i] = P.re[i] * s + P.im[i] * c;
  }
  return ifft2(shifted, nx, ny);
}

test("the fixture is a genuinely asymmetric, multi-slice problem", (t) => {
  t.diagnostic(`${nx}x${ny} grid, ${nSlices} slices of ${dz} A, probe at ${fixture.position}`);
  assert.ok(nx !== ny, "a square grid would hide a transpose");
  assert.ok(nSlices > 1, "a single slice would hide the propagator");
});


test("the antialias aperture matches abtem's", (t) => {
  // A 2/3 disc measured against the coarser sampling, with a cosine taper.
  const a = antialiasAperture(gpts, sampling);
  const inside = a.filter((v) => v === 1).length / a.length;
  t.diagnostic(`${(inside * 100).toFixed(1)}% of the grid passes unattenuated`);
  assert.ok(Math.abs(inside - 0.328) < 0.01, `passband fraction ${inside}`);
  for (const v of a) assert.ok(v >= 0 && v <= 1);
});

test("the propagator advances the wave, and reversing dz undoes it", () => {
  // The sign is pinned by the exit-wave test below; this guards the inverse.
  const forward = fresnelPropagator(gpts, sampling, energy, dz, {antialias: false});
  const back = fresnelPropagator(gpts, sampling, energy, -dz, {antialias: false});
  for (let i = 0; i < forward.re.length; i++) {
    const re = forward.re[i] * back.re[i] - forward.im[i] * back.im[i];
    const im = forward.re[i] * back.im[i] + forward.im[i] * back.re[i];
    assert.ok(Math.hypot(re - 1, im) < 1e-9, `round trip at ${i}: ${re}, ${im}`);
  }
});

for (const [label, expected] of Object.entries(fixture.cases)) {
  const defocus = expected.defocus;

  test(`the entrance wave matches abtem (${label}, defocus ${defocus} A)`, (t) => {
    const psi = entranceWave(defocus);
    const ref = {re: decode(expected.entrance_re), im: decode(expected.entrance_im)};
    const error = relativeError(psi, ref);
    t.diagnostic(`relative error ${error.toExponential(2)}`);
    assert.ok(error < 2e-5, `entrance wave differs by ${error}`);
  });

  test(`the exit wave matches abtem (${label}, defocus ${defocus} A)`, (t) => {
    const psi = multislice(entranceWave(defocus), transmissions(),
                           fresnelPropagator(gpts, sampling, energy, dz), nx, ny);
    const ref = {re: decode(expected.exit_re), im: decode(expected.exit_im)};
    const error = relativeError(psi, ref);
    t.diagnostic(`relative error ${error.toExponential(2)}`);
    assert.ok(error < 2e-4, `exit wave differs by ${error}`);
  });

  test(`the diffraction pattern matches abtem (${label}, defocus ${defocus} A)`, (t) => {
    const psi = multislice(entranceWave(defocus), transmissions(),
                           fresnelPropagator(gpts, sampling, energy, dz), nx, ny);
    const amplitude = abs(fft2({re: psi.re, im: psi.im}, nx, ny));
    // abtem returns diffraction patterns centred; ours come out of `fft2`
    // corner-centred, which is where the detector masks want them.
    const I = fftshift(Float64Array.from(amplitude, (a) => a * a), nx, ny);
    const error = realRelativeError(I, decode(expected.dp));
    t.diagnostic(`relative error ${error.toExponential(2)}`);
    assert.ok(error < 2e-4, `diffraction pattern differs by ${error}`);
  });
}

test("multisliceSpectrum is the transform of the exit wave", (t) => {
  // Equal up to the final propagator, which is unit-magnitude inside the
  // antialiasing aperture — so the intensities must agree exactly there.
  const expected = fixture.cases.focused;
  const T = transmissions();
  const P = fresnelPropagator(gpts, sampling, energy, dz);
  const spectrum = multisliceSpectrum(entranceWave(expected.defocus), T, P, nx, ny);
  const viaExit = fft2(
    (() => {
      const psi = multislice(entranceWave(expected.defocus), T, P, nx, ny);
      return {re: psi.re, im: psi.im};
    })(),
    nx, ny
  );

  const aperture = antialiasAperture(gpts, sampling);
  let worst = 0, scale = 0;
  for (let i = 0; i < spectrum.re.length; i++) {
    if (aperture[i] < 1) continue;   // the propagator's own aperture lives here
    const a = spectrum.re[i] ** 2 + spectrum.im[i] ** 2;
    const b = viaExit.re[i] ** 2 + viaExit.im[i] ** 2;
    worst = Math.max(worst, Math.abs(a - b));
    scale = Math.max(scale, b);
  }
  t.diagnostic(`max |dI| = ${(worst / scale).toExponential(2)} of the peak`);
  assert.ok(worst / scale < 1e-12, `intensities differ by ${worst / scale}`);
});

test("the probe's waist sits -C10 downstream of the plane it was built on", (t) => {
  // The sign that decides which way defocus moves the cross-over. It is not
  // visible in a single pattern — a probe defocused the wrong way looks exactly
  // like one defocused the right way — so it has to be measured, by propagating
  // through vacuum and finding the narrowest plane.
  const n = 128, dx = 0.25;
  const grid = [n, n], spacing = [dx, dx];
  const width = (P, z) => {
    const prop = fresnelPropagator(grid, spacing, energy, z, {antialias: false});
    const shifted = {re: new Float64Array(n * n), im: new Float64Array(n * n)};
    for (let i = 0; i < n * n; i++) {
      shifted.re[i] = P.re[i] * prop.re[i] - P.im[i] * prop.im[i];
      shifted.im[i] = P.re[i] * prop.im[i] + P.im[i] * prop.re[i];
    }
    const psi = ifft2(shifted, n, n);
    let m = 0, total = 0;
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        const I = psi.re[i * n + j] ** 2 + psi.im[i * n + j] ** 2;
        const di = Math.min(i, n - i) * dx, dj = Math.min(j, n - j) * dx;
        m += I * (di * di + dj * dj);
        total += I;
      }
    }
    return Math.sqrt(m / total);
  };

  for (const C10 of [-24, 0, 24]) {
    const P = complexProbe({
      gpts: grid, sampling: spacing, energy,
      semiangleCutoff: semiangle, aberrations: {C10}
    });
    let best = Infinity, at = null;
    for (let z = -40; z <= 40; z += 1) {
      const w = width(P, z);
      if (w < best) { best = w; at = z; }
    }
    t.diagnostic(`C10 = ${C10} Å -> waist ${at} Å downstream`);
    assert.ok(Math.abs(at - -C10) <= 1, `waist at ${at}, expected ${-C10}`);
  }
});
