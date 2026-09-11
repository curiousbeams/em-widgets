import test from "node:test";
import assert from "node:assert/strict";
import {readFileSync} from "node:fs";
import {fileURLToPath} from "node:url";

import {electronWavelength, interactionParameter} from "../kit/units.js";
import {spatialFrequencies, polarCoordinates, angularSpatialFrequencies, fftfreq} from "../kit/grid.js";
import {chi, softAperture, hardAperture, scherzerDefocus} from "../kit/optics.js";
import {gradient2d} from "../kit/image.js";

const num = JSON.parse(
  readFileSync(fileURLToPath(new URL("./fixtures/numerics.json", import.meta.url)), "utf8")
);

function maxDiff(a, b) {
  let worst = 0, at = -1;
  for (let i = 0; i < b.length; i++) {
    const d = Math.abs(a[i] - b[i]);
    if (d > worst) { worst = d; at = i; }
  }
  return {worst, at};
}

test("electronWavelength matches ctf/utils.py", (t) => {
  let worst = 0;
  for (const [E, want] of num.wavelength) {
    worst = Math.max(worst, Math.abs(electronWavelength(E) - want));
  }
  t.diagnostic(`max |diff| = ${worst.toExponential(2)} Angstrom`);
  assert.ok(worst < 1e-15);
});

test("interactionParameter is consistent with the wavelength", () => {
  // sanity: sigma decreases with energy, and 300 kV is ~0.00065 1/(V A)
  assert.ok(interactionParameter(300e3) < interactionParameter(80e3));
  assert.ok(Math.abs(interactionParameter(300e3) - 6.5e-4) < 1e-4);
});

test("fftfreq matches numpy for even and odd n", () => {
  assert.deepEqual([...fftfreq(4, 1)], [0, 0.25, -0.5, -0.25]);
  assert.deepEqual([...fftfreq(5, 1)].map((v) => +v.toFixed(10)), [0, 0.2, 0.4, -0.4, -0.2]);
});

test("spatialFrequencies and polarCoordinates match aberration_utils.py", (t) => {
  const g = num.grid;
  const {kx, ky} = spatialFrequencies(g.gpts, g.sampling);
  const {k, phi} = polarCoordinates(kx, ky);
  // The Python casts to float32, so compare at float32 precision.
  for (const [name, got, want] of [["kx", kx, g.kx], ["ky", ky, g.ky], ["k", k, g.k], ["phi", phi, g.phi]]) {
    const {worst, at} = maxDiff(got, want);
    t.diagnostic(`${name}: max |diff| = ${worst.toExponential(2)}`);
    assert.ok(worst < 1e-6, `${name}: ${worst} at ${at}`);
  }
});

test("chi matches aberration_utils.py:aberration_surface to 6th order", (t) => {
  const f = num.chi;
  const {alpha, phi} = angularSpatialFrequencies(f.gpts, f.sampling, f.wavelength);
  const got = chi(alpha, phi, f.wavelength, f.coefs);
  const {worst, at} = maxDiff(got, f.values);
  const scale = Math.max(...f.values.map(Math.abs));
  t.diagnostic(`max |diff| = ${worst.toExponential(2)} on values up to ${scale.toExponential(2)} (relative ${(worst / scale).toExponential(2)})`);
  // float32 inputs in the Python reference limit the achievable agreement
  assert.ok(worst / scale < 1e-6, `relative ${worst / scale} at ${at}`);
});

test("apertures match aberration_utils.py", (t) => {
  const f = num.aperture;
  const g = num.chi;
  const {alpha, phi} = angularSpatialFrequencies(g.gpts, g.sampling, g.wavelength);
  const soft = softAperture(alpha, phi, f.semiangle, f.angular_sampling);
  const hard = hardAperture(alpha, f.semiangle);
  const ds = maxDiff(soft, f.soft);
  const dh = maxDiff(hard, f.hard);
  t.diagnostic(`soft: ${ds.worst.toExponential(2)}, hard: ${dh.worst.toExponential(2)}`);
  assert.ok(ds.worst < 1e-6, `soft ${ds.worst} at ${ds.at}`);
  assert.equal(dh.worst, 0, `hard ${dh.worst} at ${dh.at}`);
});

test("scherzerDefocus has the expected sign and magnitude", () => {
  const lambda = electronWavelength(300e3);
  const Cs = 1e4;
  const df = scherzerDefocus(Cs, lambda);
  assert.ok(df < 0, "Scherzer defocus is underfocus for positive Cs");
  assert.ok(Math.abs(df + Math.sqrt(1.5 * Cs * lambda)) < 1e-12);
});

test("gradient2d reproduces an analytic gradient", () => {
  // f(x, y) = 3x + 5y on a 6x7 grid: gradient is exact for both interior and edges
  const nx = 6, ny = 7, dx = 0.3, dy = 0.7;
  const f = new Float64Array(nx * ny);
  for (let ix = 0; ix < nx; ix++)
    for (let iy = 0; iy < ny; iy++) f[ix * ny + iy] = 3 * (ix * dx) + 5 * (iy * dy);

  const {dRow, dCol} = gradient2d(f, nx, ny, dx, dy);
  for (let i = 0; i < f.length; i++) {
    assert.ok(Math.abs(dRow[i] - 3) < 1e-10, `dRow[${i}] = ${dRow[i]}`);
    assert.ok(Math.abs(dCol[i] - 5) < 1e-10, `dCol[${i}] = ${dCol[i]}`);
  }
});
