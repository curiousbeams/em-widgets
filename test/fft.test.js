import test from "node:test";
import assert from "node:assert/strict";
import {readFileSync} from "node:fs";
import {fileURLToPath} from "node:url";

import {fft2, ifft2, fftshift, ifftshift, fft1d, ifft1d} from "../kit/fft.js";

const fixtures = JSON.parse(
  readFileSync(fileURLToPath(new URL("./fixtures/numerics.json", import.meta.url)), "utf8")
);

function assertClose(actual, expected, tol, what) {
  assert.equal(actual.length, expected.length, `${what}: length`);
  let worst = 0, at = -1;
  for (let i = 0; i < expected.length; i++) {
    const d = Math.abs(actual[i] - expected[i]);
    if (d > worst) { worst = d; at = i; }
  }
  assert.ok(
    worst <= tol,
    `${what}: max |diff| = ${worst.toExponential(3)} at index ${at} ` +
      `(got ${actual[at]}, want ${expected[at]}), tol ${tol}`
  );
  return worst;
}

test("fft2 matches numpy.fft.fft2", async (t) => {
  for (const [name, f] of Object.entries(fixtures.fft2)) {
    const [nx, ny] = f.shape;
    const a = {re: Float64Array.from(f.re), im: Float64Array.from(f.im)};
    fft2(a, nx, ny);
    const wr = assertClose(a.re, f.Fre, 1e-9, `${name} real`);
    const wi = assertClose(a.im, f.Fim, 1e-9, `${name} imag`);
    t.diagnostic(`${name} (${nx}x${ny}): max err ${Math.max(wr, wi).toExponential(2)}`);
  }
});

test("ifft2 inverts fft2", () => {
  for (const [, f] of Object.entries(fixtures.fft2)) {
    const [nx, ny] = f.shape;
    const a = {re: Float64Array.from(f.re), im: Float64Array.from(f.im)};
    ifft2(fft2({re: Float64Array.from(a.re), im: Float64Array.from(a.im)}, nx, ny), nx, ny);
  }
  // explicit round trip on a non-power-of-two grid
  const nx = 6, ny = 10;
  const re = Float64Array.from(fixtures.fft2.nonpow2.re);
  const im = Float64Array.from(fixtures.fft2.nonpow2.im);
  const a = {re: Float64Array.from(re), im: Float64Array.from(im)};
  ifft2(fft2(a, nx, ny), nx, ny);
  assertClose(a.re, [...re], 1e-10, "roundtrip real");
  assertClose(a.im, [...im], 1e-10, "roundtrip imag");
});

test("fft1d handles prime lengths via Bluestein", () => {
  const n = 13;
  const re = new Float64Array(n), im = new Float64Array(n);
  for (let i = 0; i < n; i++) { re[i] = Math.sin(i * 1.7); im[i] = Math.cos(i * 0.3); }
  const re0 = Float64Array.from(re), im0 = Float64Array.from(im);

  // brute-force DFT reference
  const dr = new Float64Array(n), di = new Float64Array(n);
  for (let k = 0; k < n; k++) {
    for (let j = 0; j < n; j++) {
      const a = (-2 * Math.PI * k * j) / n;
      dr[k] += re0[j] * Math.cos(a) - im0[j] * Math.sin(a);
      di[k] += re0[j] * Math.sin(a) + im0[j] * Math.cos(a);
    }
  }
  fft1d(re, im);
  assertClose(re, [...dr], 1e-10, "prime real");
  assertClose(im, [...di], 1e-10, "prime imag");

  ifft1d(re, im);
  assertClose(re, [...re0], 1e-10, "prime roundtrip real");
  assertClose(im, [...im0], 1e-10, "prime roundtrip imag");
});

test("fftshift / ifftshift round-trip on even and odd grids", () => {
  for (const [nx, ny] of [[4, 6], [5, 7], [4, 7]]) {
    const a = Float64Array.from({length: nx * ny}, (_, i) => i);
    const back = ifftshift(fftshift(a, nx, ny), nx, ny);
    assertClose(back, [...a], 0, `shift ${nx}x${ny}`);
  }
});
