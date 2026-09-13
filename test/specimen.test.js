import test from "node:test";
import assert from "node:assert/strict";
import {readFileSync} from "node:fs";
import {fileURLToPath} from "node:url";

import {
  LOBATO, projectedScatteringFactor, splatDensity, projectedPotential
} from "../kit/specimen.js";
import {integrateGradient} from "../kit/image.js";

const read = (name) =>
  JSON.parse(readFileSync(fileURLToPath(new URL(`./fixtures/${name}`, import.meta.url)), "utf8"));

function decode(b64) {
  const buffer = Buffer.from(b64, "base64");
  return new Float32Array(buffer.buffer, buffer.byteOffset, buffer.length / 4);
}

function compare(actual, expected) {
  let worst = 0, scale = 0, at = -1;
  for (let i = 0; i < expected.length; i++) {
    scale = Math.max(scale, Math.abs(expected[i]));
    const d = Math.abs(actual[i] - expected[i]);
    if (d > worst) { worst = d; at = i; }
  }
  return {worst, scale, relative: worst / (scale || 1), at};
}

test("the Lobato projected scattering factor matches abtem's values", () => {
  // Spot values printed from abtem's own `projected_scattering_factor`.
  assert.ok(Math.abs(projectedScatteringFactor(0, LOBATO.Au) - 506.344452) < 1e-4);
  assert.ok(Math.abs(projectedScatteringFactor(1, LOBATO.Au) - 146.768692) < 1e-4);
  assert.ok(Math.abs(projectedScatteringFactor(4, LOBATO.Au) - 57.911648) < 1e-4);
  // Monotonically decreasing, and positive, for both elements.
  for (const symbol of ["Au", "C"]) {
    let previous = Infinity;
    for (let k = 0; k <= 4; k += 0.25) {
      const f = projectedScatteringFactor(k * k, LOBATO[symbol]);
      assert.ok(f > 0, `${symbol} f(${k}) = ${f} should be positive`);
      assert.ok(f < previous, `${symbol} f should decrease with k`);
      previous = f;
    }
  }
});

test("splatDensity conserves total weight and wraps at the edges", () => {
  const grid = {gpts: [8, 6], sampling: [0.5, 0.5]};
  // A position outside the cell must wrap rather than be dropped.
  const d = splatDensity([[1.25, 0.75], [-0.3, 5.9]], grid);
  let sum = 0;
  for (const v of d) sum += v;
  assert.ok(Math.abs(sum - 2) < 1e-12, `total weight ${sum}, expected 2`);
});

test("projectedPotential matches abtem for sub-pixel atom positions", (t) => {
  const f = read("abtem-potential.json");
  const groups = [...new Set(f.symbols)].map((symbol) => ({
    symbol,
    positions: f.positions.filter((_, i) => f.symbols[i] === symbol).map((p) => [p[0], p[1]])
  }));
  const V = projectedPotential(groups, {gpts: f.gpts, sampling: f.sampling});
  const {worst, scale, relative} = compare(V, decode(f.V));
  t.diagnostic(`max |V| = ${scale.toFixed(1)} V·Å, max |diff| = ${worst.toExponential(2)} (relative ${relative.toExponential(2)})`);
  // The fixture is stored as float32, which sets the achievable agreement.
  assert.ok(relative < 1e-5, `relative error ${relative}`);
});

test("projectedPotential is rectangular-grid correct", () => {
  // The fixture is 60 x 50 on a 12 x 10 Å cell, so a transposed implementation
  // cannot pass the test above — but assert the shape explicitly too.
  const f = read("abtem-potential.json");
  const groups = [{symbol: "Au", positions: [[1.0, 2.0]]}];
  const V = projectedPotential(groups, {gpts: f.gpts, sampling: f.sampling});
  assert.equal(V.length, f.gpts[0] * f.gpts[1]);
  // The peak must land at (x/dx, y/dy) = (5, 10) in [ix * ny + iy] order.
  let best = -Infinity, at = -1;
  for (let i = 0; i < V.length; i++) if (V[i] > best) { best = V[i]; at = i; }
  assert.equal(Math.floor(at / f.gpts[1]), 5);
  assert.equal(at % f.gpts[1], 10);
});

test("the sinc correction is what makes the potential exact", () => {
  // Without abtem's sinc deconvolution the peaks come out systematically low —
  // about 1.7% at 0.2 Å sampling. This guards against quietly dropping it.
  const f = read("abtem-potential.json");
  const groups = [...new Set(f.symbols)].map((symbol) => ({
    symbol,
    positions: f.positions.filter((_, i) => f.symbols[i] === symbol).map((p) => [p[0], p[1]])
  }));
  const V = projectedPotential(groups, {gpts: f.gpts, sampling: f.sampling});
  const ref = decode(f.V);
  let peakJs = -Infinity, peakRef = -Infinity;
  for (let i = 0; i < V.length; i++) {
    peakJs = Math.max(peakJs, V[i]);
    peakRef = Math.max(peakRef, ref[i]);
  }
  assert.ok(Math.abs(peakJs / peakRef - 1) < 1e-4, `peak ratio ${peakJs / peakRef}`);
});

test("integrateGradient matches abtem's _integrate_gradient_2d", (t) => {
  const f = read("abtem-icom.json");
  const [nx, ny] = f.shape;
  const T = integrateGradient(decode(f.gx), decode(f.gy), nx, ny, f.sampling);
  const {worst, scale, relative} = compare(T, decode(f.T));
  t.diagnostic(`range 0..${scale.toFixed(4)}, max |diff| = ${worst.toExponential(2)} (relative ${relative.toExponential(2)})`);
  assert.ok(relative < 1e-6, `relative error ${relative}`);
});

test("integrateGradient inverts an analytic gradient", () => {
  // T = sin(2πx/L): its gradient integrates back to itself, up to the offset
  // that `integrateGradient` removes.
  const nx = 32, ny = 32, sx = 0.5, sy = 0.5;
  const T0 = new Float64Array(nx * ny);
  const g0 = new Float64Array(nx * ny);
  const g1 = new Float64Array(nx * ny);
  const L = nx * sx;
  for (let ix = 0; ix < nx; ix++) {
    for (let iy = 0; iy < ny; iy++) {
      const x = ix * sx;
      T0[ix * ny + iy] = Math.sin((2 * Math.PI * x) / L);
      g0[ix * ny + iy] = ((2 * Math.PI) / L) * Math.cos((2 * Math.PI * x) / L);
    }
  }
  const T = integrateGradient(g0, g1, nx, ny, [sx, sy]);
  let min = Infinity;
  for (const v of T0) min = Math.min(min, v);
  let worst = 0;
  for (let i = 0; i < T.length; i++) worst = Math.max(worst, Math.abs(T[i] - (T0[i] - min)));
  assert.ok(worst < 1e-10, `max |diff| ${worst}`);
});
