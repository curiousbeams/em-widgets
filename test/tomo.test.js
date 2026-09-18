import test from "node:test";
import assert from "node:assert/strict";

import {projectVolume, backProject, tiltSeries} from "../kit/tomo.js";
import {mulberry32} from "../kit/specimen.js";

const N = 48;

/** A disc of radius `r` centred in the slice, with a soft edge. */
function disc(n, r, value = 1) {
  const volume = new Float64Array(n * n);
  const centre = (n - 1) / 2;
  for (let ix = 0; ix < n; ix++) {
    for (let iy = 0; iy < n; iy++) {
      const d = Math.hypot(ix - centre, iy - centre);
      volume[ix * n + iy] = value * Math.min(1, Math.max(0, r + 0.5 - d));
    }
  }
  return volume;
}

test("back-projection is the exact transpose of projection", () => {
  // The identity every iterative reconstruction is built on: <A x, y> must equal
  // <x, Aᵀ y>. Nothing else catches an adjoint that is merely nearly right, and
  // a nearly-right adjoint converges to a plausible wrong answer.
  const random = mulberry32(3);
  const x = Float64Array.from({length: N * N}, () => random() - 0.5);
  for (const theta of [0, 0.3, Math.PI / 4, 1.9, -0.7]) {
    const y = Float64Array.from({length: N}, () => random() - 0.5);
    const Ax = projectVolume(x, N, theta);
    const Aty = backProject(y, N, theta, new Float64Array(N * N));

    let left = 0;
    for (let i = 0; i < N; i++) left += Ax[i] * y[i];
    let right = 0;
    for (let i = 0; i < N * N; i++) right += x[i] * Aty[i];

    const error = Math.abs(left - right) / Math.max(Math.abs(left), 1e-12);
    assert.ok(error < 1e-12, `θ = ${theta}: ${left} vs ${right} (${error})`);
  }
});

test("at zero the projection is a sum across each row", () => {
  const random = mulberry32(5);
  const volume = Float64Array.from({length: N * N}, () => random());
  const line = projectVolume(volume, N, 0);
  for (let ix = 2; ix < N - 2; ix++) {
    let sum = 0;
    for (let iy = 0; iy < N; iy++) sum += volume[ix * N + iy];
    // The rays sample between samples, so this is a quadrature error rather than
    // an identity; a couple of percent on a random field is what that costs.
    assert.ok(Math.abs(line[ix] - sum) / sum < 0.03, `row ${ix}: ${line[ix]} vs ${sum}`);
  }
});

test("a disc projects to the chord through it", () => {
  const r = 14;
  const volume = disc(N, r);
  const line = projectVolume(volume, N, 0.7, undefined, 4 * N);
  const centre = (N - 1) / 2;
  let worst = 0;
  for (let d = 0; d < N; d++) {
    const s = Math.abs(d - centre);
    if (s > r - 2) continue;                  // skip the soft edge
    const chord = 2 * Math.sqrt(r * r - s * s);
    worst = Math.max(worst, Math.abs(line[d] - chord) / chord);
  }
  assert.ok(worst < 0.03, `worst relative error ${worst}`);
});

test("a projection from behind is the same line, reversed", () => {
  const random = mulberry32(7);
  const volume = Float64Array.from({length: N * N}, () => random());
  const front = projectVolume(volume, N, 0.4);
  const behind = projectVolume(volume, N, 0.4 + Math.PI);
  let worst = 0;
  for (let d = 2; d < N - 2; d++) {
    worst = Math.max(worst, Math.abs(front[d] - behind[N - 1 - d]));
  }
  assert.ok(worst < 1e-9, `worst difference ${worst}`);
});

test("the tilt series covers half a turn, or a stated wedge of it", () => {
  const full = tiltSeries(6);
  assert.equal(full.length, 6);
  assert.ok(Math.abs(full[1] - full[0] - Math.PI / 6) < 1e-12);
  assert.ok(full[0] >= -Math.PI / 2 && full.at(-1) < Math.PI / 2);

  const wedge = tiltSeries(6, Math.PI / 3);
  assert.ok(Math.abs(wedge.at(-1) - wedge[0] - (Math.PI / 3) * (5 / 6)) < 1e-12);
});

test("back-projecting a full tilt series reproduces the object's shape", () => {
  // Unfiltered back-projection is blurred by 1/|q|, so the test is on the shape
  // rather than the values: the disc has to come back as the brightest thing in
  // the middle, with its edge in the right place.
  const volume = disc(N, 10);
  const angles = tiltSeries(60);
  const out = new Float64Array(N * N);
  const line = new Float64Array(N);
  for (const theta of angles) {
    projectVolume(volume, N, theta, line);
    backProject(line, N, theta, out, 1 / angles.length);
  }
  const centre = (N - 1) / 2;
  const at = (dx, dy) => out[Math.round(centre + dx) * N + Math.round(centre + dy)];
  assert.ok(at(0, 0) > at(9, 0), "the middle is brighter than the rim");
  assert.ok(at(9, 0) > at(20, 0), "the rim is brighter than the outside");
});
