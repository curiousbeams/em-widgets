import test from "node:test";
import assert from "node:assert/strict";

import {
  makeView, project, unproject, depthSort, depthRange,
  signedArea, emitSubpath, ellipsePoints,
  quatFromAxisAngle, quatMultiply, quatToMatrix, applyMatrix, quatNormalize, trackball,
  drawPlane, drawPlaneBanded, beamRadiusAt
} from "../kit/scene3d.js";

/** A canvas 2D context stub that records the calls we care about. */
function mockContext() {
  const calls = [];
  const ctx = {
    canvas: {width: 400, height: 400},
    transforms: [],
    images: [],
    save: () => calls.push(["save"]),
    restore: () => calls.push(["restore"]),
    beginPath: () => calls.push(["beginPath"]),
    closePath: () => calls.push(["closePath"]),
    moveTo: (x, y) => calls.push(["moveTo", x, y]),
    lineTo: (x, y) => calls.push(["lineTo", x, y]),
    arc: (...a) => calls.push(["arc", ...a]),
    fill: (rule) => calls.push(["fill", rule]),
    stroke: () => calls.push(["stroke"]),
    fillRect: () => calls.push(["fillRect"]),
    clearRect: () => calls.push(["clearRect"]),
    // The kit composes with `transform` so it never clobbers a canvas's
    // device-pixel scale. The mock's base transform is the identity, so a
    // composed matrix is simply the matrix — record both the same way.
    transform(...t) { ctx.transforms.push(t); },
    setTransform(...t) { ctx.transforms.push(t); },
    drawImage(...a) { ctx.images.push(a); },
    createRadialGradient: () => ({addColorStop() {}}),
    calls
  };
  return ctx;
}

const image = {width: 64, height: 64};

test("project and unproject round-trip on the sample plane", (t) => {
  const view = makeView(0, -26, 20, 200, 200, 120);
  let worst = 0;
  for (const x of [-5, -1, 0, 2.5, 7]) {
    for (const y of [-5, -1, 0, 2.5, 7]) {
      const p = project(x, y, 0, view);
      const back = unproject(p.sx, p.sy, view);
      worst = Math.max(worst, Math.hypot(back.x - x, back.y - y));
    }
  }
  t.diagnostic(`max round-trip error ${worst.toExponential(2)} scene units`);
  assert.ok(worst < 1e-9, `round-trip error ${worst}`);
});

test("project and unproject round-trip with a rotated camera", () => {
  const view = makeView(35, -18, 14, 100, 120, 90);
  const p = project(3, -2, 0, view);
  const back = unproject(p.sx, p.sy, view);
  assert.ok(Math.hypot(back.x - 3, back.y + 2) < 1e-9);
});

test("orthographic mode ignores depth", () => {
  const view = makeView(0, -26, 20, 0, 0, 0);
  const near = project(1, -10, 0, view);
  const far = project(1, 10, 0, view);
  // With no perspective, equal x maps to equal screen x regardless of depth.
  assert.ok(Math.abs(near.sx - far.sx) < 1e-12);
});

test("perspective enlarges points nearer the camera", () => {
  // `depth` increases TOWARD the camera, so +y is the near side here.
  const view = makeView(0, -26, 20, 0, 0, 120);
  const far = project(1, -10, 0, view);
  const near = project(1, 10, 0, view);
  assert.ok(near.depth > far.depth, "depth must increase toward the camera");
  assert.ok(Math.abs(near.sx) > Math.abs(far.sx), "the nearer point should project larger");
});

test("z maps to screen-y so scenes read as a column", () => {
  const view = makeView(0, -26, 20, 0, 0, 0);
  const top = project(0, 0, 10, view);
  const bottom = project(0, 0, -10, view);
  assert.ok(top.sy < bottom.sy, "higher z must be higher on screen");
});

test("depthSort orders back to front and depthRange normalises", () => {
  const items = [{depth: 5}, {depth: -3}, {depth: 1}];
  assert.deepEqual(depthSort(items).map((d) => d.depth), [-3, 1, 5]);
  const {t} = depthRange(items);
  assert.equal(t(-3), 0);
  assert.equal(t(5), 1);
});

test("signedArea detects winding, and emitSubpath normalises it", () => {
  const cw = [{sx: 0, sy: 0}, {sx: 1, sy: 0}, {sx: 1, sy: 1}, {sx: 0, sy: 1}];
  const ccw = cw.slice().reverse();
  assert.ok(signedArea(cw) * signedArea(ccw) < 0, "reversing must flip the sign");

  // Emitting both with the same reference sign must yield the same vertex order.
  const a = mockContext();
  const b = mockContext();
  emitSubpath(a, cw, 1);
  emitSubpath(b, ccw, 1);
  assert.deepEqual(a.calls, b.calls, "both windings should emit identically");
});

test("ellipsePoints closes a full circle of the requested size", () => {
  const view = makeView(0, -26, 20, 0, 0, 0);
  const pts = ellipsePoints(view, 0, 0, 0, 2, 32);
  assert.equal(pts.length, 32);
  const xs = pts.map((p) => p.sx);
  // A circle of radius 2 at zoom 20 spans 80 screen px horizontally.
  assert.ok(Math.abs((Math.max(...xs) - Math.min(...xs)) - 80) < 1e-9);
});

test("drawPlaneBanded agrees with drawPlane when the plane is face-on", () => {
  // Looking straight down the optic axis, every point of the z=0 plane has the
  // same depth, so there is no foreshortening and the banded transform must
  // reduce exactly to the single affine one.
  const view = makeView(0, -90, 20, 100, 100, 120);
  const flat = mockContext();
  const banded = mockContext();
  drawPlane(flat, view, image, 0, 0, 0, 3);
  drawPlaneBanded(banded, view, image, 0, 0, 0, 3, {bands: 8});

  const [a] = flat.transforms;
  for (const t of banded.transforms.slice(0, 8)) {
    for (let i = 0; i < 4; i++) {
      assert.ok(Math.abs(t[i] - a[i]) < 1e-9, `band basis ${i}: ${t[i]} vs ${a[i]}`);
    }
  }
});

test("drawPlaneBanded draws every band with a slight source overlap", () => {
  const view = makeView(0, -26, 20, 100, 100, 120);
  const ctx = mockContext();
  drawPlaneBanded(ctx, view, image, 0, 0, 0, 3, {bands: 8});
  assert.equal(ctx.images.length, 8);
  const heights = ctx.images.map((a) => a[4]); // source height argument
  assert.ok(heights.slice(0, -1).every((h) => h > 64 / 8), "inner bands overlap");
  assert.equal(heights.at(-1), 64 / 8, "the last band must not overrun the image");
});

test("quaternion helpers stay orthonormal under composition", () => {
  let q = [1, 0, 0, 0];
  for (let i = 0; i < 200; i++) q = trackball(q, 7, -3, 300);
  assert.ok(Math.abs(Math.hypot(...q) - 1) < 1e-12, "drift in quaternion norm");

  const m = quatToMatrix(q);
  for (const v of [[1, 0, 0], [0, 1, 0], [0, 0, 1]]) {
    const r = applyMatrix(m, v);
    assert.ok(Math.abs(Math.hypot(...r) - 1) < 1e-12, "rotation must preserve length");
  }
});

test("a quaternion rotation of 90 degrees about z maps x to y", () => {
  const q = quatNormalize(quatFromAxisAngle([0, 0, 1], Math.PI / 2));
  const r = applyMatrix(quatToMatrix(q), [1, 0, 0]);
  assert.ok(Math.hypot(r[0] - 0, r[1] - 1, r[2] - 0) < 1e-12, `got ${r}`);
});

test("quatMultiply composes in the expected order", () => {
  const a = quatFromAxisAngle([0, 0, 1], Math.PI / 2);
  const combined = quatMultiply(a, a);
  const r = applyMatrix(quatToMatrix(combined), [1, 0, 0]);
  assert.ok(Math.hypot(r[0] + 1, r[1], r[2]) < 1e-12, `two 90° turns should give -x, got ${r}`);
});

test("beam radius vanishes at the crossover and matches at the detector", () => {
  const geom = {crossoverZ: 4, zDetector: -26, detectorRadius: 8};
  assert.equal(beamRadiusAt(4, geom), 0);
  assert.ok(Math.abs(beamRadiusAt(-26, geom) - 8) < 1e-12);
  // Focused on the sample: zero footprint there, by construction.
  assert.equal(beamRadiusAt(0, {...geom, crossoverZ: 0}), 0);
  // Defocused: the footprint grows linearly with the crossover height.
  const a = beamRadiusAt(0, {...geom, crossoverZ: 2});
  const b = beamRadiusAt(0, {...geom, crossoverZ: 4});
  assert.ok(b > a && a > 0, `expected growth, got ${a} then ${b}`);
});
