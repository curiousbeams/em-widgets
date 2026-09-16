import test from "node:test";
import assert from "node:assert/strict";

import {
  generalisedProjection, productProjection, centroid, namedParameters,
  line, polarCurve, residual, iteration
} from "../kit/projections.js";

const norm = ([x, y]) => Math.hypot(x, y);

// The two lines the notebook uses: y = x/2 and y = x/8, crossing at the origin.
const RED = {direction: [2, 1], extent: 1.5};
const BLUE = {direction: [8, 1], extent: 0.4};

// The same lines with a stretch cut out of each, which is what makes them
// non-convex and the algorithms differ.
const RED_GAPPED = {...RED, excluded: [[1 / 8, 1 / 4]]};
const BLUE_GAPPED = {...BLUE, excluded: [[1 / 32, 1 / 16]]};

// ---------------------------------------------------------------------------
// Sets
// ---------------------------------------------------------------------------

test("projecting onto a line lands on it, at right angles", () => {
  const set = line(RED);
  for (const point of [[1, 0], [0, 1], [-2, 3], [0.3, -0.7]]) {
    const [px, py] = set.project(point);
    // On the line: the projection is a multiple of the direction.
    assert.ok(Math.abs(px * 1 - py * 2) < 1e-12, `(${px}, ${py}) is off the line`);
    // At right angles: the residual is perpendicular to it.
    const dot = (point[0] - px) * 2 + (point[1] - py) * 1;
    assert.ok(Math.abs(dot) < 1e-12, `not perpendicular, dot = ${dot}`);
  }
});

test("a point already on a set stays put", () => {
  for (const set of [line(RED), line(RED_GAPPED), polarCurve({samples: 512})]) {
    const seed = set.kind === "curve" ? set.points[37] : [0.6, 0.3];
    const once = set.project(seed);
    const twice = set.project(once);
    assert.ok(norm([twice[0] - once[0], twice[1] - once[1]]) < 1e-12,
      `${set.name} is not idempotent`);
  }
});

test("a gap in a line pushes points to whichever end is nearer", () => {
  const set = line(RED_GAPPED);
  const at = (t) => [t * 2, t * 1];
  // Parameter 3/16 is the middle of the excluded stretch; either side of it
  // should snap to the opposite end.
  const below = set.project(at(0.18));
  const above = set.project(at(0.2));
  assert.ok(Math.abs(below[0] - 2 * (1 / 8)) < 1e-12, `expected the near end, got ${below}`);
  assert.ok(Math.abs(above[0] - 2 * (1 / 4)) < 1e-12, `expected the far end, got ${above}`);

  // And outside the gap nothing changes.
  const outside = set.project(at(0.5));
  assert.ok(Math.abs(outside[0] - 1) < 1e-12);
});

test("the drawn path of a gapped line has a hole in it", () => {
  const whole = line(RED);
  const gapped = line(RED_GAPPED);
  assert.equal(whole.path.length, 1, "an unbroken line is one stretch");
  assert.equal(gapped.path.length, 2, "a gapped one is two");
  // The hole is where the excluded interval says it is.
  assert.ok(Math.abs(gapped.path[0][1][0] - 2 * (1 / 8)) < 1e-12);
  assert.ok(Math.abs(gapped.path[1][0][0] - 2 * (1 / 4)) < 1e-12);
});

test("the polar curve closes and its projection is one of its samples", () => {
  const set = polarCurve({petals: 4, amplitude: 0.25, samples: 256, scale: 0.5});
  const ring = set.path[0];
  assert.equal(ring.length, 257, "the drawn path repeats its first point to close");
  assert.deepEqual(ring[0], ring.at(-1));

  const landed = set.project([3, -2]);
  assert.ok(set.points.some(([x, y]) => x === landed[0] && y === landed[1]),
    "the projection should be a sample of the curve");
});

// ---------------------------------------------------------------------------
// The iteration
// ---------------------------------------------------------------------------

test("alternating projections finds where two lines cross", (t) => {
  // The easy case, and the one everything else is measured against: two convex
  // sets, so every algorithm converges to the intersection.
  const sets = [line(RED), line(BLUE)];
  const {start, step} = iteration(sets, namedParameters().AP);
  let state = start([1.2, 0.9]);
  for (let i = 0; i < 400; i++) state = step(state);

  t.diagnostic(`after 400 steps: (${state.point.map((v) => v.toExponential(1))})`);
  assert.ok(norm(state.point) < 1e-6, `expected the origin, got ${state.point}`);
  assert.ok(residual(state.point, sets) < 1e-6);
});

test("every named algorithm solves the convex case", (t) => {
  const sets = [line(RED), line(BLUE)];
  for (const [name, parameters] of Object.entries(namedParameters())) {
    const {start, step} = iteration(sets, parameters);
    let state = start([1.2, 0.9]);
    let best = Infinity;
    for (let i = 0; i < 2000; i++) {
      state = step(state);
      best = Math.min(best, residual(state.point, sets));
    }
    t.diagnostic(`${name}: best residual ${best.toExponential(2)}`);
    assert.ok(best < 1e-4, `${name} stalled at ${best}`);
  }
});

test("alternating projections stalls in a non-convex trap where the difference map does not", (t) => {
  // The claim the widget makes, and the reason anyone uses anything other than
  // AP. Both lines have a stretch removed. Started so that its first few steps
  // carry it into the gap, AP settles onto a pair of gap ends and cannot leave;
  // DM and RRR keep moving and reach the intersection at the origin.
  const sets = [line(RED_GAPPED), line(BLUE_GAPPED)];
  const seed = [0.42, 0.2];

  const run = (parameters, steps) => {
    const {start, step} = iteration(sets, parameters);
    let state = start(seed);
    let best = Infinity;
    const tail = [];
    for (let i = 0; i < steps; i++) {
      state = step(state);
      const r = residual(state.point, sets);
      best = Math.min(best, r);
      if (i >= steps - 40) tail.push(r);
    }
    return {best, tail};
  };

  const named = namedParameters();
  const ap = run(named.AP, 2000);
  const dm = run(named.DM, 2000);
  const rrr = run(named.RRR, 2000);

  t.diagnostic(`AP  best residual ${ap.best.toExponential(2)}`);
  t.diagnostic(`DM  best residual ${dm.best.toExponential(2)}`);
  t.diagnostic(`RRR best residual ${rrr.best.toExponential(2)}`);

  assert.ok(ap.best > 1e-3, `AP should not solve this, but reached ${ap.best}`);
  // And it is stuck, not merely slow: its last forty residuals do not improve.
  const settled = Math.max(...ap.tail) - Math.min(...ap.tail);
  assert.ok(settled < 1e-9, `AP is still moving, spread ${settled}`);

  assert.ok(dm.best < 1e-6, `DM should escape, but reached ${dm.best}`);
  assert.ok(rrr.best < 1e-6, `RRR should escape, but reached ${rrr.best}`);
});

// ---------------------------------------------------------------------------
// More than two sets
// ---------------------------------------------------------------------------

test("the product-space lift keeps its copies together once they agree", (t) => {
  // The lift is not the two-set iteration in disguise. Its second projection is
  // onto the diagonal — the subspace where all the copies are equal — so
  // averaging plays the part that the second set plays with two sets. What it
  // does guarantee is that copies which agree, and agree with all the sets,
  // stay put.
  const set = line(RED);
  const onSet = set.project([0.9, -0.4]);
  for (const parameters of Object.values(namedParameters())) {
    const lifted = productProjection([set, set, set], parameters);
    const {next} = lifted([onSet, onSet, onSet]);
    for (const copy of next) {
      assert.ok(Math.hypot(copy[0] - onSet[0], copy[1] - onSet[1]) < 1e-12,
        `moved off a point that was already a solution: ${copy}`);
    }
  }

  // And three copies of one set solve that set, since its intersection with
  // itself is itself.
  const {start, step} = iteration([set, set, set], namedParameters().DM);
  let state = start([1.3, -0.8]);
  for (let i = 0; i < 500; i++) state = step(state);
  t.diagnostic(`three copies of one line: residual ${residual(state.point, [set]).toExponential(2)}`);
  assert.ok(residual(state.point, [set]) < 1e-9);
});

test("centroid averages, and the diagonal is a fixed point of it", () => {
  assert.deepEqual(centroid([[0, 0], [2, 4], [4, 2]]), [2, 2]);
  assert.deepEqual(centroid([[1, -1], [1, -1]]), [1, -1]);
});

test("three sets converge to a point on all three", (t) => {
  // Two gapped lines and a closed curve. The curve is sampled, so the residual
  // cannot fall below the spacing between its samples; the test asks for a
  // point on all three to within that, not to machine precision.
  const curve = polarCurve({
    petals: 4, amplitude: 0.25, samples: 2048,
    centre: [0.42, 0.2], scale: 0.42
  });
  const sets = [line(RED_GAPPED), line(BLUE_GAPPED), curve];
  const {start, step} = iteration(sets, namedParameters().RRR);

  let state = start([0.5, 0.25]);
  let best = Infinity;
  for (let i = 0; i < 4000; i++) {
    state = step(state);
    best = Math.min(best, residual(state.point, sets));
  }
  const spacing = Math.hypot(
    curve.points[1][0] - curve.points[0][0],
    curve.points[1][1] - curve.points[0][1]
  );
  t.diagnostic(`best residual ${best.toExponential(2)}, curve sample spacing ${spacing.toExponential(2)}`);
  // A few sample spacings, not machine precision: a sampled curve cannot place
  // a point more accurately than the gap between its samples.
  assert.ok(best < 10 * spacing, `expected within a few samples, got ${best}`);
});

test("an iteration started on the solution stays there", () => {
  const sets = [line(RED), line(BLUE)];
  for (const parameters of Object.values(namedParameters())) {
    const {start, step} = iteration(sets, parameters);
    let state = start([0, 0]);
    for (let i = 0; i < 50; i++) state = step(state);
    assert.ok(norm(state.point) < 1e-12, `moved off the intersection: ${state.point}`);
  }
});
