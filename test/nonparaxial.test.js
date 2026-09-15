import test from "node:test";
import assert from "node:assert/strict";

import {
  ETA, schiskeField, glaserField, integrateRay, principalRays, sampleAt
} from "../kit/paraxial.js";
import {
  laplaceExpansion, equationsOfMotion, traceRay3D, coneOfRays,
  larmorAngle, axialCrossing, spotDiagram
} from "../kit/nonparaxial.js";

// The lecture's own parameters, as in `test/paraxial.test.js`.
const SCHISKE = {a: 0.025, k: 1 / Math.SQRT2, phi0: 5};
const SCHISKE_SPAN = {z0: -0.5, z1: 1.0, steps: 4000};
const GLASER = {a: 0.0075, B0: 0.01, phi0: 1000};
const GLASER_SPAN = {z0: -0.1, z1: 0.1, steps: 4000};

// ---------------------------------------------------------------------------
// The analytic derivatives the expansion is built on
// ---------------------------------------------------------------------------

test("the nth-derivative closed form reproduces the written-out first two", () => {
  // `schiskeField` carries dphi and d2phi written out by hand from the chapter.
  // `d(n, z)` comes from the pole expansion instead, so agreeing at n = 1 and
  // n = 2 is a check on the closed form at the only orders where an independent
  // expression exists.
  const f = schiskeField(SCHISKE);
  for (let z = -0.15; z <= 0.15; z += 0.01) {
    assert.ok(Math.abs(f.d(0, z) - f.phi(z)) < 1e-12 * SCHISKE.phi0);
    assert.ok(Math.abs(f.d(1, z) - f.dphi(z)) < 1e-9 * Math.abs(f.dphi(0.01)));
    assert.ok(Math.abs(f.d(2, z) - f.d2phi(z)) < 1e-9 * Math.abs(f.d2phi(0)));
  }
});

test("high derivatives match central differences of the one below", (t) => {
  // Each derivative is checked against a difference of its predecessor rather
  // than of the profile, so the cancellation never gets worse than one order.
  const f = glaserField(GLASER);
  const h = 1e-5;
  let worst = 0;
  for (let n = 1; n <= 7; n++) {
    let scale = 0;
    for (let z = -0.02; z <= 0.02; z += 0.002) scale = Math.max(scale, Math.abs(f.d(n, z)));
    for (let z = -0.02; z <= 0.02; z += 0.002) {
      const differenced = (f.d(n - 1, z + h) - f.d(n - 1, z - h)) / (2 * h);
      worst = Math.max(worst, Math.abs(differenced - f.d(n, z)) / scale);
    }
  }
  t.diagnostic(`worst relative error over n = 1..7: ${worst.toExponential(1)}`);
  // The bound is set by the central difference, not by the derivative. Its
  // truncation error is h^2/6 times the next derivative up, and at n = 7 that
  // next derivative is about 1e6 larger, which lands right at 1e-5.
  assert.ok(worst < 1e-4, `off by ${worst}`);
});

test("the expansion's divergence is exactly its first unpaired term", (t) => {
  // Every term of div F cancels between the transverse and longitudinal sums
  // except the last one the truncation orphans:
  //
  //     div F = (-1)^N / (2^(2N) (N!)^2) r^(2N) F^(2N+1)(z)
  //
  // Asserting that exact expression checks each coefficient individually, which
  // asserting "div F is small" would not.
  const f = glaserField(GLASER);
  const factorial = (n) => {
    let out = 1;
    for (let i = 2; i <= n; i++) out *= i;
    return out;
  };
  const h = 5e-7;

  for (const order of [1, 2, 3, 4]) {
    const expand = laplaceExpansion((n, z) => f.d(n, z), order);
    // Away from z = 0, where every odd derivative vanishes by symmetry and both
    // sides are roundoff.
    for (const [x, y, z] of [[2e-3, -1e-3, 0.004], [-1e-3, 2e-3, -0.01]]) {
      const divergence =
        ((expand(x + h, y, z)[0] - expand(x - h, y, z)[0]) +
         (expand(x, y + h, z)[1] - expand(x, y - h, z)[1]) +
         (expand(x, y, z + h)[2] - expand(x, y, z - h)[2])) / (2 * h);
      const r2 = x * x + y * y;
      const predicted =
        ((-1) ** order / (2 ** (2 * order) * factorial(order) ** 2)) *
        r2 ** order * f.d(2 * order + 1, z);
      const ratio = divergence / predicted;
      if (order === 4) t.diagnostic(`order 4: div ${divergence.toExponential(3)}, predicted ${predicted.toExponential(3)}`);
      // A percent, because by order 4 the orphaned term has shrunk to within a
      // few orders of magnitude of what a central difference can resolve here.
      assert.ok(Math.abs(ratio - 1) < 1e-2,
        `order ${order} at (${x}, ${y}, ${z}): ratio ${ratio}`);
    }
  }

  // And the residual shrinks steeply with order, which is what makes a
  // truncated expansion usable at all.
  const residual = (order) => {
    const expand = laplaceExpansion((n, z) => f.d(n, z), order);
    const [x, y, z] = [2e-3, -1e-3, 0.004];
    return Math.abs(
      ((expand(x + h, y, z)[0] - expand(x - h, y, z)[0]) +
       (expand(x, y + h, z)[1] - expand(x, y - h, z)[1]) +
       (expand(x, y, z + h)[2] - expand(x, y, z - h)[2])) / (2 * h) / f.d(1, z)
    );
  };
  t.diagnostic(`residual vs order: ${[1, 2, 3, 4].map((n) => residual(n).toExponential(1)).join(" -> ")}`);
  assert.ok(residual(4) < residual(1) / 100, "the residual should fall steeply with order");
});

// ---------------------------------------------------------------------------
// Order 1 must reduce to the paraxial result
// ---------------------------------------------------------------------------

test("an electrostatic lens at order 1 is exactly the paraxial equation", (t) => {
  // For a purely electrostatic field the order-1 expansion has no cubic term at
  // all, so this holds at any radius rather than only in the limit. Getting a
  // coefficient or a sign wrong in `equationsOfMotion` breaks it immediately.
  const field = schiskeField(SCHISKE);
  const system = equationsOfMotion({electrostatic: field, phi0: SCHISKE.phi0, order: 1});

  const paraxial = integrateRay(field, {...SCHISKE_SPAN, r0: 1e-3, rp0: 0});
  const traced = traceRay3D(system, {...SCHISKE_SPAN, x0: 1e-3, xp0: 0});

  let worst = 0;
  let scale = 0;
  for (let i = 0; i < traced.z.length; i++) {
    worst = Math.max(worst, Math.abs(traced.x[i] - paraxial.r[i]));
    scale = Math.max(scale, Math.abs(paraxial.r[i]));
  }
  t.diagnostic(`max |3D - paraxial| = ${(worst / scale).toExponential(2)} of full scale`);
  assert.ok(worst / scale < 1e-12, `off by ${worst / scale}`);

  // y is untouched by an x-only ray, which is what says the two axes are not
  // accidentally coupled.
  for (let i = 0; i < traced.z.length; i++) assert.equal(traced.y[i], 0);
});

test("a magnetic lens at order 1 is the paraxial solution, rotated by the Larmor angle", (t) => {
  // Substituting u = x + iy = w exp(i theta) into the coupled equations removes
  // the coupling and leaves the paraxial radial equation for w. So a traced ray
  // must be a paraxial solution turned through theta(z).
  //
  // The initial conditions need care. A ray entering with x = r0 and no slope
  // has u'(z0) = 0, which means w'(z0) = -i Omega(z0) r0 rather than zero: the
  // rotating frame is already turning when the ray arrives. So w starts with a
  // real part r0 and no slope, and an imaginary part 0 with slope
  // -Omega(z0) r0. The paraxial equation is real, so those two evolve as the
  // primary and marginal rays independently.
  const field = glaserField(GLASER);
  const system = equationsOfMotion({magnetic: field, phi0: GLASER.phi0, order: 1});
  const {theta} = larmorAngle(field, GLASER_SPAN);
  const {g, h} = principalRays(field, GLASER_SPAN);
  const omega0 = (ETA * field.B(GLASER_SPAN.z0)) / (2 * Math.sqrt(GLASER.phi0));

  const errorAt = (r0) => {
    const traced = traceRay3D(system, {...GLASER_SPAN, x0: r0, xp0: 0});
    let worst = 0;
    let scale = 0;
    for (let i = 0; i < traced.z.length; i++) {
      const wRe = r0 * g.r[i];
      const wIm = -omega0 * r0 * h.r[i];
      const expectedX = wRe * Math.cos(theta[i]) - wIm * Math.sin(theta[i]);
      const expectedY = wRe * Math.sin(theta[i]) + wIm * Math.cos(theta[i]);
      worst = Math.max(worst, Math.hypot(traced.x[i] - expectedX, traced.y[i] - expectedY));
      scale = Math.max(scale, Math.hypot(wRe, wIm));
    }
    return worst / scale;
  };

  // Exact only in the limit, because order 1 keeps the r^2 term in B_z. So the
  // error must fall as the square of the radius.
  const coarse = errorAt(2e-4);
  const fine = errorAt(1e-4);
  t.diagnostic(`relative error ${coarse.toExponential(2)} at 0.2 mm, ${fine.toExponential(2)} at 0.1 mm, ratio ${(coarse / fine).toFixed(2)}`);
  assert.ok(coarse / fine > 3.5 && coarse / fine < 4.5,
    `error fell ${(coarse / fine).toFixed(2)}x, expected 4x for a quadratic residual`);
  assert.ok(fine < 1e-4, `still off by ${fine} at 0.1 mm`);
});

test("the paraxial option removes the residual entirely", (t) => {
  // With the longitudinal field linearised, the coupled system IS the paraxial
  // equation in a rotating frame, so the same comparison should now hold to
  // integrator precision at any radius.
  const field = glaserField(GLASER);
  const system = equationsOfMotion({
    magnetic: field, phi0: GLASER.phi0, order: 1, paraxial: true
  });
  const {theta} = larmorAngle(field, GLASER_SPAN);
  const {g, h} = principalRays(field, GLASER_SPAN);
  const omega0 = (ETA * field.B(GLASER_SPAN.z0)) / (2 * Math.sqrt(GLASER.phi0));

  const r0 = 2e-3; // ten times the radius the test above needed
  const traced = traceRay3D(system, {...GLASER_SPAN, x0: r0, xp0: 0});
  let worst = 0;
  let scale = 0;
  for (let i = 0; i < traced.z.length; i++) {
    const wRe = r0 * g.r[i];
    const wIm = -omega0 * r0 * h.r[i];
    worst = Math.max(worst, Math.hypot(
      traced.x[i] - (wRe * Math.cos(theta[i]) - wIm * Math.sin(theta[i])),
      traced.y[i] - (wRe * Math.sin(theta[i]) + wIm * Math.cos(theta[i]))
    ));
    scale = Math.max(scale, Math.hypot(wRe, wIm));
  }
  t.diagnostic(`relative error at 2 mm: ${(worst / scale).toExponential(2)}`);
  assert.ok(worst / scale < 1e-9, `off by ${worst / scale}`);
});

test("the Larmor angle matches its closed form for the Glaser field", (t) => {
  // theta = integral of eta B / 2 sqrt(phi0); for B0/(1 + (z/a)^2) that is
  // eta B0 a / (2 sqrt(phi0)) * arctan(z/a).
  const field = glaserField(GLASER);
  const {z, theta} = larmorAngle(field, GLASER_SPAN);
  const {a, B0, phi0} = GLASER;
  const k = (ETA * B0 * a) / (2 * Math.sqrt(phi0));
  const exact = (zz) => k * (Math.atan(zz / a) - Math.atan(GLASER_SPAN.z0 / a));

  let worst = 0;
  for (let i = 0; i < z.length; i++) worst = Math.max(worst, Math.abs(theta[i] - exact(z[i])));
  const total = exact(GLASER_SPAN.z1);
  t.diagnostic(`total rotation ${((total * 180) / Math.PI).toFixed(2)} deg, max error ${worst.toExponential(2)} rad`);
  assert.ok(worst < 1e-8, `off by ${worst} rad`);
  assert.ok(total > 0.5, "the lecture's lens should turn the image by a good fraction of a radian");
});

// ---------------------------------------------------------------------------
// Order 3 — spherical aberration
// ---------------------------------------------------------------------------

test("outer rays cross the axis before inner ones, and the paraxial system has no such shift", (t) => {
  // Spherical aberration. The crossing moves towards the lens by a term in the
  // square of the entry angle, so the shift at 20 mrad should be four times the
  // shift at 10 mrad. The paraxial system is the control: its crossing must not
  // depend on the entry angle at all.
  const field = glaserField(GLASER);
  const start = {z0: -0.1, z1: 0.1, steps: 8000, x0: 0, y0: 0};

  const crossingAt = (options, angle) => {
    const system = equationsOfMotion({magnetic: field, phi0: GLASER.phi0, ...options});
    return axialCrossing(traceRay3D(system, {...start, xp0: angle}), {after: 0}).z;
  };

  const flatInner = crossingAt({order: 1, paraxial: true}, 2e-3);
  const flatOuter = crossingAt({order: 1, paraxial: true}, 2e-2);
  t.diagnostic(`paraxial: ${(flatInner * 1e3).toFixed(4)} mm vs ${(flatOuter * 1e3).toFixed(4)} mm`);
  assert.ok(Math.abs(flatInner - flatOuter) < 1e-6,
    "a linearised lens must focus every angle at the same plane");

  const inner = crossingAt({order: 3}, 2e-3);
  const mid = crossingAt({order: 3}, 1e-2);
  const outer = crossingAt({order: 3}, 2e-2);
  t.diagnostic(`order 3: ${(inner * 1e3).toFixed(4)} mm, ${(mid * 1e3).toFixed(4)} mm, ${(outer * 1e3).toFixed(4)} mm`);
  t.diagnostic(`shift at 20 mrad: ${((outer - inner) * 1e6).toFixed(1)} um`);
  assert.ok(outer < inner, "the outer ray should cross first");

  const ratio = (inner - outer) / (inner - mid);
  t.diagnostic(`shift ratio between 20 mrad and 10 mrad: ${ratio.toFixed(2)}`);
  assert.ok(ratio > 3.5 && ratio < 4.5, `expected ~4, got ${ratio}`);
});

test("most of a magnetic lens's spherical aberration is already there at order 1", (t) => {
  // Worth asserting because it is easy to assume otherwise. For a magnetic lens
  // the cubic term comes from B_z's r^2 B'' term multiplying a slope, and that
  // term is present in the expansion from order 1. Going to order 3 changes the
  // crossing by only a few percent.
  //
  // An electrostatic lens behaves the other way round: its equations never use
  // E_z, so order 1 carries no cubic term at all and order 3 is where
  // aberration first appears.
  const field = glaserField(GLASER);
  const start = {z0: -0.1, z1: 0.1, steps: 8000, xp0: 2e-2};
  const crossingAt = (options) => axialCrossing(
    traceRay3D(equationsOfMotion({magnetic: field, phi0: GLASER.phi0, ...options}), start),
    {after: 0}
  ).z;

  const flat = crossingAt({order: 1, paraxial: true});
  const first = crossingAt({order: 1});
  const third = crossingAt({order: 3});
  const fourth = crossingAt({order: 4});
  const fifth = crossingAt({order: 5});

  const firstShare = (flat - first) / (flat - third);
  t.diagnostic(`paraxial ${(flat * 1e3).toFixed(4)} mm, order 1 ${(first * 1e3).toFixed(4)} mm, order 3 ${(third * 1e3).toFixed(4)} mm`);
  t.diagnostic(`order 1 captures ${(firstShare * 100).toFixed(1)}% of the order-3 shift`);
  assert.ok(firstShare > 0.9, `order 1 captured only ${firstShare}`);

  // The expansion has converged by order 3, which is why the widget stops
  // there: each further order changes the crossing by far less than the last.
  t.diagnostic(`order 3->4: ${((fourth - third) * 1e9).toFixed(2)} nm, order 4->5: ${((fifth - fourth) * 1e9).toFixed(2)} nm`);
  assert.ok(Math.abs(fourth - third) < 1e-4 * Math.abs(third - flat),
    "order 4 should barely move the crossing");
  assert.ok(Math.abs(fifth - fourth) < Math.abs(fourth - third) / 10,
    "the expansion should keep converging");

  const schiske = schiskeField(SCHISKE);
  const schiskeStart = {...SCHISKE_SPAN, steps: 8000, xp0: 2e-2};
  const schiskeAt = (order) => axialCrossing(
    traceRay3D(equationsOfMotion({electrostatic: schiske, phi0: SCHISKE.phi0, order}), schiskeStart),
    {after: 0}
  ).z;
  const s1 = schiskeAt(1);
  const s3 = schiskeAt(3);
  t.diagnostic(`Schiske: order 1 ${(s1 * 1e3).toFixed(3)} mm, order 3 ${(s3 * 1e3).toFixed(3)} mm`);
  assert.ok(s3 < s1, "the electrostatic lens should gain its aberration at order 3");
});

test("a cone of rays lands in a disc whose size grows with the shell", (t) => {
  const field = glaserField(GLASER);
  const system = equationsOfMotion({magnetic: field, phi0: GLASER.phi0, order: 3});
  const starts = coneOfRays({azimuthal: 12, shells: 3, minAngle: 5e-3, maxAngle: 2e-2});
  assert.equal(starts.length, 36);

  const rays = starts.map((s) => ({
    ...traceRay3D(system, {z0: -0.1, z1: 0.1, steps: 4000, ...s}),
    shell: s.shell, angle: s.angle
  }));

  // The Gaussian image plane, from a ray so close to the axis that aberration
  // cannot reach it.
  const gaussian = axialCrossing(
    traceRay3D(equationsOfMotion({magnetic: field, phi0: GLASER.phi0, order: 1}),
      {z0: -0.1, z1: 0.1, steps: 4000, xp0: 1e-4}),
    {after: 0}
  ).z;

  const spots = spotDiagram(rays, gaussian);
  const radiusOf = (shell) => {
    const r = spots.filter((s) => s.shell === shell).map((s) => Math.hypot(s.x, s.y));
    return r.reduce((a, b) => a + b, 0) / r.length;
  };
  const [r0, r1, r2] = [radiusOf(0), radiusOf(1), radiusOf(2)];
  t.diagnostic(`spot radii at z = ${gaussian.toFixed(4)} m: ${(r0 * 1e6).toFixed(2)}, ${(r1 * 1e6).toFixed(2)}, ${(r2 * 1e6).toFixed(2)} um`);
  assert.ok(r0 < r1 && r1 < r2, "the disc must grow with entry angle");

  // Each shell should be a ring, so the spread within a shell is small next to
  // its radius. Anything else means the rays are not being rotated together.
  const shell2 = spots.filter((s) => s.shell === 2).map((s) => Math.hypot(s.x, s.y));
  const spread = Math.max(...shell2) - Math.min(...shell2);
  assert.ok(spread < 0.02 * r2, `outer shell is not a ring: spread ${spread} vs radius ${r2}`);
});

test("the coupled integrator is fourth order", (t) => {
  const field = glaserField(GLASER);
  const system = equationsOfMotion({magnetic: field, phi0: GLASER.phi0, order: 3});
  const at = (steps) => {
    const r = traceRay3D(system, {z0: -0.1, z1: 0.1, steps, xp0: 1e-2, yp0: 5e-3});
    return [r.x.at(-1), r.y.at(-1)];
  };
  const [ex, ey] = at(64000);
  const coarse = Math.hypot(at(250)[0] - ex, at(250)[1] - ey);
  const fine = Math.hypot(at(500)[0] - ex, at(500)[1] - ey);
  t.diagnostic(`error ${coarse.toExponential(2)} -> ${fine.toExponential(2)}, ratio ${(coarse / fine).toFixed(1)}`);
  assert.ok(coarse / fine > 10, `convergence ratio ${coarse / fine}, expected ~16`);
});

test("a field-free span leaves rays straight in both coordinates", () => {
  const system = equationsOfMotion({magnetic: {d: () => 0}, phi0: 1000, order: 3});
  const {z, x, y} = traceRay3D(system, {z0: 0, z1: 1, steps: 100, xp0: 0.01, yp0: -0.02});
  for (let i = 0; i < z.length; i++) {
    assert.ok(Math.abs(x[i] - 0.01 * z[i]) < 1e-14);
    assert.ok(Math.abs(y[i] + 0.02 * z[i]) < 1e-14);
  }
});

test("sampling a bundle at a plane agrees with the trajectory there", () => {
  // `spotDiagram` interpolates; checked against a grid point, where it must
  // return the sample itself.
  const field = schiskeField(SCHISKE);
  const system = equationsOfMotion({electrostatic: field, phi0: SCHISKE.phi0, order: 3});
  const ray = traceRay3D(system, {...SCHISKE_SPAN, xp0: 3e-3, yp0: 1e-3});
  const i = 1234;
  const [spot] = spotDiagram([ray], ray.z[i]);
  assert.ok(Math.abs(spot.x - ray.x[i]) < 1e-15);
  assert.ok(Math.abs(spot.y - ray.y[i]) < 1e-15);

  // And halfway between two samples it should sit between them.
  const mid = (ray.z[i] + ray.z[i + 1]) / 2;
  const [between] = spotDiagram([ray], mid);
  const lo = Math.min(ray.x[i], ray.x[i + 1]);
  const hi = Math.max(ray.x[i], ray.x[i + 1]);
  assert.ok(between.x >= lo && between.x <= hi);
});

test("sampleAt and the 3D tracer agree for a purely radial electrostatic ray", () => {
  // Ties the two modules together at one number, so a future change to either
  // interpolator cannot drift them apart unnoticed.
  const field = schiskeField(SCHISKE);
  const system = equationsOfMotion({electrostatic: field, phi0: SCHISKE.phi0, order: 1});
  const paraxial = integrateRay(field, {...SCHISKE_SPAN, r0: 1e-3, rp0: 0});
  const traced = traceRay3D(system, {...SCHISKE_SPAN, x0: 1e-3, xp0: 0});
  const [spot] = spotDiagram([traced], 0.42);
  assert.ok(Math.abs(spot.x - sampleAt(paraxial, 0.42)) < 1e-15);
});
