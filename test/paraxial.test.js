import test from "node:test";
import assert from "node:assert/strict";

import {
  ETA, schiskeField, glaserField, paraxialCoefficients,
  integrateRay, principalRays, firstZero, imagingProperties, thinLensPower
} from "../kit/paraxial.js";

// The lecture's own parameters, so the widget's defaults are the ones tested.
// `w4.02.paraxial-odes.md`, the two `paraxialSolution` blocks.
const SCHISKE = {a: 0.025, k: 1 / Math.SQRT2, phi0: 5};
const SCHISKE_SPAN = {z0: -0.5, z1: 1.0, steps: 4000};
const GLASER = {a: 0.0075, B0: 0.01, phi0: 1000};
const GLASER_SPAN = {z0: -0.1, z1: 0.1, steps: 4000};

test("eta is sqrt(e/2m)", () => {
  assert.ok(Math.abs(ETA / 2.9654e5 - 1) < 1e-3, `eta = ${ETA}`);
});

test("the Schiske derivatives are the analytic ones", (t) => {
  const f = schiskeField(SCHISKE);
  const zs = [];
  for (let z = -0.2; z <= 0.2; z += 0.005) zs.push(z);

  // Normalised against the largest value over the sweep, not pointwise: phi'
  // passes through zero at the lens centre, where a pointwise relative error is
  // a ratio of two roundoff-sized numbers and means nothing.
  const scale1 = Math.max(...zs.map((z) => Math.abs(f.dphi(z))));
  const scale2 = Math.max(...zs.map((z) => Math.abs(f.d2phi(z))));

  const h = 1e-6;
  let worstFirst = 0, worstSecond = 0;
  for (const z of zs) {
    const d1 = (f.phi(z + h) - f.phi(z - h)) / (2 * h);
    const d2 = (f.phi(z + h) - 2 * f.phi(z) + f.phi(z - h)) / (h * h);
    worstFirst = Math.max(worstFirst, Math.abs(d1 - f.dphi(z)) / scale1);
    worstSecond = Math.max(worstSecond, Math.abs(d2 - f.d2phi(z)) / scale2);
  }
  t.diagnostic(`worst error, as a fraction of full scale: phi' ${worstFirst.toExponential(1)}, phi'' ${worstSecond.toExponential(1)}`);
  assert.ok(worstFirst < 1e-8, `phi' off by ${worstFirst}`);
  // Second differences lose about half the mantissa to cancellation, so this
  // bound is set by the finite difference, not by the derivative.
  assert.ok(worstSecond < 1e-4, `phi'' off by ${worstSecond}`);
});

test("a constant restoring field integrates to a cosine, exactly", (t) => {
  // B constant => restoring is constant => r'' + K r = 0, whose solution with
  // r(0)=1, r'(0)=0 is cos(sqrt(K) z). An exact reference with no derivation.
  const B0 = 0.01, phi0 = 1000;
  const field = {phi0, B: () => B0};
  const K = (ETA * ETA * B0 * B0) / (4 * phi0);
  const {z, r} = integrateRay(field, {z0: 0, z1: 0.05, r0: 1, rp0: 0, steps: 2000});
  let worst = 0;
  for (let i = 0; i < z.length; i++) {
    worst = Math.max(worst, Math.abs(r[i] - Math.cos(Math.sqrt(K) * z[i])));
  }
  t.diagnostic(`max |RK4 - cos| = ${worst.toExponential(2)}`);
  assert.ok(worst < 1e-12, `off by ${worst}`);
});

test("the integrator is fourth order", (t) => {
  // Halving the step must cut the error about sixteenfold. This is what
  // separates a correct RK4 from one whose stages are subtly wrong — such a
  // thing still converges, just at a lower order.
  const field = glaserField(GLASER);
  const reference = integrateRay(field, {...GLASER_SPAN, steps: 64000, r0: 1, rp0: 0});
  const endOf = (steps) =>
    integrateRay(field, {...GLASER_SPAN, steps, r0: 1, rp0: 0}).r.at(-1);
  const exact = reference.r.at(-1);

  const coarse = Math.abs(endOf(250) - exact);
  const fine = Math.abs(endOf(500) - exact);
  const ratio = coarse / fine;
  t.diagnostic(`error ${coarse.toExponential(2)} -> ${fine.toExponential(2)}, ratio ${ratio.toFixed(1)}`);
  assert.ok(ratio > 10, `convergence ratio ${ratio}, expected ~16`);
});

/**
 * The Glaser field's exact solution.
 *
 * Substituting z = a cot(theta) and r = u / sin(theta) turns
 *     r'' + (eta^2 B^2 / 4 phi0) r = 0,   B = B0/(1 + (z/a)^2)
 * into the harmonic equation u'' + omega^2 u = 0 with
 *     omega^2 = 1 + eta^2 B0^2 a^2 / (4 phi0),
 * because the field's z-dependence is exactly absorbed by the substitution.
 * This is the classical result for the bell-shaped field, and it gives the
 * integrator something to be checked against that is not another integrator.
 */
function glaserExact({a, B0, phi0}, {z0, r0, rp0}) {
  const K = (ETA * ETA * B0 * B0) / (4 * phi0);
  const omega = Math.sqrt(1 + K * a * a);
  const theta = (z) => Math.atan2(a, z);          // cot(theta) = z/a, theta in (0, pi)

  const t0 = theta(z0);
  const u0 = r0 * Math.sin(t0);
  // r_theta = r_z dz/dtheta = -a r_z / sin^2(theta); u = r sin(theta).
  const du0 = (-a * rp0) / Math.sin(t0) + r0 * Math.cos(t0);

  return (z) => {
    const t = theta(z);
    const d = t - t0;
    const u = u0 * Math.cos(omega * d) + (du0 / omega) * Math.sin(omega * d);
    return u / Math.sin(t);
  };
}

test("the Glaser solution matches its closed form", (t) => {
  const field = glaserField(GLASER);
  for (const [r0, rp0, name] of [[1, 0, "primary"], [0, 1, "marginal"]]) {
    const {z, r} = integrateRay(field, {...GLASER_SPAN, r0, rp0});
    const exact = glaserExact(GLASER, {z0: GLASER_SPAN.z0, r0, rp0});
    let worst = 0, scale = 0;
    for (let i = 0; i < z.length; i++) {
      worst = Math.max(worst, Math.abs(r[i] - exact(z[i])));
      scale = Math.max(scale, Math.abs(exact(z[i])));
    }
    t.diagnostic(`${name} ray: max |diff| = ${(worst / scale).toExponential(2)} of full scale`);
    assert.ok(worst / scale < 1e-9, `${name} ray off by ${worst / scale}`);
  }
});

test("the focal plane matches the closed form's zero crossing", (t) => {
  const field = glaserField(GLASER);
  const {focalPlane, imagePlane, magnification} = imagingProperties(field, GLASER_SPAN);
  assert.ok(focalPlane !== null, "the primary ray should cross the axis");

  // Bisect the analytic solution for the same crossing.
  const exact = glaserExact(GLASER, {z0: GLASER_SPAN.z0, r0: 1, rp0: 0});
  let lo = 0, hi = GLASER_SPAN.z1;
  assert.ok(exact(lo) * exact(hi) < 0, "the analytic ray should bracket a zero");
  for (let i = 0; i < 80; i++) {
    const mid = (lo + hi) / 2;
    if (exact(lo) * exact(mid) <= 0) hi = mid; else lo = mid;
  }
  const analytic = (lo + hi) / 2;
  t.diagnostic(`focal plane ${focalPlane.toFixed(8)} m vs analytic ${analytic.toFixed(8)} m`);
  t.diagnostic(`image plane ${imagePlane?.toFixed(5)} m, magnification ${magnification?.toFixed(4)}`);

  // A step is 5e-5 m here; the Hermite refinement must do far better than that.
  assert.ok(Math.abs(focalPlane - analytic) < 1e-9, `off by ${Math.abs(focalPlane - analytic)}`);
});

test("a weak lens converges to the thin-lens focal length", (t) => {
  // 1/f = integral of the restoring coefficient, which for the Glaser field is
  // eta^2 B0^2 pi a / (8 phi0). Exact only when the ray barely moves crossing
  // the field, so the agreement must improve as the lens weakens — and this is
  // the number `raytrace.js` takes as given.
  const a = 0.0075, phi0 = 1000;
  const span = {z0: -0.5, z1: 8.0, steps: 40000};
  let previous = Infinity;
  for (const B0 of [4e-3, 2e-3, 1e-3]) {
    const field = glaserField({a, B0, phi0});
    const closedForm = (ETA * ETA * B0 * B0 * Math.PI * a) / (8 * phi0);
    assert.ok(Math.abs(thinLensPower(field, span) / closedForm - 1) < 1e-6,
      "the integrated power should match its closed form");

    const {focalPlane} = imagingProperties(field, span);
    const error = Math.abs(1 / focalPlane / closedForm - 1);
    t.diagnostic(`B0 = ${B0} T: f = ${focalPlane.toFixed(4)} m, thin lens ${(1 / closedForm).toFixed(4)} m, off by ${(error * 100).toFixed(2)}%`);
    assert.ok(error < previous, "weakening the lens should improve the agreement");
    previous = error;
  }
  assert.ok(previous < 0.02, `weakest lens still off by ${previous}`);
});

test("the Schiske lens images, at the lecture's own parameters", (t) => {
  const field = schiskeField(SCHISKE);
  const {focalPlane, imagePlane, magnification} = imagingProperties(field, SCHISKE_SPAN);
  t.diagnostic(`focal ${focalPlane?.toFixed(4)} m, image ${imagePlane?.toFixed(4)} m, M = ${magnification?.toFixed(4)}`);
  assert.ok(focalPlane !== null && focalPlane > 0, `focal plane ${focalPlane}`);
  assert.ok(imagePlane !== null && imagePlane > focalPlane,
    "the image plane must sit beyond the focal plane for a real object");
  assert.ok(magnification < 0, "a real image through a single lens is inverted");
});

test("the principal rays are what their initial conditions say", () => {
  const field = schiskeField(SCHISKE);
  const {g, h} = principalRays(field, SCHISKE_SPAN);
  assert.equal(g.r[0], 1);
  assert.equal(g.rp[0], 0);
  assert.equal(h.r[0], 0);
  assert.equal(h.rp[0], 1);
});

test("a field-free span leaves rays straight", () => {
  // No lens at all: the primary ray stays put and the marginal ray is a line of
  // unit slope. Catches a stray sign or a restoring term that is never zero.
  const field = {phi0: 1000, B: () => 0};
  const {g, h} = principalRays(field, {z0: 0, z1: 1, steps: 100});
  for (let i = 0; i < g.r.length; i++) {
    assert.ok(Math.abs(g.r[i] - 1) < 1e-12);
    assert.ok(Math.abs(h.r[i] - h.z[i]) < 1e-12);
  }
});

test("paraxialCoefficients drops the damping term for a magnetic lens", () => {
  // phi is constant there, so phi'/2phi vanishes — the chapter writes that case
  // as a separate equation, and conflating them would add a spurious drag.
  const magnetic = paraxialCoefficients(glaserField(GLASER));
  assert.equal(magnetic(0.003).damping, 0);
  assert.ok(magnetic(0.003).restoring > 0);

  const electrostatic = paraxialCoefficients(schiskeField(SCHISKE));
  assert.ok(Math.abs(electrostatic(0.01).damping) > 0, "an Einzel lens does damp");
});
