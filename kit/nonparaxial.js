// Non-paraxial ray tracing in Cartesian 3D.
//
// `paraxial.js` uses the cylindrical symmetry of a round lens to reduce the
// problem to one radial coordinate. That hides two things. Larmor rotation
// becomes an effective restoring force rather than an actual rotation, and
// there is no radius at which the linearisation is allowed to fail, so there
// are no aberrations.
//
// This module keeps x and y separate and expands the off-axis field from the
// axial one to a chosen order. Order 1 is the paraxial case. Order 3 is where
// spherical aberration appears.
//
// Ported from the Wolfram source in
// `lecture-ap3402-electron-microscopy-instrumentation/w4.04.non-paraxial-odes.md`.
//
// UNITS ARE SI, as in `paraxial.js` — metres, volts, tesla.

import {ETA} from "./paraxial.js";

// ---------------------------------------------------------------------------
// Cartesian Laplace expansion
// ---------------------------------------------------------------------------

/**
 * Expand an axisymmetric field off the axis from its axial profile `F(z)`,
 * satisfying Laplace's equation term by term:
 *
 *     F_z = sum_{n=0..N}  (-1)^n / (2^(2n) (n!)^2)      r^(2n)   F^(2n)(z)
 *     F_x = sum_{n=1..N}  (-1)^n / (2^(2n-1) (n-1)! n!) x r^(2n-2) F^(2n-1)(z)
 *     F_y = same as F_x with y in place of x
 *
 * The axial profile enters only through its derivatives, so the caller passes a
 * function of `(n, z)` rather than the profile itself.
 *
 * Truncating at N leaves exactly one unpaired term in the divergence,
 *
 *     div F = (-1)^N / (2^(2N) (N!)^2) r^(2N) F^(2N+1)(z),
 *
 * since every lower-order term cancels between the transverse and longitudinal
 * sums. `test/nonparaxial.test.js` asserts that identity, which checks the
 * coefficients rather than just their consequences.
 *
 * The transverse and longitudinal sums can be truncated at different orders.
 * Even at `N = 1` the longitudinal component carries an `r^2` term, which the
 * Wolfram source keeps and which becomes a cubic term once it multiplies a
 * slope in the equations of motion. Passing a lower `longitudinalOrder` drops
 * it, giving the strictly linearised field.
 *
 * @param {(n: number, z: number) => number} derivative nth derivative of F
 * @param {number} order N for the transverse components
 * @param {number} [longitudinalOrder=order] N for the longitudinal component
 * @returns {(x: number, y: number, z: number) => [number, number, number]}
 */
export function laplaceExpansion(derivative, order, longitudinalOrder = order) {
  const factorial = (n) => {
    let out = 1;
    for (let i = 2; i <= n; i++) out *= i;
    return out;
  };
  const transverse = [];
  const longitudinal = [];
  for (let n = 1; n <= order; n++) {
    transverse[n] =
      (-1) ** n / (2 ** (2 * n - 1) * factorial(n - 1) * factorial(n));
  }
  for (let n = 0; n <= longitudinalOrder; n++) {
    longitudinal[n] = (-1) ** n / (2 ** (2 * n) * factorial(n) ** 2);
  }

  return (x, y, z) => {
    const r2 = x * x + y * y;
    let fx = 0;
    let fy = 0;
    let fz = 0;
    for (let n = 1; n <= order; n++) {
      const term = transverse[n] * r2 ** (n - 1) * derivative(2 * n - 1, z);
      fx += x * term;
      fy += y * term;
    }
    for (let n = 0; n <= longitudinalOrder; n++) {
      fz += longitudinal[n] * r2 ** n * derivative(2 * n, z);
    }
    return [fx, fy, fz];
  };
}

/**
 * Wrap a field so repeated derivative evaluations at the same `z` are computed
 * once.
 *
 * Tracing a bundle is the same RK4 grid over and over: every ray is evaluated
 * at the same `z0 + i h` and `z0 + (i + 1/2) h`, and at order 3 each evaluation
 * needs eight derivatives. Caching on `z` turns the cost from "per ray" into
 * "per grid point", which for a bundle of fifty rays is most of the work.
 *
 * The key is the float `z` itself, which is safe here because the tracer
 * computes those values the same way each time rather than accumulating them.
 *
 * @param {object} field a field with `d(n, z)`
 * @returns {object} the same field with a memoised `d`
 */
export function memoiseDerivatives(field) {
  const cache = new Map();
  return {
    ...field,
    d(n, z) {
      let row = cache.get(z);
      if (row === undefined) {
        row = [];
        cache.set(z, row);
      }
      let value = row[n];
      if (value === undefined) {
        value = field.d(n, z);
        row[n] = value;
      }
      return value;
    }
  };
}

// ---------------------------------------------------------------------------
// Equations of motion
// ---------------------------------------------------------------------------

/**
 * The coupled ray equations in the paraxial longitudinal gauge:
 *
 *     x'' + (phi'/2phi) x' = -E_x / 2phi + (eta/sqrt(phi0)) (B_y - B_z y')
 *     y'' + (phi'/2phi) y' = -E_y / 2phi - (eta/sqrt(phi0)) (B_x - B_z x')
 *
 * with `E = -grad phi` expanded from the axial `phi'`, and `B` expanded from
 * the axial `B`. The electrostatic damping term is present only when there is
 * an electrostatic field, since a purely magnetic lens has `phi = phi0`
 * everywhere.
 *
 * `paraxial: true` linearises the system, by truncating the longitudinal field
 * one order below the transverse one. For a magnetic lens this is the
 * difference between the chapter's coupled paraxial equations and its full
 * ones: `B_z` then contributes only `B(z) y'`, which is linear, instead of also
 * `-r^2 B''(z) y' / 4`, which is cubic and is most of the spherical aberration.
 * It has no effect on an electrostatic lens, whose equations never use `E_z`.
 *
 * @param {object} options
 * @param {object} [options.electrostatic] a field with `phi(z)` and `d(n, z)`
 * @param {object} [options.magnetic] a field with `d(n, z)`
 * @param {number} options.phi0 accelerating potential [V]
 * @param {number} [options.order=3] expansion order
 * @param {boolean} [options.paraxial=false] linearise, as described above
 * @returns {(z, x, xp, y, yp) => [number, number]} the two second derivatives
 */
export function equationsOfMotion({
  electrostatic = null, magnetic = null, phi0, order = 3, paraxial = false
}) {
  const longitudinal = paraxial ? order - 1 : order;
  // The electric field is expanded from phi', so the nth derivative of the
  // profile being expanded is the (n+1)th derivative of the potential.
  const expandE = electrostatic
    ? laplaceExpansion((n, z) => electrostatic.d(n + 1, z), order, longitudinal)
    : null;
  const expandB = magnetic
    ? laplaceExpansion((n, z) => magnetic.d(n, z), order, longitudinal)
    : null;
  const scale = ETA / Math.sqrt(phi0);

  return (z, x, xp, y, yp) => {
    let xpp = 0;
    let ypp = 0;

    if (expandE) {
      const [gx, gy] = expandE(x, y, z);
      const phi = electrostatic.phi(z);
      const damping = electrostatic.d(1, z) / (2 * phi);
      // E = -grad phi, so the expansion's transverse components flip sign; the
      // equation then divides by 2phi. The two sign flips leave `+gx / 2phi`.
      xpp += gx / (2 * phi) - damping * xp;
      ypp += gy / (2 * phi) - damping * yp;
    }

    if (expandB) {
      const [bx, by, bz] = expandB(x, y, z);
      xpp += scale * (by - bz * yp);
      ypp -= scale * (bx - bz * xp);
    }

    return [xpp, ypp];
  };
}

/**
 * Integrate one ray through the coupled system with classical RK4 on the state
 * `[x, x', y, y']`.
 *
 * @param {(z, x, xp, y, yp) => [number, number]} system from {@link equationsOfMotion}
 * @param {object} options
 * @param {number} options.z0 start of the span [m]
 * @param {number} options.z1 end of the span [m]
 * @param {number} [options.x0=0] entry position [m]
 * @param {number} [options.y0=0]
 * @param {number} [options.xp0=0] entry slope
 * @param {number} [options.yp0=0]
 * @param {number} [options.steps=2000]
 * @returns {{z, x, y, xp, yp}} five Float64Arrays of length `steps + 1`
 */
export function traceRay3D(system, {
  z0, z1, x0 = 0, y0 = 0, xp0 = 0, yp0 = 0, steps = 2000
}) {
  const h = (z1 - z0) / steps;
  const n = steps + 1;
  const zs = new Float64Array(n);
  const xs = new Float64Array(n);
  const ys = new Float64Array(n);
  const xps = new Float64Array(n);
  const yps = new Float64Array(n);

  // state = [x, x', y, y']; derivative = [x', x'', y', y'']
  const slope = (z, s) => {
    const [xpp, ypp] = system(z, s[0], s[1], s[2], s[3]);
    return [s[1], xpp, s[3], ypp];
  };
  const step = (s, k, f) => [s[0] + f * k[0], s[1] + f * k[1], s[2] + f * k[2], s[3] + f * k[3]];

  let z = z0;
  let s = [x0, xp0, y0, yp0];
  zs[0] = z; xs[0] = s[0]; xps[0] = s[1]; ys[0] = s[2]; yps[0] = s[3];

  for (let i = 1; i <= steps; i++) {
    const k1 = slope(z, s);
    const k2 = slope(z + h / 2, step(s, k1, h / 2));
    const k3 = slope(z + h / 2, step(s, k2, h / 2));
    const k4 = slope(z + h, step(s, k3, h));
    for (let j = 0; j < 4; j++) {
      s[j] += (h / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
    }
    z = z0 + i * h;
    zs[i] = z; xs[i] = s[0]; xps[i] = s[1]; ys[i] = s[2]; yps[i] = s[3];
  }

  return {z: zs, x: xs, y: ys, xp: xps, yp: yps};
}

/**
 * Starting conditions for a hollow cone of rays, as concentric shells of
 * increasing entry slope.
 *
 * Shells matter more than a filled cone here. Spherical aberration scales as
 * the cube of the entry angle, so a ray's crossing point depends on which shell
 * it is on, and colouring by shell is what makes that visible.
 *
 * @param {object} options
 * @param {number} [options.azimuthal=16] rays per shell
 * @param {number} [options.shells=3]
 * @param {number} [options.minAngle=1e-3] slope of the innermost shell
 * @param {number} [options.maxAngle=1e-2] slope of the outermost shell
 * @param {number} [options.x0=0] common entry position [m]
 * @param {number} [options.y0=0]
 * @returns {Array<{x0, y0, xp0, yp0, shell: number, angle: number}>}
 */
export function coneOfRays({
  azimuthal = 16, shells = 3, minAngle = 1e-3, maxAngle = 1e-2, x0 = 0, y0 = 0
} = {}) {
  const out = [];
  for (let s = 0; s < shells; s++) {
    const angle = shells === 1
      ? maxAngle
      : minAngle + ((maxAngle - minAngle) * s) / (shells - 1);
    for (let i = 0; i < azimuthal; i++) {
      const phi = (2 * Math.PI * i) / azimuthal;
      out.push({
        x0, y0,
        xp0: angle * Math.cos(phi),
        yp0: angle * Math.sin(phi),
        shell: s, angle
      });
    }
  }
  return out;
}

// ---------------------------------------------------------------------------
// Larmor rotation
// ---------------------------------------------------------------------------

/**
 * The Larmor rotation angle accumulated along the axis,
 *
 *     theta(z) = integral of eta B(z) / (2 sqrt(phi0)) dz
 *
 * A magnetic lens focuses in a frame rotating by this angle. Substituting
 * `u = x + iy = w exp(i theta)` into the coupled equations removes the coupling
 * and leaves the paraxial radial equation for `w`, which is what
 * `test/nonparaxial.test.js` checks the traced rays against.
 *
 * Simpson's rule on each interval, using its midpoint, so the result is
 * cumulative and fourth-order. Trapezoid is second-order, which over a field as
 * peaked as the Glaser profile leaves an error of about a microradian — small,
 * but the rotation is drawn, so it is worth not having.
 *
 * @returns {{z: Float64Array, theta: Float64Array}}
 */
export function larmorAngle(field, {z0, z1, steps = 2000}) {
  const h = (z1 - z0) / steps;
  const zs = new Float64Array(steps + 1);
  const theta = new Float64Array(steps + 1);
  const omega = (z) => (ETA * field.B(z)) / (2 * Math.sqrt(field.phi0));
  zs[0] = z0;
  for (let i = 1; i <= steps; i++) {
    const za = z0 + (i - 1) * h;
    const zb = z0 + i * h;
    zs[i] = zb;
    theta[i] = theta[i - 1] + (h / 6) * (omega(za) + 4 * omega((za + zb) / 2) + omega(zb));
  }
  return {z: zs, theta};
}

// ---------------------------------------------------------------------------
// Reading aberration off a bundle
// ---------------------------------------------------------------------------

/**
 * Where a traced ray crosses the axis, to better than one step.
 *
 * "Crosses the axis" means the radius `sqrt(x^2 + y^2)` reaches a local
 * minimum, because an aberrated skew ray misses the axis rather than touching
 * it. The minimum is located by fitting a parabola to the three samples around
 * the smallest one.
 *
 * @param {object} ray from {@link traceRay3D}
 * @param {object} [options]
 * @param {number} [options.after=-Infinity] ignore minima before here
 * @returns {{z: number, radius: number}|null}
 */
export function axialCrossing({z, x, y}, {after = -Infinity} = {}) {
  let best = -1;
  let bestR2 = Infinity;
  for (let i = 1; i < z.length - 1; i++) {
    if (z[i] <= after) continue;
    const r2 = x[i] * x[i] + y[i] * y[i];
    if (r2 < bestR2) {
      bestR2 = r2;
      best = i;
    }
  }
  if (best <= 0 || best >= z.length - 1) return null;

  const at = (i) => Math.hypot(x[i], y[i]);
  const [a, b, c] = [at(best - 1), at(best), at(best + 1)];
  const denom = a - 2 * b + c;
  // A flat or non-convex triple means no well-defined minimum; fall back to the
  // sample itself rather than dividing by something near zero.
  const shift = Math.abs(denom) < 1e-30 ? 0 : (0.5 * (a - c)) / denom;
  const h = z[1] - z[0];
  return {
    z: z[best] + Math.max(-1, Math.min(1, shift)) * h,
    radius: b - 0.25 * (a - c) * shift
  };
}

/**
 * Transverse positions of a bundle of rays at one axial plane, for a spot
 * diagram.
 *
 * @param {Array<object>} rays traced rays
 * @param {number} at axial position [m]
 * @returns {Array<{x: number, y: number, shell: number, angle: number}>}
 */
export function spotDiagram(rays, at) {
  return rays.map((ray) => {
    const {z, x, y, xp, yp} = ray;
    const n = z.length;
    const h = z[1] - z[0];
    const i = Math.min(n - 2, Math.max(0, Math.floor((at - z[0]) / h)));
    const t = Math.max(0, Math.min(1, (at - z[i]) / h));
    const hermite = (v, vp) => {
      const t2 = t * t;
      const t3 = t2 * t;
      return (
        (2 * t3 - 3 * t2 + 1) * v[i] +
        (t3 - 2 * t2 + t) * h * vp[i] +
        (-2 * t3 + 3 * t2) * v[i + 1] +
        (t3 - t2) * h * vp[i + 1]
      );
    };
    return {
      x: hermite(x, xp),
      y: hermite(y, yp),
      shell: ray.shell,
      angle: ray.angle
    };
  });
}
