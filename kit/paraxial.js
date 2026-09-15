// Paraxial electron optics from the field up: analytic axial fields, the
// linearised ray equation, and the imaging properties that fall out of it.
//
// `raytrace.js` takes focal lengths as given and pushes 2x2 matrices around.
// Here a real field profile is integrated to produce that focal length, so the
// two can be checked against each other in the weak-lens limit (see
// `thinLensPower`).
//
// Ported from the Wolfram source in
// `lecture-ap3402-electron-microscopy-instrumentation/w4.02.paraxial-odes.md`,
// which is a static page of Mathematica with no interactive counterpart.
//
// UNITS ARE SI HERE — metres, volts, tesla — not the Angstroms the rest of the
// kit uses. The chapter specifies lenses in metres and fields in tesla, and
// `eta` is an SI constant, so converting would put a factor in every line for
// no gain. Nothing in this module talks to the Angstrom-based modules.

import {E_CHARGE, M_E} from "./units.js";

/** `sqrt(|e| / 2m)`, the constant the ray equation carries [C^0.5 kg^-0.5]. */
export const ETA = Math.sqrt(E_CHARGE / (2 * M_E));

// ---------------------------------------------------------------------------
// Analytic axial fields
// ---------------------------------------------------------------------------

/**
 * The nth derivative of `1/(a^2 + z^2)`, in closed form.
 *
 * Both analytic fields below are that function times a constant, and the
 * Cartesian Laplace expansion needs derivatives up to order `2N + 1` — seven of
 * them at `N = 3`. Finite differences that high lose all their significant
 * digits, so the derivative is taken in the complex plane instead.
 *
 * Writing `1/(a^2 + z^2) = (1/2ia)[1/(z - ia) - 1/(z + ia)]` and differentiating
 * each simple pole gives
 *
 *     f^(n)(z) = (-1)^(n+1) n! sin((n+1) psi) / (a rho^(n+1))
 *
 * with `rho = |z - ia|` and `psi = arg(z - ia)`. Exact at every order and the
 * same cost at every order.
 */
function poleDerivative(n, z, a) {
  const rho = Math.hypot(z, a);
  const psi = Math.atan2(-a, z);
  let factorial = 1;
  for (let i = 2; i <= n; i++) factorial *= i;
  return (
    ((n % 2 === 0 ? -1 : 1) * factorial * Math.sin((n + 1) * psi)) /
    (a * rho ** (n + 1))
  );
}

/**
 * Schiske potential — the on-axis field of an Einzel (three-electrode
 * electrostatic) lens:
 *
 *     phi(z) = phi0 (1 - k^2 / (1 + (z/a)^2))
 *
 * The electron decelerates through the middle electrode and re-accelerates,
 * focusing with no net change of energy. `k` sets the strength, `a` the axial
 * extent.
 *
 * Derivatives are analytic rather than differenced. They appear inside the ray
 * equation divided by `phi`, so a differencing error there becomes a
 * focal-length error. `test/paraxial.test.js` checks them against central
 * differences.
 *
 * @param {object} options
 * @param {number} options.a lens half-width [m]
 * @param {number} options.k strength, 0 to 1 (the potential dips to phi0(1-k^2))
 * @param {number} options.phi0 accelerating potential far from the lens [V]
 */
export function schiskeField({a, k, phi0}) {
  const k2 = k * k;
  const s = (z) => 1 + (z / a) ** 2;
  return {
    kind: "electrostatic",
    a, k, phi0,
    phi: (z) => phi0 * (1 - k2 / s(z)),
    dphi: (z) => (phi0 * k2 * ((2 * z) / (a * a))) / s(z) ** 2,
    d2phi: (z) =>
      phi0 * k2 * (2 / (a * a * s(z) ** 2) - (8 * z * z) / (a ** 4 * s(z) ** 3)),
    // phi(z) = phi0 - phi0 k^2 a^2 / (a^2 + z^2), so every derivative past the
    // zeroth is a scaled pole derivative. `test/nonparaxial.test.js` checks
    // n = 1 and n = 2 against the two written out above.
    d: (n, z) =>
      (n === 0 ? phi0 : 0) - phi0 * k2 * a * a * poleDerivative(n, z, a)
  };
}

/**
 * Glaser bell-shaped field — the on-axis field of a round magnetic lens:
 *
 *     B(z) = B0 / (1 + (z/a)^2)
 *
 * Not an exact coil geometry, but it reproduces the focusing of a real
 * pole-piece lens closely enough that analytical aberration theory is built on
 * it. `a` is the effective half-width.
 *
 * @param {object} options
 * @param {number} options.a lens half-width [m]
 * @param {number} options.B0 peak axial field [T]
 * @param {number} options.phi0 accelerating potential [V]
 */
export function glaserField({a, B0, phi0}) {
  return {
    kind: "magnetic",
    a, B0, phi0,
    B: (z) => B0 / (1 + (z / a) ** 2),
    /** nth derivative of B, for the Cartesian expansion. */
    d: (n, z) => B0 * a * a * poleDerivative(n, z, a)
  };
}

/** Superpose an electrostatic and a magnetic field on the same axis. */
export function combinedField(electrostatic, magnetic) {
  return {
    kind: "combined",
    phi0: electrostatic.phi0,
    phi: electrostatic.phi,
    dphi: electrostatic.dphi,
    d2phi: electrostatic.d2phi,
    B: magnetic.B
  };
}

// ---------------------------------------------------------------------------
// The linearised ray equation
// ---------------------------------------------------------------------------

/**
 * The two coefficients of the paraxial ray equation
 *
 *     r'' + (phi'/2phi) r' + (phi''/4phi + eta^2 B^2 / 4phi0) r = 0
 *
 * as a function of z. A purely magnetic field has no damping term, because
 * `phi` is then the constant `phi0`.
 *
 * @returns {(z: number) => {damping: number, restoring: number}}
 */
export function paraxialCoefficients(field) {
  const {phi0} = field;
  return (z) => {
    let damping = 0;
    let restoring = 0;
    if (field.phi) {
      const p = field.phi(z);
      damping = field.dphi(z) / (2 * p);
      restoring += field.d2phi(z) / (4 * p);
    }
    if (field.B) {
      const b = field.B(z);
      restoring += (ETA * ETA * b * b) / (4 * phi0);
    }
    return {damping, restoring};
  };
}

/**
 * Integrate one ray through the field with classical RK4.
 *
 * @param {object} field from {@link schiskeField} / {@link glaserField}
 * @param {object} options
 * @param {number} options.z0 start of the span [m]
 * @param {number} options.z1 end of the span [m]
 * @param {number} [options.r0=0] transverse displacement at `z0` [m]
 * @param {number} [options.rp0=0] slope at `z0`
 * @param {number} [options.steps=2000]
 * @returns {{z: Float64Array, r: Float64Array, rp: Float64Array}} the whole
 *   trajectory, slopes included — they are what lets a zero crossing be
 *   located to better than a step (see {@link firstZero}).
 */
export function integrateRay(field, {z0, z1, r0 = 0, rp0 = 0, steps = 2000}) {
  const coefficients = paraxialCoefficients(field);
  const h = (z1 - z0) / steps;

  // y = [r, r']; y' = [r', -damping r' - restoring r]
  const slope = (z, r, rp) => {
    const {damping, restoring} = coefficients(z);
    return [rp, -damping * rp - restoring * r];
  };

  const zs = new Float64Array(steps + 1);
  const rs = new Float64Array(steps + 1);
  const rps = new Float64Array(steps + 1);
  let z = z0, r = r0, rp = rp0;
  zs[0] = z; rs[0] = r; rps[0] = rp;

  for (let i = 1; i <= steps; i++) {
    const [k1r, k1v] = slope(z, r, rp);
    const [k2r, k2v] = slope(z + h / 2, r + (h / 2) * k1r, rp + (h / 2) * k1v);
    const [k3r, k3v] = slope(z + h / 2, r + (h / 2) * k2r, rp + (h / 2) * k2v);
    const [k4r, k4v] = slope(z + h, r + h * k3r, rp + h * k3v);
    r += (h / 6) * (k1r + 2 * k2r + 2 * k3r + k4r);
    rp += (h / 6) * (k1v + 2 * k2v + 2 * k3v + k4v);
    z = z0 + i * h;
    zs[i] = z; rs[i] = r; rps[i] = rp;
  }
  return {z: zs, r: rs, rp: rps};
}

/**
 * The two principal rays, which span every paraxial trajectory:
 *
 *   g — the **primary** ray, `g(z0) = 1, g'(z0) = 0`: enters off-axis and
 *       parallel, so it carries transverse position through the lens.
 *   h — the **marginal** ray, `h(z0) = 0, h'(z0) = 1`: leaves the axis at an
 *       angle, so it carries angle.
 *
 * Any ray is then `r(z) = r0 g(z) + r0' h(z)`.
 */
export function principalRays(field, span) {
  return {
    g: integrateRay(field, {...span, r0: 1, rp0: 0}),
    h: integrateRay(field, {...span, r0: 0, rp0: 1})
  };
}

/**
 * The first zero crossing after `after`, to better than one step.
 *
 * The trajectory carries its own slopes, so the segment containing the crossing
 * can be interpolated as a cubic Hermite rather than a straight line, which
 * keeps the fourth-order accuracy the integrator has. The crossing position
 * becomes the focal length, so linear interpolation here would lose most of
 * that accuracy.
 *
 * @returns {number|null} z of the crossing, or null if the ray never crosses
 */
export function firstZero({z, r, rp}, {after = -Infinity} = {}) {
  for (let i = 0; i < r.length - 1; i++) {
    if (z[i + 1] <= after) continue;
    if (r[i] === 0 && z[i] > after) return z[i];
    if (r[i] * r[i + 1] >= 0) continue;

    const h = z[i + 1] - z[i];
    const hermite = (t) => {
      const t2 = t * t, t3 = t2 * t;
      return (
        (2 * t3 - 3 * t2 + 1) * r[i] +
        (t3 - 2 * t2 + t) * h * rp[i] +
        (-2 * t3 + 3 * t2) * r[i + 1] +
        (t3 - t2) * h * rp[i + 1]
      );
    };

    let lo = 0, hi = 1;
    const sign = Math.sign(hermite(0));
    for (let n = 0; n < 60; n++) {
      const mid = (lo + hi) / 2;
      if (Math.sign(hermite(mid)) === sign) lo = mid;
      else hi = mid;
    }
    return z[i] + ((lo + hi) / 2) * h;
  }
  return null;
}

/**
 * Where the lens focuses and images, and by how much it magnifies.
 *
 *   focalPlane — where a ray entering parallel to the axis crosses it, `g = 0`
 *   imagePlane — where rays from an on-axis point meet again, `h = 0`
 *   magnification — `g` at the image plane, relative to `g` at the object
 *
 * This is the JS counterpart of the chapter's
 * `WhenEvent[r[z] == 0 && z > 0, Sow[z]]`.
 *
 * `objectPlane` defaults to the start of the span, which is the textbook case.
 * Passing it separately lets the parallel ray be integrated across a wider
 * window than the object sits in — useful when the two are drawn together and
 * the frame should not move every time the object does. The magnification is
 * then normalised by `g` at the object rather than assuming `g = 1` there, so
 * it stays right even with the object inside the field's tail.
 *
 * @param {object} field
 * @param {object} span `{z0, z1, steps}` — the window `g` is integrated over
 * @param {object} [options]
 * @param {number} [options.after=0] only count crossings past here
 * @param {number} [options.objectPlane] where `h` starts; defaults to `span.z0`
 */
export function imagingProperties(field, span, {after = 0, objectPlane} = {}) {
  const origin = objectPlane ?? span.z0;
  const g = integrateRay(field, {...span, r0: 1, rp0: 0});
  const h = integrateRay(field, {...span, z0: origin, r0: 0, rp0: 1});

  const focalPlane = firstZero(g, {after});
  const imagePlane = firstZero(h, {after});

  let magnification = null;
  if (imagePlane !== null) {
    const atObject = sampleAt(g, origin);
    if (atObject !== 0) magnification = sampleAt(g, imagePlane) / atObject;
  }

  return {g, h, focalPlane, imagePlane, magnification};
}

/** Cubic-Hermite sample of a trajectory at an arbitrary z. */
export function sampleAt({z, r, rp}, at) {
  const n = z.length;
  if (at <= z[0]) return r[0];
  if (at >= z[n - 1]) return r[n - 1];
  const h = z[1] - z[0];
  const i = Math.min(n - 2, Math.max(0, Math.floor((at - z[0]) / h)));
  const t = (at - z[i]) / h;
  const t2 = t * t, t3 = t2 * t;
  return (
    (2 * t3 - 3 * t2 + 1) * r[i] +
    (t3 - 2 * t2 + t) * h * rp[i] +
    (-2 * t3 + 3 * t2) * r[i + 1] +
    (t3 - t2) * h * rp[i + 1]
  );
}

/**
 * The thin-lens focusing power, `1/f = integral of the restoring coefficient`.
 *
 * Exact only in the weak-lens limit, where the ray barely moves while crossing
 * the field. This is the focal length `raytrace.js` takes as given, computed
 * here from the field that produces it. `test/paraxial.test.js` asserts the
 * integrated ODE converges to it as the lens weakens.
 *
 * For the Glaser field it has a closed form, `eta^2 B0^2 pi a / (8 phi0)`.
 *
 * @returns {number} 1/f [1/m]
 */
export function thinLensPower(field, {z0, z1, steps = 20000}) {
  const coefficients = paraxialCoefficients(field);
  const h = (z1 - z0) / steps;
  let sum = 0;
  for (let i = 0; i <= steps; i++) {
    const weight = i === 0 || i === steps ? 0.5 : 1;
    sum += weight * coefficients(z0 + i * h).restoring;
  }
  return sum * h;
}
