// Parallel-beam tomography: the projection of a volume, and its adjoint.
//
// Two functions, and the second one is the transpose of the first. That matters
// more than it sounds: every iterative reconstruction in this file's orbit —
// SIRT, SART, and the joint ptychographic-tomographic solvers — is built from a
// forward operator and its adjoint, and if the pair is not an exact transpose
// the iteration converges to the wrong thing while looking healthy. The test
// suite checks `<A x, y> = <x, Aᵀ y>` to machine precision, which is the only
// way to know.
//
// Everything is two-dimensional here: the volume is a slice and a projection is
// a line, which is the whole of tomography with one dimension taken out and the
// arithmetic small enough to watch.
//
// Conventions follow the rest of the kit: arrays are flat and row-major,
// `[ix * n + iy]`, with `ix` the vertical axis. At an angle of zero the rays run
// along `iy`, so the projection is a sum across each row and is indexed by `ix`.

/**
 * Line integrals through a square slice, at one angle.
 *
 * Rays run parallel, one per sample of `out`, and each is integrated by stepping
 * along it and sampling the volume bilinearly. The detector shares the volume's
 * sampling and its centre.
 *
 * @param {ArrayLike<number>} volume n*n
 * @param {number} n
 * @param {number} theta radians
 * @param {Float64Array} [out] n
 * @param {number} [steps] samples along each ray; defaults to n
 * @returns {Float64Array} n
 */
export function projectVolume(volume, n, theta, out = new Float64Array(n), steps = n) {
  out.fill(0);
  const centre = (n - 1) / 2;
  const cos = Math.cos(theta);
  const sin = Math.sin(theta);
  // The ray direction, and the detector direction across it.
  const step = n / steps;
  for (let d = 0; d < n; d++) {
    const s = d - centre;
    let sum = 0;
    for (let k = 0; k < steps; k++) {
      const t = (k + 0.5) * step - n / 2;
      // A point at distance `t` along the ray, offset `s` across it.
      const x = centre + s * cos - t * sin;
      const y = centre + s * sin + t * cos;
      sum += sample(volume, n, x, y);
    }
    out[d] = sum * step;
  }
  return out;
}

/**
 * The transpose of {@link projectVolume}: smear a projection back over the slice.
 *
 * Accumulates into `out` rather than replacing it, so a whole set of angles can
 * be summed in one pass.
 *
 * @param {ArrayLike<number>} line n
 * @param {number} n
 * @param {number} theta radians
 * @param {Float64Array} out n*n, accumulated into
 * @param {number} [scale=1]
 * @param {number} [steps] must match the forward projection's
 * @returns {Float64Array} out
 */
export function backProject(line, n, theta, out, scale = 1, steps = n) {
  const centre = (n - 1) / 2;
  const cos = Math.cos(theta);
  const sin = Math.sin(theta);
  const step = n / steps;
  for (let d = 0; d < n; d++) {
    const s = d - centre;
    const value = line[d] * scale * step;
    if (value === 0) continue;
    for (let k = 0; k < steps; k++) {
      const t = (k + 0.5) * step - n / 2;
      const x = centre + s * cos - t * sin;
      const y = centre + s * sin + t * cos;
      splat(out, n, x, y, value);
    }
  }
  return out;
}

/**
 * Angles spread over half a turn, which is all a parallel-beam measurement has:
 * projecting from behind gives the same line reversed.
 *
 * @param {number} count
 * @param {number} [range=Math.PI] narrow it to model a missing wedge
 * @returns {Float64Array}
 */
export function tiltSeries(count, range = Math.PI) {
  const angles = new Float64Array(count);
  for (let i = 0; i < count; i++) {
    angles[i] = -range / 2 + (range * i) / count;
  }
  return angles;
}

/** Bilinear sample, zero outside. */
function sample(volume, n, x, y) {
  const x0 = Math.floor(x);
  const y0 = Math.floor(y);
  const fx = x - x0;
  const fy = y - y0;
  let sum = 0;
  for (let dx = 0; dx <= 1; dx++) {
    const ix = x0 + dx;
    if (ix < 0 || ix >= n) continue;
    const wx = dx ? fx : 1 - fx;
    for (let dy = 0; dy <= 1; dy++) {
      const iy = y0 + dy;
      if (iy < 0 || iy >= n) continue;
      sum += wx * (dy ? fy : 1 - fy) * volume[ix * n + iy];
    }
  }
  return sum;
}

/** The transpose of {@link sample}: spread a value over the same four points. */
function splat(volume, n, x, y, value) {
  const x0 = Math.floor(x);
  const y0 = Math.floor(y);
  const fx = x - x0;
  const fy = y - y0;
  for (let dx = 0; dx <= 1; dx++) {
    const ix = x0 + dx;
    if (ix < 0 || ix >= n) continue;
    const wx = dx ? fx : 1 - fx;
    for (let dy = 0; dy <= 1; dy++) {
      const iy = y0 + dy;
      if (iy < 0 || iy >= n) continue;
      volume[ix * n + iy] += wx * (dy ? fy : 1 - fy) * value;
    }
  }
}
