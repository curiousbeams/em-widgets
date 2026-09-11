// Fast Fourier transforms for arbitrary sizes.
//
// Radix-2 Cooley-Tukey for power-of-two lengths, Bluestein's chirp-z algorithm
// otherwise — the lab's arrays are routinely 244x242, 192x192 or 96x96, so
// non-power-of-two support is not optional.
//
// Complex data is a pair of Float64Arrays, `{re, im}`, row-major for 2D
// (`[ix * ny + iy]`). Normalisation matches numpy: `fft` is unnormalised,
// `ifft` divides by n. Output is corner-centered, like numpy — use `fftshift`
// only when you are about to draw.

/** Allocate a zeroed complex array of length n. */
export function complex(n) {
  return {re: new Float64Array(n), im: new Float64Array(n)};
}

/** Wrap a real array as complex (copies). */
export function toComplex(re) {
  return {re: Float64Array.from(re), im: new Float64Array(re.length)};
}

/** Elementwise magnitude. */
export function abs({re, im}) {
  const out = new Float64Array(re.length);
  for (let i = 0; i < re.length; i++) out[i] = Math.hypot(re[i], im[i]);
  return out;
}

/** Elementwise phase [rad]. */
export function angle({re, im}) {
  const out = new Float64Array(re.length);
  for (let i = 0; i < re.length; i++) out[i] = Math.atan2(im[i], re[i]);
  return out;
}

/** Elementwise product, into a new array (or `out` if given). */
export function multiply(a, b, out = complex(a.re.length)) {
  for (let i = 0; i < a.re.length; i++) {
    const ar = a.re[i], ai = a.im[i], br = b.re[i], bi = b.im[i];
    out.re[i] = ar * br - ai * bi;
    out.im[i] = ar * bi + ai * br;
  }
  return out;
}

/** `exp(i * phase)` for a real phase array. */
export function expi(phase) {
  const n = phase.length;
  const out = complex(n);
  for (let i = 0; i < n; i++) {
    out.re[i] = Math.cos(phase[i]);
    out.im[i] = Math.sin(phase[i]);
  }
  return out;
}

const isPow2 = (n) => n > 0 && (n & (n - 1)) === 0;

// ---------------------------------------------------------------------------
// 1D
// ---------------------------------------------------------------------------

/**
 * In-place 1D FFT of a whole Float64Array pair.
 * @param {boolean} [inverse=false] if true, computes the unscaled inverse
 *   transform (conjugate convention); `ifft1d` applies the 1/n.
 */
export function fft1d(re, im, inverse = false) {
  const n = re.length;
  if (n <= 1) return;
  if (inverse) {
    for (let i = 0; i < n; i++) im[i] = -im[i];
    transform(re, im, n, 0, 1);
    for (let i = 0; i < n; i++) im[i] = -im[i];
  } else {
    transform(re, im, n, 0, 1);
  }
}

/** In-place normalised 1D inverse FFT. */
export function ifft1d(re, im) {
  const n = re.length;
  fft1d(re, im, true);
  for (let i = 0; i < n; i++) {
    re[i] /= n;
    im[i] /= n;
  }
}

/**
 * Forward transform of a strided segment: elements at `offset + j * stride`
 * for j in [0, n). Dispatches to radix-2 or Bluestein.
 */
function transform(re, im, n, offset, stride) {
  if (n <= 1) return;
  if (isPow2(n)) radix2(re, im, n, offset, stride);
  else bluestein(re, im, n, offset, stride);
}

// Twiddle factors and the bit-reversal permutation depend only on the length,
// and a 2D transform runs hundreds of same-length passes — so compute them once.
const radix2Cache = new Map();

function radix2Plan(n) {
  let plan = radix2Cache.get(n);
  if (plan) return plan;

  const half = n >> 1;
  const cos = new Float64Array(half);
  const sin = new Float64Array(half);
  for (let k = 0; k < half; k++) {
    const a = (-2 * Math.PI * k) / n;
    cos[k] = Math.cos(a);
    sin[k] = Math.sin(a);
  }

  const levels = Math.log2(n) | 0;
  const rev = new Uint32Array(n);
  for (let i = 0; i < n; i++) {
    let j = 0;
    for (let b = 0; b < levels; b++) j = (j << 1) | ((i >>> b) & 1);
    rev[i] = j;
  }

  plan = {cos, sin, rev};
  radix2Cache.set(n, plan);
  return plan;
}

/** Iterative in-place radix-2 Cooley-Tukey, decimation in time. */
function radix2(re, im, n, offset, stride) {
  const {cos, sin, rev} = radix2Plan(n);

  for (let i = 0; i < n; i++) {
    const j = rev[i];
    if (j > i) {
      const a = offset + i * stride;
      const c = offset + j * stride;
      let t = re[a]; re[a] = re[c]; re[c] = t;
      t = im[a]; im[a] = im[c]; im[c] = t;
    }
  }

  for (let size = 2; size <= n; size <<= 1) {
    const half = size >> 1;
    const twiddleStep = n / size; // index into the length-n twiddle table
    for (let i = 0; i < n; i += size) {
      for (let j = 0, k = 0; j < half; j++, k += twiddleStep) {
        const wr = cos[k];
        const wi = sin[k];
        const a = offset + (i + j) * stride;
        const b = offset + (i + j + half) * stride;
        const xr = re[b] * wr - im[b] * wi;
        const xi = re[b] * wi + im[b] * wr;
        re[b] = re[a] - xr; im[b] = im[a] - xi;
        re[a] += xr;        im[a] += xi;
      }
    }
  }
}

// Everything Bluestein needs for a given length: the chirp, and — crucially —
// the already-transformed convolution kernel, which depends only on n. A 2D
// transform of a 244 x 242 grid makes 486 Bluestein calls, so rebuilding and
// re-transforming that kernel each time was most of the cost.
const bluesteinCache = new Map();

function bluesteinPlan(n) {
  let plan = bluesteinCache.get(n);
  if (plan) return plan;

  let m = 1;
  while (m < 2 * n + 1) m <<= 1;

  const cos = new Float64Array(n);
  const sin = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    // (i*i % 2n) keeps the argument small, so the trig stays accurate for large n
    const a = (Math.PI * ((i * i) % (2 * n))) / n;
    cos[i] = Math.cos(a);
    sin[i] = Math.sin(a);
  }

  // The convolution kernel, transformed once and reused for every call.
  const kr = new Float64Array(m);
  const ki = new Float64Array(m);
  kr[0] = cos[0];
  ki[0] = sin[0];
  for (let i = 1; i < n; i++) {
    kr[i] = kr[m - i] = cos[i];
    ki[i] = ki[m - i] = sin[i];
  }
  radix2(kr, ki, m, 0, 1);

  plan = {m, cos, sin, kr, ki, ar: new Float64Array(m), ai: new Float64Array(m)};
  bluesteinCache.set(n, plan);
  return plan;
}

/**
 * Bluestein's algorithm: an arbitrary-length DFT expressed as a power-of-two
 * convolution. Used for any length that is not a power of two.
 *
 * The scratch buffers live on the plan, so a 2D transform allocates nothing.
 * This is not re-entrant, which is fine — it is only ever called synchronously
 * from `transform`.
 */
function bluestein(re, im, n, offset, stride) {
  const {m, cos, sin, kr, ki, ar, ai} = bluesteinPlan(n);

  for (let i = 0; i < n; i++) {
    const x = re[offset + i * stride];
    const y = im[offset + i * stride];
    ar[i] = x * cos[i] + y * sin[i];
    ai[i] = y * cos[i] - x * sin[i];
  }
  ar.fill(0, n);
  ai.fill(0, n);

  radix2(ar, ai, m, 0, 1);
  for (let i = 0; i < m; i++) {
    const tr = ar[i] * kr[i] - ai[i] * ki[i];
    ai[i] = ar[i] * ki[i] + ai[i] * kr[i];
    ar[i] = tr;
  }

  // inverse transform of the product, via conjugation
  for (let i = 0; i < m; i++) ai[i] = -ai[i];
  radix2(ar, ai, m, 0, 1);

  for (let i = 0; i < n; i++) {
    const xr = ar[i] / m;
    const xi = ai[i] / -m;
    re[offset + i * stride] = xr * cos[i] + xi * sin[i];
    im[offset + i * stride] = xi * cos[i] - xr * sin[i];
  }
}

// ---------------------------------------------------------------------------
// 2D
// ---------------------------------------------------------------------------

/**
 * 2D forward FFT, matching `numpy.fft.fft2`. Operates in place on `{re, im}`
 * of length nx*ny (row-major).
 */
export function fft2({re, im}, nx, ny) {
  for (let ix = 0; ix < nx; ix++) transform(re, im, ny, ix * ny, 1); // rows
  for (let iy = 0; iy < ny; iy++) transform(re, im, nx, iy, ny);     // columns
  return {re, im};
}

/** 2D inverse FFT, matching `numpy.fft.ifft2` (includes the 1/(nx*ny)). */
export function ifft2({re, im}, nx, ny) {
  const n = nx * ny;
  for (let i = 0; i < n; i++) im[i] = -im[i];
  fft2({re, im}, nx, ny);
  for (let i = 0; i < n; i++) {
    re[i] /= n;
    im[i] /= -n;
  }
  return {re, im};
}

/** Non-mutating variants. */
export const fft2Of = (a, nx, ny) => fft2({re: Float64Array.from(a.re), im: Float64Array.from(a.im)}, nx, ny);
export const ifft2Of = (a, nx, ny) => ifft2({re: Float64Array.from(a.re), im: Float64Array.from(a.im)}, nx, ny);

// ---------------------------------------------------------------------------
// Shifts
// ---------------------------------------------------------------------------

/** Move the zero frequency to the centre (numpy.fft.fftshift) for a real 2D array. */
export function fftshift(a, nx, ny) {
  return rollBoth(a, nx, ny, nx >> 1, ny >> 1);
}

/** Inverse of {@link fftshift}. */
export function ifftshift(a, nx, ny) {
  return rollBoth(a, nx, ny, nx - (nx >> 1), ny - (ny >> 1));
}

function rollBoth(a, nx, ny, dx, dy) {
  const out = new (a.constructor)(a.length);
  for (let ix = 0; ix < nx; ix++) {
    const sx = (ix + dx) % nx;
    for (let iy = 0; iy < ny; iy++) {
      out[sx * ny + ((iy + dy) % ny)] = a[ix * ny + iy];
    }
  }
  return out;
}

/** {@link fftshift} for a complex array. */
export const fftshiftComplex = ({re, im}, nx, ny) => ({
  re: fftshift(re, nx, ny),
  im: fftshift(im, nx, ny)
});

/**
 * Fourier translation operator `exp(-2 pi i (kx dx + ky dy))`, for shifting an
 * array by a sub-pixel amount without interpolation.
 *
 * Multiply this into the Fourier transform, then inverse transform.
 */
export function fourierShift(kx, ky, dx, dy) {
  const n = kx.length;
  const out = complex(n);
  for (let i = 0; i < n; i++) {
    const a = -2 * Math.PI * (kx[i] * dx + ky[i] * dy);
    out.re[i] = Math.cos(a);
    out.im[i] = Math.sin(a);
  }
  return out;
}
