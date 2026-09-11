// Reciprocal-space coordinate grids.
//
// Ported from `aberration_utils.py:spatial_frequencies` / `polar_coordinates`,
// the single most duplicated block in the lab's Python notebooks.
//
// CONVENTION, inherited from the Python: grids are **corner-centered**
// (numpy.fft.fftfreq ordering). Only shift to a centered layout at display time,
// with `fftshift` from ./fft.js. Arrays are row-major, indexed `[ix * ny + iy]`,
// matching numpy's `indexing="ij"` meshgrid.

/**
 * Discrete Fourier transform sample frequencies, in numpy's fftfreq order:
 * `[0, 1, ..., n/2-1, -n/2, ..., -1] / (n * d)`.
 *
 * @param {number} n number of samples
 * @param {number} [d=1] sample spacing
 * @returns {Float64Array}
 */
export function fftfreq(n, d = 1) {
  const out = new Float64Array(n);
  const scale = 1 / (n * d);
  const half = Math.floor((n - 1) / 2) + 1;
  for (let i = 0; i < half; i++) out[i] = i * scale;
  for (let i = half; i < n; i++) out[i] = (i - n) * scale;
  return out;
}

/**
 * 2D spatial frequency grids kx, ky [1/Angstrom].
 *
 * @param {[number, number]} gpts grid points (nx, ny)
 * @param {[number, number]} sampling real-space sampling (sx, sy) [Angstrom]
 * @returns {{kx: Float64Array, ky: Float64Array, nx: number, ny: number}} flattened row-major grids
 */
export function spatialFrequencies([nx, ny], [sx, sy]) {
  const kxv = fftfreq(nx, sx);
  const kyv = fftfreq(ny, sy);
  const kx = new Float64Array(nx * ny);
  const ky = new Float64Array(nx * ny);
  for (let ix = 0; ix < nx; ix++) {
    for (let iy = 0; iy < ny; iy++) {
      kx[ix * ny + iy] = kxv[ix];
      ky[ix * ny + iy] = kyv[iy];
    }
  }
  return {kx, ky, nx, ny};
}

/**
 * Polar coordinates of a Cartesian frequency grid.
 *
 * @returns {{k: Float64Array, phi: Float64Array}} radius and azimuth [rad]
 */
export function polarCoordinates(kx, ky) {
  const n = kx.length;
  const k = new Float64Array(n);
  const phi = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    k[i] = Math.hypot(kx[i], ky[i]);
    phi[i] = Math.atan2(ky[i], kx[i]);
  }
  return {k, phi};
}

/**
 * Scattering angles alpha = k * lambda [rad] and azimuth phi, the natural
 * coordinates for the aberration function.
 *
 * @param {[number, number]} gpts
 * @param {[number, number]} sampling [Angstrom]
 * @param {number} wavelength [Angstrom] — see units.electronWavelength
 */
export function angularSpatialFrequencies(gpts, sampling, wavelength) {
  const {kx, ky, nx, ny} = spatialFrequencies(gpts, sampling);
  const {k, phi} = polarCoordinates(kx, ky);
  const alpha = new Float64Array(k.length);
  for (let i = 0; i < k.length; i++) alpha[i] = k[i] * wavelength;
  return {alpha, phi, k, kx, ky, nx, ny};
}

/**
 * Angular sampling (the scattering-angle step per pixel) in mrad, as expected by
 * `optics.softAperture`.
 *
 * @returns {[number, number]}
 */
export function angularSampling([nx, ny], [sx, sy], wavelength) {
  return [(wavelength / (nx * sx)) * 1e3, (wavelength / (ny * sy)) * 1e3];
}

/**
 * A checkerboard test pattern, used by the geometric-aberration widget to make
 * ray displacement visible.
 *
 * Ported from `aberration_utils.py:make_checkerboard`.
 */
export function checkerboard([nx, ny], nBlocks = 16) {
  const block = Math.floor(nx / nBlocks);
  const out = new Float64Array(nx * ny);
  const shift = block >> 1;
  for (let ix = 0; ix < nx; ix++) {
    for (let iy = 0; iy < ny; iy++) {
      // np.roll by block//2 along both axes
      const sx = (ix - shift + nx) % nx;
      const sy = (iy - shift + ny) % ny;
      const value = (Math.floor(sx / block) + Math.floor(sy / block)) % 2 === 0 ? 1 : 0;
      out[ix * ny + iy] = value;
    }
  }
  return out;
}
