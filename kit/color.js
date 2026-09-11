// Colour: perceptually-uniform complex-domain colouring, and the colormaps the
// lab's figures use.
//
// `complexToRGB` is the single most reused primitive in the Python notebooks
// (~15 widgets) and had three separate implementations: `ctf/visualize.py`,
// `py4DSTEM.visualize.Complex2RGB`, and quantem's `array_to_rgba`. All three
// compute the same thing:
//
//     J = amplitude * 61.5            (luminance, capped at the monotonic chroma cutoff)
//     C = min(98 * J / 123, 110)      (chroma)
//     h = degrees(phase) + 180        (hue)
//
// and then convert JCh -> sRGB. In Python that last step is
// `colorspacious.cspace_convert(JCh, "JCh", "sRGB1")`, i.e. an inverse CIECAM02
// under colorspacious's sRGB viewing conditions followed by the sRGB transfer
// function. That conversion is reimplemented below and is checked against
// colorspacious 1.1.2 in test/color.test.js.

import {histogramScaling} from "./image.js";
import {COLORMAP_DATA} from "./colormap-data.js";

// ---------------------------------------------------------------------------
// Inverse CIECAM02 (colorspacious's CIECAM02Space.sRGB viewing conditions)
// ---------------------------------------------------------------------------

// Precomputed appearance constants for XYZ100_w = D65, Y_b = 20, L_A = 64/pi/5,
// average surround. Printed directly from colorspacious rather than re-derived,
// so the JS cannot drift from the Python.
const VC = {
  c: 0.69,
  N_c: 1.0,
  F_L: 0.27313053667320736,
  n: 0.2,
  z: 1.9272135954999579,
  N_bb: 1.0003040045593807,
  N_cb: 1.0003040045593807,
  A_w: 25.515986669780446,
  D_RGB: [1.0444367755899502, 0.9715728131498355, 0.9332914123252806]
};

const M_CAT02 = [
  [0.7328, 0.4296, -0.1624],
  [-0.7036, 1.6975, 0.0061],
  [0.003, 0.0136, 0.9834]
];
const M_HPE = [
  [0.38971, 0.68898, -0.07868],
  [-0.22981, 1.1834, 0.04641],
  [0.0, 0.0, 1.0]
];

function matMul(A, B) {
  return A.map((row) => B[0].map((_, j) => row.reduce((s, v, k) => s + v * B[k][j], 0)));
}
function matInv3([[a, b, c], [d, e, f], [g, h, i]]) {
  const det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
  return [
    [(e * i - f * h) / det, (c * h - b * i) / det, (b * f - c * e) / det],
    [(f * g - d * i) / det, (a * i - c * g) / det, (c * d - a * f) / det],
    [(d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det]
  ];
}
const matVec = (M, v) => M.map((row) => row[0] * v[0] + row[1] * v[1] + row[2] * v[2]);

const M_CAT02_inv = matInv3(M_CAT02);
const M_CAT02_M_HPE_inv = matMul(M_CAT02, matInv3(M_HPE));

// colorspacious rounds the XYZ100/100 -> linear sRGB matrix to 4 decimals;
// match it exactly so the two implementations agree to float precision.
const M_XYZ_TO_LINEAR = [
  [3.2406, -1.5372, -0.4986],
  [-0.9689, 1.8758, 0.0415],
  [0.0557, -0.204, 1.057]
];

/** The sRGB transfer function (its linear segment is extended to negatives). */
function srgbGamma(v) {
  return v <= 0.0031308 ? 12.92 * v : 1.055 * Math.pow(v, 1 / 2.4) - 0.055;
}

/**
 * One CIECAM02 JCh triplet to (possibly out-of-gamut) sRGB.
 *
 * @param {number} J lightness, 0..100
 * @param {number} C chroma
 * @param {number} h hue angle [degrees]
 * @returns {[number, number, number]} sRGB, nominally 0..1 but NOT clipped
 */
export function jchToSrgb(J, C, h) {
  const {c, N_c, F_L, n, z, N_bb, N_cb, A_w, D_RGB} = VC;

  // J = 0 is black. Guarding here also avoids the 0/0 in `t` below, which
  // numpy silently produces as NaN.
  if (!(J > 0)) return [0, 0, 0];

  // Step 2
  const t = Math.pow(C / (Math.sqrt(J / 100) * Math.pow(1.64 - Math.pow(0.29, n), 0.73)), 1 / 0.9);
  const hRad = (h * Math.PI) / 180;
  const e_t = 0.25 * (Math.cos(hRad + 2) + 3.8);
  const A = A_w * Math.pow(J / 100, 1 / (c * z));

  const oneOverT = t === 0 ? Infinity : 1 / t;
  const p_1 = (50000 / 13) * N_c * N_cb * e_t * oneOverT;
  const p_2 = A / N_bb + 0.305;
  const p_3 = 21 / 20;

  // Step 3 — branch on |sin| vs |cos| to dodge the divide-by-zero
  const sinH = Math.sin(hRad);
  const cosH = Math.cos(hRad);
  const num = p_2 * (2 + p_3) * (460 / 1403);
  const denom2 = (2 + p_3) * (220 / 1403);
  const denom3 = -27 / 1403 + p_3 * (6300 / 1403);

  let a, b;
  if (Math.abs(sinH) >= Math.abs(cosH)) {
    b = num / (p_1 / sinH + denom2 * (cosH / sinH) + denom3);
    a = (b * cosH) / sinH;
  } else {
    a = num / (p_1 / cosH + denom2 + denom3 * (sinH / cosH));
    b = (a * sinH) / cosH;
  }

  // Step 4
  const p2ab = [p_2, a, b];
  const RGBprime_a = matVec(
    [
      [460 / 1403, 451 / 1403, 288 / 1403],
      [460 / 1403, -891 / 1403, -261 / 1403],
      [460 / 1403, -220 / 1403, -6300 / 1403]
    ],
    p2ab
  );

  // Step 5 — invert the nonlinear response compression
  const RGBprime = RGBprime_a.map((v) => {
    const d = v - 0.1;
    const abs = Math.abs(d);
    return Math.sign(d) * (100 / F_L) * Math.pow((27.13 * abs) / (400 - abs), 1 / 0.42);
  });

  // Steps 6-8
  const RGB_C = matVec(M_CAT02_M_HPE_inv, RGBprime);
  const RGB = [RGB_C[0] / D_RGB[0], RGB_C[1] / D_RGB[1], RGB_C[2] / D_RGB[2]];
  const XYZ100 = matVec(M_CAT02_inv, RGB);

  const linear = matVec(M_XYZ_TO_LINEAR, [XYZ100[0] / 100, XYZ100[1] / 100, XYZ100[2] / 100]);
  return [srgbGamma(linear[0]), srgbGamma(linear[1]), srgbGamma(linear[2])];
}

// ---------------------------------------------------------------------------
// Complex domain colouring
// ---------------------------------------------------------------------------

/**
 * Convert a complex array to a perceptually-uniform RGBA image.
 *
 * Amplitude drives lightness, phase drives hue — the same mapping as
 * `ctf/visualize.py:complex_to_rgb` and `py4DSTEM.visualize.Complex2RGB`.
 *
 * @param {{re: Float64Array, im: Float64Array}} array
 * @param {object} [options]
 * @param {number} [options.vmin=0.02] lower quantile for amplitude scaling
 * @param {number} [options.vmax=0.98] upper quantile
 * @param {number} [options.power=1] raise the scaled amplitude to this power
 * @param {number} [options.chromaBoost=1] multiply chroma (py4DSTEM's chroma_boost)
 * @returns {Uint8ClampedArray} RGBA, 4 bytes per pixel, ready for ImageData
 */
export function complexToRGB(array, {vmin = 0.02, vmax = 0.98, power = 1, chromaBoost = 1} = {}) {
  const {re, im} = array;
  const n = re.length;

  const magnitude = new Float64Array(n);
  for (let i = 0; i < n; i++) magnitude[i] = Math.hypot(re[i], im[i]);

  const amp = histogramScaling(magnitude, {vmin, vmax, normalize: true});
  const out = new Uint8ClampedArray(n * 4);

  for (let i = 0; i < n; i++) {
    let a = Math.min(Math.max(amp[i], 1e-16), 1);
    if (power !== 1) a = Math.pow(a, power);
    const J = a * 61.5;
    const C = Math.min((98 * J) / 123, 110) * chromaBoost;
    const h = (Math.atan2(im[i], re[i]) * 180) / Math.PI + 180;
    const [r, g, bl] = jchToSrgb(J, C, h);
    out[i * 4] = clamp01(r) * 255;
    out[i * 4 + 1] = clamp01(g) * 255;
    out[i * 4 + 2] = clamp01(bl) * 255;
    out[i * 4 + 3] = 255;
  }
  return out;
}

const clamp01 = (v) => (v < 0 ? 0 : v > 1 ? 1 : v);

/**
 * The hue wheel used as a legend for {@link complexToRGB}.
 *
 * @param {number} size pixel width/height of the (square) wheel
 * @returns {Uint8ClampedArray} RGBA; pixels outside the disc are transparent
 */
export function phaseWheel(size = 128) {
  const out = new Uint8ClampedArray(size * size * 4);
  const r0 = size / 2;
  for (let y = 0; y < size; y++) {
    for (let x = 0; x < size; x++) {
      const dx = x - r0 + 0.5;
      const dy = y - r0 + 0.5;
      const r = Math.hypot(dx, dy) / r0;
      const i = (y * size + x) * 4;
      if (r > 1) continue;
      const J = r * 61.5;
      const C = Math.min((98 * J) / 123, 110);
      const h = (Math.atan2(dy, dx) * 180) / Math.PI + 180;
      const [cr, cg, cb] = jchToSrgb(J, C, h);
      out[i] = clamp01(cr) * 255;
      out[i + 1] = clamp01(cg) * 255;
      out[i + 2] = clamp01(cb) * 255;
      out[i + 3] = 255;
    }
  }
  return out;
}

// ---------------------------------------------------------------------------
// Colormaps
// ---------------------------------------------------------------------------

function decodeColormap(b64) {
  const bin = typeof atob === "function" ? atob(b64) : Buffer.from(b64, "base64").toString("binary");
  const out = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
  return out;
}

const cache = new Map();

/**
 * A 256-entry RGB lookup table, as a flat Uint8Array of length 768.
 *
 * Available: magma, gray, twilight, RdBu, PuOr, PiYG, eclipse.
 * These match matplotlib / cmasher so JS figures line up with the Python ones.
 */
export function colormap(name) {
  if (cache.has(name)) return cache.get(name);
  let table;
  if (name === "gray" || name === "grey") {
    table = new Uint8Array(768);
    for (let i = 0; i < 256; i++) table[i * 3] = table[i * 3 + 1] = table[i * 3 + 2] = i;
  } else {
    const data = COLORMAP_DATA[name];
    if (!data) {
      throw new Error(`Unknown colormap "${name}". Available: ${listColormaps().join(", ")}.`);
    }
    table = decodeColormap(data);
  }
  cache.set(name, table);
  return table;
}

/** Names accepted by {@link colormap}. */
export function listColormaps() {
  return ["gray", ...Object.keys(COLORMAP_DATA)];
}

/**
 * Map a real array to RGBA using a named colormap.
 *
 * @param {ArrayLike<number>} array
 * @param {object} [options]
 * @param {string} [options.colormap="magma"]
 * @param {number} [options.vmin] value mapped to 0; defaults to the array min
 * @param {number} [options.vmax] value mapped to 1; defaults to the array max
 * @param {number} [options.power=1] applied to the normalised value
 * @returns {Uint8ClampedArray} RGBA
 */
export function applyColormap(array, {colormap: name = "magma", vmin, vmax, power = 1} = {}) {
  const table = colormap(name);
  const n = array.length;

  if (vmin === undefined || vmax === undefined) {
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < n; i++) {
      const v = array[i];
      if (Number.isNaN(v)) continue;
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    }
    vmin ??= lo;
    vmax ??= hi;
  }
  const span = vmax - vmin || 1;

  const out = new Uint8ClampedArray(n * 4);
  for (let i = 0; i < n; i++) {
    let t = (array[i] - vmin) / span;
    t = t < 0 ? 0 : t > 1 ? 1 : t;
    if (power !== 1) t = Math.pow(t, power);
    const j = Math.min(255, Math.round(t * 255)) * 3;
    out[i * 4] = table[j];
    out[i * 4 + 1] = table[j + 1];
    out[i * 4 + 2] = table[j + 2];
    out[i * 4 + 3] = 255;
  }
  return out;
}
