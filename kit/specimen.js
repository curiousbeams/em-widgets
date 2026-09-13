// Building a specimen in the browser: atoms, and the projected potential they
// produce.
//
// This replaces loading a precomputed potential. Generating it here means the
// field of view, the particles and the support are parameters rather than fixed
// assets — and it is cheap: splatting a few thousand atoms costs single-digit
// milliseconds.
//
// The potential follows abtem's default exactly — the Lobato parametrisation with
// `projection="infinite"`, i.e. a Fourier-space superposition:
//
//     V(r) = IFFT[ f_proj(|k|) · Σ_j exp(-2πi k·r_j) ] / (dx·dy)
//
// Done in Fourier space rather than by summing real-space atom kernels, for two
// reasons: it is what abtem actually does, so the result can be checked against
// it; and carbon's Lobato coefficients alternate sign with magnitudes of ~30000,
// which makes a naive real-space sum badly conditioned while still looking
// plausible.

import {complex, fft2, ifft2} from "./fft.js";
import {fftfreq} from "./grid.js";

/**
 * Lobato scattering-factor parameters, pre-scaled the way abtem scales them
 * (`parametrizations/__init__.py:540`):
 *
 *     a = π² A / B^1.5 / kappa,   b = 2π / sqrt(B)
 *
 * Hardcoded in scaled form rather than re-deriving `kappa` — which depends on
 * ASE's unit system — so there is one less way to be silently wrong.
 * Extracted from abtem 1.0.9.
 */
export const LOBATO = {
  Au: {
    Z: 79,
    a: [61.0250781385, 875.7999798205, 4305.3050058947, 6464.688669934, 2748.7804540293],
    b: [2.6737416448, 5.3484625576, 15.5996661716, 66.1639107877, 322.6273884016]
  },
  C: {
    Z: 6,
    a: [15611.1645535928, -29746.6186271721, 31465.2296764499, -17242.9598187154, 314.7179248225],
    b: [4.0379763786, 4.1381731713, 4.3899526516, 4.5186105612, 22.6580872319]
  }
};

/**
 * The projected scattering factor, `lobato.py:101`.
 *
 * @param {number} k2 |k|² in 1/Å²
 * @param {{a: number[], b: number[]}} p scaled parameters
 */
export function projectedScatteringFactor(k2, p) {
  const K = 4 * Math.PI * Math.PI * k2;
  let sum = 0;
  for (let i = 0; i < 5; i++) {
    const b = p.b[i];
    const d = K + b * b;
    sum += p.a[i] / b / d + (p.a[i] * b) / (d * d);
  }
  return 8 * Math.PI * sum;
}

/**
 * Splat atoms onto a grid with bilinear weights, wrapping at the edges.
 *
 * The result is a sampled delta comb, so its DFT is the structure factor
 * `Σ_j exp(-2πi k·r_j)`. Bilinear rather than nearest keeps sub-pixel atom
 * positions meaningful.
 */
export function splatDensity(positions, {gpts: [nx, ny], sampling: [sx, sy]}) {
  const out = new Float64Array(nx * ny);
  for (const [x, y] of positions) {
    const fx = x / sx;
    const fy = y / sy;
    const i0 = Math.floor(fx);
    const j0 = Math.floor(fy);
    const tx = fx - i0;
    const ty = fy - j0;
    const ia = ((i0 % nx) + nx) % nx;
    const ib = ((i0 + 1) % nx + nx) % nx;
    const ja = ((j0 % ny) + ny) % ny;
    const jb = ((j0 + 1) % ny + ny) % ny;
    out[ia * ny + ja] += (1 - tx) * (1 - ty);
    out[ib * ny + ja] += tx * (1 - ty);
    out[ia * ny + jb] += (1 - tx) * ty;
    out[ib * ny + jb] += tx * ty;
  }
  return out;
}

/**
 * Projected electrostatic potential of a set of atoms, in V·Å.
 *
 * Matches `abtem.Potential(atoms, gpts, slice_thickness=cell_z,
 * parametrization="lobato")` for a single slice.
 *
 * @param {Array<{symbol: string, positions: Array<[number, number]>}>} groups
 *   atoms grouped by element — one FFT pair per element, not per atom
 * @param {{gpts: [number, number], sampling: [number, number]}} grid
 * @returns {Float64Array} row-major `[ix * ny + iy]`, matching the kit convention
 *   and abtem's `(gpts_x, gpts_y)` array order
 */
export function projectedPotential(groups, {gpts, sampling}) {
  const [nx, ny] = gpts;
  const [sx, sy] = sampling;

  const kx = fftfreq(nx, sx);
  const ky = fftfreq(ny, sy);

  const total = complex(nx * ny);
  for (const group of groups) {
    const params = LOBATO[group.symbol];
    if (!params) throw new Error(`No Lobato parameters for "${group.symbol}".`);
    if (!group.positions.length) continue;

    const density = splatDensity(group.positions, {gpts, sampling});
    const spectrum = fft2({re: density, im: new Float64Array(nx * ny)}, nx, ny);

    for (let ix = 0; ix < nx; ix++) {
      for (let iy = 0; iy < ny; iy++) {
        const i = ix * ny + iy;
        const f = projectedScatteringFactor(kx[ix] * kx[ix] + ky[iy] * ky[iy], params);
        const g = f * sincCorrection(kx[ix] * sx, ky[iy] * sy);
        total.re[i] += spectrum.re[i] * g;
        total.im[i] += spectrum.im[i] * g;
      }
    }
  }

  ifft2(total, nx, ny);
  // The discrete inverse transform carries 1/(nx·ny); the continuous one needs
  // the reciprocal-space volume element, which leaves a 1/(dx·dy). That factor
  // is abtem's `dk2` inside its sinc, kept here for clarity.
  const scale = 1 / (sx * sy);
  const out = new Float64Array(nx * ny);
  for (let i = 0; i < out.length; i++) out[i] = total.re[i] * scale;
  return out;
}

/**
 * Undo the smoothing that splatting a delta across neighbouring pixels causes.
 *
 * Bilinear splatting convolves each atom with a pixel-scale kernel, which damps
 * high spatial frequencies — 3.5% by k = 2.25 1/Å at 0.2 Å sampling, which shows
 * up as atom peaks that are systematically too low.
 *
 * This is abtem's correction verbatim (`integrals.py:299`): `sin(k)/k` of the
 * *combined* magnitude `|(kx·dx, ky·dy)|`. Note it is `sin(k)/k`, not the
 * unnormalised `sin(πk)/(πk)` that the exact bilinear kernel transform would
 * give — matching abtem matters more here than deriving it independently.
 */
function sincCorrection(kxdx, kydy) {
  const k = Math.hypot(kxdx, kydy);
  if (k === 0) return 1;
  return k / Math.sin(k);
}

// ---------------------------------------------------------------------------
// Atomic models
// ---------------------------------------------------------------------------

/** Deterministic RNG, so a specimen is reproducible from its seed. */
export function mulberry32(seed) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** Rotate a vector about an axis by `angle` radians (Rodrigues). */
function rotateAbout([x, y, z], axis, angle) {
  const n = Math.hypot(...axis) || 1;
  const [ux, uy, uz] = axis.map((v) => v / n);
  const c = Math.cos(angle), s = Math.sin(angle), t = 1 - c;
  return [
    (t * ux * ux + c) * x + (t * ux * uy - s * uz) * y + (t * ux * uz + s * uy) * z,
    (t * ux * uy + s * uz) * x + (t * uy * uy + c) * y + (t * uy * uz - s * ux) * z,
    (t * ux * uz - s * uy) * x + (t * uy * uz + s * ux) * y + (t * uz * uz + c) * z
  ];
}

/** The rotation taking `from` onto `to`, as used to put a zone axis along z. */
function alignRotation(from, to) {
  const nf = Math.hypot(...from), nt = Math.hypot(...to);
  const f = from.map((v) => v / nf);
  const t = to.map((v) => v / nt);
  const axis = [
    f[1] * t[2] - f[2] * t[1],
    f[2] * t[0] - f[0] * t[2],
    f[0] * t[1] - f[1] * t[0]
  ];
  const dot = f[0] * t[0] + f[1] * t[1] + f[2] * t[2];
  if (Math.hypot(...axis) < 1e-12) return dot > 0 ? null : {axis: [1, 0, 0], angle: Math.PI};
  return {axis, angle: Math.acos(Math.max(-1, Math.min(1, dot)))};
}

/**
 * An FCC nanoparticle, cut by families of crystal planes.
 *
 * Parameterised by **diameter** rather than by ASE's layer counts, because the
 * diameter is what a reader wants on a slider. The `ratios` are relative plane
 * distances in the spirit of a Wulff construction — `sim.py`'s decahedron uses
 * {100}:1.25, {111}:1.0, {110}:1.05 — and are normalised so the largest sets the
 * requested diameter.
 *
 * @param {object} options
 * @param {number} [options.diameter=22] Å, across the widest direction
 * @param {number} [options.latticeConstant=4.08] Å (gold)
 * @param {Array<[number,number,number]>} [options.surfaces]
 * @param {number[]} [options.ratios] relative truncation per surface family
 * @param {[number,number,number]} [options.zoneAxis] crystal direction along the beam
 * @param {number} [options.twist=0] rotation about the beam, radians
 * @param {[number,number,number]} [options.centre=[0,0,0]]
 * @returns {Array<[number, number, number]>} atom positions
 */
/** Every distinct permutation and sign of a Miller index — the family {hkl}. */
function family(hkl) {
  const seen = new Map();
  const order = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
  for (const p of order) {
    for (const sx of [1, -1]) for (const sy of [1, -1]) for (const sz of [1, -1]) {
      const n = [hkl[p[0]] * sx, hkl[p[1]] * sy, hkl[p[2]] * sz];
      seen.set(n.join(","), n);
    }
  }
  return [...seen.values()];
}

export function fccCluster({
  diameter = 22,
  latticeConstant = 4.08,
  surfaces = [[1, 0, 0], [1, 1, 0], [1, 1, 1]],
  ratios = [1.25, 1.05, 1.0],
  zoneAxis = null,
  twist = 0,
  centre = [0, 0, 0]
} = {}) {
  const scale = diameter / 2 / Math.max(...ratios);

  const planes = [];
  for (let s = 0; s < surfaces.length; s++) {
    const hkl = surfaces[s];
    const norm = Math.hypot(...hkl);
    const distance = ratios[s] * scale;
    // A surface *family* {hkl} is every permutation of the indices as well as
    // every sign: {100} is six planes, {110} twelve, {111} eight. Taking only the
    // signs leaves {100} bounding ±x alone, and the particle comes out as a
    // double pyramid capped by {111} rather than a truncated octahedron.
    for (const n of family(hkl)) {
      planes.push({n: n.map((v) => v / norm), distance});
    }
  }

  const reach = Math.ceil(Math.max(...planes.map((p) => p.distance)) / latticeConstant) + 1;
  const basis = [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]];

  let positions = [];
  for (let i = -reach; i <= reach; i++) {
    for (let j = -reach; j <= reach; j++) {
      for (let k = -reach; k <= reach; k++) {
        for (const b of basis) {
          const p = [
            (i + b[0]) * latticeConstant,
            (j + b[1]) * latticeConstant,
            (k + b[2]) * latticeConstant
          ];
          let inside = true;
          for (const plane of planes) {
            if (p[0] * plane.n[0] + p[1] * plane.n[1] + p[2] * plane.n[2] > plane.distance + 1e-9) {
              inside = false;
              break;
            }
          }
          if (inside) positions.push(p);
        }
      }
    }
  }

  if (zoneAxis) {
    const r = alignRotation(zoneAxis, [0, 0, 1]);
    if (r) positions = positions.map((p) => rotateAbout(p, r.axis, r.angle));
  }
  if (twist) positions = positions.map((p) => rotateAbout(p, [0, 0, 1], twist));
  return positions.map((p) => [p[0] + centre[0], p[1] + centre[1], p[2] + centre[2]]);
}

/**
 * An amorphous carbon support.
 *
 * Following `sim.py:43`, a diamond lattice jittered by half a bond length —
 * which gives a diffuse ring and no Bragg reflections, the point being that the
 * support is invisible to anything working from detected reflections.
 *
 * One deliberate departure: `sim.py` folds four diamond cells into an 8 Å slab
 * (`z % 8 + 33`), which doubles the density to 0.346 atoms/Å³, about twice
 * diamond and three times real amorphous carbon. Here the atom count comes from
 * a physical areal density instead.
 *
 * @param {object} options
 * @param {[number, number]} options.extent field of view, Å
 * @param {number} [options.thickness=8] Å
 * @param {number} [options.arealDensity=1.1] C atoms per Å² — ~2 g/cm³ at 8 Å
 * @param {number} [options.z=0] centre of the support along the beam
 * @param {[number, number]} [options.origin=[0, 0]] in-plane offset, for laying
 *   extra film outside the simulated field of view
 * @param {number} [options.seed=10]
 */
export function amorphousSupport({
  extent,
  thickness = 8,
  arealDensity = 1.1,
  z = 0,
  origin = [0, 0],
  seed = 10
} = {}) {
  const [Lx, Ly] = extent;
  const a = 3.567;                       // diamond lattice constant
  const jitter = 0.5 * 1.54;             // half a C–C bond, as in sim.py
  const random = mulberry32(seed);
  const gauss = () => {
    // Box–Muller; the pair's second value is discarded for simplicity.
    const u = Math.max(random(), Number.MIN_VALUE);
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * random());
  };

  const basis = [
    [0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0],
    [0.25, 0.25, 0.25], [0.25, 0.75, 0.75], [0.75, 0.25, 0.75], [0.75, 0.75, 0.25]
  ];
  const nx = Math.ceil(Lx / a);
  const ny = Math.ceil(Ly / a);
  const nz = Math.max(1, Math.round(thickness / a));

  const candidates = [];
  for (let i = 0; i < nx; i++) {
    for (let j = 0; j < ny; j++) {
      for (let k = 0; k < nz; k++) {
        for (const b of basis) {
          candidates.push([
            (i + b[0]) * a + gauss() * jitter,
            (j + b[1]) * a + gauss() * jitter,
            (k + b[2]) * a + gauss() * jitter
          ]);
        }
      }
    }
  }

  // Thin to the requested areal density rather than keeping the lattice count.
  const wanted = Math.round(arealDensity * Lx * Ly);
  const keep = Math.min(1, wanted / candidates.length);
  const out = [];
  for (const [x, y, zz] of candidates) {
    if (random() > keep) continue;
    out.push([
      (((x % Lx) + Lx) % Lx) + origin[0],
      (((y % Ly) + Ly) % Ly) + origin[1],
      z + (((zz % thickness) + thickness) % thickness) - thickness / 2
    ]);
  }
  return out;
}
