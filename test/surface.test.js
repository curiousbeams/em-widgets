import test from "node:test";
import assert from "node:assert/strict";
import {readFileSync} from "node:fs";

import {complex} from "../kit/fft.js";
import {
  reciprocalVectors, latticeOrigins, splatAtoms, beamWindow, diffractionPattern,
  illuminatedPattern, superstructureSpots, spotSum
} from "../kit/surface.js";

// The real Si(111)-(7x7) cell, as `tools/si111-7x7.py` writes it.
const SI = JSON.parse(readFileSync(new URL("../data/si111-7x7.json", import.meta.url)));

const N = 192;
const EXTENT = 160;          // Angstrom; about six 7x7 cells across
const BEAM = 12;

/** Where a cell's content actually sits: its origin plus half a diagonal. */
const centreOf = (origin) => ({
  x: origin.x + (SI.a1[0] + SI.a2[0]) / 2,
  y: origin.y + (SI.a1[1] + SI.a2[1]) / 2
});

/** Tile one of the two bases over the field, choosing per cell. */
function field(choose) {
  const atoms = [];
  for (const origin of latticeOrigins(SI.a1, SI.a2, EXTENT, 60)) {
    const basis = choose(centreOf(origin)) ? SI.das : SI.bulk;
    for (const atom of basis) atoms.push([origin.x + atom[0], origin.y + atom[1], atom[2]]);
  }
  return splatAtoms({n: N, extent: EXTENT, atoms});
}

function patternAt(density, x, y) {
  const window = beamWindow({n: N, extent: EXTENT, x, y, sigma: BEAM});
  return diffractionPattern(density, window, N, complex(N * N));
}

const {b1, b2} = reciprocalVectors(SI.a1, SI.a2);
const SPOTS = superstructureSpots(b1, b2, {multiple: 7});

/** The share of the intensity in the spots the reconstruction alone produces. */
function reconstructionSignal(pattern) {
  const fractional = spotSum(pattern, N, EXTENT, SPOTS.fractional);
  const integer = spotSum(pattern, N, EXTENT, SPOTS.integer);
  return fractional / (fractional + integer);
}

test("reciprocal vectors are dual to the lattice", () => {
  const a1 = SI.a1;
  const a2 = SI.a2;
  assert.ok(Math.abs(a1[0] * b1[0] + a1[1] * b1[1] - 1) < 1e-12);
  assert.ok(Math.abs(a2[0] * b2[0] + a2[1] * b2[1] - 1) < 1e-12);
  assert.ok(Math.abs(a1[0] * b2[0] + a1[1] * b2[1]) < 1e-12);
  assert.ok(Math.abs(a2[0] * b1[0] + a2[1] * b1[1]) < 1e-12);
  // A 7 x 7 cell 26.9 A across: its own spots are at 1/26.9, the substrate's
  // seven times further out.
  assert.ok(Math.abs(Math.hypot(...b1) - 1 / 26.877 / Math.sin(Math.PI / 3)) < 1e-4);
});

test("the lattice covers the field, and only just", () => {
  const origins = latticeOrigins(SI.a1, SI.a2, EXTENT, 30);
  const inside = origins.filter((o) => o.x >= 0 && o.x < EXTENT && o.y >= 0 && o.y < EXTENT);
  // Six cells of 26.877 A across a 160 A field, so of order six squared inside.
  assert.ok(inside.length >= 25 && inside.length <= 45, `${inside.length} origins inside`);
  // Every corner of the field is covered by some cell.
  for (const [x, y] of [[0, 0], [EXTENT, 0], [0, EXTENT], [EXTENT, EXTENT]]) {
    const near = origins.some((o) => Math.hypot(o.x - x, o.y - y) < 27);
    assert.ok(near, `no cell near (${x}, ${y})`);
  }
});

test("the seven-fold cell puts six fractional spots between every substrate pair", () => {
  // Along one axis: indices 1..6 are fractional and 7 is the substrate's own.
  assert.equal(SPOTS.integer.length, 6);            // the first substrate ring
  assert.ok(SPOTS.fractional.length > 100, `${SPOTS.fractional.length} fractional spots`);
  // The innermost ring is left out: it sits under the direct beam's own skirt.
  const closest = Math.min(...SPOTS.fractional.map((q) => Math.hypot(...q)));
  assert.ok(closest > 1.4 * Math.hypot(...b1), `closest fractional spot at ${closest}`);
  // To a part in ten thousand, which is where the cell's own file rounds off.
  for (const q of SPOTS.integer) {
    assert.ok(Math.abs(Math.hypot(...q) / (7 * Math.hypot(...b1)) - 1) < 1e-4);
  }
});

test("a bulk-terminated surface has no fractional spots; the reconstruction does", () => {
  const bulk = field(() => false);
  const reconstructed = field(() => true);
  const middle = EXTENT / 2;

  const bulkSignal = reconstructionSignal(patternAt(bulk, middle, middle));
  const dasSignal = reconstructionSignal(patternAt(reconstructed, middle, middle));

  assert.ok(bulkSignal < 0.05, `bulk should be quiet between the orders, got ${bulkSignal}`);
  assert.ok(dasSignal > 0.8, `the reconstruction should be loud there, got ${dasSignal}`);
});

test("a focused beam reads back which side of a boundary it is on", () => {
  // Reconstructed on one side of the field, bulk-terminated on the other.
  const density = field((centre) => centre.x < EXTENT / 2);
  const profile = [];
  for (let x = 20; x <= 140; x += 5) {
    profile.push(reconstructionSignal(patternAt(density, x, EXTENT / 2)));
  }
  assert.ok(profile[0] > 0.8, `reconstructed side, got ${profile[0]}`);
  assert.ok(profile.at(-1) < 0.15, `bulk side, got ${profile.at(-1)}`);
  // And the boundary is resolved rather than falling between two samples: the
  // beam is wider than the scan step, so something has to read as half.
  assert.ok(profile.some((v) => v > 0.2 && v < 0.8),
    `no partly-covered position in ${profile.map((v) => v.toFixed(2)).join(" ")}`);
});

test("a beam too wide to resolve the boundary averages across it", () => {
  const density = field((centre) => centre.x < EXTENT / 2);
  const wide = beamWindow({n: N, extent: EXTENT, x: EXTENT * 0.2, y: EXTENT / 2, sigma: 45});
  const spread = reconstructionSignal(diffractionPattern(density, wide, N, complex(N * N)));
  const focused = reconstructionSignal(patternAt(density, EXTENT * 0.2, EXTENT / 2));
  assert.ok(spread < focused, `a wide beam should dilute the signal: ${spread} vs ${focused}`);
});

test("transforming the illuminated patch says the same as transforming the field", () => {
  const density = field((centre) => centre.x < EXTENT / 2);
  for (const [x, y] of [[40, 80], [80, 80], [120, 60]]) {
    const whole = reconstructionSignal(patternAt(density, x, y));

    const patch = illuminatedPattern({n: N, extent: EXTENT, density, x, y, sigma: BEAM, patch: 128});
    // The spots are the same spots; only the sampling they land on is coarser,
    // so the disks around them have to be smaller than the gap between them.
    const share = (spots) => spotSum(patch.pattern, patch.n, patch.extent, spots,
      {radius: 1.5, background: 3});
    const fractional = share(SPOTS.fractional);
    const cropped = fractional / (fractional + share(SPOTS.integer));

    assert.ok(Math.abs(cropped - whole) < 0.06,
      `patch ${cropped} vs whole field ${whole} at (${x}, ${y})`);
  }
});

test("the patch follows the beam, and slides in at the edge rather than reading outside", () => {
  const density = field(() => true);
  const sampling = EXTENT / N;
  // A beam closer to the edge than half a patch: the patch stops at the field.
  const corner = illuminatedPattern({
    n: N, extent: EXTENT, density, x: 10, y: 10, sigma: BEAM, patch: 128
  });
  assert.equal(corner.n, 128);
  assert.ok(Math.abs(corner.extent - 128 * sampling) < 1e-9);
  let total = 0;
  for (const value of corner.pattern) total += value;
  assert.ok(total > 0, "a patch at the corner still holds surface");
});

test("the beam window sums to one wherever it sits", () => {
  const window = new Float64Array(N * N);
  for (const [x, y] of [[80, 80], [30, 120], [140, 40]]) {
    beamWindow({n: N, extent: EXTENT, x, y, sigma: BEAM}, window);
    const total = window.reduce((sum, value) => sum + value, 0);
    assert.ok(Math.abs(total - 1) < 1e-12, `total ${total} at (${x}, ${y})`);
  }
});

test("the pattern obeys Parseval against the illuminated patch", () => {
  const density = field(() => true);
  const window = beamWindow({n: N, extent: EXTENT, x: 70, y: 90, sigma: BEAM});
  const pattern = diffractionPattern(density, window, N, complex(N * N));

  let real = 0;
  for (let i = 0; i < density.length; i++) real += (density[i] * window[i]) ** 2;
  let reciprocal = 0;
  for (let i = 0; i < pattern.length; i++) reciprocal += pattern[i];

  const error = Math.abs(reciprocal / (N * N) - real) / real;
  assert.ok(error < 1e-12, `relative Parseval error ${error}`);
});

test("atoms are weighted by depth, and dropped outside the field", () => {
  const one = (depth) => splatAtoms({
    n: 64, extent: 64, atoms: [[32, 32, depth]], sigma: 1, escapeDepth: 3
  }).reduce((sum, value) => sum + value, 0);
  assert.ok(Math.abs(one(3) / one(0) - Math.E ** -1) < 1e-9);

  const flat = splatAtoms({
    n: 64, extent: 64, atoms: [[32, 32, 3]], sigma: 1, escapeDepth: Infinity
  }).reduce((sum, value) => sum + value, 0);
  assert.ok(Math.abs(flat / one(0) - 1) < 1e-9, "an infinite escape depth weights alike");

  const outside = splatAtoms({n: 64, extent: 64, atoms: [[200, 200, 0]], sigma: 1});
  assert.equal(outside.reduce((sum, value) => sum + value, 0), 0);
});

test("a lattice with parallel vectors is refused rather than approximated", () => {
  assert.throws(() => reciprocalVectors([1, 0], [2, 0]), /parallel/);
});
