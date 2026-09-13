import test from "node:test";
import assert from "node:assert/strict";

import {angularSpatialFrequencies} from "../kit/grid.js";
import {electronWavelength} from "../kit/units.js";
import {
  annularMask, segmentedDetector, integrateSignals, integrateOne,
  centreOfMassFromSegments, centreOfMass, segmentIndex, integrateByIndex
} from "../kit/detectors.js";

const gpts = [96, 96];
const sampling = [0.2, 0.2];
const wavelength = electronWavelength(80e3);
const angular = angularSpatialFrequencies(gpts, sampling, wavelength);
const alphaMax = (wavelength / (2 * sampling[0])) * 1e3; // mrad

test("the angular grid spans a sensible range", (t) => {
  t.diagnostic(`alpha max = ${alphaMax.toFixed(1)} mrad at 80 kV, 0.2 Å sampling`);
  assert.ok(alphaMax > 50 && alphaMax < 200);
});

test("annularMask selects exactly the requested band", () => {
  const mask = annularMask(angular, 20, 40);
  for (let i = 0; i < mask.length; i++) {
    const mrad = angular.alpha[i] * 1e3;
    const expected = mrad >= 20 && mrad < 40 ? 1 : 0;
    assert.equal(mask[i], expected);
  }
});

test("the 16 segments tile the detector exactly, with no overlap or gap", () => {
  const det = segmentedDetector(angular, {inner: 10, outer: 60, rings: 2});
  assert.equal(det.masks.length, 16);

  const whole = annularMask(angular, 10, 60);
  const covered = new Float64Array(whole.length);
  for (const mask of det.masks) {
    for (let i = 0; i < mask.length; i++) covered[i] += mask[i];
  }
  for (let i = 0; i < whole.length; i++) {
    assert.ok(covered[i] <= 1, `pixel ${i} is in ${covered[i]} segments — segments overlap`);
    assert.equal(covered[i], whole[i], `pixel ${i}: covered ${covered[i]}, expected ${whole[i]}`);
  }
});

test("every segment gets a share of the detector", () => {
  const det = segmentedDetector(angular, {inner: 10, outer: 60, rings: 2});
  for (let s = 0; s < det.area.length; s++) {
    assert.ok(det.area[s] > 0, `segment ${s} is empty`);
  }
});

test("sector areas are equal once pixelation is not the limit", (t) => {
  // On the widget's own grid the inner disc is ~11 px in radius, so a sector
  // holds only ~50 pixels and pixelation dominates. Check the geometry on a grid
  // fine enough for a real asymmetry to show — this is what caught the sector
  // boundaries sitting on the grid diagonals.
  const fine = angularSpatialFrequencies([512, 512], sampling, wavelength);
  const det = segmentedDetector(fine, {inner: 10, outer: 60, rings: 2});
  for (const zone of [0, 1]) {
    const areas = det.area.slice(zone * 8, zone * 8 + 8);
    const mean = areas.reduce((a, b) => a + b, 0) / 8;
    const spread = Math.max(...areas.map((a) => Math.abs(a / mean - 1)));
    t.diagnostic(`zone ${zone}: ${mean.toFixed(0)} px/sector, spread ${(spread * 100).toFixed(1)}%`);
    assert.ok(spread < 0.02, `sector areas uneven in zone ${zone}: ${areas}`);
  }
});

test("segment centroids point outward, one per sector direction", () => {
  const det = segmentedDetector(angular, {inner: 10, outer: 60, rings: 2});
  const seen = new Set();
  for (let s = 0; s < 8; s++) {
    const c = det.centroids[s];
    const angle = Math.atan2(c.ky, c.kx);
    const sector = Math.round((((angle / (2 * Math.PI)) % 1) + 1) % 1 * 8) % 8;
    seen.add(sector);
    assert.ok(Math.hypot(c.kx, c.ky) > 0, `segment ${s} centroid is at the origin`);
  }
  assert.equal(seen.size, 8, "the eight disc sectors should point in eight directions");
});

/** A uniform disc of radius `r` mrad, shifted by (dx, dy) mrad. */
function shiftedDisc(radiusMrad, dxMrad, dyMrad) {
  const I = new Float64Array(angular.alpha.length);
  for (let i = 0; i < I.length; i++) {
    const kx = angular.kx[i] * wavelength * 1e3;   // mrad
    const ky = angular.ky[i] * wavelength * 1e3;
    I[i] = Math.hypot(kx - dxMrad, ky - dyMrad) < radiusMrad ? 1 : 0;
  }
  return I;
}

test("a pixelated CoM recovers the shift of a uniform disc", (t) => {
  for (const [dx, dy] of [[0, 0], [4, 0], [-3, 5]]) {
    const I = shiftedDisc(20, dx, dy);
    const com = centreOfMass(I, angular);
    const mx = com.kx * wavelength * 1e3;
    const my = com.ky * wavelength * 1e3;
    t.diagnostic(`shift (${dx}, ${dy}) -> CoM (${mx.toFixed(2)}, ${my.toFixed(2)}) mrad`);
    assert.ok(Math.hypot(mx - dx, my - dy) < 0.6, `got (${mx}, ${my}), expected (${dx}, ${dy})`);
  }
});

test("a segmented CoM follows the shift, more coarsely than a pixelated one", (t) => {
  const det = segmentedDetector(angular, {inner: 10, outer: 60, rings: 2});
  let worstSegmented = 0, worstPixelated = 0;
  for (const [dx, dy] of [[4, 0], [-3, 5], [0, -6]]) {
    const I = shiftedDisc(20, dx, dy);

    const seg = centreOfMassFromSegments(integrateSignals(I, det.masks), det.centroids);
    const pix = centreOfMass(I, angular, annularMask(angular, 10, 60));

    const segErr = Math.hypot(seg.kx * wavelength * 1e3 - dx, seg.ky * wavelength * 1e3 - dy);
    const pixErr = Math.hypot(pix.kx * wavelength * 1e3 - dx, pix.ky * wavelength * 1e3 - dy);
    worstSegmented = Math.max(worstSegmented, segErr);
    worstPixelated = Math.max(worstPixelated, pixErr);

    // The segmented estimate must at least point the right way.
    const dot = seg.kx * dx + seg.ky * dy;
    assert.ok(dot > 0, `segmented CoM points the wrong way for shift (${dx}, ${dy})`);
  }
  t.diagnostic(`worst error: segmented ${worstSegmented.toFixed(2)} mrad, pixelated ${worstPixelated.toFixed(2)} mrad`);
  assert.ok(worstPixelated < worstSegmented, "a pixelated detector should be the more accurate of the two");
});

test("integrateOne agrees with integrateSignals", () => {
  const det = segmentedDetector(angular, {inner: 10, outer: 60, rings: 2});
  const I = shiftedDisc(20, 2, -1);
  const all = integrateSignals(I, det.masks);
  for (let s = 0; s < det.masks.length; s++) {
    assert.ok(Math.abs(integrateOne(I, det.masks[s]) - all[s]) < 1e-12);
  }
});

test("segmentIndex and integrateByIndex agree with the mask-by-mask version", () => {
  const det = segmentedDetector(angular, {inner: 10, outer: 60, rings: 2});
  const I = shiftedDisc(20, 2, -1);
  const index = segmentIndex(det.masks, angular.alpha.length);

  // Every pixel is claimed by at most one segment, and by the right one.
  for (let i = 0; i < index.length; i++) {
    let owner = -1;
    for (let m = 0; m < det.masks.length; m++) if (det.masks[m][i]) owner = m;
    assert.equal(index[i], owner, `pixel ${i}`);
  }

  const slow = integrateSignals(I, det.masks);
  const fast = integrateByIndex(I, index, det.masks.length);
  for (let s = 0; s < slow.length; s++) {
    assert.ok(Math.abs(slow[s] - fast[s]) < 1e-9, `segment ${s}: ${slow[s]} vs ${fast[s]}`);
  }
});

test("integrateByIndex reuses its output buffer without accumulating", () => {
  const det = segmentedDetector(angular, {inner: 10, outer: 60, rings: 2});
  const index = segmentIndex(det.masks, angular.alpha.length);
  const out = new Float64Array(det.masks.length);
  const I = shiftedDisc(20, 0, 0);
  const first = Float64Array.from(integrateByIndex(I, index, det.masks.length, out));
  const second = integrateByIndex(I, index, det.masks.length, out);
  for (let s = 0; s < first.length; s++) assert.equal(first[s], second[s]);
});
