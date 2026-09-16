import test from "node:test";
import assert from "node:assert/strict";

import {
  transferMatrix, applyTransferMatrix, postprocessOpticalSequence,
  traceRaySequence, traceRays, traceParallel, traceFrom, statesAt,
  bundleRadiusAt, bundleSpreadAt, bundleAt, rayHeightAt, pairedSeparationAt, findCrossovers, findPlanes,
  COLUMNS, COLUMN_MODES, buildColumn, columnDefaults, reverseColumn,
  semElectroOpticalComponents
} from "../kit/raytrace.js";

// ---------------------------------------------------------------------------
// The matrices
// ---------------------------------------------------------------------------

test("a thin lens images according to 1/u + 1/v = 1/f", (t) => {
  // The one identity the whole module rests on. A ray from a point at distance
  // u in front of a lens of focal length f must return to the axis at v.
  const f = 0.25;
  const u = 0.75;
  const v = 1 / (1 / f - 1 / u);
  const components = [{type: "lens", name: "lens", axial_position: u, focal_length: f}];
  const rays = traceRays(components, 0, u + v + 0.5, 0.2, 8);

  for (const trajectory of rays) {
    const atImage = trajectory.find((p) => Math.abs(p.z - (u + v)) < 1e-9)
      ?? interpolate(trajectory, u + v);
    assert.ok(Math.abs(atImage) < 1e-12, `ray misses the image plane by ${atImage}`);
  }
  t.diagnostic(`u = ${u} m, f = ${f} m, image at ${v.toFixed(4)} m`);
});

function interpolate(trajectory, z) {
  for (let i = 0; i < trajectory.length - 1; i++) {
    const a = trajectory[i];
    const b = trajectory[i + 1];
    if (z >= a.z && z <= b.z) {
      const span = b.z - a.z;
      return span === 0 ? b.r : a.r + ((z - a.z) / span) * (b.r - a.r);
    }
  }
  return NaN;
}

test("propagation and lens matrices compose the way optics says", () => {
  // A lens at its own focal distance turns a point source into a parallel beam.
  const f = 0.3;
  const M = transferMatrix.lens(f);
  const [r, theta] = applyTransferMatrix(
    M, applyTransferMatrix(transferMatrix.propagation(f), [0, 0.1])
  );
  assert.ok(Math.abs(r - 0.03) < 1e-15, `height ${r}`);
  assert.ok(Math.abs(theta) < 1e-15, `should be collimated, slope is ${theta}`);
});

// ---------------------------------------------------------------------------
// Sequencing
// ---------------------------------------------------------------------------

test("the sequence interleaves propagation and drops what is above the start", () => {
  const components = [
    {type: "lens", axial_position: 0.5, focal_length: 0.1},
    {type: "lens", axial_position: 1.5, focal_length: 0.1},
    {type: "sample", axial_position: 2.0}
  ];
  const full = postprocessOpticalSequence(components, 0, 3);
  assert.deepEqual(full.map((e) => e.type),
    ["propagation", "lens", "propagation", "lens", "propagation", "sample", "propagation"]);

  // Launching from below the first lens must not apply it. A ray leaving the
  // specimen has not been through the condensers, and before this filter it
  // was — silently, because the component was still in the list.
  const partial = postprocessOpticalSequence(components, 1.0, 3);
  assert.deepEqual(partial.map((e) => e.type),
    ["propagation", "lens", "propagation", "sample", "propagation"]);
  assert.equal(partial[0].dz, 0.5);
});

test("a trajectory launched part-way down starts at that z", () => {
  const sequence = postprocessOpticalSequence([], 1.2, 2.0);
  const trajectory = traceRaySequence([0, 0.1], sequence, 1.2);
  assert.equal(trajectory[0].z, 1.2);
  assert.ok(Math.abs(trajectory.at(-1).z - 2.0) < 1e-12);
  assert.ok(Math.abs(trajectory.at(-1).r - 0.08) < 1e-12);
});

test("an aperture stops the rays outside it and passes the rest", () => {
  const components = [{type: "aperture", axial_position: 1, aperture_radius: 0.05}];
  const rays = traceRays(components, 0, 2, 0.2, 40);
  const blocked = rays.filter((t) => t.at(-1).z < 2 - 1e-9);
  const passed = rays.filter((t) => t.at(-1).z >= 2 - 1e-9);
  assert.ok(blocked.length > 0 && passed.length > 0, "the aperture should do something");
  // Every ray that survived must be inside the radius where the aperture is.
  for (const trajectory of passed) {
    assert.ok(Math.abs(interpolate(trajectory, 1)) <= 0.05 + 1e-12);
  }
  for (const trajectory of blocked) {
    assert.ok(Math.abs(trajectory.at(-1).r) > 0.05 - 1e-12);
  }
});

test("a deflector changes the slope and not the position", () => {
  const components = [{type: "deflector", axial_position: 1, deflection: 0.05}];
  const sequence = postprocessOpticalSequence(components, 0, 2);
  const trajectory = traceRaySequence([0, 0], sequence, 0);
  const atDeflector = trajectory.filter((p) => Math.abs(p.z - 1) < 1e-12);
  for (const point of atDeflector) assert.equal(point.r, 0);
  // and then travels off axis at the deflected slope
  assert.ok(Math.abs(trajectory.at(-1).r - 0.05) < 1e-12, `ended at ${trajectory.at(-1).r}`);
});

test("two opposed deflectors shift the beam without tilting it", () => {
  // The double-deflector identity: a kick and an equal, opposite kick a
  // distance L later leave the ray parallel to where it started, displaced.
  const L = 0.4;
  const kick = 0.05;
  const components = [
    {type: "deflector", axial_position: 0.5, deflection: kick},
    {type: "deflector", axial_position: 0.5 + L, deflection: -kick}
  ];
  const sequence = postprocessOpticalSequence(components, 0, 2);
  const trajectory = traceRaySequence([0, 0], sequence, 0);
  const shift = kick * L;
  assert.ok(Math.abs(trajectory.at(-1).r - shift) < 1e-12);
  // Parallel again: the last two samples are at the same height.
  const last = trajectory.at(-1);
  const beforeLast = trajectory.at(-2);
  assert.ok(Math.abs(last.r - beforeLast.r) < 1e-12, "the beam should be parallel again");
});

// ---------------------------------------------------------------------------
// The two launchers
// ---------------------------------------------------------------------------

test("a parallel bundle stays parallel in free space and focuses at f", () => {
  const f = 0.4;
  const components = [{type: "lens", axial_position: 0.6, focal_length: f}];
  const rays = traceParallel(components, 0, 2, 0.1, 16);
  // Collimated on the way in.
  for (const trajectory of rays) {
    assert.ok(Math.abs(interpolate(trajectory, 0.1) - trajectory[0].r) < 1e-12);
  }
  // A collimated bundle focuses at the back focal plane, whatever its width.
  for (const trajectory of rays) {
    assert.ok(Math.abs(interpolate(trajectory, 0.6 + f)) < 1e-12);
  }
});

test("tilting a parallel bundle moves its focus off axis by f tan(tilt)", (t) => {
  const f = 0.4;
  const tilt = 0.05;
  const components = [{type: "lens", axial_position: 0.6, focal_length: f}];
  const rays = traceParallel(components, 0, 2, 0.1, 16, {tilt});
  const expected = f * Math.tan(tilt);
  for (const trajectory of rays) {
    const r = interpolate(trajectory, 0.6 + f);
    assert.ok(Math.abs(r - expected) < 1e-12, `focus at ${r}, expected ${expected}`);
  }
  t.diagnostic(`tilt ${tilt} rad -> focus ${(expected * 1e3).toFixed(2)} mm off axis`);
});

test("an off-axis point source images off axis by the magnification", () => {
  const f = 0.25;
  const u = 0.75;
  const v = 1 / (1 / f - 1 / u);
  const height = 0.02;
  const components = [{type: "lens", axial_position: u, focal_length: f}];
  const rays = traceRays(components, 0, u + v + 0.2, 0.15, 12, {height});
  const expected = -height * (v / u); // inverted, magnified by v/u
  for (const trajectory of rays) {
    const r = interpolate(trajectory, u + v);
    assert.ok(Math.abs(r - expected) < 1e-12, `landed at ${r}, expected ${expected}`);
  }
});

// ---------------------------------------------------------------------------
// Reading the bundle
// ---------------------------------------------------------------------------

test("findCrossovers locates a known focus", (t) => {
  const f = 0.25;
  const u = 0.75;
  const v = 1 / (1 / f - 1 / u);
  const components = [{type: "lens", axial_position: u, focal_length: f}];
  const rays = traceRays(components, 0, u + v + 0.5, 0.2, 24);
  const found = findCrossovers(rays, {zStart: 0.05, zEnd: u + v + 0.5});
  t.diagnostic(`found ${found.map((c) => c.z.toFixed(4)).join(", ")}, expected ${(u + v).toFixed(4)}`);
  assert.equal(found.length, 1, `expected one crossover, got ${found.length}`);
  assert.ok(Math.abs(found[0].z - (u + v)) < 2e-3, `off by ${found[0].z - (u + v)}`);
});

test("bundleRadiusAt is zero at a focus and grows linearly away from it", () => {
  const rays = traceRays([], 0, 1, 0.1, 20);
  assert.ok(bundleRadiusAt(rays, 0) < 1e-12);
  const a = bundleRadiusAt(rays, 0.25);
  const b = bundleRadiusAt(rays, 0.5);
  assert.ok(Math.abs(b / a - 2) < 1e-9, `free-space growth should be linear, got ${b / a}`);
});

// ---------------------------------------------------------------------------
// Columns as data
// ---------------------------------------------------------------------------

test("every mode names a column that exists", () => {
  for (const [mode, column] of Object.entries(COLUMN_MODES)) {
    assert.ok(COLUMNS[column], `mode ${mode} names a missing column ${column}`);
  }
  // And every S/TEM mode override names a real element.
  for (const column of Object.values(COLUMNS)) {
    const keys = new Set(column.elements.map((e) => e.key).filter(Boolean));
    for (const [mode, overrides] of Object.entries(column.modes ?? {})) {
      for (const key of Object.keys(overrides)) {
        assert.ok(keys.has(key), `mode ${mode} overrides unknown element ${key}`);
      }
    }
  }
});

test("defaults sit inside their slider ranges", () => {
  for (const [name, column] of Object.entries(COLUMNS)) {
    const modes = Object.keys(column.modes ?? {});
    for (const mode of modes.length ? modes : [undefined]) {
      const values = columnDefaults(column, mode);
      for (const element of column.elements) {
        if (element.key === undefined) continue;
        const [lo, hi] = element.range;
        const v = values[element.key];
        assert.ok(v >= lo && v <= hi,
          `${name}/${mode}: ${element.key} = ${v} outside [${lo}, ${hi}]`);
      }
    }
  }
});

test("buildColumn reproduces the hand-written SEM column", () => {
  // The table has to agree with the builder the shipped notebook uses, or
  // swapping the notebook over would quietly change the picture.
  const column = COLUMNS.SEM;
  const values = columnDefaults(column, "SEM");
  const built = buildColumn(column, values);
  const original = semElectroOpticalComponents(
    [values.c1, values.c2, values.objective],
    [values.spray, values.objAperture]
  );
  assert.equal(built.length, original.length);
  for (let i = 0; i < built.length; i++) {
    assert.equal(built[i].type, original[i].type, `component ${i} type`);
    assert.equal(built[i].name, original[i].name, `component ${i} name`);
    assert.equal(built[i].axial_position, original[i].axial_position, `component ${i} z`);
    assert.equal(built[i].focal_length, original[i].focal_length, `component ${i} f`);
    assert.equal(built[i].aperture_radius, original[i].aperture_radius, `component ${i} radius`);
  }
});

test("the objective aperture follows the back focal plane", () => {
  const column = COLUMNS["S/TEM"];
  const values = columnDefaults(column, "TEM imaging");
  const find = (v) => buildColumn(column, v).find((c) => c.name === "objective aperture");
  assert.equal(find(values).axial_position, 2.28 + values.objPost);
  const refocused = {...values, objPost: 0.3};
  assert.equal(find(refocused).axial_position, 2.28 + 0.3);
});

test("the three S/TEM modes differ only in their defaults", () => {
  const column = COLUMNS["S/TEM"];
  const shapes = Object.keys(column.modes).map((mode) =>
    buildColumn(column, columnDefaults(column, mode))
      .map((c) => `${c.type}:${c.name}`)
      .join("|")
  );
  assert.equal(new Set(shapes).size, 1, "the modes should share their components");

  // And they do actually differ, or the presets are pointless.
  const values = Object.keys(column.modes).map((m) => JSON.stringify(columnDefaults(column, m)));
  assert.equal(new Set(values).size, 3, "the three modes should have distinct defaults");
});

test("TEM imaging and TEM diffraction put different planes on the detector", (t) => {
  // The claim the widget makes. Both modes illuminate the specimen the same
  // way; the intermediate lens decides whether the detector sees an image of
  // the specimen or of the back focal plane. So the illumination bundle — which
  // is focused at the back focal plane — must come to a focus at the detector
  // in diffraction mode and not in imaging mode.
  const column = COLUMNS["S/TEM"];
  const detector = column.zEnd;

  const radiusAtDetector = (mode) => {
    const values = columnDefaults(column, mode);
    const components = buildColumn(column, values);
    const rays = traceRays(components, 0, detector, column.semiangle, 200);
    const widest = Math.max(
      ...Array.from({length: 200}, (_, i) =>
        bundleRadiusAt(rays, (detector * i) / 199))
    );
    return bundleRadiusAt(rays, detector) / widest;
  };

  const diffraction = radiusAtDetector("TEM diffraction");
  const imaging = radiusAtDetector("TEM imaging");
  t.diagnostic(`bundle at the detector: diffraction ${(diffraction * 100).toFixed(1)}%, imaging ${(imaging * 100).toFixed(1)}% of its widest`);
  assert.ok(diffraction < imaging / 3,
    `diffraction mode should focus the source at the detector (${diffraction} vs ${imaging})`);
});

test("an off-axis object's image is found where the rays converge, not on the axis", (t) => {
  // The bundle from an off-axis specimen point converges to an off-axis image.
  // Measuring the spread about the optic axis finds no minimum there at all,
  // because the whole bundle is far from the axis; measuring it about the
  // bundle's own centre finds it exactly.
  const f = 0.15;
  const u = 0.225;
  const v = 1 / (1 / f - 1 / u);
  const height = 0.0124;
  const lensZ = 2.28;
  const components = [{type: "lens", axial_position: lensZ, focal_length: f}];
  const rays = traceRays(components, lensZ - u, lensZ + v + 0.4, 0.035, 32, {height});

  const imageZ = lensZ + v;
  const expected = -height * (v / u);
  t.diagnostic(`image at z = ${imageZ.toFixed(3)} m, ${(expected * 1e3).toFixed(2)} mm off axis`);

  assert.ok(bundleSpreadAt(rays, imageZ) < 1e-12, "the rays should converge there");
  assert.ok(Math.abs(bundleRadiusAt(rays, imageZ) - Math.abs(expected)) < 1e-9,
    "and be far from the axis while they do");

  const found = findCrossovers(rays, {zStart: lensZ - u + 0.02, zEnd: lensZ + v + 0.4});
  assert.equal(found.length, 1, `expected one image, got ${found.length}`);
  assert.ok(Math.abs(found[0].z - imageZ) < 2e-3, `off by ${found[0].z - imageZ}`);
});

test("spread and radius agree for a fan that is symmetric about the axis", () => {
  // Which is why using spread everywhere is safe: it generalises the on-axis
  // case rather than replacing it.
  const rays = traceRays([{type: "lens", axial_position: 0.5, focal_length: 0.2}], 0, 1.5, 0.1, 21);
  for (const z of [0.2, 0.5, 0.8, 1.2]) {
    assert.ok(Math.abs(bundleSpreadAt(rays, z) - bundleRadiusAt(rays, z)) < 1e-12, `at z = ${z}`);
  }
});

test("pairing survives an aperture clipping one beam", (t) => {
  // The bug this guards: an objective aperture takes the rays furthest off
  // axis, and a tilted beam is not clipped symmetrically, so the surviving
  // rays' *mean* shifts sideways. Comparing the two beams by their means then
  // finds the image plane in the wrong place — in the widget it moved onto the
  // diffraction planes as soon as the scattering angle clipped the beam.
  const column = COLUMNS["S/TEM"];
  const sample = column.elements.find((e) => e.kind === "sample").z;
  const values = columnDefaults(column, "TEM imaging");
  const components = buildColumn(column, values);

  const illumination = traceRays(components, 0, column.zEnd, column.semiangle, 192);
  const states = statesAt(illumination, sample);

  const imagesAt = (mrad) => {
    const tilt = Math.tan(mrad * 1e-3);
    const direct = traceFrom(components, sample, column.zEnd, states);
    const diffracted = traceFrom(components, sample, column.zEnd,
      states.map(([r, theta]) => [r, theta + tilt]));
    const surviving = diffracted.filter((r) => r.at(-1).z >= column.zEnd - 1e-9).length;
    const planes = findPlanes((z) => {
      const {separation, spread, count} = pairedSeparationAt(direct, diffracted, z);
      return count === 0 ? Infinity : separation / Math.max(spread, 1e-9);
    }, {zStart: sample + 0.08, zEnd: column.zEnd, tolerance: 0.05});
    return {surviving, total: diffracted.length, planes: planes.map((p) => p.z)};
  };

  const clear = imagesAt(30);
  const clipped = imagesAt(66);
  t.diagnostic(`30 mrad: ${clear.surviving}/${clear.total} rays, images at ${clear.planes.map((z) => z.toFixed(3))}`);
  t.diagnostic(`66 mrad: ${clipped.surviving}/${clipped.total} rays, images at ${clipped.planes.map((z) => z.toFixed(3))}`);

  assert.ok(clipped.surviving > 0 && clipped.surviving < clipped.total,
    "66 mrad should clip the beam partially, or this tests nothing");
  assert.equal(clipped.planes.length, clear.planes.length, "same number of image planes");
  for (let i = 0; i < clear.planes.length; i++) {
    assert.ok(Math.abs(clipped.planes[i] - clear.planes[i]) < 5e-3,
      `image plane ${i} moved from ${clear.planes[i]} to ${clipped.planes[i]} under clipping`);
  }
  // And the last of them is the detector, which is what imaging mode is for.
  assert.ok(Math.abs(clear.planes.at(-1) - column.zEnd) < 5e-3);
});

test("pairedSeparationAt counts only rays both bundles still have", () => {
  const components = [{type: "aperture", axial_position: 1, aperture_radius: 0.05}];
  const a = traceRays(components, 0, 2, 0.1, 20);
  const b = traceRays(components, 0, 2, 0.1, 20, {tilt: 0.04});
  const {count} = pairedSeparationAt(a, b, 2);
  const aliveA = a.filter((t) => t.at(-1).z >= 2 - 1e-9).length;
  const aliveB = b.filter((t) => t.at(-1).z >= 2 - 1e-9).length;
  assert.ok(count <= Math.min(aliveA, aliveB), `counted ${count} of ${aliveA}/${aliveB}`);
  assert.ok(count > 0, "some rays should survive both");

  // Identical bundles have no separation, however many rays an aperture took.
  const {separation} = pairedSeparationAt(a, a, 2);
  assert.equal(separation, 0);

  // And a plane no ray reaches reports no pairing rather than a coincidence.
  assert.equal(pairedSeparationAt(a, b, 5).count, 0);
  assert.equal(pairedSeparationAt(a, b, 5).separation, Infinity);
});

// ---------------------------------------------------------------------------
// Reciprocity
// ---------------------------------------------------------------------------

test("a ray reversed through a reversed column retraces its own path", (t) => {
  // The theorem the reciprocity widget is about. Send a ray down a column; take
  // where and how steeply it comes out; send it back in from the far end with
  // its slope negated, through the column with its ends swapped. It must
  // retrace the same path, mirrored.
  const zStart = 0;
  const zEnd = 2;
  const components = [
    {type: "lens", name: "condenser", axial_position: 0.4, focal_length: 0.3},
    {type: "deflector", name: "tilt", axial_position: 0.8, deflection: 0.06},
    {type: "sample", name: "specimen", axial_position: 1.0},
    {type: "lens", name: "objective", axial_position: 1.4, focal_length: 0.25}
  ];

  const forward = traceRaySequence(
    [0.01, 0.05], postprocessOpticalSequence(components, zStart, zEnd), zStart
  );
  const exit = forward.at(-1);

  const mirrored = reverseColumn(components, {zStart, zEnd});
  const backward = traceRaySequence(
    [exit.r, -exit.theta], postprocessOpticalSequence(mirrored, zStart, zEnd), zStart
  );

  // The backward ray at z must sit where the forward ray was at zStart+zEnd-z.
  let worst = 0;
  let scale = 0;
  for (let i = 0; i <= 200; i++) {
    const z = zStart + ((zEnd - zStart) * i) / 200;
    const there = rayHeightAt(backward, z);
    const here = rayHeightAt(forward, zStart + zEnd - z);
    if (there === null || here === null) continue;
    worst = Math.max(worst, Math.abs(there - here));
    scale = Math.max(scale, Math.abs(here));
  }
  t.diagnostic(`max mismatch ${worst.toExponential(2)} m over a path reaching ${scale.toFixed(4)} m`);
  assert.ok(worst < 1e-12, `the reversed ray is off by ${worst}`);

  // And it comes back out where it went in, which is the statement in the form
  // people quote: swap source and detector and nothing changes.
  assert.ok(Math.abs(backward.at(-1).r - forward[0].r) < 1e-12);
  assert.ok(Math.abs(backward.at(-1).theta + forward[0].theta) < 1e-12);
});

test("reversing a column twice gives the column back", () => {
  const span = {zStart: 0, zEnd: 2};
  const components = [
    {type: "lens", axial_position: 0.4, focal_length: 0.3},
    {type: "deflector", axial_position: 0.8, deflection: 0.06},
    {type: "aperture", axial_position: 1.4, aperture_radius: 0.05}
  ];
  const twice = reverseColumn(reverseColumn(components, span), span);
  // Positions compared with a tolerance rather than exactly: 2 - 0.4 is 1.6 and
  // 2 - 1.6 is 0.3999999999999999, so an exact round trip is not something
  // subtraction offers.
  assert.equal(twice.length, components.length);
  for (let i = 0; i < components.length; i++) {
    assert.ok(Math.abs(twice[i].axial_position - components[i].axial_position) < 1e-12);
    assert.equal(twice[i].type, components[i].type);
    assert.equal(twice[i].focal_length, components[i].focal_length);
    assert.equal(twice[i].deflection, components[i].deflection);
    assert.equal(twice[i].aperture_radius, components[i].aperture_radius);
  }
});

test("a TEM and a bright-field STEM are the same column reversed", (t) => {
  // The widget's claim, on the columns it actually draws. Parallel illumination
  // onto a specimen with the objective and aperture below it, reversed, becomes
  // an aperture and objective above the specimen with a point detector below.
  const span = {zStart: 0, zEnd: 1.5};
  const focal = 0.205;
  const tem = [
    {type: "sample", name: "specimen", axial_position: 0.75},
    {type: "lens", name: "objective", axial_position: 1.0, focal_length: focal},
    // At the back focal plane, which is where an objective aperture belongs.
    {type: "aperture", name: "objective aperture", axial_position: 1.0 + focal,
     aperture_radius: 0.0875}
  ];
  const stem = reverseColumn(tem, span);

  const at = (list, name) => list.find((c) => c.name === name).axial_position;
  t.diagnostic(`specimen ${at(tem, "specimen")} -> ${at(stem, "specimen")}, ` +
    `objective ${at(tem, "objective")} -> ${at(stem, "objective")}`);
  assert.equal(at(stem, "specimen"), 0.75, "the specimen sits midway, so it does not move");
  assert.equal(at(stem, "objective"), 0.5);
  assert.ok(at(stem, "objective aperture") < at(stem, "objective"),
    "the aperture should now come before the objective");

  // A parallel beam into the TEM focuses at the objective's back focal plane.
  // Reversed, a beam through the aperture leaves the specimen parallel.
  const parallel = traceParallel(tem, span.zStart, span.zEnd, 0.06, 24);
  const bfp = at(tem, "objective") + focal;
  assert.ok(bundleSpreadAt(parallel, bfp) < 1e-12,
    "parallel illumination focuses at the objective aperture");

  // And the content of reciprocity on these two columns: every ray of that
  // illumination, sent back in from the detector end with its slope negated,
  // retraces its own path through the STEM column.
  let worst = 0;
  for (const ray of parallel) {
    const exit = ray.at(-1);
    const back = traceRaySequence(
      [exit.r, -exit.theta], postprocessOpticalSequence(stem, span.zStart, span.zEnd), span.zStart
    );
    for (let i = 0; i <= 60; i++) {
      const z = span.zStart + ((span.zEnd - span.zStart) * i) / 60;
      const there = rayHeightAt(back, z);
      const here = rayHeightAt(ray, span.zStart + span.zEnd - z);
      if (there === null || here === null) continue;
      worst = Math.max(worst, Math.abs(there - here));
    }
  }
  t.diagnostic(`worst retrace mismatch across the bundle: ${worst.toExponential(2)} m`);
  assert.ok(worst < 1e-12, `off by ${worst}`);
});
