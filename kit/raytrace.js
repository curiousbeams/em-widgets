// Paraxial ray tracing by transfer matrices.
//
// Lifted verbatim from the pinned helper cells of
// `observable-notebooks/sem-ray-diagram.html` so that notebook (and the AP3402
// w4.* ray-diagram pages, which currently ship static Mathematica output) can
// share one implementation.
//
// A ray is `[r, theta]`: off-axis distance and slope. Optical components are
// plain objects with an `axial_position` and a `type` of "lens", "aperture" or
// "sample"; everything is in consistent length units (the SEM notebook uses metres).

/** 2x2 paraxial transfer matrices. */
export const transferMatrix = {
  /** Free-space propagation over a distance d. */
  propagation: (d) => [
    [1, d],
    [0, 1]
  ],
  /** A thin lens of focal length f. */
  lens: (f) => [
    [1, 0],
    [-1 / f, 1]
  ]
};

/** Apply a 2x2 matrix to a ray. */
export const applyTransferMatrix = (M, ray) => [
  M[0][0] * ray[0] + M[0][1] * ray[1],
  M[1][0] * ray[0] + M[1][1] * ray[1]
];

/**
 * Sort components by axial position and interleave free-space propagation steps,
 * so the result can be walked in a single pass.
 */
export function postprocessOpticalSequence(components, zStart, zEnd) {
  const sorted = [...components].sort((a, b) => a.axial_position - b.axial_position);
  const sequence = [];
  let zPrev = zStart;
  for (const component of sorted) {
    if (component.axial_position > zPrev) {
      sequence.push({type: "propagation", dz: component.axial_position - zPrev});
    }
    sequence.push(component);
    zPrev = component.axial_position;
  }
  if (zPrev < zEnd) sequence.push({type: "propagation", dz: zEnd - zPrev});
  return sequence;
}

/**
 * Trace a single ray through a prepared sequence.
 *
 * Returns the trajectory as `{z, r}` samples. A ray blocked by an aperture stops
 * there, which is what makes the aperture sliders visually meaningful.
 */
export function traceRaySequence(ray0, sequence) {
  let ray = [...ray0];
  let z = 0;
  const trajectory = [{z, r: ray[0]}];

  for (const element of sequence) {
    switch (element.type) {
      case "propagation":
        z += element.dz;
        ray = applyTransferMatrix(transferMatrix.propagation(element.dz), ray);
        trajectory.push({z, r: ray[0]});
        break;
      case "lens":
        ray = applyTransferMatrix(transferMatrix.lens(element.focal_length), ray);
        trajectory.push({z, r: ray[0]});
        break;
      case "aperture":
        if (Math.abs(ray[0]) > element.aperture_radius) return trajectory; // blocked
        break;
    }
  }
  return trajectory;
}

/**
 * Trace a fan of rays leaving the origin over +/- `semiangle`.
 *
 * @param {object[]} components
 * @param {number} zStart
 * @param {number} zEnd
 * @param {number} semiangle half-angle of the fan [rad]
 * @param {number} [numRays=128]
 * @returns {Array<Array<{z: number, r: number}>>}
 */
export function traceRays(components, zStart, zEnd, semiangle, numRays = 128) {
  const sequence = postprocessOpticalSequence(components, zStart, zEnd);
  const angles = Array.from(
    {length: numRays + 1},
    (_, i) => -semiangle + (2 * semiangle * i) / numRays
  );
  return angles.map((alpha) => traceRaySequence([0, Math.tan(alpha)], sequence));
}

/**
 * The two-condenser SEM column used by the sem-ray-diagram notebook.
 *
 * @param {[number, number, number]} focalLengths C1, C2 and objective focal lengths
 * @param {[number, number]} apertureRadii spray and objective aperture radii
 */
export const semElectroOpticalComponents = ([f1, f2, f3], [ap1, ap2]) => [
  {type: "lens", name: "1st condenser lens", axial_position: 0.2, focal_length: f1},
  {type: "aperture", name: "spray aperture", axial_position: 0.5, aperture_radius: ap1},
  {type: "lens", name: "2nd condenser lens", axial_position: 0.7, focal_length: f2},
  {type: "aperture", name: "objective aperture", axial_position: 1.2, aperture_radius: ap2},
  {type: "lens", name: "objective lens", axial_position: 1.3, focal_length: f3},
  {type: "sample", name: "sample", axial_position: 1.5}
];
