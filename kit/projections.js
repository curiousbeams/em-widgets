// Projection-set algorithms: the family of iterations that solve a feasibility
// problem by repeatedly projecting onto each constraint.
//
// Phase retrieval is one of these. The measured intensities are one set, the
// support and reality constraints another, and the reconstruction is a point in
// the intersection. Everything interesting about the algorithms — why
// alternating projections stalls, why the difference map escapes — is already
// visible in two dimensions with sets you can draw, which is what this module
// is for.
//
// Ported from `@gvarnavi/projection-set-algorithms`. The original hard-codes
// each set's nearest-point map as a branch on inequalities in x and y; here a
// set knows its own geometry and the projection follows from it.

/**
 * The generalised projection of Combettes, written as one iteration:
 *
 *     x <- (1 - a - b) x + a P1[x] + b P2[c P1[x] + (1 - c) x]
 *
 * Each named algorithm is a choice of `(a, b, c)` — see {@link namedParameters}.
 *
 * @param {object} first a set with a `project` method
 * @param {object} second
 * @param {[number, number, number]} parameters `[a, b, c]`
 * @returns {(point: number[]) => {first: number[], second: number[], next: number[]}}
 */
export function generalisedProjection(first, second, [a, b, c]) {
  const keep = 1 - a - b;
  const reflect = 1 - c;
  return ([x, y]) => {
    const [ax, ay] = first.project([x, y]);
    const [bx, by] = second.project([c * ax + reflect * x, c * ay + reflect * y]);
    return {
      first: [ax, ay],
      second: [bx, by],
      // The point the second projection was actually applied to, which is not
      // the iterate unless c is 1. A drawing that joins the iterate straight to
      // the second projection is drawing a step the algorithm never took.
      reflected: [c * ax + reflect * x, c * ay + reflect * y],
      next: [keep * x + a * ax + b * bx, keep * y + a * ay + b * by]
    };
  };
}

/**
 * The same iteration for more than two sets, lifted into a product space.
 *
 * With N sets the iterate becomes N copies of the point. Projecting each copy
 * onto its own set plays the role of the first projection, and averaging the
 * copies back together plays the role of the second — the copies agree exactly
 * when the point is in every set at once.
 */
export function productProjection(sets, [a, b, c]) {
  const keep = 1 - a - b;
  const reflect = 1 - c;
  return (copies) => {
    const projected = copies.map(([x, y], i) => sets[i].project([x, y]));
    const reflected = projected.map(([px, py], i) => [
      c * px + reflect * copies[i][0],
      c * py + reflect * copies[i][1]
    ]);
    const [mx, my] = centroid(reflected);
    return {
      projected,
      averaged: [mx, my],
      next: copies.map(([x, y], i) => [
        keep * x + a * projected[i][0] + b * mx,
        keep * y + a * projected[i][1] + b * my
      ])
    };
  };
}

/** The mean of a list of points. */
export function centroid(points) {
  let x = 0;
  let y = 0;
  for (const [px, py] of points) {
    x += px;
    y += py;
  }
  return [x / points.length, y / points.length];
}

/**
 * The named algorithms, as `(a, b, c)`.
 *
 * @param {number} [gamma=0.875] relaxation, used by RRR and RAAR
 */
export function namedParameters(gamma = 0.875) {
  return {
    AP: [0, 1, 1],
    DM: [-1, 1, 2],
    RRR: [-gamma, gamma, 2],
    RAAR: [1 - 2 * gamma, gamma, 2]
  };
}

// ---------------------------------------------------------------------------
// Sets
// ---------------------------------------------------------------------------

/**
 * A line through the origin, optionally with stretches removed.
 *
 * `excluded` is a list of `[from, to]` intervals in the line's own parameter,
 * where the point at parameter `t` is `t * direction`. Removing one makes the
 * set non-convex, and non-convexity is the entire reason these algorithms
 * differ from each other: a point whose nearest neighbour falls in a gap gets
 * pushed to whichever end of the gap is closer, and alternating projections can
 * then get stuck bouncing between two such ends forever.
 *
 * @param {object} options
 * @param {[number, number]} options.direction
 * @param {Array<[number, number]>} [options.excluded=[]]
 * @param {number} [options.extent=1] how far along the parameter to draw
 */
export function line({direction: [dx, dy], excluded = [], extent = 1, name = "line"}) {
  const norm2 = dx * dx + dy * dy;
  const at = (t) => [t * dx, t * dy];

  const project = ([x, y]) => {
    let t = (x * dx + y * dy) / norm2;
    for (const [from, to] of excluded) {
      if (t > from && t < to) t = t - from < to - t ? from : to;
    }
    return at(t);
  };

  // Drawing: the parameter range minus the gaps, as one polyline per stretch.
  const path = (() => {
    const cuts = [-extent, ...excluded.flat(), extent].sort((p, q) => p - q);
    const out = [];
    for (let i = 0; i < cuts.length - 1; i += 2) {
      out.push([at(cuts[i]), at(cuts[i + 1])]);
    }
    return out;
  })();

  return {name, kind: "line", project, path, closed: false};
}

/**
 * A closed polar curve, sampled, with the nearest sample as its projection.
 *
 *     r(theta) = 1 + amplitude * sin(petals * theta)
 *
 * Brute force over the samples rather than anything cleverer: a thousand
 * distances is nothing, and the alternative is solving a trigonometric
 * polynomial for every step of every iteration.
 *
 * Sampling makes the projection piecewise constant, so the iteration cannot
 * converge below the sample spacing. That is a property of the set, not an
 * error in the algorithm, and the widget's residual plot bottoms out there.
 */
export function polarCurve({
  petals = 4, amplitude = 0.25, samples = 1024,
  centre = [0, 0], scale = 1, name = "curve"
} = {}) {
  const points = new Array(samples);
  for (let i = 0; i < samples; i++) {
    const theta = (2 * Math.PI * i) / samples;
    const r = 1 + amplitude * Math.sin(petals * theta);
    points[i] = [
      centre[0] + scale * r * Math.cos(theta),
      centre[1] + scale * r * Math.sin(theta)
    ];
  }

  const project = ([x, y]) => {
    let best = 0;
    let bestDistance = Infinity;
    for (let i = 0; i < samples; i++) {
      const dx = points[i][0] - x;
      const dy = points[i][1] - y;
      const d = dx * dx + dy * dy;
      if (d < bestDistance) {
        bestDistance = d;
        best = i;
      }
    }
    return points[best];
  };

  return {name, kind: "curve", project, path: [[...points, points[0]]], closed: true, points};
}

// ---------------------------------------------------------------------------
// Running an iteration
// ---------------------------------------------------------------------------

/**
 * How far a point is from satisfying every constraint at once.
 *
 * The rms distance from the point to each set. Zero exactly when the point lies
 * in the intersection, and it is what separates the algorithms: in a non-convex
 * trap alternating projections settles at a residual it cannot reduce, while
 * the difference map keeps moving and eventually finds a way through. Neither
 * fact is visible in a plot of the trajectory alone, which is why the widget
 * draws this alongside.
 */
export function residual(point, sets) {
  let sum = 0;
  for (const set of sets) {
    const [px, py] = set.project(point);
    sum += (px - point[0]) ** 2 + (py - point[1]) ** 2;
  }
  return Math.sqrt(sum / sets.length);
}

/**
 * An iterator over whichever form the algorithm takes for this many sets.
 *
 * Two sets use {@link generalisedProjection} directly; more are lifted into the
 * product space. Both are hidden behind the same `step`, so a caller holds one
 * opaque state and reads `point` out of it.
 *
 * @param {object[]} sets
 * @param {[number, number, number]} parameters
 * @returns {{start: (p: number[]) => object, step: (s: object) => object}}
 */
export function iteration(sets, parameters) {
  if (sets.length === 2) {
    const advance = generalisedProjection(sets[0], sets[1], parameters);
    return {
      start: (point) => ({
        point,
        projections: sets.map((s) => s.project(point)),
        rays: sets.map((s) => ({from: point, to: s.project(point)}))
      }),
      step: (state) => {
        const {next, first, second, reflected} = advance(state.point);
        return {
          point: next,
          projections: [first, second],
          // One ray per set, each from whatever was projected onto it. With two
          // sets that is the iterate and then the reflected point, which for
          // alternating projections is the textbook zig-zag between them.
          rays: [{from: state.point, to: first}, {from: reflected, to: second}]
        };
      }
    };
  }

  const advance = productProjection(sets, parameters);
  return {
    start: (point) => ({
      copies: sets.map(() => [...point]),
      point,
      projections: sets.map((s) => s.project(point)),
      rays: sets.map((s) => ({from: point, to: s.project(point)}))
    }),
    step: (state) => {
      const {next, projected} = advance(state.copies);
      return {
        copies: next,
        point: centroid(next),
        projections: projected,
        // Here each set gets its own copy of the iterate, so the rays start
        // from different places. Drawing them all from the average would hide
        // the fact that the copies have drifted apart, which is the whole
        // mechanism by which more than two sets are handled.
        rays: projected.map((to, i) => ({from: state.copies[i], to}))
      };
    }
  };
}
