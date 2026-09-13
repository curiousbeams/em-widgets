// Animation helpers.
//
// This module is deliberately small, because the Observable runtime already does
// the hard part: its scheduler *is* `requestAnimationFrame`, and it pulls
// generator cells with `.next()` on every tick. So a generator cell advances once
// per frame and every dependent cell re-runs — no rAF plumbing, no manual
// invalidation.
//
//     const t = frames();                 // a cell: advances ~60x a second
//     const angle = t * 0.01;             // any cell that reads it re-runs
//
// ONE TRAP, worth internalising: a generator cell that *depends* on a reactive
// input is destroyed and recreated when that input changes. Any state it was
// accumulating — an eased value, a phase — resets. So anything that must ease
// across an input change has to own that input itself, which is why
// `hoverEase` attaches its own listeners rather than taking a boolean.

/**
 * A frame counter. Use as a generator cell to drive an animation.
 *
 * @param {object} [options]
 * @param {number} [options.fps] cap the rate; omit to run at display refresh.
 *   The references throttle to ~30 fps on desktop and ~19 on mobile, since the
 *   scenes cost more to draw than they gain from 60.
 * @param {boolean} [options.playing=true] when false, yields 0 forever — so a
 *   play/pause toggle is just `frames({playing: isPlaying})`.
 * @param {number} [options.start=0]
 * @yields {number} frame index
 */
export function* frames({fps, playing = true, start = 0} = {}) {
  if (!playing) {
    while (true) yield start;
  }
  const interval = fps ? 1000 / fps : 0;
  let i = start;
  let last = -Infinity;
  while (true) {
    const now = performance.now();
    // Throttle inside the loop rather than with a timer: the runtime is already
    // calling us once per animation frame, so we just skip the ones we don't want.
    if (now - last >= interval) {
      last = now;
      yield i++;
    } else {
      yield i;
    }
  }
}

/**
 * Elapsed seconds since the cell started. Convenient for `Math.cos(t * 0.5)`
 * style idle motion.
 */
export function* seconds({fps, playing = true} = {}) {
  const t0 = performance.now();
  for (const _ of frames({fps, playing})) {
    yield (performance.now() - t0) / 1000;
  }
}

/** Exponential smoothing, the references' easing. `k` 0.07 ≈ 200 ms at 60 fps. */
export function ease(value, target, k = 0.07) {
  return value + (target - value) * k;
}

/**
 * A 0→1 value that eases in while the pointer is over `el` and out when it
 * leaves. Use it to crossfade between an idle animation and reader control:
 *
 *     const hov = hoverEase(canvas);
 *     const defocus = (1 - hov) * idleSweep + hov * readerDefocus;
 *
 * Owns its own listeners so the easing survives — see the note at the top.
 *
 * @param {HTMLElement} el
 * @param {object} [options]
 * @param {number} [options.k=0.07] easing rate
 * @param {boolean} [options.pointerOnly=true] ignore touch, where "hover" has no
 *   meaning and the value would stick at 1 after a tap
 */
export function* hoverEase(el, {k = 0.07, pointerOnly = true} = {}) {
  let target = 0;
  let value = 0;
  const enter = (e) => {
    if (pointerOnly && e.pointerType === "touch") return;
    target = 1;
  };
  const leave = () => (target = 0);
  el.addEventListener("pointerenter", enter);
  el.addEventListener("pointerleave", leave);
  try {
    while (true) {
      value = ease(value, target, k);
      yield value;
    }
  } finally {
    // The runtime calls .return() when the cell is disposed or recomputed.
    el.removeEventListener("pointerenter", enter);
    el.removeEventListener("pointerleave", leave);
  }
}

/**
 * Hover easing as a plain mutable holder rather than a generator cell.
 *
 * Prefer this over {@link hoverEase} in a scene. A generator yields once per
 * animation frame, and *every* cell that reads it re-runs at that rate — so a
 * render cell depending on both a 30 fps frame counter and a 60 fps hover value
 * ends up rendering at 60. Reading `.value` off a holder instead leaves the
 * render rate entirely to the frame counter.
 *
 * The easing runs on its own cheap rAF loop.
 *
 * @param {HTMLElement} el
 * @returns {{value: number, stop(): void}} `value` eases 0 → 1 on hover
 */
export function hoverTracker(el, {k = 0.07, pointerOnly = true} = {}) {
  let target = 0;
  const state = {value: 0, stop: () => {}};

  const enter = (e) => {
    if (pointerOnly && e.pointerType === "touch") return;
    target = 1;
  };
  const leave = () => (target = 0);
  el.addEventListener("pointerenter", enter);
  el.addEventListener("pointerleave", leave);

  let raf = requestAnimationFrame(function loop() {
    state.value = ease(state.value, target, k);
    raf = requestAnimationFrame(loop);
  });

  state.stop = () => {
    cancelAnimationFrame(raf);
    el.removeEventListener("pointerenter", enter);
    el.removeEventListener("pointerleave", leave);
  };
  return state;
}

/**
 * Whether `el` is on screen, so an animation can stop when scrolled away.
 *
 * Gate a frame counter with it: `frames({playing: visible})`.
 */
export async function* whenVisible(el, {threshold = 0.05} = {}) {
  if (typeof IntersectionObserver === "undefined") {
    yield true;
    return;
  }

  // Yields ONLY when visibility actually changes, by awaiting between yields.
  //
  // This matters more than it looks. A generator that yields on every pull
  // produces a new value every animation frame, and the runtime then rebuilds
  // every cell that depends on it — including other generators, which restart
  // from scratch. Gating `frames({playing: visible})` on a per-frame generator
  // therefore freezes it at its first value, because it is recreated before it
  // can ever advance.
  //
  // Start optimistic, too: an observer that never fires, or an element still
  // detached when the cell runs (the normal case, since Observable attaches a
  // cell's node only after it returns), would otherwise freeze the widget for
  // good. At worst this animates slightly off-screen.
  let visible = true;
  let wake;
  let changed = new Promise((resolve) => (wake = resolve));

  const io = new IntersectionObserver((entries) => {
    const next = entries[entries.length - 1].isIntersecting;
    if (next === visible) return;
    visible = next;
    const resolve = wake;
    changed = new Promise((r) => (wake = r));
    resolve();
  }, {threshold});
  io.observe(el);

  try {
    while (true) {
      yield visible;
      await changed;
    }
  } finally {
    io.disconnect();
  }
}

/**
 * A 0→1→0 triangle wave with the given period, for sweeping a parameter back and
 * forth (the reference widgets scan a probe this way on a 6 s period).
 */
export function* pingPong(periodMs = 6000, options = {}) {
  const t0 = performance.now();
  for (const _ of frames(options)) {
    const phase = ((performance.now() - t0) % periodMs) / periodMs;
    yield phase < 0.5 ? phase * 2 : 2 - phase * 2;
  }
}

/** A smooth 0→1→0 sweep — `pingPong` without the corners. */
export function* sweep(periodMs = 6000, options = {}) {
  const t0 = performance.now();
  for (const _ of frames(options)) {
    const phase = ((performance.now() - t0) % periodMs) / periodMs;
    yield 0.5 - 0.5 * Math.cos(phase * 2 * Math.PI);
  }
}

/**
 * Adaptive detail: measure how long a frame took and trade quality for
 * smoothness, the trick `diffraction-sim.js` uses while dragging.
 *
 * ```js
 * const quality = qualityBudget({min: 32, max: 256});
 * // each frame:
 * const t0 = performance.now();
 * render(quality.value);
 * quality.report(performance.now() - t0);
 * ```
 *
 * @param {object} [options]
 * @param {number} [options.min=32] floor on the detail level
 * @param {number} [options.max=256] ceiling
 * @param {number} [options.slowMs=90] above this, cut detail hard
 * @param {number} [options.fastMs=30] below this, creep back up
 */
export function qualityBudget({min = 32, max = 256, slowMs = 90, fastMs = 30} = {}) {
  let value = max;
  return {
    get value() {
      return value;
    },
    report(elapsedMs) {
      if (elapsedMs > slowMs) value = Math.max(min, Math.round(value * 0.7));
      else if (elapsedMs < fastMs) value = Math.min(max, Math.round(value * 1.15) + 1);
      return value;
    },
    reset() {
      value = max;
    }
  };
}
