// @curiousbeams/em-kit — shared client-side primitives for electron-microscopy
// teaching widgets.
//
// Import the whole kit:
//
//     import * as em from "npm:@curiousbeams/em-kit";
//
// or, from a notebook in this repo (works in Observable Desktop and when the
// notebook is embedded in a MyST page, because the wrapper resolves relative
// imports against the notebook's own URL):
//
//     import {chi, complexProbe} from "../kit/index.js";
//
// Conventions
// -----------
//  * Lengths are Angstroms, energies eV, angles radians (except where a name
//    says mrad).
//  * 2D arrays are flat and row-major, indexed `[ix * ny + iy]` — numpy's
//    `indexing="ij"`. `ix` is the vertical axis, as in `imshow`.
//  * Complex arrays are `{re, im}` pairs of Float64Arrays.
//  * Reciprocal-space grids are corner-centered (fftfreq order). Only apply
//    `fftshift` when you are about to draw.
//
// Every numeric primitive here is checked against the lab's Python
// (numpy / py4DSTEM / colorspacious) in test/ — run `npm test`.

export * from "./units.js";
export * from "./grid.js";
export * from "./fft.js";
export * from "./optics.js";
export * from "./image.js";
export * from "./color.js";
export * from "./canvas.js";
export * from "./raytrace.js";
export * from "./scene3d.js";
export * from "./anim.js";
export * from "./ui.js";
export * from "./specimen.js";
export * from "./detectors.js";
export * from "./paraxial.js";
export * from "./nonparaxial.js";
export * from "./projections.js";
export * from "./ptycho.js";
export * from "./ptycho-sim.js";
export * from "./epie.js";
