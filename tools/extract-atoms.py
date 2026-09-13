#!/usr/bin/env python
"""Extract atom positions from a sliced projected-potential stack.

A widget that draws atoms should draw them exactly where the simulation's atoms
are, so pull the positions out of the potential itself rather than rebuilding a
lattice from an assumed lattice constant. Peaks are refined to sub-pixel
precision by centre of mass in a small periodic window.

Example:

    .venv/bin/python tools/extract-atoms.py \\
        notebooks/data/FCC-slab-potential-7x244x242-float32.npy \\
        --shape 7 244 242 --sampling 0.1 --dz 2.857142857142857 \\
        --out ~/Documents/em-widgets/data/fcc-slab-atoms.json
"""

import argparse
import json

import numpy as np
from scipy import ndimage


def refine(slab, row, col, sampling, window=7):
    """Centre-of-mass refinement of one peak, wrapping at the edges."""
    nx, ny = slab.shape
    half = window // 2
    rows = (np.arange(row - half, row + half + 1) % nx)[:, None]
    cols = (np.arange(col - half, col + half + 1) % ny)[None, :]
    patch = slab[rows, cols]
    patch = np.clip(patch - patch.min(), 0, None)
    total = patch.sum()
    if total <= 0:
        return row * sampling, col * sampling
    dr = (patch.sum(1) * np.arange(-half, half + 1)).sum() / total
    dc = (patch.sum(0) * np.arange(-half, half + 1)).sum() / total
    return (row + dr) * sampling, (col + dc) * sampling


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input")
    ap.add_argument("--shape", nargs=3, type=int, required=True, help="slices nx ny")
    ap.add_argument("--sampling", type=float, default=0.1, help="Angstrom per pixel")
    ap.add_argument("--dz", type=float, required=True, help="slice thickness, Angstrom")
    ap.add_argument("--threshold", type=float, default=0.15,
                    help="peak height as a fraction of the slice maximum")
    ap.add_argument("--min-separation", type=int, default=15,
                    help="peak separation in pixels")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    nz, nx, ny = args.shape
    potential = np.fromfile(args.input, dtype=np.float32).reshape(nz, nx, ny)

    atoms = []
    for s in range(nz):
        slab = potential[s]
        peaks = np.argwhere(
            (slab == ndimage.maximum_filter(slab, size=args.min_separation, mode="wrap"))
            & (slab > slab.max() * args.threshold)
        )
        for row, col in peaks:
            x, y = refine(slab, row, col, args.sampling)
            # z measured from the middle of the slab, so the specimen is centred
            atoms.append([round(float(x), 4), round(float(y), 4),
                          round(float((s - (nz - 1) / 2) * args.dz), 4)])

    xs = np.array([a[0] for a in atoms])
    ys = np.array([a[1] for a in atoms])
    extent = [nx * args.sampling, ny * args.sampling]

    # In-plane nearest neighbour within one layer, with periodic wrap — reported
    # so the caller can sanity-check which material this actually is.
    first = np.array([[a[0], a[1]] for a in atoms if abs(a[2] - atoms[0][2]) < 1e-9])
    d = first[:, None, :] - first[None, :, :]
    d -= np.array(extent) * np.round(d / np.array(extent))
    dist = np.linalg.norm(d, axis=-1)
    np.fill_diagonal(dist, np.inf)
    nn = float(np.median(dist.min(1)))

    out = {
        "note": "atom positions refined from the projected-potential peaks",
        "source": args.input,
        "extent": extent,
        "sampling": args.sampling,
        "dz": args.dz,
        "n_layers": nz,
        "n_atoms": len(atoms),
        "in_plane_nn": round(nn, 4),
        "lattice_constant_fcc111": round(nn * np.sqrt(2), 4),
        "positions": atoms,
    }
    with open(args.out, "w") as f:
        json.dump(out, f)

    print(f"{len(atoms)} atoms over {nz} layers -> {args.out}")
    print(f"  in-plane nearest neighbour {nn:.3f} A")
    print(f"  implied FCC lattice constant a = nn*sqrt(2) = {nn * np.sqrt(2):.3f} A")
    print(f"  interlayer spacing dz = {args.dz:.3f} A (FCC(111) expects a/sqrt(3) = "
          f"{nn * np.sqrt(2) / np.sqrt(3):.3f} A)")


if __name__ == "__main__":
    main()
