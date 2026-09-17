#!/usr/bin/env python
"""Build the Si(111)-(7x7) DAS reconstruction and write it out for the widgets.

This is the construction from the ERC interview material
(`grants/ERC-StG-2026/interview/silicon-reconstruction.ipynb`), which builds the
dimer-adatom-stacking-fault model in the four steps it is usually explained in:

    0. a bulk-terminated bilayer, repeated 7 x 7
    1. remove the atom at the cell corner        -> the corner hole
    2. displace the atoms along the cell edges   -> the dimers
    3. mirror one half of the cell               -> the stacking fault
    4. add twelve adatoms                        -> the adatoms

Written out as JSON: the two primitive vectors of the 7 x 7 cell, and two bases
on that same cell — the reconstruction and the bulk-terminated surface it forms
on. A widget tiles whichever of the two belongs at each cell.

Needs ASE, which lives in the ERC project's environment:

    /Users/gvarnavides/Documents/grants/ERC-StG-2026/interview/.venv/bin/python \\
        tools/si111-7x7.py --out data/si111-7x7.json
"""

import argparse
import json

import ase.build
import numpy as np


def bilayer():
    """One bulk-terminated Si(111) bilayer, rotated so its cell sits on x."""
    bulk = ase.build.bulk("Si", "diamond")
    si_111 = ase.build.surface(bulk, (1, 1, 1), layers=3, periodic=True)
    si_111.wrap()
    vacuum = si_111.cell[2, 2]
    slab = ase.build.surface(bulk, (1, 1, 1), layers=2, vacuum=vacuum)
    del slab[(slab.positions[:, 2] < 11) | (slab.positions[:, 2] > 14.5)]
    # The notebook marks one of the two sublattices so it can pick adatom sites
    # off it later. Kept, for the same reason.
    slab.numbers[0] = 15
    slab.rotate(-30, "z", rotate_cell=True)
    return slab


def reconstruct(step0):
    """Steps 1 to 4, on a 7 x 7 supercell of the bilayer."""
    layer_spacing = float(np.diff(np.unique(step0.positions[:, 2].round(6)))[0])

    # 1. The corner hole.
    step1 = step0.copy()
    scaled = np.round(step1.get_scaled_positions(), decimals=4)
    del step1[(scaled[:, 0] == 0) & (scaled[:, 1] == 0)]

    # 2. The dimers, along the three cell edges.
    step2 = step1.copy()
    scaled = step2.get_scaled_positions()
    rounded = np.round(scaled, decimals=8)
    along_a = np.where(rounded[:, 0] == 0)[0]
    along_b = np.where(rounded[:, 1] == 0)[0]
    along_diagonal = np.where((rounded[:, 0] + rounded[:, 1]) == 1)[0]
    d = 1 / 7 / 8
    scaled[along_a[::2], 1] += d
    scaled[along_a[1::2], 1] -= d
    scaled[along_b[::2], 0] += d
    scaled[along_b[1::2], 0] -= d
    scaled[along_diagonal[::2], 0] += d
    scaled[along_diagonal[::2], 1] -= d
    scaled[along_diagonal[1::2], 0] -= d
    scaled[along_diagonal[1::2], 1] += d
    step2.set_scaled_positions(scaled)

    # 3. The stacking fault: one half of the cell is the mirror of the other.
    step3 = step2.copy()
    rounded = np.round(step3.get_scaled_positions(), decimals=8)
    mirror_line = step3[(rounded[:, 0] + rounded[:, 1]) == 1]
    left = step3[(rounded[:, 0] + rounded[:, 1]) < 1]
    right = left.copy()
    right.positions[:, 0] = left.cell[:2, :2].sum(0)[0] - right.positions[:, 0]
    right.wrap()
    rounded = np.round(right.get_scaled_positions(), decimals=8)
    del right[(rounded[:, 0] + rounded[:, 1]) < 1]
    step3 = left + mirror_line + right

    # 4. Twelve adatoms, six on each half, one layer above the surface.
    step4 = step3.copy()
    rounded = np.round(step4.get_scaled_positions(), decimals=8)
    left = step4[(rounded[:, 0] + rounded[:, 1]) < 1]
    marked = left[left.numbers == 15]
    left_adatoms = marked[[6, 8, 10, 17, 19, 24]]
    right_adatoms = left_adatoms.copy()
    right_adatoms.positions[:, 0] = (
        right_adatoms.cell[:2, :2].sum(0)[0] - right_adatoms.positions[:, 0]
    )
    adatoms = left_adatoms + right_adatoms
    adatoms.numbers[:] = 13
    adatoms.positions[:, 2] += 2 * layer_spacing
    return step4 + adatoms


def basis(atoms, top):
    """Positions relative to the cell origin, with depth measured from the top.

    Each atom carries which half of the cell it belongs to, 0 or 1. The 7 x 7
    cell is a rhombus and its two triangular halves are the faulted and unfaulted
    ones, so the diagonal between them is a direction the surface itself uses: an
    island that ends there ends on a straight edge, where one that has to stop at
    whole rhombi comes out stepped. The dimer row lies exactly on the diagonal and
    is counted with the first half, so that two halves meeting there do not place
    it twice.
    """
    scaled = atoms.get_scaled_positions()
    return [
        [round(float(x), 4), round(float(y), 4), round(float(top - z), 4),
         0 if (s[0] + s[1]) <= 1 + 1e-6 else 1]
        for (x, y, z), s in zip(atoms.positions, scaled)
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    step0 = bilayer().repeat((7, 7, 1))
    step4 = reconstruct(step0)
    top = float(max(step0.positions[:, 2].max(), step4.positions[:, 2].max()))
    cell = step4.cell[:2, :2]

    data = {
        "comment": "Si(111)-(7x7) DAS reconstruction and the bulk-terminated "
                   "bilayer it forms on, sharing one 7x7 cell. Angstrom. "
                   "Each atom is [x, y, depth below the topmost atom, which "
                   "half of the cell it is in].",
        "a1": [round(float(cell[0, 0]), 4), round(float(cell[0, 1]), 4)],
        "a2": [round(float(cell[1, 0]), 4), round(float(cell[1, 1]), 4)],
        "das": basis(step4, top),
        "bulk": basis(step0, top)
    }
    with open(args.out, "w") as file:
        json.dump(data, file)
    halves = [sum(1 for a in data["das"] if a[3] == h) for h in (0, 1)]
    print(f"{args.out}: {len(data['das'])} atoms reconstructed ({halves[0]} + "
          f"{halves[1]} by half), {len(data['bulk'])} bulk-terminated, "
          f"cell {np.linalg.norm(cell[0]):.3f} x {np.linalg.norm(cell[1]):.3f} A")


if __name__ == "__main__":
    main()
