#!/usr/bin/env python
"""Regenerate the reference fixtures the JS numerics are tested against.

The values come from the lab's own Python, so `npm test` verifies that the
JS port still agrees with the notebooks it was ported from.

Run this from a checkout of one of the science repos, using its venv, so that
py4DSTEM / colorspacious / the local `ctf` and `aberration_utils` modules are
importable. For example:

    cd ~/Documents/myst-sites/presentation-20260812-streaming-ptychography
    .venv/bin/python ~/Documents/em-widgets/tools/generate-fixtures.py \
        --lecture ~/Documents/myst-sites/lecture-ap3402-electron-microscopy-instrumentation \
        --out ~/Documents/em-widgets/test/fixtures
"""

import argparse
import base64
import json
import os
import sys

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--notebooks", default="notebooks",
                        help="directory holding the `ctf` package (default: ./notebooks)")
    parser.add_argument("--lecture", required=True,
                        help="path to the lecture-ap3402 repo (for aberration_utils.py)")
    parser.add_argument("--out", required=True, help="output directory for the fixtures")
    args = parser.parse_args()

    sys.path.insert(0, args.notebooks)
    sys.path.insert(0, os.path.join(args.lecture, "notebooks"))

    from ctf.utils import electron_wavelength_angstrom, radially_average_ctf
    from ctf.visualize import histogram_scaling
    import aberration_utils as au

    os.makedirs(args.out, exist_ok=True)
    write_jch(args.out)
    write_numerics(args.out, electron_wavelength_angstrom, radially_average_ctf,
                   histogram_scaling, au)
    write_colormaps(args.out)
    write_abtem_probe(args.out)


def write_jch(out: str) -> None:
    """JCh -> sRGB1, the colour step inside complex_to_rgb."""
    from colorspacious import cspace_convert

    amps = np.linspace(0.0, 1.0, 11)
    phases = np.linspace(-np.pi, np.pi, 13)
    A, P = np.meshgrid(amps, phases, indexing="ij")
    amp = A.clip(1e-16, 1)
    J = amp * 61.5
    C = np.minimum(98 * J / 123, 110)
    h = np.rad2deg(P) + 180
    JCh = np.stack((J, C, h), axis=-1)
    rgb = cspace_convert(JCh, "JCh", "sRGB1")

    path = os.path.join(out, "jch-to-srgb.json")
    with open(path, "w") as f:
        json.dump({
            "note": "colorspacious cspace_convert(JCh,'JCh','sRGB1'), UNCLIPPED",
            "amps": amps.tolist(),
            "phases": phases.tolist(),
            "JCh": JCh.reshape(-1, 3).tolist(),
            "rgb": rgb.reshape(-1, 3).tolist(),
        }, f)
    print("wrote", path)


def write_numerics(out, electron_wavelength_angstrom, radially_average_ctf,
                   histogram_scaling, au) -> None:
    data = {}

    data["wavelength"] = [
        [float(E), electron_wavelength_angstrom(E)]
        for E in (20e3, 60e3, 80e3, 100e3, 200e3, 300e3)
    ]

    gpts, sampling, energy = (16, 16), (0.2, 0.25), 80e3
    kxa, kya = au.spatial_frequencies(gpts, sampling)
    k, phi = au.polar_coordinates(kxa, kya)
    lam = au.electron_wavelength(energy)
    alpha = k * lam

    coefs = {
        "C10": -50.0, "C12": 12.0, "phi12": 0.3,
        "C21": 300.0, "phi21": -0.7, "C23": 150.0, "phi23": 1.1,
        "C30": 10000.0, "C32": 500.0, "phi32": 0.2, "C34": 250.0, "phi34": -1.3,
        "C41": 50.0, "phi41": 0.4, "C43": 40.0, "phi43": 0.1, "C45": 30.0, "phi45": -0.2,
        "C50": 1e5, "C52": 900.0, "phi52": 0.6, "C54": 800.0, "phi54": -0.9,
        "C56": 700.0, "phi56": 1.4,
    }
    data["chi"] = {
        "gpts": list(gpts), "sampling": list(sampling), "energy": energy,
        "wavelength": float(lam), "coefs": coefs,
        "values": au.aberration_surface(alpha, phi, lam, coefs).ravel().tolist(),
    }
    data["grid"] = {
        "gpts": list(gpts), "sampling": list(sampling),
        "kx": kxa.ravel().tolist(), "ky": kya.ravel().tolist(),
        "k": k.ravel().tolist(), "phi": phi.ravel().tolist(),
    }

    ang = (float(np.abs(kxa).max() * lam * 1e3 / (gpts[0] // 2)),
           float(np.abs(kya).max() * lam * 1e3 / (gpts[1] // 2)))
    data["aperture"] = {
        "semiangle": 20.0, "angular_sampling": list(ang),
        "soft": au.soft_aperture_function(alpha, phi, 20.0, ang).ravel().tolist(),
        "hard": au.hard_aperture_function(alpha, 20.0).ravel().tolist(),
    }

    arr = np.random.default_rng(0).normal(size=(8, 8)) * 3 + 1
    data["histogramScaling"] = {
        "input": arr.ravel().tolist(),
        "output": histogram_scaling(arr).ravel().tolist(),
    }

    ctf_arr = np.abs(np.random.default_rng(1).normal(size=(16, 20)))
    kb, ib = radially_average_ctf(ctf_arr, (0.2, 0.25))
    data["radialAverage"] = {
        "shape": [16, 20], "sampling": [0.2, 0.25],
        "input": ctf_arr.ravel().tolist(), "k": kb.tolist(), "I": ib.tolist(),
    }

    rng = np.random.default_rng(2)
    data["fft2"] = {}
    for name, (nx, ny) in {"pow2": (8, 8), "nonpow2": (6, 10)}.items():
        a = rng.normal(size=(nx, ny)) + 1j * rng.normal(size=(nx, ny))
        F = np.fft.fft2(a)
        data["fft2"][name] = {
            "shape": [nx, ny],
            "re": a.real.ravel().tolist(), "im": a.imag.ravel().tolist(),
            "Fre": F.real.ravel().tolist(), "Fim": F.imag.ravel().tolist(),
        }

    path = os.path.join(out, "numerics.json")
    with open(path, "w") as f:
        json.dump(data, f)
    print("wrote", path)


def write_abtem_probe(out: str) -> None:
    """Reference probe from abtem — the end-to-end FFT-convention test.

    Requires abtem, so run this from a repo whose venv has it (the workshop or
    lecture repos do). Skipped with a warning otherwise.
    """
    try:
        from abtem.transfer import CTF
    except ImportError:
        print("warning: abtem not installed, skipping the probe fixture")
        return

    gpts, sampling, energy = (64, 64), (0.25, 0.25), 300e3
    # abtem's `defocus` is -C10; everything else maps straight to Krivanek Cnm.
    ctf = CTF(
        semiangle_cutoff=20,
        sampling=sampling,
        gpts=gpts,
        energy=energy,
        aberration_coefficients=dict(
            defocus=100.0,
            astigmatism=50.0,
            astigmatism_angle=np.deg2rad(45),
            coma=1e4,
            coma_angle=np.pi,
        ),
    )
    fourier = np.asarray(ctf.to_diffraction_patterns(gpts=ctf.gpts).array)
    real = np.asarray(ctf.to_point_spread_functions(gpts=ctf.gpts, extent=ctf.extent).array)

    def b64(a):
        return base64.b64encode(np.ascontiguousarray(a, dtype=np.float32).tobytes()).decode()

    path = os.path.join(out, "abtem-probe.json")
    with open(path, "w") as f:
        json.dump({
            "note": "abtem CTF; note abtem's defocus == -C10.",
            "gpts": list(gpts), "sampling": list(sampling), "energy": energy,
            "semiangle": 20.0,
            "coefficients": {
                "C10": -100.0, "C12": 50.0, "phi12": float(np.deg2rad(45)),
                "C21": 1e4, "phi21": float(np.pi),
            },
            "wavelength": float(ctf.wavelength),
            "angular_sampling": list(ctf.angular_sampling),
            "fourier_re": b64(fourier.real), "fourier_im": b64(fourier.imag),
            "real_re": b64(real.real), "real_im": b64(real.imag),
        }, f)
    print("wrote", path)


def write_colormaps(out: str) -> None:
    """Rewrite kit/colormap-data.js from matplotlib and cmasher."""
    import matplotlib

    names = {"magma": "magma", "twilight": "twilight",
             "RdBu": "RdBu", "PuOr": "PuOr", "PiYG": "PiYG"}
    maps = {}
    for js_name, mpl_name in names.items():
        colors = matplotlib.colormaps[mpl_name](np.linspace(0, 1, 256))[:, :3]
        maps[js_name] = base64.b64encode(
            np.round(colors * 255).astype(np.uint8).tobytes()).decode()
    try:
        import cmasher  # noqa: F401  (registers the cmr.* colormaps)
        colors = matplotlib.colormaps["cmr.eclipse"](np.linspace(0, 1, 256))[:, :3]
        maps["eclipse"] = base64.b64encode(
            np.round(colors * 255).astype(np.uint8).tobytes()).decode()
    except ImportError:
        print("warning: cmasher not installed, skipping the eclipse colormap")

    path = os.path.join(os.path.dirname(out.rstrip("/")), "..", "kit", "colormap-data.js")
    path = os.path.normpath(path)
    lines = [
        "// Colormap lookup tables: 256 RGB triplets each, base64-encoded bytes.",
        "// Generated from matplotlib / cmasher to match the lab's Python figures exactly.",
        "// Regenerate with: tools/generate-fixtures.py",
        "",
        "export const COLORMAP_DATA = {",
    ]
    for name, encoded in maps.items():
        lines.append(f'  {name}:\n    "{encoded}",')
    lines.append("};")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("wrote", path)


if __name__ == "__main__":
    main()
