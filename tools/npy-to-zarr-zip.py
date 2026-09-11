#!/usr/bin/env python
"""Convert .npy / .h5 arrays into a zarr ZipStore for the browser to read lazily.

Why a zip of zarr rather than a raw binary blob:

  * one file to host and version, instead of thousands of chunk files;
  * the browser fetches only the chunks a widget actually touches, via HTTP
    range requests, so a 20 MB dataset costs a few hundred KB to start;
  * dtype, shape and chunking travel with the data, so the widget does not
    hard-code them.

Two rules matter for browser use:

  1. The zip itself must be STORED, not deflated. Zarr chunks are already
     compressed, and a deflated zip cannot be range-read entry by entry.
  2. Pick a codec the browser can actually decode. `zstd` and `blosc` need a
     wasm build of numcodecs; `gzip` is decoded natively. Default to zstd and
     fall back with --codec gzip if you hit trouble.

Example:

    python tools/npy-to-zarr-zip.py \
        ~/…/notebooks/data/FCC-slab-potential-7x244x242-float32.npy \
        --shape 7 244 242 --dtype float32 \
        --chunks 1 244 242 \
        --out data/fcc-slab-potential.zarr.zip
"""

import argparse
import os
import zipfile

import numpy as np


def load_array(path: str, shape, dtype: str) -> np.ndarray:
    if path.endswith(".npy"):
        try:
            return np.load(path)
        except ValueError:
            # Several of the lab's "npy" files are raw buffers written with
            # ndarray.tofile(), so they need an explicit shape and dtype.
            if not shape:
                raise SystemExit(
                    f"{path} is a raw buffer, not a real .npy — pass --shape and --dtype"
                )
            return np.fromfile(path, dtype=dtype).reshape(shape)
    if path.endswith((".h5", ".hdf5")):
        import h5py

        with h5py.File(path, "r") as f:
            keys = list(f.keys())
            raise SystemExit(f"pass --dataset for an HDF5 file; top-level keys: {keys}")
    raise SystemExit(f"unsupported input: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input")
    parser.add_argument("--out", required=True, help="output .zarr.zip path")
    parser.add_argument("--shape", nargs="*", type=int, default=None,
                        help="shape, for raw buffers saved with tofile()")
    parser.add_argument("--dtype", default="float32")
    parser.add_argument("--chunks", nargs="*", type=int, default=None,
                        help="chunk shape; defaults to one chunk per leading index")
    parser.add_argument("--codec", default="zstd", choices=["zstd", "gzip", "none"])
    parser.add_argument("--cast", default=None,
                        help="cast to this dtype before writing, e.g. float32 or uint8")
    parser.add_argument("--zarr-format", type=int, default=3, choices=[2, 3])
    args = parser.parse_args()

    import zarr
    from zarr.codecs import BytesCodec, GzipCodec, ZstdCodec

    array = load_array(args.input, args.shape, args.dtype)
    if args.cast:
        array = array.astype(args.cast)

    chunks = tuple(args.chunks) if args.chunks else default_chunks(array.shape)
    codecs = [BytesCodec()]
    if args.codec == "zstd":
        codecs.append(ZstdCodec(level=5))
    elif args.codec == "gzip":
        codecs.append(GzipCodec(level=5))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    if os.path.exists(args.out):
        os.remove(args.out)

    # ZIP_STORED so each chunk stays individually range-readable.
    store = zarr.storage.ZipStore(args.out, mode="w", compression=zipfile.ZIP_STORED)
    root = zarr.create_group(store=store, zarr_format=args.zarr_format)
    z = root.create_array(
        name="data",
        shape=array.shape,
        chunks=chunks,
        dtype=array.dtype,
        compressors=codecs[1:] if len(codecs) > 1 else None,
    )
    z[:] = array
    store.close()

    raw = array.nbytes
    packed = os.path.getsize(args.out)
    print(f"{args.input}")
    print(f"  shape {array.shape} {array.dtype}, chunks {chunks}, codec {args.codec}")
    print(f"  {raw/1048576:.2f} MB raw -> {packed/1048576:.2f} MB zipped "
          f"({100*packed/raw:.0f}%)  ->  {args.out}")


def default_chunks(shape):
    """One chunk per leading index for a stack; a single chunk for a 2D array."""
    if len(shape) >= 3:
        return (1,) + tuple(shape[1:])
    return tuple(shape)


if __name__ == "__main__":
    main()
