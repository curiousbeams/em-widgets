// Loading large arrays into the browser.
//
// The lab's datasets range from 150 KB to 300 MB. Anything past a few MB should
// not be a single blob the page downloads before it can draw anything, so we
// store arrays as zarr in a ZipStore: one file to host, but chunked inside, and
// read lazily over HTTP range requests.
//
// Write the files with tools/npy-to-zarr-zip.py. Two rules there matter here:
// the zip must be STORED (not deflated) so entries stay range-readable, and the
// zarr codec must be one the browser can decode.

const ZARRITA = "https://cdn.jsdelivr.net/npm/zarrita@0.7.5/+esm";
const ZARRITA_ZIP = "https://cdn.jsdelivr.net/npm/@zarrita/storage@0.2.0/dist/src/zip.js/+esm";

/**
 * Apply the host's URL rewrites to a data URL.
 *
 * Notebooks reference data by absolute published URL (Observable Desktop cannot
 * see anything outside a notebook's own folder). `rewriteImports` on the widget
 * handles ES imports, but a plain `fetch` needs the same treatment — otherwise a
 * local checkout silently loads production data, or 404s on data that has not
 * been published yet.
 *
 * Outside the widget the global is unset and this is the identity, so notebooks
 * using it still work unchanged in Observable Desktop.
 */
export function resolveAsset(url) {
  const rewrites = globalThis.__emWidgetsRewrite;
  if (rewrites) {
    for (const [from, to] of Object.entries(rewrites)) {
      if (String(url).startsWith(from)) return to + String(url).slice(from.length);
    }
  }
  return url;
}

let modules;

/** Load zarrita once, lazily — widgets that never touch data never pay for it. */
async function zarrita() {
  return (modules ??= Promise.all([import(ZARRITA), import(ZARRITA_ZIP)]).then(
    ([zarr, zip]) => ({zarr, ZipFileStore: zip.default ?? zip.ZipFileStore})
  ));
}

/**
 * Open a zarr array stored inside a `.zarr.zip`.
 *
 * @param {string} url location of the zip, absolute or relative to the notebook
 * @param {string} [path="data"] array path within the zarr group
 * @returns {Promise<object>} a zarrita array: `.shape`, `.dtype`, `.chunks`
 *
 * @example
 * const potential = await openZarrZip("../data/fcc-slab-potential.zarr.zip");
 * const slice = await readSlice(potential, 0);   // just that chunk is fetched
 */
export async function openZarrZip(url, path = "data") {
  const {zarr, ZipFileStore} = await zarrita();
  const store = ZipFileStore.fromUrl(new URL(resolveAsset(url), import.meta.url).href);
  return zarr.open.v3(zarr.root(store).resolve(path), {kind: "array"}).catch(() =>
    zarr.open.v2(zarr.root(store).resolve(path), {kind: "array"})
  );
}

/**
 * Read a whole zarr array into memory.
 *
 * @returns {Promise<{data: TypedArray, shape: number[]}>}
 */
export async function readAll(array) {
  const {zarr} = await zarrita();
  const chunk = await zarr.get(array);
  return {data: chunk.data, shape: chunk.shape};
}

/**
 * Read one index along the leading axis — e.g. one slice of a potential stack,
 * or one diffraction pattern out of a 4D dataset.
 *
 * With one chunk per leading index (the default from npy-to-zarr-zip.py) this
 * fetches exactly one chunk.
 *
 * @param {object} array from {@link openZarrZip}
 * @param {number} index
 * @returns {Promise<{data: TypedArray, shape: number[]}>}
 */
export async function readSlice(array, index) {
  const {zarr} = await zarrita();
  const selection = [index, ...Array(array.shape.length - 1).fill(null)];
  const chunk = await zarr.get(array, selection);
  return {data: chunk.data, shape: chunk.shape};
}

/**
 * Read an arbitrary selection, using zarrita's slicing.
 *
 * @example
 * const {slice} = await zarrSlice();
 * await read(array, [0, slice(0, 64), slice(0, 64)]);
 */
export async function read(array, selection) {
  const {zarr} = await zarrita();
  const chunk = await zarr.get(array, selection);
  return {data: chunk.data, shape: chunk.shape};
}

/** zarrita's `slice` helper, for building selections. */
export async function zarrSlice() {
  const {zarr} = await zarrita();
  return {slice: zarr.slice};
}

/**
 * Convert a typed array to Float64Array, which is what the rest of the kit
 * expects. Cheap for anything up to a few million elements.
 */
export function toFloat64(data) {
  return data instanceof Float64Array ? data : Float64Array.from(data);
}
