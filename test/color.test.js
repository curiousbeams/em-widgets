import test from "node:test";
import assert from "node:assert/strict";
import {readFileSync} from "node:fs";
import {fileURLToPath} from "node:url";

import {jchToSrgb, colormap, applyColormap, listColormaps} from "../kit/color.js";
import {histogramScaling, radialAverage} from "../kit/image.js";

const read = (name) =>
  JSON.parse(readFileSync(fileURLToPath(new URL(`./fixtures/${name}`, import.meta.url)), "utf8"));

const jch = read("jch-to-srgb.json");
const num = read("numerics.json");

test("jchToSrgb matches colorspacious cspace_convert(JCh -> sRGB1)", (t) => {
  let worst = 0, at = -1;
  for (let i = 0; i < jch.JCh.length; i++) {
    const [J, C, h] = jch.JCh[i];
    const got = jchToSrgb(J, C, h);
    for (let c = 0; c < 3; c++) {
      const d = Math.abs(got[c] - jch.rgb[i][c]);
      if (d > worst) { worst = d; at = i; }
    }
  }
  t.diagnostic(`${jch.JCh.length} samples, max |diff| = ${worst.toExponential(3)}`);
  assert.ok(
    worst < 1e-9,
    `max |diff| ${worst.toExponential(3)} at sample ${at} ` +
      `(JCh ${JSON.stringify(jch.JCh[at])}: got ${JSON.stringify(jchToSrgb(...jch.JCh[at]))}, ` +
      `want ${JSON.stringify(jch.rgb[at])})`
  );
});

test("jchToSrgb is well-behaved at the achromatic axis", () => {
  for (const h of [0, 90, 180, 270, 359.9]) {
    const rgb = jchToSrgb(0, 0, h);
    for (const v of rgb) assert.ok(Number.isFinite(v), `J=0 C=0 h=${h} gave ${rgb}`);
  }
});

test("histogramScaling matches ctf/visualize.py", (t) => {
  const got = histogramScaling(num.histogramScaling.input);
  const want = num.histogramScaling.output;
  let worst = 0;
  for (let i = 0; i < want.length; i++) worst = Math.max(worst, Math.abs(got[i] - want[i]));
  t.diagnostic(`max |diff| = ${worst.toExponential(3)}`);
  assert.ok(worst < 1e-12, `max |diff| ${worst}`);
});

test("radialAverage matches ctf/utils.py:radially_average_ctf", (t) => {
  const f = num.radialAverage;
  const [nx, ny] = f.shape;
  const {k, I} = radialAverage(f.input, nx, ny, f.sampling);

  assert.equal(k.length, f.k.length, `bin count: got ${k.length}, want ${f.k.length}`);
  let worstK = 0, worstI = 0;
  for (let i = 0; i < f.k.length; i++) {
    worstK = Math.max(worstK, Math.abs(k[i] - f.k[i]));
    worstI = Math.max(worstI, Math.abs(I[i] - f.I[i]));
  }
  t.diagnostic(`${k.length} bins, max |dk| = ${worstK.toExponential(2)}, max |dI| = ${worstI.toExponential(2)}`);
  assert.ok(worstK < 1e-12 && worstI < 1e-12, `dk ${worstK}, dI ${worstI}`);
});

test("colormaps decode to 256 RGB triplets", () => {
  for (const name of listColormaps()) {
    const table = colormap(name);
    assert.equal(table.length, 768, `${name} length`);
  }
  assert.throws(() => colormap("viridis"), /Unknown colormap/);
});

test("applyColormap spans the table and writes opaque RGBA", () => {
  const rgba = applyColormap([0, 0.5, 1], {colormap: "magma"});
  const table = colormap("magma");
  assert.deepEqual([...rgba.slice(0, 3)], [...table.slice(0, 3)]);
  assert.deepEqual([...rgba.slice(8, 11)], [...table.slice(255 * 3, 255 * 3 + 3)]);
  assert.equal(rgba[3], 255);
});
