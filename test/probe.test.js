import test from "node:test";
import assert from "node:assert/strict";
import {readFileSync} from "node:fs";
import {fileURLToPath} from "node:url";

import {electronWavelength} from "../kit/units.js";
import {complexProbe} from "../kit/optics.js";
import {ifft2, fftshift} from "../kit/fft.js";

// Reference arrays generated with abtem's CTF — see tools/generate-fixtures.py.
// This is the end-to-end convention test: fftfreq ordering, the sign of chi, the
// ifft2 normalisation and the fftshift all have to agree for it to pass.
const f = JSON.parse(
  readFileSync(fileURLToPath(new URL("./fixtures/abtem-probe.json", import.meta.url)), "utf8")
);

function decode(b64) {
  const buffer = Buffer.from(b64, "base64");
  return new Float32Array(buffer.buffer, buffer.byteOffset, buffer.length / 4);
}

function l2(re, im) {
  let sum = 0;
  for (let i = 0; i < re.length; i++) sum += re[i] * re[i] + im[i] * im[i];
  return Math.sqrt(sum);
}

/** Largest complex difference after normalising both arrays to unit L2 norm. */
function maxComplexDiff(ar, ai, br, bi) {
  const na = l2(ar, ai);
  const nb = l2(br, bi);
  let worst = 0;
  for (let i = 0; i < ar.length; i++) {
    worst = Math.max(worst, Math.hypot(ar[i] / na - br[i] / nb, ai[i] / na - bi[i] / nb));
  }
  return worst;
}

const [nx, ny] = f.gpts;
const reference = {
  fourierRe: decode(f.fourier_re),
  fourierIm: decode(f.fourier_im),
  realRe: decode(f.real_re),
  realIm: decode(f.real_im)
};

// abtem's `defocus` is -C10; f.coefficients is already in Krivanek Cnm form.
const probeOptions = {
  gpts: f.gpts,
  sampling: f.sampling,
  energy: f.energy,
  semiangleCutoff: f.semiangle,
  aberrations: f.coefficients
};

test("electronWavelength agrees with abtem to the precision of their constants", (t) => {
  // The kit deliberately uses the lab's Python constants (ctf/utils.py,
  // aberration_utils.py), which are rounded; abtem uses CODATA. The two differ
  // by ~1e-7 relative, which is far below anything that matters here — but it is
  // why this tolerance is not 1e-12.
  const mine = electronWavelength(f.energy);
  const relative = Math.abs(mine - f.wavelength) / f.wavelength;
  t.diagnostic(`kit ${mine.toPrecision(12)} vs abtem ${f.wavelength.toPrecision(12)} (relative ${relative.toExponential(2)})`);
  assert.ok(relative < 1e-6, `relative difference ${relative}`);
});

test("complexProbe matches abtem's Fourier-space probe", (t) => {
  const probe = complexProbe(probeOptions);
  const worst = maxComplexDiff(
    fftshift(probe.re, nx, ny),
    fftshift(probe.im, nx, ny),
    reference.fourierRe,
    reference.fourierIm
  );
  t.diagnostic(`max |Δ| = ${worst.toExponential(2)} (fixture stored as float32)`);
  assert.ok(worst < 1e-6, `max |Δ| ${worst}`);
});

test("ifft2 of the probe matches abtem's real-space probe", (t) => {
  const probe = complexProbe(probeOptions);
  const psi = ifft2(probe, nx, ny);
  const worst = maxComplexDiff(
    fftshift(psi.re, nx, ny),
    fftshift(psi.im, nx, ny),
    reference.realRe,
    reference.realIm
  );
  t.diagnostic(`max |Δ| = ${worst.toExponential(2)}`);
  assert.ok(worst < 1e-6, `max |Δ| ${worst}`);
});

test("the sign of chi is exp(-i chi), not exp(+i chi)", () => {
  // Flipping every aberration coefficient conjugates the probe, which is what an
  // exp(+i chi) convention would produce. It must NOT match abtem.
  const flipped = Object.fromEntries(
    Object.entries(f.coefficients).map(([k, v]) => [k, k.startsWith("C") ? -v : v])
  );
  const probe = complexProbe({...probeOptions, aberrations: flipped});
  const worst = maxComplexDiff(
    fftshift(probe.re, nx, ny),
    fftshift(probe.im, nx, ny),
    reference.fourierRe,
    reference.fourierIm
  );
  assert.ok(worst > 1e-2, `the wrong sign convention should disagree, got ${worst}`);
});

test("a 128x128 probe recomputes fast enough for a slider drag", (t) => {
  const options = {...probeOptions, gpts: [128, 128], sampling: [0.125, 0.125]};
  const start = performance.now();
  const runs = 20;
  for (let i = 0; i < runs; i++) {
    const probe = complexProbe({...options, aberrations: {...f.coefficients, C10: -100 - i}});
    ifft2(probe, 128, 128);
  }
  const ms = (performance.now() - start) / runs;
  t.diagnostic(`${ms.toFixed(2)} ms per probe + ifft2 at 128x128`);
  assert.ok(ms < 50, `too slow for interactive use: ${ms} ms`);
});
