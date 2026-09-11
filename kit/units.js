// Electron-optical unit conversions.
//
// Ported from `ctf/utils.py:electron_wavelength_angstrom` and
// `aberration_utils.py:electron_wavelength`, which are byte-identical and were
// duplicated in at least four places across the lab's notebooks.
//
// Energies are in eV, lengths in Angstroms.

/** Electron rest mass [kg]. */
export const M_E = 9.109383e-31;
/** Elementary charge [C]. */
export const E_CHARGE = 1.602177e-19;
/** Speed of light [m/s]. */
export const C_LIGHT = 299792458;
/** Planck constant [J s]. */
export const H_PLANCK = 6.62607e-34;

/**
 * Relativistic electron wavelength.
 *
 * @param {number} energy accelerating voltage [eV] (e.g. 300e3 for 300 kV)
 * @returns {number} wavelength [Angstrom]
 */
export function electronWavelength(energy) {
  const lambda =
    H_PLANCK /
    Math.sqrt(2 * M_E * E_CHARGE * energy) /
    Math.sqrt(1 + (E_CHARGE * energy) / 2 / M_E / (C_LIGHT * C_LIGHT));
  return lambda * 1e10; // m -> Angstrom
}

/**
 * Interaction parameter sigma, relating projected potential to phase shift:
 * `phase = sigma * V_projected`.
 *
 * Matches `py4DSTEM.process.utils.electron_interaction_parameter` and
 * `abtem.core.energy.energy2sigma`.
 *
 * @param {number} energy accelerating voltage [eV]
 * @returns {number} sigma [1 / (Angstrom eV)]
 */
export function interactionParameter(energy) {
  const lambda = electronWavelength(energy);
  const mc2 = M_E * C_LIGHT * C_LIGHT;
  return ((2 * Math.PI) / lambda / energy) * (mc2 + E_CHARGE * energy) / (2 * mc2 + E_CHARGE * energy);
}

/** Alias for {@link interactionParameter}, matching abtem's spelling. */
export const energy2sigma = interactionParameter;

/** Alias for {@link electronWavelength}, matching abtem's spelling. */
export const energy2wavelength = electronWavelength;
