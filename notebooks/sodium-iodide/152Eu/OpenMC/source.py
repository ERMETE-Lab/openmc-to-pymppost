"""
Eu-152 photon source for OpenMC (and matching MCNP SI/SP cards).

Decay data: LNHB/DDEP recommended values, "Table de Radionucleides - Eu-152",
V.R. Vanin, R.M. de Castro, E. Browne (2004), http://www.lnhb.fr/nuclides/Eu-152_tables.pdf

Lines included:
  * every gamma line with emission probability >= 0.1 per 100 decays
    (unresolved lines at the same energy are summed);
  * optionally the K X-rays of Sm and Gd (>= 0.1 per 100 decays).
L X-rays (5.6-7.8 keV), the 511 keV annihilation line (0.054 %) and the
weaker gamma lines are neglected.

One photon is emitted per source history: gamma-gamma and gamma-X
cascades are NOT correlated (no true-coincidence summing). The MCNP SDEF
built from mcnp_cards() samples the source in the same way.
"""

import numpy as np
import openmc

# ================================ CONFIGURATION ================================
INCLUDE_XRAYS = True                # Sm/Gd K X-rays (~40-49 keV)
SOURCE_POS    = (0.0, 0.0, -1.0)    # [cm], 1 cm before the crystal front face
BEAM_DIR      = (0.0, 0.0, 1.0)     # monodirectional beam along +z
# ===============================================================================

# (energy [keV], photons per 100 decays) - gamma lines, DDEP 2004
GAMMA_LINES = [
    (121.7817, 28.41), (244.6974, 7.55),  (295.9387, 0.442), (329.425, 0.129),
    (344.2785, 26.59), (367.7891, 0.862), (411.1165, 2.238), (416.048, 0.1090),
    (443.965,  3.12),  # Sm 0.32 + 2.80, same energy
    (488.6792, 0.4139), (503.474, 0.1533), (563.990, 0.457), (566.442, 0.131),
    (586.265,  0.462), (656.489, 0.1437), (674.675, 0.1871),  # Sm 0.170 + Gd 0.0171
    (678.623,  0.470), (688.670, 0.841),
    (719.349,  0.327),  # Sm 0.268 + 0.059, same energy
    (764.900,  0.190), (778.9045, 12.97), (810.451, 0.317), (841.574, 0.163),
    (867.380,  4.243), (919.337, 0.429), (926.317, 0.273), (963.390, 0.1341),
    (964.079, 14.50),  (1005.272, 0.665), (1084.0, 0.244), (1085.837, 10.13),
    (1089.737, 1.73),  (1109.174, 0.186), (1112.076, 13.41), (1212.948, 1.416),
    (1249.938, 0.186), (1292.778, 0.104), (1299.142, 1.633), (1408.013, 20.85),
    (1457.643, 0.498), (1528.103, 0.281),
]

# (energy [keV], photons per 100 decays) - K X-rays, DDEP 2004.
# Kbeta groups are represented by their strongest component.
XRAY_LINES = [
    (39.5229, 20.8),  (40.1186, 37.7),   # Sm K-alpha2, K-alpha1
    (45.413,  11.78), (46.575,  3.04),   # Sm K'beta1, K'beta2 groups
    (42.3093, 0.243), (42.9967, 0.437),  # Gd K-alpha2, K-alpha1
    (48.697,  0.138),                    # Gd K'beta1 group
]

lines = GAMMA_LINES + (XRAY_LINES if INCLUDE_XRAYS else [])
energies_eV = np.array([e for e, _ in lines]) * 1e3
yields      = np.array([p for _, p in lines]) / 100.0      # photons per decay

PHOTONS_PER_DECAY = yields.sum()     # multiply per-source-photon results by this
probabilities     = yields / PHOTONS_PER_DECAY             # sampling probabilities

source = openmc.IndependentSource(
    space=openmc.stats.Point(SOURCE_POS),
    angle=openmc.stats.Monodirectional(BEAM_DIR),
    energy=openmc.stats.Discrete(energies_eV, probabilities),
    particle="photon",
)


def mcnp_cards(dist=1, width=78):
    """SI/SP cards with exactly the same lines (energies in MeV)."""
    def block(card, values, fmt):
        out, line = [], card
        for v in values:
            item = " " + format(v, fmt)
            if len(line) + len(item) > width:
                out.append(line)
                line = "    "
            line += item
        out.append(line)
        return "\n".join(out)
    return (block(f"si{dist} L", energies_eV / 1e6, ".7f") + "\n" +
            block(f"sp{dist}  ", yields, ".6e"))


if __name__ == "__main__":
    print(f"Lines: {len(lines)} | photons per decay represented: {PHOTONS_PER_DECAY:.4f}")
    print(f"Fraction of lines above 1.022 MeV (pair production possible): "
          f"{probabilities[energies_eV > 1.022e6].sum():.3f}")
    print()
    print(mcnp_cards())