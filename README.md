# 🔬 OpenMC → pyMPPost

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
<!-- TODO: Zenodo DOI and arXiv badges, once available -->

**Detector response simulations with OpenMC and pyMPPost.**
This repository converts the particle tracks written by [OpenMC](https://openmc.org) into an
MCNPX-PoliMi collision file, so that the detector response post-processor
[pyMPPost](https://gitlab.eecs.umich.edu/umich-dnng/pymppost) can be used with OpenMC,
without any modification of OpenMC or pyMPPost.

```
OpenMC (tracks.h5) ──► TrackAnalyzer ──► collision file (MCNPX-PoliMi format) ──► pyMPPost ──► light output, pulses
```

---

## ✨ Why

- **pyMPPost** computes realistic detector responses (light output with quenching,
  energy and time resolution, pulse pile-up within a pulse generation time) from the
  collision-by-collision output of **MCNPX-PoliMi**.
- MCNPX-PoliMi is distributed by RSICC and is subject to export control.
  **OpenMC** is open source, but it does not write an MCNPX-PoliMi collision file.
- OpenMC can write the full sequence of particle states of every history (**tracks**).
  `TrackAnalyzer` reconstructs the interactions in the detector from these tracks,
  including the **target nuclide of neutron collisions**, which the track file does not store.
- It works with **any OpenMC version that writes tracks**, including the development branch
  with photonuclear physics
  ([GuySten/openmc, `photonuclear-physics`](https://github.com/GuySten/openmc/tree/photonuclear-physics)),
  where the collision tracking feature of the official release is not available.

---

## ⚙️ How it works

OpenMC stores a state (position, direction, energy, time, weight, cell, material)
after every event of every particle. A **collision** is detected as a change of
direction between two consecutive states inside the detector material.

### 🟣 Neutrons (`NeutronTrackAnalyzer`)
- Only **elastic scattering** is written (PoliMi code `-99`).
- The target nuclide is chosen among a list of **candidates** (e.g. ¹H and C for an organic
  scintillator) using two-body elastic kinematics: the measured E′/E is compared with the
  value predicted at the measured scattering angle, and the candidate with the smallest
  residual is kept if the residual is below a tolerance (default 0.05).
- Not written: collisions below the incident-energy threshold (default 10 eV), thermal
  up-scattering, near-forward scattering, and collisions not compatible with elastic kinematics
  (e.g. inelastic scattering).

### 🟡 Photons (`PhotonTrackAnalyzer`)
- Electrons are not transported (thick-target bremsstrahlung, `electron_treatment = 'ttb'`).
- Direction change with energy loss → **Compton** (1); without energy loss → coherent (2, not
  written by default); track ending with E = 0 → **photoelectric** (3), or **pair production** (5)
  if two 511 keV photons are born at the end point.
- Deposited energy = incident energy − outgoing energy − energy of the secondary photons
  born at the same point (they are tracked separately and would otherwise be counted twice).

Neutrons and photons are always analysed in **separate runs and separate collision files**.

For both particles, collisions that deposit less than a given energy, or that occur outside the
selected cells, can be left out of the collision file, as with the `RPOL` and `IPOL` cards of
MCNPX-PoliMi (see [Same options as MCNPX-PoliMi](#-same-options-as-mcnpx-polimi)).

---

## 📦 Requirements

| What | Notes |
|---|---|
| 🐍 Python ≥ 3.11 | tested with 3.11 and 3.14 |
| `numpy`, `h5py` | used by `TrackAnalyzer` |
| `openmc` (Python API) | the same OpenMC build used for the transport |
| `pandas`, `matplotlib`, `joblib`, `tomli` | used by the notebooks |
| [pyMPPost](https://gitlab.eecs.umich.edu/umich-dnng/pymppost) | not included; install it following its own instructions |
| Nuclear data in OpenMC HDF5 format | e.g. ENDF/B-VII.1 from [openmc.org](https://openmc.org/data/) |

`src/collision_track_check.py` additionally requires the official OpenMC release, version 0.15.3 or later.

```bash
git clone https://github.com/ERMETE-Lab/openmc-to-pymppost.git
cd openmc-to-pymppost
pip install numpy h5py pandas matplotlib joblib tomli
```

---

## 🚀 How to use it

### Step 1 — Run OpenMC with track output 🏃

Any OpenMC model works. A few settings matter for the coupling:

```python
settings.photon_transport = True          # needed for photon cases
settings.electron_treatment = "ttb"       # electrons deposited locally
settings.max_tracks = n_batches * n_particles   # ⚠️ write EVERY history (default is much lower)

openmc.run(tracks=True)                   # writes tracks.h5
```

- ⚠️ `max_tracks` must be at least the number of simulated histories, otherwise the
  response is normalized to the wrong number of source particles.
- ⚠️ Run in **analog mode** (no variance reduction): all particles must have unit weight.
- ⚠️ Avoid **reflective boundaries** touching the detector: a reflection looks like a collision.
- Give the detector material a fixed ID (e.g. `material_id=99`): you will need it in Step 2.

**Official release or photonuclear branch?** The `0_OpenMC_run.ipynb` notebooks have a switch:

```python
PHOTONUCLEAR = True    # True: photonuclear-physics branch; False: official OpenMC release
```

With `True`, the notebook loads a cross-section library that includes photonuclear data and
sets `settings.photonuclear_physics = True` (and stops with an error if the installed OpenMC does
not support it). With `False`, it uses a standard library and no photonuclear settings.

### Step 2 — Write the collision file 📝

```python
import sys
sys.path.append("path/to/OpenMC-to-pyMPPost/src")
from TrackAnalyzer import build_analyzer

# Describe the detector: candidate nuclides with atomic weight ratio (A) and ZAID
detector = {
    "neutron_targets": {"H-1":  {"A": 0.99917, "zaid": 1001},
                        "C-12": {"A": 11.8969, "zaid": 6000}},
    "photon_zaid": 6000,          # single representative ZAID written for photon events
}

analyzer = build_analyzer("neutron", detector,
                          Ethreshold=10.0,       # incident-energy threshold, eV (use 1e3 for photons)
                          rpol=(0.01, 100),      # as the RPOL card of MCNPX-PoliMi, MeV
                          max_residual=0.05,     # neutrons only
                          polimi_file_path="collision_file.d")
analyzer.process_file("tracks.h5",
                      target_materials=[99],     # detector material ID(s)
                      target_cells=[20])         # as IPOL(7)-(8...), OpenMC cell ID(s)
analyzer.summary()
```

For photons: `build_analyzer("photon", detector, Ethreshold=1e3, rpol=(0, 0), polimi_file_path=...)`.

#### 🔁 Same options as MCNPX-PoliMi

When the OpenMC results are compared with MCNPX-PoliMi, use the same values as in the
MCNPX-PoliMi input file, so that the two collision files contain the same kind of events.

| MCNPX-PoliMi card | Meaning | `TrackAnalyzer` | Notebook variable |
|---|---|---|---|
| `RPOL(1)`, `RPOL(2)` | thresholds on the energy **deposited** in a collision by neutrons and photons [MeV]; collisions that deposit less are not written; 0 writes every energy release | `build_analyzer(..., rpol=(RPOL1, RPOL2))` (or `Edep_threshold` in eV in the analyzer constructors) | `RPOL` |
| `IPOL(7)`, `IPOL(8...)` | number and list of the cells for which collisions are written | `process_file(..., target_cells=[...])` with the **OpenMC** cell IDs; `IPOL(7)` is the length of the list; `None` = every cell | `IPOL_CELLS` |
| other `IPOL` entries | source definition, fission multiplicity, photon production | no counterpart: the source is defined in the OpenMC model | — |
| — | threshold on the **incident** energy of the particle [eV] | `Ethreshold` | `E_THRESHOLD_EV` |

- `Ethreshold` and `rpol` are different cuts. `Ethreshold` removes the collisions of particles
  below a given energy: 10 eV for neutrons, below which the target can no longer be considered
  at rest, and the OpenMC photon cutoff (1 keV by default) for photons. `rpol` removes the
  collisions that deposit little energy, as MCNPX-PoliMi does.
- ⚠️ Do not use `Ethreshold` to reproduce `RPOL`: with photons, a 10 keV incident-energy cut
  removes the energy deposited by low-energy photons (e.g. fluorescence X-rays) and distorts
  the photopeak.
- With `rpol=None` (or `(0, 0)`) and `target_cells=None` every reconstructed collision in the
  selected materials is written.
- `summary()` reports how many events were skipped by the `rpol` threshold.

#### 📓 Configuration of `1_write_collision_file.ipynb`

All the settings are in the configuration cell, except the path to `src/`. The values used for
the cases of the paper are given in the last column.

| Variable | Meaning | How to set it (values used in the paper) |
|---|---|---|
| `sys.path.append(...)` (2nd cell) | folder that contains `TrackAnalyzer.py` | path to `src/` of this repository |
| `PARTICLE` | particle written in the collision file | `"neutron"` or `"photon"`, never both: run the notebook twice for both (`"neutron"` for stilbene and ⁴He, `"photon"` for NaI) |
| `DETECTOR` | key of `DETECTOR_LIBRARY` that describes the detector | `"stilbene"`, `"helium4"`, `"nai"` |
| `DETECTOR_MAT_IDS` | OpenMC material ID(s) of the detector; collisions in other materials are not written | the `material_id` given in `0_OpenMC_run.ipynb` (`[99]`) |
| `BASE_DIR` | root folder of this repository on your machine | your path; results are written in `BASE_DIR/results/CASE` |
| `CASE` | case folder, `<detector>/<source>` | e.g. `"trans-stilbene/252Cf"` |
| `TRACK_DIR`, `TRACK_FILE` | folder and name of the OpenMC track file | where `0_OpenMC_run.ipynb` saved `tracks.h5` |
| `E_THRESHOLD_EV` | incident-energy threshold per particle [eV] | `{"neutron": 1e1, "photon": 1e3}`; for photons keep it ≥ the OpenMC photon cutoff |
| `RPOL` | (neutron, photon) thresholds on the deposited energy [MeV], as `RPOL(1)`, `RPOL(2)` | copy the `rpol` card of the MCNPX-PoliMi input: `(0.01, 100)` for stilbene and ⁴He, `(0, 0)` for NaI; `(0, 0)` = write every energy release |
| `IPOL_CELLS` | OpenMC cell ID(s) for which collisions are written, as `IPOL(7)`-`(8...)` | `[20]`; must match `detector_cells` in the pyMPPost `input_config.toml`; `None` = every cell |
| `MAX_RESIDUAL` | neutrons: tolerance on the difference between predicted and measured E′/E | `0.05`; smaller values reject more collisions as not elastic |
| `WRITE_COHERENT` | photons: write coherent scattering events (zero deposited energy) | `False`; they do not change the light output |
| `N_INSPECT` | number of initial histories printed event by event, with the reason why each event is written or skipped | `10`; `0` for a silent run |
| `DETECTOR_LIBRARY` | candidate nuclides (`A` = atomic weight ratio of the nuclear data, ZAID written in the file) and representative photon ZAID of each detector | see *Choosing the candidates* below; add an entry for a new detector |

The collision file is written to `BASE_DIR/results/CASE/OpenMC_collision_file_<PARTICLE>_<DETECTOR>_<N>E5.d`,
where `N` is the number of histories found in the track file. Give this path to pyMPPost in Step 3.

💡 **Choosing the candidates**
- List **every nuclide** of the detector on which collisions can occur: a collision on a
  nuclide that is not listed cannot be identified.
- `A` is the **atomic weight ratio** (nucleus mass / neutron mass) of the nuclear data used
  in the transport, *not* the atomic mass in u (e.g. ⁴He: A = 3.968, not 4.0026). You can read
  it from the data file:
  ```python
  import openmc.data
  openmc.data.IncidentNeutron.from_hdf5("He4.h5").atomic_weight_ratio
  ```
- The ZAID is what pyMPPost uses to pick the light-output function (e.g. `1001` for proton
  recoils, `6000` for carbon recoils): it must match the material card of pyMPPost.

Use `n_inspect=10` in `process_file` to print, for the first histories, every event with the
reason why it was written or skipped. `describe_history("tracks.h5", index)` lists all the
particles of one history.

### Step 3 — Run pyMPPost 💡

Point pyMPPost to the collision file in `input_config.toml`:

```toml
["i/o"]
polimi_det_in = "/path/to/collision_file.d"
output_root   = "./pyMPPost_output_OpenMC/my_case.mpp"

[detectors]
NPS            = 4e5          # number of simulated source particles
detector_cells = [[20]]       # OpenMC cell ID(s) of the detector
detector_mats  = [0]
material_list  = ["material_card.toml"]
```

Then run `pyMPPost input_config.toml`, or use `run_pyMPPost_OpenMC/pyMPPost_analyzer.ipynb`,
which runs pyMPPost, reads the pulses and saves a light-output histogram (`.joblib`)
normalized per source particle.

### Step 4 — Analyse 📊

`notebooks/compare_LO_OpenMC_MCNP.ipynb` loads the histograms of OpenMC and MCNPX-PoliMi and
produces the comparison figures in `images/` (distributions, ratio OpenMC/MCNPX-PoliMi and totals).

---

## 🔁 Reproducing the results of the paper

Each case has its own folder `notebooks/<detector>/<source>/`. Run the notebooks in order:

| # | Notebook | What it does |
|---|---|---|
| 0 | `OpenMC/0_OpenMC_run.ipynb` | builds the model and runs OpenMC with tracks |
| 1 | `OpenMC/1_write_collision_file.ipynb` | `tracks.h5` → collision file |
| 2 | `OpenMC/2_plot_tracks.ipynb` | (optional) plots a few histories |
| 3 | `OpenMC/run_pyMPPost_OpenMC/pyMPPost_analyzer.ipynb` | pyMPPost on the OpenMC collision file → histogram |
| 4 | `MCNPX-PoliMi/pyMPPost_analyzer_MCNP.ipynb` | pyMPPost on the MCNPX-PoliMi collision files → histogram |
| 5 | `notebooks/compare_LO_OpenMC_MCNP.ipynb` | comparison of all cases |

📁 **Paths.** Notebooks and TOML files contain the absolute paths used by the authors
(`/Users/...`, `/Volumes/...`): replace them with your own before running.

💾 **Not included** (too large, regenerable): OpenMC track files (`*.h5`), collision files
(`*.d`, `*.d.neutron`, `*.d.photon`) and pyMPPost pulse files (`*_All_Pulses`, `*_Pulse_Hist`).
The light-output histograms (`results/**/*.joblib`) are included, so the figures can be
reproduced directly with Step 4.

☢️ **MCNPX-PoliMi.** The input files used for the reference simulations are in
`notebooks/<detector>/<source>/MCNPX-PoliMi/*.i`. MCNPX-PoliMi itself must be obtained from
RSICC (package CCC-791).

---

## ✅ Verification

The coupling was verified against MCNPX-PoliMi + pyMPPost on six cases:

| Particle | Detector | Source | Histories |
|---|---|---|---|
| neutrons | trans-stilbene Ø5.08×5.08 cm | ²⁵²Cf (Watt) | 4·10⁵ |
| neutrons | trans-stilbene Ø5.08×5.08 cm | DD (2.45 MeV) | 4·10⁵ |
| neutrons | ⁴He, 190 bar, Ø4.95×20 cm | ²⁵²Cf (Watt) | 4·10⁵ |
| neutrons | ⁴He, 190 bar, Ø4.95×20 cm | DD (2.45 MeV) | 4·10⁵ |
| photons | NaI(Tl) 3″×3″ | ¹³⁷Cs | 10⁵ |
| photons | NaI(Tl) 3″×3″ | ¹⁵²Eu | 10⁵ |

The total number of pulses agrees within 1% in all cases. The neutron reconstruction was also
checked event by event against the collision tracking feature of OpenMC 0.15.3
(`src/collision_track_check.py`, `notebooks/trans-stilbene/252Cf/OpenMC/2_collision_track_check.ipynb`):
99.98% of the written events and 99.82% of the written energy are assigned to the correct nuclide.
See the paper for details.

---

## 📄 Collision file format

Fixed-width columns, as in MCNPX-PoliMi:

| Column | Content | Unit |
|---|---|---|
| 1 | history number | |
| 2 | particle number within the history | |
| 3 | particle type (1 neutron, 2 photon) | |
| 4 | interaction type (−99 elastic; 1, 2, 3, 5 for photons) | |
| 5 | target ZAID | |
| 6 | cell number | |
| 7 | deposited energy | MeV |
| 8 | time | shake |
| 9–11 | position x, y, z | cm |
| 12 | generation number (always 0) | |
| 13 | weight | |
| 14 | number of written events of the particle before this one | |
| 15 | code (always 0) | |
| 16 | incident energy | MeV |

Differences from MCNPX-PoliMi: column 14 counts all written events (MCNPX-PoliMi counts
elastic or Compton scatterings only); column 2 is the index of the particle track in the
OpenMC history; for photons a single representative ZAID is written. None of these
affects the light-output computation.

---

## ⚠️ Limitations

- Only elastic neutron scattering is written. Inelastic collisions are mostly rejected, but
  about 6% of the inelastic collisions on carbon are accepted as elastic scattering on
  hydrogen (relevant only above the carbon inelastic threshold, about 4.8 MeV).
- Photon code 4 (capture with double fluorescence) cannot be separated from code 3.
- Pair production is implemented but has not been specifically verified.
- Photonuclear reactions in the detector material are not identified.
- Electrons are not transported (thick-target bremsstrahlung approximation).
- A direction change at a reflective boundary is indistinguishable from a collision.
- Track files are large: about 1.2 kB per history for the ²⁵²Cf–stilbene case.

---

## 🗂️ Repository structure

```
src/
  TrackAnalyzer.py                tracks.h5 → MCNPX-PoliMi collision file
  collision_track_check.py        event-by-event check against OpenMC collision tracking
notebooks/
  <detector>/<source>/            trans-stilbene/{252Cf, DD}, helium4/{252Cf, DD},
                                  sodium-iodide/{137Cs, 152Eu}
    OpenMC/
      0_OpenMC_run.ipynb          model and OpenMC run (writes tracks.h5)
      1_write_collision_file.ipynb
      2_plot_tracks.ipynb
      geometry.xml, materials.xml, settings.xml
      run_pyMPPost_OpenMC/        pyMPPost input files and pyMPPost_analyzer.ipynb
    MCNPX-PoliMi/
      *.i                         MCNPX-PoliMi input file
      input_config.toml, material_card.toml
      pyMPPost_analyzer_MCNP.ipynb
  sodium-iodide/152Eu/OpenMC/source.py      ¹⁵²Eu source (DDEP data)
  trans-stilbene/252Cf/OpenMC/2_collision_track_check.ipynb
  compare_LO_OpenMC_MCNP.ipynb    OpenMC vs MCNPX-PoliMi comparison and figures
data/birks/stilbene_dEdx.txt      proton stopping power table for the stilbene material card
results/<detector>/<source>/      light-output histograms (*.joblib)
images/                           comparison figures (PDF)
```

---

## 📚 Citation

If you use this code, please cite: <!-- TODO: paper reference (arXiv / journal) and Zenodo DOI -->

## 📜 License

This project is released under the [MIT License](LICENSE).
OpenMC and pyMPPost are not distributed with this repository and are covered by their own licenses.

## 🙏 Acknowledgements

pyMPPost and MPPost are developed by the Detection for Nuclear Nonproliferation Group (DNNG),
University of Michigan. MCNPX-PoliMi is developed by Politecnico di Milano and the University of
Michigan and distributed by RSICC.
The proton stopping power table for trans-stilbene (`data/birks/stilbene_dEdx.txt`) was provided
by the DNNG, University of Michigan.
<!-- TODO: funding -->
