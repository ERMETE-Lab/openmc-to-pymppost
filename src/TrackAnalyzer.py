"""
TrackAnalyzer.py
================

Convert OpenMC particle tracks (``tracks.h5``) into an MCNPX-PoliMi style
collision file that can be used as input by pyMPPost.

pyMPPost can be downloaded from:
    https://gitlab.eecs.umich.edu/umich-dnng/pymppost

Workflow
--------
1. Run OpenMC with track output enabled (e.g. ``openmc.run(tracks=True)``
   together with ``settings.max_tracks``).
2. Build ONE analyzer for the particle of interest. Neutrons and photons are
   always analysed separately, never in the same collision file.
3. Call ``process_file()`` on the tracks file to write the collision file.

Interactions written to the collision file
------------------------------------------
Neutrons (particle type 1):
    -99  elastic scattering. The target nuclide is identified from two-body
         elastic kinematics (energy ratio vs. scattering angle).

Photons (particle type 2):
     1   Compton scattering: the photon survives with lower energy and a new
         direction.
     2   coherent scattering: the photon survives with the same energy and a
         new direction (zero deposited energy, optional).
     3   capture with fluorescence (photoelectric): the photon track ends with
         E = 0 in its last state.
     5   pair production: the photon track ends with E = 0 and two 511 keV
         annihilation photons are born at the end point (not yet validated).
    Code 4 (capture with double fluorescence) cannot be separated from code 3
    using track data alone, so it is never written.

Implementation notes
--------------------
* OpenMC stores every track state AFTER the corresponding event. A collision
  detected between states i and i+1 therefore happens at state i+1: position,
  time, cell and material are read at i+1, while the incident energy and the
  weight are read at i.
* Track states do not store the reaction type (MT). Interactions are
  reconstructed from energy/direction changes and, for photons, from the
  secondary particles born at the interaction site.
* A direction change at a reflective boundary is indistinguishable from a
  collision. Avoid reflective surfaces touching the detector material.

Correspondence with the MCNPX-PoliMi input cards
------------------------------------------------
* RPOL(1), RPOL(2): thresholds on the energy released in a collision, for
  neutrons and photons [MeV]. Collisions that deposit less are not written.
  Here: ``build_analyzer(..., rpol=(RPOL1, RPOL2))``, or ``Edep_threshold``
  [eV] in the analyzer constructors. 0 writes every energy release.
* IPOL(7), IPOL(8...): number and list of the cells for which collisions are
  written. Here: ``process_file(..., target_cells=[...])`` with the OpenMC
  cell IDs (IPOL(7) is simply the length of the list). None = every cell.
* The other IPOL entries (source definition, fission multiplicity, photon
  production options) have no counterpart: the source is defined in OpenMC.
"""

import numpy as np
import h5py
import openmc

# OpenMC particle type codes: enum values (OpenMC <= 0.15.x) and PDG numbers
# (OpenMC >= 0.16, where openmc.ParticleType wraps the PDG code)
PARTICLE_NAMES = {0: "neutron", 1: "photon", 2: "electron", 3: "positron",
                  2112: "neutron", 22: "photon", 11: "electron", -11: "positron"}

# Electron rest mass energy [eV]
ELECTRON_MASS_EV = 0.51099895e6

# Photon interaction codes used in the collision file
PHOTON_COMPTON = 1
PHOTON_COHERENT = 2
PHOTON_CAPTURE = 3
PHOTON_CAPTURE_DOUBLE = 4   # not reconstructable from tracks, kept for reference
PHOTON_PAIR = 5

# Neutron interaction codes used in the collision file
NEUTRON_ELASTIC = -99


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
def particle_name(p):
    """Return the lowercase name of an OpenMC particle type (enum, int or str)."""
    if hasattr(p, "name"):
        return p.name.lower()
    if isinstance(p, (int, np.integer)):
        return PARTICLE_NAMES.get(int(p), "unknown")
    return str(p).lower()


def polimi_line(history, particle, ptype, itype, zaid, cell, edep_MeV, time_shake,
                x, y, z, wgt, n_scat, einc_MeV, gen=0, code=0):
    """One line of an MCNPX-PoliMi collision file (fixed-width columns)."""
    return (f"{history:>10}{particle:>5}{ptype:>3}{itype:>5}{zaid:>6}{cell:>11}"
            f"{edep_MeV:>11.6f}{time_shake:>17.2f}{x:>9.2f}{y:>9.2f}{z:>9.2f}"
            f"{gen:>5}{wgt:>11.3E}{n_scat:>6}{code:>5}{einc_MeV:>11.3E}\n")


def _track_sort_key(dset_name):
    """Sort 'track_<batch>_<generation>_<particle>' datasets like OpenMC does."""
    _, batch, gen, particle = dset_name.split("_")
    return int(batch), int(gen), int(particle)


def _vec(field, i):
    """Return entry i of a structured (x, y, z) field as a plain numpy vector."""
    try:
        return np.array([field[i]["x"], field[i]["y"], field[i]["z"]], dtype=float)
    except (ValueError, IndexError, KeyError, TypeError):
        return np.asarray(field[i], dtype=float)


def describe_history(track_file, index=0):
    """
    Print every particle track of one history (type, number of states, birth
    energy and position). Useful to check which secondaries OpenMC writes.
    """
    with h5py.File(track_file, "r") as fh:
        names = sorted((k for k in fh.keys() if k.startswith("track_")),
                       key=_track_sort_key)
        track = openmc.Track(fh[names[index]])
        print(f"History #{index} ({names[index]}): "
              f"{len(track.particle_tracks)} particle tracks")
        for k, pt in enumerate(track.particle_tracks):
            s = pt.states.flatten()
            r0 = _vec(s["r"], 0)
            print(f"  P #{k:<3} {particle_name(pt.particle):<9} states={len(s):<4} "
                  f"E0={s['E'][0]:.4e} eV  r0=({r0[0]:.4f}, {r0[1]:.4f}, {r0[2]:.4f})  "
                  f"E_end={s['E'][-1]:.4e} eV")


# ---------------------------------------------------------------------------
# Base class
# ---------------------------------------------------------------------------
class TrackAnalyzer:
    """
    Common machinery: HDF5 reading, particle/material filtering, collision
    file writing and bookkeeping. Subclasses only implement ``_find_events``
    (and optionally ``_history_context``).
    """

    PARTICLE = None       # "neutron" or "photon", set by subclasses
    POLIMI_PTYPE = None   # particle type code in the collision file

    def __init__(self, Ethreshold, Edep_threshold=0.0, polimi_file_path=None,
                 verbose=False, dir_atol=1e-6):
        self.Ethreshold = Ethreshold      # events with incident energy below this [eV] are skipped
        self.Edep_threshold = Edep_threshold  # events depositing less than this [eV] are skipped (RPOL)
        self.verbose = verbose
        self.dir_atol = dir_atol          # tolerance used to detect a direction change

        # Bookkeeping
        self.counts = {}                  # written events per label (nuclide or interaction)
        self.n_tracks = 0                 # analysed particle tracks
        self.n_tracks_no_events = 0       # tracks without any written event
        self.n_below_edep = 0             # events skipped by Edep_threshold (RPOL)

        self.polimi_file = open(polimi_file_path, "w") if polimi_file_path else None

    # --- file handling ------------------------------------------------------
    def close(self):
        """Flush and close the collision file."""
        if self.polimi_file:
            self.polimi_file.flush()
            self.polimi_file.close()
            self.polimi_file = None

    def __del__(self):
        self.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    # --- hooks for subclasses -------------------------------------------------
    def _history_context(self, track):
        """Per-history data shared by all tracks of a history (None by default)."""
        return None

    def _find_events(self, states, context=None):
        """Yield event dictionaries (see ``_new_event``) in chronological order."""
        raise NotImplementedError

    @staticmethod
    def _new_event(j, E_in, E_out, wgt, itype):
        """Create an event record. ``j`` is the state index where the event happens."""
        return {
            "j": j,                            # state index of the interaction
            "E_in": float(E_in),               # incident energy [eV]
            "E_out": float(E_out),             # outgoing energy of the same particle [eV]
            "edep": max(float(E_in - E_out), 0.0),  # deposited energy [eV]
            "wgt": float(wgt),                 # particle weight before the event
            "itype": itype,                    # interaction code in the collision file
            "zaid": 0,                         # target ZAID in the collision file
            "label": None,                     # name used for counting / printing
            "error": None,                     # reason why the event is not written
            "info": "",                        # extra text for verbose output
        }

    # --- main loop ------------------------------------------------------------
    def process_file(self, track_file, target_materials=None, target_cells=None,
                     n_inspect=0, max_histories=None):
        """
        Analyse every history in ``track_file`` and write the collision file.

        target_materials : list of OpenMC material IDs of the detector(s);
                           None keeps events in every material.
        target_cells     : list of OpenMC cell IDs for which collisions are
                           written (as IPOL(7)-(8...) in MCNPX-PoliMi);
                           None keeps events in every cell.
        n_inspect        : number of initial histories printed in verbose mode.
        max_histories    : optional limit on the number of histories (debug).
        Returns the number of processed histories.
        """
        with h5py.File(track_file, "r") as fh:
            names = sorted((k for k in fh.keys() if k.startswith("track_")),
                           key=_track_sort_key)
            if max_histories is not None:
                names = names[:max_histories]

            for h, name in enumerate(names):
                # Lazy load: only this history is kept in memory
                track = openmc.Track(fh[name])
                self.verbose = h < n_inspect
                if self.verbose:
                    print(f"\n========== HISTORY #{h} ==========")

                context = self._history_context(track)
                for p_idx, pt in enumerate(track.particle_tracks):
                    if particle_name(pt.particle) != self.PARTICLE:
                        continue
                    self.analyze_track(pt.states.flatten(), history_id=h,
                                       particle_id=p_idx,
                                       target_materials=target_materials,
                                       target_cells=target_cells,
                                       context=context)
        self.close()
        return len(names)

    def analyze_track(self, states, history_id=0, particle_id=0,
                      target_materials=None, target_cells=None, context=None,
                      return_results=False):
        """
        Analyse a single particle track and write its valid events.
        Event dictionaries are returned only if ``return_results=True``.
        """
        tag = f"[H #{history_id} | P #{particle_id}] "
        results = [] if return_results else None
        n_written = 0   # events already written for this particle (PoliMi scatter number)

        for ev in self._find_events(states, context):
            j = ev["j"]

            # 1. Material filter, evaluated where the interaction happens
            if target_materials is not None and states["material_id"][j] not in target_materials:
                continue
            # 2. Cell filter (IPOL(7)-(8...) in MCNPX-PoliMi)
            if target_cells is not None and states["cell_id"][j] not in target_cells:
                continue

            # 3. Energy threshold on the incident energy
            if ev["error"] is None and ev["E_in"] < self.Ethreshold:
                ev["error"] = (f"Incident energy below threshold "
                               f"({ev['E_in']:.2e} < {self.Ethreshold:.2e} eV).")

            # 4. Threshold on the deposited energy (RPOL in MCNPX-PoliMi)
            if ev["error"] is None and ev["edep"] < self.Edep_threshold:
                ev["error"] = (f"Deposited energy below threshold "
                               f"({ev['edep']:.2e} < {self.Edep_threshold:.2e} eV, RPOL).")
                self.n_below_edep += 1

            # 5. Write and count
            if ev["error"] is None:
                self._write_polimi_event(ev, states, history_id, particle_id, n_written)
                n_written += 1
                self.counts[ev["label"]] = self.counts.get(ev["label"], 0) + 1

            if self.verbose:
                self._print_event(ev, tag)
            if return_results:
                results.append(ev)

        self.n_tracks += 1
        if n_written == 0:
            self.n_tracks_no_events += 1
            if self.verbose:
                where = f" in materials {target_materials}" if target_materials is not None else ""
                where += f" in cells {target_cells}" if target_cells is not None else ""
                print(f"  → {tag}No valid events{where}.")
        if self.verbose:
            print(f"  → {tag}End of track: {len(states)} states, "
                  f"final energy {states['E'][-1]:.2e} eV")
        return results

    # --- output -------------------------------------------------------------
    def _write_polimi_event(self, ev, states, history_id, particle_id, n_scat):
        """Write one event using the fixed-width MCNPX-PoliMi column layout."""
        if self.polimi_file is None:
            return

        j = ev["j"]                          # interaction state (see module notes)
        x, y, z = _vec(states["r"], j)
        cell = int(states["cell_id"][j])
        time = states["time"][j] * 1e8       # s -> shakes
        edep = ev["edep"] / 1e6              # eV -> MeV
        einc = ev["E_in"] / 1e6              # eV -> MeV

        # generation and PoliMi code columns are not tracked (written as 0)
        self.polimi_file.write(polimi_line(
            history_id + 1, particle_id + 1, self.POLIMI_PTYPE, ev["itype"],
            ev["zaid"], cell, edep, time, x, y, z, ev["wgt"], n_scat, einc))

    def _print_event(self, ev, tag):
        """Verbose, human-readable description of one event."""
        status = f"✓ {ev['label']}" if ev["error"] is None else "✗ SKIPPED"
        print(f"  {tag}Event at state {ev['j']} (code {ev['itype']}): {status}")
        print(f"    Energy: {ev['E_in']:.3e} eV → {ev['E_out']:.3e} eV | "
              f"Edep = {ev['edep']:.3e} eV")
        if ev["info"]:
            print(ev["info"])
        if ev["error"]:
            print(f"    Reason: {ev['error']}")
        print("  " + "-" * 48)

    def summary(self):
        """Print a short summary of the analysis."""
        print(f"Particle            : {self.PARTICLE}")
        print(f"Tracks analysed     : {self.n_tracks}")
        print(f"Tracks w/o events   : {self.n_tracks_no_events}")
        print(f"Edep threshold      : {self.Edep_threshold / 1e6:g} MeV (RPOL), "
              f"{self.n_below_edep} events skipped")
        print(f"Written events      : {sum(self.counts.values())}")
        for label, n in sorted(self.counts.items(), key=lambda x: -x[1]):
            print(f"  {label:<18}: {n}")


# ---------------------------------------------------------------------------
# Neutrons
# ---------------------------------------------------------------------------
class NeutronTrackAnalyzer(TrackAnalyzer):
    """
    Elastic neutron scattering. The target nuclide is chosen among
    ``targets`` by comparing the measured E'/E with the two-body elastic
    kinematics prediction at the measured scattering angle.

    targets : {name: {"A": atomic weight ratio, "zaid": ZAID in the output}}
    """

    PARTICLE = "neutron"
    POLIMI_PTYPE = 1

    def __init__(self, targets, Ethreshold=1e1, max_residual=0.05, **kwargs):
        super().__init__(Ethreshold=Ethreshold, **kwargs)
        self.targets = targets
        self.max_residual = max_residual   # max |predicted - measured| on E'/E

    def _find_events(self, states, context=None):
        u = states["u"]
        for i in range(len(states["E"]) - 1):
            u0, u1 = _vec(u, i), _vec(u, i + 1)
            if np.allclose(u0, u1, atol=self.dir_atol):
                continue   # no direction change: surface crossing, not a collision
            yield self._identify_nuclide(states, i, u0, u1)

    def _identify_nuclide(self, states, i, u0, u1):
        """Identify the target nuclide for the collision between states i and i+1."""
        j = i + 1
        E_in, E_out = states["E"][i], states["E"][j]
        ev = self._new_event(j, E_in, E_out, states["wgt"][i], NEUTRON_ELASTIC)

        cos_theta = float(np.clip(np.dot(u0, u1), -1.0, 1.0))
        sin2_theta = 1.0 - cos_theta**2
        ratio = E_out / E_in

        if ratio > 1.0 + 1e-5:
            ev["error"] = "Thermal up-scattering (E' > E): kinematics not applicable."
            return ev
        if ratio > 1.0 - 1e-4:
            ev["error"] = "Near-forward scatter: energy loss too small to discriminate."
            return ev

        residuals = {}
        for name, info in self.targets.items():
            A = info["A"]
            # Hydrogen-like targets (A <= 1) cannot back-scatter a neutron
            if A < 1.01 and cos_theta < -1e-9:
                continue
            # Minimum reachable E'/E for this target
            alpha = ((A - 1) / (A + 1))**2
            if ratio < alpha - 1e-5:
                continue
            inner = A**2 - sin2_theta
            if inner < 0:
                continue
            predicted = ((cos_theta + np.sqrt(inner)) / (A + 1))**2
            residuals[name] = abs(predicted - ratio)

        if not residuals:
            ev["error"] = "No target compatible with elastic kinematics (likely inelastic)."
            return ev

        ranked = sorted(residuals.items(), key=lambda x: x[1])
        best, best_res = ranked[0]
        ev["info"] = (f"    cos θ = {cos_theta:.3f} | E'/E = {ratio:.3f}\n" +
                      "\n".join(f"    {n:<8} residual {r:.2e}{'  ✅' if n == best else ''}"
                                for n, r in ranked))

        if best_res > self.max_residual:
            ev["error"] = (f"Best match {best} exceeds residual tolerance "
                           f"({best_res:.2e} > {self.max_residual}); likely inelastic.")
            return ev

        ev["label"] = best
        ev["zaid"] = self.targets[best]["zaid"]
        ev["residual"] = best_res
        return ev


# Backward compatibility with the original class name
KinematicTrackAnalyzer = NeutronTrackAnalyzer


# ---------------------------------------------------------------------------
# Photons
# ---------------------------------------------------------------------------
class PhotonTrackAnalyzer(TrackAnalyzer):
    """
    Photon interactions reconstructed from tracks:

    * direction change, same energy                -> coherent scattering (2)
    * direction change, lower energy               -> Compton scattering  (1)
    * track ends with E = 0 and two 511 keV photons
      born at the end point                        -> pair production     (5)
    * track ends with E = 0 otherwise              -> photoelectric capture (3)
    * track ends with E > 0                        -> leakage / energy cutoff,
                                                      nothing written

    With the thick-target bremsstrahlung approximation (TTB, OpenMC default)
    electrons and positrons are not transported and do not appear in the
    tracks: their energy is deposited locally and only the bremsstrahlung,
    fluorescence and annihilation photons are tracked as new particles.

    Deposited energy at a site = incident energy - outgoing photon energy
    - energy of every secondary photon born at that site (fluorescence,
    bremsstrahlung, annihilation photons). Those secondary photons are
    tracked separately and may produce their own events in the detector,
    so they must not be counted twice.

    target_zaid : ZAID written in the nucleus column. The target of a photon
                  interaction cannot be identified from tracks, so a single
                  representative value is used for the detector material.

    Note: pair production (code 5) has not been validated yet against a
    source above 1.022 MeV.
    """

    PARTICLE = "photon"
    POLIMI_PTYPE = 2

    LABELS = {PHOTON_COMPTON: "Compton", PHOTON_COHERENT: "Coherent",
              PHOTON_CAPTURE: "Photoelectric", PHOTON_PAIR: "Pair production"}

    def __init__(self, target_zaid, Ethreshold=1e3, write_coherent=False,
                 site_atol=1e-6, **kwargs):
        super().__init__(Ethreshold=Ethreshold, **kwargs)
        self.target_zaid = target_zaid
        self.write_coherent = write_coherent   # coherent events deposit no energy
        self.site_atol = site_atol             # [cm] tolerance to match birth positions

    def _history_context(self, track):
        """Birth type, position and energy of every particle in the history."""
        types, pos, energy = [], [], []
        for pt in track.particle_tracks:
            s = pt.states.flatten()
            if len(s) == 0:
                continue
            types.append(particle_name(pt.particle))
            pos.append(_vec(s["r"], 0))
            energy.append(float(s["E"][0]))
        return {"types": np.array(types, dtype=object),
                "pos": np.array(pos, dtype=float).reshape(-1, 3),
                "E": np.array(energy, dtype=float)}

    def _secondary_photons_at(self, context, site):
        """Birth energies [eV] of the photons born at ``site``."""
        if context is None or len(context["E"]) == 0:
            return np.zeros(0)
        born = np.all(np.abs(context["pos"] - site) <= self.site_atol, axis=1)
        return context["E"][born & (context["types"] == "photon")]

    def _find_events(self, states, context=None):
        E, u, r = states["E"], states["u"], states["r"]
        n = len(E)

        # --- how does the track end? -----------------------------------------
        # OpenMC sets E = 0 in the last state when a photon is absorbed
        # (photoelectric or pair production); leaking photons keep their energy.
        absorbed = n >= 2 and E[n - 1] == 0.0 and E[n - 2] > 0.0

        # --- scattering events along the track -------------------------------
        # When the photon is absorbed, the last transition is the absorption
        # itself and is handled separately below.
        n_transitions = n - 2 if absorbed else n - 1
        for i in range(n_transitions):
            u0, u1 = _vec(u, i), _vec(u, i + 1)
            if np.allclose(u0, u1, atol=self.dir_atol):
                continue   # surface crossing
            j = i + 1

            if np.isclose(E[j], E[i], rtol=1e-12, atol=0.0):
                ev = self._new_event(j, E[i], E[j], states["wgt"][i], PHOTON_COHERENT)
                ev["edep"] = 0.0
                if not self.write_coherent:
                    ev["error"] = "Coherent scattering (no energy deposited, not written)."
            elif E[j] > E[i]:
                ev = self._new_event(j, E[i], E[j], states["wgt"][i], PHOTON_COMPTON)
                ev["error"] = "Photon energy increased: not a physical photon collision."
            else:
                ev = self._new_event(j, E[i], E[j], states["wgt"][i], PHOTON_COMPTON)
                e_sec = self._secondary_photons_at(context, _vec(r, j)).sum()
                ev["edep"] = max(E[i] - E[j] - e_sec, 0.0)
                ev["info"] = f"    Secondary photons born at site: {e_sec:.3e} eV"

            ev["zaid"] = self.target_zaid
            ev["label"] = self.LABELS[ev["itype"]]
            yield ev

        # --- absorption at the end of the track ------------------------------
        if absorbed:
            j = n - 1
            E_in = E[n - 2]                    # energy during the last flight
            sec = self._secondary_photons_at(context, _vec(r, j))

            # Pair production: only possible above 2 m_e c^2; with TTB the
            # positron annihilates on the spot, so two 511 keV photons are
            # born at the end point
            n_511 = np.sum(np.isclose(sec, ELECTRON_MASS_EV, rtol=1e-3))
            is_pair = E_in > 2.0 * ELECTRON_MASS_EV and n_511 >= 2
            itype = PHOTON_PAIR if is_pair else PHOTON_CAPTURE

            ev = self._new_event(j, E_in, 0.0, states["wgt"][n - 2], itype)
            edep = E_in - sec.sum()
            if is_pair:
                # The rest mass of the pair never becomes deposited energy
                edep = min(edep, E_in - 2.0 * ELECTRON_MASS_EV)
            ev["edep"] = max(edep, 0.0)
            ev["zaid"] = self.target_zaid
            ev["label"] = self.LABELS[itype]
            ev["info"] = f"    Secondary photons born at site: {sec.sum():.3e} eV"
            yield ev


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------
def build_analyzer(particle, detector, rpol=None, **kwargs):
    """
    Create the right analyzer from a detector description, e.g.

        detector = {"neutron_targets": {"H-1": {"A": 0.99917, "zaid": 1001}, ...},
                    "photon_zaid": 6000}

    rpol : (neutron, photon) thresholds on the deposited energy [MeV], as the
           RPOL card of MCNPX-PoliMi; the one of ``particle`` is used.
           None = no threshold (or ``Edep_threshold`` [eV] if given).
    Extra keyword arguments are forwarded to the analyzer constructor.
    """
    particle = particle.lower()
    if rpol is not None:
        if "Edep_threshold" in kwargs:
            raise ValueError("Give either rpol [MeV] or Edep_threshold [eV], not both.")
        kwargs["Edep_threshold"] = 1e6 * rpol[0 if particle == "neutron" else 1]
    if particle == "neutron":
        return NeutronTrackAnalyzer(targets=detector["neutron_targets"], **kwargs)
    if particle == "photon":
        return PhotonTrackAnalyzer(target_zaid=detector["photon_zaid"], **kwargs)
    raise ValueError(f"Unsupported particle '{particle}': use 'neutron' or 'photon'.")