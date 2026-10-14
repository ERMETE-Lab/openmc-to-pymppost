"""
collision_track_check.py
========================

Event-by-event verification of ``NeutronTrackAnalyzer`` against the
``collision_track`` feature of OpenMC (Saliba et al., Nucl. Sci. Eng., 2026,
doi:10.1080/00295639.2026.2675195; available from OpenMC 0.15.3).

``collision_track`` records, for every collision, the true reaction (MT), the
target nuclide and the energy lost by the particle. Running ONE OpenMC
simulation with both ``tracks=True`` and ``settings.collision_track`` enabled
gives two descriptions of the same collisions:

* truth          : collision_track.h5 (MT, nuclide, dE)
* reconstruction : tracks.h5 -> NeutronTrackAnalyzer (kinematic identification)

The two are matched by history (source particle id) and collision position.

The module also provides ``write_polimi_from_collision_track``, which writes an
MCNPX-PoliMi collision file directly from collision_track data, so that
pyMPPost can be run on the "true" elastic events and compared with the
reconstructed ones at the light-output level.

Conventions (checked against the OpenMC source, collision_track.cpp):
* E   = energy before the collision [eV];  dE = E_before - E_after [eV]
* r   = collision position [cm];           time [s]
* parent_id = id of the source particle (same number as in 'track_<b>_<g>_<id>')
* nuclide_id = PDG number (OpenMC >= 0.16) or ZAM = 10000 Z + 10 A + m (0.15.x)
* particle   = PDG number (OpenMC >= 0.16) or enum value (0.15.x)
"""

import h5py
import numpy as np
import pandas as pd
import openmc

from TrackAnalyzer import _track_sort_key, _vec, particle_name, polimi_line

try:
    from openmc.data import REACTION_NAME
except ImportError:  # pragma: no cover
    REACTION_NAME = {}

NEUTRON_CODES = {0, 2112}   # enum (<= 0.15.x) and PDG (>= 0.16)

# Short categories for the rejection messages of NeutronTrackAnalyzer
REJECTION_REASONS = [
    ("Incident energy below threshold", "below threshold"),
    ("Thermal up-scattering", "up-scattering"),
    ("Near-forward scatter", "near-forward"),
    ("No target compatible", "no compatible target"),
    ("Best match", "residual > tolerance"),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def decode_nuclide(code):
    """Return (Z, A) from a PDG nuclear code or a ZAM number (A = 0: natural element)."""
    code = int(code)
    if code >= 1_000_000_000:                 # PDG: 10LZZZAAAI
        return (code // 10_000) % 1000, (code // 10) % 1000
    if code > 0:                              # ZAM: 10000 Z + 10 A + m
        return code // 10_000, (code // 10) % 1000
    return 0, 0


def nuclide_label(Z, A):
    if Z == 0:
        return "none"
    symbol = openmc.data.ATOMIC_SYMBOL.get(int(Z), f"Z{Z}")
    return f"{symbol}-nat" if A == 0 else f"{symbol}-{A}"


def reaction_label(mt):
    return REACTION_NAME.get(int(mt), f"MT{int(mt)}")


def rejection_category(message):
    if message is None:
        return None
    for prefix, short in REJECTION_REASONS:
        if message.startswith(prefix):
            return short
    return "other"


# ---------------------------------------------------------------------------
# Truth: collision_track file
# ---------------------------------------------------------------------------
def load_truth(collision_track_file, material_ids=None):
    """Neutron collisions recorded by collision_track, one row per collision."""
    data = openmc.read_collision_track_file(collision_track_file)

    particle = np.asarray(data["particle"]).astype(int)
    keep = np.isin(particle, list(NEUTRON_CODES))
    if material_ids is not None:
        keep &= np.isin(data["material_id"], material_ids)
    d = data[keep]

    ZA = np.array([decode_nuclide(c) for c in d["nuclide_id"]]).reshape(-1, 2)
    truth = pd.DataFrame({
        "history": d["parent_id"].astype(np.int64),
        "x": d["r"]["x"], "y": d["r"]["y"], "z": d["r"]["z"],
        "time": d["time"],
        "E_in": d["E"], "dE": d["dE"],
        "mt": d["event_mt"].astype(int),
        "Z": ZA[:, 0], "A": ZA[:, 1],
        "cell": d["cell_id"], "wgt": d["wgt"],
    })
    truth["nuclide"] = [nuclide_label(z, a) for z, a in zip(truth["Z"], truth["A"])]
    truth["reaction"] = [reaction_label(mt) for mt in truth["mt"]]
    truth["true_label"] = truth["nuclide"] + " " + truth["reaction"]
    return truth


# ---------------------------------------------------------------------------
# Reconstruction: tracks file + NeutronTrackAnalyzer
# ---------------------------------------------------------------------------
def load_reconstructed(track_file, analyzer, material_ids=None, max_histories=None):
    """
    Every collision found by the analyzer (written or rejected), one row each.
    ``analyzer`` must be a NeutronTrackAnalyzer built WITHOUT an output file.
    """
    rows = []
    with h5py.File(track_file, "r") as fh:
        names = sorted((k for k in fh if k.startswith("track_")), key=_track_sort_key)
        if max_histories is not None:
            names = names[:max_histories]

        for name in names:
            history = _track_sort_key(name)[2]          # source particle id
            track = openmc.Track(fh[name])
            for p_idx, pt in enumerate(track.particle_tracks):
                if particle_name(pt.particle) != "neutron":
                    continue
                states = pt.states.flatten()
                events = analyzer.analyze_track(
                    states, history_id=history, particle_id=p_idx,
                    target_materials=material_ids, return_results=True)
                for ev in events:
                    x, y, z = _vec(states["r"], ev["j"])
                    rows.append({
                        "history": history, "particle_idx": p_idx,
                        "x": x, "y": y, "z": z,
                        "time_reco": float(states["time"][ev["j"]]),
                        "E_in_reco": ev["E_in"], "edep_reco": ev["edep"],
                        "written": ev["error"] is None,
                        "reco_label": ev["label"] if ev["error"] is None else None,
                        "reject_reason": rejection_category(ev["error"]),
                        "residual": ev.get("residual", np.nan),
                    })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Matching and summary
# ---------------------------------------------------------------------------
def match_events(truth, reco, decimals=7):
    """
    Outer join of truth and reconstruction on (history, position).
    Positions are rounded to ``decimals`` cm (default 1 nm): both files store the
    same double, rounding only protects against I/O round-off.
    """
    def keyed(df):
        df = df.copy()
        for c in ("x", "y", "z"):
            df[f"k{c}"] = df[c].round(decimals)
        return df

    t, r = keyed(truth), keyed(reco)
    keys = ["history", "kx", "ky", "kz"]
    for label, df in (("truth", t), ("reconstruction", r)):
        n_dup = df.duplicated(keys).sum()
        if n_dup:
            print(f"Warning: {n_dup} duplicated keys in {label}; increase `decimals`.")

    m = t.merge(r.drop(columns=["x", "y", "z"]), on=keys, how="outer", indicator=True)
    m["status"] = m["_merge"].map({"both": "matched", "left_only": "truth only",
                                   "right_only": "reconstruction only"})

    def outcome(row):
        if row["status"] == "truth only":
            return "not seen (no direction change)"
        if row["status"] == "reconstruction only":
            return "spurious"
        if row["written"]:
            return f"written as {row['reco_label']}"
        return f"rejected: {row['reject_reason']}"

    m["outcome"] = m.apply(outcome, axis=1)
    m["true_label"] = m["true_label"].fillna("no true collision")
    return m.drop(columns=["_merge"])


def summarize(m, energy_unit=1e6):
    """
    Return three tables:
      counts : number of collisions, true reaction (rows) vs outcome (columns)
      energy : true energy lost by the neutron [MeV], same layout
      edep   : agreement of deposited energy for correctly identified events
    """
    counts = pd.crosstab(m["true_label"], m["outcome"], margins=True)
    energy = pd.crosstab(m["true_label"], m["outcome"], values=m["dE"] / energy_unit,
                         aggfunc="sum", margins=True).fillna(0.0)

    ok = m[(m["status"] == "matched") & m["written"].fillna(False).astype(bool)]
    diff = (ok["edep_reco"] - ok["dE"]).abs()
    edep = pd.Series({
        "matched written events": len(ok),
        "max |edep_reco - dE| [eV]": diff.max() if len(ok) else np.nan,
        "max relative difference": (diff / ok["dE"].where(ok["dE"] > 0)).max() if len(ok) else np.nan,
        "max |t_reco - t_true| [s]": (ok["time_reco"] - ok["time"]).abs().max() if len(ok) else np.nan,
    })
    return counts, energy, edep


def identification_check(m, expected):
    """
    Fraction of true elastic collisions on each nuclide that are written with
    the right label. ``expected`` maps the collision_track nuclide label to the
    analyzer label, e.g. {"H-1": "H-1", "C-nat": "C-12"}.
    """
    el = m[(m["mt"] == 2) & (m["status"] == "matched")]
    rows = []
    for true_nuc, reco_label in expected.items():
        sub = el[el["nuclide"] == true_nuc]
        n = len(sub)
        right = (sub["reco_label"] == reco_label).sum()
        wrong = (sub["written"].astype(bool) & (sub["reco_label"] != reco_label)).sum()
        rows.append({"nuclide": true_nuc, "true elastic (seen)": n,
                     "correct": right, "wrong nuclide": wrong,
                     "rejected": n - right - wrong,
                     "correct [%]": 100 * right / n if n else np.nan,
                     "dE rejected / dE total [%]":
                         100 * sub.loc[~sub["written"].astype(bool), "dE"].sum() / sub["dE"].sum()
                         if n else np.nan})
    return pd.DataFrame(rows).set_index("nuclide")


# ---------------------------------------------------------------------------
# collision_track -> MCNPX-PoliMi collision file
# ---------------------------------------------------------------------------
def write_polimi_from_collision_track(collision_track_file, output_path,
                                      material_ids=None, zaid_map=None,
                                      mts=(2,), Ethreshold=10.0):
    """
    Write an MCNPX-PoliMi collision file from collision_track data.

    Only reactions in ``mts`` are written (default: elastic, PoliMi code -99),
    i.e. the same event class written by NeutronTrackAnalyzer. Writing other
    reactions would require their PoliMi interaction codes.

    zaid_map : {Z: ZAID written in the file}, default {1: 1001, 6: 6000}
    Returns the number of written lines.
    """
    polimi_code = {2: -99}
    missing = set(mts) - set(polimi_code)
    if missing:
        raise ValueError(f"No PoliMi code defined for MT {sorted(missing)}.")
    zaid_map = zaid_map or {1: 1001, 6: 6000}

    t = load_truth(collision_track_file, material_ids)
    t = t[t["mt"].isin(mts) & (t["E_in"] >= Ethreshold)]
    t = t.sort_values(["history", "time"], kind="stable")   # bank order is not chronological

    n = 0
    with open(output_path, "w") as f:
        for history, grp in t.groupby("history", sort=True):
            for n_scat, row in enumerate(grp.itertuples(index=False)):
                f.write(polimi_line(
                    int(history), 1, 1, polimi_code[row.mt], zaid_map.get(int(row.Z), 0),
                    int(row.cell), row.dE / 1e6, row.time * 1e8,
                    row.x, row.y, row.z, row.wgt, n_scat, row.E_in / 1e6))
                n += 1
    return n
