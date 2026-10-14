"""Figure 3 and its Source Data, generated together from the public measurements.

Every plotted value is recomputed here from the released measurements with the scoring code's own estimators and must agree
with the frozen evaluation outputs before anything is drawn:
  b  per-target squared errors of the owner-only and partner-only predictions (partner_external.train and .score, single-role
     profiles fitted on the analysis rows of each array) against record_transfer_losses.csv (Set A) and target_losses.csv
     (Set B); pooled fractions, D, per-array D and the whole-array intervals (partner_external.summarise) against
     reserve_validation.json and cross_library.json. Set A's interval comes from the generator state of the frozen run, which
     first drew 2,000 molecule resamples for the amide ownership estimate; those draws are replayed (anatomy.molecule_weights).
     The 2,000 whole-array draws that follow, and the 2,000 owning-building-block draws after them (the interval reported in
     Supplementary Note 10), are written out with the Source Data, so both intervals can be recomputed from them and the
     per-target errors alone.
  c  the same targets with each history cut to two partners (the Phase 2 rule) against phase2_setAB_k2.csv and
     phase2_summary.json.
Public inputs, and nothing else with a measurement:
  05_DATA/Release/00_RAW_IMMUTABLE/DATASET_analysis.csv.gz, A1 (177,308 rows; process-development CMC rows dropped on reading):
    the analysis arrays, which fix each molecule's history
  05_DATA/Validation_Set_A/never_analysed_reserve.csv.gz, SA1 (333 pairs, 29,516 rows; the other reserve pairs are dropped on
    reading, by the release rule: record-transfer targets and product-forming amide pairs)
  05_DATA/Validation_Set_B/never_opened_cross_library.csv.gz, SB1 (535 pairs, 30,556 rows)
Each array's owner role is read from the frozen anatomy output (Anatomy_arrays.csv).
Writes 02_FIGURES/Main/Figure3.{pdf,svg,png}, 04_SOURCE_DATA/Figure3/Figure3_SourceData.xlsx,
SetA_whole_array_bootstrap_draws.csv, SetA_owner_block_bootstrap_draws.csv and Figure3_SourceData_provenance.json, and
08_REVIEW_AND_QC/figure3/figure3_plot_readback.json.
    python 06_CODE/Figures/figure3.py             all outputs in place
    python 06_CODE/Figures/figure3.py --out DIR   all outputs into DIR, for comparison with the delivered files
"""
import argparse
import hashlib
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                                              # noqa: E402
from matplotlib.colors import to_hex                                         # noqa: E402
from matplotlib.patches import Rectangle                                     # noqa: E402
from matplotlib.text import Text                                             # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path[:0] = [str(ROOT / "06_CODE/Core/00_common"), str(ROOT / "06_CODE/Core/00_common/ctb"),
                str(ROOT / "06_CODE/Validation/06_validation"), str(ROOT / "06_CODE/Validation/11_partner_asymmetry")]
import fig_style as fs                                                       # noqa: E402
from figure1 import public_units, write_xlsx                                 # noqa: E402
import anatomy                                                               # noqa: E402
import history                                                               # noqa: E402
import partner_external as pe                                                # noqa: E402

DATA = ROOT / "05_DATA"
RAW = DATA / "Release/00_RAW_IMMUTABLE/DATASET_analysis.csv.gz"
SET_A, SET_B = DATA / "Validation_Set_A/never_analysed_reserve.csv.gz", DATA / "Validation_Set_B/never_opened_cross_library.csv.gz"
APPROVED = {"record_transfer_losses.csv": DATA / "Results/reserve_validation/record_transfer_losses.csv",
            "reserve_validation.json": DATA / "Results/reserve_validation/reserve_validation.json",
            "target_losses.csv": DATA / "Results/cross_library/target_losses.csv",
            "cross_library.json": DATA / "Results/cross_library/cross_library.json",
            "phase2_setAB_k2.csv": DATA / "Results/directional/phase2_setAB_k2.csv",
            "phase2_summary.json": DATA / "Results/directional/phase2_summary.json",
            "Anatomy_arrays.csv": DATA / "Results/display_inputs/Anatomy_arrays.csv"}
OUT = dict(pdf=ROOT / "02_FIGURES/Main/Figure3.pdf", svg=ROOT / "02_FIGURES/Main/Figure3.svg",
           png=ROOT / "02_FIGURES/Main/Figure3.png", xlsx=ROOT / "04_SOURCE_DATA/Figure3/Figure3_SourceData.xlsx",
           draws=ROOT / "04_SOURCE_DATA/Figure3/SetA_whole_array_bootstrap_draws.csv",
           owner_draws=ROOT / "04_SOURCE_DATA/Figure3/SetA_owner_block_bootstrap_draws.csv",
           provenance=ROOT / "04_SOURCE_DATA/Figure3/Figure3_SourceData_provenance.json",
           readback=ROOT / "08_REVIEW_AND_QC/figure3/figure3_plot_readback.json")
N_SA1, N_SB1 = 29516, 30556
SEED_A, SEED_B, SEED_K2, NBOOT, K2_MAX = 20261004, 20261007, 20261004, 2000, 45
TOL = 1e-9
PNG_DPI = 1200
TITLE = dict(a="The owner was fixed from the analysis arrays",
             b="The owner's history predicts out-of-array pairs better",
             c="The advantage remains with two partners on each side")
TRANSFORM = dict(
    a="Design only, no measurement. Targets: pairs of Set A (SA1) and Set B (SB1) whose two molecules each have at least two "
      "distinct partners in the analysis rows (A1) of their array and that have at least eight measured cells.",
    b="Per array, condition means s and slope-adjusted patterns of the analysis rows (A1); a molecule's profile in a role = sum "
      "of its pairs' patterns / (number of observed contributions + 4), for molecules with at least two distinct partners. A "
      "target's pattern e = residual of its yields on an intercept and s over the conditions it covers; each profile is "
      "projected the same way. SSE = sum of squared errors of e against the projected owner (or partner) profile, SST = sum of "
      "e squared. Fraction explained = 1 - sum SSE / sum SST over a set's targets; D = (sum SSE partner - sum SSE owner) / "
      "sum SST; per-array D the same within an array. 95% interval: 2,000 draws resampling whole arrays (Set A: the frozen "
      "evaluation's draws, listed in SetA_whole_array_bootstrap_draws.csv; Set B: seed 20261007). Set B: a condition measured "
      "twice for a pair is dropped.",
    c="As panel b, with each side's profile built from exactly two of its molecule's analysis pairs: all two-pair subsets when "
      "there are at most 45, otherwise the 45 whose sorted partner identities have the smallest SHA-256; squared errors are "
      "averaged over subsets before pooling. Intervals: 2,000 whole-array draws, seed 20261004 for each set. History depth = "
      "the number of analysis pairs of the molecule. Post hoc sensitivity, specified after both sets had been scored.")
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
sha_text = lambda t: hashlib.sha256(t.encode("utf-8")).hexdigest()
rel = lambda p: Path(p).resolve().relative_to(ROOT).as_posix()


# ====================================================================== values
def public_sets():
    """SA1 and SB1 as the release defines them; SA1 by its rule (record-transfer targets and product-forming amide pairs)."""
    SA = pd.read_csv(SET_A, low_memory=False)
    key = SA.screen + "|" + SA.subpair
    tl = pd.read_csv(APPROVED["record_transfer_losses.csv"])
    best = pd.to_numeric(SA["yield"], errors="coerce").groupby(key).max()
    keep = set(tl.screen + "|" + tl.subpair) | {k for k, v in best.items() if k.startswith("Acylation|1|") and v >= anatomy.FORMS}
    SA = SA[key.isin(keep)].assign(**{"yield": lambda d: pd.to_numeric(d["yield"], errors="coerce")})
    SB = pd.read_csv(SET_B, low_memory=False)
    SB["yield"] = pd.to_numeric(SB["yield"], errors="coerce")
    assert len(SA) == N_SA1 and len(SB) == N_SB1, (len(SA), len(SB))
    return SA, SB


def models(A1, screens):
    """partner_external's single-role profiles on the analysis rows of each array (A1), with the patterns behind them."""
    out = {}
    for key in screens:
        G = A1[A1.screen == key]
        P = G.pivot(index="subpair", columns="condition_id", values="yield")
        meta = G.drop_duplicates("subpair").set_index("subpair").loc[P.index]
        can = {s: pe.canonical(s) for s in set(meta.sub_1_smiles) | set(meta.sub_2_smiles)}
        Y, s1, s2 = P.to_numpy(float), meta.sub_1_smiles.map(can).to_numpy(), meta.sub_2_smiles.map(can).to_numpy()
        s, prof = pe.train(Y, s1, s2)
        E = pe.slope_residuals(Y, s)
        out[key] = dict(s=s, prof=prof, conds=list(P.columns), s1=s1, s2=s2, E=E, Mk=~np.isnan(E), pairs=len(Y), rows=len(G))
    return out


def subsets_k2(ids):
    c = list(itertools.combinations(range(len(ids)), 2))
    return c if len(c) <= K2_MAX else sorted(c, key=lambda ij: sha_text("\t".join(sorted((ids[ij[0]], ids[ij[1]])))))[:K2_MAX]


def losses(S, M, owner, dedup, set_name):
    rows = []
    for (key, sp), g in S.groupby(["screen", "subpair"], sort=False):
        m = M[key]
        a, b = pe.canonical(g.sub_1_smiles.iloc[0]), pe.canonical(g.sub_2_smiles.iloc[0])
        if a not in m["prof"][1] or b not in m["prof"][2]:
            continue
        z = g.set_index("condition_id")["yield"]
        if dedup:
            z = z[~z.index.duplicated(keep=False)]
        y = np.array([z.get(c, np.nan) for c in m["conds"]], float)
        r = pe.score(y, m["s"], m["prof"], a, b, owner[key])
        if not r:
            continue
        o = ~np.isnan(y) & ~np.isnan(m["s"])
        X = np.c_[np.ones(o.sum()), m["s"][o]]
        proj = lambda v: v - X @ np.linalg.lstsq(X, v, rcond=None)[0]
        e = proj(y[o])
        k2 = {}
        for side, role in (("e", owner[key]), ("n", 3 - owner[key])):
            keys, other = (m["s1"], m["s2"]) if role == 1 else (m["s2"], m["s1"])
            rws = np.flatnonzero(keys == (a if role == 1 else b))
            sse = [float(np.sum((e - proj((np.nansum(m["E"][rws[[i, j]]], 0) / (m["Mk"][rws[[i, j]]].sum(0) + pe.SHRINK))[o])) ** 2))
                   for i, j in subsets_k2(list(other[rws]))]
            k2[f"sse_{side}_k2"], k2[f"subsets_{side}"], k2[f"history_{side}"] = float(np.mean(sse)), len(sse), len(rws)
        rows.append(dict(set=set_name, screen=key, subpair=sp, electrophile=a if owner[key] == 1 else b, **r, **k2))
    return pd.DataFrame(rows)


def check_losses(L, approved, cols):
    J = L.merge(approved, on=["screen", "subpair"], how="outer", suffixes=("", "_approved"), indicator=True)
    assert (J._merge == "both").all() and len(J) == len(approved), J._merge.value_counts().to_dict()
    d = max(float((J[c] - J[c + "_approved"]).abs().max()) for c in cols)
    assert d < TOL and (J.electrophile == J.electrophile_approved).all(), d
    return d


def summarise_k2(L, seed):
    """Phase 2's summary of the two-partner losses: pooled D, whole-array interval, per-array D."""
    L = L.assign(gain=L.sse_n_k2 - L.sse_e_k2)
    A = L.groupby("screen")[["gain", "sst"]].sum()
    ia = np.random.default_rng(seed).integers(0, len(A), (NBOOT, len(A)))
    num, den = A.gain.to_numpy()[ia].sum(1), A.sst.to_numpy()[ia].sum(1)
    whole = np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)
    per = (A.gain / A.sst).to_dict()
    return dict(targets=len(L), arrays=len(A), D=float(L.gain.sum() / L.sst.sum()),
                r2_owner=float(1 - L.sse_e_k2.sum() / L.sst.sum()), r2_partner=float(1 - L.sse_n_k2.sum() / L.sst.sum()),
                ci_whole_array=[float(np.nanpercentile(whole, 2.5)), float(np.nanpercentile(whole, 97.5))],
                arrays_positive=int(sum(v > 0 for v in per.values())), per_array={k: float(v) for k, v in per.items()})


def frozen_run_state(SA):
    """The generator of the frozen Set A run, after the amide ownership estimate's 2,000 molecule resamples."""
    rng = np.random.default_rng(SEED_A)
    for key, A in history.partner_arrays(SA).items():
        keep = np.nanmax(A["M"], 1) >= anatomy.FORMS
        if keep.sum() >= 50:                                                  # the estimate's eligibility rule
            for _ in range(NBOOT):
                anatomy.molecule_weights(A["a"][keep], A["b"][keep], rng)
    return rng


def frozen_draws(LA, S, rng):
    """The frozen Set A run's draws, taken from the state it is given in partner_external.summarise's order: first the whole
    arrays, then, array by array, the owning building blocks within each array. Two tables with one row per draw: the arrays
    drawn, or in each array the owner molecules drawn (';' between them), and the draw's D. The 2.5th and 97.5th percentiles
    of D must be the frozen intervals."""
    L = LA.assign(gain=LA.sse_n - LA.sse_e)
    G = L.groupby("screen")[["gain", "sst"]].sum()                                     # summarise's units, in its order
    ia = rng.integers(0, len(G), (NBOOT, len(G)))
    d = G.gain.to_numpy()[ia].sum(1) / G.sst.to_numpy()[ia].sum(1)
    assert np.array_equal(np.percentile(d, [2.5, 97.5]), S["ci_whole_array"]), (np.percentile(d, [2.5, 97.5]), S["ci_whole_array"])
    whole = pd.DataFrame(G.index.to_numpy()[ia], columns=[f"screen_{j + 1}" for j in range(len(G))])
    num, den, owners = np.zeros(NBOOT), np.zeros(NBOOT), {}
    for key, g in L.groupby("screen"):
        c = g.groupby("electrophile")[["gain", "sst"]].sum()
        assert not any(";" in m for m in c.index)
        ic = rng.integers(0, len(c), (NBOOT, len(c)))
        num += c.gain.to_numpy()[ic].sum(1); den += c.sst.to_numpy()[ic].sum(1)
        owners[key] = [";".join(r) for r in c.index.to_numpy()[ic]]
    e = num / den
    assert np.array_equal(np.percentile(e, [2.5, 97.5]), S["ci_electrophile_cluster"]), (np.percentile(e, [2.5, 97.5]),
                                                                                       S["ci_electrophile_cluster"])
    block = pd.DataFrame(owners)
    for t in (whole, block):
        t.insert(0, "draw", np.arange(1, NBOOT + 1))
    return whole.assign(D=d), block.assign(D=e)


def compute(A1, SA, SB):
    owner = pd.read_csv(APPROVED["Anatomy_arrays.csv"]).set_index("screen").owner_role.astype(int).to_dict()
    M = models(A1, sorted(set(SA.screen) | set(SB.screen)))
    LA, LB = losses(SA, M, owner, False, "Set A"), losses(SB, M, owner, True, "Set B")
    cols = ["sse_e", "sse_n", "sst"]
    dA = check_losses(LA, pd.read_csv(APPROVED["record_transfer_losses.csv"]), cols)
    dB = check_losses(LB, pd.read_csv(APPROVED["target_losses.csv"]), cols)
    k2_ref = pd.read_csv(APPROVED["phase2_setAB_k2.csv"])
    k2cols = ["sse_e_k2", "sse_n_k2", "subsets_e", "subsets_n", "history_e", "history_n"]
    dK = max(check_losses(L, k2_ref[k2_ref.set == n].drop(columns="set"), cols + k2cols) for n, L in (("Set A", LA), ("Set B", LB)))
    SAs = pe.summarise(LA.drop(columns="set"), NBOOT, frozen_run_state(SA))
    SBs = pe.summarise(LB.drop(columns="set"), NBOOT, np.random.default_rng(SEED_B))
    for mine, ref in ((SAs, json.loads(APPROVED["reserve_validation.json"].read_text(encoding="utf-8"))["e2"]),
                      (SBs, json.loads(APPROVED["cross_library.json"].read_text(encoding="utf-8")))):
        for k in ("targets", "arrays", "arrays_positive"):
            assert mine[k] == ref[k], k
        for k in ("D", "r2_electrophile", "r2_partner"):
            assert abs(mine[k] - ref[k]) < TOL, k
        for k in ("ci_whole_array", "ci_electrophile_cluster"):
            assert np.allclose(mine[k], ref[k], atol=1e-12, rtol=0), k
        assert mine["per_array"].keys() == ref["per_array"].keys() and all(abs(mine["per_array"][k] - v) < TOL for k, v in ref["per_array"].items())
    draws, owner_draws = frozen_draws(LA, SAs, frozen_run_state(SA))
    P2 =json.loads(APPROVED["phase2_summary.json"].read_text(encoding="utf-8"))["7.1"]
    K = {}
    for n, L in (("Set A", LA), ("Set B", LB)):
        K[n] = summarise_k2(L, SEED_K2)
        ref = P2[n]
        assert (abs(K[n]["D"] - ref["k2"]["D"]) < TOL and np.allclose(K[n]["ci_whole_array"], ref["k2"]["ci_whole_array"], atol=1e-12, rtol=0)
                and K[n]["arrays_positive"] == ref["k2"]["arrays_positive"] and abs(K[n]["r2_owner"] - ref["r2_owner_k2"]) < TOL
                and abs(K[n]["r2_partner"] - ref["r2_partner_k2"]) < TOL
                and all(abs(K[n]["per_array"][k] - v) < TOL for k, v in ref["k2"]["per_array"].items()))
        K[n]["owner_history_deeper"] = int((L.history_e > L.history_n).sum())
        K[n]["partner_history_deeper"] = int((L.history_n > L.history_e).sum())
        assert (K[n]["owner_history_deeper"], K[n]["partner_history_deeper"]) == (ref["targets_owner_history_deeper"],
                                                                                 ref["targets_partner_history_deeper"])
    union = len(set(LA.screen) | set(LB.screen))
    design = dict(analysis_arrays=int(A1.screen.nunique()), set_A=dict(pairs=len(LA), arrays=int(LA.screen.nunique())),
                  set_B=dict(pairs=len(LB), arrays=int(LB.screen.nunique())), union_pairs=len(LA) + len(LB), union_systems=union)
    assert (design["set_A"], design["set_B"], design["union_pairs"], union) == (dict(pairs=309, arrays=6), dict(pairs=535, arrays=13), 844, 14)
    info = dict(max_abs_diff=dict(set_A_losses=dA, set_B_losses=dB, two_partner_losses=dK),
                analysis_rows_used={k: v["rows"] for k, v in M.items()})
    return LA, LB, dict(A=SAs, B=SBs), K, design, info, draws, owner_draws


# ====================================================================== drawing
W, H = fs.WIDTH_MM, 88.0
XB = 60.0                                   # left edge of panels b and c
BAR0, BAR1, D0, D1 = 72.0, 113.0, 131.0, 175.0
ROWS = dict(b=(9.5, 35.0), c=(53.5, 79.0))  # axes top and bottom of b and c
YB = 44.0                                   # top of panel c
CAT = ("Set A", "Set B")
DX_LIM, DX_DOTS = (-0.4, 2.6), (0.0, 1.4)    # D axes: x range, centre of each set's dots


def offsets(vals, dy=0.017, dx=0.07):
    """Deterministic swarm for the per-array dots: each dot, in order of value, takes the smallest sideways offset at which it
    keeps one dot width (dx in x units, dy in y units) from every dot already placed."""
    placed, out = [], {}
    for i in np.argsort(vals, kind="stable"):
        for off in [0.0] + [s_ * k * dx for k in range(1, 20) for s_ in (1, -1)]:
            if all(abs(vals[i] - pv) >= dy or abs(off - po) >= dx for pv, po in placed):
                placed.append((vals[i], off)); out[i] = off
                break
    return np.array([out[i] for i in range(len(vals))])


def draw(design, S, K):
    fs.apply(); fs.check_font()
    plt.rcParams["svg.hashsalt"] = "figure3"
    fig = plt.figure(figsize=(W * fs.MM, H * fs.MM))
    cv = fig.add_axes([0, 0, 1, 1]); cv.set_xlim(0, W); cv.set_ylim(H, 0); cv.axis("off"); cv.set_gid("canvas")
    put = lambda x, y, s, gid=None, **kw: cv.text(x, y, s, gid=gid, **kw)
    for L, x, y in (("a", 0, 0.6), ("b", XB - 2.0, 0.6), ("c", XB - 2.0, YB)):
        fs.letter(fig, x, y, L, H).set_gid(f"letter|{L}")
        fs.title(fig, x + 4.2, y + 0.35, TITLE[L], H).set_gid(f"title|{L}")
    small = dict(fontsize=fs.PT["small"], va="center")
    tiny = dict(fontsize=fs.PT["tiny"], va="center")
    owner, partner, faint, muted = fs.C["owner"], fs.C["partner"], fs.C["faint"], fs.C["muted"]

    # ---------------------------------------------------------------- a: design (cells are a diagram, not measurements)
    c, gx, ncol = 2.6, 4.5, 11
    cell = lambda x, y, colr, gid: cv.add_patch(Rectangle((x, y), c - 0.4, c - 0.4, facecolor=colr, edgecolor="none", gid=gid))
    put(gx, 10.2, "Analysis arrays", "a|array_label", ha="left", fontweight="bold", **small)
    for i in range(6):                                                         # rows: pairs; columns: conditions
        for j in range(ncol):
            cell(gx + j * c, 12.0 + i * c, fs.C["noise"], f"a|array_cell|{i}|{j}")
    put(gx, 30.6, f"{design['analysis_arrays']} arrays: each array's owner and every molecule's", "a|array_note|1", ha="left", **tiny)
    put(gx, 33.0, "history are fixed here", "a|array_note|2", ha="left", **tiny)
    for i in range(2):                                                         # Set A: out-of-array pairs on at least 90% of the library
        for j in range(ncol):
            cell(gx + j * c, 38.0 + i * c, muted, f"a|setA_cell|{i}|{j}")
    put(gx, 46.6, f"Set A: {design['set_A']['pairs']} reserve pairs", "a|setA|1", ha="left", fontweight="bold", **small)
    put(gx, 49.2, f"{design['set_A']['arrays']} arrays, ≥90% of each condition library", "a|setA|2", ha="left", **tiny)
    for i, cols in ((0, (0, 1, 3, 4, 7, 9)), (1, (2, 3, 5, 6, 8))):           # Set B: same library, partial panels
        for j in cols:
            cell(gx + j * c, 54.0 + i * c, faint, f"a|setB_cell|{i}|{j}")
    put(gx, 62.6, f"Set B: {design['set_B']['pairs']} pairs", "a|setB|1", ha="left", fontweight="bold", **small)
    put(gx, 65.2, f"{design['set_B']['arrays']} arrays, same libraries, outside the dense arrays", "a|setB|2", ha="left", **tiny)
    put(gx, 72.6, f"Together {design['union_pairs']} pairs outside the analysed arrays", "a|union|1", ha="left", **small)
    put(gx, 75.4, f"in {design['union_systems']} reaction systems", "a|union|2", ha="left", **small)

    # ---------------------------------------------------------------- b and c: owner against partner
    def panel(tag, R2, Dset, top_bottom):
        y0, y1 = top_bottom
        bx = fs.axes_mm(fig, BAR0, y0, BAR1 - BAR0, y1 - y0, H); bx.set_gid(f"{tag}|axes|bars")
        bx.set_xlim(-0.55, 1.55); bx.set_ylim(0, 0.25)
        for k, name in enumerate(CAT):
            for side, dx, colr in (("owner", -0.18, owner), ("partner", 0.18, partner)):
                v = R2[name][side]
                bx.add_patch(Rectangle((k + dx - 0.16, 0), 0.32, v, facecolor=colr, edgecolor="none", lw=0, gid=f"{tag}|bar|{name}|{side}"))
                bx.text(k + dx, v + 0.006, f"{v:.3f}", ha="center", va="bottom", fontsize=fs.PT["tiny"], gid=f"{tag}|bar_value|{name}|{side}")
        bx.set_xticks([0, 1], CAT); bx.set_yticks([0, 0.1, 0.2])
        bx.set_ylabel("Fraction explained", labelpad=1.5)
        for sp in ("right", "top"):
            bx.spines[sp].set_visible(False)
        bx.patch.set_visible(False)
        dx_ = fs.axes_mm(fig, D0, y0, D1 - D0, y1 - y0, H); dx_.set_gid(f"{tag}|axes|D")
        dx_.set_xlim(*DX_LIM); dx_.set_ylim(-0.1, 0.32)
        dx_.axhline(0, color=faint, lw=0.4, zorder=0, gid=f"{tag}|zero")
        for k, name in enumerate(CAT):
            x0 = DX_DOTS[k]
            per = Dset[name]["per_array"]
            keys = sorted(per)
            vals = np.array([per[q] for q in keys])
            for q, v, off in zip(keys, vals, offsets(vals)):
                dx_.plot([x0 + off], [v], marker="o", ms=2.6, mew=0, mfc=owner if v > 0 else faint, ls="none", zorder=3,
                         gid=f"{tag}|array|{name}|{q}")
            d, ci = Dset[name]["D"], Dset[name]["ci_whole_array"]
            dx_.plot([x0 + 0.30, x0 + 0.30], ci, color=fs.C["text"], lw=0.8, zorder=4, gid=f"{tag}|ci|{name}")
            dx_.plot([x0 + 0.20, x0 + 0.40], [d, d], color=fs.C["text"], lw=1.4, solid_capstyle="butt", zorder=4, gid=f"{tag}|D|{name}")
            dx_.text(x0 + 0.47, d + 0.018, f"{d:.3f}", ha="left", va="center", fontsize=fs.PT["tiny"], gid=f"{tag}|D_value|{name}",
                     clip_on=False)
            dx_.text(x0 + 0.47, d - 0.018, f"{ci[0]:.3f}–{ci[1]:.3f}", ha="left", va="center", fontsize=fs.PT["tiny"], color=muted,
                     gid=f"{tag}|D_ci|{name}", clip_on=False)
            dx_.text(x0, 0.305, f"{Dset[name]['arrays_positive']} of {Dset[name]['arrays']} > 0", ha="center", va="center",
                     fontsize=fs.PT["tiny"], color=muted, gid=f"{tag}|positive|{name}")
        dx_.set_xticks([x + 0.15 for x in DX_DOTS], CAT); dx_.set_yticks([-0.1, 0, 0.1, 0.2, 0.3])
        dx_.set_ylabel("Owner − partner, D", labelpad=1.5)
        for sp in ("right", "top"):
            dx_.spines[sp].set_visible(False)
        dx_.patch.set_visible(False)
        return bx, dx_

    r2b = {n: dict(owner=S[k]["r2_electrophile"], partner=S[k]["r2_partner"]) for n, k in zip(CAT, "AB")}
    panel("b", r2b, {n: S[k] for n, k in zip(CAT, "AB")}, ROWS["b"])
    r2c = {n: dict(owner=K[n]["r2_owner"], partner=K[n]["r2_partner"]) for n in CAT}
    panel("c", r2c, K, ROWS["c"])
    # one line, at the height of panel b's key row; the history-depth counts are given in the legend and the Results
    put(XB + 2.2, YB + 5.8, "Matched-depth sensitivity; two partners per side", "c|note", ha="left", color=muted, **tiny)
    for j,(s, colr) in enumerate((("owner's history", owner), ("partner's history", partner))):
        x = BAR0 + 2.0 + 21.0 * j
        cv.add_patch(Rectangle((x, 5.8 - 0.7), 1.4, 1.4, facecolor=colr, edgecolor="none", gid=f"b|key|{j}"))
        put(x + 2.0, 5.8, s, f"b|key_label|{j}", ha="left", **tiny)
    for j, (s, colr) in enumerate((("array, owner ahead", owner), ("array, partner ahead", faint))):
        x = D0 + 1.0 + 21.0 * j
        cv.plot([x + 0.7], [5.8], marker="o", ms=2.6, mew=0, mfc=colr, ls="none", gid=f"b|dot_key|{j}")
        put(x + 2.0, 5.8, s, f"b|dot_key_label|{j}", ha="left", **tiny)
    return fig


# ====================================================================== read back what was drawn
def readback(fig):
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    px_mm = 25.4 / fig.dpi

    def bbox(a):
        e = a.get_window_extent(rend)
        return [round(e.x0 * px_mm, 3), round(H - e.y1 * px_mm, 3), round(e.x1 * px_mm, 3), round(H - e.y0 * px_mm, 3)]
    out = dict(width_mm=W, height_mm=H, rectangles=[], markers=[], lines=[], texts=[])
    for ax in fig.axes:
        for p in ax.patches:
            if p.get_gid() and isinstance(p, Rectangle):
                out["rectangles"].append(dict(gid=p.get_gid(), x=p.get_x(), y=p.get_y(), width=p.get_width(), height=p.get_height(),
                                              color=to_hex(p.get_facecolor()), bbox_mm=bbox(p)))
        for ln in ax.lines:
            if ln.get_gid():
                rec = dict(gid=ln.get_gid(), x=np.asarray(ln.get_xdata(), float).tolist(), y=np.asarray(ln.get_ydata(), float).tolist(),
                           color=to_hex(ln.get_color()), bbox_mm=bbox(ln))
                if ln.get_marker() not in (None, "None", ""):
                    rec.update(marker=ln.get_marker(), face=to_hex(ln.get_markerfacecolor()))
                    out["markers"].append(rec)
                else:
                    out["lines"].append(rec)
    for t in fig.findobj(Text):
        if t.get_visible() and t.get_text().strip():
            out["texts"].append(dict(gid=t.get_gid(), text=t.get_text(), size=t.get_fontsize(), bbox_mm=bbox(t)))
    return out


# ====================================================================== Source Data
def source_data(LA, LB, S, K, design):
    T = pd.concat([LA, LB], ignore_index=True).rename(columns={
        "electrophile": "owner_molecule", "sse_e": "sse_owner", "sse_n": "sse_partner", "sse_e_k2": "sse_owner_two_partners",
        "sse_n_k2": "sse_partner_two_partners", "subsets_e": "owner_subsets", "subsets_n": "partner_subsets",
        "history_e": "owner_history_pairs", "history_n": "partner_history_pairs"})
    T.insert(2, "array", T.screen.map(fs.ARRAY_LABEL))
    a = pd.DataFrame([dict(set="Set A", targets=design["set_A"]["pairs"], arrays=design["set_A"]["arrays"],
                           array_names="; ".join(sorted(LA.screen.map(fs.ARRAY_LABEL).unique())), coverage="at least 90% of the condition library for every pair"),
                      dict(set="Set B", targets=design["set_B"]["pairs"], arrays=design["set_B"]["arrays"],
                           array_names="; ".join(sorted(LB.screen.map(fs.ARRAY_LABEL).unique())),
                           coverage="same libraries, partial panels outside the dense arrays"),
                      dict(set="together", targets=design["union_pairs"], arrays=design["union_systems"]),
                      dict(set="analysis arrays (histories, owners)", arrays=design["analysis_arrays"])])
    rows_b, rows_c = [], []
    for n, k in zip(CAT, "AB"):
        s = S[k]
        rows_b += [dict(set=n, row_type="pooled", owner_fraction=s["r2_electrophile"], partner_fraction=s["r2_partner"], D=s["D"],
                        D_ci_lower=s["ci_whole_array"][0], D_ci_upper=s["ci_whole_array"][1], targets=s["targets"],
                        arrays_positive=s["arrays_positive"], arrays=s["arrays"])]
        rows_b += [dict(set=n, row_type="array", screen=q, array=fs.ARRAY_LABEL[q], D=v, targets=int(((LA if k == "A" else LB).screen == q).sum()))
                   for q, v in sorted(s["per_array"].items())]
        kk = K[n]
        rows_c += [dict(set=n, row_type="pooled", owner_fraction=kk["r2_owner"], partner_fraction=kk["r2_partner"], D=kk["D"],
                        D_ci_lower=kk["ci_whole_array"][0], D_ci_upper=kk["ci_whole_array"][1], targets=kk["targets"],
                        arrays_positive=kk["arrays_positive"], arrays=kk["arrays"],
                        targets_partner_history_deeper=kk["partner_history_deeper"], targets_owner_history_deeper=kk["owner_history_deeper"])]
        rows_c += [dict(set=n, row_type="array", screen=q, array=fs.ARRAY_LABEL[q], D=v) for q, v in sorted(kk["per_array"].items())]
    ints = ("targets", "arrays_positive", "arrays", "targets_partner_history_deeper", "targets_owner_history_deeper")
    b = pd.DataFrame(rows_b).astype({c_: "Int64" for c_ in ints if c_ in pd.DataFrame(rows_b)})
    c = pd.DataFrame(rows_c).astype({c_: "Int64" for c_ in ints if c_ in pd.DataFrame(rows_c)})
    readme = pd.DataFrame([
        dict(sheet="Fig3a", panel="a", rows="one per set, then the union and the analysis arrays", units="counts", transformation=TRANSFORM["a"]),
        dict(sheet="Fig3bc_targets", panel="b, c", rows="one per scored target (309 + 535)", units="squared errors of yield "
             "fractions; history in analysis pairs", transformation=TRANSFORM["b"] + " " + TRANSFORM["c"]),
        dict(sheet="Fig3b", panel="b", rows="per set: the pooled row (bars, D, interval, count), then one per array (dots)",
             units="fractions of the pattern; D as a fraction", transformation=TRANSFORM["b"]),
        dict(sheet="Fig3c", panel="c", rows="per set: the pooled row, then one per array", units="fractions; D as a fraction",
             transformation=TRANSFORM["c"]),
        dict(sheet="SetA_whole_array_bootstrap_draws.csv (separate file)", panel="b",
             rows="one per whole-array draw of Set A (2,000): the six arrays drawn, with repeats, and the draw's D",
             units="screen identifiers as in Fig3bc_targets; D as a fraction",
             transformation="The frozen evaluation's own draws. D of a draw = sum over the drawn arrays of (sse_partner - "
                            "sse_owner) / the same sum of sst, from the Set A rows of Fig3bc_targets. The 95% interval of panel b "
                            "(Set A) = the 2.5th and 97.5th percentiles of the 2,000 values (linear interpolation)."),
        dict(sheet="SetA_owner_block_bootstrap_draws.csv (separate file)", panel="Supplementary Note 10",
             rows="one per owning-building-block draw of Set A (2,000): for each array, the owner molecules drawn with repeats "
                  "within it, and the draw's D",
             units="one column per array (screen identifier); owner molecules as in owner_molecule of Fig3bc_targets, separated "
                   "by ';'; D as a fraction",
             transformation="The frozen evaluation's own draws, which follow its whole-array draws in the same stream. D of a "
                            "draw = sum over the arrays and the owners drawn in each of (sse_partner - sse_owner) of that owner's "
                            "Set A targets in Fig3bc_targets / the same sum of sst. The owning-building-block interval of Set A "
                            "(Supplementary Note 10) = the 2.5th and 97.5th percentiles of the 2,000 values (linear "
                            "interpolation)."),
        dict(sheet="(all)", panel="-", rows="-", units="-",
             transformation="Measurements: units A1 (the analysis arrays), SA1 (Set A) and SB1 (Set B) of the data release. "
                            "Owner roles: Anatomy_arrays.csv. "
                            "Input hashes and the agreement with the frozen outputs: Figure3_SourceData_provenance.json.")])
    return dict(README=readme, Fig3a=a.astype({"targets": "Int64", "arrays": "Int64"}), Fig3bc_targets=T, Fig3b=b, Fig3c=c)


# ====================================================================== main
def main(out=None):
    paths = {k: (Path(out) / v.name if out else v) for k, v in OUT.items()}
    for p in paths.values():
        p.parent.mkdir(parents=True, exist_ok=True)
    A1, _ = public_units()
    SA, SB = public_sets()
    LA, LB, S, K, design, info, draws, owner_draws = compute(A1, SA, SB)
    fig = draw(design, S, K)
    rb = readback(fig)
    fig.savefig(paths["pdf"], metadata={"Creator": None, "Producer": None, "CreationDate": None})
    fig.savefig(paths["svg"], metadata={"Creator": None, "Date": None})
    fig.savefig(paths["png"], dpi=PNG_DPI, metadata={"Software": None})
    plt.close(fig)
    paths["readback"].write_text(json.dumps(rb, indent=1, ensure_ascii=False), encoding="utf-8")
    write_xlsx(source_data(LA, LB, S, K, design), paths["xlsx"])
    for k, t in (("draws", draws), ("owner_draws", owner_draws)):
        t.to_csv(paths[k], index=False, lineterminator="\n")
        assert pd.read_csv(paths[k], float_precision="round_trip").equals(t), f"{paths[k].name} does not read back exactly"
    prov = dict(
        figure="Figure 3 | Ownership predicts combinations beyond the analysed arrays", generated_by=rel(__file__),
        inputs=dict(measurements={n: dict(file=rel(p), sha256=sha(p)) for n, p in (("A1", RAW), ("SA1", SET_A), ("SB1", SET_B))},
                    rows=dict(A1=int(len(A1)), SA1=int(len(SA)), SB1=int(len(SB)), analysis_rows_used=info["analysis_rows_used"]),
                    approved_outputs={n: dict(file=rel(p), sha256=sha(p)) for n, p in APPROVED.items()}),
        panels=dict(a=dict(design=design, transformation=TRANSFORM["a"]),
                    b=dict(set_A={k: v for k, v in S["A"].items()}, set_B={k: v for k, v in S["B"].items()}, transformation=TRANSFORM["b"]),
                    c=dict(K, transformation=TRANSFORM["c"])),
        agreement_with_frozen_outputs=info["max_abs_diff"],
        outputs={k: dict(file=rel(OUT[k]), sha256=sha(p)) for k, p in paths.items() if k != "provenance"})
    paths["provenance"].write_text(json.dumps(prov, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"Figure 3: A1 {len(A1)}, SA1 {len(SA)}, SB1 {len(SB)} rows; losses max |diff| "
          f"{max(info['max_abs_diff'].values()):.1e}; wrote {', '.join(p.name for p in paths.values())}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out")
    main(ap.parse_args().out)
