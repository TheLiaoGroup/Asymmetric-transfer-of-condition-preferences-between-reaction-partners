"""Figure 2 and its Source Data, generated together from the public measurements.

Every plotted value is recomputed here from the released measurements with the analysis code's own estimators and must
agree with the approved analysis outputs to 1e-9 before anything is drawn:
  a, b  T of each molecule in each array (the Phase 1 estimator, on the condition patterns of anatomy.pattern) against
        05_DATA/Results/directional/phase1_T.csv; the molecule differences against phase1_D.csv; means, the 95% interval
        and the counts against phase1_summary.json
  c     the two SNAr role shares against Anatomy_arrays.csv, and the class subset that selects the N-H structures
        (anatomy.class_pattern on amplitude's nucleophile classes) against Anatomy.json
The structures follow 07_METHODS_AND_PROVENANCE/FIGURE_STRUCTURES_PROMOTION_READY.csv; their selection rules are applied
again to the released rows and must give the same molecules.

Public inputs, and nothing else with a measurement:
  05_DATA/Release/00_RAW_IMMUTABLE/DATASET_analysis.csv.gz, A1: rows of the 15 arrays with a numeric yield (177,308;
    Hiyama C46 excluded). Rows of the process-development (CMC) screens are dropped as the file is read.
  05_DATA/Analysis/directional_cohorts/DESIGN_COUNTS.csv, DC1: the frozen cohorts (identities and counts, no yields).
Writes 02_FIGURES/Main/Figure2.{pdf,svg,png}, 04_SOURCE_DATA/Figure2/Figure2_SourceData.xlsx and
Figure2_SourceData_provenance.json, and 08_REVIEW_AND_QC/figure2/figure2_plot_readback.json (what each panel draws,
read back from the drawn artists).
    python 06_CODE/Figures/figure2.py             all outputs in place
    python 06_CODE/Figures/figure2.py --out DIR   all outputs into DIR, for comparison with the delivered files
"""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                                              # noqa: E402
from matplotlib.colors import to_hex                                         # noqa: E402
from matplotlib.patches import Polygon, Rectangle                            # noqa: E402
from matplotlib.text import Text                                             # noqa: E402
from rdkit import Chem, RDLogger                                             # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "06_CODE/Core/00_common/ctb"))
import fig_style as fs                                                       # noqa: E402
import structure_art as sa                                                   # noqa: E402
from figure1 import matrix, public_units, write_xlsx                         # noqa: E402
from anatomy import FORMS, Shares, class_pattern, pattern                    # noqa: E402
import amplitude                                                             # noqa: E402  (after anatomy, which sets the path)
from history import ROLE                                                     # noqa: E402

RDLogger.DisableLog("rdApp.*")
RAW = ROOT / "05_DATA/Release/00_RAW_IMMUTABLE/DATASET_analysis.csv.gz"
COHORTS = ROOT / "05_DATA/Analysis/directional_cohorts/DESIGN_COUNTS.csv"
RES, DI = ROOT / "05_DATA/Results/directional", ROOT / "05_DATA/Results/display_inputs"
APPROVED = {"phase1_T.csv": RES / "phase1_T.csv", "phase1_D.csv": RES / "phase1_D.csv",
            "phase1_summary.json": RES / "phase1_summary.json", "Anatomy.json": DI / "Anatomy.json",
            "Anatomy_arrays.csv": DI / "Anatomy_arrays.csv"}
CHEMDRAW = ROOT / "02_FIGURES/ChemDraw"
SPEC = ROOT / "07_METHODS_AND_PROVENANCE/FIGURE_STRUCTURES_PROMOTION_READY.csv"
OUT = dict(pdf=ROOT / "02_FIGURES/Main/Figure2.pdf", svg=ROOT / "02_FIGURES/Main/Figure2.svg",
           png=ROOT / "02_FIGURES/Main/Figure2.png", xlsx=ROOT / "04_SOURCE_DATA/Figure2/Figure2_SourceData.xlsx",
           provenance=ROOT / "04_SOURCE_DATA/Figure2/Figure2_SourceData_provenance.json",
           readback=ROOT / "08_REVIEW_AND_QC/figure2/figure2_plot_readback.json")
AMIDE, BH1, BH2, SNAR = "Acylation|1", "Buchwald-Hartwig|1", "Buchwald-Hartwig|2", "SNAr|1"
SEED_BOOT, NB, SEED_CLASS = 20261004, 10000, 20261008
TOL = 1e-9
LABEL_PT = 7.0          # atom-label size of every structure in the figure, as in Figure 1
PNG_DPI = 1200
NAME = {"Cc1ccc(N)cc1": "4-methylaniline", "Cc1ccc(N)nc1": "5-methylpyridin-2-amine", "Cc1ccccc1N": "2-methylaniline",
        "COc1ccc(N)cc1": "4-methoxyaniline", "CSc1ccc(N)cc1": "4-(methylthio)aniline", "CSc1cccc(N)c1": "3-(methylthio)aniline",
        "Cc1ccc(N)c(C)c1": "2,4-dimethylaniline", "Nc1ccc2nccnc2c1": "quinoxalin-6-amine",
        "Cc1cc(F)ccc1C#N": "4-fluoro-2-methylbenzonitrile", "Cc1ccc(=O)[nH]c1": "5-methylpyridin-2(1H)-one",
        "c1cc2[nH]ccc2cn1": "1H-pyrrolo[3,2-c]pyridine", "C1CCNC1": "pyrrolidine", "CC1CCCNC1": "3-methylpiperidine"}
TITLE = dict(a="The same molecule carries more where its role owns",
             b="The same amines carry less in both C–N libraries",
             c="In SNAr the aryl halide owns")
TRANSFORM = dict(
    a="A1 rows of the arrays named in each row. In an array a pair is product-forming if its highest yield over the "
      "conditions is at least 0.10. Condition pattern = each product-forming pair's yields minus its least-squares line on "
      "the condition means (Methods). T(molecule, array) = mean, over ordered pairs of the molecule's product-forming pairs "
      "in that array, of the mean product of the two patterns over the conditions measured for both, divided by the array's "
      "mean squared pattern. Contexts from DESIGN_COUNTS.csv (frozen cohort G1, arrays whose owner holds at least twice its "
      "partner's share): owner = the array in which the molecule's role owns; non-owner = the arrays in which it does not; "
      "a context with two arrays takes the mean of their T. Difference = T(owner) - T(non-owner). 95% interval: 10,000 "
      "resamples of the 23 molecules (numpy default_rng, seed 20261004), 2.5th and 97.5th percentiles of the mean.",
    b="T as in panel a, for the eight amines measured in all three arrays (triple_context in DESIGN_COUNTS.csv), numbered "
      "M1-M8 by heavy atoms, then SMILES. C-N mean = (T in the ligand library + T in the metal library) / 2; R = T(amide) - C-N mean. The two C-N "
      "libraries share no condition identifier.",
    c="SNAr|1 A1 rows, product-forming pairs, patterns as in panel a. Shares: role shares of Figure 1d (Methods) for "
      "SNAr|1, as fractions of the pattern variance. Structures: class = amplitude.nucleophile class of the N-H nucleophile "
      "(role 1); the two most frequent classes with at least three pairs form the class subset (36 pairs), and the two "
      "molecules of each class with the fewest heavy atoms, then SMILES, are drawn.")
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
rel = lambda p: Path(p).resolve().relative_to(ROOT).as_posix()
_can = {}


def can(s):
    if s not in _can:
        m = Chem.MolFromSmiles(str(s))
        _can[s] = Chem.MolToSmiles(m, isomericSmiles=True) if m is not None else str(s)
    return _can[s]


heavy = lambda s: Chem.MolFromSmiles(s).GetNumHeavyAtoms()


# ====================================================================== values
def condition_patterns(A1, key):
    P, a, b = matrix(A1[A1.screen == key])
    M = P.to_numpy(float)
    act = np.nanmax(M, 1) >= FORMS
    return dict(P=P, M=M, act=act, E=pattern(M[act])[0], a=a[act], b=b[act], pairs=list(P.index[act]), rows=int((A1.screen == key).sum()))


def molecule_T(E, keys):
    """T for every molecule of one role (the Phase 1 estimator): mean pattern product over its ordered pairs of pairs
    sharing a condition, over the array's mean pattern variance."""
    Mk = np.isfinite(E).astype(float); E0 = np.nan_to_num(E); den = Mk @ Mk.T
    C = np.divide(E0 @ E0.T, den, out=np.zeros_like(den), where=den > 0)
    V = (np.nansum(E ** 2, 1) / Mk.sum(1)).mean()
    out = {}
    for u in np.unique(keys):
        idx = np.flatnonzero(keys == u)
        if len(idx) < 2:
            continue
        ok = (den[np.ix_(idx, idx)] > 0) & ~np.eye(len(idx), dtype=bool)
        if ok.sum():
            out[u] = (float(C[np.ix_(idx, idx)][ok].sum() / ok.sum() / V), int(ok.sum()))
    return out


def boot_ci(x, rng):
    x = np.asarray(x, float)
    return [float(v) for v in np.percentile(x[rng.integers(0, len(x), (NB, len(x)))].mean(1), [2.5, 97.5])]


def t_values(A1, DC):
    arrays = sorted(set(DC[(DC.set == "G1") & DC.array_clear_owner].array) | set(DC[DC.triple_context].array) | {SNAR})
    X, rows = {}, []
    for key in arrays:
        X[key] = condition_patterns(A1, key)
        for r, keys in ((1, X[key]["a"]), (2, X[key]["b"])):
            rows += [dict(molecule=can(u), array=key, role=r, T=t, n_ordered_pair_pairs=n) for u, (t, n) in molecule_T(X[key]["E"], keys).items()]
    TT = pd.DataFrame(rows)
    ref = pd.read_csv(APPROVED["phase1_T.csv"])
    ref = ref[ref.array.isin(arrays)]
    J = TT.merge(ref, on=["molecule", "array", "role"], how="outer", suffixes=("", "_approved"), indicator=True)
    assert (J._merge == "both").all() and len(J) == len(ref), J._merge.value_counts().to_dict()
    dT = float((J["T"] - J["T_approved"]).abs().max())
    assert dT < TOL and (J.n_ordered_pair_pairs == J.n_ordered_pair_pairs_approved).all(), dT
    return TT, X, dict(arrays=arrays, molecule_array_role_values=len(TT), max_abs_diff_T=dT)


def panel_a(TT, DC, S):
    Tmap = {(r.molecule, r.array, r.role): r.T for r in TT.itertuples()}
    prim = DC[(DC.set == "G1") & DC.array_clear_owner]
    rows, arr_rows = [], []
    for m, h in prim.groupby("molecule"):
        o, n = h[h.context == "owner"], h[h.context == "non-owner"]
        to = np.mean([Tmap[(m, r.array, r.role)] for r in o.itertuples()])
        tn = np.mean([Tmap[(m, r.array, r.role)] for r in n.itertuples()])
        rows.append(dict(molecule=m, stratum=h.stratum.iloc[0], owner_arrays="+".join(sorted(o.array)),
                         nonowner_arrays="+".join(sorted(n.array)), T_O=to, T_N=tn, D=to - tn))
        arr_rows += [dict(molecule=m, stratum=h.stratum.iloc[0], array=r.array, array_name=fs.ARRAY_LABEL[r.array], role=int(r.role),
                          role_name=r.role_label, context=r.context, T=Tmap[(m, r.array, r.role)], product_forming_partners=int(r.n_pf_partners),
                          ordered_pairs_of_pairs=int(r.n_ordered_pair_pairs)) for r in h.sort_values(["context", "array"], ascending=[False, True]).itertuples()]
    PM = pd.DataFrame(rows)                                                  # molecule order, as the frozen analysis resamples it
    assert len(PM) == 23
    ci = boot_ci(PM.D, np.random.default_rng(SEED_BOOT))
    am, ac = PM[PM.stratum == "amine"], PM[PM.stratum == "acid"]
    ref = pd.read_csv(APPROVED["phase1_D.csv"])
    mine = PM.sort_values(["stratum", "D"], ascending=[False, False]).reset_index(drop=True)
    assert mine.molecule.tolist() == ref.molecule.tolist() and (mine.owner_arrays == ref.owner_arrays).all() \
        and (mine.nonowner_arrays == ref.nonowner_arrays).all()
    dD = float(np.abs(mine[["T_O", "T_N", "D"]].to_numpy() - ref[["T_O", "T_N", "D"]].to_numpy()).max())
    assert dD < TOL, dD
    A, B, Dsc = S["A"], S["B"], S["descriptives"]
    stats = dict(n=len(PM), mean_D=float(PM.D.mean()), ci=ci, positive=int((PM.D > 0).sum()), mean_T_O=float(PM.T_O.mean()),
                 mean_T_N=float(PM.T_N.mean()),
                 amine=dict(n=len(am), mean_D=float(am.D.mean()), positive=int((am.D > 0).sum()), mean_T_O=float(am.T_O.mean()),
                            mean_T_N=float(am.T_N.mean())),
                 acid=dict(n=len(ac), mean_D=float(ac.D.mean()), positive=int((ac.D > 0).sum()), mean_T_O=float(ac.T_O.mean()),
                           mean_T_N=float(ac.T_N.mean())))
    assert (abs(stats["mean_D"] - A["mean"]) < TOL and np.allclose(ci, A["ci"], atol=1e-12, rtol=0) and stats["positive"] == A["positive"]
            and stats["n"] == A["n"] and abs(stats["amine"]["mean_D"] - B["amine_mean"]) < TOL
            and stats["amine"]["positive"] == B["amine_positive"] and abs(stats["acid"]["mean_D"] - B["acid_mean"]) < TOL
            and abs(stats["mean_T_O"] - Dsc["T_O_mean"]) < TOL and abs(stats["mean_T_N"] - Dsc["T_N_mean"]) < TOL)
    plot = pd.concat([PM[PM.stratum == s].sort_values(["D", "molecule"], ascending=[False, True]) for s in ("amine", "acid")],
                     ignore_index=True)
    plot.insert(0, "plot_row", np.arange(1, len(plot) + 1))
    return plot, pd.DataFrame(arr_rows), stats, dict(max_abs_diff_D=dD)


def panel_b(TT, DC, S, spec):
    Tmap = {(r.molecule, r.array, r.role): r.T for r in TT.itertuples()}
    tri = sorted(DC[DC.triple_context].molecule.unique(), key=lambda m: (heavy(m), m))
    labels = {r["smiles"]: r["label"] for r in spec if r["panel"] == "b"}
    assert [labels.get(m) for m in tri] == [f"M{i}" for i in range(1, 9)], "M1-M8 differ from the structure list"
    partners = DC.set_index(["molecule", "array"]).n_pf_partners
    rows = []
    for m in tri:
        ta, t1, t2 = Tmap[(m, AMIDE, 1)], Tmap[(m, BH1, 1)], Tmap[(m, BH2, 1)]
        rows.append(dict(label=labels[m], name=NAME[m], smiles=m, heavy_atoms=heavy(m), T_amide=ta, T_BH1=t1, T_BH2=t2,
                         CN_mean=(t1 + t2) / 2, R=ta - (t1 + t2) / 2, amide_above_CN_mean=bool(ta - (t1 + t2) / 2 > 0),
                         amide_above_both=bool(ta > t1 and ta > t2), partners_amide=int(partners[(m, AMIDE)]),
                         partners_BH1=int(partners[(m, BH1)]), partners_BH2=int(partners[(m, BH2)])))
    R = pd.DataFrame(rows)
    C = S["C"]
    ref = {r["molecule"]: r for r in C["rows"]}
    assert all(abs(round(getattr(r, k), 4) - ref[r.smiles][j]) < 1e-12 for r in R.itertuples()
               for k, j in (("T_amide", "T_amide"), ("T_BH1", "T_BH1"), ("T_BH2", "T_BH2"), ("R", "R")))
    stats = dict(mean_T=dict(amide=float(R.T_amide.mean()), BH1=float(R.T_BH1.mean()), BH2=float(R.T_BH2.mean())),
                 amide_minus_BH1=float((R.T_amide - R.T_BH1).mean()), amide_minus_BH2=float((R.T_amide - R.T_BH2).mean()),
                 mean_R=float(R.R.mean()), R_positive=int((R.R > 0).sum()), above_both=int(R.amide_above_both.sum()),
                 BH2_partners=sorted({int(v) for v in R.partners_BH2}), BH2_range=[float(R.T_BH2.min()), float(R.T_BH2.max())])
    assert (abs(stats["amide_minus_BH1"] - C["mean_amide_minus_BH1"]) < TOL and abs(stats["amide_minus_BH2"] - C["mean_amide_minus_BH2"]) < TOL
            and abs(stats["mean_R"] - C["mean_R"]) < TOL and stats["R_positive"] == C["n_R_positive"] == 6 and stats["above_both"] == 5
            and stats["BH2_partners"] == [2, 3])
    return R, stats


def panel_c(X, J, ARR, spec):
    x = X[SNAR]
    s = Shares(x["E"], x["a"], x["b"]).compute(np.ones(len(x["E"])))
    ref = ARR.set_index("screen").loc[SNAR]
    assert abs(s["1"] - ref.share1) < TOL and abs(s["2"] - ref.share2) < TOL
    assert ROLE[SNAR][0][0].startswith("N") and "halide" in ROLE[SNAR][1][0]
    classify = dict(amplitude.CLASS_ARRAYS)[SNAR]
    raw = [classify(m) for m in x["a"]]
    counts = pd.Series(raw).value_counts()
    top = [c for c in counts.index if c != "other" and counts[c] >= 3][:2]
    lab = [c if c in top else None for c in raw]
    ce = class_pattern(x["E"], x["a"], lab, np.random.default_rng(SEED_CLASS))
    approved = next(e for e in J["substrate_class"] if e["screen"] == SNAR)
    assert (abs(ce["explained"] - approved["explained"]) < TOL and abs(ce["null_median"] - approved["null_median"]) < TOL
            and abs(ce["null_p95"] - approved["null_p95"]) < TOL and ce["p"] == approved["p"] and ce["pairs"] == approved["pairs"] == 36
            and ce["molecules"] == approved["molecules"] and top == approved["classes"])
    subset = pd.DataFrame([dict(pair=p, n_h_nucleophile=a, nucleophile_class=c, aryl_halide=b)
                           for p, a, b, c in zip(x["pairs"], x["a"], x["b"], lab) if c is not None])
    # structures: the electrophile-side example and the two smallest molecules of each class (FIGURE_STRUCTURES rules)
    P = x["P"]; M = x["M"]
    full = [i for i in range(len(P)) if x["act"].tolist()[i] and np.isfinite(M[i]).all()]
    by = {}
    for i in full:
        nuc, hal = P.index[i].split(" >> ")
        by.setdefault(hal, []).append(nuc)
    best = min(by, key=lambda h: (-len(by[h]), heavy(h), h))
    c_spec = [r for r in spec if r["panel"] == "c"]
    assert M.shape[1] == 47 and can(best) == c_spec[0]["smiles"] and len(by[best]) == 4 \
        and {classify(n) for n in by[best]} == {"azole/lactam N–H"}
    smallest = {c: sorted(set(subset[subset.nucleophile_class == c].n_h_nucleophile), key=lambda m: (heavy(m), m))[:2] for c in top}
    shown = {"azole/lactam N–H": [r["smiles"] for r in c_spec if r["label"].startswith("azole")],
             "aliphatic amine": [r["smiles"] for r in c_spec if r["label"].startswith("aliphatic")]}
    assert all([can(m) for m in smallest[c]] == shown[c] for c in top)
    stats = dict(share_nucleophile=float(s["1"]), share_aryl_halide=float(s["2"]), pairs=ce["pairs"], molecules=ce["molecules"],
                 classes=top, product_forming=int(x["act"].sum()), pairs_all=len(P), conditions=int(M.shape[1]),
                 example=dict(smiles=can(best), fully_measured_product_forming_partners=len(by[best]), partners=sorted(by[best])))
    return subset, stats


# ====================================================================== drawing
W, H = fs.WIDTH_MM, 145.0                  # 6 mm less than with the two summary lines under titles a and b
AX0, AX1 = 3.0, 85.0                       # panel a axes, mm
A_YROW0, A_PITCH = 15.6, 3.2               # first amine row centre, row pitch
B0, B1, BY0, BY1 = 101.0, 163.0, 11.0, 72.0  # panel b axes
B_XLIM, B_YLIM = (-0.3, 2.3), (-0.5, 0.72)
STRIP_TOP, STRIP_X0, STRIP_X1 = 80.5, 92.5, 179.5
YC = 114.0                                 # top of panel c
XR, XS = 78.0, 122.0                       # panel c: left edge of the exact-molecule region, of its share axes
MS = 2.7                                   # marker size, pt


def a_rows():
    """Row centres (mm) of panel a: 19 amines, their mean, 4 acids, their mean."""
    am = [A_YROW0 + A_PITCH * i for i in range(19)]
    am_mean = am[-1] + 4.2
    ac0 = am_mean + 11.2
    ac = [ac0 + A_PITCH * i for i in range(4)]
    return am, am_mean, ac, ac[-1] + 4.2


def repel(values, gap):
    """Label positions near the values, at least gap apart, order kept (deterministic)."""
    order = np.argsort(values, kind="stable")
    pos = np.array(values, float)[order].copy()
    for _ in range(200):
        moved = False
        for i in range(1, len(pos)):
            if pos[i] - pos[i - 1] < gap:
                push = (gap - (pos[i] - pos[i - 1])) / 2
                pos[i - 1] -= push; pos[i] += push; moved = True
        if not moved:
            break
    out = np.empty_like(pos); out[order] = pos
    return out


def draw(PA, stats_a, R, stats_b, CL, stats_c, art):
    fs.apply(); fs.check_font()
    plt.rcParams["svg.hashsalt"] = "figure2"
    fig = plt.figure(figsize=(W * fs.MM, H * fs.MM))
    cv = fig.add_axes([0, 0, 1, 1]); cv.set_xlim(0, W); cv.set_ylim(H, 0); cv.axis("off"); cv.set_gid("canvas")
    put = lambda x, y, s, gid=None, **kw: cv.text(x, y, s, gid=gid, **kw)
    for L, x, y in (("a", 0, 0.6), ("b", 92.0, 0.6), ("c", 0, YC)):           # 0.6 mm: no glyph touches the top edge
        fs.letter(fig, x, y, L, H).set_gid(f"letter|{L}")
        fs.title(fig, x + 4.2, y + 0.35, TITLE[L], H).set_gid(f"title|{L}")
    small = dict(fontsize=fs.PT["small"], va="center")
    tiny = dict(fontsize=fs.PT["tiny"], va="center")
    owner, faint, muted = fs.C["owner"], fs.C["faint"], fs.C["muted"]

    # ---------------------------------------------------------------- a: the same molecule in two reactions
    sa_ = stats_a
    for gid, x, filled, s in (("owner", AX0 + 0.9, True, "where the role owns"), ("nonowner", AX0 + 29.0, False, "where it does not")):
        cv.plot([x], [5.9], marker="o", ms=MS, mew=0.6, mec=owner, mfc=owner if filled else "white", ls="none", gid=f"a|key|{gid}")
        put(x + 1.4, 5.9, s, f"a|key_label|{gid}", ha="left", **tiny)
    am_y, am_mean_y, ac_y, ac_mean_y = a_rows()
    bottom = ac_mean_y + 2.2
    ax = fs.axes_mm(fig, AX0, 8.5, AX1 - AX0, bottom - 8.5, H); ax.set_gid("a|axes")
    ax.set_xlim(-0.5, 0.8); ax.set_ylim(bottom, 8.5)
    for s, y0, y1 in (("amine", am_y[0] - 1.6, am_mean_y + 1.6), ("acid", ac_y[0] - 1.6, ac_mean_y + 1.6)):
        ax.plot([0, 0], [y0, y1], color=faint, lw=0.4, zorder=0, gid=f"a|zero|{s}")      # under the rows only, not the headers
    ys = am_y + ac_y
    for i, r in PA.iterrows():
        y, shape = ys[i], ("o" if r.stratum == "amine" else "s")
        ax.plot([r.T_N, r.T_O], [y, y], color=faint, lw=0.7, solid_capstyle="butt", zorder=1, gid=f"a|link|{r.plot_row}")
        ax.plot([r.T_N], [y], marker=shape, ms=MS, mew=0.6, mec=owner, mfc="white", ls="none", zorder=3, gid=f"a|nonowner|{r.plot_row}")
        ax.plot([r.T_O], [y], marker=shape, ms=MS, mew=0.6, mec=owner, mfc=owner, ls="none", zorder=3, gid=f"a|owner|{r.plot_row}")
    blocks = (("amine", "o", am_y[0] - 6.0, am_mean_y, f"{sa_['amine']['n']} amines: mean difference {sa_['amine']['mean_D']:.3f}, larger "
               f"for {sa_['amine']['positive']} of {sa_['amine']['n']}",
               "owner in amide coupling; non-owner in C–N coupling, SNAr or Sonogashira coupling"),
              ("acid", "s", ac_y[0] - 6.0, ac_mean_y, f"{sa_['acid']['n']} carboxylic acids: mean difference {sa_['acid']['mean_D']:.3f}, "
               f"larger for {sa_['acid']['positive']} of {sa_['acid']['n']}", "owner in esterification; non-owner in amide coupling"))
    for s, shape, yh, ym, l1, l2 in blocks:
        cv.plot([AX0 + 0.9], [yh], marker=shape, ms=MS, mew=0.6, mec=owner, mfc=owner, ls="none", gid=f"a|header_marker|{s}")
        put(AX0 + 2.3, yh, l1, f"a|header|{s}|1", ha="left", **small)
        put(AX0 + 2.3, yh + 2.6, l2, f"a|header|{s}|2", ha="left", color=muted, **tiny)
        ax.plot([-0.5, 0.8], [ym - 2.1, ym - 2.1], color=fs.C["grid"], lw=0.4, zorder=0, gid=f"a|rule|{s}")
        st = sa_[s]
        ax.plot([st["mean_T_N"], st["mean_T_O"]], [ym, ym], color=muted, lw=1.3, solid_capstyle="butt", zorder=1, gid=f"a|mean_link|{s}")
        ax.plot([st["mean_T_N"]], [ym], marker=shape, ms=MS + 1.2, mew=0.8, mec=owner, mfc="white", ls="none", zorder=3,
                gid=f"a|mean_nonowner|{s}")
        ax.plot([st["mean_T_O"]], [ym], marker=shape, ms=MS + 1.2, mew=0.8, mec=owner, mfc=owner, ls="none", zorder=3,
                gid=f"a|mean_owner|{s}")
        ax.text(-0.485, ym, "mean", gid=f"a|mean_label|{s}", ha="left", color=muted, **tiny)
    ax.set_xticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8]); ax.set_yticks([])
    ax.set_xlabel("T, share of the array's condition-pattern variance carried across partners", labelpad=1.5)
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    ax.patch.set_visible(False)

    # ---------------------------------------------------------------- b: eight amines in three libraries
    sb = stats_b
    bx = fs.axes_mm(fig, B0, BY0, B1 - B0, BY1 - BY0, H); bx.set_gid("b|axes")
    bx.set_xlim(*B_XLIM); bx.set_ylim(*B_YLIM)
    col_mm = lambda k: B0 + (k - B_XLIM[0]) / (B_XLIM[1] - B_XLIM[0]) * (B1 - B0)
    bx.axhline(0, color=faint, lw=0.4, zorder=0, gid="b|zero")
    for k in range(3):
        bx.axvline(k, color=fs.C["grid"], lw=0.4, zorder=0, gid=f"b|column_line|{k}")
    accent = fs.C["accent"]
    colour = {r.label: (accent if r.R < 0 else fs.C["partner"]) for r in R.itertuples()}
    for r in R.itertuples():
        bx.plot([0, 1, 2], [r.T_amide, r.T_BH1, r.T_BH2], color=colour[r.label], lw=0.8, marker="o", ms=2.2, mew=0, zorder=2,
                gid=f"b|line|{r.label}")
    mt = sb["mean_T"]
    bx.plot([0, 1, 2], [mt["amide"], mt["BH1"], mt["BH2"]], color=owner, lw=1.6, marker="o", ms=3.0, mew=0, zorder=4, gid="b|mean")
    bx.text(-0.07, mt["amide"], "mean", ha="right", color=owner, gid="b|mean_label", **tiny)
    ylab = repel(R.T_BH2.to_numpy(), 0.047)
    for r, yl in zip(R.itertuples(), ylab):
        bx.plot([2.03, 2.11], [r.T_BH2, yl], color=colour[r.label], lw=0.4, zorder=1, clip_on=False, gid=f"b|leader|{r.label}")
        bx.text(2.13, yl, r.label, ha="left", color=accent if r.R < 0 else fs.C["text"], clip_on=False, gid=f"b|label|{r.label}", **tiny)
    for k, (name, own) in enumerate(((fs.ARRAY_LABEL[AMIDE], "amine owns"), (fs.ARRAY_LABEL[BH1], "aryl halide owns"),
                                     (fs.ARRAY_LABEL[BH2], "aryl halide owns"))):     # the array names of Figure 1
        put(col_mm(k), BY1 + 2.0, name, f"b|column|{k}", ha="center", **small)
        put(col_mm(k), BY0 - 1.6, own, f"b|owner_status|{k}", ha="center", color=owner if k == 0 else muted, **tiny)
    yb = BY0 - 4.6
    cv.plot([col_mm(1), col_mm(1), col_mm(2), col_mm(2)], [yb + 0.8, yb, yb, yb + 0.8], color=muted, lw=0.4, gid="b|bracket")
    put((col_mm(1) + col_mm(2)) / 2, yb - 1.3, "no condition in common", "b|bracket_label", ha="center", color=muted, **tiny)
    bx.set_xticks([]); bx.set_yticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6])
    bx.set_ylabel("T (as in a)", labelpad=1.5)
    for sp in ("right", "top", "bottom"):
        bx.spines[sp].set_visible(False)
    bx.patch.set_visible(False)
    # structure strip: tops aligned, equal gaps
    wid = {r.label: sa.size_mm(art[r.label], LABEL_PT) for r in R.itertuples()}
    gap = (STRIP_X1 - STRIP_X0 - sum(w for w, _ in wid.values())) / (len(R) - 1)
    x = STRIP_X0
    hmax = max(h for _, h in wid.values())
    for r in R.itertuples():
        w, h = wid[r.label]
        sa.draw(cv, art[r.label], x, STRIP_TOP, LABEL_PT, gid=f"b|struct|{r.label}")
        put(x + w / 2, STRIP_TOP + hmax + 2.0, r.label, f"b|strip_label|{r.label}", ha="center",
            color=accent if r.R < 0 else fs.C["text"], **small)
        x += w + gap

    # ---------------------------------------------------------------- c: the two SNAr roles
    sc = stats_c
    put(4.2, YC + 5.4, "The library varies only the base or acid and the solvent.", "c|subtitle", ha="left",
        color=muted, **small)
    put(4.2, YC + 9.0, "N–H nucleophiles", "c|region|class", ha="left", fontweight="bold", **small)
    put(XR + 0.5, YC + 9.0, "Exact-molecule ownership", "c|region|exact", ha="left", fontweight="bold", **small)
    cv.plot([XR - 4.0, XR - 4.0], [YC + 7.5, H - 1.5], color=fs.C["grid"], lw=0.6, gid="c|divider")
    top_s = YC + 11.8
    groups = (("azole/lactam N–H", ("azolelactam1", "azolelactam2"), 4.2), ("aliphatic amine", ("aliphatic1", "aliphatic2"), 37.5))
    for cls, names, x0 in groups:
        x = x0
        for n in names:
            w, h = sa.size_mm(art[n], LABEL_PT)
            sa.draw(cv, art[n], x, top_s, LABEL_PT, gid=f"c|struct|{n}")
            x += w + 1.6
        put((x0 + x - 1.6) / 2, H - 2.6, cls, f"c|class_label|{cls}", ha="center", **small)
    # exact-molecule shares
    w, h = sa.size_mm(art["aryl"], LABEL_PT)
    sa.draw(cv, art["aryl"], XR + 1.0, top_s + 1.0, LABEL_PT, gid="c|struct|aryl")
    put(XR + 1.0 + w / 2, top_s + h + 3.4, "example aryl halide", "c|example_label", ha="center", color=muted, **tiny)
    ex = fs.axes_mm(fig, XS, top_s + 0.5, 172.0 - XS, 8.0, H); ex.set_gid("c|axes|exact")
    ex.set_xlim(0, 0.5); ex.set_ylim(2, 0)
    for j, (k, v, colr, lbl) in enumerate((("aryl halide", sc["share_aryl_halide"], owner, "aryl halide"),
                                            ("N-H nucleophile", sc["share_nucleophile"], fs.C["partner"], "N–H nucleophile"))):
        ex.add_patch(Rectangle((0, j + 0.18), v, 0.64, facecolor=colr, edgecolor="none", lw=0, gid=f"c|share|{k}"))
        ex.text(max(v, 0) + 0.012, j + 0.5, f"{v:.3f}", ha="left", gid=f"c|share_value|{k}", clip_on=False, **tiny)
        ex.text(-0.014, j + 0.5, lbl, ha="right", gid=f"c|share_label|{k}", clip_on=False, **tiny)
    ex.set_xticks([0, 0.2, 0.4]); ex.set_yticks([])
    ex.set_xlabel("Share carried by the exact molecule", labelpad=1.0, fontsize=fs.PT["tiny"])
    ex.tick_params(axis="x", labelsize=fs.PT["tiny"])
    for sp in ("left", "right", "top"):
        ex.spines[sp].set_visible(False)
    ex.patch.set_visible(False)
    return fig


# ====================================================================== read back what was drawn
def readback(fig):
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    px_mm = 25.4 / fig.dpi

    def bbox(a):
        e = a.get_window_extent(rend)
        return [round(e.x0 * px_mm, 3), round(H - e.y1 * px_mm, 3), round(e.x1 * px_mm, 3), round(H - e.y0 * px_mm, 3)]
    out = dict(width_mm=W, height_mm=H, rectangles=[], markers=[], lines=[], texts=[], polygons={})
    for ax in fig.axes:
        for p in ax.patches:
            g = p.get_gid()
            if not g:
                continue
            if isinstance(p, Polygon):
                out["polygons"][g.rsplit("|", 1)[0]] = out["polygons"].get(g.rsplit("|", 1)[0], 0) + 1
            elif isinstance(p, Rectangle):
                out["rectangles"].append(dict(gid=g, x=p.get_x(), y=p.get_y(), width=p.get_width(), height=p.get_height(),
                                              color=to_hex(p.get_facecolor()), bbox_mm=bbox(p)))
        for ln in ax.lines:
            if ln.get_gid():
                rec = dict(gid=ln.get_gid(), x=np.asarray(ln.get_xdata(), float).tolist(), y=np.asarray(ln.get_ydata(), float).tolist(),
                           color=to_hex(ln.get_color()), bbox_mm=bbox(ln))
                if ln.get_marker() not in (None, "None", ""):
                    rec.update(marker=ln.get_marker(), face=to_hex(ln.get_markerfacecolor()), size=ln.get_markersize())
                    out["markers"].append(rec)
                else:
                    out["lines"].append(rec)
    for t in fig.findobj(Text):
        if t.get_visible() and t.get_text().strip():
            out["texts"].append(dict(gid=t.get_gid(), text=t.get_text(), size=t.get_fontsize(), bbox_mm=bbox(t)))
    return out


# ====================================================================== Source Data
def source_data(PA, PT, sa_, R, sb, CL, sc, spec):
    disp = lambda s: " + ".join(fs.ARRAY_LABEL[k] for k in s.split("+"))
    a = PA.assign(owner_context=PA.owner_arrays.map(disp), nonowner_context=PA.nonowner_arrays.map(disp),
                  larger_in_owner_context=PA.D > 0)
    a = a.rename(columns={"stratum": "class", "molecule": "smiles", "T_O": "T_owner_context", "T_N": "T_nonowner_context",
                          "D": "difference"})[["plot_row", "class", "smiles", "owner_context", "T_owner_context", "nonowner_context",
                                               "T_nonowner_context", "difference", "larger_in_owner_context"]]
    a.insert(0, "row_type", "molecule")
    summ = [dict(row_type="mean over the 23 molecules", T_owner_context=sa_["mean_T_O"], T_nonowner_context=sa_["mean_T_N"],
                 difference=sa_["mean_D"]),
            dict(row_type="95% interval of the mean difference, lower", difference=sa_["ci"][0]),
            dict(row_type="95% interval of the mean difference, upper", difference=sa_["ci"][1]),
            dict(row_type="molecules larger in the owner context", count=sa_["positive"], of=sa_["n"])]
    for s, nm in (("amine", "amines"), ("acid", "carboxylic acids")):
        st = sa_[s]
        summ += [dict(row_type=f"mean over the {st['n']} {nm} (mean row)", **{"class": s}, T_owner_context=st["mean_T_O"],
                      T_nonowner_context=st["mean_T_N"], difference=st["mean_D"]),
                 dict(row_type=f"{nm} larger in the owner context", **{"class": s}, count=st["positive"], of=st["n"])]
    a = pd.concat([a, pd.DataFrame(summ)], ignore_index=True).astype({"plot_row": "Int64", "count": "Int64", "of": "Int64"})
    t = PT.rename(columns={"molecule": "smiles", "stratum": "class"})
    order = dict(zip(PA.molecule, PA.plot_row))
    t.insert(0, "plot_row", t.smiles.map(order))
    t = t.sort_values(["plot_row", "context", "array"], ascending=[True, False, True], ignore_index=True)
    b = R.rename(columns={"label": "plot_label"})[["plot_label", "name", "smiles", "heavy_atoms", "T_amide", "T_BH1", "T_BH2", "CN_mean",
                                                   "R", "amide_above_CN_mean", "amide_above_both", "partners_amide", "partners_BH1",
                                                   "partners_BH2"]]
    b.insert(0, "row_type", "amine")
    b["line_colour"] = ["highlighted (R < 0)" if v < 0 else "standard" for v in R.R]
    bs = pd.DataFrame([dict(row_type="mean over the 8 amines (mean line)", T_amide=sb["mean_T"]["amide"], T_BH1=sb["mean_T"]["BH1"],
                            T_BH2=sb["mean_T"]["BH2"], R=sb["mean_R"]),
                       dict(row_type="mean T(amide) - T(C–N ligand library)", R=sb["amide_minus_BH1"]),
                       dict(row_type="mean T(amide) - T(C–N metal library)", R=sb["amide_minus_BH2"]),
                       dict(row_type="amines with R > 0", count=sb["R_positive"], of=8),
                       dict(row_type="amines above both C-N libraries", count=sb["above_both"], of=8)])
    b = pd.concat([b, bs], ignore_index=True).astype({"heavy_atoms": "Int64", "partners_amide": "Int64", "partners_BH1": "Int64",
                                                      "partners_BH2": "Int64", "count": "Int64", "of": "Int64"})
    c_shares = pd.DataFrame([dict(bar="aryl halide", role=ROLE[SNAR][1][0], share=sc["share_aryl_halide"], array=SNAR,
                                  product_forming_pairs=sc["product_forming"]),
                             dict(bar="N-H nucleophile", role=ROLE[SNAR][0][0], share=sc["share_nucleophile"], array=SNAR,
                                  product_forming_pairs=sc["product_forming"])])
    names = sa_names(spec)
    st = pd.DataFrame([dict(panel=r["panel"], label=r["label"], drawing=names[i], name=NAME[r["smiles"]], smiles=r["smiles"],
                            role=r["role"], rule=r["rule"], chemdraw_file=f"{names[i]}.cdxml") for i, r in enumerate(spec)])
    readme = pd.DataFrame([
        dict(sheet="Fig2a", panel="a", rows="one per molecule (23), in plot order; then the means, interval and counts printed",
             units="T as a fraction of the array's pattern variance", transformation=TRANSFORM["a"]),
        dict(sheet="Fig2a_T", panel="a", rows="one per molecule and array (the T values behind each context)",
             units="T as a fraction", transformation=TRANSFORM["a"]),
        dict(sheet="Fig2b", panel="b", rows="one per amine (8, M1-M8); then the mean line and the printed differences and counts",
             units="T as a fraction", transformation=TRANSFORM["b"]),
        dict(sheet="Fig2c_class_pairs", panel="c", rows="one per pair of the class subset (36), from which the N-H structures "
             "are chosen", units="-",
             transformation=TRANSFORM["c"]),
        dict(sheet="Fig2c_shares", panel="c", rows="the two bars", units="fraction of the pattern variance",
             transformation=TRANSFORM["c"]),
        dict(sheet="Structures", panel="b, c", rows="one per structure drawn (13)", units="-",
             transformation="Molecules and selection rules of FIGURE_STRUCTURES_PROMOTION_READY.csv, drawn by ChemDraw "
                            "(02_FIGURES/ChemDraw); the rules are re-applied to the released rows by the generator."),
        dict(sheet="(all)", panel="-", rows="-", units="-",
             transformation="Measurements: DATASET_analysis.csv.gz of the data release, rows of the 15 recurring-molecule arrays "
                            "with a numeric yield (Hiyama C46 excluded); cohorts: DESIGN_COUNTS.csv (identities and counts, no "
                            "yields). Input hashes, counts and the agreement with the analysis outputs: "
                            "Figure2_SourceData_provenance.json."),
    ])
    return dict(README=readme, Fig2a=a, Fig2a_T=t, Fig2b=b, Fig2c_class_pairs=CL, Fig2c_shares=c_shares,
                Structures=st)


def sa_names(spec):
    import chemdraw_structures                                              # file names only; no ChemDraw needed
    return chemdraw_structures.drawing_names(spec)


# ====================================================================== main
def main(out=None):
    paths = {k: (Path(out) / v.name if out else v) for k, v in OUT.items()}
    for p in paths.values():
        p.parent.mkdir(parents=True, exist_ok=True)
    S = json.loads(APPROVED["phase1_summary.json"].read_text(encoding="utf-8"))
    J = json.loads(APPROVED["Anatomy.json"].read_text(encoding="utf-8"))
    ARR = pd.read_csv(APPROVED["Anatomy_arrays.csv"])
    DC = pd.read_csv(COHORTS)
    spec = [r for r in csv.DictReader(open(SPEC, encoding="utf-8")) if r["figure"] == "2"]
    names = sa_names(spec)
    man = {r["name"]: r for r in csv.DictReader(open(CHEMDRAW / "structures_manifest.csv", encoding="utf-8"))}
    art = {}
    for r, n in zip(spec, names):
        m = man[n]
        assert m["smiles"] == r["smiles"], n
        a = sa.load(CHEMDRAW / f"{n}.svg")
        assert a["sha256"] == m["svg_sha256"], f"{n}: SVG differs from the ChemDraw manifest"
        art[r["label"] if r["panel"] == "b" else n.split("_")[-1]] = a
    A1, _ = public_units()
    TT, X, t_info = t_values(A1, DC)
    PA, PT, stats_a, a_info = panel_a(TT, DC, S)
    R, stats_b = panel_b(TT, DC, S, spec)
    CL, stats_c = panel_c(X, J, ARR, spec)

    fig = draw(PA, stats_a, R, stats_b, CL, stats_c, art)
    rb = readback(fig)
    fig.savefig(paths["pdf"], metadata={"Creator": None, "Producer": None, "CreationDate": None})
    fig.savefig(paths["svg"], metadata={"Creator": None, "Date": None})
    fig.savefig(paths["png"], dpi=PNG_DPI, metadata={"Software": None})
    plt.close(fig)
    paths["readback"].write_text(json.dumps(rb, indent=1, ensure_ascii=False), encoding="utf-8")
    write_xlsx(source_data(PA, PT, stats_a, R, stats_b, CL, stats_c, spec), paths["xlsx"])
    prov = dict(
        figure="Figure 2 | Reaction context reassigns molecular ownership", generated_by=rel(__file__),
        inputs=dict(
            measurements=dict(file=rel(RAW), sha256=sha(RAW), A1=dict(rows=int(len(A1)), selection="screens of the 15 arrays (no "
                              "process-development CMC screen), numeric yield, not Hiyama|1 C46"),
                              rows_used={k: v["rows"] for k, v in X.items()}),
            cohorts=dict(file=rel(COHORTS), sha256=sha(COHORTS)),
            approved_outputs={n: dict(file=rel(p), sha256=sha(p)) for n, p in APPROVED.items()},
            structures={n: dict(file=rel(CHEMDRAW / f"{n}.svg"), sha256=sha(CHEMDRAW / f"{n}.svg"),
                                cdxml_sha256=sha(CHEMDRAW / f"{n}.cdxml")) for n in names},
            structure_spec=dict(file=rel(SPEC), sha256=sha(SPEC))),
        panels=dict(a=dict(t_info, **a_info, records=len(PA), stats=stats_a, transformation=TRANSFORM["a"]),
                    b=dict(records=len(R), stats=stats_b, transformation=TRANSFORM["b"]),
                    c=dict(records=dict(class_pairs=len(CL)), stats=stats_c, transformation=TRANSFORM["c"])),
        outputs={k: dict(file=rel(OUT[k]), sha256=sha(p)) for k, p in paths.items() if k != "provenance"})
    paths["provenance"].write_text(json.dumps(prov, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"Figure 2: A1 {len(A1)} rows; T max |diff| {t_info['max_abs_diff_T']:.1e}; D max |diff| {a_info['max_abs_diff_D']:.1e}; "
          f"wrote {', '.join(p.name for p in paths.values())}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out")
    main(ap.parse_args().out)
