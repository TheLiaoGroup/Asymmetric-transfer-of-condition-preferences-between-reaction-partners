"""Record-level provenance of the 244,436 released measurements under the versioned provenance rule PROV-1+A2.

Implements 07_METHODS_AND_PROVENANCE/provenance_public/PROV-1+A2_PUBLIC_RULE.md, the final rule only. Each released
measurement (P2) is matched exactly and one to one against the dataset of the companion manuscript (P1) and against the
published amide-coupling collection of Zhang et al. (2025). Writes P2_RECORD_PROVENANCE.csv (one row per released
measurement) and P1_P2_ZHANG_OVERLAP_COUNTS_PUBLIC.json to 07_METHODS_AND_PROVENANCE/provenance_public/, or to --out.
    python 06_CODE/Release/provenance_overlap_public.py --p1 <companion data package> --p4 <yield-prediction data package> [--zhang <all_HTE_with_condition.csv>] [--out DIR]
--zhang takes the file of the pinned commit, plain or gzip-compressed; it defaults to the project's gzip copy. The run stops
unless the bytes are the commit's.
Without --p1 (the companion package is not part of this release) nothing is rewritten: the released records are keyed as the
rule keys them and compared row for row with the released P2_RECORD_PROVENANCE.csv, and, if the Zhang file is supplied, the
Zhang flags are recomputed and compared; the companion flags and the counts that involve the companion dataset need its package.
--p4 names the data package of the yield-prediction study whose released dataset contains every measurement released here
(ref. 31 of the manuscript): a folder holding its released records, data/p4_records.csv.gz, as published in its GitHub repository.
It is read with --p1, and adds p4_row_uid: the record's row_uid in that released file. Only the records whose screen is a library
of this release are used. That identifier is not part of PROV-1+A2. Its key is the library, the condition label, the two partners
as RDKit canonical isomeric SMILES, taken as an unordered pair, and the yield rounded as PROV-1+A2 rounds it; records are paired one
to one, the k-th record of a key in released order with the k-th of that key in the P4 file, in file order.
"""
import argparse
import ast
import gzip
import hashlib
import io
import itertools
import json
import platform
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem, RDLogger

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "05_DATA/Release/00_RAW_IMMUTABLE"
DATA = ROOT / "05_DATA"
OUT = ROOT / "07_METHODS_AND_PROVENANCE/provenance_public"
RULE = OUT / "PROV-1+A2_PUBLIC_RULE.md"
VERSION = "PROV-1+A2"
ZHANG = dict(repository="github.com/aichemeco/amide_coupling", commit="f6fdbe283fc24af53389aab8dcba5e46386ac9c0",
             path="data/all_HTE_with_condition.csv", git_blob="778db5176ed6a1b3b11cbb80da4eacbd40ea74c5", bytes=17657403,
             sha256="5085e1ae50ee26173bd170a4b10a99c4e3101a46b706d115e425aff17c091b41",
             citation="Zhang, C. et al. Chem. Sci. 16, 11809-11822 (2025). doi:10.1039/d5sc03364k")
UNITS = dict(A1=177308, A3=74, SA1=36498, SB1=30556)
N = dict(P2=244436, P1=165131, Zhang=47015)
P4_RECORDS = "data/p4_records.csv.gz"                  # the released records of ref. 31, as in its GitHub repository
REGIONS = ("P1 only", "P2 only", "Zhang only", "P1 ∩ P2 only", "P1 ∩ Zhang only", "P2 ∩ Zhang only", "P1 ∩ P2 ∩ Zhang")
NAMES = [f"Reagent{i}_Standardized Name" for i in range(1, 7)]
SEP = "\x1f"
RDLogger.DisableLog("rdApp.*")
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
_can, _mol = {}, {}


def canon(s):
    if not isinstance(s, str):
        return None
    if s not in _can:
        m = Chem.MolFromSmiles(s)
        _mol[s], _can[s] = m, (Chem.MolToSmiles(m, isomericSmiles=True) if m is not None else None)
    return _can[s]


def yield_key(t, two_stage=True):
    """The stored yield rounded to nine and then six decimals; two_stage=False rounds once, for the guard only."""
    try:
        v = float(t)
    except (TypeError, ValueError):
        return None
    if v != v:
        return None
    return round(round(v, 9) if two_stage else v, 6) + 0.0


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def key(lib, label, s1, s2, y):
    return None if any(x is None for x in (lib, label, s1, s2, y)) else SEP.join((lib, label, s1, s2, repr(y)))


def condition_table(path):
    C = pd.read_csv(path, dtype=str)
    assert C.groupby("rxn_class").reaction_type.nunique().max() == 1 and C.groupby(["reaction_type", "condition_id"]).size().max() == 1
    C["names"] = [tuple(x for x in r if isinstance(x, str) and x.strip()) for r in C[NAMES].itertuples(index=False)]
    return C.set_index(["reaction_type", "condition_id"]), dict(C.drop_duplicates("rxn_class")[["rxn_class", "reaction_type"]].to_numpy())


def role_table():
    src = (ROOT / "06_CODE/Core/00_common/ctb/history.py").read_text(encoding="utf-8")
    return next(ast.literal_eval(n.value) for n in ast.parse(src).body
                if isinstance(n, ast.Assign) and any(getattr(t, "id", "") == "ROLE" for t in n.targets))


def share(smiles, screens, ROLE, role):
    """Share of the distinct molecules matching the role's SMARTS under any of the screens."""
    pats = [Chem.MolFromSmarts(ROLE[s][role][1]) for s in screens]
    ms = [_mol[s] for s in set(smiles) if canon(s) is not None]
    return float(np.mean([any(m.HasSubstructMatch(p) for p in pats) for m in ms])) if ms else 0.0


# ------------------------------------------------------------------ the three datasets
def released(T, lib_of):
    """P2: the released measurements, unit by unit (A1, A3, SA1, SB1) in released order."""
    A = pd.read_csv(RAW / "DATASET_analysis.csv.gz", dtype=str)
    y = pd.to_numeric(A["yield"], errors="coerce")
    process = A.screen.str.contains("CMC", regex=False)                  # process-development arrays, not released
    dup = (A.screen == "Hiyama|1") & (A.condition_id == "C46")
    SA = pd.read_csv(DATA / "Validation_Set_A/never_analysed_reserve.csv.gz", dtype=str)
    sk = SA.screen + "|" + SA.subpair
    tl = pd.read_csv(DATA / "Results/reserve_validation/record_transfer_losses.csv", dtype=str)
    best = pd.to_numeric(SA["yield"], errors="coerce").groupby(sk).max()
    dt = pd.read_csv(DATA / "Results/reserve_validation/decision_targets.csv", dtype=str)
    keep = (set(tl.screen + "|" + tl.subpair) | set(dt.screen + "|" + dt.subpair)
            | {k for k, v in best.items() if k.startswith("Acylation|1|") and v >= 0.10})
    units = dict(A1=A[y.notna() & ~dup & ~process], A3=A[dup], SA1=SA[sk.isin(keep)],
                 SB1=pd.read_csv(DATA / "Validation_Set_B/never_opened_cross_library.csv.gz", dtype=str))
    for u, d in units.items():
        assert len(d) == UNITS[u], (u, len(d))
    P = pd.concat([d.assign(unit=u) for u, d in units.items()], ignore_index=True)
    assert len(P) == N["P2"] and pd.to_numeric(P["yield"], errors="coerce").notna().all()
    lib = P.rxn_class.map(lib_of)
    rows = T.reindex(list(zip(lib, P.condition_id)))
    ok = rows.names.notna().to_numpy()
    for a, b in [(f"reagent{i}_id", f"Reagent{i}_id") for i in range(1, 7)] + [("solvent_id", "solvent_id")]:
        ok &= P[a].fillna("").to_numpy() == rows[b].fillna("").to_numpy()
    for a in ("temp_id", "time_id"):
        u, v = pd.to_numeric(P[a], errors="coerce").to_numpy(), pd.to_numeric(rows[a], errors="coerce").to_numpy()
        ok &= (u == v) | (np.isnan(u) & np.isnan(v))
    P["lib"], P["label"], P["ok"] = lib.to_numpy(), P.condition_id, ok
    P["s1"], P["s2"] = P.sub_1_smiles.map(canon), P.sub_2_smiles.map(canon)
    return P


def companion(p1, T, lib_of, ROLE, screens_of):
    """P1: the companion dataset; conditions resolved by its own compositions, partners oriented per library."""
    W = pd.read_csv(p1 / "data/experimental_wells.csv", dtype=str)
    assert len(W) == N["P1"]
    comp = {}
    for r in pd.read_csv(p1 / "metadata/base_condition_roles.csv", dtype=str).itertuples():
        comp[(r.rxn_class, r.condition_id)] = (tuple(r.reagents.split("; ")) if isinstance(r.reagents, str) and r.reagents.strip() else (),
                                               r.solvent, num(r.temperature_C), num(r.time_h))
    for r in pd.read_csv(p1 / "metadata/extended_conditions.csv", dtype=str).to_dict("records"):
        comp[(r["rxn_class"], r["condition_id"])] = (tuple(x for x in (r[n] for n in NAMES) if isinstance(x, str) and x.strip()),
                                                     r["Standardized Name"], num(r["temp_id"]), num(r["time_id"]))
    lib, ok = [], []
    for rc, lab in zip(W.rxn_class, W.condition_id):
        L = lib_of.get(rc)
        lib.append(L)
        ok.append(L is not None and (L, lab) in T.index and comp.get((rc, lab)) == (
            T.loc[(L, lab)].names, T.loc[(L, lab)]["Standardized Name"], num(T.loc[(L, lab)].temp_id), num(T.loc[(L, lab)].time_id)))
    W["lib"], W["label"], W["ok"] = lib, W.condition_id, ok
    cE, cN = W.electrophile_smiles.map(canon), W.nucleophile_smiles.map(canon)
    first = {}                                     # library -> True if the electrophile is sub_1, False if sub_2
    for L, g in W[W.lib.notna()].groupby("lib"):
        if L in screens_of:
            e1 = min(share(g.electrophile_smiles, screens_of[L], ROLE, 0), share(g.nucleophile_smiles, screens_of[L], ROLE, 1))
            e2 = min(share(g.electrophile_smiles, screens_of[L], ROLE, 1), share(g.nucleophile_smiles, screens_of[L], ROLE, 0))
            if abs(e1 - e2) >= 0.5:
                first[L] = e1 > e2
    W["s1"] = [None if L not in first else (e if first[L] else n) for L, e, n in zip(W.lib, cE, cN)]
    W["s2"] = [None if L not in first else (n if first[L] else e) for L, e, n in zip(W.lib, cE, cN)]
    files = ("data/experimental_wells.csv", "metadata/base_condition_roles.csv", "metadata/extended_conditions.csv")
    return W, {f: sha(p1 / f) for f in files}


def p4_dataset(p4, lib_of):
    """P4: the numeric records of the yield-prediction study's released file, in file order, with their row_uid; only the screens
    that are a library of this release."""
    cols = ["screen", "row_uid", "condition_id", "sub_1_smiles", "sub_2_smiles", "yield_fraction"]
    D = pd.read_csv(p4 / P4_RECORDS, dtype=str, keep_default_na=False, usecols=cols)
    assert D.row_uid.is_unique, "a repeated row_uid in the P4 file"
    D = D[D.screen.isin(lib_of.keys()) & (D.yield_fraction != "")].reset_index(drop=True)
    D["lib"] = D.screen.map(lib_of)
    return D, {P4_RECORDS: sha(p4 / P4_RECORDS)}


def unordered_key(lib, label, s1, s2, y):
    return None if any(x is None or x != x for x in (lib, label, s1, s2, y)) else SEP.join((lib, label, *sorted((s1, s2)), repr(y)))


def p4_uids(P, D):
    """For each released record, the row_uid of its one-to-one partner in the P4 dataset, or ''."""
    dk = [unordered_key(l, c, canon(a), canon(b), yield_key(y)) for l, c, a, b, y in
          zip(D.lib, D.condition_id, D.sub_1_smiles, D.sub_2_smiles, D.yield_fraction)]
    pool = {}
    for k, u in zip(dk, D.row_uid):
        if k is not None:
            pool.setdefault(k, []).append(u)
    seen, out = Counter(), []
    for k in (unordered_key(l, c, a, b, yield_key(y)) for l, c, a, b, y in zip(P.lib, P.label, P.s1, P.s2, P["yield"])):
        u = pool.get(k, []) if k is not None else []
        out.append(u[seen[k]] if seen[k] < len(u) else "")
        seen[k] += 1
    return out


def zhang(path, T, lib_of, ROLE, screens_of):
    """Zhang: the commit's file; conditions resolved through the code bridge, partners in stored order."""
    raw = Path(path).read_bytes()
    raw = gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw
    if hashlib.sha256(raw).hexdigest() != ZHANG["sha256"]:
        raise SystemExit("STOP: the Zhang file is not the file of the pinned commit")
    Z = pd.read_csv(io.BytesIO(raw), dtype=str)
    assert len(Z) == N["Zhang"]
    L = lib_of["Acylation"]
    Z["label"] = "C" + Z.condition_id
    labs = [(r.label, r) for r in Z.drop_duplicates("condition_id").itertuples() if (L, r.label) in T.index]
    roles = ("Activation_ID", "Additive_ID", "Base_ID")
    cand = {}
    for lab, r in labs:
        names, solvent = set(T.loc[(L, lab)].names), T.loc[(L, lab)]["Standardized Name"]
        for k, s in [((c, getattr(r, c)), names) for c in roles] + [(("solvent_id", r.solvent_id), {solvent})]:
            cand[k] = cand[k] & s if k in cand else set(s)
    assign = {}
    while True:                                    # elimination within each component role; no candidate left means "none"
        changed = False
        for k, s in cand.items():
            if k not in assign and len(s) <= 1:
                assign[k], changed = (next(iter(s)) if s else None), True
        for k in cand:
            if k not in assign:
                cand[k] = cand[k] - {v for kk, v in assign.items() if v is not None and kk[0] == k[0]}
        if not changed:
            break

    def agrees(lab, r):
        names = list(T.loc[(L, lab)].names)
        opts = [{assign[k]} if k in assign else cand[k] for k in ((c, getattr(r, c)) for c in roles) if not (k in assign and assign[k] is None)]
        return (len(opts) == len(names) and any(all(n in o for n, o in zip(p, opts)) for p in itertools.permutations(names))
                and assign.get(("solvent_id", r.solvent_id)) == T.loc[(L, lab)]["Standardized Name"])
    agree = {lab for lab, r in labs if agrees(lab, r)}
    stored_order = share(Z.sub_1_smiles, screens_of[L], ROLE, 0) >= 0.95 and share(Z.sub_2_smiles, screens_of[L], ROLE, 1) >= 0.95
    Z["lib"], Z["ok"] = L, Z.label.isin(agree) & stored_order
    Z["s1"], Z["s2"] = Z.sub_1_smiles.map(canon), Z.sub_2_smiles.map(canon)
    return Z


# ------------------------------------------------------------------ counting
def regions(kP1, kP2, kZ):
    """The seven regions, key by key."""
    c1, c2, cz = (Counter(k for k in ks if k) for ks in (kP1, kP2, kZ))
    r = dict.fromkeys(REGIONS, 0)
    for k in set(c1) | set(c2) | set(cz):
        a, b, z = c1.get(k, 0), c2.get(k, 0), cz.get(k, 0)
        m3 = min(a, b, z)
        r["P1 ∩ P2 ∩ Zhang"] += m3
        r["P1 ∩ P2 only"] += min(a, b) - m3
        r["P1 ∩ Zhang only"] += min(a, z) - m3
        r["P2 ∩ Zhang only"] += min(b, z) - m3
        r["P1 only"] += max(0, a - max(b, z))
        r["P2 only"] += max(0, b - max(a, z))
        r["Zhang only"] += max(0, z - max(a, b))
    for name, ks in zip(("P1 only", "P2 only", "Zhang only"), (kP1, kP2, kZ)):
        r[name] += sum(1 for k in ks if not k)
    return r


def regions_by_position(kP1, kP2, kZ):
    """The same regions position by position: the j-th record of a key lies in every dataset holding at least j of them."""
    c1, c2, cz = (Counter(k for k in ks if k) for ks in (kP1, kP2, kZ))
    name = {(1, 0, 0): "P1 only", (0, 1, 0): "P2 only", (0, 0, 1): "Zhang only", (1, 1, 0): "P1 ∩ P2 only",
            (1, 0, 1): "P1 ∩ Zhang only", (0, 1, 1): "P2 ∩ Zhang only", (1, 1, 1): "P1 ∩ P2 ∩ Zhang"}
    r = dict.fromkeys(REGIONS, 0)
    for k in set(c1) | set(c2) | set(cz):
        a, b, z = c1.get(k, 0), c2.get(k, 0), cz.get(k, 0)
        for j in range(1, max(a, b, z) + 1):
            r[name[(int(j <= a), int(j <= b), int(j <= z))]] += 1
    for nm, ks in zip(("P1 only", "P2 only", "Zhang only"), (kP1, kP2, kZ)):
        r[nm] += sum(1 for k in ks if not k)
    return r


def flags(X, Y):
    """For each record of X, whether it is among the first n(Y) records of its key, X in released order."""
    rank = X.groupby("key", sort=False).cumcount()
    return (X.key.notna() & (rank < X.key.map(Y.key.dropna().value_counts()).fillna(0))).to_numpy()


def guard(datasets):
    """Records of two datasets joined by the two-stage rounding but not by a single rounding: stored yields <= 1e-9 apart."""
    A = pd.concat([pd.DataFrame({"ds": n, "k": D.key, "k1": D.key_single, "y": pd.to_numeric(D["yield"], errors="coerce")})
                   for n, D in datasets.items()], ignore_index=True).dropna(subset=["k"])
    worst, pairs = 0.0, 0
    for _, g in A[A.groupby("k").k1.transform("nunique") > 1].groupby("k"):
        for x, z in itertools.combinations(list(g.itertuples(index=False)), 2):
            if x.ds != z.ds and x.k1 != z.k1:
                pairs, worst = pairs + 1, max(worst, abs(x.y - z.y))
    if worst > 1e-9:
        raise SystemExit(f"STOP: the two-stage rounding joins stored yields {worst:.3g} apart")
    return dict(record_pairs_joined_by_the_two_stage_rounding=pairs, max_abs_stored_yield_difference=worst, threshold=1e-9)


def main(p1, p4, out=OUT, zhang_file=RAW / "public_amide_reference.csv.gz", conditions=RAW / "DATASET_conditions.csv.gz"):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    T, lib_of = condition_table(conditions)
    ROLE = role_table()
    screens_of = {}
    for s in ROLE:
        screens_of.setdefault(lib_of[s.split("|")[0]], []).append(s)
    P = released(T, lib_of)
    W, p1_hashes = companion(Path(p1), T, lib_of, ROLE, screens_of)
    Z = zhang(zhang_file, T, lib_of, ROLE, screens_of)
    for D in (P, W, Z):
        for col, two in (("key", True), ("key_single", False)):
            D[col] = [key(l, c, a, b, yield_key(t, two)) if ok else None
                      for l, c, a, b, t, ok in zip(D.lib, D.label, D.s1, D.s2, D["yield"], D.ok)]
    g = guard(dict(P2=P, P1=W, Zhang=Z))
    r = regions(list(W.key), list(P.key), list(Z.key))
    assert r == regions_by_position(list(W.key), list(P.key), list(Z.key)), "the two region computations disagree"
    for ds in ("P1", "P2", "Zhang"):
        assert sum(v for k, v in r.items() if ds in k) == N[ds], ds
    counts = {"P2": N["P2"], "P2 ∩ P1": r["P1 ∩ P2 only"] + r["P1 ∩ P2 ∩ Zhang"], "P2 ∩ Zhang": r["P2 ∩ Zhang only"] + r["P1 ∩ P2 ∩ Zhang"],
              "P2 ∩ P1 ∩ Zhang": r["P1 ∩ P2 ∩ Zhang"], "P2 ∩ (P1 ∪ Zhang)": N["P2"] - r["P2 only"], "P2 \\ (P1 ∪ Zhang)": r["P2 only"]}
    fc, fz = flags(P, W), flags(P, Z)
    cls = np.where(fc & fz, "both", np.where(fc, "p1_only", np.where(fz, "zhang_only", "neither")))
    t = Counter(cls)
    assert int(fc.sum()) == counts["P2 ∩ P1"] and int(fz.sum()) == counts["P2 ∩ Zhang"] and t["both"] == counts["P2 ∩ P1 ∩ Zhang"]
    assert N["P2"] - t["neither"] == counts["P2 ∩ (P1 ∪ Zhang)"] and t["neither"] == counts["P2 \\ (P1 ∪ Zhang)"]
    rec = P[["unit", "screen", "condition_id", "sub_1_smiles", "sub_2_smiles", "yield"]].copy()
    rec["matches_zhang2025"], rec["matches_companion_p1"], rec["overlap_class"], rec["match_rule_version"] = fz, fc, cls, VERSION
    D4, p4_hashes = p4_dataset(Path(p4), lib_of)
    rec["p4_row_uid"] = p4_uids(P, D4)
    n4 = int((rec.p4_row_uid != "").sum())
    assert rec.p4_row_uid[rec.p4_row_uid != ""].is_unique
    rec.to_csv(out / "P2_RECORD_PROVENANCE.csv", index=False, lineterminator="\n")
    inputs = ["05_DATA/Release/00_RAW_IMMUTABLE/DATASET_analysis.csv.gz", "05_DATA/Validation_Set_A/never_analysed_reserve.csv.gz",
              "05_DATA/Validation_Set_B/never_opened_cross_library.csv.gz", "05_DATA/Results/reserve_validation/record_transfer_losses.csv",
              "05_DATA/Results/reserve_validation/decision_targets.csv"]
    result = dict(
        rule=dict(version=VERSION, text=RULE.relative_to(ROOT).as_posix(), text_sha256=sha(RULE),
                  history="PROV-1 was frozen before counting; Amendment 2 was adopted after the first PROV-1 run to correct two "
                          "documented mechanical matching defects; the final counts use PROV-1+A2"),
        environment=dict(python=platform.python_version(), pandas=pd.__version__, numpy=np.__version__, rdkit=rdkit.__version__),
        inputs=dict(P2={f: sha(ROOT / f) for f in inputs},
                    condition_dictionary={"05_DATA/Release/00_RAW_IMMUTABLE/DATASET_conditions.csv.gz": sha(conditions)},
                    role_table={"06_CODE/Core/00_common/ctb/history.py": sha(ROOT / "06_CODE/Core/00_common/ctb/history.py")},
                    P1=p1_hashes, Zhang=ZHANG, P4=p4_hashes),
        records=N, venn_regions=r, p2_disclosure_counts=counts,
        p4=dict(identifier="p4_row_uid, not part of PROV-1+A2: library, condition label, unordered canonical partners, yield as "
                           "PROV-1+A2 rounds it; one to one in released order",
                p4_records_read=len(D4), p4_screens_read=sorted(D4.screen.unique().tolist()),
                p2_with_p4_row_uid=n4, p2_without=N["P2"] - n4),
        cross_check=dict(regions_key_by_key_equal_regions_position_by_position=True,
                         per_record_flags_reproduce_every_p2_count=True, yield_guard=g))
    (out / "P1_P2_ZHANG_OVERLAP_COUNTS_PUBLIC.json").write_text(json.dumps(result, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(dict(venn_regions=r, p2_disclosure_counts=counts, yield_guard=g, p2_with_p4_row_uid=n4), ensure_ascii=False, indent=1))
    return result


def without_companion(zhang_file, released_file=OUT / "P2_RECORD_PROVENANCE.csv", conditions=RAW / "DATASET_conditions.csv.gz"):
    """Without the companion package: the released records, keyed by the rule, against the released flags file; with the Zhang
    file, the Zhang flags recomputed. Writes nothing."""
    T, lib_of = condition_table(conditions)
    ROLE = role_table()
    screens_of = {}
    for s in ROLE:
        screens_of.setdefault(lib_of[s.split("|")[0]], []).append(s)
    P = released(T, lib_of)
    R = pd.read_csv(released_file, dtype=str)
    cols = ["unit", "screen", "condition_id", "sub_1_smiles", "sub_2_smiles", "yield"]
    rows_equal = len(R) == len(P) and bool((R[cols].fillna("").to_numpy() == P[cols].fillna("").to_numpy()).all())
    out = dict(records=len(P), rows_equal_to_released_flags_file=rows_equal, conditions_resolved=int(P.ok.sum()),
               companion="not supplied: matches_companion_p1 and the counts involving the companion dataset are not recomputed")
    zf = Path(zhang_file) if zhang_file else None
    if zf is not None and zf.is_file():
        Z = zhang(zf, T, lib_of, ROLE, screens_of)
        for D in (P, Z):
            D["key"] = [key(l, c, a, b, yield_key(t)) if ok else None for l, c, a, b, t, ok in zip(D.lib, D.label, D.s1, D.s2, D["yield"], D.ok)]
        fz = flags(P, Z)
        out.update(zhang_flags_equal=bool((R.matches_zhang2025.map({"True": True, "False": False}).to_numpy() == fz).all()),
                   p2_matching_zhang=int(fz.sum()))
    else:
        out["zhang"] = "not supplied: matches_zhang2025 is not recomputed"
    ok = rows_equal and out["conditions_resolved"] == N["P2"] and out.get("zhang_flags_equal", True)
    out["result"] = "PASS" if ok else "FAIL"
    print(json.dumps(out, ensure_ascii=False, indent=1))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--p1", help="folder of the companion manuscript's data package (optional; without it nothing is rewritten)")
    ap.add_argument("--p4", help="folder of the yield-prediction study's data package (required with --p1)")
    ap.add_argument("--zhang", default=str(RAW / "public_amide_reference.csv.gz"), help="the commit's all_HTE_with_condition.csv, plain or gzip")
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args()
    if a.p1:
        if not a.p4:
            raise SystemExit("--p1 rewrites the flags, which include p4_row_uid: give --p4 as well")
        main(a.p1, a.p4, a.out, a.zhang)
    else:
        raise SystemExit(0 if without_companion(a.zhang)["result"] == "PASS" else 1)
