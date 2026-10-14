# Data and code for "Reaction-condition knowledge has a molecular owner"

This package holds the measurements analysed in the article, the files that fix how they were analysed (frozen analysis plans,
and the folds, draws and groups the analyses used), the code that produced the figures, and the figures and Source Data
themselves. The figures, Extended Data figures and Source Data can be rebuilt from this package alone.

## Contents

| Folder or file | What it holds |
|---|---|
| `data/` | the 244,436 released yield measurements in four files; the design mask; the condition table; the reagent table |
| `metadata/` | the outcome-free manifest of Set A and its release record; the frozen same-molecule cohorts (`DESIGN_COUNTS.csv`); the publication use of every released pair |
| `provenance/` | record-level provenance of the released measurements: matches under the versioned matching rule PROV-1+A2 and each record's identifier in the dataset it is cited from; the rule, the counts and the matching code |
| `frozen_plans/` | the frozen analysis plans behind the reported results, with their SHA-256 records (`frozen_plans/PUBLIC_FROZEN_PLAN_README.md`), and beside them the provenance correction of the two out-of-array evaluation sets with its per-pair Set B flags |
| `derived/` | the analysis tables made from `data/` (task construction), with `DERIVED_FILES.json` |
| `results/` | the retained results the figures are checked against, with `RESULTS_MANIFEST.csv`; `TEXT_CLAIMS.csv`, the statements of the text recomputed by the rebuild |
| `source_data/` | the Source Data of Figures 1-4 and Extended Data Figs 1-10, with the released folds, draws and groups |
| `figures/` | Figures 1-4 and Extended Data Figs 1-10 (PDF, SVG, PNG at 1,200 dpi); `structures/`, the ChemDraw structures of Figs 1 and 2 |
| `code/` | the analysis and figure code; `code/release/`, the code of the two entry points below |
| `environment/` | the Python and package versions |
| `rebuild_all.py`, `verify_release.py` | the two entry points |
| `DATA_MANIFEST.csv`, `CHECKSUMS.sha256` | every file of the package with its size, SHA-256 and role |
| `LICENSE_DATA_CC-BY-4.0.txt`, `LICENSE_CODE_MIT.txt` | the licences: CC BY 4.0 for the data and documents, MIT for the code |
| `.gitattributes` | Git settings: no line-ending conversion, so a clone holds the bytes listed in `CHECKSUMS.sha256` |

## The measurements: 244,436 numeric yields in four files

| Unit | File | Records | Pairs | What they are |
|---|---|---|---|---|
| A1 | `data/measurements_A1_analysis_arrays.csv.gz` | 177,308 | 2,210 | the fifteen analysis arrays (the second copy of the duplicated Hiyama condition is in A3) |
| A3 | `data/measurements_A3_hiyama_duplicate.csv.gz` | 74 | 74 | Hiyama condition C46, measured as a duplicate of C43; used for repeatability and noise calibration |
| SA1 | `data/measurements_SA1_set_A.csv.gz` | 36,498 | 420 | the pairs of Set A, the first out-of-array evaluation set, that enter the record-transfer test (309), the three-condition decision (397) or the ownership re-estimate (217 product-forming amide pairs) |
| SB1 | `data/measurements_SB1_set_B.csv.gz` | 30,556 | 535 | every pair of Set B, the second out-of-array evaluation set |
| | | **244,436** | | |

Each record has a `public_record_id` of the form `P2-<unit>-<seven-digit number>`, numbered in file order and unique across
the four files. The Hiyama duplicate is a separate measurement and keeps its own identifier; records with the same chemistry are never
merged. The other columns are the unit, the array (`screen`), the reaction class, the plate (array measurements only), the
condition identifier, the substrate pair and its two substrates as SMILES, the yield as recorded (a fraction between 0 and
1) and the inherited flag `check_passed` (90 records carry `False`; the flag imposes no exclusion).

`data/design_mask_no_yield.csv` lists the 384 designed cells of the fifteen arrays that carry no numeric yield. They are
design information, not measurements: they have no identifier and are not among the 244,436.

`data/conditions.csv` gives the 1,192 conditions named by the records with their full compositions (reagents with CAS
number, equivalents, name and identifier; solvent; temperature; time); a condition is identified by `rxn_class` and
`condition_id`. `data/reagents.csv` gives the 401 reagent entries referenced in the reagent slots of those conditions;
solvents are named in the condition table itself.

For Set A, `metadata/set_A_frozen_manifest.csv` lists all 442 frozen pairs without outcomes, and
`metadata/set_A_release.json` gives the SHA-256 of the frozen data file, the rule that selects the released pairs and the
SHA-256 of the released subset.

Set A and Set B are the article's two out-of-array evaluation sets: pairs outside the fifteen analysis arrays, excluded from
model fitting and scored under frozen complete-set plans. `frozen_plans/OUT_OF_ARRAY_EVALUATION_PROVENANCE_CORRECTION.md`
documents the earlier outcome access to parts of both sets in the direct project lineage and corrects the historical plans;
`frozen_plans/SET_B_DOCUMENTED_PRIOR_SCORING.csv` flags the 124 Set B pairs concerned.

## Provenance of the records

`provenance/P2_RECORD_PROVENANCE.csv` has one row per released measurement, in the order of the four files (unit A1 first),
so its row for record *i* of a unit is `P2-<unit>-<i>`. For each record it says whether the record matches, under the
versioned rule PROV-1+A2 (`provenance/PROV-1+A2_PUBLIC_RULE.md`), a record of the previously published amide-coupling
collection of Zhang et al. and a record of the dataset released with a companion manuscript; the totals are in
`provenance/P1_P2_ZHANG_OVERLAP_COUNTS_PUBLIC.json`. Every released measurement is a record of the dataset released with a study
of yield prediction across reaction screens (ref. 31 of the manuscript); the last column, `p4_row_uid`, gives the record's
`row_uid` in that dataset. That column is matched on the library, the condition label, the two molecules as an unordered pair
and the yield, one to one (`provenance/provenance_overlap_public.py`).

A match is an exact record match after harmonising the reaction library, the molecular identities, the condition identity and
the yield representation. It does not establish that two records come from the same physical experiment. No first disclosure is
claimed for any released measurement.

Neither compared dataset is redistributed here. The Zhang et al. collection is pinned by repository, commit, path and hash in
the rule and the counts. The companion manuscript's dataset is released with that manuscript and is not copied into this
package. To re-derive the flags, give both:

    python rebuild_all.py --p1 <folder of the companion data package> --zhang <all_HTE_with_condition.csv of the pinned commit>

Without them, the rebuild checks the released flags against the released records only, and rebuilds the figures all the same.

## Verify only, or rebuild

    python verify_release.py

checks the package as it stands, without rebuilding or writing anything: the hashes of every file, the 244,436 records and
their identifiers, the design mask, the condition and reagent tables, the release rule of Set A, the
alignment of the provenance flags, the derived tables, the results manifest, the publication use of every measurement, the
statements of the text and the frozen-plan records.

    python rebuild_all.py

rebuilds everything from this package (about 15 minutes with 6 processes on the computer used for the study; `--jobs` sets the
number). It makes the derived
tables again from `data/`, runs the fourteen figure generators, the re-fit of Extended Data Fig. 10, the owner labels fixed without the scored pairs
(`code/Checks/nested_owner.py`) the owner margins, equal-role null and SNAr class proxy (`code/Checks/role_null_checks.py`) and the per-array tables of the
Supplementary Information (`code/Checks/si_tables_checks.py`), recomputes the
statements of the text and the publication use of every measurement, and compares every rebuilt figure, Source Data file and
provenance record with the released one. Every Python process of the rebuild runs under a file-access guard: it may read this
package, the Python installation and the system fonts, and write only in `_rebuild/`; any other access stops the rebuild. The
outcome is written to `_rebuild/REBUILD_REPORT.json`; `_rebuild/` is not part of the release.

On the reference environment (`environment/`) every rebuilt figure and Source Data file is byte-identical to the released
one. A provenance record (`*_SourceData_provenance.json`) is identical apart from its list of inputs, which names the files
the rebuild read: the derived tables and filtered results of this package where the original build read internal files, and
the released folds and draws where it read the internal fold registry. `rebuild_all.py` accepts exactly these differences
and no other.

## How the code finds its inputs

The code under `code/` is the code that produced the figures, unchanged, and it addresses its inputs by the paths of the
original working layout. `rebuild_all.py` therefore builds a work tree in `_rebuild/tree/` from copies of the package files
at those paths (`code/release/package_lib.py`, `stage_map`), and runs the generators in public mode
(`P2_RELEASE_MODE=public`), in which the hidden state (the folds of Fig. 4a and Extended Data Figs 7 and 8, the draws of
Fig. 4a and Extended Data Fig. 7, the molecule groups, condition folds and pair halves of Extended Data Fig. 10) is read from
the released Source Data instead of the internal fold registry. The shared modules also hold functions for analyses that are
not part of the article; the rebuild does not run them.

## Licences

The source code is released under the MIT licence (`LICENSE_CODE_MIT.txt`, Copyright (c) 2026 Kuangbiao Liao): `code/`,
`rebuild_all.py`, `verify_release.py`, `provenance/provenance_overlap_public.py` and the frozen scoring code in
`frozen_plans/` (`reserve_validation/reserve_validation.py`, `cross_library/cross_library_validation.py`). Everything else (the
measurements, tables, metadata, provenance records, frozen plans and their correction, derived tables, results, Source Data,
figures, structures and documentation) is released under the Creative Commons Attribution 4.0 International licence
(`LICENSE_DATA_CC-BY-4.0.txt`). Neither licence covers the amide-coupling collection of Zhang et al. or the dataset of the
companion manuscript; neither is redistributed here.

## Fonts

No font file is included: Arial is not redistributed. `environment/fonts.json` records the versions and SHA-256 of the Arial
files the released figures were drawn with (regular, Version 6.90; bold and italic, Version 7.00). The numeric results and the
Source Data do not depend on fonts; byte-identical figure files (PDF, SVG, PNG) require the recorded Arial files.
