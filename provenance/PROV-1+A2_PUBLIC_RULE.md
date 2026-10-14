# Provenance matching rule PROV-1+A2

This is a versioned rule. PROV-1 was frozen before any count was computed. Amendment 1 was made before the first count, and Amendment 2 was adopted after the first PROV-1 run (section 6). The final counts use PROV-1+A2.

The rule decides which of the 244,436 measurements released with this manuscript exactly match records of two other datasets: the dataset analysed in a companion manuscript, and the published amide-coupling collection of Zhang *et al.* (2025). The counts serve disclosure only. No analysis, figure or release scope depends on them.

`06_CODE/Release/provenance_overlap_public.py` implements this rule and nothing else. Its results are the two files next to this one: `P1_P2_ZHANG_OVERLAP_COUNTS_PUBLIC.json` (counts) and `P2_RECORD_PROVENANCE.csv` (one row per released measurement).

## 1. Datasets

| Name | Records | Content |
|---|---:|---|
| P2 | 244,436 | The released measurements of this manuscript, in four release units: A1, 177,308 (the fifteen arrays without the duplicated Hiyama condition); A3, 74 (the duplicated Hiyama condition); SA1, 36,498 (the 420 released pairs of Set A, the first out-of-array evaluation set); SB1, 30,556 (the 535 pairs of Set B, the second). |
| P1 | 165,131 | The dataset of the companion manuscript: `data/experimental_wells.csv` of its data package. Condition compositions come from `metadata/base_condition_roles.csv` and `metadata/extended_conditions.csv`. The hashes of the three files are recorded with the counts. |
| Zhang | 47,015 | `data/all_HTE_with_condition.csv` of github.com/aichemeco/amide_coupling at commit `f6fdbe283fc24af53389aab8dcba5e46386ac9c0`: Git blob `778db5176ed6a1b3b11cbb80da4eacbd40ea74c5`, 17,657,403 bytes, SHA-256 `5085e1ae50ee26173bd170a4b10a99c4e3101a46b706d115e425aff17c091b41`. The script stops unless the file it reads, decompressed if gzip-compressed, has this SHA-256. |

The condition dictionary is the condition table of the release, 1,192 conditions, whose reaction types are the libraries below. The project copy, `05_DATA/Release/00_RAW_IMMUTABLE/DATASET_conditions.csv.gz`, also holds conditions used only by unreleased records; with or without them, every count and flag is the same. The roles of the two partners in each screen are given by the table `ROLE` in `06_CODE/Core/00_common/ctb/history.py`.

## 2. Record key

A record's key has five parts: library, condition label, partner 1, partner 2 and yield. Two records match only when all five are equal. A record without a key matches nothing.

1. **Library.** The reaction type (`reaction_type` of the condition dictionary) to which the record's condition resolves under section 3.
2. **Condition label.** The record's condition label within that library, under section 3.
3. **Partners.** The two partners as RDKit canonical isomeric SMILES, `Chem.MolToSmiles(Chem.MolFromSmiles(s), isomericSmiles=True)`, placed in the library's role order under section 4. A SMILES that RDKit cannot parse gives no key. The RDKit version is recorded with the counts.
4. **Yield.** The stored yield, a fraction of 1 in all three datasets, parsed as a floating-point number. It is rounded with Python `round` to nine decimal places and then to six, and negative zero is written as zero.

Rounding to nine places first joins floating-point representations of one value that would fall on opposite sides of a six-place boundary, such as 0.2700134999999999 and 0.2700135.

**Guard.** Take every pair of records from two datasets whose keys agree under this rule but differ under a single rounding to six places. Their stored yields must differ by at most 10⁻⁹, or the run stops.

## 3. Condition resolution

- **P2.** The library is the reaction type paired with the record's `rxn_class` in the condition dictionary, and the label is `condition_id`. The record's own reagent slots 1–6, `solvent_id`, `temp_id` and `time_id` must equal that dictionary row; empty equals empty, and numbers are compared as numbers. Otherwise the record has no key.
- **P1.** The library is the reaction type of the well's `rxn_class`, and the label is `condition_id`; the label must exist in that library. The companion package's own composition for (`rxn_class`, label) must equal the dictionary row: the reagent names in slot order, the solvent name, the temperature and the time.
  - Base wells take these from `base_condition_roles.csv` (`reagents`, `solvent`, `temperature_C`, `time_h`).
  - The second Buchwald–Hartwig campaign takes them from `extended_conditions.csv` (the reagent name slots, `Standardized Name`, `temp_id`, `time_id`).
  - Otherwise the well has no key.
- **Zhang.** The library is the reaction type of `rxn_class` "Acylation", and the label is `'C'` + `condition_id`. The label must exist in that library, and its composition must agree through the code bridge below. Otherwise the record has no key.

**Code bridge.** Zhang records name their components by codes (`Activation_ID`, `Additive_ID`, `Base_ID`, `solvent_id`). Each code is mapped to a dictionary name by elimination over the labels present in both datasets:

1. A code's candidates are the intersection, over every label that uses it, of that label's dictionary names. These are the reagent names for activation, additive and base codes, and the solvent name for solvent codes.
2. Repeat until nothing changes:
   - a code with exactly one candidate is assigned it, and a code with none is assigned "none";
   - each assigned name is then removed from the candidates of the unassigned codes of the same component role (activation, additive, base or solvent).
3. A label agrees when two conditions hold:
   - its codes that are not "none" can be paired one to one with its dictionary reagent names, each code with a name among its candidates. Codes left with several candidates are allowed whenever such a pairing exists;
   - its solvent code is assigned the dictionary solvent.

## 4. Role order

The key places the partners in the library's order (`sub_1`, `sub_2`), whose roles in each screen are given by `ROLE`.

- **P2.** Records are stored in this order.
- **Zhang.** Records are taken in stored order (`sub_1_smiles`, `sub_2_smiles`). This holds only if two shares both reach 95%: the distinct `sub_1` molecules matching the role-1 SMARTS of the amide library, and the distinct `sub_2` molecules matching its role-2 SMARTS. Otherwise no Zhang record has a key.
- **P1.** The companion dataset stores an electrophile and a nucleophile, and their order is decided per library:
  1. Each of the two assignments of (electrophile, nucleophile) to (`sub_1`, `sub_2`) gets a score. For each side, take the share of its distinct molecules that match the `ROLE` SMARTS of the role assigned to it, under any screen of the library. The score is the smaller of the two shares.
  2. The higher-scoring assignment is used if it exceeds the other by at least 0.5.
  3. Otherwise, and in libraries without a `ROLE` entry, the library's wells have no key.

## 5. Counting: one to one

For each key, let n(P2), n(P1) and n(Zhang) be the numbers of records carrying it. Records are paired one to one, never by set membership: the j-th record of a key in one dataset is paired with the j-th record of that key in every other dataset holding at least j of them.

Each position j therefore falls in exactly one of seven regions: one of the three datasets only, one of the three pairs of datasets only, or all three. A record without a key falls in its own dataset's single region. Every count is a sum of regions:

| Count | Regions | Formula |
|---|---|---|
| P2 ∩ P1 | P1 ∩ P2 only + all three | Σ min(n(P2), n(P1)) |
| P2 ∩ Zhang | P2 ∩ Zhang only + all three | Σ min(n(P2), n(Zhang)) |
| P2 ∩ P1 ∩ Zhang | all three | Σ min(n(P2), n(P1), n(Zhang)) |
| P2 ∩ (P1 ∪ Zhang) | every P2 region except P2 only | Σ min(n(P2), max(n(P1), n(Zhang))) |
| P2 \ (P1 ∪ Zhang) | P2 only | \|P2\| − \|P2 ∩ (P1 ∪ Zhang)\| |

**Per-record flags.** Within a key, P2 records are taken in release order: by unit (A1, A3, SA1, SB1), then by row.
- The first min(n(P2), n(P1)) records of a key are flagged `matches_companion_p1`.
- The first min(n(P2), n(Zhang)) records of a key are flagged `matches_zhang2025`.
- `overlap_class` is `both`, `p1_only`, `zhang_only` or `neither`, and `match_rule_version` is `PROV-1+A2`.

The flags reproduce every P2 count above.

A match means that two records are identical under this rule. It does not establish that they come from one physical experiment, and a record without a match is not thereby a first-time experiment.

## 6. Version history

- **PROV-1** was frozen on 2026-10-05, before any count was computed. Before freezing, only the datasets' formats and condition dictionaries had been inspected; no yields were compared between datasets.
- **Amendment 1** was made before the first count. A check of the rule's mechanics, which compared no yields, found two defects:
  - every code-bridge candidate set contained "none", so no reagent code could reach a single candidate;
  - a floor of 0.95 on the orientation score rejected the companion Suzuki library (scores 0.886 and 0.000), whose orientation is not in doubt.

  Sections 3 and 4 give the amended text.
- **Amendment 2** was adopted after the first PROV-1 run. It is a post-first-run correction of two documented mechanical matching defects; it was not specified before counting.
  - **Cross-role elimination.** Elimination in the code bridge crossed component roles. DMAP is an additive in some amide conditions and the base in C46 (BOPCl; DMAP), so both C46 codes were mapped to BOPCl, and the 462 Zhang records under C46 had no key. Yet the collection's own condition SMILES for C46 lists BOPCl, DMAP and DMF, exactly as the dictionary does. Elimination now stays within a component role (section 3).
  - **Rounding boundary.** A single rounding to six places separated floating-point representations of one value. The yield is now rounded to nine places first, under the guard of section 2.
- **Internal record.** It is available to editors and reviewers on request and comprises:
  - the PROV-1 text, SHA-256 `7ca45c28de481539776bc28f48fa564e58d731a14d0c95fa724982e99a599c1d`;
  - the rule with Amendments 1 and 2, SHA-256 `0907e12e838e053a3cc7555ffb2cd7fc308ad5b72bbc24406abc3153dc1d3e7b`;
  - the full counts record, SHA-256 `d4348e20a8c77072454f0538fb9a3e212e409c24cd0abb86555bea5191b52b25`, including the PROV-1 counts.
