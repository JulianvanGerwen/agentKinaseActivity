# DATA — provenance of every dataset used

Nothing under `data/`, `external/` or `cache/` is committed (all git-ignored). Checksums are sha256 of the files as received on 2026-10-27. Re-verify with `sha256sum data/* data/Klaeger/*`.

| File | Source | Version / date | sha256 | License | Commit? |
|---|---|---|---|---|---|
| `data/00_guideToFiles.txt` | provided by user | 2026-10-27 | `8ec5f995…594e4f` | unknown — provided by user | no |
| `data/01_phosphoproteomics_conditions.csv` | user-provided compendium of human phosphoproteomics (1937 conditions, 178 PMIDs). The `source` column says "This study" (1437), "qPTM" (461; qPTM database, Yu et al. 2023 NAR) or "Ochoa2016" (39; Ochoa et al. 2016 MSB) | unknown; received 2026-10-27 | `1257eb67…0af690` | unknown — provided by user (the underlying studies each have their own terms) | no |
| `data/02_phosphosite_information.csv` | same compendium, site metadata (peptide_UID to a 15-mer window) | as above | `56cab92d…ac8bfc` | unknown — provided by user | no |
| `data/03_phosphosite_quantifications.csv` | same compendium, long table of log2FC (17,032,346 rows) | as above | `8877144f…3198d7` | unknown — provided by user | no |
| `data/Kinase_Substrate_Dataset` | PhosphoSitePlus (Cell Signaling Technology), www.phosphosite.org | header date code `051925` (19 May 2025) | `bf057b8b…ed3c4c` | CC BY-NC-SA 3.0 (cite Hornbeck et al. 2015 NAR) | no |
| `data/Phosphorylation_site_dataset` | PhosphoSitePlus | `051925` | `9d79ab45…3d49e6` | CC BY-NC-SA 3.0 | no |
| `data/Regulatory_sites` | PhosphoSitePlus | `051925` | `1e654c31…c974c7` | CC BY-NC-SA 3.0 | no |
| `data/kinhib_inhib_comb__1.csv` | user-provided; inhibitor × kinase 0/1 matrix for Hijazi et al. 2020 Nat Biotechnol (PMID 31959955). Per the paper: DiscoverX >50% inhibition at 200 nM for 20 compounds, plus vendor/literature data (incl. Klaeger 2017) for 40 | unknown | `b0217e79…1ee545` | unknown — provided by user | no |
| `data/Klaeger/aan4368_Table_S1.xlsx` | Klaeger et al. 2017 Science 358:eaan4368, Table S1 (inhibitor annotation) | publication supplement (2017) | `f993057c…39f433` | journal supplementary material; terms per Science/AAAS — treat as unknown, do not redistribute | no |
| `data/Klaeger/aan4368_Table_S2.xlsx` | Klaeger 2017, Table S2 (kinobeads dose–response targets) | 2017 | `d28b91e6…db2c32` | as above | no |
| `data/Klaeger/Table_S3.xlsx` | Klaeger 2017, Table S3 (drug matrix of apparent Kd in nM) | 2017 | `57b399b0…b6ecf5` | as above | no |
| `data/Klaeger/desktop.ini` | Google Drive artefact, not data | – | `663a74a5…2b606f` | – | no |
| `data/Paper_Hijazi.pdf` | Hijazi et al. 2020 Nat Biotechnol 38:493 | publisher PDF | `a795c281…b61ee0` | publisher copyright | no |
| `data/Paper_KSEA.pdf` | Casado et al. 2013 Sci Signal (KSEA) | publisher PDF | `fb2a5291…0d8268` | publisher copyright | no |
| `data/Paper_JumpSEM.pdf` | JUMPsem paper | publisher PDF | `62d96d6a…731350` | publisher copyright | no |
| `external/hgnc/hgnc_complete_set.txt` | HGNC complete set, downloaded 2026-10-27 from `https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt` | latest `date_modified` 2026-10-15 | `12363916…2c2f91` | HGNC data are freely available (no restrictions) | no (re-downloadable) |

Full checksums:
```
8ec5f995235877e85f0092caedb80d0b77f9f9707802be2ab50885f4eb594e4f  data/00_guideToFiles.txt
1257eb673d1b6db532773f4d6a54a75f2960a7f4434ccb67b0239834170af690  data/01_phosphoproteomics_conditions.csv
56cab92d90ea1760302169674552161bee405d946ff32c8221935f483dac8bfc  data/02_phosphosite_information.csv
8877144f25cc95b4625e9c55e05d33ffe2ef11af07b1297bf27791d13d3198d7  data/03_phosphosite_quantifications.csv
bf057b8bb347881f8d0ac82d705adcc0be6f346f8c72882078fbd440a2ed3c4c  data/Kinase_Substrate_Dataset
a795c281c909114a354146c44f802aef7142fba6da468293b52acf7b67b61ee0  data/Paper_Hijazi.pdf
62d96d6a049d878ac75cf51a8e037fc7840af70eb6677876b561588aee731350  data/Paper_JumpSEM.pdf
fb2a5291cf09cc3123155004a53742f9afda1bd39ec86b0fc416ece2620d8268  data/Paper_KSEA.pdf
9d79ab45cbcf9f4813db1afcfe0520f3202e579b75e119c8455ef930443d49e6  data/Phosphorylation_site_dataset
1e654c316bee94dcf2713544796ed2dfbfdcb1ebd42b72db89166bf6d6c974c7  data/Regulatory_sites
b0217e7962ba936f18a7e4d060d849d10d0b38d2f0e04d3aa730d0c0861ee545  data/kinhib_inhib_comb__1.csv
57b399b0855c911b04f4060976873e61ad17421555a6df81b655458b6cb6ecf5  data/Klaeger/Table_S3.xlsx
f993057ccf9e41621614ba9ee8acede93fc490da8cf741a770f3977c2739f433  data/Klaeger/aan4368_Table_S1.xlsx
d28b91e62e78e5e011b60da27672875621fef5cdabbea793ac9cce4b98db2c32  data/Klaeger/aan4368_Table_S2.xlsx
663a74a582a69d0dba5f4086981374e4aa8bf6d01dcab479c2e9ff0c8f2b606f  data/Klaeger/desktop.ini
12363916be36f8f6f91879b0faade4eeaac973722de7136ec61329b77f2c2f91  external/hgnc/hgnc_complete_set.txt
```

Notes
- PSP files are latin-1 encoded and start with 3 header lines (date code, license, blank). Parse them with `scripts/prep_data.py`.
- Derived caches are in `cache/` (regenerate with `python scripts/prep_data.py`, about 30 s). Small derived tables are in `results/` and may be committed. Some of them hold PSP-derived rows (for example `results/data_profile/psp_kinase_coverage.tsv` holds counts only). Do not commit raw PSP rows: CC BY-NC-SA requires attribution and share-alike.
- Manual mappings are in `results/ground_truth/drug_targets_manual.tsv` and `condition_rules_manual.tsv`, each with a note column.
