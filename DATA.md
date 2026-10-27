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

## External benchmark data (copied 2026-10-27 from kinaseSubMain, never committed; checksums in results/external_checksums.txt)
- `external/benchmarkin/` — benchmarKIN Zenodo 10.5281/zenodo.12566560 (Müller-Dott et al. 2025, Nat Commun 16:4771), CC BY 4.0. Originally downloaded 2026-10-26 by kinaseSubMain (see its DATA.md for per-file provenance).
- `external/benchmarkin_github/` — tyrosineData.rda / tyrosineMeta.rda from github.com/saezlab/benchmarKIN main@2aa6e25c4e5c, GPL-3.
- `external/benchmarkin_py/` — Python (parquet/json) conversions of the above produced by kinaseSubMain `src/ksp/downstream/bench_data.py` (derived; same licences). The `r2_clean`/`r2_nophospho` dirs were NOT copied.

## KS-library sources (added 2026-10-27 by the library builder; copied READ-ONLY from kinaseSubMain = `KSM/` = /cluster/project/beltrao/JvG/Projects/hackathonKinaseSubs/kinaseSubMain; never committed)
Used only through `scripts/build_libraries.py` → `cache/libraries/*.parquet`; report `results/libraries/LIBRARY_REPORT.md`; per-build input hashes in `results/libraries/build_manifest.json`.

| File | Source | Version / date | sha256 | Licence |
|---|---|---|---|---|
| `external/omnipath/omnipath_enzsub.tsv` | OmniPath enzyme–substrate web service (`https://omnipathdb.org/enzsub?genesymbols=yes&fields=sources,references&format=tsv&license=academic`), downloaded by KSM 2026-10-26 | OmniPath web service 2026-10-26 | `57fa74706ef9c54ad23fe12e7390ba11a0694a8168534416d183fedc070e4f5d` | per-resource; academic licence set (includes non-commercial resources: PhosphoSite CC BY-NC-SA 3.0, HPRD academic, SIGNOR CC BY-SA 4.0, phospho.ELM academic, PhosphoNetworks, Li2012, KEA, dbPTM; MIMP/ProtMapper not used). Non-commercial research use only. |
| `external/iptmnet/{ptm,score,protein,readme}.txt` | iPTMnet 6.2 bulk download `https://research.bioinformatics.udel.edu/iptmnet_data/files/current/` (KSM, 2026-10-26) | release 6.2 | ptm.txt `d5b2ed7c7138bc77f6f712feb038e18ddc4c791e96823b286a680e1982eff20f`; score.txt `e7dfa494099ab9d8cdff54c727aba65251067a041b931846d8656c622563fe63`; protein.txt `381b90d7a9d5b1653c33bc3076db62ba9ffb3130c89f453276f310b2413e9309`; readme `03526d7031f5270804b26605a67325a6e0a979d30fc62fbf6b896e0e2760ef2c` | **Resolved: CC BY-NC-SA 4.0.** The readme (also the live copy, checked 2026-10-27) states the download is provided under CC BY-NC-SA 4.0 and, in a footer, CC BY 4.0. We apply the more restrictive NC-SA 4.0 (same non-commercial/share-alike terms as PSP, which our libraries already carry); cite iPTMnet and the original source DBs; derived libraries not committed. |
| `external/sugiyama2019/supp_PMC6642169.zip` | Sugiyama et al. 2019 Sci Rep 9:10503, "Large-scale discovery of substrates of the human kinome", PMID 31324866, doi 10.1038/s41598-019-46385-4; Europe PMC supplementary files (KSM, 2026-10-26). MOESM3 = Table S2 (198,536 in-vitro kinase–substrate–site rows) | journal suppl. 2019 | `178aee0f5d6ed76d044083d23c8a5f2d40611655979be325c510b690a178ba88` | CC BY 4.0 (open-access article) |
| `external/ksm_index/{kinase_list.tsv,kinase_synonyms.tsv,kinases.parquet,sites.parquet,universe.parquet}` | KSM `results/labels/` and `cache/index/`, `cache/sites/universe.parquet` (kinase universe U = include==True, kinase name resolver, score-matrix row/column index, site universe with compendium mapping `comp_uid`) | KSM 2026-10-26 | kinase_list `12fd39e611e0d335f08afea083156d3f62bdbe201fdd22991b0b6d356ddd857b`; kinase_synonyms `f198604c2ff5ff5e4b7bab15b506f29dd07cc622ea4246f282dc0e6cfd738d31`; kinases `e278ec9b4147d04dd33911f27b76e54e5528c1b1125c5140c588a1d1c49f4637`; sites `a280d3bb1a843c965b148b2067715567c9a515e9f3e18ab486bb97ed0d052975`; universe `17f9c6d38f5fbedd307c4a6ffad8e3a5d5cf022e5e12f7f8b9442657a5eb0ce7` | derived by KSM from UniProt (CC BY 4.0), HGNC, Manning/KinHub, PSP (CC BY-NC-SA 3.0) |
| `external/uniprot/UP000005640_{annotations.tsv.gz,reviewed_isoforms.fasta.gz}`, `MANIFEST.json` | UniProt human reference proteome (KSM `external/uniprot/`) — entry-name→accession and residue check for Sugiyama | UniProt 2026_03 | annotations `bb1ab565fc7b3cbee92cdf271f36980d5483a14fe172d8ccd5faa4f7d8c0c6dd`; fasta `63011ee75fbfd0c92f35412fb660cfe920980b0c96276f3f941de91039e8ece1` | CC BY 4.0 |
| `external/ksm_pred/r2np/uobs_predictions.parquet` (+ `manifest.json` `53230cef…`) | KSM R2 **no-phospho**, DEV-only labels, seed 0: `KSM/cache/runs/r2_final_nophospho_uobs__20261026-213950__s0/uobs_predictions.parquet` (groups F1,F0,F2a,F3,F4,F5,F7; no F6) | 2026-10-26 | `1edbd2691e3c2d638f449e1a423e97088c6d1299fc791a1eaf10f5d2aca2c052` | KSM internal (derived from PSP CC BY-NC-SA 3.0 labels) — non-commercial |
| `external/ksm_pred/r3/uobs_predictions.parquet` (+ `README.json` `ae5403ce…`) | KSM R3 deploy (5 seeds, DEV+TEST PSP labels, ESM-2 windows; no perturbation features): `KSM/cache/release/r3_deploy_20261027-005816/uobs_predictions.parquet` | 2026-10-27 | `2f62988cf798f15ace6cd9925dce77e45ca2c6cf297549203734d256e9d002ba` | as above |
| `external/ksm_pred/cantley/cantley_unprimed_score.npy` (+ `.json` `05c82894…`, `cantley_meta.json` `6a81dddb…`) | Kinase Library / Cantley atlas (Johnson 2023 Nature; Yaron-Barir 2024 Nature) unprimed published-definition scores, reproduced by KSM from kinase-library 1.8.0 matrices: `KSM/cache/features/motif/cantley_unprimed_score.npy` (503 × 350,736 float32; label_derived=False, not imputed) | kinase-library 1.8.0 | `e5597327729ab3605ca1649d510d89fbeffb0af1c416cd1fe96c6527bd4da2dc` | CC BY-NC-SA 3.0 (kinase-library package) |
| `cache/pubmed/pubmed_meta.jsonl` | PubMed E-utilities efetch (title, abstract, MeSH) for all PMIDs cited by OmniPath/iPTMnet phospho edges; fetched once by `kinact.library.pubmed_metadata` on 2026-10-27; drives `results/libraries/excluded_pmids.tsv` | 2026-10-27 | regenerated | NLM terms (public) |

Verified (`kinact.library.verify_predicted_file`, tests/test_libraries.py): the three predicted files match the audited sha256 and none matches the banned KSM outputs (release_scores `cceb32c7…`, r2_deploy `10e25be1…`, r2r3_release outputs `98f348ed…`/`61e20f51…`, clean_bench `ff5bfe80…`).
**Not used / dropped:** KSM `external/ikipdb_zenodo_{20037245,21469258}` ("iKiP-DB": cell-line atlas Cscores + Table7 iKiP∪PSP gmt; mislabelled, provenance inseparable — r2_data N3); SELPHI2, PTMsigDB, KSM release/clean_bench (banned, DESIGN §2); Phosformer (B17); NetworKIN predictions; OmniPath MIMP/ProtMapper-only edges; iPTMnet RLIMS-P text mining.

## Ground-truth inputs (copied/downloaded 2026-10-27 by the GT builder; never committed)
| File | Source | Version | sha256 | License |
|---|---|---|---|---|
| `external/kinasesubmain/kinase_list.tsv` | copied from kinaseSubMain `results/labels/kinase_list.tsv` (READ-ONLY sister project; 503 include==True kinases, Manning group/family) | KSM 2026-10-26 | `12fd39e611e0d335f08afea083156d3f62bdbe201fdd22991b0b6d356ddd857b` | derived from UniProt (CC BY 4.0), KinHub/Manning 2002 (academic use) |
| `external/reinecke2024/41589_2023_1459_MOESM2_ESM.xlsx` | Reinecke et al. 2024 Nat Chem Biol 20:577 suppl. Table S1 (compound annotation); copied from kinaseSubMain `external/reinecke2024/` (originally https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41589-023-01459-3/MediaObjects/) | 2024 suppl. | `9c13b2c95e2cd31de46b74fc9ecc25a2f3c72e0b101d0720bd6e917c3cecf5bb` | publisher supplementary; cite; do not redistribute |
| `external/reinecke2024/41589_2023_1459_MOESM3_ESM.xlsx` | Reinecke 2024 kinobeads drug matrix | 2024 suppl. | `aa62105218abe10c6c1eb3584a13a87769c4dd1d96cd9cc496cf6c2743d36606` | as above |
| `external/davis2011/41587_2011_BFnbt1990_MOESM{3,4,5}_ESM.xls` | Davis et al. 2011 Nat Biotechnol 29:1046 suppl. (constructs, compound synonyms, Kd matrix); copied from kinaseSubMain `external/davis2011/` | 2011 suppl. | MOESM3 `69c2957c1fa1b5ba032721031af2fffe1eb896acd0723829e4500da8f75779c9`; MOESM4 `f8b01f6bfa47adb61b6c5a6de2ae81281990da89b703848efec4f2e37488fae8`; MOESM5 `e47821b2e959ff22aa6d1613b236a63fb530f0da19cfa8d45408598cf27498d7` | publisher supplementary; cite; do not redistribute |
| `external/hijazi2020_si/41587_2019_391_MOESM3_ESM.xlsm` | Hijazi et al. 2020 Nat Biotechnol 38:493, Supplementary Dataset 1 ("Kinase-inhibitor selectivity": DiscoverX % control + Klaeger Kd); downloaded from https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41587-019-0391-9/MediaObjects/41587_2019_391_MOESM3_ESM.xlsm | 2020 suppl. | `e5931daafa4ef9a7b7ba26cdb1b96bddea62f14145f1c83b59b40ef3e2c6ea8f` | publisher supplementary; cite; do not redistribute |

Notes: the Hijazi SI (Datasets 1–2, SI PDF) and PRIDE PXD015943 (raw/mzid only) contain no per-drug concentrations,
so Hijazi doses use the pre-registered 10 µM default (paper: 1–10 µM, 1 h). Klaeger 2017 is read from `data/Klaeger/`.
The "kinoplex" files in kinaseSubMain (`external/kinoplex_*`) are a structural S/T/Y atlas without drug–kinase
affinities and are not used for masking. Python dependency added for Davis .xls parsing: `xlrd==2.0.2` (pip, plm_env).
