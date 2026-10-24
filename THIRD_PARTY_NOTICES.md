# Third-party notices and attribution

The following notices cover material actually retained in this candidate. Original author contributions use the grants specified in [FILE_LICENSES.json](FILE_LICENSES.json) and [LICENSE_SCOPE.md](LICENSE_SCOPE.md). This repository is a collection of separately scoped materials; different licenses alone do not make the collection incompatible. No blanket license is applied to all files.

## FDA DILIrank

FDA Liver Toxicity Knowledge Base; DILIrank 1.0 (Chen et al., 2016) and DILIrank 2.0 (Olubamiwa et al., 2025). The historical workbook preserved on 2026-10-07 contains separate `version 1` and `version 2` sheets. Included N/C/bridge CSV labels and the N/bridge-DILI NPZ label arrays come from those records. Changes in the study: fixed binary label mapping, pre-existing identity qualification, cohort selection and saved score pairing. This release only selects columns; it does not change labels.

The [FDA website policy](https://www.fda.gov/about-fda/about-website/website-policies) permits reuse of FDA content unless otherwise noted and recommends source credit, date and a source link. The source is the [DILIrank dataset](https://www.fda.gov/science-research/liver-toxicity-knowledge-base-ltkb/drug-induced-liver-injury-rank-dilirank-20-dataset). The page asks downloaders to notify FDA; no notification was sent. The FDA policy is not applied to unrelated content hosted by PubChem.

## ACS supplementary datasets, distinct from the article

Marina Garcia de Lomana, Domenico Gadaleta, Marian Raschke, Robert Fricke and Floriane Montanari, “Predicting Liver-Related In Vitro Endpoints with Machine Learning to Support Early Detection of Drug-Induced Liver Injury,” supplementary dataset versions 1, posted 2025-03-11:

- [SI002, DOI 10.1021/acs.chemrestox.4c00453.s002](https://doi.org/10.1021/acs.chemrestox.4c00453.s002): historical auxiliary-model training source; original matrix and trained models are not bundled.
- [SI004, DOI 10.1021/acs.chemrestox.4c00453.s004](https://doi.org/10.1021/acs.chemrestox.4c00453.s004): experimental BSEP labels retained in `BSEP_predictions.csv`, `bridge_membership.csv`, `E_BSEP_bootstrap_indices.npz`, `E_no_anyA_bootstrap_indices.npz` and `bridge_indices_BSEP.npz` under `release_candidate/legacy_fixed_replay/data/`.

The official Figshare metadata identifies **CC BY-NC 4.0** for both attachments: [SI002 metadata](https://api.figshare.com/v2/articles/28571061), [SI004 metadata](https://api.figshare.com/v2/articles/28571067). The current API was read successfully on 2026-10-23. Preserve attribution, the [license link](https://creativecommons.org/licenses/by-nc/4.0/), indication of changes, and noncommercial conditions; no endorsement is implied. Included adaptations are selected original labels paired to the study's historical members, original sampling positions and author predictions. No label threshold or population was changed in this release.

The separate article's recorded CC BY-NC-ND 4.0 status is not substituted for the attachments' license. The article text and the full original attachments are not bundled. The listed noncommercial grants for author-generated BSEP/bridge summaries are an explicit author choice, not an automatic claim that every prediction trained using SI002 inherits NC conditions. Purely original DILI summaries and other listed material use their separately stated grants.

## Open TG-GATEs

Toxicogenomics Project and Toxicogenomics Informatics Project under CC Attribution-Share Alike 2.1 Japan.

The [Open TG-GATEs database](https://dbarchive.biosciencedbc.jp/en/open-tggates/desc.html), archived study inputs acquired 2026-10-08, supplies the rat ALT measurements underlying `NAM_predictions.csv` and the three NAM expected-summary CSVs. Historical changes: matched-condition response transformation and parent-level evaluation, paired with saved predictions; this release removes only the unused model-path column. Raw animal tables and CEL files are not included. Cite Igarashi et al. (2015), [doi:10.1093/nar/gku955](https://doi.org/10.1093/nar/gku955), and the database name and URL.

The applicable [Standard License is CC BY-SA 2.1 Japan](https://creativecommons.org/licenses/by-sa/2.1/jp/legalcode). The source [license page](https://dbarchive.biosciencedbc.jp/en/open-tggates/lic.html) requires its Additional License to accompany derivatives. That notice is reproduced unchanged in [sources/OPEN_TGGATES_ADDITIONAL_LICENSE.txt](sources/OPEN_TGGATES_ADDITIONAL_LICENSE.txt). TG-derived data remain under those conditions; they are not MIT. No assumption is made that removing compound names would extinguish source conditions.

## Fonts embedded in current figure PDFs

`figures/Figure6_NHANES_Q3.pdf` contains embedded Arial and Arial Bold subsets. The original font terms remain applicable; CC BY 4.0 covers the original figure artwork, with the separate icon conditions retained. See [figure_patch/rights/LICENSE_SCOPE.md](figure_patch/rights/LICENSE_SCOPE.md). No standalone font files are distributed.

`figures/FigureS12_CBC_distributions.pdf` retains three DejaVu Sans and DejaVu Sans Bold subset font resources. Bitstream Vera copyright/permission, public-domain DejaVu changes, and the bundled Arev notice are reproduced in [sources/LICENSE_DEJAVU.txt](sources/LICENSE_DEJAVU.txt). The plotting process subsets fonts; this project has not edited their glyph designs. See the [DejaVu project](https://dejavu-fonts.github.io/) for the upstream family. The PDF bytes and embedded font programs are unchanged in this preparation.

## Source references and install-only dependencies

Hallmark H gene sets **2025.1.Hs**, 50 sets, were used in historical expression processing; no GMT/member matrix or GSEA executable is bundled. Attribution: Molecular Signatures Database, Broad Institute, MIT and the Regents of the University of California; see [MSigDB terms](https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp) (CC BY 4.0 for the applicable releases, with collection-specific exceptions) and the retained [version record](sources/Hallmark_source_and_license.txt). Current MSigDB releases must not silently replace 2025.1.Hs. The KEGG MEDICUS exception concerns another collection, not these Hallmark H sets.

PubChem annotation/structure caches are not bundled. [NCBI policies](https://www.ncbi.nlm.nih.gov/home/about/policies/) distinguish NCBI government information from contributor-controlled content. The excluded chemical representations in the enriched historical tables have not been granted a blanket PubChem license. Stable study IDs and the original grouping equivalence relations remain in the minimal input collection.

NHANES original survey/mortality files and the official R reader are not bundled. Official provenance and access conditions remain in [NHANES source notes](NHANES_Q3_Project/sources/SOURCES.md); included aggregate references are generated study results.

Python, NumPy, pandas, SciPy, scikit-learn, Matplotlib, Pillow, joblib and other packages in the dependency lock are installed separately. Historical R code calls affy, frma, hgu133plus2cdf, hgu133plus2.db, AnnotationDbi and jsonlite; these packages are not vendored. No source copy from the cited Bayer/liver-mtnn repository, official NHANES R reader, or third-party model weights was identified in the included code. Importing a dependency does not replace its license or establish author ownership of the calling code. The author has now confirmed ownership of the original code. Existing upstream code and dependency rights remain separate from that confirmation.

## Expanded saved-result attachments

The same source attributions apply to retained source-derived fields in scientific_SI. The reader index and recovery record identify each table and its projection. BSEP and bridge labels retain SI004 conditions; original DILI labels retain FDA provenance. NAM and measurement-derived ALT results preserve the applicable TG-GATEs conditions and Additional License. S40/S41 saved memberships, similarities and score summaries are delivered without chemical strings. Neither column selection nor numerical aggregation is treated as a new blanket permission. The current author grants are specified at file/contribution scope in FILE_LICENSES.json and LICENSE_SCOPE.md.
