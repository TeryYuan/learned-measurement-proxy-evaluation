# File-level license scope

Yuan Cao owns the original code and original contributions in this collection. Copyright (c) 2026 Yuan Cao. The author affiliation is Independent Researcher.

The exact machine-readable allocation is [FILE_LICENSES.json](FILE_LICENSES.json); the [reader index](source_index.tsv) also states each file's scope. An author's license applies only to the author's original contribution. Source labels, source measurements, covered adaptations and embedded third-party software retain their existing terms.

| Material | Author's license | Retained conditions |
|---|---|---|
| Original Python/R code, tests, tools, executable configuration and dependency specifications listed in the map | MIT | Installed dependencies and any upstream notices retain their own terms |
| Original documentation, metadata, clinical aggregate results and figure artwork | CC BY 4.0 | Excludes third-party quotations, legal notices and embedded font software |
| Original DILI results using FDA labels, original grouping/sampling records, and author-generated similarity/coverage summaries | CC BY 4.0 where listed | FDA and chemical-source provenance remains; no grant over source facts or omitted structures |
| Retained SI002/SI004 source labels or covered adaptations, and specifically listed BSEP/bridge summaries | CC BY-NC 4.0 for listed author contributions | Preserve the attachments' CC BY-NC 4.0 source conditions, attribution, DOI, license link and change notice |
| TG-GATEs measurements, transformed outcomes, covered derivatives, and listed ALT/semisynthetic result contributions | CC BY-SA 2.1 Japan for listed author contributions | Preserve applicable TG-GATEs terms and the unchanged Additional License |
| Existing third-party license notices and embedded DejaVu fonts | Their own terms | The author does not relicense them |

Using NC for selected BSEP summaries and compatible SA for selected ALT and semisynthetic summaries is an explicit licensing choice for the author's contributions. It is not a claim that every statistical result or model prediction automatically inherits a training dataset's license. The file map distinguishes actual source fields from author-generated metrics and metadata.

In particular, `all_routes_absolute_and_increment.csv` contains DILI S23 result summaries, and `model_information_blocks.csv` contains original method descriptions; neither is treated as a mixed TG measurement table. `NAM_context.json` and `NAM_indices.csv` are original sampling metadata. Seven chemical-coverage attachments contain saved similarities and coverage information without source structures or measured BSEP labels. Conversely, the scaffold-stratified input table retains BSEP labels and its source conditions.

The MIT text is in [LICENSE](LICENSE). The CC application notices and canonical legal-code links are in [LICENSES](LICENSES/README.md). [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) gives attribution and transformations. Different licenses in this collection do not grant users a choice to ignore the conditions attached to material they reuse.

Raw NHANES files and participant-level records, original third-party workbooks/CEL/GMT files, enriched chemical caches, fitted model objects, unpublished manuscript drafts and private author materials are outside this distribution. The complete historical liver training pipeline is outside the validated reproduction scope.

GitHub and Zenodo identifiers are recorded only after they exist. [RELEASE_STATUS.json](RELEASE_STATUS.json) gives the current publication state.
