# Data availability

## Why the molecule corpus is not in this package

The full SemanticSMILES table (`SemanticSMILES.csv`, on the order of ~5M SMILES) aggregates many third-party chemistry datasets. Some upstream sources allow redistribution; others are freely downloadable for research but do **not** state a clear permissive data licence for redistribution.

This public package therefore distributes **metainformation and build code**, not the aggregated SMILES corpus or per-source molecule dumps. Researchers can obtain sources from the listed URLs and rebuild locally where applicable law and upstream terms allow.

## Where to find source locations

| Resource | Role |
|----------|------|
| [`metadata/supplementary_dataset_references.csv`](metadata/supplementary_dataset_references.csv) | Dataset names, source links, citations |
| [`metadata/supplementary_dataset_overview.csv`](metadata/supplementary_dataset_overview.csv) | Dataset keys, approximate sizes, years |
| [`metadata/supplementary_to_the_data_sources_section.csv`](metadata/supplementary_to_the_data_sources_section.csv) | Descriptions of major curated chemical resources (PubChem, ChEMBL, etc.) |

## How to rebuild

1. **Download sources** using the reference and overview CSVs. Respect each source’s terms.
2. **Lay out sources** under a workspace directory as `workspace/<DATASET_KEY>/` (see [`code/README.md`](code/README.md)).
3. **Configure** [`code/paths.yaml`](code/paths.yaml) (`workspace`, `workspace_data`, `workspace_summary`). Copy [`code/.env.example`](code/.env.example) → `code/.env` if using LLM steps.
4. **Run notebooks** `step01` … `step05` in [`code/notebooks/`](code/notebooks/).

## Licence stance for rebuilt data

- Any locally rebuilt molecule tables remain subject to upstream source licence terms.
- This package does not redistribute the aggregated corpus or per-source dumps.
- Code remains under MIT ([LICENSE](LICENSE)).
