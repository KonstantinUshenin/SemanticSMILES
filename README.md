# SemanticSMILES — public metadata and code package

This repository describes the **SemanticSMILES** project: curated SMILES strings linked to chemistry dataset keys and semantic tags.

**This release ships metadata and build code only.** The molecule corpus (`SemanticSMILES.csv`, ~5M SMILES) is **not** redistributed here because many upstream sources lack a clear redistributable data licence. See [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md) for how to locate sources and rebuild locally.

## Contents

```
SemanticSMILES_public/
├── README.md              # this file
├── DATA_AVAILABILITY.md   # why the corpus is absent; how to rebuild
├── LICENSE                # MIT for code; notes on metadata / rebuilt data
├── metadata/              # dataset overview, references, data-source notes
└── code/                  # notebooks to build SemanticSMILES from local sources
    ├── paths.yaml
    ├── .env.example
    ├── README.md
    └── notebooks/         # steps 01–05 + helpers
```

### `metadata/`

| File | Contents |
|------|----------|
| `supplementary_dataset_overview.csv` | Dataset keys, sizes, years (716 entries) |
| `supplementary_dataset_references.csv` | Dataset names, source links, citations |
| `supplementary_to_the_data_sources_section.csv` | Prose descriptions of major chemical resources |

### `code/`

Jupyter notebooks that turn a local folder of source chemistry datasets into curated SMILES tables, tags, and summaries. Configure paths in `paths.yaml`, place each source under `workspace/<DATASET_KEY>/`, then run steps 01–05. Details: [`code/README.md`](code/README.md).

## Typical rebuild flow

1. Use `metadata/supplementary_dataset_references.csv` (and the overview CSV) to download sources you are allowed to use.
2. Configure and run [`code/`](code/) (see [`code/README.md`](code/README.md)).
3. Keep rebuilt molecule tables local; respect each upstream source’s redistribution terms.

## Licence

- **Code** (`code/`): MIT — see [LICENSE](LICENSE).
- **Metadata** in `metadata/`: attribute the SemanticSMILES project (CC BY 4.0 for curated project metadata).
- **Rebuilt molecule tables** inherit upstream source terms.
