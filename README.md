# World Marine Habitats

Code and decision record for testing global marine habitat classifications with species occurrence data, and for deriving new, fauna-validated habitat maps from ocean environmental gradients.

The project asks two questions:

1. **Do existing top-down classifications separate the fauna?** It tests the Ecological Benthic Units (EBU; Harris et al. 2026) and the Ecological Marine Units (EMU; Sayre et al. 2017) against OBIS occurrence records, globally (level 1, 25 km) and in eight regional windows (level 2, 6.25–25 km).
2. **Which environmental gradients and thresholds structure the fauna, and what habitat map follows from them?** It uses Gradient Forest and GDM on GLORYS12 ocean reanalysis fields, builds new habitats that are independent of EBU and EMU, and finally unites the top-down and bottom-up views in one map.

The repository accompanies a scientific article (in preparation). The article text, figures and the internal validation report are not part of this repository.

> **Status:** work in progress. Code is published now; output data (habitat maps, thresholds, tables) will be released with a tagged version once the results are final (see [`data/README.md`](data/README.md)).

## Repository layout

```
world_marine_habitats/
├── README.md               this file
├── LICENSE                 MIT licence (code)
├── LICENSE-DATA.md         CC BY 4.0 (output data)
├── CITATION.cff            how to cite
├── .Renviron.example       environment variables used by the scripts
├── scripts/                complete R pipeline, numbered in run order
├── docs/
│   ├── workflow.md         what each script does, inputs and outputs, by phase
│   ├── decision_log.md     every methodological decision (Gates 1–9), rationale and key results
│   ├── glossary.md         Portuguese identifiers used in code and outputs, with English meaning
│   └── environment.md      software versions and hardware
└── data/                   output data (released with the final version)
```

## Approach in brief

| Stage | Question | Main methods | Scripts |
|---|---|---|---|
| Phases 0–7 | Data acquisition, curation and certainty tables | Local OBIS snapshot, cascade of quality filters, compartment (benthic/pelagic) assignment, EBU/EMU class per record | `00`–`74` |
| Stage 1 | Do EBU/EMU classes separate the fauna? | PERMANOVA on √β<sub>sim</sub> within spatial blocks, ω² (effect size) with equal-effort resampling and power; IndVal | `81`, `82` |
| Stage 2 | Where does the classification fail, which classes are equivalent? | Leave-one-block-out classification, pairwise equivalence with consensus across taxonomic groups, bottom-up synthesis map | `83`–`85` |
| Stage 3 | Which gradients and thresholds matter? | Gradient Forest and GDM on GLORYS12 predictors, robust thresholds, basin covariate, transferability across basins | `86`–`89`, `88c` |
| Stage 3.5 | New habitats from gradients | GF-transformed environment, k-means by depth band, merges guided by the fauna, threshold rules (classification trees) | `89b`, `89c` |
| Gate 9 | Final map | Top-down × bottom-up candidates, merges only on fauna equivalence, protected features, minimum area | `91` |

Principles that run through the analysis:

- **Effect size over p-values.** With tens of thousands of sampling units almost everything is "significant"; decisions use partial ω² (and R²) with equal-effort resampling and power, not p alone.
- **Precaution.** Classes are merged only when the fauna shows equivalence with adequate power; untested classes stay separate and are shown hatched ("stand-by") as potential habitats.
- **Benthos static, pelagic dynamic.** Benthic habitats use annual summaries (the fauna showed negligible seasonal change); pelagic habitats are defined by month and depth.

## Data sources

| Source | Use | Licence |
|---|---|---|
| [OBIS](https://obis.org) full occurrence snapshot (AWS Open Data, `s3://obis-open-data/occurrence/`) | Fauna records | Per dataset (CC0, CC BY, CC BY-NC); dataset licences are recorded by `24_conjuntos.R` |
| EBU — Ecological Benthic Units (Harris et al. 2026; Esri feature service `EBU_Publication_WFL1`) | Top-down benthic classes | CC BY 4.0 |
| EMU — Ecological Marine Units v1.0 (Sayre et al. 2017; Esri netCDF `emu_v5_20200914.nc`) | Top-down pelagic classes (3-D) | CC BY 4.0 |
| [Copernicus Marine](https://marine.copernicus.eu) GLORYS12V1 monthly climatology 1993–2016 (`GLOBAL_MULTIYEAR_PHY_001_030`) | Environmental predictors | Copernicus Marine licence (free, attribution) |
| WoRMS classification (as distributed with OBIS) | Taxonomy | CC BY 4.0 |

Raw and curated input data are **not** redistributed here; the scripts download and rebuild them.

## Reproducing the analysis

Requirements: R ≥ 4.5 with the packages listed in [`docs/environment.md`](docs/environment.md); Rtools (Windows) to build `gradientForest` from source (`scripts/80b_gradientforest.R`); about 150 GB of free disk (the OBIS snapshot alone is about 90 GB); 64 GB RAM recommended. Heavy steps run in parallel on the CPU (`parallel`, PSOCK clusters).

1. Copy `.Renviron.example` to `.Renviron` (or set the variables in your session) and point `WMH_ROOT` to a working folder. The scripts create `bruto/` (raw), `camadas/` (layers), `intermediario/` (intermediate), `entrada/` (curated lookup tables) and `saida/` (outputs) inside it.
2. Copernicus Marine downloads need your own account; credentials are read from a local file outside the repository (see `scripts/funcoes_copernicus.R`). Never commit credentials.
3. Run the scripts in numeric order (`Rscript scripts/00_ambiente.R`, …). Some scripts take arguments, for example `83_etapa2.R nivel=1` or `82_etapa1_n2.R bloco=100 unidade=genero`; see the header of each script and [`docs/workflow.md`](docs/workflow.md).
4. Steps that stopped at a decision Gate were run with the options recorded in [`docs/decision_log.md`](docs/decision_log.md).

Code identifiers, file names and output column names are in Portuguese (the working language of the project); comments and documentation are in English. [`docs/glossary.md`](docs/glossary.md) maps the identifiers to English.

## How to cite

See [`CITATION.cff`](CITATION.cff). Please also cite OBIS, the EBU and EMU publications and Copernicus Marine when using the outputs.

## Licence

- Code: [MIT](LICENSE).
- Output data in `data/`: [CC BY 4.0](LICENSE-DATA.md).

Input data keep their own licences (see the table above).
