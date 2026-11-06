# Workflow

Scripts are numbered in run order. Every script reads and writes inside the working folder `WMH_ROOT` (see `.Renviron.example`); folder names are explained in [`glossary.md`](glossary.md). Each script starts with a header that states its purpose, method, inputs and outputs. Decisions taken at each Gate are in [`decision_log.md`](decision_log.md).

## Phase 0 — environment and raw data

| Script | Purpose | Main outputs |
|---|---|---|
| `00_ambiente.R` | Record R and package versions | `saida/ambiente.csv` |
| `01_baixar_obis.R` | Local copy of the full OBIS occurrence snapshot (AWS Open Data) | `bruto/occurrence/*.parquet` |
| `02_verificar_obis.R` | Verify the copy against the manifest (size; MD5 = ETag for single-part files) | `bruto/verificacao.csv` |

## Phase 1 — grids and habitat layers

| Script | Purpose | Main outputs |
|---|---|---|
| `funcoes_grade.R` | Nested equal-area EASE-Grid 2.0 (EPSG:6933), 25 km family | (functions) |
| `10_grades.R` | Document and test the grids | `camadas/grades/` |
| `11_ebu_baixar.R` | Download the EBU polygons (Esri feature service) | `camadas/ebu/` |
| `12_ebu_celulas.R` | EBU group, region and unit per grid cell, with purity | `camadas/ebu/ebu_celula_25km.parquet` |
| `13_emu_celulas.R` | EMU class per cell in each of the 102 depth levels (3-D) | `camadas/emu/emu_25km.tif` |
| `14_glorys_amplitude.R` | Seasonal amplitude of GLORYS12 fields, used to score level-2 windows | `camadas/glorys/` |
| `funcoes_copernicus.R` | Load Copernicus Marine credentials from a local file into R options (never printed) | (function) |

## Phases 2–3 — extraction, profiling and scoping

| Script | Purpose |
|---|---|
| `21_extrair_nivel1.R` | Extract all records of the level-1 focal taxa from the OBIS snapshot (no filtering yet) |
| `22_perfil_nivel1.R` | Profile of the extraction: field completeness, dates, coordinate precision, flags, basisOfRecord, licences, dataset dominance |
| `21b_extrair_nivel1_b.R` | Option B: extraction of the small pelagic fish and micronekton families, apart from the approved level 1 |
| `23_taxa.R` | Taxon table from the WoRMS classification shipped with OBIS |
| `24_conjuntos.R` | Dataset metadata, type, scope and standardised licence (OBIS API) |
| `31_sondagem.R` | Scoping: purity, effort, testable classes and boundary pairs at 25, 50 and 100 km |
| `34_janelas.R` | Candidate level-2 windows (500 × 500 km), scored by effort, contrasts and dynamics |

## Phases 4–6 — curation and certainty tables (level 1)

| Script | Purpose |
|---|---|
| `40_curadoria_n1.R` | Cascade of quality filters 4.1–4.10 (one quarantine reason per record) |
| `50_tratamento_n1.R` | Habitat class of every retained record: EBU (benthic) and 3-D EMU (pelagic), with purity |
| `60_tabelas_n1.R` | Certainty tables (levels A, B, …) per record and validation summaries |

## Phase 7 — level 2 (regional windows)

| Script | Purpose |
|---|---|
| `71_extrair_nivel2.R` | Extract focal and extended taxa in the 8 windows |
| `71b_aparelho_emof.R` | Sampling gear from the eMoF extension, for records without `samplingProtocol` |
| `72_proposta_compartimento.R` | Proposed compartment and ecological group per family (WoRMS functional group, FishBase), approved at a Gate |
| `72b_aparelho_conjunto.R` | Sampling gear inferred from dataset title/abstract |
| `73_curadoria_n2.R` | Level-2 curation, position precision and compartment per record |
| `74_tabelas_n2.R` | Classes at 25, 12.5 and 6.25 km, resolution choice, validation and testability per window |

## Phase 8 — analyses

### Stage 1: do EBU and EMU separate the fauna?

| Script | Purpose |
|---|---|
| `80_pacotes.R`, `80b_gradientforest.R` | Install and smoke-test vegan, indicspecies, gdm, gradientForest |
| `funcoes_permanova.R` | Fast PERMANOVA (same statistic as `vegan::adonis2`, class after block, permutations within blocks), partial ω², equal-effort resampling with power, pairwise tests |
| `81_etapa1_n1.R` | Level 1: PERMANOVA on √β<sub>sim</sub> (10 events drawn per unit, 5 draws), equal-effort ω² and power, IndVal.g; argument `unidade=genero` for the genus sensitivity |
| `82_etapa1_n2.R` | Level 2 within the 8 windows (blocks of 100 km; `bloco=50` sensitivity) |

### Stage 2: diagnosis, equivalences and bottom-up synthesis

| Script | Purpose |
|---|---|
| `funcoes_etapa2.R` | Unit preparation, leave-one-block-out classification by fauna centroid, discontinuity search (PAM, silhouette) |
| `83_etapa2.R` | Where classification fails (`nivel=1` or `nivel=2`) |
| `84_equivalencias.R` | Distinct / equivalent pairs per taxonomic group (family or genus), consensus across groups, conservative merges |
| `85_sintese.R` | Final bottom-up list (distinct on its own / partly / not distinct / stand-by) and global maps |
| `funcoes_mapa.R` | Class maps with stand-by classes in their own colour and hatched |

### Stage 3: gradients, thresholds and new habitats

| Script | Purpose |
|---|---|
| `86_glorys_campos.R` | Monthly GLORYS12 fields on the 25 km grid (T and S at 14 depths, bottom T/S/current, surface current, ice, MLD, depth) |
| `86b_bacias.R` | Ocean basins (rule-based, ocean-mask connectivity, Panama barrier) used as a biogeographic covariate |
| `87_preditores.R` | Sampling units (≥ 1993, with month) and predictors; response = frequency of occurrence |
| `funcoes_gf.R` | Thresholds from Gradient Forest cumulative-importance curves |
| `88_gf_gdm.R` | Gradient Forest and GDM per taxonomic group (benthic without bottom salinity since Gate 8) |
| `88c_transferencia.R` | Do thresholds transfer across basins? (local vs leave-basin-out GF) |
| `89_ligacao_habitats.R` | Environmental profile of habitats, GF distance between distinct and equivalent pairs (AUC), robust thresholds |
| `89b_habitats_novos.R` | New habitats independent of EBU/EMU: GF-transformed environment, k-means by depth band (benthic) or month × depth (pelagic), merges guided by the fauna, threshold rules (rpart), validation (ω²) and agreement (ARI) |
| `89c_bentos_sazonal.R` | Does the benthos need dynamic habitats? (intra-annual variation, habitat stability, fauna between quarters) |

### Gate 9: united map

| Script | Purpose |
|---|---|
| `91_mapa_unido_bentos.R` | Final benthic map: depth band × new habitat × EBU region candidates on one tree whose dissimilarity is a false-merge rate (fauna-equivalent = 0; distinct, other band or other protected feature = never; untested = share of fauna-distinct pairs at a smaller GF distance); base cut at the 25th percentile of equivalent pairs and nested tolerance levels; protected features (seamounts, ridges); minimum area; evidence attributes; acceptance test; criterion 9 applies the literature splits of `92c` |
| `91b_mapa_unido_pelagico.R` | Final pelagic map with the same rules on the dynamic grid (cell × month × depth): vertical zone × pelagic habitat (HP, updated with level 2) × EMU class; fauna of levels 1 and 2; persistence (resident / seasonal / occasional) |

### Pelagic at level 2

| Script | Purpose |
|---|---|
| `93_pelagico_n2.R` | Level-2 pelagic units (cell × quarter × depth layer), Gradient Forest per ecological group, thresholds vs level 1, level-1 pelagic clusters tested with the level-2 fauna |

### Gate 10: published thresholds

| Script | Purpose |
|---|---|
| `92_limiares_publicados.R` | L1 standardise the published-threshold base (group, predictor, study box, basin, independence heuristic, sampling-stratum flag); L2 compare with the model robust and basin thresholds |
| `92b_enriquecimento.R` | L2b tolerance-free test: does the published interval concentrate more faunal turnover on the GF curve than random intervals of the same width? Gate 10 use and weight; stratum diagnostic of the pelagic z |
| `92c_redefinicao.R` | L3 undetected published depth / bottom-temperature thresholds tested as new class boundaries with the fauna inside the study box |
| `92d_afericao.R` | L4 in-situ vs GLORYS depth in the study boxes, published values moved to the GLORYS scale, L2 and L2b recomputed; variables in the literature missing from the model |
| `92f_tabela_validacao.R` | Worksheet for checking each published threshold against the full text (DOIs in `entrada/limiares_referencias.csv`) |
| `92e_figura_limiares.R` | L5 model vs published thresholds per variable and study, with threshold ids (key table in `saida/figura_limiares_tabela.csv`) |

`94_comparacao_mapas.R` draws, in one common style, the original classification, the classification tested by the fauna, the new habitats and the united map, for the benthos and for the pelagic.

The internal HTML validation report (`90_relatorio.R`) is not part of this repository.

## Stage 4 (planned)

Monthly dynamic pelagic maps, validation with blocks and years left out, area of applicability, and a wider level-1 pelagic set (small pelagic fish, micronekton, gelatinous plankton).
