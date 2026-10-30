# Glossary of Portuguese identifiers

Script names, function and variable names, folder names, output file names and output column names were kept in Portuguese, the working language of the project, so that the code that produced the results is exactly the code published here. Comments and documentation are in English. This table maps the recurring terms.

## Folders (inside `WMH_ROOT`)

| Identifier | Meaning |
|---|---|
| `bruto/` | raw data (OBIS snapshot, API caches) |
| `camadas/` | spatial layers (grids, EBU, EMU, GLORYS12 fields, basins) |
| `entrada/` | curated lookup tables approved at decision Gates (taxon → compartment, regional windows, sampling gear by dataset) |
| `intermediario/` | intermediate files (extractions, model objects) |
| `saida/` | outputs (tables, rasters, maps) |

## Script and output name parts

| Term | Meaning |
|---|---|
| `baixar` | download |
| `verificar` | verify |
| `grades` / `grade` | grids / grid |
| `celulas` | cells |
| `extrair` | extract |
| `perfil` | profile (summary) |
| `conjuntos` | datasets |
| `sondagem` | scoping survey (preliminary analysis of testability) |
| `janelas` | windows (level-2 regional windows, 500 × 500 km) |
| `curadoria` | curation (quality filters) |
| `tratamento` | processing (class assignment) |
| `tabelas` | tables |
| `aparelho` | sampling gear |
| `proposta_compartimento` | proposed compartment (benthic/pelagic) per taxon |
| `pacotes` | packages |
| `etapa1`, `etapa2`, `etapa3` | stage 1, 2, 3 |
| `n1`, `n2`, `nivel1`, `nivel2` | level 1 (global, 25 km), level 2 (regional windows) |
| `equivalencias` | equivalences |
| `sintese` | synthesis (bottom-up final list and map) |
| `campos` | (environmental) fields |
| `bacias` | ocean basins |
| `preditores` | predictors |
| `ligacao_habitats` | link between environment and habitats |
| `habitats_novos` | new habitats (derived from gradients, independent of EBU/EMU) |
| `bentos_sazonal` | benthic seasonality test |
| `transferencia` | transferability (across basins) |
| `mapa_unido` | united map (top-down × bottom-up) |
| `relatorio` | report |
| `funcoes_*` | shared function files |
| `legenda` | legend |
| `limiares` (`limiares_robustos`) | thresholds (robust thresholds) |
| `importancia` | importance |
| `curvas` | curves |
| `pares` | pairs |
| `fusoes` | merges |
| `regras` | rules |
| `fidelidade` | fidelity (accuracy of the rules on held-out cells) |
| `validacao` | validation |
| `concordancia` | agreement (adjusted Rand index) |
| `avaliacao` | evaluation |
| `representacao` | representation |
| `comsfundo` / `semsfundo` | with / without bottom salinity |

## Frequent variables and columns

| Identifier | Meaning |
|---|---|
| `raiz` | root folder (`WMH_ROOT`) |
| `bentos`, `pelagico` | benthos (benthic), pelagic |
| `comp` | compartment (benthic or pelagic) |
| `grupo` | taxonomic group (stage 1–3) or EBU geomorphic group (`ebu_grupo`) |
| `regiao` (`ebu_regiao`) | EBU benthic region (57) |
| `unidade` | sampling unit (cell, or cell × EMU [× quarter] for the pelagic) |
| `evento` | sampling event |
| `bloco` | spatial block (2,500 km at level 1; 100 km at level 2) |
| `classe` | habitat class |
| `pureza` | purity (share of the cell occupied by its dominant class) |
| `nivel_certeza` | certainty level (A, B, …) |
| `faixa` | depth band (0–200, 200–2,000, > 2,000 m) |
| `bacia`, `bac_1` … `bac_8` | ocean basin and its 0/1 indicator variables |
| `limiar` | threshold |
| `separa` | separates (pair is distinct: ω² ≥ 0.01 and power ≥ 0.8) |
| `distinto`, `equivalente`, `depende do grupo` | distinct, equivalent, depends on the taxonomic group (pair consensus) |
| `poder_eq` | power under equal effort |
| `omega2`, `omega2_eq` | partial ω² (all data / equal effort) |
| `distingue-se por si só`, `distingue-se em parte`, `não se distingue`, `inconclusivo` | distinct on its own, partly distinct, not distinct, inconclusive |
| `stand-by: sem par` / `sem dados` | stand-by: no comparable pair / no data |
| `standby` | stand-by (untested; drawn hatched on maps) |
| `celulas` | number of cells |
| `prof` | depth (GLORYS12 bottom depth) |
| `tfundo_*`, `sfundo_*`, `vfundo_*` | bottom temperature, bottom salinity, bottom current speed (`_med` mean, `_min`, `_max`) |
| `sst_med`, `sst_amp` | mean and annual range of sea-surface temperature |
| `mld`, `mld_max` | mixed-layer depth (monthly / annual maximum) |
| `gelo`, `gelo_meses` | sea-ice concentration / months with ice ≥ 0.15 |
| `z`, `t_z`, `s_z`, `z_mld` | depth of the pelagic unit, temperature and salinity at z, log2 ratio of z to the mixed layer |
| `vsup` | surface current speed |
| `tri` (`T1`…`T4`) | quarter (Jan–Mar … Oct–Dec) |
| `HB`*n*, `HP`*n*, `HF`*n*, `U`*n* | new benthic habitat, new pelagic habitat, final united benthic class, union of merged classes (bottom-up synthesis) |
| `protegida` | protected feature (EBU groups C seamounts and D ridges/rifts) |
| `evidencia` | evidence type of a final class (`distinta pela fauna` distinct by fauna, `fundida pela fauna` merged by fauna, `testada, sem par conclusivo` tested without a conclusive pair, `não testada` untested) |
| `sazonalidade` (`forte`, `moderada`, `fraca`) | seasonality of shelf bottom temperature (strong, moderate, weak) |
