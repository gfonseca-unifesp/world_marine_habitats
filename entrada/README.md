# Curated input tables (`entrada/`)

Small lookup tables that encode decisions approved by the PI at decision Gates (see [`../docs/decision_log.md`](../docs/decision_log.md)). The scripts read them from `WMH_ROOT/entrada/`; copy this folder there before running the pipeline. Values are in Portuguese because the scripts match them literally (see [`../docs/glossary.md`](../docs/glossary.md)).

| File | Produced / used by | Content |
|---|---|---|
| `regioes_nivel2.csv` | `34_janelas.R` (selection), Gate 1 item G4; read by `71_extrair_nivel2.R` and later level-2 scripts | The 8 level-2 windows (500 × 500 km): order, window id, name, centre, latitude zone (`faixa_lat`), reason for selection (`motivo`), bounding box (`lon_min`, `lat_min`, `lon_max`, `lat_max`, WGS84) |
| `taxon_compartimento_proposta.csv` | `72_proposta_compartimento.R` | Proposed compartment and ecological group per family, from WoRMS functional groups and FishBase, before review |
| `taxon_compartimento.csv` | approved version of the table above (stop G5, 2,312 families); read by `73_curadoria_n2.R` | Per family: broad group, taxonomy, number of records and species, WoRMS/FishBase evidence counts (`bentos`, `pelagico`, `misto`, `fb_*`), final `compartimento` (`bentos` benthic, `pelagico` pelagic, `ambiguo` decided per record, `fora` excluded), its source, `grupo_ecologico` (ecological group) and its source |
| `aparelho_conjunto.csv` | `72b_aparelho_conjunto.R` (Gate 3 item 4); read by `73_curadoria_n2.R` | Sampling gear inferred per OBIS dataset from its title or abstract: `dataset_id`, title, inferred compartment of the gear (`aparelho_conjunto`: `bentos`, `pelagico` or empty) and the reason (`fonte`, e.g. `tipo de conjunto sem aparelho` dataset type says nothing about gear, `sem termos` no gear terms found) |
