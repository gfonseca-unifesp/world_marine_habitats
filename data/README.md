# Output data

This folder will hold the output data of the project, released under [CC BY 4.0](../LICENSE-DATA.md) with a tagged version of the repository once the results are final (after the Gate 9 united map and stage 4, the dynamic monthly maps).

Planned content (file and column names will be in English, with a data dictionary):

| Folder | Content |
|---|---|
| `habitats/` | Habitat rasters on the EASE-Grid 2.0 25 km grid (EPSG:6933): bottom-up synthesis of EBU and EMU, new benthic and pelagic habitats (pelagic by month × depth), the united benthic map; class legends and attributes |
| `thresholds/` | Gradient Forest importance, cumulative-importance curves (resampled), thresholds per taxonomic group, robust thresholds and their transferability across basins |
| `tests/` | Stage 1–3 test results: PERMANOVA effect sizes, pairwise distinctness and consensus, classification diagnostics, new-habitat validation, seasonality test |
| `units/` | Sampling units per taxonomic group (cell, compartment, block, basin, environmental predictors) without the underlying occurrence records |

Raw OBIS records, curated record tables and third-party layers (EBU, EMU, GLORYS12) are not redistributed; the scripts rebuild them from the original sources.
