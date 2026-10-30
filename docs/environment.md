# Computing environment

Versions used to produce the results (October 2026).

| Item | Version |
|---|---|
| R | 4.5.2 (2025-10-31, ucrt) |
| OS | Windows 11 |
| Rtools | 4.5 (to build `extendedForest` / `gradientForest` from source) |
| arrow | 25.0.1 |
| dplyr | 1.2.0 |
| tidyr | 1.3.2 |
| terra | 1.9.50 |
| sf | 1.1.0 |
| data.table | 1.18.2.1 |
| vegan | 2.7.2 |
| permute | 0.9.10 |
| indicspecies | 1.8.0 |
| cluster | 2.1.8.1 |
| rpart | 4.1.24 |
| gdm | 1.6.0.7 |
| gradientForest | 0.1.37 (R-Forge, GitHub mirror `r-forge/gradientforest`) |
| extendedForest | 1.6.2 (R-Forge) |
| CopernicusMarine | 0.4.9 |
| jsonlite | 2.0.0 |
| base64enc | 0.1.6 |

`scripts/00_ambiente.R` records the versions of the data-handling packages, and `scripts/80_pacotes.R` installs and smoke-tests the analysis packages (`80b_gradientforest.R` builds `gradientForest` when R-Forge binaries are not available).

## Hardware

Intel Core i7-14700 (28 threads), 64 GB RAM. Heavy steps (PERMANOVA with permutations, equal-effort resampling, Gradient Forest, pairwise tests) run in parallel on the CPU with `parallel` PSOCK clusters of 7–14 workers. A full run of stages 1–3 takes several hours; the OBIS download and curation (phases 0–7) take about a day.

## Reproducibility notes

- Random steps use fixed seeds written in the scripts (for example 1001 for the draw of 10 events per unit).
- `permanova_rapida()` in `scripts/funcoes_permanova.R` reproduces the `vegan::adonis2` statistic for large *n* (about 65× faster); it was checked against `adonis2` and `betadisper`.
- Some scripts resume from intermediate files and skip work already done; delete the corresponding files in `intermediario/` to force a full rerun.
