# 60b_tabela_n1b.R  (option B, Gate 11: certainty table of the new pelagic groups, same columns as 60_tabelas_n1.R)
# Output: saida/tabela_certa_pelagico_n1b.parquet
suppressPackageStartupMessages({library(arrow); library(dplyr)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
ds <- open_dataset(file.path(raiz, "intermediario/n1b_tratado"), unify_schemas = TRUE)
comuns <- c("record_id", "evento", "dataset_id", "lat", "lon", "coordUncertainty", "decimais", "year", "mes", "dia_conhecido",
            "depth_min", "depth_max", "bathymetry", "samplingProtocol", "scientificName", "taxonRank_u", "aphiaid",
            "family_u", "genus_u", "unidade_tax", "taxon_unidade", "grupo", "cell_25", "cell_50", "cell_100", "bloco_2500",
            "fatia", "nivel_certeza", "data_relida", "basis_citacao")
renomeia <- c(event_id = "evento", dataset = "dataset_id", coord_uncertainty_m = "coordUncertainty", depth_min_m = "depth_min",
              depth_max_m = "depth_max", bathymetry_m = "bathymetry", gear = "samplingProtocol", taxon = "scientificName",
              rank = "taxonRank_u", aphia_id = "aphiaid", family = "family_u", genus = "genus_u", month = "mes",
              cell_id = "cell_25", bloco_espacial = "bloco_2500")
pel <- ds |> filter(compartimento == "pelagico", !so_esforco, situacao_classe == "ok") |> collect()
pel <- pel |> select(all_of(c(comuns, "emu", "emu_rotulo", "pureza_emu", "prof_padrao_conjunto"))) |> rename(all_of(renomeia))
write_parquet(pel, file.path(raiz, "saida/tabela_certa_pelagico_n1b.parquet"))
print(pel |> group_by(grupo) |> summarise(registros = n(), desde_1993 = sum(year >= 1993, na.rm = TRUE), eventos = n_distinct(event_id), conjuntos = n_distinct(dataset), generos = n_distinct(taxon_unidade)) |> as.data.frame())
cat(format(Sys.time()), "FIM\n")
