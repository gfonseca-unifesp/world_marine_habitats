# 60c_tabela_n1c.R  (Gate 12, 13/11/2026: certainty tables of the test-fauna groups of 21c, same columns as 60_tabelas_n1.R)
# Outputs: saida/tabela_certa_bentos_n1c.parquet, saida/tabela_certa_pelagico_n1c.parquet
suppressPackageStartupMessages({library(arrow); library(dplyr)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
ds <- open_dataset(file.path(raiz, "intermediario/n1c_tratado"), unify_schemas = TRUE)
comuns <- c("record_id", "evento", "dataset_id", "lat", "lon", "coordUncertainty", "decimais", "year", "mes", "dia_conhecido",
            "depth_min", "depth_max", "bathymetry", "samplingProtocol", "scientificName", "taxonRank_u", "aphiaid",
            "family_u", "genus_u", "unidade_tax", "taxon_unidade", "grupo", "cell_25", "cell_50", "cell_100", "bloco_2500",
            "fatia", "nivel_certeza", "data_relida", "basis_citacao")
renomeia <- c(event_id = "evento", dataset = "dataset_id", coord_uncertainty_m = "coordUncertainty", depth_min_m = "depth_min",
              depth_max_m = "depth_max", bathymetry_m = "bathymetry", gear = "samplingProtocol", taxon = "scientificName",
              rank = "taxonRank_u", aphia_id = "aphiaid", family = "family_u", genus = "genus_u", month = "mes",
              cell_id = "cell_25", bloco_espacial = "bloco_2500")
bent <- ds |> filter(compartimento == "bentos", !so_esforco, situacao_classe == "ok") |> collect() |>
  select(all_of(c(comuns, "ebu_grupo", "pureza_grupo", "ebu_regiao", "pureza_regiao", "ebu_unidade", "pureza_unidade"))) |> rename(all_of(renomeia))
pel <- ds |> filter(compartimento == "pelagico", !so_esforco, situacao_classe == "ok") |> collect() |>
  select(all_of(c(comuns, "emu", "emu_rotulo", "pureza_emu", "prof_padrao_conjunto"))) |> rename(all_of(renomeia))
write_parquet(bent, file.path(raiz, "saida/tabela_certa_bentos_n1c.parquet")); write_parquet(pel, file.path(raiz, "saida/tabela_certa_pelagico_n1c.parquet"))
print(bind_rows(bent, pel) |> group_by(grupo) |> summarise(registros = n(), desde_1993 = sum(year >= 1993, na.rm = TRUE), eventos = n_distinct(event_id),
  conjuntos = n_distinct(dataset), generos = n_distinct(taxon_unidade)) |> as.data.frame())
cat(format(Sys.time()), "FIM\n")
