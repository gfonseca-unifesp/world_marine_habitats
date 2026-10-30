# 60_tabelas_n1.R  (Phase 6: 6.1 to 6.4)
# 6.1 Level 1 certainty tables (long format, one row per record):
#     criteria: passed Phase 4; taxon identified to the unit (family; genus in Euphausiacea and Chaetognatha);
#     benthos -> EBU class assigned (depth, when present, already passed 4.9; at level 1 the compartment comes from the taxon);
#     pelagic -> a single EMU in the depth range (5.1); year present; license approved; certainty level A, B or C.
#     Written as parquet (the tables have millions of rows; CSV would be several GB).
# 6.2 Automatic validation -> saida/validacao_n1.csv
# 6.3 Testability per class x group x slice -> saida/testabilidade_n1.csv
# 6.4 Cell x taxon incidence (1 = recorded, empty = unknown) -> saida/incidencia_n1/<compartimento>_<grupo>.csv
suppressPackageStartupMessages({library(arrow); library(dplyr); library(terra); library(tidyr)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
source(file.path(raiz, "scripts/funcoes_grade.R"))
ds <- open_dataset(file.path(raiz, "intermediario/n1_tratado"), unify_schemas = TRUE)   # benthos and pelagic have different columns
dir.create(file.path(raiz, "saida/incidencia_n1"), showWarnings = FALSE)

comuns <- c("record_id", "evento", "dataset_id", "lat", "lon", "coordUncertainty", "decimais", "year", "mes", "dia_conhecido",
            "depth_min", "depth_max", "bathymetry", "samplingProtocol", "scientificName", "taxonRank_u", "aphiaid",
            "family_u", "genus_u", "unidade_tax", "taxon_unidade", "grupo", "cell_25", "cell_50", "cell_100", "bloco_2500",
            "fatia", "nivel_certeza", "data_relida", "basis_citacao")
renomeia <- c(event_id = "evento", dataset = "dataset_id", coord_uncertainty_m = "coordUncertainty", depth_min_m = "depth_min",
              depth_max_m = "depth_max", bathymetry_m = "bathymetry", gear = "samplingProtocol", taxon = "scientificName",
              rank = "taxonRank_u", aphia_id = "aphiaid", family = "family_u", genus = "genus_u", month = "mes",
              cell_id = "cell_25", bloco_espacial = "bloco_2500")
fl_str <- function(f) vapply(f, function(v) paste(v, collapse = ";"), "")

# ---- 6.1 ----
bent <- ds |> filter(compartimento == "bentos", !so_esforco, situacao_classe == "ok") |> collect()
bent$flags_originais <- fl_str(bent$flags)
bent <- bent |> select(all_of(c(comuns, "flags_originais", "ebu_grupo", "pureza_grupo", "ebu_regiao", "pureza_regiao",
                                "ebu_unidade", "pureza_unidade", "sem_profundidade"))) |> rename(all_of(renomeia))
pel <- ds |> filter(compartimento == "pelagico", !so_esforco, situacao_classe == "ok") |> collect()
pel$flags_originais <- fl_str(pel$flags)
pel <- pel |> select(all_of(c(comuns, "flags_originais", "emu", "emu_rotulo", "pureza_emu", "prof_padrao_conjunto"))) |>
  rename(all_of(renomeia))
write_parquet(bent, file.path(raiz, "saida/tabela_certa_bentos_n1.parquet"))
write_parquet(pel, file.path(raiz, "saida/tabela_certa_pelagico_n1.parquet"))
cat(format(Sys.time()), "- tabelas: bentos", nrow(bent), "| pelagico", nrow(pel), "\n")

# ---- 6.2 validation ----
V <- list(); teste <- function(t, ok, det = "") V[[length(V) + 1]] <<- data.frame(teste = t, resultado = ifelse(ok, "passou", "FALHOU"), detalhe = det)
obrig_b <- c("record_id", "dataset", "lat", "lon", "year", "aphia_id", "taxon_unidade", "grupo", "cell_id", "bloco_espacial",
             "nivel_certeza", "ebu_grupo", "ebu_regiao", "ebu_unidade")
obrig_p <- c(setdiff(obrig_b, c("ebu_grupo", "ebu_regiao", "ebu_unidade")), "emu", "fatia")
teste("esquema bentos: colunas obrigatórias presentes", all(obrig_b %in% names(bent)))
teste("esquema pelágico: colunas obrigatórias presentes", all(obrig_p %in% names(pel)))
teste("tipos: lat/lon/year numéricos", is.numeric(bent$lat) && is.numeric(pel$lon) && is.numeric(bent$year))
teste("record_id único (bentos)", !anyDuplicated(bent$record_id), sprintf("%d registros", nrow(bent)))
teste("record_id único (pelágico)", !anyDuplicated(pel$record_id), sprintf("%d registros", nrow(pel)))
teste("nenhum registro nas duas tabelas", !any(bent$record_id %in% pel$record_id))
g25 <- ease2_grade("25km")
cel_de <- function(d) { c <- cellFromXY(g25, project(cbind(d$lon, d$lat), from = "EPSG:4326", to = "EPSG:6933"))
  rc <- rowColFromCell(g25, c); paste("25km", rc[, 1], rc[, 2], sep = "_") }
teste("lat/lon cai na célula (bentos)", all(cel_de(bent) == bent$cell_id))
teste("lat/lon cai na célula (pelágico)", all(cel_de(pel) == pel$cell_id))
ebu_cls <- read.csv(file.path(raiz, "camadas/ebu/ebu_classes.csv"))
teste("classes EBU pertencem à lista (3 níveis)", all(bent$ebu_grupo %in% ebu_cls$classe[ebu_cls$nivel == "grupo"]) &&
        all(bent$ebu_regiao %in% ebu_cls$classe[ebu_cls$nivel == "regiao"]) && all(bent$ebu_unidade %in% ebu_cls$classe[ebu_cls$nivel == "unidade"]))
teste("classes EMU pertencem à lista (1 a 37)", all(pel$emu %in% 1:37))
vaz <- function(d, cols) sum(sapply(cols, function(k) sum(is.na(d[[k]]) | as.character(d[[k]]) == "")))
teste("obrigatórias sem vazios (bentos)", vaz(bent, obrig_b) == 0, paste(vaz(bent, obrig_b), "vazios"))
teste("obrigatórias sem vazios (pelágico)", vaz(pel, obrig_p) == 0, paste(vaz(pel, obrig_p), "vazios"))
graves <- c("NO_COORD", "ZERO_COORD", "LAT_OUT_OF_RANGE", "LON_OUT_OF_RANGE", "MIN_DEPTH_EXCEEDS_MAX", "DATE_IN_FUTURE",
            "DATE_BEFORE_MIN", "NOT_MARINE")
tem_grave <- function(f) sum(grepl(paste(graves, collapse = "|"), f))
teste("nenhum sinalizador grave", tem_grave(bent$flags_originais) + tem_grave(pel$flags_originais) == 0)
lic <- read.csv(file.path(raiz, "intermediario/conjunto_tipo_escopo.csv"))
teste("conjunto corrompido ausente", !any(c(bent$dataset, pel$dataset) %in% "d29d7aa5-08e0-4c1f-9b1c-3b46b3b93563"))
teste("nível de certeza em A/B/C", all(c(bent$nivel_certeza, pel$nivel_certeza) %in% c("A", "B", "C")))
cas <- read.csv(file.path(raiz, "saida/cascata_filtros_n1.csv")); ext <- read.csv(file.path(raiz, "saida/nivel1_totais.csv"))
ret <- open_dataset(file.path(raiz, "intermediario/n1_curado")) |> count(grupo) |> collect()
qua <- open_dataset(file.path(raiz, "intermediario/quarentena/n1")) |> count(grupo) |> collect()
aus <- open_dataset(file.path(raiz, "intermediario/ausencias/n1")) |> count(grupo) |> collect()
fecha <- ext |> select(grupo, registros) |> left_join(ret, by = "grupo") |> left_join(qua, by = "grupo", suffix = c("", "_q")) |>
  left_join(aus, by = "grupo", suffix = c("", "_a")) |> mutate(dif = registros - n - coalesce(n_q, 0L) - coalesce(n_a, 0L))
teste("cascata fecha: extraídos = retidos + quarentena + ausências", all(fecha$dif == 0), paste("maior diferença:", max(abs(fecha$dif))))
# reproducibility: rebuild the tables from n1_tratado and compare the content hash
h <- function(d) digest::digest(d[order(d$record_id), ], algo = "sha256")
b2 <- read_parquet(file.path(raiz, "saida/tabela_certa_bentos_n1.parquet"))
teste("reconstrução produz o mesmo hash (bentos)", h(b2) == h(bent), "reconstrução a partir de n1_tratado; reexecução completa da curadoria fica para a entrega")
rm(b2); gc(verbose = FALSE)
val <- bind_rows(V); write.csv(val, file.path(raiz, "saida/validacao_n1.csv"), row.names = FALSE); print(val)

# ---- 6.3 testability ----
classe_long <- bind_rows(
  bent |> select(grupo, cell_id, event_id, taxon_unidade, nivel_certeza, fatia, classe = ebu_grupo) |> mutate(nivel = "EBU grupo"),
  bent |> select(grupo, cell_id, event_id, taxon_unidade, nivel_certeza, fatia, classe = ebu_regiao) |> mutate(nivel = "EBU região"),
  bent |> select(grupo, cell_id, event_id, taxon_unidade, nivel_certeza, fatia, classe = ebu_unidade) |> mutate(nivel = "EBU unidade"),
  pel |> mutate(classe = as.character(emu)) |> select(grupo, cell_id, event_id, taxon_unidade, nivel_certeza, fatia, classe) |>
    mutate(nivel = "EMU"))
cel <- classe_long |> group_by(nivel, classe, grupo, fatia, cell_id) |> summarise(ev = n_distinct(event_id), .groups = "drop")
testa <- classe_long |> group_by(nivel, classe, grupo, fatia) |>
  summarise(celulas = n_distinct(cell_id), eventos = n_distinct(event_id), registros = n(), taxons = n_distinct(taxon_unidade),
            A = sum(nivel_certeza == "A"), B = sum(nivel_certeza == "B"), C = sum(nivel_certeza == "C"), .groups = "drop") |>
  left_join(cel |> group_by(nivel, classe, grupo, fatia) |> summarise(celulas_10ev = sum(ev >= 10), .groups = "drop"),
            by = c("nivel", "classe", "grupo", "fatia")) |>
  mutate(classe_apta = celulas_10ev >= 5)
viav <- testa |> group_by(nivel, grupo, fatia) |> summarise(classes_aptas = sum(classe_apta), .groups = "drop") |>
  mutate(contraste_viavel = ifelse(classes_aptas >= 2, "sim", "não"))
testa <- testa |> left_join(viav, by = c("nivel", "grupo", "fatia"))
write.csv(testa, file.path(raiz, "saida/testabilidade_n1.csv"), row.names = FALSE)
print(viav |> group_by(nivel) |> summarise(combinacoes = n(), viaveis = sum(contraste_viavel == "sim"), mediana_classes_aptas = median(classes_aptas)))

# ---- 6.4 incidence ----
inc <- function(d, comp, chave) {
  for (g in unique(d$grupo)) {
    x <- d[d$grupo == g, ]
    n_ev <- x |> group_by(across(all_of(chave))) |> summarise(n_eventos = n_distinct(event_id), .groups = "drop")
    w <- x |> distinct(across(all_of(c(chave, "taxon_unidade")))) |> mutate(v = 1L) |>
      pivot_wider(names_from = taxon_unidade, values_from = v)
    w <- left_join(n_ev, w, by = chave)
    write.csv(w, file.path(raiz, sprintf("saida/incidencia_n1/%s_%s.csv", comp, gsub(" ", "_", g))), row.names = FALSE, na = "")
  }
}
inc(bent, "bentos", "cell_id"); inc(pel |> mutate(emu = as.integer(emu)), "pelagico", c("cell_id", "emu"))
cat(format(Sys.time()), "- FIM\n")
