# 50_tratamento_n1.R  (Phase 5: 5.1 to 5.4)
# 5.1 Class of each retained record, at the level 1 resolution (25 km, Gate 1):
#     benthos  -> EBU of the cell at the three levels (group, region, unit) with the purities;
#     pelagic -> EMU from the native 3D (D4.1): the depth range of the record (min to max) is crossed with the EMU
#                 levels in the cell; a single EMU -> assigned; two or more -> pelagico_multi_emu; no depth ->
#                 documented default depth of the dataset (D4.3) or outside the EMU test (counts toward effort).
# 5.2 Level 1 time slice: pelagic in quarters (from the month); benthos static. Certainty level already comes from Phase 4.
# 5.3 Effort per cell x group x slice (x EMU in pelagic). Species lists (type d) do not count toward effort (G9).
# 5.4 Scope per dataset: focal groups that each dataset covers.
# Outputs: intermediario/n1_tratado/ (by group), saida/pelagico_multi_emu.parquet, saida/esforco_nivel1.parquet,
#         saida/conjunto_grupo.csv, saida/tratamento_n1_resumo.csv
suppressPackageStartupMessages({library(arrow); library(dplyr); library(terra)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
source(file.path(raiz, "scripts/funcoes_grade.R"))
unlink(file.path(raiz, "intermediario/n1_tratado"), recursive = TRUE)
ds  <- open_dataset(file.path(raiz, "intermediario/n1_curado"))
cj  <- read.csv(file.path(raiz, "intermediario/conjunto_tipo_escopo.csv"))
lista <- cj$dataset_id[cj$tipo == "d_lista_especies"]
prof_padrao <- setNames(cj$prof_padrao_m, cj$dataset_id)
g25 <- ease2_grade("25km")
ebu <- read_parquet(file.path(raiz, "camadas/ebu/ebu_celula_25km.parquet"))
niv <- read.csv(file.path(raiz, "camadas/emu/emu_niveis.csv"))                 # bands from the bottom (5500 m) to the surface (0 m)
emu_c <- values(rast(file.path(raiz, "camadas/emu/emu_25km.tif")))
emu_p <- values(rast(file.path(raiz, "camadas/emu/emu_pureza_25km.tif")))
emu_cls <- read.csv(file.path(raiz, "camadas/emu/emu_classes.csv"))
banda_de <- function(z) { z <- pmin(pmax(z, 0), 5500); vapply(z, function(v) which.min(abs(niv$profundidade_m - v)), 1L) }

grupos <- ds |> distinct(grupo, compartimento) |> collect()
multi <- list(); esf <- list(); resumo <- list(); escopo <- list()
for (i in seq_len(nrow(grupos))) {
  g <- grupos$grupo[i]; comp <- grupos$compartimento[i]
  x <- ds |> filter(grupo == g) |> collect()
  # pelagic without month (level C): slice "sem_mes", used only in stage 1 (T4)
  x$fatia <- if (comp == "pelagico") ifelse(is.na(x$mes), "sem_mes", paste0("T", (x$mes - 1) %/% 3 + 1)) else "estatico"
  if (comp == "bentos") {
    k <- match(x$cell_25, ebu$cell_id)
    x$ebu_grupo <- ebu$ebu_grupo[k]; x$ebu_regiao <- ebu$ebu_regiao[k]; x$ebu_unidade <- ebu$ebu_unidade[k]
    x$pureza_grupo <- ebu$pureza_grupo[k]; x$pureza_regiao <- ebu$pureza_regiao[k]; x$pureza_unidade <- ebu$pureza_unidade[k]
    x$situacao_classe <- ifelse(is.na(x$ebu_grupo), "sem_classe", "ok")
  } else {
    zmin <- pmin(x$depth_min, x$depth_max, x$depth, na.rm = TRUE); zmax <- pmax(x$depth_min, x$depth_max, x$depth, na.rm = TRUE)
    pp <- prof_padrao[x$dataset_id]; usa_pp <- is.na(zmin) & !is.na(pp)
    zmin[usa_pp] <- pp[usa_pp]; zmax[usa_pp] <- pp[usa_pp]
    x$prof_padrao_conjunto <- usa_pp
    rc <- do.call(rbind, strsplit(x$cell_25, "_")); cel <- cellFromRowCol(g25, as.integer(rc[, 2]), as.integer(rc[, 3]))
    tem <- !is.na(zmin)
    blo <- rep(NA_integer_, nrow(x)); bhi <- blo
    blo[tem] <- banda_de(zmax[tem]); bhi[tem] <- banda_de(zmin[tem])      # lower band = deeper
    comb <- unique(data.frame(cel = cel[tem], blo = blo[tem], bhi = bhi[tem]))
    res_c <- t(vapply(seq_len(nrow(comb)), function(j) {
      v <- emu_c[comb$cel[j], comb$blo[j]:comb$bhi[j]]; p <- emu_p[comb$cel[j], comb$blo[j]:comb$bhi[j]]
      u <- unique(v[!is.na(v)])
      c(n = length(u), emu = if (length(u) == 1) u else NA, pur = if (length(u) == 1) min(p[!is.na(v)]) / 100 else NA,
        agua = sum(!is.na(v)) / length(v))
    }, c(n = 0, emu = 0, pur = 0, agua = 0)))
    key <- paste(cel, blo, bhi); kc <- paste(comb$cel, comb$blo, comb$bhi); m <- match(key, kc)
    nemu <- res_c[m, "n"]; x$emu <- res_c[m, "emu"]; x$pureza_emu <- res_c[m, "pur"]; agua <- res_c[m, "agua"]
    x$situacao_classe <- ifelse(!tem, "sem_profundidade", ifelse(nemu == 0, "fora_coluna_emu",
                          ifelse(nemu > 1, "multi_emu", ifelse(agua < 1, "parcial_fora_coluna", "ok"))))
    x$situacao_classe[x$situacao_classe == "parcial_fora_coluna"] <- "ok"      # part of the range below the EMU bottom: the water EMU is kept
    x$emu_rotulo <- emu_cls$rotulo[match(x$emu, emu_cls$Cluster37)]
    mm <- x[x$situacao_classe == "multi_emu", c("record_id", "dataset_id", "cell_25", "depth_min", "depth_max", "aphiaid")]
    mm$emus <- vapply(which(x$situacao_classe == "multi_emu"), function(j) {
      v <- emu_c[cel[j], blo[j]:bhi[j]]; paste(unique(v[!is.na(v)]), collapse = ",") }, "")
    multi[[g]] <- mm |> mutate(grupo = g)
  }
  write_dataset(x |> mutate(grupo = g), file.path(raiz, "intermediario/n1_tratado"), partitioning = "grupo", format = "parquet")
  # effort (5.3)
  e <- x |> filter(!dataset_id %in% lista)
  chave <- if (comp == "pelagico") c("cell_25", "fatia", "emu") else c("cell_25", "fatia")
  esf[[g]] <- e |> group_by(across(all_of(chave))) |>
    summarise(eventos = n_distinct(evento), registros = n(), conjuntos = n_distinct(dataset_id),
              ano_min = min(year), ano_max = max(year), anos_1993 = n_distinct(year[year >= 1993]),
              aparelhos = paste(head(names(sort(table(samplingProtocol[!is.na(samplingProtocol) & samplingProtocol != ""]),
                                                decreasing = TRUE)), 3), collapse = " | "), .groups = "drop") |>
    mutate(grupo = g, compartimento = comp, emu = if (comp == "pelagico") emu else NA_real_)
  escopo[[g]] <- x |> distinct(dataset_id) |> mutate(grupo = g)
  resumo[[g]] <- x |> count(situacao_classe) |> mutate(grupo = g, compartimento = comp, pct = round(100 * n / sum(n), 1))
  cat(format(Sys.time()), "-", g, ":", nrow(x), "registros;", paste(resumo[[g]]$situacao_classe, resumo[[g]]$n, sep = "=", collapse = " "), "\n")
  rm(x, e); gc(verbose = FALSE)
}
write_parquet(bind_rows(multi), file.path(raiz, "saida/pelagico_multi_emu.parquet"))
write_parquet(bind_rows(esf), file.path(raiz, "saida/esforco_nivel1.parquet"))
write.csv(bind_rows(resumo), file.path(raiz, "saida/tratamento_n1_resumo.csv"), row.names = FALSE)
cg <- bind_rows(escopo) |> group_by(dataset_id) |> summarise(grupos = paste(sort(grupo), collapse = "; "), n_grupos = n())
cg$titulo <- cj$titulo[match(cg$dataset_id, cj$dataset_id)]; cg$tipo <- cj$tipo[match(cg$dataset_id, cj$dataset_id)]
write.csv(cg, file.path(raiz, "saida/conjunto_grupo.csv"), row.names = FALSE)
cat("conjuntos com registros retidos:", nrow(cg), "\n")
