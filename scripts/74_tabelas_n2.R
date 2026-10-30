# 74_tabelas_n2.R  (task 7.4: tables, resolution, validation and testability per region)
# Classes at the three level 2 resolutions (25, 12.5 and 6.25 km):
#   bentos   -> EBU (group, region, unit) and purities, per cell;
#   pelagico -> EMU by native 3D (D4.1; same as 50_tratamento_n1.R). The EMU (27 km) uses the 12.5 km grid for
#               12.5 and 6.25 km (finer adds no EMU information).
# Level 2 time slices: pelagic in climatological month (01 to 12; "sem_mes" only in stage 1); benthos in quarter
#   (T1 to T4; "sem_mes") and static.
# Resolution rule (3.2) per region x class level x ecological group: the finest at which (a) the median purity
#   in cells with data is >= 0.8 and (b) >= 2 classes with >= 5 pure cells of >= 10 events. Records only enter
#   the resolution their position allows (res_max, 7.3).
# Outputs: saida/tabela_certa_bentos_n2.parquet, saida/tabela_certa_pelagico_n2.parquet, saida/resolucao_n2.csv,
#         saida/testabilidade_n2.csv, saida/validacao_n2.csv
suppressPackageStartupMessages({library(arrow); library(dplyr); library(terra)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
source(file.path(raiz, "scripts/funcoes_grade.R"))
ds <- open_dataset(file.path(raiz, "intermediario/n2_curado"), unify_schemas = TRUE)
RES <- c("25km", "12.5km", "6.25km"); ORD <- c("25km" = 1, "12.5km" = 2, "6.25km" = 3)   # ORD: larger = finer
ebu <- lapply(RES, function(r) read_parquet(file.path(raiz, sprintf("camadas/ebu/ebu_celula_%s.parquet", r)))); names(ebu) <- RES
niv <- read.csv(file.path(raiz, "camadas/emu/emu_niveis.csv"))
emu <- list("25km" = list(c = values(rast(file.path(raiz, "camadas/emu/emu_25km.tif"))), p = values(rast(file.path(raiz, "camadas/emu/emu_pureza_25km.tif")))),
            "12.5km" = list(c = values(rast(file.path(raiz, "camadas/emu/emu_12.5km.tif"))), p = values(rast(file.path(raiz, "camadas/emu/emu_pureza_12.5km.tif")))))
cj <- read.csv(file.path(raiz, "intermediario/conjunto_tipo_escopo.csv")); prof_padrao <- setNames(cj$prof_padrao_m, cj$dataset_id)
banda_de <- function(z) { z <- pmin(pmax(z, 0), 5500); vapply(z, function(v) which.min(abs(niv$profundidade_m - v)), 1L) }
rc_de <- function(id) { p <- do.call(rbind, strsplit(id, "_")); cbind(as.integer(p[, 2]), as.integer(p[, 3])) }

emu_3d <- function(x, res) {
  g <- ease2_grade(res); E <- emu[[res]]
  zmin <- pmin(x$depth_min, x$depth_max, x$depth, na.rm = TRUE); zmax <- pmax(x$depth_min, x$depth_max, x$depth, na.rm = TRUE)
  pp <- prof_padrao[x$dataset_id]; u <- is.na(zmin) & !is.na(pp); zmin[u] <- pp[u]; zmax[u] <- pp[u]
  rc <- rc_de(x[[paste0("cell_", sub("km", "", res))]]); cel <- cellFromRowCol(g, rc[, 1], rc[, 2]); tem <- !is.na(zmin)
  blo <- rep(NA_integer_, nrow(x)); bhi <- blo; blo[tem] <- banda_de(zmax[tem]); bhi[tem] <- banda_de(zmin[tem])
  comb <- unique(data.frame(cel = cel[tem], blo = blo[tem], bhi = bhi[tem]))
  rr <- t(vapply(seq_len(nrow(comb)), function(j) { v <- E$c[comb$cel[j], comb$blo[j]:comb$bhi[j]]; p <- E$p[comb$cel[j], comb$blo[j]:comb$bhi[j]]
    u <- unique(v[!is.na(v)]); c(n = length(u), emu = if (length(u) == 1) u else NA, pur = if (length(u) == 1) min(p[!is.na(v)]) / 100 else NA) },
    c(n = 0, emu = 0, pur = 0)))
  m <- match(paste(cel, blo, bhi), paste(comb$cel, comb$blo, comb$bhi))
  data.frame(emu = rr[m, "emu"], pureza = rr[m, "pur"], situacao = ifelse(!tem, "sem_profundidade", ifelse(rr[m, "n"] == 0, "fora_coluna_emu",
             ifelse(rr[m, "n"] > 1, "multi_emu", "ok"))))
}

regioes <- ds |> distinct(regiao) |> collect() |> pull()
B <- list(); P <- list()
for (rg in regioes) {
  x <- ds |> filter(regiao == rg, !so_esforco) |> collect()
  b <- x[x$compartimento == "bentos", ]; p <- x[x$compartimento == "pelagico", ]
  for (r in RES) {
    k <- match(b[[paste0("cell_", sub("km", "", r))]], ebu[[r]]$cell_id); e <- ebu[[r]]
    sfx <- sub("km", "", r)
    b[[paste0("ebu_grupo_", sfx)]] <- e$ebu_grupo[k]; b[[paste0("pureza_grupo_", sfx)]] <- e$pureza_grupo[k]
    b[[paste0("ebu_regiao_", sfx)]] <- e$ebu_regiao[k]; b[[paste0("pureza_regiao_", sfx)]] <- e$pureza_regiao[k]
    b[[paste0("ebu_unidade_", sfx)]] <- e$ebu_unidade[k]; b[[paste0("pureza_unidade_", sfx)]] <- e$pureza_unidade[k]
  }
  if (nrow(p)) for (r in c("25km", "12.5km")) { em <- emu_3d(p, r); sfx <- sub("km", "", r)
    p[[paste0("emu_", sfx)]] <- em$emu; p[[paste0("pureza_emu_", sfx)]] <- em$pureza; p[[paste0("situacao_emu_", sfx)]] <- em$situacao }
  if (nrow(p)) { p$emu_6.25 <- p$emu_12.5; p$pureza_emu_6.25 <- p$pureza_emu_12.5; p$situacao_emu_6.25 <- p$situacao_emu_12.5 }
  b$fatia_trimestre <- ifelse(is.na(b$mes), "sem_mes", paste0("T", (b$mes - 1) %/% 3 + 1)); b$fatia_estatica <- "estatico"
  if (nrow(p)) p$fatia_mes <- ifelse(is.na(p$mes), "sem_mes", sprintf("M%02d", p$mes))
  B[[rg]] <- b |> select(-any_of(c("flags"))); P[[rg]] <- p |> select(-any_of(c("flags")))
  cat(format(Sys.time()), "-", rg, ": bentos", nrow(b), "| pelagico", nrow(p), "\n"); rm(x, b, p); gc(verbose = FALSE)
}
bent <- bind_rows(B); pel <- bind_rows(P)
bent <- bent[!is.na(bent$ebu_grupo_25), ]                        # benthos without an EBU class even at 25 km is excluded
write_parquet(bent, file.path(raiz, "saida/tabela_certa_bentos_n2.parquet")); write_parquet(pel, file.path(raiz, "saida/tabela_certa_pelagico_n2.parquet"))

# ---- resolution rule (3.2) and testability, per region x level x ecological group x slice ----
avalia <- function(d, classe, pureza, res, fatia) {
  d <- d[ORD[d$res_max] >= ORD[res] & !is.na(d[[classe]]), ]
  if (!nrow(d)) return(NULL)
  cel <- paste0("cell_", sub("km", "", res))
  cc <- d |> group_by(regiao, grupo_ecologico, fatia = .data[[fatia]], cell = .data[[cel]]) |>
    summarise(classe = first(.data[[classe]]), pureza = first(.data[[pureza]]), ev = n_distinct(evento), taxons = n_distinct(taxon_unidade), .groups = "drop")
  # eligible class: >= 5 pure cells (>= 0.8) with >= 10 events
  apt <- cc |> filter(pureza >= 0.8, ev >= 10) |> count(regiao, grupo_ecologico, fatia, classe) |> filter(n >= 5) |>
    count(regiao, grupo_ecologico, fatia, name = "classes_aptas")
  cc |> group_by(regiao, grupo_ecologico, fatia) |>
    summarise(celulas = n(), mediana_pureza = median(pureza, na.rm = TRUE), classes_com_dados = n_distinct(classe),
              eventos = sum(ev), taxons_max = max(taxons), .groups = "drop") |>
    left_join(apt, by = c("regiao", "grupo_ecologico", "fatia")) |>
    mutate(classes_aptas = coalesce(classes_aptas, 0L), res = res, nivel = sub("_[0-9.]+$", "", classe))
}
linhas <- list()
for (r in RES) { s <- sub("km", "", r)
  for (nv in c("grupo", "regiao", "unidade")) linhas[[length(linhas) + 1]] <- bind_rows(
    avalia(bent, paste0("ebu_", nv, "_", s), paste0("pureza_", nv, "_", s), r, "fatia_estatica"),
    avalia(bent, paste0("ebu_", nv, "_", s), paste0("pureza_", nv, "_", s), r, "fatia_trimestre"))
  linhas[[length(linhas) + 1]] <- bind_rows(avalia(pel |> filter(.data[[paste0("situacao_emu_", s)]] == "ok"), paste0("emu_", s), paste0("pureza_emu_", s), r, "fatia_mes"))
}
tst <- bind_rows(linhas) |> mutate(passa_a = mediana_pureza >= 0.8, passa_b = classes_aptas >= 2, contraste_viavel = ifelse(passa_a & passa_b, "sim", "não"))
esc <- tst |> filter(passa_a, passa_b) |> group_by(regiao, nivel, grupo_ecologico, fatia) |> slice_max(ORD[res], n = 1, with_ties = FALSE) |>
  select(regiao, nivel, grupo_ecologico, fatia, res_escolhida = res)
tst <- tst |> left_join(esc, by = c("regiao", "nivel", "grupo_ecologico", "fatia"))
write.csv(tst, file.path(raiz, "saida/testabilidade_n2.csv"), row.names = FALSE)
write.csv(esc, file.path(raiz, "saida/resolucao_n2.csv"), row.names = FALSE)

# ---- validation ----
V <- list(); teste <- function(t, ok, det = "") V[[length(V) + 1]] <<- data.frame(teste = t, resultado = ifelse(ok, "passou", "FALHOU"), detalhe = det)
teste("record_id único (bentos)", !anyDuplicated(bent$record_id), paste(nrow(bent), "registros"))
teste("record_id único (pelágico)", !anyDuplicated(pel$record_id), paste(nrow(pel), "registros"))
teste("nenhum registro nas duas tabelas", !any(bent$record_id %in% pel$record_id))
reg <- read.csv(file.path(raiz, "entrada/regioes_nivel2.csv")); W <- 80                     # 80 cells of 6.25 km per window
jan <- function(d) { rc <- rc_de(d$cell_6.25); paste((rc[, 1] - 1) %/% W + 1, (rc[, 2] - 1) %/% W + 1, sep = "_") }
teste("todo registro cai na janela da sua região", all(jan(bent) == reg$janela[match(bent$regiao, reg$nome)]) &&
        all(jan(pel) == reg$janela[match(pel$regiao, reg$nome)]))
g6 <- ease2_grade("6.25km"); cel6 <- function(d) { c <- cellFromXY(g6, project(cbind(d$lon, d$lat), from = "EPSG:4326", to = "EPSG:6933")); rc <- rowColFromCell(g6, c); paste("6.25km", rc[, 1], rc[, 2], sep = "_") }
teste("lat/lon cai na célula de 6,25 km", all(cel6(bent) == bent$cell_6.25) && all(cel6(pel) == pel$cell_6.25))
teste("compartimento coerente com a tabela", all(bent$compartimento == "bentos") && all(pel$compartimento == "pelagico"))
teste("grupo ecológico preenchido", !any(is.na(c(bent$grupo_ecologico, pel$grupo_ecologico))))
teste("resolução máxima informada", !any(is.na(c(bent$res_max, pel$res_max))))
teste("nenhum táxon do Anexo I da UNCLOS", !any(c(bent$genus, pel$genus) %in% c("Thunnus", "Katsuwonus", "Xiphias", "Makaira", "Istiophorus")))
teste("nenhum Molidae", !any(c(bent$family, pel$family) %in% "Molidae"))
cas2 <- read.csv(file.path(raiz, "saida/cascata_filtros_n2.csv"))
esp <- sum(cas2$restam[cas2$filtro == "7.3 compartimento"])
# count from each file's metadata (parts may have different types in empty columns)
linhas_em <- function(d) sum(vapply(list.files(d, recursive = TRUE, full.names = TRUE, pattern = "parquet$"),
                                    function(f) ParquetFileReader$create(f)$num_rows, 0))
gravados <- linhas_em(file.path(raiz, "intermediario/n2_curado"))
teste("registros gravados = restantes na cascata", gravados == esp, paste(gravados, "gravados;", esp, "na cascata"))
ext <- linhas_em(file.path(raiz, "intermediario/nivel2_extracao"))
qn <- linhas_em(file.path(raiz, "intermediario/quarentena/n2")); an <- linhas_em(file.path(raiz, "intermediario/ausencias/n2"))
teste("cascata fecha: extraídos = retidos + quarentena + ausências", ext == gravados + qn + an, paste(ext, "=", gravados, "+", qn, "+", an))
val <- bind_rows(V); write.csv(val, file.path(raiz, "saida/validacao_n2.csv"), row.names = FALSE); print(val)
print(tst |> filter(fatia %in% c("estatico")) |> group_by(regiao, nivel) |>
        summarise(grupos = n_distinct(grupo_ecologico), viaveis = n_distinct(grupo_ecologico[contraste_viavel == "sim"]), .groups = "drop") |> as.data.frame())
