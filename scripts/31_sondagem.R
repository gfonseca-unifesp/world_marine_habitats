# 31_sondagem.R  (tasks 3.1 to 3.3)
# Level 1 survey (scoping analysis) on the PRESENCES that pass the cheap filters (valid coordinate, year present,
# cell with class). Numbers are an upper bound: the full curation has not run yet.
# For 25, 50 and 100 km, by group: purity, effort (distinct events), testable classes, boundary pairs,
# temporal coverage, time series and dataset dominance. Applies the resolution choice rule (E3.3)
# and estimates the autocorrelation range of the composition (3.3).
# Output: intermediario/sondagem_celulas.parquet (cell x group x res), saida/sondagem_resolucao.csv,
#        saida/sondagem_curva_esforco.csv, saida/sondagem_autocorrelacao.csv
suppressPackageStartupMessages({library(arrow); library(dplyr); library(terra)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
source(file.path(raiz, "scripts/funcoes_grade.R"))
set.seed(42)

ds  <- open_dataset(file.path(raiz, "intermediario/nivel1_extracao"))
cj  <- read.csv(file.path(raiz, "intermediario/conjunto_tipo_escopo.csv"))
grp <- read.csv(file.path(raiz, "saida/nivel1_definicao_grupos.csv")) |> distinct(grupo, compartimento)
g25 <- ease2_grade("25km")
RES <- c("25km" = 1, "50km" = 2, "100km" = 4)          # factor applied to the 25 km row/column

# ---- class layers per cell ----
ebu <- lapply(names(RES), function(r) read_parquet(file.path(raiz, sprintf("camadas/ebu/ebu_celula_%s.parquet", r))))
names(ebu) <- names(RES)
niv <- read.csv(file.path(raiz, "camadas/emu/emu_niveis.csv"))         # band -> depth (bands from the bottom to the surface)
emu_m <- lapply(names(RES), function(r) list(c = values(rast(file.path(raiz, sprintf("camadas/emu/emu_%s.tif", r)))),
                                              p = values(rast(file.path(raiz, sprintf("camadas/emu/emu_pureza_%s.tif", r))))))
names(emu_m) <- names(RES)
banda_de <- function(z) { z <- pmin(pmax(z, 0), 5500); sapply(z, function(v) which.min(abs(niv$profundidade_m - v))) }

ms <- function(x) as.POSIXct(as.numeric(x) / 1000, origin = "1970-01-01", tz = "UTC")
prof_padrao <- setNames(cj$prof_padrao_m, cj$dataset_id)

celulas <- list(); resumo <- list(); curva <- list(); autoc <- list()
for (g in grp$grupo) {
  comp <- grp$compartimento[grp$grupo == g]
  x <- ds |> filter(grupo == g, !(absence %in% TRUE)) |>
    select(dataset_id, lon, lat, date_start, date_end, year, depth_min, depth_max, depth, eventID, family) |> collect()
  n0 <- nrow(x)
  x <- x |> filter(!is.na(lon), !is.na(lat), abs(lat) <= 90, abs(lon) <= 180, !(lon == 0 & lat == 0), !is.na(year))
  # 25 km cell and derived cells
  xy <- project(cbind(x$lon, x$lat), from = "EPSG:4326", to = "EPSG:6933")
  cel <- cellFromXY(g25, xy); ok <- !is.na(cel); x <- x[ok, ]; cel <- cel[ok]
  rc <- rowColFromCell(g25, cel); x$r25 <- rc[, 1]; x$c25 <- rc[, 2]
  # event: eventID; without it, dataset + date + coordinate
  x$evento <- ifelse(!is.na(x$eventID) & x$eventID != "", paste(x$dataset_id, x$eventID),
                     paste(x$dataset_id, x$date_start, round(x$lon, 4), round(x$lat, 4)))
  di <- ms(x$date_start); dfim <- ms(x$date_end)
  mes_ok <- !is.na(di) & format(di, "%Y-%m") == format(dfim, "%Y-%m")
  x$mes <- ifelse(mes_ok, as.integer(format(di, "%m")), NA); x$estacao <- (x$mes - 1) %/% 3 + 1
  if (comp == "pelagico") {
    z <- ifelse(!is.na(x$depth_min) & !is.na(x$depth_max), (x$depth_min + x$depth_max) / 2,
         ifelse(!is.na(x$depth), x$depth, ifelse(!is.na(x$depth_min), x$depth_min, x$depth_max)))
    z <- ifelse(is.na(z), prof_padrao[x$dataset_id], z)
    x$banda <- ifelse(is.na(z), NA, banda_de(ifelse(is.na(z), 0, z)))
  }
  for (r in names(RES)) {
    f <- RES[[r]]; x$row <- (x$r25 - 1) %/% f + 1; x$col <- (x$c25 - 1) %/% f + 1
    gr <- ease2_grade(r); x$cel <- cellFromRowCol(gr, x$row, x$col); x$cell_id <- paste(r, x$row, x$col, sep = "_")
    if (comp == "bentos") {
      e <- ebu[[r]]; k <- match(x$cell_id, e$cell_id)
      niveis <- list(grupo = c("ebu_grupo", "pureza_grupo"), regiao = c("ebu_regiao", "pureza_regiao"),
                     unidade = c("ebu_unidade", "pureza_unidade"))
      cls <- lapply(niveis, function(n) data.frame(classe = e[[n[1]]][k], pureza = e[[n[2]]][k]))
    } else {
      k <- cbind(x$cel, x$banda); valido <- !is.na(x$banda)
      cl <- rep(NA_integer_, nrow(x)); pu <- rep(NA_real_, nrow(x))
      cl[valido] <- emu_m[[r]]$c[k[valido, , drop = FALSE]]; pu[valido] <- emu_m[[r]]$p[k[valido, , drop = FALSE]] / 100
      cls <- list(emu = data.frame(classe = as.character(cl), pureza = pu))
    }
    for (nv in names(cls)) {
      y <- cbind(x[, c("dataset_id", "cell_id", "row", "col", "evento", "year", "mes", "estacao", "family")], cls[[nv]])
      y <- y[!is.na(y$classe), ]
      if (!nrow(y)) next
      cc <- y |> group_by(cell_id, row, col, classe) |>
        summarise(pureza = first(pureza), eventos = n_distinct(evento), registros = n(), conjuntos = n_distinct(dataset_id),
                  meses = n_distinct(mes, na.rm = TRUE), estacoes = n_distinct(estacao[year >= 1993], na.rm = TRUE),
                  anos93 = n_distinct(year[year >= 1993]), familias = n_distinct(family[!is.na(family) & family != ""]),
                  .groups = "drop") |>
        group_by(cell_id) |> slice_max(eventos, n = 1, with_ties = FALSE) |> ungroup() |>   # (pelagic: most sampled stratum)
        mutate(serie = estacoes >= 3 & anos93 >= 5, pura = pureza >= 0.8)
      # boundary pairs (4-neighbourhood) between cells with data of different classes
      key <- paste(cc$row, cc$col); cl_de <- setNames(cc$classe, key); ev_de <- setNames(cc$eventos, key)
      viz <- rbind(data.frame(a = key, b = paste(cc$row, cc$col + 1)), data.frame(a = key, b = paste(cc$row + 1, cc$col)))
      viz <- viz[viz$b %in% key, ]
      dif <- cl_de[viz$a] != cl_de[viz$b]
      par1 <- sum(dif); par5 <- sum(dif & ev_de[viz$a] >= 5 & ev_de[viz$b] >= 5)
      # testable classes: >= 5 pure cells with >= 10 events
      cl_ok <- cc |> filter(pura, eventos >= 10) |> count(classe) |> filter(n >= 5)
      top <- sort(table(y$dataset_id), decreasing = TRUE)
      resumo[[length(resumo) + 1]] <- data.frame(
        grupo = g, compartimento = comp, nivel = nv, res = r, celulas = nrow(cc),
        mediana_pureza = round(median(cc$pureza), 3), pct_puras = round(100 * mean(cc$pura), 1),
        cel_ev1 = sum(cc$eventos >= 1), cel_ev3 = sum(cc$eventos >= 3), cel_ev5 = sum(cc$eventos >= 5),
        cel_ev10 = sum(cc$eventos >= 10), classes_com_dados = n_distinct(cc$classe), classes_testaveis = nrow(cl_ok),
        pares_fronteira = par1, pares_fronteira_5ev = par5, celulas_serie = sum(cc$serie),
        pct_maior_conjunto = round(100 * as.numeric(top[1]) / nrow(y), 1),
        passa_a = median(cc$pureza) >= 0.8, passa_b = nrow(cl_ok) >= 2)
      celulas[[length(celulas) + 1]] <- cc |> mutate(grupo = g, compartimento = comp, nivel = nv, res = r)
      if (r == "25km") curva[[length(curva) + 1]] <- data.frame(grupo = g, nivel = nv, min_eventos = c(1, 2, 3, 5, 10, 20, 50),
        celulas = sapply(c(1, 2, 3, 5, 10, 20, 50), function(m) sum(cc$eventos >= m)))
    }
  }
  # composition autocorrelation (25 km, families, cells with >= 10 events): Jaccard x distance
  a <- x |> filter(!is.na(family), family != "") |> mutate(cell_id = paste(r25, c25)) |>
    group_by(cell_id, r25, c25) |> filter(n_distinct(evento) >= 10) |> ungroup()
  if (n_distinct(a$cell_id) >= 50) {
    cs <- sample(unique(a$cell_id), min(1500, n_distinct(a$cell_id)))
    fam <- split(a$family[a$cell_id %in% cs], a$cell_id[a$cell_id %in% cs]); fam <- lapply(fam, unique)
    pos <- a |> filter(cell_id %in% cs) |> distinct(cell_id, r25, c25)
    pos <- pos[match(names(fam), pos$cell_id), ]
    xyc <- project(xyFromCell(g25, cellFromRowCol(g25, pos$r25, pos$c25)), from = "EPSG:6933", to = "EPSG:4326")
    pr <- t(combn(length(fam), 2)); if (nrow(pr) > 3e5) pr <- pr[sample(nrow(pr), 3e5), ]
    dk <- geosphere_dist <- terra::distance(xyc[pr[, 1], , drop = FALSE], xyc[pr[, 2], , drop = FALSE], lonlat = TRUE, pairwise = TRUE) / 1000
    jac <- mapply(function(i, j) { u <- length(union(fam[[i]], fam[[j]])); 1 - length(intersect(fam[[i]], fam[[j]])) / u }, pr[, 1], pr[, 2])
    br <- c(0, 25, 50, 100, 150, 200, 300, 400, 500, 750, 1000, 1500, 2000, 3000, 5000, 20000)
    bin <- cut(dk, br, include.lowest = TRUE)
    v <- data.frame(grupo = g, faixa_km = levels(bin), dist_media = tapply(dk, bin, mean), jaccard = tapply(jac, bin, mean),
                    pares = as.integer(table(bin)))
    patamar <- mean(jac[dk > 3000]); limiar <- 0.95 * patamar
    alc <- v$dist_media[which(v$jaccard >= limiar & v$pares >= 30)[1]]
    v$alcance_km <- round(alc); v$patamar <- round(patamar, 3)
    autoc[[g]] <- v
  }
  cat(format(Sys.time()), "-", g, ": presencas", n0, "-> usadas", nrow(x), "\n")
  rm(x); gc(verbose = FALSE)
}
write_parquet(bind_rows(celulas), file.path(raiz, "intermediario/sondagem_celulas.parquet"))
res <- bind_rows(resumo)
ord <- c("25km" = 1, "50km" = 2, "100km" = 3)
escolha <- res |> filter(passa_a & passa_b) |> group_by(grupo, nivel) |> slice_min(ord[res], n = 1) |>
  select(grupo, nivel, res_escolhida = res)
res <- res |> left_join(escolha, by = c("grupo", "nivel"))
write.csv(res, file.path(raiz, "saida/sondagem_resolucao.csv"), row.names = FALSE)
write.csv(bind_rows(curva), file.path(raiz, "saida/sondagem_curva_esforco.csv"), row.names = FALSE)
write.csv(bind_rows(autoc), file.path(raiz, "saida/sondagem_autocorrelacao.csv"), row.names = FALSE)
print(res |> filter(res == "25km") |> select(grupo, nivel, celulas, mediana_pureza, cel_ev10, classes_testaveis, celulas_serie, pct_maior_conjunto))
print(distinct(res, grupo, nivel, res_escolhida))
print(bind_rows(autoc) |> distinct(grupo, alcance_km, patamar))
