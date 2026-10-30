# 87_preditores.R  (Phase 8, stage 3.1: sampling units and environmental predictors for level 1)
# Guide v3 (5.2): stages 3 and 4 use only records from 1993 onward with at least the month (certainty A and B).
# Gate 7 (30/10/2026): GF/GDM response = frequency of occurrence (fraction of ALL events in the unit in which the taxon
#   appears, 0 to 1); X (10 sampled events) is kept for the fauna tests (beta_sim). Ocean basin (86b_bacias.R) enters
#   as a covariate (0/1 variables bac_<codigo>). prof = GLORYS depth (prof_glorys) in training and prediction,
#   so that the model uses the same variable that exists in every map cell.
# Units (>= 10 events; X: 10 sampled events, seed 1001, as in stage 1):
#   bentos   = 25 km cell (static; annual predictors);
#   pelagico = cell x EMU x quarter (quarter predictors at the unit's depth), so that the relationships
#              hold month by month in stage 4.
# The class (EBU group/region with purity; EMU) is carried along only to link the results to habitats (3.4); it does not filter units.
# Predictors (camadas/glorys/mes/glorys_<mm>_25km.tif):
#   bentos: prof (GLORYS depth), tfundo_med/min/max, sfundo_med, vfundo_med/max, sst_med, sst_amp,
#           mld_max, gelo_meses (months with concentration >= 0.15)
#   pelagico: z (median of the mid-point of the records' depth band), t_z and s_z (temperature and salinity at z, quarter mean,
#           interpolated between the 14 standard levels), sst, sss, mld (quarter mean), z_mld = log2((z + 1)/(mld + 1)),
#           prof (GLORYS depth), vsup, gelo; both with bac_1 ... bac_8
# Outputs: intermediario/etapa3_n1/<comp>_<grupo>.rds (list: meta + X + predictors), saida/etapa3_n1_unidades.csv (summary)
suppressPackageStartupMessages({library(arrow); library(dplyr); library(terra); library(tidyr)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); source(file.path(raiz, "scripts/funcoes_grade.R"))
dir_t <- file.path(raiz, "intermediario/etapa3_n1"); dir.create(dir_t, showWarnings = FALSE)
EV_MIN <- 10; Z_PADRAO <- c(0, 10, 25, 50, 100, 200, 300, 500, 1000, 1500, 2000, 3000, 4000, 5000)
TRI <- list(T1 = 1:3, T2 = 4:6, T3 = 7:9, T4 = 10:12)
g25 <- ease2_grade("25km")
M <- lapply(1:12, function(m) rast(file.path(raiz, sprintf("camadas/glorys/mes/glorys_%02d_25km.tif", m))))
cel_de <- function(id) { p <- do.call(rbind, strsplit(id, "_")); cellFromRowCol(g25, as.integer(p[, 2]), as.integer(p[, 3])) }
BAC <- values(rast(file.path(raiz, "camadas/bacias_25km.tif")))[, 1]
dummies <- function(cel) { b <- BAC[cel]; m <- sapply(1:8, function(k) as.integer(b %in% k)); colnames(m) <- paste0("bac_", 1:8); as.data.frame(m) }
frequencia <- function(d, un) {   # fraction of the unit's events with the taxon
  f <- d |> distinct(unidade, evento, taxon_unidade) |> count(unidade, taxon_unidade) |> left_join(un |> select(unidade, ev), by = "unidade") |>
    mutate(f = n / ev) |> select(-n, -ev) |> pivot_wider(names_from = taxon_unidade, values_from = f, values_fill = 0) |> arrange(unidade)
  F <- as.matrix(f[, -1]); rownames(F) <- f$unidade; F }
val <- function(cel, camada) sapply(1:12, function(m) values(M[[m]][[camada]])[cel, 1])   # cells x 12 months

incid <- function(d) {   # d: unidade, evento, taxon_unidade -> incidence matrix of the 10-event draw
  set.seed(1001)
  sor <- d |> distinct(unidade, evento) |> group_by(unidade) |> slice_sample(n = EV_MIN) |> ungroup()
  m <- d |> semi_join(sor, by = c("unidade", "evento")) |> distinct(unidade, taxon_unidade) |> mutate(v = 1L) |>
    pivot_wider(names_from = taxon_unidade, values_from = v, values_fill = 0L) |> arrange(unidade)
  X <- as.matrix(m[, -1]); rownames(X) <- m$unidade; X[, colSums(X) > 0, drop = FALSE]
}

resumo <- list()
# ---- benthos ----
b <- open_dataset(file.path(raiz, "saida/tabela_certa_bentos_n1.parquet")) |> filter(year >= 1993, nivel_certeza %in% c("A", "B")) |>
  select(grupo, cell_id, event_id, taxon_unidade, bathymetry_m, bloco_espacial, ebu_grupo, pureza_grupo, ebu_regiao, pureza_regiao) |> collect() |>
  filter(!is.na(taxon_unidade), taxon_unidade != "")
for (g in sort(unique(b$grupo))) {
  d <- b |> filter(grupo == g) |> rename(unidade = cell_id, evento = event_id)
  un <- d |> group_by(unidade) |> summarise(ev = n_distinct(evento), prof = median(bathymetry_m, na.rm = TRUE), bloco = first(bloco_espacial),
    ebu_grupo = first(ebu_grupo), pureza_grupo = first(pureza_grupo), ebu_regiao = first(ebu_regiao), pureza_regiao = first(pureza_regiao), .groups = "drop") |>
    filter(ev >= EV_MIN) |> arrange(unidade)
  X <- incid(d |> filter(unidade %in% un$unidade)); stopifnot(identical(rownames(X), un$unidade))
  cel <- cel_de(un$unidade)
  tf <- val(cel, "tfundo"); sf <- val(cel, "sfundo"); vf <- val(cel, "vfundo"); ss <- val(cel, "thetao_0"); ml <- val(cel, "mld"); ge <- val(cel, "gelo")
  pg <- values(M[[1]][["prof_glorys"]])[cel, 1]
  E <- data.frame(prof = pg, tfundo_med = rowMeans(tf), tfundo_min = apply(tf, 1, min), tfundo_max = apply(tf, 1, max),
    sfundo_med = rowMeans(sf), vfundo_med = rowMeans(vf), vfundo_max = apply(vf, 1, max), sst_med = rowMeans(ss), sst_amp = apply(ss, 1, max) - apply(ss, 1, min),
    mld_max = apply(ml, 1, max), gelo_meses = rowSums(ge >= 0.15, na.rm = TRUE))
  E <- cbind(E, dummies(cel)); F <- frequencia(d |> filter(unidade %in% un$unidade), un); stopifnot(identical(rownames(F), un$unidade))
  ok <- complete.cases(E); xy <- xyFromCell(g25, cel)
  r <- list(comp = "bentos", grupo = g, meta = cbind(un, bacia = BAC[cel], x = xy[, 1], y = xy[, 2])[ok, ], X = X[ok, colSums(X[ok, , drop = FALSE]) > 0, drop = FALSE],
            F = F[ok, colSums(F[ok, , drop = FALSE]) > 0, drop = FALSE], E = E[ok, ])
  saveRDS(r, file.path(dir_t, paste0("bentos_", gsub(" ", "_", g), ".rds")))
  resumo[[length(resumo) + 1]] <- data.frame(comp = "bentos", grupo = g, unidades = sum(ok), sem_ambiente = sum(!ok), taxons = ncol(r$X))
  cat(format(Sys.time()), "bentos", g, ":", sum(ok), "unidades\n")
}
rm(b); gc(verbose = FALSE)

# ---- pelagic ----
p <- open_dataset(file.path(raiz, "saida/tabela_certa_pelagico_n1.parquet")) |> filter(year >= 1993, nivel_certeza %in% c("A", "B"), !is.na(month)) |>
  select(grupo, cell_id, event_id, taxon_unidade, month, depth_min_m, depth_max_m, bathymetry_m, bloco_espacial, emu, pureza_emu) |> collect() |>
  filter(!is.na(taxon_unidade), taxon_unidade != "", !is.na(emu)) |>
  mutate(tri = paste0("T", (month - 1) %/% 3 + 1), zreg = rowMeans(cbind(depth_min_m, depth_max_m), na.rm = TRUE))
interp <- function(perfil, z) vapply(seq_len(nrow(perfil)), function(i) { v <- perfil[i, ]; k <- !is.na(v)
  if (sum(k) == 0) NA_real_ else if (sum(k) == 1) v[k] else approx(Z_PADRAO[k], v[k], xout = min(z[i], max(Z_PADRAO[k])), rule = 2)$y }, 0)
for (g in sort(unique(p$grupo))) {
  d <- p |> filter(grupo == g) |> mutate(unidade = paste(cell_id, emu, tri), evento = event_id)
  un <- d |> group_by(unidade) |> summarise(ev = n_distinct(evento), cell_id = first(cell_id), emu = first(emu), pureza_emu = first(pureza_emu), tri = first(tri),
    z = median(zreg, na.rm = TRUE), prof = median(bathymetry_m, na.rm = TRUE), bloco = first(bloco_espacial), .groups = "drop") |>
    filter(ev >= EV_MIN, !is.na(z)) |> arrange(unidade)
  X <- incid(d |> filter(unidade %in% un$unidade)); stopifnot(identical(rownames(X), un$unidade))
  cel <- cel_de(un$cell_id)
  E <- data.frame(   # filled by quarter below
z = un$z, t_z = NA_real_, s_z = NA_real_, sst = NA_real_, sss = NA_real_, mld = NA_real_, prof = un$prof, vsup = NA_real_, gelo = NA_real_)
  for (tr in names(TRI)) {
    i <- which(un$tri == tr); if (!length(i)) next; ms <- TRI[[tr]]
    med <- function(camada) Reduce(`+`, lapply(ms, function(m) values(M[[m]][[camada]])[cel[i], 1])) / length(ms)
    tz <- sapply(Z_PADRAO, function(zz) med(paste0("thetao_", zz))); sz <- sapply(Z_PADRAO, function(zz) med(paste0("so_", zz)))
    if (length(i) == 1) { tz <- matrix(tz, 1); sz <- matrix(sz, 1) }
    E$t_z[i] <- interp(tz, un$z[i]); E$s_z[i] <- interp(sz, un$z[i]); E$sst[i] <- med("thetao_0"); E$sss[i] <- med("so_0")
    E$mld[i] <- med("mld"); E$vsup[i] <- med("vsup"); E$gelo[i] <- med("gelo")
  }
  E$prof <- values(M[[1]][["prof_glorys"]])[cel, 1]
  E$gelo[is.na(E$gelo)] <- 0; E$z_mld <- log2((E$z + 1) / (E$mld + 1))
  E <- cbind(E, dummies(cel)); F <- frequencia(d |> filter(unidade %in% un$unidade), un); stopifnot(identical(rownames(F), un$unidade))
  ok <- complete.cases(E); xy <- xyFromCell(g25, cel)
  r <- list(comp = "pelagico", grupo = g, meta = cbind(un, bacia = BAC[cel], x = xy[, 1], y = xy[, 2])[ok, ], X = X[ok, colSums(X[ok, , drop = FALSE]) > 0, drop = FALSE],
            F = F[ok, colSums(F[ok, , drop = FALSE]) > 0, drop = FALSE], E = E[ok, ])
  saveRDS(r, file.path(dir_t, paste0("pelagico_", gsub(" ", "_", g), ".rds")))
  resumo[[length(resumo) + 1]] <- data.frame(comp = "pelagico", grupo = g, unidades = sum(ok), sem_ambiente = sum(!ok), taxons = ncol(r$X))
  cat(format(Sys.time()), "pelagico", g, ":", sum(ok), "unidades\n")
}
res <- bind_rows(resumo); write.csv(res, file.path(raiz, "saida/etapa3_n1_unidades.csv"), row.names = FALSE); print(res)
