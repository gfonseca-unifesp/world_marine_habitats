# 15_gebco.R  (pending items 4.9 and A9; PI approval 13/11/2026: download GEBCO_2026 and run the check)
# GEBCO_2026 Grid (ice surface elevation, 15 arc-second; doi:10.5285/4f68d5c7-45eb-f999-e063-7086abc036fa; public domain) against
# the depth used by the project: prof_glorys = depth of the last valid GLORYS12 level in the 25 km cell (stepwise by construction).
#  1 GEBCO on the project grid: land (elevation >= 0) removed; projected to the 3.125 km EASE grid (median) and aggregated to 25 km
#    (median of the ocean subcells, as for the EBU in 12). Outputs camadas/gebco/gebco_25km.tif (depth > 0, fraction of ocean).
#  2 GLORYS vs GEBCO per 25 km cell of the united benthic map: offset by depth band and EBU group, and the band each cell would get
#    with GEBCO (0-200, 200-2,000, > 2,000 m; same rule as 91). Nothing approved is rebuilt.
#  3 Item 4.9: sample of benthic events (level 1 certainty table): OBIS bathymetry and sampling depth against GEBCO at the point.
#  4 L4 (92d): in the study boxes of the published depth thresholds, the offsets in situ - GEBCO and GEBCO - GLORYS.
# Outputs: saida/gebco_{faixas,ebu,troca_faixa,pontos,afericao}.csv, saida/gebco_diferenca_25km.png
suppressPackageStartupMessages({library(terra); library(dplyr); library(arrow)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); source(file.path(raiz, "scripts/funcoes_grade.R"))
terraOptions(memfrac = 0.6, progress = 0, tempdir = file.path(raiz, "intermediario"))
dg <- file.path(raiz, "camadas/gebco"); nc <- list.files(dg, pattern = "[.]nc$", full.names = TRUE)[1]
g25 <- ease2_grade("25km"); faixa <- function(p) findInterval(p, c(200, 2000), left.open = TRUE) + 1; NOME_FX <- c("0–200 m", "200–2.000 m", "> 2.000 m")
cel_de <- function(id) { p <- do.call(rbind, strsplit(id, "_")); cellFromRowCol(g25, as.integer(p[, 2]), as.integer(p[, 3])) }

# ---------------- 1 GEBCO on the 25 km grid ----------------
f25 <- file.path(dg, "gebco_25km.tif")
if (!file.exists(f25)) {
  r <- rast(nc); cat(format(Sys.time()), "GEBCO", nrow(r), "x", ncol(r), "\n")
  ro <- classify(r, cbind(0, Inf, NA), right = FALSE, filename = file.path(dg, "tmp_oceano.tif"), overwrite = TRUE, datatype = "INT2S")
  g3 <- ease2_grade("3.125km")
  p3 <- project(ro, g3, method = "med", filename = file.path(dg, "tmp_3km.tif"), overwrite = TRUE, datatype = "INT2S"); cat(format(Sys.time()), "projetado a 3,125 km\n")
  med <- aggregate(p3, 8, fun = "median", na.rm = TRUE); oc <- aggregate(!is.na(p3), 8, fun = "mean", na.rm = TRUE)
  if (!compareGeom(med, g25, stopOnError = FALSE)) { med <- resample(med, g25, "near"); oc <- resample(oc, g25, "near") }
  out <- c(-med, oc); names(out) <- c("prof_gebco", "frac_oceano"); writeRaster(out, f25, overwrite = TRUE)
  unlink(file.path(dg, c("tmp_oceano.tif", "tmp_3km.tif"))); cat(format(Sys.time()), "gebco_25km gravado\n")
}
G <- rast(f25); pgb <- values(G[["prof_gebco"]])[, 1]; foc <- values(G[["frac_oceano"]])[, 1]
pgl <- values(rast(file.path(raiz, "camadas/glorys/mes/glorys_01_25km.tif"))[["prof_glorys"]])[, 1]

# ---------------- 2 GLORYS x GEBCO in the cells of the united benthic map ----------------
ub <- values(rast(file.path(raiz, "saida/mapa_unido_bentos_25km.tif"))[["base"]])[, 1]
e25 <- read_parquet(file.path(raiz, "camadas/ebu/ebu_celula_25km.parquet")); e25$cel <- cel_de(e25$cell_id)
k <- which(!is.na(ub) & !is.na(pgl) & !is.na(pgb) & pgb > 0)
C <- data.frame(cel = k, glorys = pgl[k], gebco = pgb[k], oceano = foc[k], grupo_ebu = e25$ebu_grupo[match(k, e25$cel)]) |>
  mutate(dif = glorys - gebco, rel = dif / gebco, fx_glorys = NOME_FX[faixa(glorys)], fx_gebco = NOME_FX[faixa(gebco)])
resumo <- function(d) d |> summarise(celulas = n(), gebco_med = median(gebco), dif_med = median(dif), dif_q1 = quantile(dif, 0.25), dif_q3 = quantile(dif, 0.75),
  rel_med = median(rel), abs_rel_med = median(abs(rel)), frac_rel_10 = mean(abs(rel) > 0.10), frac_rel_25 = mean(abs(rel) > 0.25),
  spearman = cor(glorys, gebco, method = "spearman"), .groups = "drop")
fx <- bind_rows(C |> mutate(fx_glorys = "todas") |> group_by(faixa_glorys = fx_glorys) |> resumo(), C |> group_by(faixa_glorys = fx_glorys) |> resumo()) |>
  mutate(costeiras = c(mean(C$oceano < 0.5), tapply(C$oceano < 0.5, C$fx_glorys, mean)[NOME_FX]))
eb <- C |> group_by(grupo_ebu) |> resumo()
tr <- as.data.frame.matrix(table(factor(C$fx_glorys, NOME_FX), factor(C$fx_gebco, NOME_FX))) |> tibble::rownames_to_column("faixa_glorys")
tr$muda <- rowSums(tr[, NOME_FX]) - diag(as.matrix(tr[, NOME_FX])); tr$frac_muda <- tr$muda / rowSums(tr[, NOME_FX])
write.csv(fx, file.path(raiz, "saida/gebco_faixas.csv"), row.names = FALSE); write.csv(eb, file.path(raiz, "saida/gebco_ebu.csv"), row.names = FALSE)
write.csv(tr, file.path(raiz, "saida/gebco_troca_faixa.csv"), row.names = FALSE)
rd <- rast(g25); v <- rep(NA_real_, ncell(rd)); v[C$cel] <- pmax(pmin(C$rel, 0.5), -0.5); rd[] <- v
png(file.path(raiz, "saida/gebco_diferenca_25km.png"), 2000, 1100, res = 150, bg = "white")
plot(rd, col = hcl.colors(21, "Blue-Red 3"), breaks = seq(-0.5, 0.5, length.out = 22), axes = FALSE, box = FALSE, maxcell = ncell(rd),
  main = "(GLORYS − GEBCO) / GEBCO, profundidade do fundo por célula de 25 km (limitado a ±50%)"); dev.off()
cat(format(Sys.time()), "comparacao por celula feita\n")

# ---------------- 3 item 4.9: OBIS points ----------------
set.seed(49)
P <- open_dataset(file.path(raiz, "saida/tabela_certa_bentos_n1.parquet")) |> select(event_id, lon, lat, depth_min_m, depth_max_m, bathymetry_m, cell_id) |> collect() |>
  distinct(event_id, .keep_all = TRUE)
P <- P[sample(nrow(P), min(nrow(P), 300000)), ]
P$gebco_ponto <- -extract(rast(nc), cbind(P$lon, P$lat))[, 1]; P$glorys <- pgl[cel_de(P$cell_id)]
P <- P |> mutate(insitu = coalesce(depth_max_m, depth_min_m), gebco_ponto = ifelse(gebco_ponto > 0, gebco_ponto, NA))
pt <- function(x, y, nome) { k <- !is.na(x) & !is.na(y) & y > 0; d <- x[k] - y[k]
  data.frame(comparacao = nome, eventos = sum(k), dif_med = median(d), rel_med = median(d / y[k]), abs_rel_med = median(abs(d / y[k])),
    frac_rel_10 = mean(abs(d / y[k]) > 0.10), frac_rel_25 = mean(abs(d / y[k]) > 0.25), spearman = cor(x[k], y[k], method = "spearman"),
    abaixo_do_fundo = mean(d > pmax(0.1 * y[k], 50))) }
pts <- bind_rows(pt(P$bathymetry_m, P$gebco_ponto, "batimetria OBIS − GEBCO no ponto"), pt(P$insitu, P$gebco_ponto, "profundidade de coleta − GEBCO no ponto"),
  pt(P$insitu, P$glorys, "profundidade de coleta − GLORYS na célula"), pt(P$glorys, P$gebco_ponto, "GLORYS na célula − GEBCO no ponto"))
write.csv(pts, file.path(raiz, "saida/gebco_pontos.csv"), row.names = FALSE)
cat(format(Sys.time()), "pontos feitos\n")

# ---------------- 4 L4: study boxes of the published depth thresholds ----------------
lp <- read.csv(file.path(raiz, "saida/limiares_publicados_padronizados.csv"), encoding = "UTF-8")
ids <- read.csv(file.path(raiz, "saida/limiares_afericao.csv"), encoding = "UTF-8") |> filter(!is.na(eventos), eventos >= 30) |> pull(id)
B <- open_dataset(file.path(raiz, "saida/tabela_certa_bentos_n1.parquet")) |> select(event_id, lon, lat, depth_min_m, depth_max_m, cell_id) |> collect() |> distinct(event_id, .keep_all = TRUE)
af <- bind_rows(lapply(ids, function(i) { x <- lp[lp$id == i, ][1, ]; pol <- vect(x$caixa_wkt, crs = "EPSG:4326"); e <- ext(pol)
  d <- B |> filter(lat >= e$ymin, lat <= e$ymax); d <- d[is.related(vect(cbind(d$lon, d$lat), crs = "EPSG:4326"), pol, "intersects"), ]
  d <- d |> mutate(insitu = coalesce(depth_max_m, depth_min_m), glorys = pgl[cel_de(cell_id)]) |> filter(!is.na(insitu), insitu > 0, !is.na(glorys))
  d$gebco <- -extract(rast(nc), cbind(d$lon, d$lat))[, 1]; d <- d |> filter(gebco > 0)
  data.frame(id = i, referencia = x$referencia, eventos = nrow(d), insitu_menos_glorys = median(d$insitu - d$glorys), insitu_menos_gebco = median(d$insitu - d$gebco),
    gebco_menos_glorys = median(d$gebco - d$glorys), rel_gebco_glorys = median((d$gebco - d$glorys) / d$gebco)) }))
write.csv(af, file.path(raiz, "saida/gebco_afericao.csv"), row.names = FALSE)
options(width = 200); r3 <- function(d) as.data.frame(d |> mutate(across(where(is.double), ~ round(.x, 3))))
print(r3(fx)); print(r3(eb)); print(tr); print(r3(pts)); print(r3(af)); cat(format(Sys.time()), "FIM\n")
