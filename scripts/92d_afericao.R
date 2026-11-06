# 92d_afericao.R  (GUIA_LIMIARES_PUBLICADOS.md, step L4: published thresholds measured in situ vs. the global model scale)
# 1 For each published threshold with a study box and a model variable: the in-situ value of that variable in the box is compared
#   with GLORYS in the same cells. Only depth has an in-situ measurement in our data: benthic sampling depth of the OBIS events
#   (depth_max_m, else depth_min_m) and the OBIS bathymetry at the point (bathymetry_m, GEBCO-based), against the GLORYS depth of
#   the 25 km cell. Offset = median(in situ - GLORYS) over events. The published value is moved to the GLORYS scale as the median
#   GLORYS depth of the events sampled within +-10% (at least +-50 m) of it (>= 10 events; else not converted). Marginal quantile
#   mapping was tried and rejected: the boxes mix shelf and slope events, so marginal percentiles do not pair the same places.
#   Pelagic z is itself a sampling depth (no conversion); temperature, salinity, MLD and ice have no in-situ values in our data.
# 2 With the converted values, the L2 situation (same tolerance as 92) and the L2b enrichment ratio (same rule as 92b) are
#   recomputed; both versions are reported.
# 3 (GEBCO pending item, A9) the OBIS bathymetry offset is reported next to the sampling-depth offset.
# 4 Variables in the literature that are not in the model, with the number of studies and thresholds.
# Outputs: saida/limiares_afericao.csv (one row per published threshold), saida/limiares_afericao_grupos.csv (per threshold x
#          project group, both versions), saida/limiares_variaveis_candidatas.csv
suppressPackageStartupMessages({library(dplyr); library(terra); library(arrow); library(tidyr)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); source(file.path(raiz, "scripts/funcoes_grade.R"))
TOL <- 0.10; WMIN <- 0.05; NNULL <- 999; EV_MIN <- 30
g25 <- ease2_grade("25km"); pg <- values(rast(file.path(raiz, "camadas/glorys/mes/glorys_01_25km.tif"))[["prof_glorys"]])[, 1]
lp <- read.csv(file.path(raiz, "saida/limiares_publicados_padronizados.csv"), encoding = "UTF-8")
cmp <- read.csv(file.path(raiz, "saida/limiares_comparacao.csv"), encoding = "UTF-8")
en <- read.csv(file.path(raiz, "saida/limiares_enriquecimento.csv"), encoding = "UTF-8")
cel_de <- function(id) { p <- do.call(rbind, strsplit(id, "_")); cellFromRowCol(g25, as.integer(p[, 2]), as.integer(p[, 3])) }

# ---------------- 1 in situ vs GLORYS in the study boxes ----------------
fonte_de <- function(x) if (is.na(x$variavel_modelo)) "variável fora do modelo" else if (is.na(x$caixa_wkt) || x$caixa_wkt == "") "sem caixa (global ou bacia)" else
  if (x$compartimento == "bentos" && x$variavel_modelo == "prof") "profundidade de coleta (OBIS)" else if (x$variavel_modelo == "z") "z já é profundidade de coleta" else "sem medida in situ nos dados"
lp$afericao <- vapply(seq_len(nrow(lp)), function(i) fonte_de(lp[i, ]), "")
quant_map <- function(v, ins, glo) { if (is.na(v)) return(NA_real_)   # conditional median: GLORYS depth of the events sampled near v
  k <- abs(ins - v) <= max(0.10 * v, 50); if (sum(k) < 10) NA_real_ else median(glo[k]) }
af <- bind_rows(lapply(which(lp$afericao == "profundidade de coleta (OBIS)"), function(i) { x <- lp[i, ]; pol <- vect(x$caixa_wkt, crs = "EPSG:4326"); e <- ext(pol)
  d <- open_dataset(file.path(raiz, "saida/tabela_certa_bentos_n1.parquet")) |> filter(lat >= e$ymin, lat <= e$ymax) |>
    select(event_id, lon, lat, depth_min_m, depth_max_m, bathymetry_m, cell_id) |> collect()
  d <- d[is.related(vect(cbind(d$lon, d$lat), crs = "EPSG:4326"), pol, "intersects"), ]
  d <- d |> mutate(insitu = coalesce(depth_max_m, depth_min_m)) |> filter(!is.na(insitu), insitu > 0) |> distinct(event_id, .keep_all = TRUE) |>
    mutate(glorys = pg[cel_de(cell_id)]) |> filter(!is.na(glorys))
  if (nrow(d) < EV_MIN) return(data.frame(id = x$id, eventos = nrow(d)))
  data.frame(id = x$id, eventos = nrow(d), insitu_med = median(d$insitu), glorys_med = median(d$glorys), desloc_med = median(d$insitu - d$glorys),
    desloc_rel = median((d$insitu - d$glorys) / d$glorys), desloc_iqr = IQR(d$insitu - d$glorys), cor_insitu_glorys = cor(d$insitu, d$glorys, method = "spearman"),
    batim_obis_desloc_med = median(d$bathymetry_m - d$glorys, na.rm = TRUE), batim_obis_n = sum(!is.na(d$bathymetry_m)),
    limiar_min_glorys = quant_map(x$limiar_min, d$insitu, d$glorys), limiar_max_glorys = quant_map(x$limiar_max, d$insitu, d$glorys)) }))
afericao <- lp |> select(id, referencia, compartimento, variavel_modelo, limiar_min, limiar_max, bacia, afericao) |> left_join(af, by = "id")
cat(format(Sys.time()), "aferidos:", sum(!is.na(afericao$limiar_min_glorys)), "de", sum(lp$afericao == "profundidade de coleta (OBIS)"), "\n")

# ---------------- 2 L2 and L2b again on the GLORYS scale ----------------
rob <- read.csv(file.path(raiz, "saida/etapa3_n1_limiares_robustos.csv")) |> mutate(id_rob = row_number()) |> separate_rows(lista, sep = ", ") |>
  mutate(grupo = sub(" \\((bentos|pelagico)\\)$", "", lista), comp = sub("^.*\\((bentos|pelagico)\\)$", "\\1", lista))
compara <- function(a, b, de, ate) { if (!length(de)) return("nao_detectado")   # as in 92
  if (any(a * (1 - TOL) <= ate & b * (1 + TOL) >= de)) "confirmado" else if (any(a * (1 - 2 * TOL) <= ate & b * (1 + 2 * TOL) >= de)) "proximo" else "nao_detectado" }
gc_ds <- open_dataset(file.path(raiz, "saida/etapa3_n1_gf_curvas.csv"), format = "csv")
enriq <- function(x, y, a, b) {   # as in 92b
  o <- order(x); x <- x[o]; y <- y[o]; f <- function(v) approx(x, y, v, rule = 2, ties = max)$y
  lo <- min(x); hi <- max(x); w <- b - a; if (w <= 0) { w <- WMIN * (hi - lo); a <- a - w / 2; b <- b + w / 2 }
  if (b < lo || a > hi || w >= hi - lo) return(c(razao = NA, p = NA))
  g <- f(b) - f(a); s <- runif(NNULL, lo, hi - w); nul <- f(s + w) - f(s); c(razao = g / mean(nul), p = (1 + sum(nul >= g)) / (NNULL + 1)) }
set.seed(2026)
conv <- afericao |> filter(!is.na(limiar_min_glorys), !is.na(limiar_max_glorys))
gr <- cmp |> filter(id %in% conv$id) |> select(id, grupo_projeto, variavel_modelo, situacao) |> left_join(conv |> select(id, limiar_min_glorys, limiar_max_glorys), by = "id") |>
  left_join(en |> select(id, grupo_projeto, razao, p), by = c("id", "grupo_projeto"))
gr <- bind_rows(lapply(seq_len(nrow(gr)), function(i) { li <- gr[i, ]; g_ <- li$grupo_projeto; v_ <- li$variavel_modelo
  r <- rob |> filter(comp == "bentos", grupo == g_, preditor == v_)
  cv <- gc_ds |> filter(comp == "bentos", grupo == g_, preditor == v_) |> select(x, y) |> collect()
  s <- if (nrow(cv) && (li$limiar_max_glorys < min(cv$x) || li$limiar_min_glorys > max(cv$x))) "fora_do_alcance" else compara(li$limiar_min_glorys, li$limiar_max_glorys, r$de, r$ate)
  if (li$situacao == "tipo_incompativel") s <- "tipo_incompativel"
  e <- if (nrow(cv) >= 5) enriq(cv$x, cv$y, li$limiar_min_glorys, li$limiar_max_glorys) else c(razao = NA, p = NA)
  cbind(li, situacao_glorys = s, razao_glorys = e[["razao"]], p_glorys = e[["p"]]) }))
gr <- gr |> mutate(confirmado_enriq = !is.na(razao) & razao >= 2 & p <= 0.05, confirmado_enriq_glorys = !is.na(razao_glorys) & razao_glorys >= 2 & p_glorys <= 0.05)

# ---------------- 4 variables outside the model ----------------
cand_v <- lp |> filter(variavel_no_modelo == "nao") |> group_by(variavel_padrao) |>
  summarise(estudos = n_distinct(referencia), limiares = n(), compartimentos = paste(sort(unique(compartimento)), collapse = ", "),
    tipos = paste(sort(unique(tipo_limiar)), collapse = ", "), referencias = paste(sort(unique(referencia)), collapse = "; "), .groups = "drop") |> arrange(desc(estudos))
write.csv(afericao, file.path(raiz, "saida/limiares_afericao.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(gr, file.path(raiz, "saida/limiares_afericao_grupos.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(cand_v, file.path(raiz, "saida/limiares_variaveis_candidatas.csv"), row.names = FALSE, fileEncoding = "UTF-8")
options(width = 220); print(afericao |> filter(!is.na(eventos)) |> select(id, limiar_min, limiar_max, eventos, insitu_med, glorys_med, desloc_med, desloc_rel, batim_obis_desloc_med,
  limiar_min_glorys, limiar_max_glorys) |> mutate(across(where(is.double), ~ round(.x, 2))) |> as.data.frame())
print(table(afericao$afericao)); print(table(original = gr$situacao, glorys = gr$situacao_glorys))
print(gr |> group_by(id) |> summarise(razao_max = max(razao, na.rm = TRUE), razao_glorys_max = max(razao_glorys, na.rm = TRUE), enriq = any(confirmado_enriq), enriq_glorys = any(confirmado_enriq_glorys)) |> as.data.frame())
print(cand_v |> select(-referencias) |> as.data.frame()); cat(format(Sys.time()), "FIM\n")
