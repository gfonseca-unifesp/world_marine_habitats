# 89c_bentos_sazonal.R  (Phase 8, stage 3; before Gate 9: does benthos need dynamic habitats?)
# The variables that separate the new benthic habitats (HB) are annual summaries of monthly variables (sst_med, sst_amp,
# mld_max, tfundo_*, vfundo_*, gelo_meses); only depth is static. Three tests:
#  A within-year variation: median monthly range of each variable per depth band, and fraction of cells in which the
#    monthly value crosses a robust benthic threshold (89_ligacao_habitats.R) over the year;
#  B HB stability: in each month, the level variables (tfundo_med, sst_med, vfundo_med) take the month's value and the
#    regime ones (min, max, amplitude, mld_max, gelo_meses) stay annual; the cell goes to the nearest group within the band
#    (same centroids and merges as 89b). Fraction of cells that change HB in at least one month;
#  C fauna: unit = cell x quarter (>= 10 events, 10 sampled), only cells with >= 2 quarters; PERMANOVA of quarter
#    within cell (sqrt(beta_sim), permutation within cell), omega2 and omega2 with equal effort (10 units per
#    quarter). Compared with the spatial omega2 of the HB (89b). Limitation: quarter is confounded with year and gear.
# Outputs: saida/etapa3_n1_sazonal_{variacao,limiares,estabilidade,fauna}.csv, saida/mapa_bentos_meses_fora.png
suppressPackageStartupMessages({library(dplyr); library(terra); library(arrow); library(parallel); library(tidyr)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); source(file.path(raiz, "scripts/funcoes_grade.R"))
dir_t <- file.path(raiz, "intermediario/etapa3_n1"); KF <- 10; FX <- rep(1:3, each = KF); NOME_FX <- c("0–200 m", "200–2.000 m", "> 2.000 m")
faixa <- function(prof) findInterval(prof, c(200, 2000), left.open = TRUE) + 1
g25 <- ease2_grade("25km"); M <- lapply(1:12, function(m) rast(file.path(raiz, sprintf("camadas/glorys/mes/glorys_%02d_25km.tif", m))))
vm <- function(camada) sapply(1:12, function(k) values(M[[k]][[camada]])[, 1])
curvas <- read.csv(file.path(raiz, "saida/etapa3_n1_gf_curvas.csv"))
transf <- function(E, cp) { cv <- curvas[curvas$comp == cp, ]; vs <- intersect(unique(cv$preditor), names(E))   # as in 89b
  sapply(vs, function(v) { s <- 0; for (g in unique(cv$grupo[cv$preditor == v])) { c1 <- cv[cv$grupo == g & cv$preditor == v, ]
    s <- s + approx(c1$x, c1$y, E[[v]], rule = 2, ties = mean)$y }; s }) }
mais_proximo <- function(T, C) { d <- sapply(seq_len(nrow(C)), function(k) colSums((t(T) - C[k, ])^2)); max.col(-d, ties.method = "first") }
atribui_b <- function(T, prof, C) { f <- faixa(prof); k <- integer(nrow(T))
  for (j in 1:3) { i <- which(f == j); if (length(i)) k[i] <- (j - 1) * KF + mais_proximo(T[i, , drop = FALSE], C[FX == j, , drop = FALSE]) }; k }

# ---- benthos grid (as in 89b) ----
e25 <- read_parquet(file.path(raiz, "camadas/ebu/ebu_celula_25km.parquet"))
p <- do.call(rbind, strsplit(e25$cell_id, "_")); cel_b <- cellFromRowCol(g25, as.integer(p[, 2]), as.integer(p[, 3]))
MES <- list(tfundo = vm("tfundo")[cel_b, ], vfundo = vm("vfundo")[cel_b, ], sst = vm("thetao_0")[cel_b, ], mld = vm("mld")[cel_b, ], gelo = vm("gelo")[cel_b, ])
MES$gelo[is.na(MES$gelo)] <- 0
GB <- data.frame(cel = cel_b, prof = values(M[[1]][["prof_glorys"]])[cel_b, 1], tfundo_med = rowMeans(MES$tfundo), tfundo_min = apply(MES$tfundo, 1, min),
  tfundo_max = apply(MES$tfundo, 1, max), vfundo_med = rowMeans(MES$vfundo), vfundo_max = apply(MES$vfundo, 1, max), sst_med = rowMeans(MES$sst),
  sst_amp = apply(MES$sst, 1, max) - apply(MES$sst, 1, min), mld_max = apply(MES$mld, 1, max), gelo_meses = rowSums(MES$gelo >= 0.15))
ok <- complete.cases(GB[, 2:11]); GB <- GB[ok, ]; MES <- lapply(MES, function(x) x[ok, , drop = FALSE]); GB$faixa <- faixa(GB$prof)
cat(format(Sys.time()), "grade:", nrow(GB), "celulas\n")

# ---- A within-year variation and threshold crossing ----
amp <- sapply(MES, function(x) apply(x, 1, max) - apply(x, 1, min))
variacao <- bind_rows(lapply(names(MES), function(v) data.frame(variavel = v, faixa = NOME_FX[GB$faixa], a = amp[, v]) |> group_by(variavel, faixa) |>
  summarise(amp_mediana = median(a), amp_p90 = quantile(a, 0.9), .groups = "drop")))
rob <- read.csv(file.path(raiz, "saida/etapa3_n1_limiares_robustos.csv")) |> filter(grepl("bentos", lista))
mapa_v <- c(tfundo_med = "tfundo", tfundo_min = "tfundo", tfundo_max = "tfundo", vfundo_med = "vfundo", vfundo_max = "vfundo", sst_med = "sst", mld_max = "mld")
rob <- rob |> filter(preditor %in% names(mapa_v))
limiares <- bind_rows(lapply(seq_len(nrow(rob)), function(k) { x <- MES[[mapa_v[[rob$preditor[k]]]]]; L <- rob$limiar[k]
  cruza <- apply(x, 1, min) < L & apply(x, 1, max) > L
  data.frame(preditor = rob$preditor[k], variavel_mensal = mapa_v[[rob$preditor[k]]], limiar = L, grupos = rob$grupos[k], faixa = NOME_FX[GB$faixa], cruza = cruza) |>
    group_by(preditor, variavel_mensal, limiar, grupos, faixa) |> summarise(celulas = n(), frac_cruza = mean(cruza), .groups = "drop") }))

# ---- B HB stability month by month ----
mod <- readRDS(file.path(dir_t, "habitats_novos_modelo.rds")); CB <- mod$CB; FB <- mod$FB; rm(mod); gc(verbose = FALSE)
hab_ano <- FB$habitat[atribui_b(transf(GB, "bentos"), GB$prof, CB)]
HM <- sapply(1:12, function(m) { E <- GB; E$tfundo_med <- MES$tfundo[, m]; E$sst_med <- MES$sst[, m]; E$vfundo_med <- MES$vfundo[, m]
  FB$habitat[atribui_b(transf(E, "bentos"), GB$prof, CB)] })
fora <- rowSums(HM != hab_ano)
estab <- data.frame(faixa = NOME_FX[GB$faixa], fora = fora) |> group_by(faixa) |>
  summarise(celulas = n(), muda_algum_mes = mean(fora > 0), muda_3_ou_mais = mean(fora >= 3), meses_fora_media = mean(fora), .groups = "drop") |>
  bind_rows(data.frame(faixa = "todas", celulas = length(fora), muda_algum_mes = mean(fora > 0), muda_3_ou_mais = mean(fora >= 3), meses_fora_media = mean(fora)))
r <- rast(g25); v <- rep(NA, ncell(r)); v[GB$cel] <- fora; r[] <- v
png(file.path(raiz, "saida/mapa_bentos_meses_fora.png"), 2000, 1080, res = 150, bg = "white"); layout(matrix(1:2), heights = c(1000, 80)); par(mar = c(0.5, 1, 2, 1))
cores <- c("#e4e4e4", hcl.colors(12, "YlOrRd", rev = TRUE))
plot(r, col = cores, breaks = c(-0.5, seq(0.5, 12.5, 1)), legend = FALSE, axes = FALSE, box = FALSE, maxcell = ncell(r),
     main = "Bentos: meses em que o habitat do mês difere do habitat anual")
par(mar = c(0, 1, 0, 1)); plot.new(); legend("top", legend = 0:12, fill = cores, horiz = TRUE, bty = "n", cex = 0.9, title = "meses fora do habitat anual")
dev.off()
cat(format(Sys.time()), "estabilidade feita\n")

# ---- C fauna: quarter within cell ----
b <- open_dataset(file.path(raiz, "saida/tabela_certa_bentos_n1.parquet")) |> filter(year >= 1993, nivel_certeza %in% c("A", "B"), !is.na(month)) |>
  select(grupo, cell_id, event_id, taxon_unidade, month) |> collect() |> filter(!is.na(taxon_unidade), taxon_unidade != "") |>
  mutate(tri = paste0("T", (month - 1) %/% 3 + 1))
teste <- function(g) {
  suppressPackageStartupMessages({library(dplyr); library(tidyr); library(vegan)}); source(file.path(raiz, "scripts/funcoes_permanova.R"))
  d <- b[b$grupo == g, ] |> mutate(unidade = paste(cell_id, tri))
  un <- d |> group_by(unidade, cell_id, tri) |> summarise(ev = n_distinct(event_id), .groups = "drop") |> filter(ev >= 10) |>
    group_by(cell_id) |> filter(n() >= 2) |> ungroup() |> arrange(unidade)
  if (length(unique(un$cell_id)) < 5) return(data.frame(grupo = g, celulas = length(unique(un$cell_id)), unidades = nrow(un)))
  set.seed(1001); sor <- d |> semi_join(un, by = "unidade") |> distinct(unidade, event_id) |> group_by(unidade) |> slice_sample(n = 10) |> ungroup()
  X <- d |> semi_join(sor, by = c("unidade", "event_id")) |> distinct(unidade, taxon_unidade) |> mutate(v = 1L) |>
    pivot_wider(names_from = taxon_unidade, values_from = v, values_fill = 0L) |> arrange(unidade)
  Xm <- as.matrix(X[, -1]); Dm <- as.matrix(sqrt(betadiver(Xm, "sim"))); un <- un[match(X$unidade, un$unidade), ]
  set.seed(7); x <- permanova_rapida(as.dist(Dm), un$tri, un$cell_id, 999)
  set.seed(8); eq <- tryCatch(permanova_esforco_igual(Dm, un$tri, un$cell_id, n_classe = 10, B = 200), error = function(e) NULL)
  data.frame(grupo = g, celulas = length(unique(un$cell_id)), unidades = nrow(un), R2 = x$R2, omega2 = x$omega2, p = x$p,
    omega2_eq = if (!is.null(eq)) eq$omega2_eq else NA, poder_eq = if (!is.null(eq)) eq$poder_eq else NA)
}
cl <- makeCluster(7); clusterExport(cl, c("raiz", "b"))
fauna <- bind_rows(parLapply(cl, sort(unique(b$grupo)), teste)); stopCluster(cl)
hv <- read.csv(file.path(raiz, "saida/etapa3_n1_habitats_novos_validacao.csv")) |> filter(comp == "bentos", classificacao == "habitats novos") |> select(grupo, omega2_espacial_HB = omega2)
fauna <- fauna |> left_join(hv, by = "grupo")
write.csv(variacao, file.path(raiz, "saida/etapa3_n1_sazonal_variacao.csv"), row.names = FALSE)
write.csv(limiares, file.path(raiz, "saida/etapa3_n1_sazonal_limiares.csv"), row.names = FALSE)
write.csv(estab, file.path(raiz, "saida/etapa3_n1_sazonal_estabilidade.csv"), row.names = FALSE)
write.csv(fauna, file.path(raiz, "saida/etapa3_n1_sazonal_fauna.csv"), row.names = FALSE)
options(width = 200); print(as.data.frame(variacao)); print(as.data.frame(limiares |> mutate(across(where(is.double), ~ signif(.x, 3)))))
print(as.data.frame(estab)); print(fauna |> mutate(across(where(is.double), ~ round(.x, 3)))); cat(format(Sys.time()), "FIM\n")
