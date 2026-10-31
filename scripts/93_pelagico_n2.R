# 93_pelagico_n2.R  (Phase 8; Gate 8 item 3 and option C of 30/10/2026: the pelagic at level 2)
# Data: saida/tabela_certa_pelagico_n2.parquet (8 regional windows), records >= 1993 with month (certainty A and B) and with a
# sampling depth (34% of the records have none, so no z). Ecological groups (grupo_ecologico) replace the level-1 taxonomic groups.
#  1 units = 25 km cell x quarter x depth layer (0-50, 50-200, 200-500, 500-1000, > 1000 m), >= 10 events; X = 10 events drawn
#    (seed 1001), F = frequency of occurrence; predictors as in 87_preditores.R (pelagic) plus the level-2 region as 0/1
#    covariates (reg_*; like the basin at level 1: in the importance, not in the curves); spatial blocks of 100 km;
#  2 Gradient Forest (as 88_gf_gdm.R, 500 trees) for groups with >= MIN_GF units; thresholds compared with the level-1 pelagic
#    robust thresholds (match = within 10% of the predictor range);
#  3 the level-1 initial pelagic clusters (k0 of 89b) tested with the level-2 fauna: each unit gets the k0 of the nearest
#    level-1 centroid (level-1 GF transformation, kmP); pairs per ecological group (pares_classes, omega2 >= 0.01 and power >= 0.8,
#    blocks of 100 km); counts added to the level-1 consensus and the 89b merge rule rerun, to see how many stand-by pelagic
#    habitats (HP) become tested.
# Outputs: saida/pelagico_n2_{unidades,gf_importancia,gf_limiares,gf_curvas,limiares_vs_n1,pares,consenso,fusoes}.csv
suppressPackageStartupMessages({library(arrow); library(dplyr); library(tidyr); library(terra); library(parallel)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); source(file.path(raiz, "scripts/funcoes_grade.R")); source(file.path(raiz, "scripts/funcoes_gf.R"))
dir_t <- file.path(raiz, "intermediario/pelagico_n2"); dir.create(dir_t, showWarnings = FALSE)
EV_MIN <- 10; MIN_GF <- 150; K0 <- 30; Z_PADRAO <- c(0, 10, 25, 50, 100, 200, 300, 500, 1000, 1500, 2000, 3000, 4000, 5000)
CAMADAS <- c(-Inf, 50, 200, 500, 1000, Inf); TRI <- list(T1 = 1:3, T2 = 4:6, T3 = 7:9, T4 = 10:12)
g25 <- ease2_grade("25km"); M <- lapply(1:12, function(m) rast(file.path(raiz, sprintf("camadas/glorys/mes/glorys_%02d_25km.tif", m))))
cel_de <- function(id) { p <- do.call(rbind, strsplit(id, "_")); cellFromRowCol(g25, as.integer(p[, 2]), as.integer(p[, 3])) }
bloco_de <- function(id) { p <- do.call(rbind, strsplit(id, "_")); paste(as.integer(p[, 2]) %/% 4, as.integer(p[, 3]) %/% 4, sep = "_") }   # 100 km
interp <- function(perfil, z) vapply(seq_len(nrow(perfil)), function(i) { v <- perfil[i, ]; k <- !is.na(v)
  if (sum(k) == 0) NA_real_ else if (sum(k) == 1) v[k] else approx(Z_PADRAO[k], v[k], xout = min(z[i], max(Z_PADRAO[k])), rule = 2)$y }, 0)

# ---------------- 1 units and predictors ----------------
p <- open_dataset(file.path(raiz, "saida/tabela_certa_pelagico_n2.parquet")) |> filter(year >= 1993, nivel_certeza %in% c("A", "B"), !is.na(mes), !is.na(depth_min) | !is.na(depth_max)) |>
  select(regiao, grupo_ecologico, cell_25, evento, taxon_unidade, mes, depth_min, depth_max) |> collect() |>
  filter(!is.na(taxon_unidade), taxon_unidade != "", !is.na(grupo_ecologico), !is.na(cell_25)) |>
  mutate(tri = paste0("T", (mes - 1) %/% 3 + 1), zreg = rowMeans(cbind(depth_min, depth_max), na.rm = TRUE),
         camada = as.character(cut(zreg, CAMADAS, labels = c("0-50", "50-200", "200-500", "500-1000", ">1000"))))
REG <- sort(unique(p$regiao))
monta <- function(g) {
  d <- p |> filter(grupo_ecologico == g) |> mutate(unidade = paste(cell_25, tri, camada))
  un <- d |> group_by(unidade) |> summarise(ev = n_distinct(evento), cell_id = first(cell_25), tri = first(tri), camada = first(camada), regiao = first(regiao),
    z = median(zreg), .groups = "drop") |> filter(ev >= EV_MIN) |> arrange(unidade)
  if (nrow(un) < 10) return(NULL)
  d <- d |> filter(unidade %in% un$unidade)
  set.seed(1001); sor <- d |> distinct(unidade, evento) |> group_by(unidade) |> slice_sample(n = EV_MIN) |> ungroup()
  X <- d |> semi_join(sor, by = c("unidade", "evento")) |> distinct(unidade, taxon_unidade) |> mutate(v = 1L) |>
    pivot_wider(names_from = taxon_unidade, values_from = v, values_fill = 0L) |> arrange(unidade)
  X <- as.matrix(X[, -1]) |> `rownames<-`(X$unidade)
  F <- d |> distinct(unidade, evento, taxon_unidade) |> count(unidade, taxon_unidade) |> left_join(un |> select(unidade, ev), by = "unidade") |>
    mutate(f = n / ev) |> select(-n, -ev) |> pivot_wider(names_from = taxon_unidade, values_from = f, values_fill = 0) |> arrange(unidade)
  F <- as.matrix(F[, -1]) |> `rownames<-`(F$unidade)
  cel <- cel_de(un$cell_id)
  E <- data.frame(z = un$z, t_z = NA_real_, s_z = NA_real_, sst = NA_real_, sss = NA_real_, mld = NA_real_, prof = values(M[[1]][["prof_glorys"]])[cel, 1], vsup = NA_real_, gelo = NA_real_)
  for (tr in names(TRI)) { i <- which(un$tri == tr); if (!length(i)) next; ms <- TRI[[tr]]
    med <- function(camada) Reduce(`+`, lapply(ms, function(m) values(M[[m]][[camada]])[cel[i], 1])) / length(ms)
    tz <- sapply(Z_PADRAO, function(zz) med(paste0("thetao_", zz))); sz <- sapply(Z_PADRAO, function(zz) med(paste0("so_", zz)))
    if (length(i) == 1) { tz <- matrix(tz, 1); sz <- matrix(sz, 1) }
    E$t_z[i] <- interp(tz, un$z[i]); E$s_z[i] <- interp(sz, un$z[i]); E$sst[i] <- med("thetao_0"); E$sss[i] <- med("so_0")
    E$mld[i] <- med("mld"); E$vsup[i] <- med("vsup"); E$gelo[i] <- med("gelo") }
  E$gelo[is.na(E$gelo)] <- 0; E$z_mld <- log2((E$z + 1) / (E$mld + 1))
  R <- sapply(seq_along(REG), function(k) as.integer(un$regiao == REG[k])); colnames(R) <- paste0("reg_", seq_along(REG)); E <- cbind(E, as.data.frame(R))
  ok <- complete.cases(E)
  list(grupo = g, meta = cbind(un, bloco = bloco_de(un$cell_id))[ok, ], X = X[ok, colSums(X[ok, , drop = FALSE]) > 0, drop = FALSE],
       F = F[ok, colSums(F[ok, , drop = FALSE]) > 0, drop = FALSE], E = E[ok, ])
}
U <- Filter(Negate(is.null), lapply(sort(unique(p$grupo_ecologico)), monta)); names(U) <- sapply(U, `[[`, "grupo")
resumo <- bind_rows(lapply(U, function(u) data.frame(grupo = u$grupo, unidades = nrow(u$meta), taxons = ncol(u$X), regioes = n_distinct(u$meta$regiao), blocos = n_distinct(u$meta$bloco))))
write.csv(resumo, file.path(raiz, "saida/pelagico_n2_unidades.csv"), row.names = FALSE); print(resumo)
for (u in U) saveRDS(u, file.path(dir_t, paste0(gsub(" ", "_", u$grupo), ".rds")))

# ---------------- 2 Gradient Forest ----------------
roda_gf <- function(u) {
  suppressPackageStartupMessages({library(gradientForest); library(dplyr)}); source(file.path(raiz, "scripts/funcoes_gf.R"))
  X <- u$F; E <- u$E; n <- nrow(X); E <- E[, !(grepl("^reg_", names(E)) & colSums(E) == 0), drop = FALSE]; E <- E[, vapply(E, function(x) sd(x) > 0, TRUE), drop = FALSE]
  prev <- colSums(X > 0); X <- X[, prev >= 10 & prev <= 0.95 * n, drop = FALSE]
  set.seed(1); dat <- cbind(E, as.data.frame(X)); names(dat) <- make.names(names(dat), unique = TRUE)
  gf <- gradientForest(dat, predictor.vars = names(E), response.vars = names(dat)[-(1:ncol(E))], ntree = 500, maxLevel = floor(log2(0.368 * n / 2)), corr.threshold = 0.5, trace = FALSE)
  imp <- importance(gf, type = "Weighted"); vs <- grep("^reg_", names(E), value = TRUE, invert = TRUE)
  cur <- lapply(vs, function(v) tryCatch(cumimp(gf, v, type = "Overall"), error = function(e) NULL)); names(cur) <- vs; cur <- cur[!vapply(cur, function(c) is.null(c) || length(c$x) < 5, TRUE)]
  list(grupo = u$grupo, imp = data.frame(grupo = u$grupo, preditor = names(imp), importancia = as.numeric(imp), especies_ajustadas = length(gf$result),
      r2_medio = mean(gf$result), unidades = n, taxons = ncol(X)),
    curvas = bind_rows(lapply(names(cur), function(v) data.frame(grupo = u$grupo, preditor = v, x = cur[[v]]$x, y = cur[[v]]$y))),
    lim = bind_rows(lapply(names(cur), function(v) { l <- limiares_curva(cur[[v]]$x, cur[[v]]$y); if (!is.null(l)) cbind(grupo = u$grupo, preditor = v, l) })))
}
UG <- Filter(function(u) nrow(u$meta) >= MIN_GF, U)
cl <- makeCluster(length(UG)); clusterExport(cl, "raiz"); G <- parLapply(cl, UG, roda_gf); stopCluster(cl)
gi <- bind_rows(lapply(G, `[[`, "imp")); gl <- bind_rows(lapply(G, `[[`, "lim"))
write.csv(gi, file.path(raiz, "saida/pelagico_n2_gf_importancia.csv"), row.names = FALSE)
write.csv(gl, file.path(raiz, "saida/pelagico_n2_gf_limiares.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(G, `[[`, "curvas")), file.path(raiz, "saida/pelagico_n2_gf_curvas.csv"), row.names = FALSE)
rob <- read.csv(file.path(raiz, "saida/etapa3_n1_limiares_robustos.csv")) |> filter(grepl("pelagico", lista))
amp <- bind_rows(lapply(UG, function(u) data.frame(grupo = u$grupo, preditor = names(u$E), amp = sapply(u$E, function(x) diff(range(x))))))
vs_n1 <- gl |> left_join(amp, by = c("grupo", "preditor")) |> rowwise() |>
  mutate(limiar_n1 = { r <- rob$limiar[rob$preditor == preditor]; if (length(r)) r[which.min(abs(r - limiar))] else NA_real_ },
         coincide_n1 = !is.na(limiar_n1) && abs(limiar_n1 - limiar) <= 0.10 * amp) |> ungroup()
write.csv(vs_n1, file.path(raiz, "saida/pelagico_n2_limiares_vs_n1.csv"), row.names = FALSE)
cat(format(Sys.time()), "GF feito\n")

# ---------------- 3 level-1 pelagic clusters tested with the level-2 fauna ----------------
curvas_p <- open_dataset(file.path(raiz, "saida/etapa3_n1_gf_curvas.csv"), format = "csv") |> filter(comp == "pelagico") |> collect() |> as.data.frame()
transf <- function(E) { vs <- intersect(unique(curvas_p$preditor), names(E))   # as in 89b (level-1 pelagic curves)
  sapply(vs, function(v) { s <- 0; for (g in unique(curvas_p$grupo[curvas_p$preditor == v])) { c1 <- curvas_p[curvas_p$grupo == g & curvas_p$preditor == v, ]
    s <- s + approx(c1$x, c1$y, E[[v]], rule = 2, ties = mean)$y }; s }) }
mais_proximo <- function(T, C) { d <- sapply(seq_len(nrow(C)), function(k) colSums((t(T) - C[k, ])^2)); max.col(-d, ties.method = "first") }
mod <- readRDS(file.path(raiz, "intermediario/etapa3_n1/habitats_novos_modelo.rds")); CP <- mod$kmP$centers; FP <- mod$FP; rm(mod); gc(verbose = FALSE)
for (k in names(U)) U[[k]]$k0 <- mais_proximo(transf(U[[k]]$E), CP)
testa <- function(u) {
  suppressPackageStartupMessages({library(vegan); library(dplyr)}); source(file.path(raiz, "scripts/funcoes_permanova.R"))
  Dm <- as.matrix(sqrt(betadiver(u$X, "sim"))); set.seed(7100)
  x <- pares_classes(Dm, as.character(u$k0), u$meta$bloco, un_min = 5, n_classe = 10, B = 100)
  if (!is.null(x)) cbind(grupo = u$grupo, x)
}
cl <- makeCluster(length(U)); clusterExport(cl, "raiz"); pares <- bind_rows(parLapply(cl, U, testa)); stopCluster(cl)
pares <- pares |> filter(!is.na(separa)); write.csv(pares, file.path(raiz, "saida/pelagico_n2_pares.csv"), row.names = FALSE)
n1 <- read.csv(file.path(raiz, "saida/etapa3_n1_habitats_novos_pares.csv")) |> filter(comp == "pelagico") |>
  transmute(a = pmin(classe_a, classe_b), b = pmax(classe_a, classe_b), grupos, distintos) |> rename(classe_a = a, classe_b = b)   # pairs were sorted as text in 89b
n2 <- pares |> mutate(a = pmin(as.integer(classe_a), as.integer(classe_b)), b = pmax(as.integer(classe_a), as.integer(classe_b))) |>
  group_by(classe_a = a, classe_b = b) |> summarise(grupos = n(), distintos = sum(separa), .groups = "drop")
cons <- full_join(n1, n2, by = c("classe_a", "classe_b"), suffix = c("_n1", "_n2")) |> mutate(across(starts_with(c("grupos", "distintos")), ~ coalesce(.x, 0L))) |>
  mutate(grupos = grupos_n1 + grupos_n2, distintos = distintos_n1 + distintos_n2, frac = distintos / grupos,
         consenso = ifelse(grupos >= 2 & frac <= 1/3, "equivalente", ifelse(frac >= 2/3, "distinto", "depende do grupo")))
write.csv(cons, file.path(raiz, "saida/pelagico_n2_consenso.csv"), row.names = FALSE)
funde <- function(e) { testados <- sort(unique(c(e$classe_a, e$classe_b)))   # 89b merge rule
  Mx <- matrix(1, K0, K0); diag(Mx) <- 0; Mx[cbind(e$classe_a, e$classe_b)] <- ifelse(e$grupos >= 2, e$frac, 1); Mx[cbind(e$classe_b, e$classe_a)] <- Mx[cbind(e$classe_a, e$classe_b)]
  g <- cutree(hclust(as.dist(Mx), method = "complete"), h = 1/3); fin <- match(g, unique(g)); te <- (1:K0) %in% testados
  data.frame(k0 = 1:K0, final = fin, testado = te, standby = !tapply(te, fin, any)[as.character(fin)]) }
F1 <- funde(cons |> filter(grupos_n1 > 0) |> transmute(classe_a, classe_b, grupos = grupos_n1, frac = distintos_n1 / grupos_n1))
F2 <- funde(cons)
fus <- data.frame(k0 = 1:K0, hp_n1 = FP$habitat, final_n1 = F1$final, standby_n1 = F1$standby, final_n1n2 = F2$final, standby_n1n2 = F2$standby, testado_n2 = (1:K0) %in% c(n2$classe_a, n2$classe_b))
write.csv(fus, file.path(raiz, "saida/pelagico_n2_fusoes.csv"), row.names = FALSE)
options(width = 200); print(gi |> group_by(grupo) |> slice_max(importancia, n = 4) |> mutate(across(where(is.double), ~ round(.x, 3))) |> as.data.frame())
print(vs_n1 |> select(grupo, preditor, limiar, limiar_n1, coincide_n1) |> mutate(across(where(is.double), ~ signif(.x, 4))) |> as.data.frame())
cat("pares testados no nivel 2:", nrow(pares), "| pares de k0 distintos:", nrow(n2), "\n")
cat("HP (nivel 1):", length(unique(F1$final)), "habitats,", sum(!tapply(F1$testado, F1$final, any)), "stand-by | com o nivel 2:", length(unique(F2$final)), "habitats,",
    sum(!tapply(F2$testado, F2$final, any)), "stand-by\n"); cat(format(Sys.time()), "FIM\n")
