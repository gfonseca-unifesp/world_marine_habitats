# 88_gf_gdm.R  (Phase 8, stages 3.2 and 3.3: Gradient Forest and GDM at level 1; M1 approved at Gate 3)
# Input: intermediario/etapa3_n1/<comp>_<grupo>.rds (87_preditores.R). One group per process (parallel).
# Gate 7 (30/10/2026): response = frequency of occurrence (F, 0 to 1); ocean basin (bac_1..8, 0/1) among the predictors,
#   to absorb the biogeographic effect: it enters the importance, but generates no curves or thresholds. The GF object is
#   saved (gfobj_*) for 3.5.
# Gate 8 (30/10/2026): bottom salinity is dropped from benthos (the Gate 7 sensitivity showed an R2 drop <= 0.01). The
#   results with it were kept in gf_comsfundo_bentos_*.rds and saida/etapa3_n1_gf_importancia_comsfundo.csv.
# Gradient Forest (Ellis et al. 2012): taxa present in >= 10 units and in <= 95% of them; 500 trees; conditional importance
#   (corr.threshold 0.5). Per predictor: importance (weighted R2), cumulative importance curve and thresholds =
#   peaks of the curve's rate of change with density >= 2x the gradient mean (at most 3 per predictor).
# GDM (Ferrier et al. 2007): Bray-Curtis on the frequencies, predictors + geographic distance; 3 subsamples of up to 1,000
#   units; deviance explained; importance = I-spline height; threshold = steepest stretch of the I-spline.
# Outputs: intermediario/etapa3_n1/gf_<id>.rds; saida/etapa3_n1_gf_importancia.csv, _gf_limiares.csv, _gf_curvas.csv,
#         _gdm.csv, _gdm_curvas.csv
suppressPackageStartupMessages({library(parallel); library(dplyr)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); dir_t <- file.path(raiz, "intermediario/etapa3_n1"); SUF <- ""
arqs <- list.files(dir_t, pattern = "^(bentos|pelagico)_.*\\.rds$", full.names = TRUE)

roda <- function(a) {
  suppressPackageStartupMessages({library(gradientForest); library(gdm); library(dplyr)}); source(file.path(raiz, "scripts/funcoes_gf.R"))
  saida <- file.path(dirname(a), paste0("gf", SUF, "_", basename(a))); if (file.exists(saida)) return(saida)
  r <- readRDS(a); X <- r$F; E <- r$E; n <- nrow(X)
  if (r$comp == "bentos") E$sfundo_med <- NULL                                    # Gate 8, item 2
  E <- E[, !(grepl("^bac_", names(E)) & colSums(E) == 0), drop = FALSE]          # basin without units is dropped
  prev <- colSums(X > 0); X <- X[, prev >= 10 & prev <= 0.95 * n, drop = FALSE]
  out <- list(comp = r$comp, grupo = r$grupo, unidades = n, taxons = ncol(X))
  # ---- Gradient Forest ----
  set.seed(1); dat <- cbind(E, as.data.frame(X)); names(dat) <- make.names(names(dat), unique = TRUE); resp <- names(dat)[-(1:ncol(E))]
  gf <- gradientForest(dat, predictor.vars = names(E), response.vars = resp, ntree = 500, maxLevel = floor(log2(0.368 * n / 2)),
                       corr.threshold = 0.5, trace = FALSE)
  saveRDS(gf, file.path(dirname(a), paste0("gfobj", SUF, "_", basename(a))))
  imp <- importance(gf, type = "Weighted"); out$gf_imp <- data.frame(preditor = names(imp), importancia = as.numeric(imp))
  out$gf_r2 <- data.frame(especies_ajustadas = length(gf$result), r2_medio = mean(gf$result), r2_total = sum(gf$result))
  curvas <- list(); lim <- list()
  for (v in grep("^bac_", names(E), value = TRUE, invert = TRUE)) {
    ci <- tryCatch(cumimp(gf, v, type = "Overall"), error = function(e) NULL); if (is.null(ci) || length(ci$x) < 5) next
    curvas[[v]] <- data.frame(preditor = v, x = ci$x, y = ci$y)
    l <- limiares_curva(ci$x, ci$y); if (!is.null(l)) lim[[v]] <- cbind(preditor = v, l)
  }
  out$gf_curvas <- bind_rows(curvas); out$gf_lim <- bind_rows(lim)
  # ---- GDM ----
  gd <- list(); gc_ <- list()
  for (s in 1:3) {
    set.seed(100 + s); i <- if (n > 1000) sort(sample(n, 1000)) else seq_len(n)
    bio <- data.frame(site = seq_along(i), x = r$meta$x[i], y = r$meta$y[i], X[i, , drop = FALSE]); bio <- bio[, c(TRUE, TRUE, TRUE, colSums(X[i, , drop = FALSE]) > 0)]
    pr <- data.frame(site = seq_along(i), x = r$meta$x[i], y = r$meta$y[i], E[i, , drop = FALSE])
    sp <- formatsitepair(bio, bioFormat = 1, dist = "bray", abundance = TRUE, siteColumn = "site", XColumn = "x", YColumn = "y", predData = pr)
    sp <- sp[is.finite(sp$distance), ]
    m <- tryCatch(gdm(sp, geo = TRUE), error = function(e) NULL); if (is.null(m)) next
    isp <- isplineExtract(m)
    alt <- apply(isp$y, 2, max)
    gd[[s]] <- data.frame(subamostra = s, preditor = names(alt), importancia = alt, deviance_explicada = m$explained)
    gc_[[s]] <- bind_rows(lapply(colnames(isp$x), function(v) { xs <- isp$x[, v]; ys <- isp$y[, v]; k <- which.max(diff(ys))
      data.frame(subamostra = s, preditor = v, x = xs, y = ys, limiar = (xs[k] + xs[k + 1]) / 2) }))
  }
  out$gdm <- bind_rows(gd); out$gdm_curvas <- bind_rows(gc_)
  saveRDS(out, saida); saida
}

cl <- makeCluster(min(12, length(arqs))); clusterExport(cl, c("roda", "raiz", "SUF"))
t0 <- Sys.time(); feitos <- unlist(parLapplyLB(cl, arqs, function(a) tryCatch(roda(a), error = function(e) paste("ERRO", basename(a), conditionMessage(e)))))
stopCluster(cl); print(feitos[grepl("^ERRO", feitos)]); cat("tempo:", format(Sys.time() - t0), "\n")

rs <- lapply(list.files(dir_t, pattern = paste0("^gf", SUF, "_(bentos|pelagico)"), full.names = TRUE), readRDS)
cab <- function(r, x) if (!is.null(x) && nrow(x)) cbind(comp = r$comp, grupo = r$grupo, x)
imp <- bind_rows(lapply(rs, function(r) cab(r, merge(r$gf_imp, data.frame(r$gf_r2, unidades = r$unidades, taxons = r$taxons)))))
write.csv(imp, file.path(raiz, paste0("saida/etapa3_n1_gf_importancia", SUF, ".csv")), row.names = FALSE)
write.csv(bind_rows(lapply(rs, function(r) cab(r, r$gf_lim))), file.path(raiz, paste0("saida/etapa3_n1_gf_limiares", SUF, ".csv")), row.names = FALSE)
write.csv(bind_rows(lapply(rs, function(r) cab(r, r$gf_curvas))), file.path(raiz, paste0("saida/etapa3_n1_gf_curvas", SUF, ".csv")), row.names = FALSE)
gdm_t <- bind_rows(lapply(rs, function(r) cab(r, r$gdm)))
write.csv(gdm_t, file.path(raiz, paste0("saida/etapa3_n1_gdm", SUF, ".csv")), row.names = FALSE)
write.csv(bind_rows(lapply(rs, function(r) cab(r, r$gdm_curvas))), file.path(raiz, paste0("saida/etapa3_n1_gdm_curvas", SUF, ".csv")), row.names = FALSE)
options(width = 200)
print(imp |> group_by(comp, grupo) |> slice_max(importancia, n = 3) |> mutate(importancia = round(importancia, 4), r2_medio = round(r2_medio, 3)) |> as.data.frame())
print(gdm_t |> group_by(comp, grupo) |> summarise(deviance = round(mean(deviance_explicada), 1), .groups = "drop") |> as.data.frame())
