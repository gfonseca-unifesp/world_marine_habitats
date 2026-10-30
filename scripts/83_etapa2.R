# 83_etapa2.R  (Phase 8, stage 2: diagnosis; argument nivel=1 or nivel=2)
# For each contrast tested in stage 1 (family; level 1 without the no-MLPA sensitivity; level 2 with 100 km blocks):
#  2c where the classification fails: each unit is classified by the nearest fauna centroid (sqrt(beta_sim)),
#     computed without the unit's block (2,500 km at level 1; 100 km at level 2). Accuracy, balanced accuracy, kappa,
#     sensitivity and precision per habitat, with which habitats each one is confused with, and the map of units.
#  2d discontinuities within each habitat: fauna subgroups (PAM, silhouette) that are spatially coherent.
#  Dispersion of each habitat: mean distance to its own centroid.
# Same unit, filters and draw 1 of stage 1. Parallel in 12 processes.
# Outputs: intermediario/etapa2_n<k>/<id>.rds; saida/etapa2_n<k>_{classif,classif_classe,confusao,descontinuidades}.csv;
#         saida/etapa2_n<k>_unidades.parquet
suppressPackageStartupMessages({library(dplyr); library(parallel); library(arrow)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
NIVEL <- if ("nivel=2" %in% commandArgs(TRUE)) 2L else 1L; NUCLEOS <- 12
dir_e1 <- file.path(raiz, if (NIVEL == 1) "intermediario/etapa1_n1" else "intermediario/etapa1_n2_b100")
dir_t <- file.path(raiz, sprintf("intermediario/etapa2_n%d", NIVEL)); dir.create(dir_t, showWarnings = FALSE)
arqs <- list.files(dir_e1, pattern = "\\.rds$", full.names = TRUE)
arqs <- arqs[vapply(arqs, function(a) { r <- readRDS(a); !is.null(r$perm) && !isTRUE(r$tarefa$sem_mlpa) }, TRUE)]
cat("nivel", NIVEL, "| contrastes:", length(arqs), "\n")

roda <- function(a) {
  saida <- file.path(dir_t, basename(a)); if (file.exists(saida)) return(basename(a))
  source(file.path(raiz, "scripts/funcoes_grade.R")); source(file.path(raiz, "scripts/funcoes_etapa2.R"))
  tk <- readRDS(a)$tarefa
  pr <- if (NIVEL == 1) prepara_n1(tk, raiz) else prepara_n2(tk, raiz, 100)
  if (is.null(pr)) { saveRDS(list(tarefa = tk), saida); return(basename(a)) }
  meta <- pr$meta; S <- as.matrix(betadiver(pr$X, "sim"))      # S = D^2, D = sqrt(beta_sim)
  cl <- classifica_lbo(S, meta$classe, meta$bloco, un_min = UN_MIN)
  mt <- metricas_classif(cl)
  # dispersion: distance to own centroid (all units of the class)
  grp <- factor(meta$classe); Z <- model.matrix(~ grp - 1); nc <- colSums(Z); Q <- S %*% Z; W <- colSums(Z * Q)
  d2 <- sweep(Q, 2, nc, "/") - matrix(W / (2 * nc^2), nrow(S), ncol(Z), byrow = TRUE)
  meta$dist_centroide <- sqrt(pmax(d2[cbind(seq_len(nrow(S)), as.integer(grp))], 0))
  meta$prevista <- cl$prevista; meta$correta <- cl$correta; meta$subgrupo <- NA_integer_
  ds <- list()
  set.seed(8000)
  for (c in levels(grp)) {
    i <- which(meta$classe == c)
    r <- descontinuidade(sqrt(S[i, i]), cbind(meta$x[i], meta$y[i]))
    if (!is.null(r)) { meta$subgrupo[i] <- r$subgrupo; ds[[c]] <- cbind(classe = c, r$resumo) }
  }
  saveRDS(list(tarefa = tk, unidades = meta, metricas = mt, descont = bind_rows(ds)), saida); basename(a)
}

cl <- makeCluster(NUCLEOS); clusterExport(cl, c("roda", "raiz", "dir_t", "NIVEL"))
t0 <- Sys.time()
feitos <- unlist(parLapplyLB(cl, arqs, function(a) tryCatch(roda(a), error = function(e) paste("ERRO", basename(a), conditionMessage(e)))))
stopCluster(cl); print(feitos[grepl("^ERRO", feitos)]); cat("tempo:", format(Sys.time() - t0), "\n")

# ---- consolidation ----
chave <- if (NIVEL == 1) c("comp", "grupo", "nivel", "fatia") else c("regiao", "nivel", "grupo_ecologico", "fatia", "res_escolhida")
rs <- lapply(file.path(dir_t, basename(arqs)), function(a) if (file.exists(a)) readRDS(a))
rs <- rs[vapply(rs, function(r) !is.null(r$metricas), TRUE)]
b <- function(r, n) r$tarefa[rep(1, n), chave, drop = FALSE]
geral <- bind_rows(lapply(rs, function(r) cbind(b(r, 1), r$metricas$geral)))
pc <- bind_rows(lapply(rs, function(r) { x <- r$metricas$por_classe
  disp <- r$unidades |> group_by(classe) |> summarise(dispersao = mean(dist_centroide), .groups = "drop")
  x <- x |> left_join(disp, by = "classe")
  if (!is.null(r$descont) && nrow(r$descont)) x <- x |> left_join(r$descont |> rename(unidades_desc = unidades), by = "classe")
  cbind(b(r, nrow(x)), x) }))
cf <- bind_rows(lapply(rs, function(r) cbind(b(r, nrow(r$metricas$confusao)), r$metricas$confusao)))
ds <- bind_rows(lapply(rs, function(r) if (!is.null(r$descont) && nrow(r$descont)) cbind(b(r, nrow(r$descont)), r$descont)))
un <- bind_rows(lapply(rs, function(r) cbind(b(r, nrow(r$unidades)), r$unidades |> select(unidade, cell, classe, bloco, lon, lat, prevista, correta, dist_centroide, subgrupo))))
pre <- sprintf("saida/etapa2_n%d_", NIVEL)
write.csv(geral, file.path(raiz, paste0(pre, "classif.csv")), row.names = FALSE)
write.csv(pc, file.path(raiz, paste0(pre, "classif_classe.csv")), row.names = FALSE)
write.csv(cf, file.path(raiz, paste0(pre, "confusao.csv")), row.names = FALSE)
write.csv(ds, file.path(raiz, paste0(pre, "descontinuidades.csv")), row.names = FALSE)
write_parquet(un, file.path(raiz, paste0(pre, "unidades.parquet")))
options(width = 220); print(geral |> mutate(across(where(is.numeric), ~ round(.x, 3))) |> head(60))
cat("\ndescontinuidades:", sum(ds$descontinuidade, na.rm = TRUE), "de", nrow(ds), "habitats avaliados\n")
