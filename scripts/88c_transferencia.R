# 88c_transferencia.R  (Phase 8, stage 3; Gate 8, item 6: do thresholds hold in all basins?)
# GF thresholds are estimated with the units of all basins together and applied to the whole globe; here we test whether
# they recur in each basin. Groups with signal (excluding Euphausiacea and Hyperiidea); basins with >= 150 units in the group.
# For each group x basin, two Gradient Forests (same configuration as 88, with 200 trees):
#  local = only the basin's units -> local thresholds (funcoes_gf.R, same rule as 88);
#  fora  = all units except those of the basin -> curves applied to the basin's units (transferability).
# Robust threshold (89_ligacao_habitats.R) x basin, for each group in the threshold's list:
#  testable  = the range sampled in the basin covers the threshold with a margin of 5% of the group's range;
#  confirmed = a local threshold of the same variable within <= 10% of the group's range.
#  Class: global (confirmed in >= 2 basins and in >= 2/3 of the testable basins), regional (confirmed in >= 1),
#  not confirmed (testable and never confirmed), not testable.
# Transferability: Spearman correlation (Mantel-style) between the transformed environmental distance and the fauna's sqrt(beta_sim)
#  among the basin's units (up to 800): curves of the out-of-basin model (transferred), the local model and the global model (88).
# Outputs: saida/etapa3_n1_transf_{ajustes,limiares_locais,testes,robustos}.csv
suppressPackageStartupMessages({library(parallel); library(dplyr)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); dir_t <- file.path(raiz, "intermediario/etapa3_n1"); dir_c <- file.path(dir_t, "transf"); dir.create(dir_c, showWarnings = FALSE)
UN_MIN <- 150; SEM_SINAL <- c("Euphausiacea", "Hyperiidea")
bac_nome <- read.csv(file.path(raiz, "camadas/bacias_legenda.csv"))
arqs <- list.files(dir_t, pattern = "^(bentos|pelagico)_.*\\.rds$", full.names = TRUE); arqs <- arqs[!grepl(paste(SEM_SINAL, collapse = "|"), arqs)]
jobs <- bind_rows(lapply(arqs, function(a) { r <- readRDS(a); tb <- table(r$meta$bacia); n <- nrow(r$meta)
  b <- as.integer(names(tb)[tb >= UN_MIN & n - tb >= UN_MIN]); if (!length(b)) return(NULL)
  expand.grid(arq = a, bacia = b, tipo = c("local", "fora"), stringsAsFactors = FALSE) |> mutate(n = ifelse(tipo == "local", as.integer(tb[as.character(bacia)]), n - as.integer(tb[as.character(bacia)]))) }))
jobs <- jobs[order(-jobs$n), ]; cat(format(Sys.time()), nrow(jobs), "ajustes\n")

ajusta <- function(j) {
  suppressPackageStartupMessages({library(gradientForest); library(vegan)}); source(file.path(raiz, "scripts/funcoes_gf.R"))
  saida <- file.path(dir_c, sprintf("%s_b%d_%s.rds", sub("\\.rds$", "", basename(j$arq)), j$bacia, j$tipo)); if (file.exists(saida)) return(saida)
  r <- readRDS(j$arq); E <- r$E; if (r$comp == "bentos") E$sfundo_med <- NULL
  na_b <- r$meta$bacia %in% j$bacia; i <- if (j$tipo == "local") which(na_b) else which(!na_b)
  X <- r$F[i, , drop = FALSE]; Ei <- E[i, , drop = FALSE]; n <- length(i)
  Ei <- Ei[, !(grepl("^bac_", names(Ei)) & (j$tipo == "local" | colSums(Ei) == 0)), drop = FALSE]   # in the local fit the basin is constant
  Ei <- Ei[, vapply(Ei, function(x) sd(x) > 0, TRUE), drop = FALSE]   # gradientForest rejects constant predictors (e.g. sea ice in a basin without ice)
  prev <- colSums(X > 0); X <- X[, prev >= 10 & prev <= 0.95 * n, drop = FALSE]
  out <- list(comp = r$comp, grupo = r$grupo, bacia = j$bacia, tipo = j$tipo, unidades = n, taxons = ncol(X), r2 = NA, ajustados = 0, rho = NA, lim = NULL)
  if (ncol(X) >= 3) {
    set.seed(1); dat <- cbind(Ei, as.data.frame(X)); names(dat) <- make.names(names(dat), unique = TRUE)
    gf <- tryCatch(gradientForest(dat, predictor.vars = names(Ei), response.vars = names(dat)[-(1:ncol(Ei))], ntree = 200,
      maxLevel = floor(log2(0.368 * n / 2)), corr.threshold = 0.5, trace = FALSE), error = function(e) NULL)
    if (!is.null(gf) && length(gf$result)) {
      out$r2 <- mean(gf$result); out$ajustados <- length(gf$result)
      vs <- grep("^bac_", names(Ei), value = TRUE, invert = TRUE)
      cur <- lapply(vs, function(v) tryCatch(cumimp(gf, v, type = "Overall"), error = function(e) NULL)); names(cur) <- vs
      cur <- cur[!vapply(cur, function(c) is.null(c) || length(c$x) < 5, TRUE)]
      out$lim <- do.call(rbind, lapply(names(cur), function(v) { l <- limiares_curva(cur[[v]]$x, cur[[v]]$y); if (!is.null(l)) cbind(preditor = v, l) }))
      ib <- which(na_b); if (length(ib) > 800) { set.seed(3); ib <- sort(sample(ib, 800)) }
      if (length(cur)) {
        T <- sapply(names(cur), function(v) approx(cur[[v]]$x, cur[[v]]$y, E[[v]][ib], rule = 2, ties = mean)$y)
        Xb <- r$X[ib, , drop = FALSE]; Xb <- Xb[, colSums(Xb) > 0, drop = FALSE]
        out$rho <- cor(as.vector(dist(T)), as.vector(sqrt(betadiver(Xb, "sim"))), method = "spearman")
      }
    }
  }
  saveRDS(out, saida); saida
}
cl <- makeCluster(12); clusterExport(cl, c("raiz", "dir_c", "ajusta"))
t0 <- Sys.time(); feitos <- unlist(parLapplyLB(cl, split(jobs, seq_len(nrow(jobs))), function(j) tryCatch(ajusta(j), error = function(e) paste("ERRO", basename(j$arq), j$bacia, j$tipo, conditionMessage(e)))))
stopCluster(cl); print(feitos[grepl("^ERRO", feitos)]); cat("tempo:", format(Sys.time() - t0), "\n")
rs <- lapply(feitos[!grepl("^ERRO", feitos)], readRDS)

# ---- reference: global model (88) on the same basin units ----
curvas <- read.csv(file.path(raiz, "saida/etapa3_n1_gf_curvas.csv"))
rho_global <- bind_rows(lapply(unique(jobs$arq), function(a) { r <- readRDS(a); cv <- curvas[curvas$grupo == r$grupo & curvas$comp == r$comp, ]
  bind_rows(lapply(unique(jobs$bacia[jobs$arq == a]), function(b) { ib <- which(r$meta$bacia == b); if (length(ib) > 800) { set.seed(3); ib <- sort(sample(ib, 800)) }
    vs <- intersect(unique(cv$preditor), names(r$E))
    T <- sapply(vs, function(v) { c1 <- cv[cv$preditor == v, ]; approx(c1$x, c1$y, r$E[[v]][ib], rule = 2, ties = mean)$y })
    Xb <- r$X[ib, , drop = FALSE]; Xb <- Xb[, colSums(Xb) > 0, drop = FALSE]
    data.frame(comp = r$comp, grupo = r$grupo, bacia = b, rho_global = cor(as.vector(dist(T)), as.vector(sqrt(vegan::betadiver(Xb, "sim"))), method = "spearman")) })) }))

aj <- bind_rows(lapply(rs, function(o) data.frame(comp = o$comp, grupo = o$grupo, bacia = o$bacia, tipo = o$tipo, unidades = o$unidades, taxons = o$taxons,
  ajustados = o$ajustados, r2 = o$r2, rho = o$rho)))
ajw <- aj |> select(comp, grupo, bacia, tipo, unidades, ajustados, r2, rho) |> tidyr::pivot_wider(names_from = tipo, values_from = c(unidades, ajustados, r2, rho)) |>
  left_join(rho_global, by = c("comp", "grupo", "bacia")) |> mutate(razao_transf = rho_fora / rho_local, nome_bacia = bac_nome$bacia[match(bacia, bac_nome$codigo)]) |> as.data.frame()
write.csv(ajw, file.path(raiz, "saida/etapa3_n1_transf_ajustes.csv"), row.names = FALSE)
ll <- bind_rows(lapply(rs, function(o) if (o$tipo == "local" && !is.null(o$lim)) cbind(comp = o$comp, grupo = o$grupo, bacia = o$bacia, o$lim)))
write.csv(ll, file.path(raiz, "saida/etapa3_n1_transf_limiares_locais.csv"), row.names = FALSE)

# ---- robust thresholds x basins ----
rob <- read.csv(file.path(raiz, "saida/etapa3_n1_limiares_robustos.csv")); rob$id <- seq_len(nrow(rob))
amp <- bind_rows(lapply(unique(jobs$arq), function(a) { r <- readRDS(a)
  bind_rows(lapply(names(r$E)[!grepl("^bac_", names(r$E))], function(v) { x <- r$E[[v]]
    bind_rows(lapply(unique(jobs$bacia[jobs$arq == a]), function(b) { xb <- x[r$meta$bacia == b]
      data.frame(comp = r$comp, grupo = r$grupo, preditor = v, bacia = b, amp = diff(range(x)), min_b = min(xb), max_b = max(xb)) })) })) }))
locais <- aj |> filter(tipo == "local") |> select(comp, grupo, bacia, r2_local = r2)
testes <- bind_rows(lapply(seq_len(nrow(rob)), function(k) {
  g <- regmatches(rob$lista[k], gregexpr("[^,]+ \\((bentos|pelagico)\\)", rob$lista[k]))[[1]]; g <- trimws(g)
  d <- data.frame(grupo = sub(" \\((bentos|pelagico)\\)$", "", g), comp = sub("^.* \\((bentos|pelagico)\\)$", "\\1", g))
  d |> inner_join(amp |> filter(preditor == rob$preditor[k]), by = c("comp", "grupo")) |> inner_join(locais, by = c("comp", "grupo", "bacia")) |>
    mutate(id = k, preditor = rob$preditor[k], limiar = rob$limiar[k], testavel = min_b <= limiar - 0.05 * amp & max_b >= limiar + 0.05 * amp,
      confirmado = testavel & vapply(seq_along(bacia), function(i) { l <- ll[ll$comp == comp[i] & ll$grupo == grupo[i] & ll$bacia == bacia[i] & ll$preditor == rob$preditor[k], ]
        nrow(l) > 0 && any(abs(l$limiar - rob$limiar[k]) <= 0.10 * amp[i]) }, logical(1)))
}))
write.csv(testes, file.path(raiz, "saida/etapa3_n1_transf_testes.csv"), row.names = FALSE)
cls <- testes |> group_by(id) |> summarise(bacias_testaveis = n_distinct(bacia[testavel]), bacias_confirmadas = n_distinct(bacia[confirmado]),
  testes = sum(testavel), confirmados = sum(confirmado),
  onde = paste(sort(unique(bac_nome$bacia[match(bacia[confirmado], bac_nome$codigo)])), collapse = ", "), .groups = "drop")
rob2 <- rob |> left_join(cls, by = "id") |> mutate(across(c(bacias_testaveis, bacias_confirmadas, testes, confirmados), ~ ifelse(is.na(.x), 0L, .x)),
  classe = case_when(bacias_testaveis == 0 ~ "não testável", bacias_confirmadas >= 2 & bacias_confirmadas / bacias_testaveis >= 2/3 ~ "global",
                     bacias_confirmadas >= 1 ~ "regional", TRUE ~ "não confirmado")) |> select(-id)
write.csv(rob2, file.path(raiz, "saida/etapa3_n1_transf_robustos.csv"), row.names = FALSE)
options(width = 200)
print(ajw |> select(comp, grupo, nome_bacia, unidades_local, r2_local, rho_local, rho_fora, rho_global, razao_transf) |> mutate(across(where(is.double), ~ round(.x, 3))))
print(rob2 |> select(preditor, limiar, grupos, bacias_testaveis, bacias_confirmadas, classe, onde) |> mutate(limiar = signif(limiar, 4)))
print(table(rob2$classe)); cat(format(Sys.time()), "FIM\n")
