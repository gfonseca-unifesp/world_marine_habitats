# funcoes_etapa2.R  (Phase 8, stage 2: class diagnostics)
# Unit preparation (same rule as stage 1: purity >= 0.8, >= 10 events, class with >= 5 units, 10 events
# drawn with the seed of draw 1), leave-one-block-out classification and search for discontinuities.
suppressPackageStartupMessages({library(arrow); library(dplyr); library(vegan); library(cluster)})
EV_MIN <- 10; UN_MIN <- 5; PUR_MIN <- 0.8; ORD <- c("25km" = 1, "12.5km" = 2, "6.25km" = 3)

# center of cell "res_linha_coluna" in EASE-Grid 2.0 (m) and in lon/lat
centro_celula <- function(cell) {
  p <- do.call(rbind, strsplit(sub(" .*$", "", cell), "_")); x <- y <- rep(NA_real_, nrow(p))
  for (r in unique(p[, 1])) { i <- p[, 1] == r; g <- ease2_grade(r)
    x[i] <- terra::xFromCol(g, as.integer(p[i, 3])); y[i] <- terra::yFromRow(g, as.integer(p[i, 2])) }
  ll <- terra::project(cbind(x, y), from = "EPSG:6933", to = "EPSG:4326")
  data.frame(x = x, y = y, lon = ll[, 1], lat = ll[, 2])
}

# units and incidence (draw 1) of a level 1 contrast (task from 81_etapa1_n1.R)
prepara_n1 <- function(tk, raiz) {
  d <- open_dataset(file.path(raiz, sprintf("saida/tabela_certa_%s_n1.parquet", tk$comp))) |> filter(grupo == !!tk$grupo)
  if (tk$comp == "bentos") {
    col <- if (tk$nivel == "EBU grupo") c("ebu_grupo", "pureza_grupo") else c("ebu_regiao", "pureza_regiao")
    d <- d |> select(cell = cell_id, evento = event_id, taxon_unidade, bloco = bloco_espacial, classe = all_of(col[1]), pureza = all_of(col[2])) |>
      collect() |> mutate(unidade = cell)
  } else {
    d <- d |> select(cell = cell_id, evento = event_id, taxon_unidade, bloco = bloco_espacial, fatia, emu, pureza = pureza_emu)
    if (tk$fatia != "estatico") { f <- tk$fatia; d <- d |> filter(fatia == f) }
    d <- d |> collect() |> mutate(classe = as.character(emu), unidade = paste(cell, emu))
  }
  incidencia(d)
}

# same, level 2 (task from 82_etapa1_n2.R), blocks of bloco_km
prepara_n2 <- function(tk, raiz, bloco_km = 100) {
  s <- sub("km", "", tk$res_escolhida); okres <- names(ORD)[ORD >= ORD[tk$res_escolhida]]
  d <- open_dataset(file.path(raiz, sprintf("saida/tabela_certa_%s_n2.parquet", tk$comp))) |>
    filter(regiao == !!tk$regiao, grupo_ecologico == !!tk$grupo_ecologico, res_max %in% okres)
  if (tk$comp == "bentos") {
    nv <- sub("ebu_", "", tk$nivel); fcol <- if (tk$fatia == "estatico") "fatia_estatica" else "fatia_trimestre"
    d <- d |> select(cell = all_of(paste0("cell_", s)), r6, c6, evento, taxon_unidade, fatia = all_of(fcol),
                     classe = all_of(paste0("ebu_", nv, "_", s)), pureza = all_of(paste0("pureza_", nv, "_", s))) |>
      filter(fatia == !!tk$fatia) |> collect() |> mutate(unidade = cell)
  } else {
    d <- d |> select(cell = all_of(paste0("cell_", s)), r6, c6, evento, taxon_unidade, fatia = fatia_mes,
                     classe = all_of(paste0("emu_", s)), pureza = all_of(paste0("pureza_emu_", s)), situacao = all_of(paste0("situacao_emu_", s))) |>
      filter(fatia == !!tk$fatia, situacao == "ok") |> collect() |> mutate(classe = as.character(classe), unidade = paste(cell, classe))
  }
  lado <- bloco_km / 6.25
  d <- d |> mutate(bloco = paste((r6 - 1) %/% lado, (c6 - 1) %/% lado, sep = "_"))
  incidencia(d)
}

incidencia <- function(d) {
  d <- d |> filter(!is.na(classe), pureza >= PUR_MIN, !is.na(taxon_unidade), taxon_unidade != "") |> mutate(classe = as.character(classe))
  un <- d |> group_by(unidade, classe, cell) |> summarise(ev = n_distinct(evento), bloco = names(which.max(table(bloco))), .groups = "drop") |>
    filter(ev >= EV_MIN)
  ca <- un |> count(classe) |> filter(n >= UN_MIN); un <- un |> filter(classe %in% ca$classe) |> arrange(unidade)
  if (nrow(ca) < 2) return(NULL)
  ev <- d |> filter(unidade %in% un$unidade) |> distinct(unidade, evento, taxon_unidade)
  set.seed(1001)
  sor <- ev |> distinct(unidade, evento) |> group_by(unidade) |> slice_sample(n = EV_MIN) |> ungroup()
  m <- ev |> semi_join(sor, by = c("unidade", "evento")) |> distinct(unidade, taxon_unidade) |> mutate(v = 1L) |>
    tidyr::pivot_wider(names_from = taxon_unidade, values_from = v, values_fill = 0L) |> arrange(unidade)
  stopifnot(identical(m$unidade, un$unidade))
  X <- as.matrix(m[, -1]); X <- X[, colSums(X) > 0, drop = FALSE]
  list(meta = cbind(un, centro_celula(un$cell)), X = X)
}

# Classification by nearest fauna centroid, leaving out the unit's block.
# S = matrix of squared dissimilarities (here beta_sim, because D = sqrt(beta_sim)). Distance to the centroid of class c
# (Anderson 2006): d2(i, c) = mean_j S_ij - mean_jk S_jk / 2, with j, k in c outside the block of i.
# A class without >= un_min units outside the block does not compete; a unit whose class does not compete is "nao classificavel" (unclassifiable).
classifica_lbo <- function(S, grupo, bloco, un_min = 5) {
  grupo <- factor(grupo); cls <- levels(grupo); Z <- model.matrix(~ grupo - 1); colnames(Z) <- cls
  Q <- S %*% Z; tot_w <- colSums(Z * Q)                      # tot_w[c] = sum of S_jk with j, k in c
  pred <- rep(NA_character_, length(grupo)); dist_propria <- rep(NA_real_, length(grupo))
  for (b in unique(bloco)) {
    ib <- which(bloco == b); Zb <- Z[ib, , drop = FALSE]
    m <- colSums(Z) - colSums(Zb)                               # members of each class outside the block
    soma_i <- Q[ib, , drop = FALSE] - S[ib, ib, drop = FALSE] %*% Zb      # sum_j S_ij, j in c outside the block
    wb <- tot_w - 2 * colSums(Zb * Q[ib, , drop = FALSE]) + colSums(Zb * (S[ib, ib, drop = FALSE] %*% Zb))
    d2 <- sweep(soma_i, 2, m, "/") - matrix(wb / (2 * m^2), length(ib), length(cls), byrow = TRUE)
    d2[, m < un_min] <- Inf
    pred[ib] <- cls[apply(d2, 1, which.min)]
    pred[ib][apply(d2, 1, function(v) all(!is.finite(v)))] <- NA
    k <- cbind(seq_along(ib), as.integer(grupo[ib])); dist_propria[ib] <- d2[k]
    pred[ib][!is.finite(d2[k])] <- NA                            # the unit's own class does not compete
  }
  data.frame(classe = as.character(grupo), prevista = pred, classificavel = !is.na(pred), correta = pred == as.character(grupo),
             dist2_propria = dist_propria)
}

metricas_classif <- function(cl) {
  c2 <- cl[cl$classificavel, ]; if (!nrow(c2)) return(NULL)
  tab <- table(factor(c2$classe, levels = sort(unique(cl$classe))), factor(c2$prevista, levels = sort(unique(cl$classe))))
  n <- sum(tab); acc <- sum(diag(tab)) / n; pe <- sum(rowSums(tab) * colSums(tab)) / n^2
  sens <- diag(tab) / pmax(rowSums(tab), 1); prec <- diag(tab) / pmax(colSums(tab), 1)
  list(geral = data.frame(classificaveis = n, nao_classificaveis = sum(!cl$classificavel), acerto = acc,
                          acerto_equilibrado = mean(sens[rowSums(tab) > 0]), acaso = pe, kappa = (acc - pe) / (1 - pe)),
       por_classe = data.frame(classe = rownames(tab), unidades = as.integer(rowSums(tab)), sensibilidade = as.numeric(sens),
                               precisao = as.numeric(prec), confunde_com = vapply(rownames(tab), function(r) {
                                 v <- tab[r, ]; v[r] <- 0; if (sum(v) == 0) "" else paste(names(sort(v[v > 0], decreasing = TRUE))[seq_len(min(3, sum(v > 0)))], collapse = ", ") }, "")),
       confusao = as.data.frame(tab, responseName = "n") |> setNames(c("classe", "prevista", "n")))
}

# Discontinuities within a class: PAM on sqrt(beta_sim) with k = 2 to 4, smallest subgroup >= 5 units and >= 10% of the class (n > 1500: medoids on a sample of 1500 and
# assignment of the remaining units to the nearest medoid). Best k by mean silhouette. Spatial coherence: fraction of units whose
# nearest geographic neighbor (in the same class) is in the same subgroup, against 999 permutations of the labels.
# "descontinuidade" = silhouette >= 0.25, coherence with p < 0.05 and excess (observed - expected) >= 0.10.
descontinuidade <- function(Dm, xy, n_min = 20, kmax = 4, nmax = 1500) {
  n <- nrow(Dm); if (n < n_min) return(NULL)
  idx <- if (n > nmax) sort(sample(n, nmax)) else seq_len(n)
  melhor <- NULL
  for (k in 2:kmax) {
    if (length(idx) <= k + 1) break
    pm <- pam(as.dist(Dm[idx, idx]), k = k, diss = TRUE)
    rot <- if (n > nmax) apply(Dm[, idx[pm$id.med], drop = FALSE], 1, which.min) else pm$clustering
    if (min(table(rot)) < max(5, 0.1 * n)) next                 # subgroups of a few isolated units do not count
    sil <- if (n > nmax) pm$silinfo$avg.width else mean(silhouette(rot, as.dist(Dm))[, 3])
    if (is.null(melhor) || sil > melhor$sil) melhor <- list(k = k, sil = sil, rot = rot)
  }
  if (is.null(melhor)) return(list(resumo = data.frame(unidades = n, k = 1L, silhueta = NA, coerencia = NA, coerencia_esperada = NA,
    p_coerencia = NA, tamanhos = as.character(n), descontinuidade = FALSE), subgrupo = rep(1L, n)))
  nn <- vapply(seq_len(n), function(i) { d <- (xy[, 1] - xy[i, 1])^2 + (xy[, 2] - xy[i, 2])^2; d[i] <- Inf; which.min(d) }, 1L)
  obs <- mean(melhor$rot == melhor$rot[nn])
  perm <- replicate(999, { r <- sample(melhor$rot); mean(r == r[nn]) })
  esp <- mean(perm); p <- (sum(perm >= obs) + 1) / 1000
  list(resumo = data.frame(unidades = n, k = melhor$k, silhueta = melhor$sil, coerencia = obs, coerencia_esperada = esp, p_coerencia = p,
                           tamanhos = paste(sort(table(melhor$rot), decreasing = TRUE), collapse = "/"),
                           descontinuidade = melhor$sil >= 0.25 && p < 0.05 && (obs - esp) >= 0.10),
       subgrupo = melhor$rot)
}
