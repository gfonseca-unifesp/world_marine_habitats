# funcoes_permanova.R  Fast PERMANOVA (same statistic as vegan::adonis2), for large n.
# G = centered Gower matrix (-1/2 D^2, doubly centered). For class after block:
#   SS_classe|bloco = tr[(W'W)^-1 W' G W], W = class indicators centered within each block (without the last column).
# The block SS does not change under permutations within blocks, so each permutation costs only G %*% W (n x n x k).
# Without block (global): W = indicators centered on the grand mean.
# Returns R2 (= SS_classe / SS_total, as in adonis2 by = "terms"), F, p and the number of permutations.
permanova_rapida <- function(D, grupo, bloco = NULL, nperm = 999, controle = NULL) {
  grupo <- droplevels(factor(grupo)); n <- length(grupo); k <- nlevels(grupo)
  A <- -0.5 * as.matrix(D)^2; G <- sweep(sweep(A, 1, rowMeans(A)), 2, colMeans(A)) + mean(A); rm(A)
  sst <- sum(diag(G))
  if (!is.finite(sst) || sst < 1e-12) return(list(R2 = NA, omega2 = NA, F = NA, p = NA, nperm = 0L, gl = c(NA, NA), R2_bloco = NA))   # no variation
  if (is.null(bloco)) { bloco <- factor(rep(1, n)) } else bloco <- droplevels(factor(bloco))
  nb <- nlevels(bloco); Bm <- if (nb > 1) model.matrix(~ bloco - 1) else matrix(1, n, 1)
  ssb <- if (nb > 1) sum(vapply(seq_len(nb), function(j) { i <- which(Bm[, j] == 1); sum(G[i, i]) / length(i) }, 0)) else 0
  ss_g <- function(g) {
    Z <- model.matrix(~ g - 1)[, -k, drop = FALSE]
    W <- Z - Bm %*% (crossprod(Bm, Z) / colSums(Bm))            # center within blocks
    WtW <- crossprod(W)
    sum(diag(MASS::ginv(WtW) %*% crossprod(W, G %*% W)))
  }
  gl_g <- qr(cbind(Bm, model.matrix(~ grupo - 1)))$rank - nb; gl_r <- n - nb - gl_g
  Fst <- function(ss) (ss / gl_g) / ((sst - ssb - ss) / gl_r)
  ss0 <- ss_g(grupo); F0 <- Fst(ss0)
  if (is.null(controle)) controle <- if (nb > 1) permute::how(blocks = bloco, nperm = nperm) else permute::how(nperm = nperm)
  P <- permute::shuffleSet(n, control = controle)
  Fp <- apply(P, 1, function(o) Fst(ss_g(grupo[o])))
  list(R2 = ss0 / sst, omega2 = omega2_parcial(F0, gl_g, n), F = F0, p = (sum(Fp >= F0 - 1e-8) + 1) / (nrow(P) + 1), nperm = nrow(P),
       gl = c(gl_g, gl_r), R2_bloco = ssb / sst)
}

# partial omega2 (Olejnik and Algina 2003): R2 corrected for the number of units and classes; comparable across contrasts
# of different sizes (plain R2 grows with the number of classes even without an effect). Can be < 0 when F < 1.
omega2_parcial <- function(F, gl, n) gl * (F - 1) / (gl * (F - 1) + n)

# Equal-effort test: in each repetition, n_classe units drawn per class (all of them, if the class has fewer),
# only blocks that still have >= 2 classes; PERMANOVA class after block, permutations within blocks.
# Summary: median of n, R2 and omega2; fraction of repetitions with p < 0.05 (power at equal effort) and median of p.
permanova_esforco_igual <- function(M, grupo, bloco, n_classe = 10, B = 200, nperm = 199) {
  grupo <- droplevels(factor(grupo)); bloco <- droplevels(factor(bloco))
  r <- t(vapply(seq_len(B), function(b) {
    idx <- unlist(lapply(levels(grupo), function(g) { i <- which(grupo == g); if (length(i) <= n_classe) i else sample(i, n_classe) }))
    mix <- tapply(as.character(grupo[idx]), as.character(bloco[idx]), function(z) length(unique(z)) >= 2)
    idx <- idx[as.character(bloco[idx]) %in% names(mix)[mix]]
    if (length(idx) < 4 || nlevels(droplevels(grupo[idx])) < 2) return(c(n = NA, R2 = NA, omega2 = NA, p = NA))
    x <- permanova_rapida(as.dist(M[idx, idx]), grupo[idx], bloco[idx], nperm)
    c(n = length(idx), R2 = x$R2, omega2 = x$omega2, p = x$p)
  }, c(n = 0, R2 = 0, omega2 = 0, p = 0)))
  data.frame(n_classe = n_classe, repeticoes = sum(!is.na(r[, "p"])), n_med = median(r[, "n"], na.rm = TRUE),
             R2_eq = median(r[, "R2"], na.rm = TRUE), omega2_eq = median(r[, "omega2"], na.rm = TRUE),
             omega2_eq_q05 = unname(quantile(r[, "omega2"], 0.05, na.rm = TRUE)), omega2_eq_q95 = unname(quantile(r[, "omega2"], 0.95, na.rm = TRUE)),
             poder_eq = mean(r[, "p"] < 0.05, na.rm = TRUE), p_eq_med = median(r[, "p"], na.rm = TRUE))
}

# Pairwise contrasts between classes (habitats), within blocks. For each pair (a, b): units of the two classes in the
# blocks that have both; >= un_min units of each. omega2 with all the data (permutations within blocks) and test with
# equal effort (n_classe units per class, B repetitions). "separa" = omega2 >= 0.01 and power >= 0.8 (Gate 4 criterion).
# pares: 2-column matrix with the pairs to test (default: all)
pares_classes <- function(M, grupo, bloco, un_min = 5, n_classe = 10, B = 100, nperm = 199, pares = NULL) {
  grupo <- as.character(grupo); bloco <- as.character(bloco); cls <- sort(unique(grupo))
  if (length(cls) < 2) return(NULL)
  pr <- if (is.null(pares)) t(combn(cls, 2)) else pares; if (!nrow(pr)) return(NULL)
  bind_rows(lapply(seq_len(nrow(pr)), function(j) {
    a <- pr[j, 1]; b <- pr[j, 2]
    comuns <- intersect(bloco[grupo == a], bloco[grupo == b])
    i <- which(grupo %in% c(a, b) & bloco %in% comuns); na <- sum(grupo[i] == a); nb <- sum(grupo[i] == b)
    base <- data.frame(classe_a = a, classe_b = b, unidades_a = na, unidades_b = nb, blocos = length(comuns))
    if (na < un_min || nb < un_min) return(cbind(base, omega2 = NA, omega2_eq = NA, poder_eq = NA, separa = NA))
    x <- permanova_rapida(as.dist(M[i, i]), grupo[i], bloco[i], nperm)
    e <- permanova_esforco_igual(M[i, i], grupo[i], bloco[i], n_classe = n_classe, B = B, nperm = nperm)
    cbind(base, omega2 = x$omega2, omega2_eq = e$omega2_eq, poder_eq = e$poder_eq,
          separa = isTRUE(x$omega2 >= 0.01) && isTRUE(e$poder_eq >= 0.8))
  }))
}