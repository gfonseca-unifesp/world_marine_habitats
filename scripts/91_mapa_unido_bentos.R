# 91_mapa_unido_bentos.R  (Phase 8; Gate 9, 30/10/2026: final benthic map uniting top-down and bottom-up)
# Principle: separate when any evidence separates; merge only when the fauna shows equivalence with power.
#  1 candidates: each cell = depth band x new habitat (HB, 89b) x EBU region;
#  2 merge by fauna: pairs of candidates tested by taxonomic group (pares_classes: omega2 >= 0.01 and power >= 0.8);
#    consensus (equivalent: not distinct in >= 2/3 of >= 2 groups) and complete-linkage clustering with a 1/3 cut; an untested pair
#    is not merged;
#  3 protected features: EBU groups C (seamounts) and D (ridges and rifts) only merge with the same group;
#  4 hierarchy: I compartment, II depth band (0-200, 200-2,000, > 2,000 m), III HB, IV EBU region; merging only within the depth band;
#  5 minimum area: class with < 10 cells, not protected and without a fauna-distinct pair, goes to the nearest neighboring class (4-cell
#    adjacency) of the same depth band in Gradient Forest space; repeat until stable;
#  6 attributes: evidence (distinct by fauna / merged by fauna / untested = hatched), protected, shelf seasonality
#    (annual amplitude of bottom temperature, Gate 8b item 2);
#  7 acceptance: omega2 of the final classification >= that of the best component (HB, EBU region) on the same units, per group;
#    all depth bands and EBU groups represented.
#  8 the pelagic goes in a separate layer (after level 2, Gate 8 item 3).
# Outputs: saida/mapa_unido_bentos_{classes,pares,avaliacao,representacao}.csv, saida/mapa_unido_bentos_25km.tif,
#         saida/mapa_unido_bentos_25km.png
suppressPackageStartupMessages({library(dplyr); library(terra); library(arrow); library(parallel)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); source(file.path(raiz, "scripts/funcoes_grade.R")); source(file.path(raiz, "scripts/funcoes_permanova.R"))
source(file.path(raiz, "scripts/funcoes_mapa.R"))
dir_t <- file.path(raiz, "intermediario/etapa3_n1"); NOME_FX <- c("0–200 m", "200–2.000 m", "> 2.000 m"); PROT <- c("C", "D"); AREA_MIN <- 10
faixa <- function(prof) findInterval(prof, c(200, 2000), left.open = TRUE) + 1
g25 <- ease2_grade("25km"); M <- lapply(1:12, function(m) rast(file.path(raiz, sprintf("camadas/glorys/mes/glorys_%02d_25km.tif", m))))
vm <- function(camada) sapply(1:12, function(k) values(M[[k]][[camada]])[, 1])
curvas <- read.csv(file.path(raiz, "saida/etapa3_n1_gf_curvas.csv"))
transf <- function(E, cp) { cv <- curvas[curvas$comp == cp, ]; vs <- intersect(unique(cv$preditor), names(E))   # as in 89b
  sapply(vs, function(v) { s <- 0; for (g in unique(cv$grupo[cv$preditor == v])) { c1 <- cv[cv$grupo == g & cv$preditor == v, ]
    s <- s + approx(c1$x, c1$y, E[[v]], rule = 2, ties = mean)$y }; s }) }

# ---------------- 1 candidates ----------------
e25 <- read_parquet(file.path(raiz, "camadas/ebu/ebu_celula_25km.parquet"))
p <- do.call(rbind, strsplit(e25$cell_id, "_")); cel_b <- cellFromRowCol(g25, as.integer(p[, 2]), as.integer(p[, 3]))
lh <- read.csv(file.path(raiz, "saida/habitats_novos_legenda.csv")) |> filter(arquivo == "bentos_25km")
hb <- values(rast(file.path(raiz, "saida/habitats_novos_bentos_25km.tif")))[cel_b, 1]
tf <- vm("tfundo")[cel_b, ]; vf <- vm("vfundo")[cel_b, ]; ss <- vm("thetao_0")[cel_b, ]; ml <- vm("mld")[cel_b, ]; ge <- vm("gelo")[cel_b, ]
C <- data.frame(cel = cel_b, prof = values(M[[1]][["prof_glorys"]])[cel_b, 1], tfundo_med = rowMeans(tf), tfundo_min = apply(tf, 1, min), tfundo_max = apply(tf, 1, max),
  vfundo_med = rowMeans(vf), vfundo_max = apply(vf, 1, max), sst_med = rowMeans(ss), sst_amp = apply(ss, 1, max) - apply(ss, 1, min),
  mld_max = apply(ml, 1, max), gelo_meses = rowSums(ge >= 0.15, na.rm = TRUE), tfundo_amp = apply(tf, 1, max) - apply(tf, 1, min),
  hb = lh$habitat[hb], regiao = e25$ebu_regiao, grupo_ebu = e25$ebu_grupo)
rm(tf, vf, ss, ml, ge); C <- C[complete.cases(C[, 2:11]) & !is.na(C$hb) & !is.na(C$regiao), ]
C$faixa <- faixa(C$prof); C$cand <- paste(C$faixa, C$hb, C$regiao, sep = "|")
TC <- transf(C, "bentos")
cand <- C |> group_by(cand) |> summarise(faixa = first(faixa), hb = first(hb), regiao = first(regiao), grupo_ebu = first(grupo_ebu), celulas = n(), .groups = "drop")
cat(format(Sys.time()), "celulas:", nrow(C), "| candidatos:", nrow(cand), "| com < ", AREA_MIN, "celulas:", sum(cand$celulas < AREA_MIN), "\n")
pode <- function(a, b) { x <- cand[match(a, cand$cand), ]; y <- cand[match(b, cand$cand), ]   # rules 3 and 4
  x$faixa == y$faixa & (!(x$grupo_ebu %in% PROT | y$grupo_ebu %in% PROT) | x$grupo_ebu == y$grupo_ebu) }

# ---------------- 2 merge by fauna ----------------
arqs <- list.files(dir_t, pattern = "^bentos_.*\\.rds$", full.names = TRUE)
prep <- lapply(arqs, function(a) { r <- readRDS(a); p <- do.call(rbind, strsplit(rownames(r$X), "_"))
  cc <- C$cand[match(cellFromRowCol(g25, as.integer(p[, 2]), as.integer(p[, 3])), C$cel)]
  n <- table(cc); cl <- names(n)[n >= 5]; pr <- if (length(cl) >= 2) t(combn(sort(cl), 2)) else matrix(character(0), 0, 2)
  pr <- pr[pode(pr[, 1], pr[, 2]), , drop = FALSE]
  list(arq = a, grupo = r$grupo, cand = cc, pares = pr) })
jobs <- bind_rows(lapply(seq_along(prep), function(k) { n <- nrow(prep[[k]]$pares); if (!n) return(NULL)
  data.frame(k = k, de = seq(1, n, by = 150)) |> mutate(ate = pmin(de + 149, n)) }))
cat(format(Sys.time()), "pares a testar:", sum(sapply(prep, function(x) nrow(x$pares))), "em", nrow(jobs), "lotes\n")
testa <- function(j) {
  suppressPackageStartupMessages({library(vegan); library(dplyr)}); source(file.path(raiz, "scripts/funcoes_permanova.R"))
  pp <- prep[[j$k]]; r <- readRDS(pp$arq); ok <- !is.na(pp$cand)
  Dm <- as.matrix(sqrt(betadiver(r$X[ok, , drop = FALSE], "sim"))); set.seed(7000 + j$de)
  x <- pares_classes(Dm, pp$cand[ok], r$meta$bloco[ok], un_min = 5, n_classe = 10, B = 100, pares = pp$pares[j$de:j$ate, , drop = FALSE])
  if (!is.null(x)) cbind(grupo = pp$grupo, x)
}
cl <- makeCluster(12); clusterExport(cl, c("raiz", "prep"))
pares <- bind_rows(parLapplyLB(cl, split(jobs, seq_len(nrow(jobs))), testa)); stopCluster(cl)
pares <- pares |> filter(!is.na(separa))
cat(format(Sys.time()), "pares testados:", nrow(pares), "\n")
consenso <- pares |> group_by(classe_a, classe_b) |> summarise(grupos = n(), distintos = sum(separa), .groups = "drop") |>
  mutate(frac = distintos / grupos, consenso = ifelse(grupos >= 2 & frac <= 1/3, "equivalente", ifelse(frac >= 2/3, "distinto", "depende do grupo")))
# only candidates that entered some equivalent pair can merge; the others stay alone (spares the full matrix)
eq <- consenso |> filter(consenso == "equivalente"); nos <- sort(unique(c(eq$classe_a, eq$classe_b)))
cand$classe <- cand$cand
if (length(nos) >= 2) {
  Mx <- matrix(1, length(nos), length(nos), dimnames = list(nos, nos)); diag(Mx) <- 0
  cs <- consenso |> filter(classe_a %in% nos, classe_b %in% nos)
  Mx[cbind(cs$classe_a, cs$classe_b)] <- ifelse(cs$grupos >= 2, cs$frac, 1); Mx[cbind(cs$classe_b, cs$classe_a)] <- Mx[cbind(cs$classe_a, cs$classe_b)]
  pd <- outer(nos, nos, pode); Mx[!pd] <- 1
  g <- cutree(hclust(as.dist(Mx), method = "complete"), h = 1/3)
  rep_ <- tapply(names(g), g, function(x) x[1]); cand$classe[match(nos, cand$cand)] <- rep_[as.character(g)]
}
cand$fundido_fauna <- ave(cand$cand, cand$classe, FUN = length) > 1
dist_c <- unique(c(consenso$classe_a[consenso$consenso == "distinto"], consenso$classe_b[consenso$consenso == "distinto"]))
testado <- unique(c(consenso$classe_a, consenso$classe_b))

# ---------------- 5 minimum area ----------------
C$classe <- cand$classe[match(C$cand, cand$cand)]
cls_info <- function() C |> group_by(classe) |> summarise(celulas = n(), faixa = first(faixa), prot = any(grupo_ebu %in% PROT), .groups = "drop")
centro <- function() { m <- rowsum(TC, C$classe) / as.vector(table(C$classe)[rownames(rowsum(TC, C$classe))]); m }
distinta <- function(cl) any(cand$cand[cand$classe == cl] %in% dist_c)
absorvidas <- list(); repeat {
  ci <- cls_info(); peq <- ci |> filter(celulas < AREA_MIN, !prot) |> arrange(celulas); peq <- peq[!vapply(peq$classe, distinta, TRUE), ]
  if (!nrow(peq)) break
  CT <- centro(); mudou <- FALSE; r <- rast(g25); vv <- rep(NA_character_, ncell(r)); vv[C$cel] <- C$classe
  for (k in seq_len(nrow(peq))) { a <- peq$classe[k]; ca <- C$cel[C$classe == a]; if (!length(ca)) next
    viz <- unique(vv[as.vector(adjacent(r, ca, "rook"))]); viz <- setdiff(viz[!is.na(viz)], a)
    viz <- viz[ci$faixa[match(viz, ci$classe)] == peq$faixa[k] & !ci$prot[match(viz, ci$classe)]]
    if (!length(viz)) next
    d <- colSums((t(CT[viz, , drop = FALSE]) - CT[a, ])^2); dest <- viz[which.min(d)]
    C$classe[C$classe == a] <- dest; vv[ca] <- dest; cand$classe[cand$classe == a] <- dest; absorvidas[[length(absorvidas) + 1]] <- data.frame(de = a, para = dest); mudou <- TRUE }
  if (!mudou) break
}
absorvidas <- bind_rows(absorvidas); cat(format(Sys.time()), "classes absorvidas por area:", nrow(absorvidas), "\n")

# ---------------- 6 attributes ----------------
fin <- C |> group_by(classe) |> summarise(celulas = n(), faixa = NOME_FX[first(faixa)], hb = paste(sort(unique(hb)), collapse = " + "),
  regioes = paste(sort(unique(regiao)), collapse = " + "), grupos_ebu = paste(sort(unique(grupo_ebu)), collapse = ""),
  protegida = any(grupo_ebu %in% PROT), tfundo_amp = median(tfundo_amp), prof = median(prof), .groups = "drop") |>
  left_join(cand |> group_by(classe) |> summarise(candidatos = n(), cand_fauna = sum(cand %in% testado), fundida_fauna = any(fundido_fauna),
    distinta_fauna = any(cand %in% dist_c), .groups = "drop"), by = "classe") |>
  mutate(evidencia = case_when(distinta_fauna ~ "distinta pela fauna", fundida_fauna ~ "fundida pela fauna", cand_fauna > 0 ~ "testada, sem par conclusivo", TRUE ~ "não testada"),
    standby = evidencia %in% c("não testada", "testada, sem par conclusivo"),
    sazonalidade = ifelse(faixa == "0–200 m", ifelse(tfundo_amp >= 5, "forte", ifelse(tfundo_amp >= 2, "moderada", "fraca")), "–")) |>
  arrange(match(faixa, NOME_FX), hb, regioes)
fin$codigo <- paste0("HF", seq_len(nrow(fin))); C$codigo <- fin$codigo[match(C$classe, fin$classe)]
cat("classes finais:", nrow(fin), "| por evidencia:", paste(names(table(fin$evidencia)), table(fin$evidencia), collapse = "; "), "\n")

# ---------------- 7 acceptance ----------------
aval <- bind_rows(lapply(prep, function(pp) { r <- readRDS(pp$arq); ok <- !is.na(pp$cand)
  cel <- C$cel[match(pp$cand[ok], C$cand)]; Dm <- as.matrix(sqrt(vegan::betadiver(r$X[ok, , drop = FALSE], "sim")))
  cls <- list("mapa unido" = C$codigo[match(pp$cand[ok], C$cand)], "habitats novos (HB)" = cand$hb[match(pp$cand[ok], cand$cand)],
              "EBU região" = cand$regiao[match(pp$cand[ok], cand$cand)])
  bind_rows(lapply(names(cls), function(nm) { g <- cls[[nm]]; k <- g %in% names(which(table(g) >= 5)); if (length(unique(g[k])) < 2) return(NULL)
    set.seed(9); x <- permanova_rapida(as.dist(Dm[k, k]), g[k], r$meta$bloco[ok][k], 199)
    data.frame(grupo = pp$grupo, classificacao = nm, classes = length(unique(g[k])), unidades = sum(k), omega2 = x$omega2) })) }))
repres <- bind_rows(data.frame(item = paste("faixa", NOME_FX), presente = NOME_FX %in% fin$faixa),
  data.frame(item = paste("grupo da EBU", sort(unique(C$grupo_ebu))), presente = sort(unique(C$grupo_ebu)) %in% unlist(strsplit(fin$grupos_ebu, ""))),
  data.frame(item = "regiões da EBU", presente = all(unique(C$regiao) %in% unlist(strsplit(fin$regioes, " \\+ ")))))

# ---------------- outputs and map ----------------
write.csv(fin, file.path(raiz, "saida/mapa_unido_bentos_classes.csv"), row.names = FALSE)
write.csv(consenso, file.path(raiz, "saida/mapa_unido_bentos_pares.csv"), row.names = FALSE)
write.csv(aval, file.path(raiz, "saida/mapa_unido_bentos_avaliacao.csv"), row.names = FALSE)
write.csv(repres, file.path(raiz, "saida/mapa_unido_bentos_representacao.csv"), row.names = FALSE)
write.csv(absorvidas, file.path(raiz, "saida/mapa_unido_bentos_absorvidas.csv"), row.names = FALSE)
# color: HB hue (level III) with light/dark variations per class; hatched = stand-by; legend per HB
corHB <- setNames(lh$cor, lh$habitat); hb1 <- sub(" \\+.*$", "", fin$hb)
fin$cor <- unlist(unname(lapply(split(seq_len(nrow(fin)), hb1), function(i) { n <- length(i); f <- if (n == 1) 0 else seq(-0.45, 0.45, length.out = n)
  setNames(vapply(f, function(x) { c0 <- col2rgb(corHB[hb1[i[1]]]) / 255; c1 <- if (x < 0) c0 * (1 + x) else c0 + (1 - c0) * x; rgb(c1[1], c1[2], c1[3]) }, ""), fin$codigo[i]) })))[fin$codigo]
r <- rast(g25); v <- rep(NA, ncell(r)); v[C$cel] <- match(C$codigo, fin$codigo); r[] <- v
writeRaster(r, file.path(raiz, "saida/mapa_unido_bentos_25km.tif"), overwrite = TRUE, datatype = "INT2U")
write.csv(fin |> mutate(valor = seq_len(n())), file.path(raiz, "saida/mapa_unido_bentos_classes.csv"), row.names = FALSE)
arq <- file.path(raiz, "saida/mapa_unido_bentos_25km.png"); hs <- lh$habitat[lh$habitat %in% hb1]; hs <- hs[order(as.integer(sub("HB", "", hs)))]
vv <- values(r)[, 1]; rc <- rowColFromCell(r, seq_len(ncell(r))); k <- which(!is.na(vv)); k <- k[fin$standby[vv[k]] & (rc[k, 1] + rc[k, 2]) %% 7 < 2]
h <- rast(r); hv <- rep(NA, ncell(r)); hv[k] <- 1; h[] <- hv
png(arq, 2000, 1160, res = 150, bg = "white"); layout(matrix(1:2), heights = c(1000, 160)); par(mar = c(0.5, 1, 2, 1))
plot(r, col = fin$cor, breaks = c(0.5, seq_len(nrow(fin)) + 0.5), legend = FALSE, axes = FALSE, box = FALSE, maxcell = ncell(r),
     main = sprintf("Mapa unido do bentos (top-down × bottom-up): %d classes", nrow(fin)))
if (length(k)) plot(h, col = COR_ACH, breaks = c(0.5, 1.5), legend = FALSE, axes = FALSE, box = FALSE, add = TRUE, maxcell = ncell(r))
par(mar = c(0, 1, 0, 1)); plot.new()
legend("top", legend = c(hs, "achurado = não testada"), fill = c(corHB[hs], COR_ACH), density = c(rep(NA, length(hs)), 40), ncol = 10, cex = 0.75, bty = "n",
       title = "matiz = habitat novo (nível III); tons = regiões da EBU dentro dele (nível IV)")
dev.off()
options(width = 200); print(aval |> mutate(omega2 = round(omega2, 3))); print(repres); print(table(fin$faixa, fin$evidencia)); cat(format(Sys.time()), "FIM\n")
