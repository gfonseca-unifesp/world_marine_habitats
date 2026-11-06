# 91_mapa_unido_bentos.R  (Phase 8; Gates 9 and 9b, 30/10/2026: final benthic map uniting top-down and bottom-up)
# Principle: separate when any evidence separates; merge when the fauna shows equivalence, or, for untested pairs, when the
# environment is as similar as in fauna-equivalent pairs, with a measured false-merge rate.
#  1 candidates: each cell = depth band x new habitat (HB, 89b) x EBU region; protected features (EBU groups C seamounts and
#    D ridges/rifts) = depth band x EBU region, not crossed with HB (Gate 9b item 3);
#  2 merges on a single tree (complete linkage) whose dissimilarity is on a "false-merge rate" scale (Gate 9b):
#      fauna-equivalent pair (consensus: not distinct in >= 2/3 of >= 2 taxonomic groups) = 0;
#      distinct or "depende do grupo" pair, different depth band, or protected feature with another group = never (2);
#      untested pair (or tested in a single group without separating) = share of fauna-distinct pairs whose Gradient Forest
#      distance is <= the pair's distance (the false-merge rate that merging it implies);
#    base level = cut at the rate of the 25th percentile of the GF distances of equivalent pairs (about 4%); tolerance levels =
#    cuts at the 50th, 75th and 90th percentiles; levels are nested by construction (same tree);
#  3 protected features only merge within the same EBU group;
#  4 hierarchy: I compartment, II depth band, III HB (or protected feature group), IV base classes, plus tolerance levels;
#  5 minimum area at the base level: class with < 10 cells, not protected and without a fauna-distinct pair, goes to the nearest
#    neighboring class (4-cell adjacency) of the same depth band in GF space; repeat until stable;
#  6 attributes: evidence (distinct by fauna / merged by fauna / merged by environment / tested without a conclusive pair /
#    untested = hatched), evidence inherited from level III (HB validated by fauna or not), protected, shelf seasonality;
#  7 acceptance: omega2 of the base classes >= that of the best component (HB, EBU region) on the same units, per group, with
#    all data and with equal effort (Gate 9b item 4); all depth bands and EBU groups represented;
#  8 the pelagic goes in a separate layer (91b_mapa_unido_pelagico.R);
#  9 external evidence (Gate 10, L3): splits decided in 92c_redefinicao.R (saida/limiares_candidatos_redefinicao.csv) are applied
#    inside the study box; new classes carry evidencia_externa, referencia and classe_mae, and are cross-hatched on the map.
# Outputs: saida/mapa_unido_bentos_{classes,pares,avaliacao,representacao,absorvidas,hierarquia,niveis}.csv,
#          saida/mapa_unido_bentos_25km.tif (bands: base, nivel3, tolerance levels), saida/mapa_unido_bentos_{25km,nivel3}.png
suppressPackageStartupMessages({library(dplyr); library(terra); library(arrow); library(parallel)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); source(file.path(raiz, "scripts/funcoes_grade.R")); source(file.path(raiz, "scripts/funcoes_permanova.R"))
source(file.path(raiz, "scripts/funcoes_mapa.R"))
dir_t <- file.path(raiz, "intermediario/etapa3_n1"); NOME_FX <- c("0–200 m", "200–2.000 m", "> 2.000 m"); PROT <- c("C", "D"); AREA_MIN <- 10
NOME_PROT <- c(C = "montes submarinos (C)", D = "dorsais e riftes (D)"); NUNCA <- 2; Q_BASE <- 0.25; Q_NIVEIS <- c(0.50, 0.75, 0.90)
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
C$faixa <- faixa(C$prof); C$prot <- C$grupo_ebu %in% PROT
C$nivel3 <- ifelse(C$prot, paste(C$faixa, C$grupo_ebu, sep = "|"), paste(C$faixa, C$hb, sep = "|"))
C$cand <- ifelse(C$prot, paste(C$faixa, "FP", C$regiao, sep = "|"), paste(C$faixa, C$hb, C$regiao, sep = "|"))
TC <- transf(C, "bentos")
cand <- C |> group_by(cand) |> summarise(faixa = first(faixa), hb = names(which.max(table(hb))), regiao = first(regiao), grupo_ebu = first(grupo_ebu),
  prot = first(prot), nivel3 = first(nivel3), celulas = n(), .groups = "drop")
CT <- rowsum(TC, C$cand)[cand$cand, , drop = FALSE] / cand$celulas
cat(format(Sys.time()), "celulas:", nrow(C), "| candidatos:", nrow(cand), "(protegidos:", sum(cand$prot), ") | com <", AREA_MIN, "celulas:", sum(cand$celulas < AREA_MIN), "\n")
pode <- function(a, b) { x <- cand[match(a, cand$cand), ]; y <- cand[match(b, cand$cand), ]   # rules 3 and 4
  x$faixa == y$faixa & (!(x$prot | y$prot) | x$grupo_ebu == y$grupo_ebu) }

# ---------------- 2 fauna tests ----------------
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
  mutate(frac = distintos / grupos, consenso = ifelse(grupos >= 2 & frac <= 1/3, "equivalente", ifelse(frac >= 2/3, "distinto", "depende do grupo")),
    fraco = grupos == 1 & distintos == 0)   # tested in a single group without separating: treated as untested

# ---------------- 2 tolerance tree ----------------
n <- nrow(cand); D <- as.matrix(dist(CT)); dimnames(D) <- list(cand$cand, cand$cand)
ia <- match(consenso$classe_a, cand$cand); ib <- match(consenso$classe_b, cand$cand)
dq <- D[cbind(ia, ib)][consenso$consenso == "equivalente"]; dd <- D[cbind(ia, ib)][consenso$consenso == "distinto"]
taxa <- ecdf(dd)                                                     # false-merge rate implied by a GF distance
Tx <- matrix(taxa(D), n, n, dimnames = dimnames(D))                  # untested pairs
v <- ifelse(consenso$fraco, taxa(D[cbind(ia, ib)]), ifelse(consenso$consenso == "equivalente", 0, NUNCA))
Tx[cbind(ia, ib)] <- v; Tx[cbind(ib, ia)] <- v
Tx[!outer(seq_len(n), seq_len(n), function(i, j) pode(cand$cand[i], cand$cand[j]))] <- NUNCA; diag(Tx) <- 0
h_base <- taxa(quantile(dq, Q_BASE)); h_niv <- taxa(quantile(dq, Q_NIVEIS))
cat(sprintf("distancia GF: equivalentes mediana %.3f, distintos mediana %.3f | erro de fusao: base %.3f; niveis %s\n",
  median(dq), median(dd), h_base, paste(sprintf("%.3f", h_niv), collapse = ", ")))
arvore <- hclust(as.dist(Tx), method = "complete")
rep_ <- tapply(cand$cand, cutree(arvore, h = h_base), function(x) x[1]); cand$classe <- unname(rep_[as.character(cutree(arvore, h = h_base))])
eqp <- consenso |> filter(consenso == "equivalente")
cand$fundido_fauna <- vapply(seq_len(n), function(i) { mem <- cand$cand[cand$classe == cand$classe[i]]
  length(mem) > 1 && any(eqp$classe_a %in% mem & eqp$classe_b %in% mem) }, TRUE)
dist_c <- unique(c(consenso$classe_a[consenso$consenso == "distinto"], consenso$classe_b[consenso$consenso == "distinto"]))
testado <- unique(c(consenso$classe_a[!consenso$fraco], consenso$classe_b[!consenso$fraco]))

# ---------------- 5 minimum area (base level) ----------------
C$classe <- cand$classe[match(C$cand, cand$cand)]
cls_info <- function() C |> group_by(classe) |> summarise(celulas = n(), faixa = first(faixa), prot = any(prot), .groups = "drop")
centro <- function() { s <- rowsum(TC, C$classe); s / as.vector(table(C$classe)[rownames(s)]) }
distinta <- function(cl) any(cand$cand[cand$classe == cl] %in% dist_c)
absorvidas <- list(); repeat {
  ci <- cls_info(); peq <- ci |> filter(celulas < AREA_MIN, !prot) |> arrange(celulas); peq <- peq[!vapply(peq$classe, distinta, TRUE), ]
  if (!nrow(peq)) break
  CTc <- centro(); mudou <- FALSE; r <- rast(g25); vv <- rep(NA_character_, ncell(r)); vv[C$cel] <- C$classe
  for (k in seq_len(nrow(peq))) { a <- peq$classe[k]; ca <- C$cel[C$classe == a]; if (!length(ca)) next
    viz <- unique(vv[as.vector(adjacent(r, ca, "rook"))]); viz <- setdiff(viz[!is.na(viz)], a)
    viz <- viz[ci$faixa[match(viz, ci$classe)] == peq$faixa[k] & !ci$prot[match(viz, ci$classe)]]
    if (!length(viz)) next
    d <- colSums((t(CTc[viz, , drop = FALSE]) - CTc[a, ])^2); dest <- viz[which.min(d)]
    C$classe[C$classe == a] <- dest; vv[ca] <- dest; cand$classe[cand$classe == a] <- dest; absorvidas[[length(absorvidas) + 1]] <- data.frame(de = a, para = dest); mudou <- TRUE }
  if (!mudou) break
}
absorvidas <- bind_rows(absorvidas); cat(format(Sys.time()), "classes absorvidas por area:", nrow(absorvidas), "\n")

# ---------------- 6 attributes ----------------
hb_valid <- setNames(!lh$standby, lh$habitat)
fin <- C |> group_by(classe) |> summarise(celulas = n(), faixa = NOME_FX[first(faixa)], nivel3 = names(which.max(table(nivel3))),
  hb = paste(sort(unique(hb)), collapse = " + "), regioes = paste(sort(unique(regiao)), collapse = " + "),
  grupos_ebu = paste(sort(unique(grupo_ebu)), collapse = ""), protegida = any(prot), tfundo_amp = median(tfundo_amp), prof = median(prof), .groups = "drop") |>
  left_join(cand |> group_by(classe) |> summarise(candidatos = n(), cand_fauna = sum(cand %in% testado), fundida_fauna = any(fundido_fauna),
    distinta_fauna = any(cand %in% dist_c), .groups = "drop"), by = "classe") |>
  mutate(evidencia = case_when(distinta_fauna ~ "distinta pela fauna", fundida_fauna ~ "fundida pela fauna",
      candidatos > 1 ~ "fundida pelo ambiente", cand_fauna > 0 ~ "testada, sem par conclusivo", TRUE ~ "não testada"),
    standby = evidencia %in% c("não testada", "testada, sem par conclusivo", "fundida pelo ambiente"),
    nivel3_rotulo = ifelse(protegida, paste(NOME_FX[as.integer(sub("\\|.*", "", nivel3))], NOME_PROT[sub(".*\\|", "", nivel3)], sep = " · "),
                           paste(NOME_FX[as.integer(sub("\\|.*", "", nivel3))], sub(".*\\|", "", nivel3), sep = " · ")),
    evidencia_nivel3 = ifelse(protegida, "feição protegida (geomorfológica)", ifelse(hb_valid[sub(".*\\|", "", nivel3)] %in% TRUE, "HB validado pela fauna", "HB em stand-by")),
    sazonalidade = ifelse(faixa == "0–200 m", ifelse(tfundo_amp >= 5, "forte", ifelse(tfundo_amp >= 2, "moderada", "fraca")), "–")) |>
  arrange(match(faixa, NOME_FX), nivel3, regioes)
fin$codigo <- paste0("HF", seq_len(nrow(fin))); C$codigo <- fin$codigo[match(C$classe, fin$classe)]
cat("classes base:", nrow(fin), "| por evidencia:", paste(names(table(fin$evidencia)), table(fin$evidencia), collapse = "; "), "\n")

# ---------------- 2 tolerance levels above the base (same tree, nested) ----------------
# class-level dissimilarity = maximum over member candidates (complete linkage), so levels stay nested above the base classes
mem <- split(match(cand$cand, rownames(Tx)), fin$codigo[match(cand$classe, fin$classe)])[fin$codigo]
nb <- length(mem); Tc <- matrix(0, nb, nb, dimnames = list(fin$codigo, fin$codigo))
for (i in seq_len(nb - 1)) for (j in (i + 1):nb) Tc[i, j] <- Tc[j, i] <- max(Tx[mem[[i]], mem[[j]]])
arv_c <- hclust(as.dist(Tc), method = "complete")
niveis <- data.frame(codigo = fin$codigo)
for (k in seq_along(h_niv)) { g <- cutree(arv_c, h = h_niv[k]); niveis[[sprintf("tol_%02.0f", 100 * h_niv[k])]] <- unname(tapply(fin$codigo, g, function(x) x[which.max(fin$celulas[match(x, fin$codigo)])])[as.character(g)]) }
hier <- data.frame(passo = seq_len(nrow(arv_c$merge)), a = arv_c$merge[, 1], b = arv_c$merge[, 2], altura = arv_c$height) |> filter(altura < NUNCA)
cat("classes por nivel de tolerancia:", paste(names(niveis)[-1], sapply(niveis[-1], function(x) length(unique(x))), collapse = "; "), "\n")

# ---------------- 9 external evidence (Gate 10; L3, 92c_redefinicao.R) ----------------
# inside the study box, the part of a base class of the tested level-III context that lies beyond the published value becomes a
# new class "<parent>-<id>" (if it has >= AREA_MIN cells and the parent keeps cells on the other side); it never merges with the
# parent at any tolerance level (its tolerance labels get the same suffix); evidence "distinta pela fauna" or "separada pela literatura"
fin$evidencia_externa <- NA_character_; fin$referencia <- NA_character_; fin$classe_mae <- NA_character_
niveis <- niveis[match(fin$codigo, niveis$codigo), ]
fr <- file.path(raiz, "saida/limiares_candidatos_redefinicao.csv")
red <- if (file.exists(fr)) read.csv(fr, encoding = "UTF-8") |> filter(decisao == "partir") else NULL
if (!is.null(red) && nrow(red)) {
  ll <- project(xyFromCell(g25, C$cel), from = crs(g25), to = "EPSG:4326")
  for (k in seq_len(nrow(red))) { rr <- red[k, ]; e <- ext(vect(rr$caixa_wkt, crs = "EPSG:4326"))
    dentro <- C$nivel3 == rr$nivel3 & ll[, 1] >= e$xmin & ll[, 1] <= e$xmax & ll[, 2] >= e$ymin & ll[, 2] <= e$ymax
    alem <- C[[rr$variavel]] >= rr$valor
    for (pai in unique(C$codigo[dentro & alem])) { i <- which(dentro & alem & C$codigo == pai)
      if (length(i) < AREA_MIN || !any(dentro & !alem & C$codigo == pai)) next
      novo <- paste(pai, rr$id, sep = "-"); C$codigo[i] <- novo; lin <- fin[fin$codigo == pai, ]
      lin$codigo <- novo; lin$classe_mae <- pai; lin$evidencia_externa <- rr$evidencia; lin$referencia <- rr$referencia
      lin$evidencia <- rr$evidencia; lin$standby <- FALSE; fin <- bind_rows(fin, lin)
      nl <- niveis[niveis$codigo == pai, ]; nl$codigo <- novo; nl[-1] <- lapply(nl[-1], function(x) paste(x, rr$id, sep = "-")); niveis <- bind_rows(niveis, nl) } }
  aj <- C |> group_by(codigo) |> summarise(celulas = n(), prof = median(prof), tfundo_amp = median(tfundo_amp), .groups = "drop")
  fin <- fin |> select(-celulas, -prof, -tfundo_amp) |> left_join(aj, by = "codigo") |> filter(!is.na(celulas))
  niveis <- niveis[match(fin$codigo, niveis$codigo), ]
  for (nm in names(niveis)[-1]) { lab <- niveis[[nm]]; falta <- !lab %in% niveis$codigo   # suffixed labels -> largest child carrying them
    if (any(falta)) { rp <- tapply(niveis$codigo[falta], lab[falta], function(x) x[which.max(fin$celulas[match(x, fin$codigo)])]); niveis[[nm]][falta] <- unname(rp[lab[falta]]) } }
  cat("criterio 9: classes novas pela literatura:", sum(!is.na(fin$classe_mae)), "(", paste(names(table(fin$evidencia_externa)), table(fin$evidencia_externa), collapse = "; "), ")\n")
}

# ---------------- 7 acceptance (all data and equal effort) ----------------
avalia <- function(pp, chave, cls) {   # cls: classifications aligned with chave (candidate codes)
  suppressPackageStartupMessages({library(vegan); library(dplyr)}); source(file.path(raiz, "scripts/funcoes_permanova.R"))
  r <- readRDS(pp$arq); ok <- !is.na(pp$cand); Dm <- as.matrix(sqrt(betadiver(r$X[ok, , drop = FALSE], "sim"))); bl <- r$meta$bloco[ok]
  bind_rows(lapply(names(cls), function(nm) { g <- cls[[nm]][match(pp$cand[ok], chave)]; k <- g %in% names(which(table(g) >= 5))
    if (length(unique(g[k])) < 2) return(NULL)
    set.seed(9); x <- permanova_rapida(as.dist(Dm[k, k]), g[k], bl[k], 199)
    set.seed(10); e <- tryCatch(permanova_esforco_igual(Dm[k, k], g[k], bl[k], n_classe = 10, B = 100), error = function(e) NULL)
    data.frame(grupo = pp$grupo, classificacao = nm, classes = length(unique(g[k])), unidades = sum(k), omega2 = x$omega2,
      omega2_eq = if (!is.null(e)) e$omega2_eq else NA, poder_eq = if (!is.null(e)) e$poder_eq else NA) })) }
cls_av <- list("mapa unido (base)" = fin$codigo[match(cand$classe, fin$classe)], "nível III" = cand$nivel3,
               "habitats novos (HB)" = cand$hb, "EBU região" = cand$regiao)
cl <- makeCluster(7); clusterExport(cl, c("raiz"))
aval <- bind_rows(parLapply(cl, prep, avalia, chave = cand$cand, cls = cls_av)); stopCluster(cl)
repres <- bind_rows(data.frame(item = paste("faixa", NOME_FX), presente = NOME_FX %in% fin$faixa),
  data.frame(item = paste("grupo da EBU", sort(unique(C$grupo_ebu))), presente = sort(unique(C$grupo_ebu)) %in% unlist(strsplit(fin$grupos_ebu, ""))),
  data.frame(item = "regiões da EBU", presente = all(unique(C$regiao) %in% unlist(strsplit(fin$regioes, " \\+ ")))))

# ---------------- colors, outputs and maps ----------------
# color: hue of level III (HB color; protected features in fixed hues) with light/dark variations per base class
corN3 <- c(setNames(lh$cor, lh$habitat), C = "#5b4a8b", D = "#8b5a2b"); chave3 <- sub(".*\\|", "", fin$nivel3)
fin$cor <- unlist(unname(lapply(split(seq_len(nrow(fin)), fin$nivel3), function(i) { m <- length(i); f <- if (m == 1) 0 else seq(-0.45, 0.45, length.out = m)
  setNames(vapply(f, function(x) { c0 <- col2rgb(corN3[chave3[i[1]]]) / 255; c1 <- if (x < 0) c0 * (1 + x) else c0 + (1 - c0) * x; rgb(c1[1], c1[2], c1[3]) }, ""), fin$codigo[i]) })))[fin$codigo]
fin <- fin |> mutate(valor = seq_len(n())) |> left_join(niveis, by = "codigo")
n3 <- sort(unique(fin$nivel3)); fin$valor_nivel3 <- match(fin$nivel3, n3)
write.csv(fin, file.path(raiz, "saida/mapa_unido_bentos_classes.csv"), row.names = FALSE)
write.csv(consenso, file.path(raiz, "saida/mapa_unido_bentos_pares.csv"), row.names = FALSE)
write.csv(aval, file.path(raiz, "saida/mapa_unido_bentos_avaliacao.csv"), row.names = FALSE)
write.csv(repres, file.path(raiz, "saida/mapa_unido_bentos_representacao.csv"), row.names = FALSE)
write.csv(absorvidas, file.path(raiz, "saida/mapa_unido_bentos_absorvidas.csv"), row.names = FALSE)
write.csv(hier, file.path(raiz, "saida/mapa_unido_bentos_hierarquia.csv"), row.names = FALSE)
write.csv(data.frame(nivel = c("base", names(niveis)[-1]), percentil_equivalentes = c(Q_BASE, Q_NIVEIS), erro_fusao = c(h_base, h_niv),
  classes = c(nrow(fin), sapply(niveis[-1], function(x) length(unique(x))))), file.path(raiz, "saida/mapa_unido_bentos_niveis.csv"), row.names = FALSE)
vcel <- fin$valor[match(C$codigo, fin$codigo)]
camada <- function(x) { r <- rast(g25); v <- rep(NA, ncell(r)); v[C$cel] <- x; r[] <- v; r }
st <- c(camada(vcel), camada(fin$valor_nivel3[vcel]), rast(lapply(names(niveis)[-1], function(nm) camada(fin$valor[match(fin[[nm]], fin$codigo)][vcel]))))
names(st) <- c("base", "nivel3", names(niveis)[-1]); writeRaster(st, file.path(raiz, "saida/mapa_unido_bentos_25km.tif"), overwrite = TRUE, datatype = "INT2U")
# base map: too many classes for a legend per class; the legend shows the level-III hues
hbs <- unique(chave3[grepl("^HB", chave3)]); hs <- c(hbs[order(as.integer(sub("HB", "", hbs)))], intersect(c("C", "D"), chave3))
rotulo3 <- ifelse(hs %in% names(NOME_PROT), NOME_PROT[hs], hs)
rb <- st[["base"]]; vb <- values(rb)[, 1]; rc <- rowColFromCell(rb, seq_len(ncell(rb))); k0 <- which(!is.na(vb)); k <- k0[fin$standby[vb[k0]] & (rc[k0, 1] + rc[k0, 2]) %% 7 < 2]
lit <- fin$evidencia == "separada pela literatura"; k <- c(k, k0[lit[vb[k0]] & ((rc[k0, 1] + rc[k0, 2]) %% 7 < 2 | (rc[k0, 1] - rc[k0, 2]) %% 7 < 2)])   # cross-hatched
h <- rast(rb); hv <- rep(NA, ncell(rb)); hv[k] <- 1; h[] <- hv
png(file.path(raiz, "saida/mapa_unido_bentos_25km.png"), 2000, 1160, res = 150, bg = "white"); layout(matrix(1:2), heights = c(1000, 160)); par(mar = c(0.5, 1, 2, 1))
plot(rb, col = fin$cor, breaks = c(0.5, seq_len(nrow(fin)) + 0.5), legend = FALSE, axes = FALSE, box = FALSE, maxcell = ncell(rb),
     main = sprintf("Mapa unido do bentos, nível base (erro de fusão ≤ %.0f%%): %d classes", 100 * h_base, nrow(fin)))
if (length(k)) plot(h, col = COR_ACH, breaks = c(0.5, 1.5), legend = FALSE, axes = FALSE, box = FALSE, add = TRUE, maxcell = ncell(rb))
par(mar = c(0, 1, 0, 1)); plot.new()
legend("top", legend = rotulo3, fill = corN3[hs], ncol = 11, cex = 0.7, bty = "n",
       title = "matiz = nível III (habitat novo ou feição protegida); tons = classes do nível base; achurado = sem evidência da fauna; cruzado = separada pela literatura")
dev.off()
cor3 <- vapply(n3, function(x) corN3[sub(".*\\|", "", x)], ""); val3 <- fin$evidencia_nivel3[match(n3, fin$nivel3)] != "HB validado pela fauna" & !grepl("\\|[CD]$", n3)
mapa_classes(st[["nivel3"]], cor3, val3, fin$nivel3_rotulo[match(n3, fin$nivel3)], "Mapa unido do bentos, nível III (faixa × habitat novo; feições protegidas)",
  file.path(raiz, "saida/mapa_unido_bentos_nivel3.png"), ncol = 4)
options(width = 200); print(aval |> mutate(across(where(is.double), ~ round(.x, 3)))); print(repres); print(table(fin$faixa, fin$evidencia)); cat(format(Sys.time()), "FIM\n")
