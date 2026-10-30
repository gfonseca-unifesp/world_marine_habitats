# 89b_habitats_novos.R  (Phase 8, stage 3.5: new habitats from the gradients, independent of EBU and EMU)
# Gate 7 (30/10/2026). Steps, separate for benthos (static) and pelagic (dynamic: month x depth):
#  1 GF transformation: each environmental predictor becomes the summed cumulative importance of the compartment's taxonomic groups
#    (curves from 88_gf_gdm.R); basin is excluded (environmental, not geographic, habitats);
#  2 grid: bentos = seafloor cells at 25 km (annual predictors); pelagico = cell x month x standard depth
#    (0, 50, 100, 200, 500, 1000, 2000 m; only where there is water at that depth);
#  3 initial clustering: k-means with K0 = 30 in the transformed space (sample for fitting; all cells assigned
#    to the nearest centroid);
#  4 fauna-guided merging: stage 3 units (X = 10 events, sqrt(beta_sim), 2,500 km blocks) receive the group of the
#    nearest centroid; pairs tested per taxonomic group (pares_classes: omega2 >= 0.01 and power >= 0.8 with equal
#    effort); consensus across groups (equivalent: not distinct in >= 2/3 of >= 2 groups; distinct: >= 2/3) and merging by
#    complete linkage with a 1/3 cut (same rules as 84_equivalencias.R). Group with no testable pair = stand-by;
#  5 threshold rules: classification tree (rpart) of the final habitats on the original predictors (fidelity);
#  6 validation: partial omega2 within the blocks of the new habitats x EBU (group, region) and EMU, per taxonomic group;
#    agreement with EBU and EMU (adjusted Rand index).
# Gate 8 (30/10/2026):
#  - benthos clustered by depth band (0-200, 200-2,000, > 2,000 m; GLORYS prof), 10 groups per band, because
#    k-means weights by area and the abyssal absorbed the tropical shelf; fauna merges within a band and between adjacent bands,
#    never between 0-200 and > 2,000 m; fauna units receive the nearest group within their band;
#  - bottom salinity dropped from benthos (already out of the GF curves);
#  - Euphausiacea and Hyperiidea (no signal in stage 1 or in the Gate 7 validation) excluded from merging and validation;
#  - maps: each stand-by with its own colour and hatched (funcoes_mapa.R).
# Outputs: saida/habitats_novos_{bentos,pelagico}_*.{csv,tif,png}; saida/etapa3_n1_habitats_novos_*.csv
suppressPackageStartupMessages({library(dplyr); library(terra); library(parallel); library(rpart); library(arrow)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); source(file.path(raiz, "scripts/funcoes_grade.R")); source(file.path(raiz, "scripts/funcoes_permanova.R"))
source(file.path(raiz, "scripts/funcoes_mapa.R"))
dir_t <- file.path(raiz, "intermediario/etapa3_n1"); K0 <- 30; ZP <- c(0, 50, 100, 200, 500, 1000, 2000)
KF <- 10; FX <- rep(1:3, each = KF); NOME_FX <- c("0–200 m", "200–2.000 m", "> 2.000 m")   # benthos group k belongs to band FX[k]
SEM_SINAL <- c("Euphausiacea", "Hyperiidea")
faixa <- function(prof) findInterval(prof, c(200, 2000), left.open = TRUE) + 1
ordena <- function(h) h[order(as.integer(sub("^H[BP]", "", h)))]
g25 <- ease2_grade("25km"); M <- lapply(1:12, function(m) rast(file.path(raiz, sprintf("camadas/glorys/mes/glorys_%02d_25km.tif", m))))
curvas <- read.csv(file.path(raiz, "saida/etapa3_n1_gf_curvas.csv"))
transf <- function(E, cp) {   # transformed matrix (sum of the curves of the compartment's groups)
  cv <- curvas[curvas$comp == cp, ]; vs <- intersect(unique(cv$preditor), names(E))
  sapply(vs, function(v) { s <- 0; for (g in unique(cv$grupo[cv$preditor == v])) { c1 <- cv[cv$grupo == g & cv$preditor == v, ]
    s <- s + approx(c1$x, c1$y, E[[v]], rule = 2, ties = mean)$y }; s })
}
mais_proximo <- function(T, C) { d <- sapply(seq_len(nrow(C)), function(k) colSums((t(T) - C[k, ])^2)); max.col(-d, ties.method = "first") }
atribui_b <- function(T, prof, C) { f <- faixa(prof); k <- integer(nrow(T))   # nearest centroid within the band
  for (j in 1:3) { i <- which(f == j); if (length(i)) k[i] <- (j - 1) * KF + mais_proximo(T[i, , drop = FALSE], C[FX == j, , drop = FALSE]) }; k }
ari <- function(a, b) { ok <- !is.na(a) & !is.na(b); t <- table(a[ok], b[ok]); n <- sum(t); s <- function(x) sum(choose(x, 2))
  e <- s(rowSums(t)) * s(colSums(t)) / choose(n, 2); (s(t) - e) / ((s(rowSums(t)) + s(colSums(t))) / 2 - e) }
vm <- function(camada) sapply(1:12, function(k) values(M[[k]][[camada]])[, 1])   # cells x 12 (all cells)

# ---------------- grids ----------------
cat(format(Sys.time()), "grade do bentos\n")
e25 <- read_parquet(file.path(raiz, "camadas/ebu/ebu_celula_25km.parquet"))
p <- do.call(rbind, strsplit(e25$cell_id, "_")); cel_b <- cellFromRowCol(g25, as.integer(p[, 2]), as.integer(p[, 3]))
tf <- vm("tfundo")[cel_b, ]; vf <- vm("vfundo")[cel_b, ]; ss <- vm("thetao_0")[cel_b, ]; ml <- vm("mld")[cel_b, ]; ge <- vm("gelo")[cel_b, ]
GB <- data.frame(cel = cel_b, prof = values(M[[1]][["prof_glorys"]])[cel_b, 1], tfundo_med = rowMeans(tf), tfundo_min = apply(tf, 1, min), tfundo_max = apply(tf, 1, max),
  vfundo_med = rowMeans(vf), vfundo_max = apply(vf, 1, max), sst_med = rowMeans(ss), sst_amp = apply(ss, 1, max) - apply(ss, 1, min),
  mld_max = apply(ml, 1, max), gelo_meses = rowSums(ge >= 0.15, na.rm = TRUE), ebu_grupo = e25$ebu_grupo, ebu_regiao = e25$ebu_regiao)
GB <- GB[complete.cases(GB[, 2:11]), ]; GB$faixa <- faixa(GB$prof); rm(tf, vf, ss, ml, ge)
cat(format(Sys.time()), "grade do pelagico\n")
pg <- values(M[[1]][["prof_glorys"]])[, 1]; emu <- rast(file.path(raiz, "camadas/emu/emu_25km.tif"))
GP <- list()
for (m in 1:12) for (z in ZP) {
  v <- values(M[[m]][[c(paste0("thetao_", z), paste0("so_", z), "thetao_0", "so_0", "mld", "vsup", "gelo")]])
  ok <- which(!is.na(v[, 1]) & !is.na(pg) & pg >= z)
  bz <- which.min(abs(as.numeric(gsub("prof_|m", "", names(emu))) - z))
  GP[[length(GP) + 1]] <- data.frame(cel = ok, mes = m, z = z, t_z = v[ok, 1], s_z = v[ok, 2], sst = v[ok, 3], sss = v[ok, 4], mld = v[ok, 5],
    prof = pg[ok], vsup = v[ok, 6], gelo = ifelse(is.na(v[ok, 7]), 0, v[ok, 7]), emu = values(emu[[bz]])[ok, 1])
}
GP <- bind_rows(GP); GP$z_mld <- log2((GP$z + 1) / (GP$mld + 1)); GP <- GP[complete.cases(GP[, c("t_z", "s_z", "sst", "sss", "mld", "vsup")]), ]
cat("celulas bentos:", nrow(GB), "(por faixa:", paste(NOME_FX, table(GB$faixa), collapse = "; "), ") | linhas pelagico (celula x mes x z):", nrow(GP), "\n")

# ---------------- initial clustering ----------------
TB <- transf(GB, "bentos"); TP <- transf(GP, "pelagico")
kmB <- lapply(1:3, function(f) { i <- which(GB$faixa == f); set.seed(10 + f); s <- if (length(i) > 150000) sample(i, 150000) else i
  kmeans(TB[s, ], KF, nstart = 5, iter.max = 100, algorithm = "MacQueen") })
CB <- do.call(rbind, lapply(kmB, `[[`, "centers"))
set.seed(12); ip <- sample(nrow(TP), min(300000, nrow(TP))); kmP <- kmeans(TP[ip, ], K0, nstart = 5, iter.max = 100, algorithm = "MacQueen")
GB$k0 <- atribui_b(TB, GB$prof, CB); GP$k0 <- unlist(lapply(split(seq_len(nrow(TP)), ceiling(seq_len(nrow(TP)) / 2e5)), function(i) mais_proximo(TP[i, , drop = FALSE], kmP$centers)))
cat(format(Sys.time()), "k-means feitos\n")

# ---------------- fauna-guided merging ----------------
arqs <- list.files(dir_t, pattern = "^(bentos|pelagico)_.*\\.rds$", full.names = TRUE)
arqs <- arqs[!grepl(paste(SEM_SINAL, collapse = "|"), arqs)]
teste_grupo <- function(a, C) {
  suppressPackageStartupMessages({library(vegan); library(dplyr)}); source(file.path(raiz, "scripts/funcoes_permanova.R"))
  r <- readRDS(a); T <- transf(r$E, r$comp); k <- if (r$comp == "bentos") atribui_b(T, r$E$prof, C$bentos) else mais_proximo(T, C$pelagico)
  M2 <- as.matrix(sqrt(betadiver(r$X, "sim"))); set.seed(7000)
  pr <- pares_classes(M2, as.character(k), r$meta$bloco, un_min = 5, n_classe = 10, B = 100)
  list(comp = r$comp, grupo = r$grupo, k = k, pares = if (!is.null(pr)) cbind(comp = r$comp, grupo = r$grupo, pr))
}
cl <- makeCluster(12); clusterExport(cl, c("raiz", "transf", "curvas", "mais_proximo", "atribui_b", "faixa", "KF", "FX"))
TG <- parLapplyLB(cl, arqs, teste_grupo, C = list(bentos = CB, pelagico = kmP$centers)); stopCluster(cl)
pares <- bind_rows(lapply(TG, `[[`, "pares")) |> filter(!is.na(separa))
cat(format(Sys.time()), "pares testados:", nrow(pares), "\n")
consenso <- pares |> group_by(comp, classe_a, classe_b) |> summarise(grupos = n(), distintos = sum(separa), .groups = "drop") |>
  mutate(frac = distintos / grupos, consenso = ifelse(grupos >= 2 & frac <= 1/3, "equivalente", ifelse(frac >= 2/3, "distinto", "depende do grupo")))
funde <- function(cp) {
  e <- consenso |> filter(comp == cp); testados <- sort(unique(c(e$classe_a, e$classe_b)))
  Mx <- matrix(1, K0, K0, dimnames = list(1:K0, 1:K0)); diag(Mx) <- 0
  Mx[cbind(e$classe_a, e$classe_b)] <- ifelse(e$grupos >= 2, e$frac, 1); Mx[cbind(e$classe_b, e$classe_a)] <- Mx[cbind(e$classe_a, e$classe_b)]
  if (cp == "bentos") Mx[abs(outer(FX, FX, "-")) > 1] <- 1                        # 0-200 m never with > 2,000 m
  cl <- cutree(hclust(as.dist(Mx), method = "complete"), h = 1/3)
  data.frame(k0 = 1:K0, final = match(cl, unique(cl)), testado = as.character(1:K0) %in% testados)
}
FB <- funde("bentos"); FP <- funde("pelagico")
nome <- function(F, pref) { F$habitat <- paste0(pref, F$final); sb <- tapply(F$testado, F$final, any); F$standby <- !sb[as.character(F$final)]; F }
FB <- nome(FB, "HB"); FP <- nome(FP, "HP"); FB$faixa <- NOME_FX[FX]; FP$faixa <- NA_character_
GB$habitat <- FB$habitat[GB$k0]; GP$habitat <- FP$habitat[GP$k0]
cat("bentos:", K0, "grupos iniciais ->", length(unique(FB$habitat)), "habitats (", sum(!tapply(FB$testado, FB$final, any)), "stand-by )\n")
cat("pelagico:", K0, "grupos iniciais ->", length(unique(FP$habitat)), "habitats (", sum(!tapply(FP$testado, FP$final, any)), "stand-by )\n")

# ---------------- threshold rules (tree) ----------------
regras <- function(G, vars, F, cp) {
  set.seed(5); i <- sample(nrow(G), min(60000, nrow(G))); d <- G[i, c(vars, "habitat")]; d$habitat <- factor(d$habitat)
  arv <- rpart(habitat ~ ., data = d, method = "class", control = rpart.control(cp = 0.002, maxdepth = 8, minbucket = 50))
  j <- head(setdiff(sample(nrow(G), min(nrow(G), 100000)), i), 20000)            # sample outside the fit
  fid <- mean(as.character(predict(arv, G[j, vars], type = "class")) == G$habitat[j])
  fr <- arv$frame; folhas <- as.integer(rownames(fr)[fr$var == "<leaf>"])
  caminhos <- path.rpart(arv, folhas, print.it = FALSE)
  tab <- data.frame(folha = folhas, habitat = levels(d$habitat)[fr$yval[fr$var == "<leaf>"]], n = fr$n[fr$var == "<leaf>"],
    pureza = apply(fr$yval2[fr$var == "<leaf>", 1 + seq_along(levels(d$habitat)) + length(levels(d$habitat)), drop = FALSE], 1, max),
    regra = vapply(caminhos, function(x) paste(x[-1], collapse = " e "), ""))
  list(arvore = arv, fidelidade = fid, regras = tab |> group_by(habitat) |> arrange(desc(n), .by_group = TRUE) |> ungroup() |> mutate(comp = cp))
}
RB <- regras(GB, c("prof", "tfundo_med", "tfundo_min", "tfundo_max", "vfundo_med", "vfundo_max", "sst_med", "sst_amp", "mld_max", "gelo_meses"), FB, "bentos")
RP <- regras(GP, c("z", "t_z", "s_z", "sst", "sss", "mld", "z_mld", "prof", "vsup", "gelo"), FP, "pelagico")
cat("fidelidade das regras: bentos", round(RB$fidelidade, 3), "| pelagico", round(RP$fidelidade, 3), "\n")

# ---------------- validation: omega2 of the new habitats x EBU/EMU ----------------
val <- bind_rows(lapply(seq_along(TG), function(j) {
  r <- readRDS(arqs[j]); cp <- r$comp; F <- if (cp == "bentos") FB else FP
  hn <- F$habitat[TG[[j]]$k]; Dm <- as.matrix(sqrt(vegan::betadiver(r$X, "sim")))
  ref <- if (cp == "bentos") list("habitats novos" = hn, "EBU grupo" = r$meta$ebu_grupo, "EBU região" = r$meta$ebu_regiao) else list("habitats novos" = hn, "EMU" = as.character(r$meta$emu))
  bind_rows(lapply(names(ref), function(nm) { g <- ref[[nm]]; ok <- !is.na(g)
    if (nm != "habitats novos") ok <- ok & r$meta[[if (nm == "EBU grupo") "pureza_grupo" else if (nm == "EBU região") "pureza_regiao" else "pureza_emu"]] >= 0.8
    gg <- g[ok]; bb <- r$meta$bloco[ok]; keep <- gg %in% names(which(table(gg) >= 5)); gg <- gg[keep]; bb <- bb[keep]; idx <- which(ok)[keep]
    if (length(unique(gg)) < 2) return(NULL)
    set.seed(9); x <- permanova_rapida(as.dist(Dm[idx, idx]), gg, bb, 199)
    data.frame(comp = cp, grupo = r$grupo, classificacao = nm, classes = length(unique(gg)), unidades = length(gg), R2 = x$R2, omega2 = x$omega2, p = x$p) }))
}))
conc <- data.frame(comp = c("bentos", "bentos", "pelagico", "pelagico"), referencia = c("EBU grupo", "EBU região", "EMU (superfície, jan)", "EMU (200 m, jul)"),
  ari = c(ari(GB$habitat, GB$ebu_grupo), ari(GB$habitat, GB$ebu_regiao),
          with(GP[GP$mes == 1 & GP$z == 0, ], ari(habitat, emu)), with(GP[GP$mes == 7 & GP$z == 200, ], ari(habitat, emu))))

# ---------------- profiles, maps and outputs ----------------
perfil <- function(G, vars, F, cp) G |> group_by(habitat) |> summarise(celulas = n(), across(all_of(vars), median), .groups = "drop") |>
  left_join(F |> group_by(habitat) |> summarise(grupos_iniciais = paste(k0, collapse = "+"), standby = first(standby),
    faixa = paste(unique(na.omit(faixa)), collapse = " + "), .groups = "drop"), by = "habitat") |> mutate(comp = cp)
PB <- perfil(GB, c("prof", "tfundo_med", "vfundo_med", "sst_med", "sst_amp", "mld_max", "gelo_meses"), FB, "bentos")
PP <- perfil(GP, c("z", "t_z", "s_z", "sst", "mld", "z_mld", "prof", "gelo"), FP, "pelagico")
write.csv(bind_rows(PB, PP), file.path(raiz, "saida/etapa3_n1_habitats_novos_perfil.csv"), row.names = FALSE)
write.csv(bind_rows(RB$regras, RP$regras), file.path(raiz, "saida/etapa3_n1_habitats_novos_regras.csv"), row.names = FALSE)
write.csv(rbind(data.frame(comp = "bentos", fidelidade = RB$fidelidade), data.frame(comp = "pelagico", fidelidade = RP$fidelidade)), file.path(raiz, "saida/etapa3_n1_habitats_novos_fidelidade.csv"), row.names = FALSE)
write.csv(val, file.path(raiz, "saida/etapa3_n1_habitats_novos_validacao.csv"), row.names = FALSE)
write.csv(conc, file.path(raiz, "saida/etapa3_n1_habitats_novos_concordancia.csv"), row.names = FALSE)
write.csv(consenso, file.path(raiz, "saida/etapa3_n1_habitats_novos_pares.csv"), row.names = FALSE)
write.csv(bind_rows(FB |> mutate(comp = "bentos"), FP |> mutate(comp = "pelagico")), file.path(raiz, "saida/etapa3_n1_habitats_novos_fusoes.csv"), row.names = FALSE)
saveRDS(list(kmB = kmB, CB = CB, kmP = kmP, FB = FB, FP = FP, arvB = RB$arvore, arvP = RP$arvore), file.path(dir_t, "habitats_novos_modelo.rds"))

mapa <- function(G, F, titulo, arq) {
  hs <- ordena(unique(F$habitat)); sb <- hs %in% F$habitat[F$standby]; cores <- unname(PAL_HAB[seq_along(hs)])
  r <- rast(g25); v <- rep(NA, ncell(r)); v[G$cel] <- match(G$habitat, hs); r[] <- v
  writeRaster(r, file.path(raiz, sprintf("saida/habitats_novos_%s.tif", arq)), overwrite = TRUE, datatype = "INT1U")
  mapa_classes(r, cores, sb, ifelse(sb, paste(hs, "(stand-by)"), hs), titulo, file.path(raiz, sprintf("saida/mapa_habitats_novos_%s.png", arq)))
  data.frame(arquivo = arq, habitat = hs, cor = cores, standby = sb)
}
leg <- bind_rows(mapa(GB, FB, "Habitats bentônicos novos (Gradient Forest + fauna), nível 1", "bentos_25km"),
  mapa(GP[GP$mes == 1 & GP$z == 0, ], FP, "Habitats pelágicos novos: superfície, janeiro", "pelagico_0m_jan"),
  mapa(GP[GP$mes == 7 & GP$z == 0, ], FP, "Habitats pelágicos novos: superfície, julho", "pelagico_0m_jul"),
  mapa(GP[GP$mes == 1 & GP$z == 200, ], FP, "Habitats pelágicos novos: 200 m, janeiro", "pelagico_200m_jan"),
  mapa(GP[GP$mes == 7 & GP$z == 200, ], FP, "Habitats pelágicos novos: 200 m, julho", "pelagico_200m_jul"))
write.csv(leg, file.path(raiz, "saida/habitats_novos_legenda.csv"), row.names = FALSE)
# all pelagic layers (month x z) for stage 4
hsP <- ordena(unique(FP$habitat)); camadas <- list()
for (m in 1:12) for (z in ZP) { r <- rast(g25); v <- rep(NA, ncell(r)); s <- GP[GP$mes == m & GP$z == z, ]
  v[s$cel] <- match(s$habitat, hsP); r[] <- v; names(r) <- sprintf("m%02d_z%d", m, z); camadas[[length(camadas) + 1]] <- r }
st <- rast(camadas)
writeRaster(st, file.path(raiz, "saida/habitats_novos_pelagico_mes_z.tif"), overwrite = TRUE, datatype = "INT1U")
options(width = 200); print(val |> mutate(across(where(is.numeric), ~ round(.x, 3)))); print(conc)
cat(format(Sys.time()), "FIM\n")
