# 91b_mapa_unido_pelagico.R  (Phase 8; Gate 9/9b rules applied to the pelagic, after level 2 and Gate 10, 06/11/2026)
# Same principle as 91_mapa_unido_bentos.R, on the dynamic pelagic grid (25 km cell x month x standard depth):
#  1 candidates: each row = vertical zone (epipelagic z < 200 m, mesopelagic 200-1,000 m, bathypelagic >= 1,000 m; the z ~ 200 m
#    threshold of levels 1 and 2 and of the literature) x pelagic habitat (HP of 89b, updated with the level-2 tests of
#    93_pelagico_n2.R: merged HP keep the name of the lowest one) x EMU class at the nearest EMU depth;
#  2 fauna tests: level-1 pelagic units (taxonomic groups; Euphausiacea and Hyperiidea out, Gate 8) and level-2 units (ecological
#    groups) each get the candidate of their own environment (nearest kmP centroid -> HP; zone from z; EMU of the cell at z);
#    pairs per group (pares_classes, omega2 >= 0.01 and power >= 0.8), consensus as in 91;
#  3 single tree (complete linkage) on the false-merge-rate scale: equivalent = 0; distinct, "depende do grupo" or other zone =
#    never (2); untested = ecdf of the GF distances of distinct pairs; base cut at the 25th percentile of equivalent pairs,
#    tolerance levels at the 50th, 75th and 90th (nested);
#  4 minimum area: class present in < 10 cells in every month x depth layer, without a distinct pair, goes to the nearest class
#    in GF space among its 4-cell neighbours in the same layer and zone; repeat until stable;
#  5 attributes: evidence (as in 91), level III = zone x HP (validated by fauna or stand-by), persistence (mean share of the 12
#    months in which a cell x depth keeps the class: resident >= 0.75, seasonal 0.25-0.75, occasional < 0.25);
#  6 acceptance: omega2 of the base classes vs. the components (HP, EMU, level III), with all data and with equal effort.
# Outputs: saida/mapa_unido_pelagico_{classes,pares,avaliacao,absorvidas,hierarquia,niveis}.csv,
#          saida/mapa_unido_pelagico_mes_z.tif (base codes, 84 layers m01_z0 ... m12_z2000), saida/mapa_unido_pelagico_*.png
suppressPackageStartupMessages({library(dplyr); library(terra); library(arrow); library(parallel)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); source(file.path(raiz, "scripts/funcoes_grade.R")); source(file.path(raiz, "scripts/funcoes_permanova.R"))
source(file.path(raiz, "scripts/funcoes_mapa.R"))
ZP <- c(0, 50, 100, 200, 500, 1000, 2000); NOME_Z <- c("epipelágico (< 200 m)", "mesopelágico (200–1.000 m)", "batipelágico (≥ 1.000 m)")
SEM_SINAL <- c("Euphausiacea", "Hyperiidea"); AREA_MIN <- 10; NUNCA <- 2; Q_BASE <- 0.25; Q_NIVEIS <- c(0.50, 0.75, 0.90)
zona <- function(z) findInterval(z, c(200, 1000)) + 1
ordena <- function(h) h[order(as.integer(sub("^H[BP]", "", h)))]
g25 <- ease2_grade("25km"); M <- lapply(1:12, function(m) rast(file.path(raiz, sprintf("camadas/glorys/mes/glorys_%02d_25km.tif", m))))
emu <- rast(file.path(raiz, "camadas/emu/emu_25km.tif")); zemu <- as.numeric(gsub("prof_|m", "", names(emu)))
banda_emu <- function(z) vapply(z, function(x) which.min(abs(zemu - x)), 1L)
curvas_p <- open_dataset(file.path(raiz, "saida/etapa3_n1_gf_curvas.csv"), format = "csv") |> filter(comp == "pelagico") |> collect() |> as.data.frame()
transf <- function(E) { vs <- intersect(unique(curvas_p$preditor), names(E))   # as in 89b (level-1 pelagic curves)
  sapply(vs, function(v) { s <- 0; for (g in unique(curvas_p$grupo[curvas_p$preditor == v])) { c1 <- curvas_p[curvas_p$grupo == g & curvas_p$preditor == v, ]
    s <- s + approx(c1$x, c1$y, E[[v]], rule = 2, ties = mean)$y }; s }) }
mais_proximo <- function(T, C) { d <- sapply(seq_len(nrow(C)), function(k) colSums((t(T) - C[k, ])^2)); max.col(-d, ties.method = "first") }

# ---------------- HP updated with level 2 ----------------
fus <- read.csv(file.path(raiz, "saida/pelagico_n2_fusoes.csv"))
nome_final <- tapply(fus$hp_n1, fus$final_n1n2, function(x) ordena(unique(x))[1])
fus$hp <- unname(nome_final[as.character(fus$final_n1n2)])                       # merged HP keep the lowest name
hp_de_antigo <- tapply(fus$hp, fus$hp_n1, first); hsP_antigo <- ordena(unique(fus$hp_n1))   # value order of habitats_novos_pelagico_mes_z.tif
hp_standby <- tapply(fus$standby_n1n2, fus$hp, first)
cat("HP atualizados:", length(hp_standby), "(stand-by:", sum(hp_standby), ")\n")

# ---------------- 1 candidates (rows = cell x month x depth) ----------------
hpt <- rast(file.path(raiz, "saida/habitats_novos_pelagico_mes_z.tif")); pg <- values(M[[1]][["prof_glorys"]])[, 1]
G <- list(); for (m in 1:12) for (z in ZP) {
  h <- values(hpt[[sprintf("m%02d_z%d", m, z)]])[, 1]; ok <- which(!is.na(h))
  v <- values(M[[m]][[c(paste0("thetao_", z), paste0("so_", z), "thetao_0", "so_0", "mld", "vsup", "gelo")]])[ok, ]
  G[[length(G) + 1]] <- data.frame(cel = ok, mes = m, z = z, t_z = v[, 1], s_z = v[, 2], sst = v[, 3], sss = v[, 4], mld = v[, 5], vsup = v[, 6],
    gelo = ifelse(is.na(v[, 7]), 0, v[, 7]), prof = pg[ok], emu = values(emu[[banda_emu(z)]])[ok, 1], hp = unname(hp_de_antigo[hsP_antigo[h[ok]]])) }
G <- bind_rows(G); G$z_mld <- log2((G$z + 1) / (G$mld + 1)); G$zona <- zona(G$z); G <- G[!is.na(G$emu), ]
G$camada <- (G$mes - 1) * length(ZP) + match(G$z, ZP)
G$nivel3 <- paste(G$zona, G$hp, sep = "|"); G$cand <- paste(G$zona, G$hp, G$emu, sep = "|")
TG <- transf(G); cat(format(Sys.time()), "linhas:", nrow(G), "\n")
cand <- G |> count(cand, name = "linhas") |> mutate(zona = as.integer(sub("\\|.*", "", cand)), hp = sub("^[0-9]\\|(HP[0-9]+)\\|.*", "\\1", cand),
  emu = as.integer(sub(".*\\|", "", cand)), nivel3 = paste(zona, hp, sep = "|"))
CT <- rowsum(TG, G$cand)[cand$cand, , drop = FALSE] / cand$linhas
cat("candidatos:", nrow(cand), "\n")
pode <- function(a, b) cand$zona[match(a, cand$cand)] == cand$zona[match(b, cand$cand)]

# ---------------- 2 fauna tests ----------------
mod <- readRDS(file.path(raiz, "intermediario/etapa3_n1/habitats_novos_modelo.rds")); CP <- mod$kmP$centers; rm(mod); gc(verbose = FALSE)
cel_de <- function(id) { p <- do.call(rbind, strsplit(id, "_")); cellFromRowCol(g25, as.integer(p[, 2]), as.integer(p[, 3])) }
arqs <- c(list.files(file.path(raiz, "intermediario/etapa3_n1"), pattern = "^pelagico_.*\\.rds$", full.names = TRUE),
          list.files(file.path(raiz, "intermediario/pelagico_n2"), pattern = "\\.rds$", full.names = TRUE))
arqs <- arqs[!grepl(paste(SEM_SINAL, collapse = "|"), arqs)]
prep <- lapply(arqs, function(a) { r <- readRDS(a); n2 <- grepl("pelagico_n2", a)
  k0 <- mais_proximo(transf(r$E), CP); hp <- fus$hp[k0]
  em <- if (!n2) r$meta$emu else { ev <- as.matrix(terra::extract(emu, cel_de(r$meta$cell_id))); ev[cbind(seq_len(nrow(ev)), banda_emu(r$E$z))] }
  cc <- paste(zona(r$E$z), hp, em, sep = "|"); cc[!cc %in% cand$cand] <- NA
  n <- table(cc); cl <- names(n)[n >= 5]; pr <- if (length(cl) >= 2) t(combn(sort(cl), 2)) else matrix(character(0), 0, 2)
  pr <- pr[pode(pr[, 1], pr[, 2]), , drop = FALSE]
  list(arq = a, grupo = paste0(r$grupo, if (n2) " (n2)" else ""), cand = cc, hp = hp, nivel3 = paste(zona(r$E$z), hp, sep = "|"), emu = em, pares = pr) })
cat(format(Sys.time()), "pares a testar:", paste(sapply(prep, `[[`, "grupo"), sapply(prep, function(x) nrow(x$pares)), collapse = "; "), "\n")
jobs <- bind_rows(lapply(seq_along(prep), function(k) { n <- nrow(prep[[k]]$pares); if (!n) return(NULL)
  data.frame(k = k, de = seq(1, n, by = 100)) |> mutate(ate = pmin(de + 99, n)) }))
testa <- function(j) {
  suppressPackageStartupMessages({library(vegan); library(dplyr)}); source(file.path(raiz, "scripts/funcoes_permanova.R"))
  pp <- prep[[j$k]]; r <- readRDS(pp$arq); ok <- !is.na(pp$cand)
  Dm <- as.matrix(sqrt(betadiver(r$X[ok, , drop = FALSE], "sim"))); set.seed(7200 + j$de)
  x <- pares_classes(Dm, pp$cand[ok], r$meta$bloco[ok], un_min = 5, n_classe = 10, B = 100, pares = pp$pares[j$de:j$ate, , drop = FALSE])
  if (!is.null(x)) cbind(grupo = pp$grupo, x)
}
cl <- makeCluster(min(14, nrow(jobs))); clusterExport(cl, c("raiz", "prep"))
pares <- bind_rows(parLapplyLB(cl, split(jobs, seq_len(nrow(jobs))), testa)); stopCluster(cl)
pares <- pares |> filter(!is.na(separa)); cat(format(Sys.time()), "pares testados:", nrow(pares), "\n")
consenso <- pares |> group_by(classe_a, classe_b) |> summarise(grupos = n(), distintos = sum(separa), .groups = "drop") |>
  mutate(frac = distintos / grupos, consenso = ifelse(grupos >= 2 & frac <= 1/3, "equivalente", ifelse(frac >= 2/3, "distinto", "depende do grupo")),
    fraco = grupos == 1 & distintos == 0)
print(table(consenso$consenso, consenso$fraco))

# ---------------- 3 tolerance tree ----------------
n <- nrow(cand); D <- as.matrix(dist(CT)); dimnames(D) <- list(cand$cand, cand$cand)
ia <- match(consenso$classe_a, cand$cand); ib <- match(consenso$classe_b, cand$cand)
dq <- D[cbind(ia, ib)][consenso$consenso == "equivalente"]; dd <- D[cbind(ia, ib)][consenso$consenso == "distinto"]
cat("pares equivalentes:", length(dq), "| distintos:", length(dd), "\n")
taxa <- ecdf(dd)
Tx <- matrix(taxa(D), n, n, dimnames = dimnames(D))
v <- ifelse(consenso$fraco, taxa(D[cbind(ia, ib)]), ifelse(consenso$consenso == "equivalente", 0, NUNCA))
Tx[cbind(ia, ib)] <- v; Tx[cbind(ib, ia)] <- v
Tx[outer(cand$zona, cand$zona, "!=")] <- NUNCA; diag(Tx) <- 0
h_base <- taxa(quantile(dq, Q_BASE)); h_niv <- taxa(quantile(dq, Q_NIVEIS))
cat(sprintf("distancia GF: equivalentes mediana %.3f, distintos mediana %.3f | erro de fusao: base %.3f; niveis %s\n",
  median(dq), median(dd), h_base, paste(sprintf("%.3f", h_niv), collapse = ", ")))
arvore <- hclust(as.dist(Tx), method = "complete"); ct <- cutree(arvore, h = h_base)
rep_ <- tapply(cand$cand, ct, function(x) x[1]); cand$classe <- unname(rep_[as.character(ct)])
eqp <- consenso |> filter(consenso == "equivalente")
cand$fundido_fauna <- vapply(seq_len(n), function(i) { mem <- cand$cand[cand$classe == cand$classe[i]]
  length(mem) > 1 && any(eqp$classe_a %in% mem & eqp$classe_b %in% mem) }, TRUE)
dist_c <- unique(c(consenso$classe_a[consenso$consenso == "distinto"], consenso$classe_b[consenso$consenso == "distinto"]))
testado <- unique(c(consenso$classe_a[!consenso$fraco], consenso$classe_b[!consenso$fraco]))

# ---------------- 4 minimum area (cells per month x depth layer) ----------------
G$classe <- cand$classe[match(G$cand, cand$cand)]; NC <- ncell(g25)
distinta <- function(cl) any(cand$cand[cand$classe == cl] %in% dist_c)
absorvidas <- list(); repeat {
  ci <- G |> count(classe, camada) |> group_by(classe) |> summarise(max_cel = max(n), .groups = "drop") |>
    mutate(zona = cand$zona[match(classe, cand$cand)])
  peq <- ci |> filter(max_cel < AREA_MIN) |> arrange(max_cel); peq <- peq[!vapply(peq$classe, distinta, TRUE), ]
  if (!nrow(peq)) break
  s <- rowsum(TG, G$classe); CTc <- s / as.vector(table(G$classe)[rownames(s)])
  chave <- integer(NC * 84); chave[(G$camada - 1) * NC + G$cel] <- match(G$classe, ci$classe); mudou <- FALSE
  for (k in seq_len(nrow(peq))) { a <- peq$classe[k]; i <- which(G$classe == a); if (!length(i)) next
    vz <- as.data.frame(adjacent(g25, unique(G$cel[i]), "rook", pairs = TRUE))
    vz <- merge(data.frame(from = G$cel[i], lay = G$camada[i]), vz, by = "from")      # neighbours in the same month x depth layer
    viz <- unique(ci$classe[chave[(vz$lay - 1) * NC + vz$to]]); viz <- setdiff(viz[!is.na(viz)], a)
    viz <- viz[ci$zona[match(viz, ci$classe)] == peq$zona[k]]
    if (!length(viz)) next
    d <- colSums((t(CTc[viz, , drop = FALSE]) - CTc[a, ])^2); dest <- viz[which.min(d)]
    G$classe[i] <- dest; chave[(G$camada[i] - 1) * NC + G$cel[i]] <- match(dest, ci$classe); cand$classe[cand$classe == a] <- dest
    absorvidas[[length(absorvidas) + 1]] <- data.frame(de = a, para = dest); mudou <- TRUE }
  if (!mudou) break
}
absorvidas <- bind_rows(absorvidas); cat(format(Sys.time()), "classes absorvidas por area:", nrow(absorvidas), "\n")

# ---------------- 5 attributes ----------------
pers <- G |> count(classe, cel, z) |> group_by(classe) |> summarise(persistencia = mean(n) / 12, .groups = "drop")
fin <- G |> group_by(classe) |> summarise(linhas = n(), celulas_max = max(tabulate(camada)), zona = NOME_Z[first(zona)], nivel3 = names(which.max(table(nivel3))),
  hp = paste(ordena(unique(hp)), collapse = " + "), emu = paste(sort(unique(emu)), collapse = " + "), meses = paste(sort(unique(mes)), collapse = ","),
  z_med = median(z), t_z = median(t_z), sst = median(sst), mld = median(mld), .groups = "drop") |>
  left_join(pers, by = "classe") |>
  left_join(cand |> group_by(classe) |> summarise(candidatos = n(), cand_fauna = sum(cand %in% testado), fundida_fauna = any(fundido_fauna),
    distinta_fauna = any(cand %in% dist_c), .groups = "drop"), by = "classe") |>
  mutate(evidencia = case_when(distinta_fauna ~ "distinta pela fauna", fundida_fauna ~ "fundida pela fauna",
      candidatos > 1 ~ "fundida pelo ambiente", cand_fauna > 0 ~ "testada, sem par conclusivo", TRUE ~ "não testada"),
    standby = evidencia %in% c("não testada", "testada, sem par conclusivo", "fundida pelo ambiente"),
    hp3 = sub(".*\\|", "", nivel3), nivel3_rotulo = paste(c("epi", "meso", "bati")[as.integer(sub("\\|.*", "", nivel3))], hp3, sep = " · "),
    evidencia_nivel3 = ifelse(hp_standby[hp3], "HP em stand-by", "HP validado pela fauna"),
    regime = ifelse(persistencia >= 0.75, "residente", ifelse(persistencia >= 0.25, "sazonal", "ocasional"))) |>
  arrange(match(zona, NOME_Z), as.integer(sub("HP", "", hp3)), desc(linhas))
fin$codigo <- paste0("PF", seq_len(nrow(fin))); G$codigo <- fin$codigo[match(G$classe, fin$classe)]
cat("classes base:", nrow(fin), "| por evidencia:", paste(names(table(fin$evidencia)), table(fin$evidencia), collapse = "; "), "\n")

# ---------------- 3 tolerance levels above the base ----------------
mem <- split(match(cand$cand, rownames(Tx)), fin$codigo[match(cand$classe, fin$classe)])[fin$codigo]
nb <- length(mem); Tc <- matrix(0, nb, nb, dimnames = list(fin$codigo, fin$codigo))
for (i in seq_len(nb - 1)) for (j in (i + 1):nb) Tc[i, j] <- Tc[j, i] <- max(Tx[mem[[i]], mem[[j]]])
arv_c <- hclust(as.dist(Tc), method = "complete"); niveis <- data.frame(codigo = fin$codigo)
for (k in seq_along(h_niv)) { g <- cutree(arv_c, h = h_niv[k]); niveis[[sprintf("tol_%02.0f", 100 * h_niv[k])]] <- unname(tapply(fin$codigo, g, function(x) x[which.max(fin$linhas[match(x, fin$codigo)])])[as.character(g)]) }
hier <- data.frame(passo = seq_len(nrow(arv_c$merge)), a = arv_c$merge[, 1], b = arv_c$merge[, 2], altura = arv_c$height) |> filter(altura < NUNCA)
cat("classes por nivel de tolerancia:", paste(names(niveis)[-1], sapply(niveis[-1], function(x) length(unique(x))), collapse = "; "), "\n")

# ---------------- 6 acceptance ----------------
avalia <- function(pp, cls) {
  suppressPackageStartupMessages({library(vegan); library(dplyr)}); source(file.path(raiz, "scripts/funcoes_permanova.R"))
  r <- readRDS(pp$arq); Dm <- as.matrix(sqrt(betadiver(r$X, "sim"))); bl <- r$meta$bloco
  bind_rows(lapply(names(cls), function(nm) { g <- cls[[nm]]; k <- !is.na(g) & g %in% names(which(table(g) >= 5))
    if (length(unique(g[k])) < 2) return(NULL)
    set.seed(9); x <- permanova_rapida(as.dist(Dm[k, k]), g[k], bl[k], 199)
    set.seed(10); e <- tryCatch(permanova_esforco_igual(Dm[k, k], g[k], bl[k], n_classe = 10, B = 100), error = function(e) NULL)
    data.frame(grupo = pp$grupo, classificacao = nm, classes = length(unique(g[k])), unidades = sum(k), omega2 = x$omega2,
      omega2_eq = if (!is.null(e)) e$omega2_eq else NA, poder_eq = if (!is.null(e)) e$poder_eq else NA) })) }
cod_cand <- setNames(fin$codigo[match(cand$classe, fin$classe)], cand$cand)
cl <- makeCluster(length(prep)); clusterExport(cl, c("raiz", "avalia"))
aval <- bind_rows(parLapply(cl, prep, function(pp, cod) avalia(pp, list("mapa unido (base)" = unname(cod[pp$cand]), "nível III" = pp$nivel3,
  "habitats novos (HP)" = pp$hp, "EMU" = as.character(pp$emu))), cod = cod_cand)); stopCluster(cl)

# ---------------- colors, outputs and maps ----------------
lh <- read.csv(file.path(raiz, "saida/habitats_novos_legenda.csv")) |> filter(grepl("^pelagico", arquivo)) |> distinct(habitat, cor)
corHP <- setNames(lh$cor, lh$habitat)
fin$cor <- unlist(unname(lapply(split(seq_len(nrow(fin)), fin$nivel3), function(i) { m <- length(i); f <- if (m == 1) 0 else seq(-0.45, 0.45, length.out = m)
  setNames(vapply(f, function(x) { c0 <- col2rgb(corHP[fin$hp3[i[1]]]) / 255; c1 <- if (x < 0) c0 * (1 + x) else c0 + (1 - c0) * x; rgb(c1[1], c1[2], c1[3]) }, ""), fin$codigo[i]) })))[fin$codigo]
fin <- fin |> mutate(valor = seq_len(n())) |> left_join(niveis, by = "codigo")
n3 <- unique(fin$nivel3); fin$valor_nivel3 <- match(fin$nivel3, n3)
write.csv(fin, file.path(raiz, "saida/mapa_unido_pelagico_classes.csv"), row.names = FALSE)
write.csv(consenso, file.path(raiz, "saida/mapa_unido_pelagico_pares.csv"), row.names = FALSE)
write.csv(aval, file.path(raiz, "saida/mapa_unido_pelagico_avaliacao.csv"), row.names = FALSE)
write.csv(absorvidas, file.path(raiz, "saida/mapa_unido_pelagico_absorvidas.csv"), row.names = FALSE)
write.csv(hier, file.path(raiz, "saida/mapa_unido_pelagico_hierarquia.csv"), row.names = FALSE)
write.csv(data.frame(nivel = c("base", names(niveis)[-1]), percentil_equivalentes = c(Q_BASE, Q_NIVEIS), erro_fusao = c(h_base, h_niv),
  classes = c(nrow(fin), sapply(niveis[-1], function(x) length(unique(x))))), file.path(raiz, "saida/mapa_unido_pelagico_niveis.csv"), row.names = FALSE)
G$valor <- fin$valor[match(G$codigo, fin$codigo)]
camada_r <- function(m, z, col = "valor") { r <- rast(g25); v <- rep(NA, ncell(r)); s <- G[G$mes == m & G$z == z, ]; v[s$cel] <- s[[col]]; r[] <- v; names(r) <- sprintf("m%02d_z%d", m, z); r }
st <- rast(lapply(seq_len(12 * length(ZP)), function(i) camada_r((i - 1) %/% length(ZP) + 1, ZP[(i - 1) %% length(ZP) + 1])))
writeRaster(st, file.path(raiz, "saida/mapa_unido_pelagico_mes_z.tif"), overwrite = TRUE, datatype = "INT2U")
hs <- ordena(unique(fin$hp3))
mapa_base <- function(m, z, arq, titulo) { rb <- st[[sprintf("m%02d_z%d", m, z)]]; vb <- values(rb)[, 1]; rc <- rowColFromCell(rb, seq_len(ncell(rb)))
  k <- which(!is.na(vb)); k <- k[fin$standby[vb[k]] & (rc[k, 1] + rc[k, 2]) %% 7 < 2]; h <- rast(rb); hv <- rep(NA, ncell(rb)); hv[k] <- 1; h[] <- hv
  png(file.path(raiz, sprintf("saida/mapa_unido_pelagico_%s.png", arq)), 2000, 1160, res = 150, bg = "white"); layout(matrix(1:2), heights = c(1000, 160)); par(mar = c(0.5, 1, 2, 1))
  plot(rb, col = fin$cor, breaks = c(0.5, seq_len(nrow(fin)) + 0.5), legend = FALSE, axes = FALSE, box = FALSE, maxcell = ncell(rb),
       main = sprintf("%s: %d classes do nível base (erro de fusão ≤ %.0f%%)", titulo, length(unique(vb[!is.na(vb)])), 100 * h_base))
  if (length(k)) plot(h, col = COR_ACH, breaks = c(0.5, 1.5), legend = FALSE, axes = FALSE, box = FALSE, add = TRUE, maxcell = ncell(rb))
  par(mar = c(0, 1, 0, 1)); plot.new()
  legend("top", legend = hs, fill = corHP[hs], ncol = 14, cex = 0.7, bty = "n",
         title = "matiz = habitat pelágico (HP); tons = classes do nível base; achurado = sem evidência da fauna"); dev.off() }
mapa_base(1, 0, "0m_jan", "Mapa unido do pelágico, superfície, janeiro"); mapa_base(7, 0, "0m_jul", "Mapa unido do pelágico, superfície, julho")
mapa_base(1, 200, "200m_jan", "Mapa unido do pelágico, 200 m, janeiro"); mapa_base(7, 200, "200m_jul", "Mapa unido do pelágico, 200 m, julho")
c3 <- corHP[sub(".*\\|", "", n3)]; s3 <- hp_standby[sub(".*\\|", "", n3)]; rot3 <- fin$nivel3_rotulo[match(n3, fin$nivel3)]
for (zz in c(0, 500)) { G$v3 <- fin$valor_nivel3[G$valor]; r3 <- camada_r(1, zz, "v3")
  mapa_classes(r3, c3, s3, rot3, sprintf("Mapa unido do pelágico, nível III (zona × HP), %d m, janeiro", zz), file.path(raiz, sprintf("saida/mapa_unido_pelagico_nivel3_%dm_jan.png", zz)), ncol = 8) }
options(width = 200); print(aval |> mutate(across(where(is.double), ~ round(.x, 3))))
print(table(fin$zona, fin$evidencia)); print(table(fin$regime)); cat(format(Sys.time()), "FIM\n")
