# 99_etapa4.R  (stage 4, supplementary material: PI decision of 06/11/2026; Gate 12 item 4, 13/11/2026)
# 1 Monthly maps of the united pelagic map (base level), surface and 200 m, 12 panels each.
# 2 Validation of the united maps with fauna left out:
#   a blocks out: each unit is classified by its nearest fauna centroid (beta_sim; Anderson 2006) among the classes of the units outside its
#     2,500 km block (classifica_lbo of funcoes_etapa2.R); accuracy, balanced accuracy and kappa, for the united map (base) and its components
#     (benthos: HB, EBU region, EBU group; pelagic: HP, EMU). Groups: level 1, option B and the test fauna c (the last never built a map);
#   b years out: units rebuilt in two periods split at the median event year of each group; units of one period are classified by the class
#     centroids of the other period (both directions), only classes with >= 5 units in the training period.
# 3 Area of applicability (Meyer and Pebesma 2021), in the GF-transformed space that built the habitats: dissimilarity index DI = distance
#   to the nearest training unit / mean distance between training units; threshold = upper whisker (Q3 + 1.5 IQR) of the DI of the training
#   units against units outside their 2,500 km block. Training = units of all groups (level 1, option B, c); benthos on the 25 km cells,
#   pelagic on every cell x month x standard depth. Area inside the AOA per class of the united maps.
# Outputs: saida/etapa4_{validacao_blocos,validacao_anos,aoa_classes,aoa_resumo}.csv, saida/etapa4_aoa_bentos_25km.tif,
#          saida/etapa4_aoa_pelagico_mes_z.tif, saida/etapa4_*.png
suppressPackageStartupMessages({library(dplyr); library(terra); library(arrow); library(parallel); library(FNN); library(tidyr)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); source(file.path(raiz, "scripts/funcoes_grade.R")); source(file.path(raiz, "scripts/funcoes_etapa2.R"))
source(file.path(raiz, "scripts/funcoes_mapa.R"))
g25 <- ease2_grade("25km"); ZP <- c(0, 50, 100, 200, 500, 1000, 2000); zona <- function(z) findInterval(z, c(200, 1000)) + 1
ordena <- function(h) h[order(as.integer(sub("^H[BP]", "", h)))]
cel_de <- function(id) { p <- do.call(rbind, strsplit(id, "_")); cellFromRowCol(g25, as.integer(p[, 2]), as.integer(p[, 3])) }
M <- lapply(1:12, function(m) rast(file.path(raiz, sprintf("camadas/glorys/mes/glorys_%02d_25km.tif", m)))); pg <- values(M[[1]][["prof_glorys"]])[, 1]
fb <- read.csv(file.path(raiz, "saida/mapa_unido_bentos_classes.csv"), encoding = "UTF-8"); ub <- rast(file.path(raiz, "saida/mapa_unido_bentos_25km.tif"))
fp <- read.csv(file.path(raiz, "saida/mapa_unido_pelagico_classes.csv"), encoding = "UTF-8"); up <- rast(file.path(raiz, "saida/mapa_unido_pelagico_mes_z.tif"))
hbv <- values(rast(file.path(raiz, "saida/habitats_novos_bentos_25km.tif")))[, 1]; lh <- read.csv(file.path(raiz, "saida/habitats_novos_legenda.csv"))
hb_nome <- lh$habitat[lh$arquivo == "bentos_25km"]; bv <- values(ub[["base"]])[, 1]
hpt <- rast(file.path(raiz, "saida/habitats_novos_pelagico_mes_z.tif")); fus <- read.csv(file.path(raiz, "saida/pelagico_b_fusoes.csv"))
nf <- tapply(fus$hp_n1, fus$final, function(x) ordena(unique(x))[1]); hp_de <- setNames(as.character(nf[as.character(fus$final)]), fus$hp_n1); hs0 <- ordena(unique(fus$hp_n1))
e25 <- read_parquet(file.path(raiz, "camadas/ebu/ebu_celula_25km.parquet")); e25$cel <- cel_de(e25$cell_id)
le <- function(st, cel, lay) { v <- rep(NA_integer_, length(cel)); for (l in unique(lay)) { i <- which(lay == l); v[i] <- values(st[[l]])[cel[i], 1] }; v }
camada_p <- function(cel, tri, z) { m <- c(T1 = 2, T2 = 5, T3 = 8, T4 = 11)[tri]; sprintf("m%02d_z%d", m, ZP[vapply(z, function(x) which.min(abs(ZP - x)), 1L)]) }
classes_b <- function(cel) { b <- bv[cel]; k <- match(cel, e25$cel)
  list("mapa unido (base)" = fb$codigo[match(b, fb$valor)], "habitats novos (HB)" = hb_nome[hbv[cel]], "EBU região" = e25$ebu_regiao[k], "EBU grupo" = e25$ebu_grupo[k]) }
classes_p <- function(cel, tri, z, emu) { lay <- camada_p(cel, tri, z)
  list("mapa unido (base)" = fp$codigo[match(le(up, cel, lay), fp$valor)], "habitats novos (HP)" = unname(hp_de[hs0[le(hpt, cel, lay)]]), "EMU" = as.character(emu)) }

# ---------------- 1 monthly maps ----------------
for (zz in c(0, 200)) { png(file.path(raiz, sprintf("saida/etapa4_mensal_pelagico_%dm.png", zz)), 2400, 1500, res = 150, bg = "white"); par(mfrow = c(4, 3), mar = c(0.3, 0.3, 1.6, 0.3))
  for (m in 1:12) plot(up[[sprintf("m%02d_z%d", m, zz)]], col = fp$cor, breaks = c(0.5, seq_len(nrow(fp)) + 0.5), legend = FALSE, axes = FALSE, box = FALSE,
    main = sprintf("%s, %d m", month.abb[m], zz), cex.main = 1, maxcell = 2e6)
  dev.off() }
cat(format(Sys.time()), "mapas mensais\n")

# ---------------- 2a blocks out ----------------
dirs <- file.path(raiz, c("intermediario/etapa3_n1", "intermediario/etapa3_n1b", "intermediario/etapa3_n1c"))
arqs <- unlist(lapply(dirs, function(d) list.files(d, pattern = "^(bentos|pelagico)_.*[.]rds$", full.names = TRUE)))
arqs <- arqs[!grepl("Euphausiacea|Hyperiidea", arqs)]
fonte <- function(a) if (grepl("n1b", a)) "opção B" else if (grepl("n1c", a)) "fauna de teste" else "nível 1"
blocos <- function(a) { suppressPackageStartupMessages({library(vegan); library(dplyr); library(terra)}); r <- readRDS(a)
  cls <- if (r$comp == "bentos") classes_b(cel_de(r$meta$unidade)) else classes_p(cel_de(r$meta$cell_id), r$meta$tri, r$E$z, r$meta$emu)
  S <- as.matrix(betadiver(r$X, "sim")); S[is.na(S)] <- 1
  bind_rows(lapply(names(cls), function(nm) { g <- cls[[nm]]; k <- !is.na(g) & g %in% names(which(table(g) >= 5)); if (length(unique(g[k])) < 2) return(NULL)
    m <- metricas_classif(classifica_lbo(S[k, k], g[k], r$meta$bloco[k])); if (is.null(m)) return(NULL)
    cbind(comp = r$comp, grupo = r$grupo, fonte = fonte(a), classificacao = nm, classes = length(unique(g[k])), m$geral) })) }
cl <- makeCluster(min(14, length(arqs))); clusterExport(cl, c("raiz", "classes_b", "classes_p", "camada_p", "le", "cel_de", "fonte", "classifica_lbo", "metricas_classif",
  "bv", "fb", "fp", "hb_nome", "hbv", "e25", "hp_de", "hs0", "ZP"))
invisible(clusterEvalQ(cl, { library(terra); source(file.path(raiz, "scripts/funcoes_grade.R")); g25 <- ease2_grade("25km"); up <- rast(file.path(raiz, "saida/mapa_unido_pelagico_mes_z.tif")); hpt <- rast(file.path(raiz, "saida/habitats_novos_pelagico_mes_z.tif")); NULL }))
vb <- bind_rows(parLapplyLB(cl, arqs, blocos)); stopCluster(cl)
write.csv(vb, file.path(raiz, "saida/etapa4_validacao_blocos.csv"), row.names = FALSE); cat(format(Sys.time()), "validacao por blocos\n")

# ---------------- 2b years out ----------------
tabs <- list(bentos = paste0("saida/tabela_certa_bentos_n1", c("", "c"), ".parquet"), pelagico = paste0("saida/tabela_certa_pelagico_n1", c("", "b", "c"), ".parquet"))
incid10 <- function(d) { set.seed(1001); s <- d |> distinct(unidade, evento) |> group_by(unidade) |> slice_sample(n = 10) |> ungroup()
  m <- d |> semi_join(s, by = c("unidade", "evento")) |> distinct(unidade, taxon_unidade) |> mutate(v = 1L) |> pivot_wider(names_from = taxon_unidade, values_from = v, values_fill = 0L)
  X <- as.matrix(m[, -1]); rownames(X) <- m$unidade; X }
classifica_tempo <- function(S, g, treino, teste) {   # nearest centroid of the training period (Anderson 2006), as classifica_lbo
  n <- table(g[treino]); cls <- names(n)[n >= 5]; if (length(cls) < 2) return(NULL)
  d2 <- sapply(cls, function(c) { j <- which(treino & g == c); rowMeans(S[teste, j, drop = FALSE]) - sum(S[j, j]) / (2 * length(j)^2) })
  if (is.null(dim(d2))) d2 <- matrix(d2, nrow = 1, dimnames = list(NULL, cls))
  data.frame(classe = g[teste], prevista = cls[apply(d2, 1, which.min)], classificavel = g[teste] %in% cls) |> mutate(prevista = ifelse(classificavel, prevista, NA), correta = prevista == classe) }
anos <- function(cp, f) { d <- open_dataset(file.path(raiz, f)) |> filter(year >= 1993, nivel_certeza %in% c("A", "B")) |>
    select(any_of(c("grupo", "cell_id", "event_id", "taxon_unidade", "year", "month", "depth_min_m", "depth_max_m", "bloco_espacial", "emu"))) |> collect() |>
    filter(!is.na(taxon_unidade), taxon_unidade != "", !grupo %in% c("Euphausiacea", "Hyperiidea"))
  if (cp == "pelagico") d <- d |> filter(!is.na(month), !is.na(emu)) |> mutate(tri = paste0("T", (month - 1) %/% 3 + 1), zreg = rowMeans(cbind(depth_min_m, depth_max_m), na.rm = TRUE))
  bind_rows(lapply(split(d, d$grupo), function(x) { corte <- median(x$year[!duplicated(x$event_id)]); x$periodo <- ifelse(x$year <= corte, "P1", "P2")
    x$unidade <- if (cp == "bentos") paste(x$cell_id, x$periodo) else paste(x$cell_id, x$emu, x$tri, x$periodo); x$evento <- x$event_id
    un <- x |> group_by(unidade) |> summarise(ev = n_distinct(evento), cell_id = first(cell_id), periodo = first(periodo), bloco = first(bloco_espacial),
      tri = if (cp == "pelagico") first(tri) else NA, emu = if (cp == "pelagico") first(emu) else NA, z = if (cp == "pelagico") median(zreg, na.rm = TRUE) else NA, .groups = "drop") |>
      filter(ev >= 10, cp == "bentos" | !is.na(z))
    if (min(table(factor(un$periodo, c("P1", "P2")))) < 20) return(NULL)
    X <- incid10(x |> filter(unidade %in% un$unidade)); un <- un[match(rownames(X), un$unidade), ]; cel <- cel_de(un$cell_id)
    cls <- if (cp == "bentos") classes_b(cel) else classes_p(cel, un$tri, un$z, un$emu)
    S <- as.matrix(vegan::betadiver(X, "sim")); S[is.na(S)] <- 1
    bind_rows(lapply(names(cls), function(nm) { g <- cls[[nm]]; bind_rows(lapply(c("P1", "P2"), function(tr) { ok <- !is.na(g)
      c1 <- classifica_tempo(S[ok, ok], g[ok], un$periodo[ok] == tr, un$periodo[ok] != tr); if (is.null(c1)) return(NULL); m <- metricas_classif(c1); if (is.null(m)) return(NULL)
      cbind(comp = cp, grupo = x$grupo[1], classificacao = nm, treino = tr, ano_corte = corte, unidades_P1 = sum(un$periodo == "P1"), unidades_P2 = sum(un$periodo == "P2"), m$geral) })) })) })) }
va <- bind_rows(lapply(names(tabs), function(cp) bind_rows(lapply(tabs[[cp]][file.exists(file.path(raiz, tabs[[cp]]))], function(f) anos(cp, f)))))
write.csv(va, file.path(raiz, "saida/etapa4_validacao_anos.csv"), row.names = FALSE); cat(format(Sys.time()), "validacao por anos\n")

# ---------------- 3 area of applicability ----------------
cv <- open_dataset(file.path(raiz, "saida/etapa3_n1_gf_curvas.csv"), format = "csv") |> collect() |> as.data.frame()
transf <- function(E, cp) { c0 <- cv[cv$comp == cp, ]; vs <- unique(c0$preditor)
  sapply(vs, function(v) { s <- 0; for (g in unique(c0$grupo[c0$preditor == v])) { k <- c0[c0$grupo == g & c0$preditor == v, ]; s <- s + approx(k$x, k$y, E[[v]], rule = 2, ties = mean)$y }; s }) }
treino <- function(cp) { a <- arqs[grepl(paste0("/", cp, "_"), arqs)]
  T <- do.call(rbind, lapply(a, function(f) { r <- readRDS(f); t0 <- transf(r$E, cp); attr(t0, "bloco") <- r$meta$bloco; cbind(as.data.frame(t0), bloco = r$meta$bloco) }))
  T <- T[complete.cases(T), ]; list(X = as.matrix(T[, setdiff(names(T), "bloco")]), bloco = T$bloco) }
limiar <- function(tr) { set.seed(3); i <- sample(nrow(tr$X), min(3000, nrow(tr$X))); dm <- mean(dist(tr$X[i, ]))
  di <- vapply(unique(tr$bloco), function(b) 0, 0); d <- rep(NA_real_, nrow(tr$X))
  for (b in unique(tr$bloco)) { k <- tr$bloco == b; if (all(k)) next; d[k] <- get.knnx(tr$X[!k, , drop = FALSE], tr$X[k, , drop = FALSE], k = 1)$nn.dist[, 1] }
  q <- quantile(d / dm, c(0.25, 0.75), na.rm = TRUE); list(dm = dm, thr = q[2] + 1.5 * diff(q), di_treino = d / dm) }
# benthos
trB <- treino("bentos"); lmB <- limiar(trB); cb <- e25$cel
vm <- function(camada) sapply(1:12, function(k) values(M[[k]][[camada]])[cb, 1])
tf <- vm("tfundo"); vf <- vm("vfundo"); ss <- vm("thetao_0"); ml <- vm("mld"); ge <- vm("gelo")
GB <- data.frame(prof = pg[cb], tfundo_med = rowMeans(tf), tfundo_min = apply(tf, 1, min), tfundo_max = apply(tf, 1, max), vfundo_med = rowMeans(vf), vfundo_max = apply(vf, 1, max),
  sst_med = rowMeans(ss), sst_amp = apply(ss, 1, max) - apply(ss, 1, min), mld_max = apply(ml, 1, max), gelo_meses = rowSums(ge >= 0.15, na.rm = TRUE)); rm(tf, vf, ss, ml, ge)
okB <- complete.cases(GB); TB <- transf(GB[okB, ], "bentos")[, colnames(trB$X)]
diB <- get.knnx(trB$X, TB, k = 1)$nn.dist[, 1] / lmB$dm
rB <- rast(g25); v <- rep(NA, ncell(rB)); v[cb[okB]] <- as.integer(diB <= lmB$thr); rB[] <- v; names(rB) <- "aoa"
writeRaster(rB, file.path(raiz, "saida/etapa4_aoa_bentos_25km.tif"), overwrite = TRUE, datatype = "INT1U")
aoaB <- data.frame(cel = cb[okB], aoa = diB <= lmB$thr, di = diB) |> mutate(codigo = fb$codigo[match(bv[cel], fb$valor)])
cat(format(Sys.time()), "AOA bentos: limiar", round(lmB$thr, 3), "| celulas dentro:", round(mean(aoaB$aoa), 3), "\n")
# pelagic, layer by layer
trP <- treino("pelagico"); lmP <- limiar(trP); emu <- rast(file.path(raiz, "camadas/emu/emu_25km.tif"))
aoaP <- list(); camadas <- list()
for (m in 1:12) for (z in ZP) { v <- values(M[[m]][[c(paste0("thetao_", z), paste0("so_", z), "thetao_0", "so_0", "mld", "vsup", "gelo")]]); ok <- which(!is.na(v[, 1]) & !is.na(pg) & pg >= z)
  E <- data.frame(z = z, t_z = v[ok, 1], s_z = v[ok, 2], sst = v[ok, 3], sss = v[ok, 4], mld = v[ok, 5], prof = pg[ok], vsup = v[ok, 6], gelo = ifelse(is.na(v[ok, 7]), 0, v[ok, 7]))
  E$z_mld <- log2((E$z + 1) / (E$mld + 1)); k <- complete.cases(E); ok <- ok[k]; E <- E[k, ]
  di <- get.knnx(trP$X, transf(E, "pelagico")[, colnames(trP$X), drop = FALSE], k = 1)$nn.dist[, 1] / lmP$dm; lay <- sprintf("m%02d_z%d", m, z)
  r <- rast(g25); vv <- rep(NA, ncell(r)); vv[ok] <- as.integer(di <= lmP$thr); r[] <- vv; names(r) <- lay; camadas[[lay]] <- r
  aoaP[[lay]] <- data.frame(mes = m, z = z, aoa = di <= lmP$thr, codigo = fp$codigo[match(values(up[[lay]])[ok, 1], fp$valor)]) |> count(mes, z, codigo, aoa) }
stP <- rast(camadas); writeRaster(stP, file.path(raiz, "saida/etapa4_aoa_pelagico_mes_z.tif"), overwrite = TRUE, datatype = "INT1U")
aoaP <- bind_rows(aoaP); cat(format(Sys.time()), "AOA pelagico: limiar", round(lmP$thr, 3), "| linhas dentro:", round(sum(aoaP$n[aoaP$aoa]) / sum(aoaP$n), 3), "\n")
cls <- bind_rows(aoaB |> filter(!is.na(codigo)) |> group_by(codigo) |> summarise(n = n(), dentro = mean(aoa), .groups = "drop") |> mutate(comp = "bentos") |>
    left_join(fb |> select(codigo, faixa_zona = faixa, evidencia), by = "codigo"),
  aoaP |> filter(!is.na(codigo)) |> group_by(codigo) |> summarise(dentro = sum(n[aoa]) / sum(n), n = sum(n), .groups = "drop") |> mutate(comp = "pelagico") |>
    left_join(fp |> select(codigo, faixa_zona = zona, evidencia), by = "codigo"))
write.csv(cls, file.path(raiz, "saida/etapa4_aoa_classes.csv"), row.names = FALSE)
res <- cls |> group_by(comp, faixa_zona) |> summarise(classes = n(), area_dentro = sum(n * dentro) / sum(n), classes_maioria_fora = sum(dentro < 0.5), .groups = "drop") |>
  bind_rows(cls |> group_by(comp, faixa_zona = paste("evidência:", evidencia)) |> summarise(classes = n(), area_dentro = sum(n * dentro) / sum(n), classes_maioria_fora = sum(dentro < 0.5), .groups = "drop"))
write.csv(bind_rows(res, data.frame(comp = c("bentos", "pelagico"), faixa_zona = "limiar DI", area_dentro = c(lmB$thr, lmP$thr))), file.path(raiz, "saida/etapa4_aoa_resumo.csv"), row.names = FALSE)
cor_aoa <- c("#d9d9d9", "#2c7fb8")
for (x in list(list(rB, "bentos", "Área de aplicabilidade, bentos"), list(stP[["m01_z0"]], "pelagico_0m_jan", "Área de aplicabilidade, pelágico, superfície, janeiro"),
               list(stP[["m01_z500"]], "pelagico_500m_jan", "Área de aplicabilidade, pelágico, 500 m, janeiro"))) {
  png(file.path(raiz, sprintf("saida/etapa4_aoa_%s.png", x[[2]])), 2000, 1100, res = 150, bg = "white"); par(mar = c(0.5, 0.5, 2, 0.5))
  plot(x[[1]], col = cor_aoa, breaks = c(-0.5, 0.5, 1.5), legend = FALSE, axes = FALSE, box = FALSE, main = x[[3]], maxcell = 2e6)
  legend("bottomleft", fill = cor_aoa, legend = c("fora (ambiente sem análogo nos dados)", "dentro"), bty = "n", cex = 0.8); dev.off() }
options(width = 220); print(vb |> select(comp, grupo, fonte, classificacao, classes, acerto, acerto_equilibrado, kappa) |> mutate(across(where(is.double), ~ round(.x, 3))) |> as.data.frame())
print(va |> select(comp, grupo, classificacao, treino, acerto, kappa) |> mutate(across(where(is.double), ~ round(.x, 3))) |> as.data.frame()); print(as.data.frame(res))
cat(format(Sys.time()), "FIM\n")
