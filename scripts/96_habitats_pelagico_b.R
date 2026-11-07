# 96_habitats_pelagico_b.R  (option B, Gate 11 item 5: parallel version of the pelagic habitats, compared with the approved one)
# Same steps as the pelagic part of 89b_habitats_novos.R: GF transformation (level-1 pelagic curves + curves of the 21b groups), k-means with
# K0 = 30 on the cell x month x depth grid, nearest-centroid assignment of the fauna units, pairs per group, consensus and merging by complete
# linkage with a 1/3 cut. Fauna: level-1 pelagic groups (without Euphausiacea and Hyperiidea, Gate 8), the 21b groups and the level-2 groups.
# Acceptance rule (Gate 11): the new version replaces the approved one only if omega2 improves. Comparison per group on the same units, with
# all data and with equal effort. Groups that built a version favour it, so level-1 groups and the new groups are reported apart.
# Nothing approved is overwritten. Outputs: saida/opcaoB_{habitats,pares,comparacao,zonas}.csv, saida/opcaoB_pelagico_mes_z.tif,
#          intermediario/etapa3_n1b/opcaoB_modelo.rds
suppressPackageStartupMessages({library(dplyr); library(terra); library(arrow); library(parallel)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); source(file.path(raiz, "scripts/funcoes_grade.R")); source(file.path(raiz, "scripts/funcoes_permanova.R"))
K0 <- 30; ZP <- c(0, 50, 100, 200, 500, 1000, 2000); SEM_SINAL <- c("Euphausiacea", "Hyperiidea"); zona <- function(z) findInterval(z, c(200, 1000)) + 1
ordena <- function(h) h[order(as.integer(sub("^HP", "", h)))]
g25 <- ease2_grade("25km"); M <- lapply(1:12, function(m) rast(file.path(raiz, sprintf("camadas/glorys/mes/glorys_%02d_25km.tif", m))))
c1 <- open_dataset(file.path(raiz, "saida/etapa3_n1_gf_curvas.csv"), format = "csv") |> filter(comp == "pelagico") |> collect() |> as.data.frame()
c2 <- read.csv(file.path(raiz, "saida/etapa3_n1_gf_curvas_b.csv"))
curvas_a <- c1; curvas_b <- bind_rows(c1, c2)
transf <- function(E, cv) { vs <- intersect(unique(cv$preditor), names(E))
  sapply(vs, function(v) { s <- 0; for (g in unique(cv$grupo[cv$preditor == v])) { k <- cv[cv$grupo == g & cv$preditor == v, ]; s <- s + approx(k$x, k$y, E[[v]], rule = 2, ties = mean)$y }; s }) }
mais_proximo <- function(T, C) { d <- sapply(seq_len(nrow(C)), function(k) colSums((t(T) - C[k, ])^2)); max.col(-d, ties.method = "first") }
# ---------------- grid ----------------
pg <- values(M[[1]][["prof_glorys"]])[, 1]; G <- list()
for (m in 1:12) for (z in ZP) { v <- values(M[[m]][[c(paste0("thetao_", z), paste0("so_", z), "thetao_0", "so_0", "mld", "vsup", "gelo")]]); ok <- which(!is.na(v[, 1]) & !is.na(pg) & pg >= z)
  G[[length(G) + 1]] <- data.frame(cel = ok, mes = m, z = z, t_z = v[ok, 1], s_z = v[ok, 2], sst = v[ok, 3], sss = v[ok, 4], mld = v[ok, 5], prof = pg[ok], vsup = v[ok, 6], gelo = ifelse(is.na(v[ok, 7]), 0, v[ok, 7])) }
G <- bind_rows(G); G$z_mld <- log2((G$z + 1) / (G$mld + 1)); G <- G[complete.cases(G[, c("t_z", "s_z", "sst", "sss", "mld", "vsup")]), ]
TP <- transf(G, curvas_b); cat(format(Sys.time()), "linhas:", nrow(G), "\n")
set.seed(12); ip <- sample(nrow(TP), min(300000, nrow(TP))); km <- kmeans(TP[ip, ], K0, nstart = 5, iter.max = 100, algorithm = "MacQueen")
G$k0 <- unlist(lapply(split(seq_len(nrow(TP)), ceiling(seq_len(nrow(TP)) / 2e5)), function(i) mais_proximo(TP[i, , drop = FALSE], km$centers)))
rm(TP); gc(verbose = FALSE); cat(format(Sys.time()), "k-means feito\n")
# ---------------- fauna ----------------
arqs <- c(list.files(file.path(raiz, "intermediario/etapa3_n1"), pattern = "^pelagico_.*\\.rds$", full.names = TRUE), list.files(file.path(raiz, "intermediario/etapa3_n1b"), pattern = "^pelagico_.*\\.rds$", full.names = TRUE),
          list.files(file.path(raiz, "intermediario/pelagico_n2"), pattern = "\\.rds$", full.names = TRUE))
arqs <- arqs[!grepl(paste(SEM_SINAL, collapse = "|"), arqs)]
mod <- readRDS(file.path(raiz, "intermediario/etapa3_n1/habitats_novos_modelo.rds")); CPa <- mod$kmP$centers; rm(mod); gc(verbose = FALSE)
fus <- read.csv(file.path(raiz, "saida/pelagico_n2_fusoes.csv")); nf <- tapply(fus$hp_n1, fus$final_n1n2, function(x) ordena(unique(x))[1]); hp_a <- as.character(nf[as.character(fus$final_n1n2)])
tag <- function(a) if (grepl("pelagico_n2", a)) " (n2)" else if (grepl("etapa3_n1b", a)) " (B)" else ""
U <- lapply(arqs, function(a) { r <- readRDS(a); list(arq = a, grupo = paste0(r$grupo, tag(a)), n1 = tag(a) == "", k = mais_proximo(transf(r$E, curvas_b), km$centers), hpa = hp_a[mais_proximo(transf(r$E, curvas_a), CPa)], z = r$E$z) })
testa <- function(u) { suppressPackageStartupMessages({library(vegan); library(dplyr)}); source(file.path(raiz, "scripts/funcoes_permanova.R")); r <- readRDS(u$arq); Dm <- as.matrix(sqrt(betadiver(r$X, "sim"))); set.seed(7400)
  x <- pares_classes(Dm, as.character(u$k), r$meta$bloco, un_min = 5, n_classe = 10, B = 100); if (!is.null(x)) cbind(grupo = u$grupo, x) }
cl <- makeCluster(min(10, length(U))); clusterExport(cl, "raiz"); pares <- bind_rows(parLapplyLB(cl, U, testa)); stopCluster(cl); pares <- pares |> filter(!is.na(separa))
cat(format(Sys.time()), "pares:", nrow(pares), "\n")
cons <- pares |> mutate(a = pmin(as.integer(classe_a), as.integer(classe_b)), b = pmax(as.integer(classe_a), as.integer(classe_b))) |> group_by(classe_a = a, classe_b = b) |> summarise(grupos = n(), distintos = sum(separa), .groups = "drop") |>
  mutate(frac = distintos / grupos, consenso = ifelse(grupos >= 2 & frac <= 1/3, "equivalente", ifelse(frac >= 2/3, "distinto", "depende do grupo")))
Mx <- matrix(1, K0, K0); diag(Mx) <- 0; Mx[cbind(cons$classe_a, cons$classe_b)] <- ifelse(cons$grupos >= 2, cons$frac, 1); Mx[cbind(cons$classe_b, cons$classe_a)] <- Mx[cbind(cons$classe_a, cons$classe_b)]
g <- cutree(hclust(as.dist(Mx), method = "complete"), h = 1/3); fin <- match(g, unique(g)); te <- (1:K0) %in% c(cons$classe_a, cons$classe_b)
FP <- data.frame(k0 = 1:K0, final = fin, testado = te, habitat = paste0("HPb", fin)); FP$standby <- !tapply(FP$testado, FP$final, any)[as.character(FP$final)]
G$habitat <- FP$habitat[G$k0]; G$zona <- zona(G$z)
zn <- G |> distinct(habitat, zona) |> left_join(FP |> distinct(habitat, standby), by = "habitat") |> count(zona, standby)
pz <- bind_rows(lapply(U, function(u) data.frame(grupo = u$grupo, zona = zona(u$z)))) |> count(grupo, zona)
# ---------------- comparison ----------------
avalia <- function(u, FP) { suppressPackageStartupMessages({library(vegan); library(dplyr)}); source(file.path(raiz, "scripts/funcoes_permanova.R")); r <- readRDS(u$arq); Dm <- as.matrix(sqrt(betadiver(r$X, "sim"))); bl <- r$meta$bloco
  cls <- list(aprovado = u$hpa, versao_B = FP$habitat[u$k]); bind_rows(lapply(names(cls), function(nm) { g <- cls[[nm]]; k <- g %in% names(which(table(g) >= 5)); if (length(unique(g[k])) < 2) return(NULL)
    set.seed(9); x <- permanova_rapida(as.dist(Dm[k, k]), g[k], bl[k], 199); set.seed(10); e <- tryCatch(permanova_esforco_igual(Dm[k, k], g[k], bl[k], n_classe = 10, B = 100), error = function(e) NULL)
    data.frame(grupo = u$grupo, nivel_1_original = u$n1, versao = nm, classes = length(unique(g[k])), unidades = sum(k), omega2 = x$omega2, omega2_eq = if (!is.null(e)) e$omega2_eq else NA, poder_eq = if (!is.null(e)) e$poder_eq else NA) })) }
cl <- makeCluster(length(U)); clusterExport(cl, c("raiz", "avalia")); aval <- bind_rows(parLapply(cl, U, function(u, FP) avalia(u, FP), FP = FP)); stopCluster(cl)
write.csv(FP, file.path(raiz, "saida/opcaoB_habitats.csv"), row.names = FALSE); write.csv(cons, file.path(raiz, "saida/opcaoB_pares.csv"), row.names = FALSE)
write.csv(aval, file.path(raiz, "saida/opcaoB_comparacao.csv"), row.names = FALSE); write.csv(pz, file.path(raiz, "saida/opcaoB_zonas.csv"), row.names = FALSE)
hs <- unique(FP$habitat); st <- rast(lapply(seq_len(84), function(i) { m <- (i - 1) %/% 7 + 1; z <- ZP[(i - 1) %% 7 + 1]; r <- rast(g25); v <- rep(NA, ncell(r)); s <- G[G$mes == m & G$z == z, ]; v[s$cel] <- match(s$habitat, hs); r[] <- v; names(r) <- sprintf("m%02d_z%d", m, z); r }))
writeRaster(st, file.path(raiz, "saida/opcaoB_pelagico_mes_z.tif"), overwrite = TRUE, datatype = "INT1U")
saveRDS(list(km = km, FP = FP), file.path(raiz, "intermediario/etapa3_n1b/opcaoB_modelo.rds"))
options(width = 180); cat("habitats:", length(hs), "| stand-by:", sum(!tapply(FP$testado, FP$final, any)), "| consenso:", paste(names(table(cons$consenso)), table(cons$consenso), collapse = "; "), "\n"); print(zn); print(pz |> filter(zona > 1))
print(aval |> mutate(across(where(is.double), ~ round(.x, 3)))); cat(format(Sys.time()), "FIM\n")
