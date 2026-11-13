# 95_teste_fauna_c.R  (Gate 12, 13/11/2026, items 1 and 2: test fauna for the approved maps)
# The groups of 21c (benthos: demersal fishes, benthic Decapoda; pelagic: pelagic Cyclopoida, Appendicularia, Pteropoda) did not take part
# in building any map, so they are an independent test. Nothing approved is rebuilt.
#  1 acceptance, as in 91/91b: omega2 (all data) and equal-effort omega2 and power of each classification on the units of each new group:
#    benthos  = united map (base), level III, new habitats (HB), EBU region, EBU group (the unit's cell);
#    pelagic  = united map (base), level III, HP (93b names), EMU. Pelagic units are cell x EMU x quarter with a depth: the class is read
#               from the monthly rasters at the middle month of the quarter and the nearest standard depth (0, 50, 100, 200, 500, 1000, 2000 m).
#  2 pairs (pares_classes: omega2 >= 0.01 and power >= 0.8): HB and HP pairs, and pairs of base classes of the united maps in the same
#    depth band (benthos) or vertical zone (pelagic). A base class whose pairs separate in a new group gains independent fauna evidence;
#    an untested class that now has a test is counted.
# Outputs: saida/teste_fauna_c_{avaliacao,pares,classes}.csv
suppressPackageStartupMessages({library(dplyr); library(terra); library(arrow); library(parallel)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); source(file.path(raiz, "scripts/funcoes_grade.R"))
g25 <- ease2_grade("25km"); ZP <- c(0, 50, 100, 200, 500, 1000, 2000); zona <- function(z) findInterval(z, c(200, 1000)) + 1
ordena <- function(h) h[order(as.integer(sub("^H[BP]", "", h)))]
cel_de <- function(id) { p <- do.call(rbind, strsplit(id, "_")); cellFromRowCol(g25, as.integer(p[, 2]), as.integer(p[, 3])) }
dir_c <- file.path(raiz, "intermediario/etapa3_n1c")

# ---------------- class of each unit ----------------
fb <- read.csv(file.path(raiz, "saida/mapa_unido_bentos_classes.csv"), encoding = "UTF-8"); ub <- rast(file.path(raiz, "saida/mapa_unido_bentos_25km.tif"))
hbv <- values(rast(file.path(raiz, "saida/habitats_novos_bentos_25km.tif")))[, 1]
lh <- read.csv(file.path(raiz, "saida/habitats_novos_legenda.csv")); hb_nome <- lh$habitat[lh$arquivo == "bentos_25km"]
bv <- values(ub[["base"]])[, 1]; fxb <- setNames(fb$faixa, fb$codigo)
fp <- read.csv(file.path(raiz, "saida/mapa_unido_pelagico_classes.csv"), encoding = "UTF-8"); up <- rast(file.path(raiz, "saida/mapa_unido_pelagico_mes_z.tif"))
hpt <- rast(file.path(raiz, "saida/habitats_novos_pelagico_mes_z.tif"))
fus <- read.csv(file.path(raiz, "saida/pelagico_b_fusoes.csv")); nf <- tapply(fus$hp_n1, fus$final, function(x) ordena(unique(x))[1])
hp_de <- setNames(as.character(nf[as.character(fus$final)]), fus$hp_n1); hs0 <- ordena(unique(fus$hp_n1))
le <- function(st, cel, lay) { v <- rep(NA_integer_, length(cel)); for (l in unique(lay)) { i <- which(lay == l); v[i] <- values(st[[l]])[cel[i], 1] }; v }
arqs <- list.files(dir_c, pattern = "^(bentos|pelagico)_.*[.]rds$", full.names = TRUE)
U <- lapply(arqs, function(a) { r <- readRDS(a); bent <- r$comp == "bentos"
  if (bent) { cel <- cel_de(r$meta$unidade); b <- bv[cel]; cod <- fb$codigo[match(b, fb$valor)]
    cls <- list("mapa unido (base)" = cod, "nível III" = fb$nivel3[match(b, fb$valor)], "habitats novos (HB)" = hb_nome[hbv[cel]],
                "EBU região" = r$meta$ebu_regiao, "EBU grupo" = r$meta$ebu_grupo); estrato <- fxb[cod]
  } else { cel <- cel_de(r$meta$cell_id); m <- c(T1 = 2, T2 = 5, T3 = 8, T4 = 11)[r$meta$tri]; z <- ZP[vapply(r$E$z, function(x) which.min(abs(ZP - x)), 1L)]
    lay <- sprintf("m%02d_z%d", m, z); b <- le(up, cel, lay); cod <- fp$codigo[match(b, fp$valor)]; h <- hp_de[hs0[le(hpt, cel, lay)]]
    cls <- list("mapa unido (base)" = cod, "nível III" = fp$nivel3[match(b, fp$valor)], "habitats novos (HP)" = unname(h), "EMU" = as.character(r$meta$emu))
    estrato <- as.character(zona(r$E$z)) }
  list(arq = a, comp = r$comp, grupo = r$grupo, cls = cls, estrato = unname(estrato)) })
cat(format(Sys.time()), "grupos:", paste(sapply(U, `[[`, "grupo"), collapse = "; "), "\n")

# ---------------- 1 acceptance ----------------
avalia <- function(u) { suppressPackageStartupMessages({library(vegan); library(dplyr)}); source(file.path(raiz, "scripts/funcoes_permanova.R"))
  r <- readRDS(u$arq); Dm <- as.matrix(sqrt(betadiver(r$X, "sim"))); bl <- r$meta$bloco
  bind_rows(lapply(names(u$cls), function(nm) { g <- u$cls[[nm]]; k <- !is.na(g) & g %in% names(which(table(g) >= 5)); if (length(unique(g[k])) < 2) return(NULL)
    set.seed(9); x <- permanova_rapida(as.dist(Dm[k, k]), g[k], bl[k], 199)
    set.seed(10); e <- tryCatch(permanova_esforco_igual(Dm[k, k], g[k], bl[k], n_classe = 10, B = 100), error = function(e) NULL)
    data.frame(comp = u$comp, grupo = u$grupo, classificacao = nm, classes = length(unique(g[k])), unidades = sum(k), omega2 = x$omega2,
      omega2_eq = if (!is.null(e)) e$omega2_eq else NA, poder_eq = if (!is.null(e)) e$poder_eq else NA) })) }
cl <- makeCluster(length(U)); clusterExport(cl, "raiz"); aval <- bind_rows(parLapply(cl, U, avalia)); stopCluster(cl)
cat(format(Sys.time()), "aceitacao feita\n")

# ---------------- 2 pairs ----------------
jobs <- bind_rows(lapply(seq_along(U), function(k) { u <- U[[k]]; out <- list()
  for (nm in c(grep("^habitats novos", names(u$cls), value = TRUE), "mapa unido (base)")) { g <- u$cls[[nm]]; n <- table(g[!is.na(g)]); cl5 <- sort(names(n)[n >= 5])
    if (length(cl5) < 2) next; pr <- t(combn(cl5, 2))
    if (nm == "mapa unido (base)") { ea <- tapply(u$estrato, g, function(x) names(which.max(table(x)))); pr <- pr[ea[pr[, 1]] == ea[pr[, 2]], , drop = FALSE] }
    if (nrow(pr)) out[[nm]] <- data.frame(k = k, cls = nm, a = pr[, 1], b = pr[, 2]) }
  bind_rows(out) }))
jobs$lote <- ave(seq_len(nrow(jobs)), jobs$k, jobs$cls, FUN = function(i) ceiling(seq_along(i) / 25))
cat("pares a testar:", nrow(jobs), "\n")
testa <- function(j) { suppressPackageStartupMessages({library(vegan); library(dplyr)}); source(file.path(raiz, "scripts/funcoes_permanova.R"))
  u <- U[[j$k[1]]]; r <- readRDS(u$arq); g <- u$cls[[j$cls[1]]]; ok <- !is.na(g); Dm <- as.matrix(sqrt(betadiver(r$X[ok, , drop = FALSE], "sim"))); set.seed(7600 + j$k[1])
  x <- pares_classes(Dm, g[ok], r$meta$bloco[ok], un_min = 5, n_classe = 10, B = 100, pares = as.matrix(j[, c("a", "b")]))
  if (!is.null(x)) cbind(comp = u$comp, grupo = u$grupo, classificacao = j$cls[1], x) }
cl <- makeCluster(14); clusterExport(cl, c("raiz", "U")); pares <- bind_rows(parLapplyLB(cl, split(jobs, list(jobs$k, jobs$cls, jobs$lote), drop = TRUE), testa)); stopCluster(cl)
pares <- pares |> filter(!is.na(separa)); cat(format(Sys.time()), "pares testados:", nrow(pares), "| separam:", sum(pares$separa), "\n")

# ---------------- base classes: what the new fauna adds ----------------
pb <- pares |> filter(classificacao == "mapa unido (base)")
ev <- bind_rows(lapply(c("bentos", "pelagico"), function(cp) { f <- if (cp == "bentos") fb else fp; p <- pb |> filter(comp == cp)
  data.frame(comp = cp, codigo = f$codigo, evidencia_antes = f$evidencia, testada_c = f$codigo %in% c(p$classe_a, p$classe_b),
             distinta_c = f$codigo %in% c(p$classe_a[p$separa], p$classe_b[p$separa])) }))
write.csv(aval, file.path(raiz, "saida/teste_fauna_c_avaliacao.csv"), row.names = FALSE); write.csv(pares, file.path(raiz, "saida/teste_fauna_c_pares.csv"), row.names = FALSE)
write.csv(ev, file.path(raiz, "saida/teste_fauna_c_classes.csv"), row.names = FALSE)
options(width = 220); print(aval |> mutate(across(where(is.double), ~ round(.x, 3))) |> as.data.frame())
print(pares |> group_by(comp, grupo, classificacao) |> summarise(pares = n(), separam = sum(separa), .groups = "drop") |> as.data.frame())
print(ev |> group_by(comp, evidencia_antes) |> summarise(classes = n(), testadas_c = sum(testada_c), distintas_c = sum(distinta_c), .groups = "drop") |> as.data.frame())
cat(format(Sys.time()), "FIM\n")
