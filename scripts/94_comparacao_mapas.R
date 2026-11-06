# 94_comparacao_mapas.R  (side-by-side maps for the validation report, in one common style without legends)
# Four panels per compartment and layer, same grid, projection and size:
#   original classification (EBU regions; EMU classes at the depth) | classification tested by the fauna (85_sintese.R) |
#   new habitats from the gradients (HB, 89b; HP, updated with level 2 in 93) | united map, base level (91, 91b).
# Colours: the original map takes the colour of its class in the tested map, with tones when classes were merged, so the two read
# together; the united map keeps the hue of its new habitat. Hatched = stand-by / no fauna evidence; cross-hatched = separated by
# the literature (criterion 9). The original map has no hatching: it was not evaluated.
# Pelagic: surface and 1,000 m (the depths of 85_sintese.R), January and July for the two dynamic panels.
# Outputs: saida/comparacao/*.png
suppressPackageStartupMessages({library(dplyr); library(terra); library(arrow)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); source(file.path(raiz, "scripts/funcoes_grade.R")); source(file.path(raiz, "scripts/funcoes_mapa.R"))
dir_s <- file.path(raiz, "saida/comparacao"); dir.create(dir_s, showWarnings = FALSE)
g25 <- ease2_grade("25km"); ordena <- function(h) h[order(as.integer(sub("^H[BP]", "", h)))]
painel <- function(v, cores, titulo, arq, hach = NULL, cruz = NULL) {   # v: values 1..length(cores) per cell
  r <- rast(g25); r[] <- v; ii <- which(!is.na(v)); rc <- rowColFromCell(r, ii); d1 <- (rc[, 1] + rc[, 2]) %% 7 < 2; d2 <- (rc[, 1] - rc[, 2]) %% 7 < 2
  k <- c(if (!is.null(hach)) ii[hach[v[ii]] & d1], if (!is.null(cruz)) ii[cruz[v[ii]] & (d1 | d2)])
  h <- rast(r); hv <- rep(NA_integer_, ncell(r)); hv[k] <- 1L; h[] <- hv
  png(file.path(dir_s, paste0(arq, ".png")), 1400, 640, res = 130, bg = "white"); par(mar = c(0.2, 0.2, 1.8, 0.2))
  plot(r, col = cores, breaks = c(0.5, seq_along(cores) + 0.5), legend = FALSE, axes = FALSE, box = FALSE, maxcell = ncell(r), main = titulo, cex.main = 1.05, mar = NA)
  if (length(k)) plot(h, col = COR_ACH, breaks = c(0.5, 1.5), legend = FALSE, axes = FALSE, box = FALSE, add = TRUE, maxcell = ncell(r), mar = NA)
  dev.off(); invisible(arq) }
tons <- function(hab, cls, leg) {   # colour of each original habitat: its tested class, with tones inside merged classes
  cor <- leg$cor[match(cls, leg$classe_final)]; cor[is.na(cor)] <- "#BBBBBB"
  for (c0 in unique(cls[duplicated(cls)])) { i <- which(cls == c0); f <- seq(-0.4, 0.4, length.out = length(i))
    cor[i] <- vapply(seq_along(i), function(j) { b <- col2rgb(cor[i[j]]) / 255; x <- f[j]; a <- if (x < 0) b * (1 + x) else b + (1 - b) * x; rgb(a[1], a[2], a[3]) }, "") }
  cor }
sint <- read.csv(file.path(raiz, "saida/sintese_habitats_n1.csv"))
testado <- function(camada, titulo, arq) { leg <- read.csv(file.path(raiz, sprintf("saida/habitats_final_%s_legenda.csv", camada)))
  v <- values(rast(file.path(raiz, sprintf("saida/habitats_final_%s.tif", camada))))[, 1]
  x <- match(v, leg$valor); painel(x, leg$cor, sprintf("%s: %d classes", titulo, length(unique(x[!is.na(x)]))), arq, hach = leg$standby); leg }

# ---------------- benthos ----------------
e25 <- read_parquet(file.path(raiz, "camadas/ebu/ebu_celula_25km.parquet")); p <- do.call(rbind, strsplit(e25$cell_id, "_"))
cel <- cellFromRowCol(g25, as.integer(p[, 2]), as.integer(p[, 3]))
leg <- testado("ebu_regiao_25km", "EBU testada pela fauna", "bentos_2_testado")
s <- sint |> filter(nivel == "EBU região"); hs <- sort(unique(na.omit(e25$ebu_regiao))); v <- rep(NA_integer_, ncell(g25)); v[cel] <- match(e25$ebu_regiao, hs)
painel(v, tons(hs, s$classe_final[match(hs, s$habitat)], leg), sprintf("EBU original: %d regiões", length(hs)), "bentos_1_original")
lh <- read.csv(file.path(raiz, "saida/habitats_novos_legenda.csv")) |> filter(arquivo == "bentos_25km")
x <- values(rast(file.path(raiz, "saida/habitats_novos_bentos_25km.tif")))[, 1]
painel(x, lh$cor, sprintf("Habitats novos (bottom-up): %d habitats", nrow(lh)), "bentos_3_novos", hach = lh$standby)
fin <- read.csv(file.path(raiz, "saida/mapa_unido_bentos_classes.csv")); x <- values(rast(file.path(raiz, "saida/mapa_unido_bentos_25km.tif"))[["base"]])[, 1]
painel(x, fin$cor, sprintf("Mapa unido, nível base: %d classes", nrow(fin)), "bentos_4_unido", hach = fin$standby, cruz = fin$evidencia == "separada pela literatura")

# ---------------- pelagic ----------------
emu <- rast(file.path(raiz, "camadas/emu/emu_25km.tif")); s <- sint |> filter(nivel == "EMU")
fus <- read.csv(file.path(raiz, "saida/pelagico_n2_fusoes.csv")); nf <- tapply(fus$hp_n1, fus$final_n1n2, function(x) ordena(unique(x))[1]); fus$hp <- as.character(nf[as.character(fus$final_n1n2)])
hp_de <- tapply(fus$hp, fus$hp_n1, function(x) x[1]); antigo <- ordena(unique(fus$hp_n1)); hsP <- ordena(unique(fus$hp)); hp_sb <- tapply(fus$standby_n1n2, fus$hp, function(x) x[1])[hsP]
lp <- read.csv(file.path(raiz, "saida/habitats_novos_legenda.csv")) |> filter(grepl("^pelagico", arquivo)) |> distinct(habitat, cor); corHP <- setNames(lp$cor, lp$habitat)[hsP]
hpt <- rast(file.path(raiz, "saida/habitats_novos_pelagico_mes_z.tif"))
mpf <- read.csv(file.path(raiz, "saida/mapa_unido_pelagico_classes.csv")); mpt <- rast(file.path(raiz, "saida/mapa_unido_pelagico_mes_z.tif"))
for (z in c(0, 1000)) { nm <- if (z == 0) "superficie" else "1000m"; rot <- if (z == 0) "superfície" else "1.000 m"
  leg <- testado(sprintf("emu_%s_25km", nm), sprintf("EMU testada pela fauna, %s", rot), sprintf("pelagico_%s_2_testado", nm))
  ev <- values(emu[[sprintf("prof_%dm", z)]])[, 1]; hs <- sort(unique(ev[!is.na(ev)]))
  painel(match(ev, hs), tons(hs, s$classe_final[match(as.character(hs), s$habitat)], leg), sprintf("EMU original, %s: %d classes", rot, length(hs)), sprintf("pelagico_%s_1_original", nm))
  for (m in c(1, 7)) { mes <- c("janeiro", "", "", "", "", "", "julho")[m]; cam <- sprintf("m%02d_z%d", m, z)
    h <- values(hpt[[cam]])[, 1]; x <- match(unname(hp_de[antigo[h]]), hsP)
    painel(x, corHP, sprintf("Habitats novos (bottom-up), %s, %s: %d habitats", rot, mes, length(unique(x[!is.na(x)]))), sprintf("pelagico_%s_3_novos_m%02d", nm, m), hach = hp_sb)
    x <- values(mpt[[cam]])[, 1]
    painel(x, mpf$cor, sprintf("Mapa unido, %s, %s: %d classes", rot, mes, length(unique(x[!is.na(x)]))), sprintf("pelagico_%s_4_unido_m%02d", nm, m), hach = mpf$standby) } }
cat(format(Sys.time()), "FIM\n")
