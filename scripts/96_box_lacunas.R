# 96_box_lacunas.R  (decision of the PI, 06/11/2026; built on 13/11/2026 after Gate 12)
# Figure 8 / Box "Gaps and how to contribute": where the united maps lack fauna evidence, how much area that is, whether the area lies
# inside the area of applicability (AOA, stage 4), and which data would close each gap.
#  - Evidence grade per class from saida/mapa_unido_{bentos,pelagico}_classes.csv. Benthic classes that the independent test fauna of
#    Gate 12 (saida/teste_fauna_c_classes.csv) found distinct get their own grade ("distinta pela fauna de teste"), without changing the
#    approved map.
#  - Area: benthos = 25 km cells; pelagic = cell x month x standard depth layers of the monthly map (space-time share).
#  - Stand-by classes (merged by environment, untested) get their own colour and hatching.
# Outputs: saida/box_lacunas_figura.png, saida/box_lacunas_faixas.csv, saida/box_lacunas_tabela.csv
suppressPackageStartupMessages({library(dplyr); library(terra)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
fb <- read.csv(file.path(raiz, "saida/mapa_unido_bentos_classes.csv"), encoding = "UTF-8")
fp <- read.csv(file.path(raiz, "saida/mapa_unido_pelagico_classes.csv"), encoding = "UTF-8")
tc <- read.csv(file.path(raiz, "saida/teste_fauna_c_classes.csv"), encoding = "UTF-8")
aoa <- read.csv(file.path(raiz, "saida/etapa4_aoa_resumo.csv"), encoding = "UTF-8")
ub <- rast(file.path(raiz, "saida/mapa_unido_bentos_25km.tif"))[["base"]]; up <- rast(file.path(raiz, "saida/mapa_unido_pelagico_mes_z.tif"))

EV <- c("distinta pela fauna", "separada pela literatura", "fundida pela fauna", "distinta pela fauna de teste", "fundida pelo ambiente", "não testada")
COR <- c("#1b7837", "#2166ac", "#a6dba0", "#35978f", "#e08214", "#8073ac")
ROT <- c("distinta pela fauna", "separada pela literatura", "fundida pela fauna (equivalente)", "distinta pela fauna de teste (Gate 12)",
         "stand-by: fundida só pelo ambiente", "stand-by: não testada")
nova <- tc$codigo[tc$comp == "bentos" & tc$distinta_c & tc$evidencia_antes %in% c("fundida pelo ambiente", "não testada")]
fb$ev <- ifelse(fb$codigo %in% nova, "distinta pela fauna de teste", fb$evidencia); fp$ev <- fp$evidencia
fb$sb <- fb$ev %in% EV[5:6]; fp$sb <- fp$ev %in% EV[5:6]

# ---------------- area per band / zone and evidence grade ----------------
vb <- values(ub)[, 1]; cb <- table(factor(vb, fb$valor)); fb$area <- as.numeric(cb[as.character(fb$valor)])
vp <- values(up); cp <- tabulate(vp[!is.na(vp)], nbins = max(fp$valor)); fp$area <- cp[fp$valor]
ZB <- c("0–200 m", "200–2.000 m", "> 2.000 m"); ZPel <- c("epipelágico (< 200 m)", "mesopelágico (200–1.000 m)", "batipelágico (≥ 1.000 m)")
tab <- bind_rows(fb |> transmute(comp = "bentos", faixa_zona = faixa, ev, sb, area), fp |> transmute(comp = "pelagico", faixa_zona = zona, ev, sb, area))
fx <- tab |> group_by(comp, faixa_zona) |> summarise(classes = n(), classes_standby = sum(sb), area_standby = sum(area[sb]) / sum(area),
  area_fauna = sum(area[ev %in% EV[c(1, 3, 4)]]) / sum(area), area_literatura = sum(area[ev == EV[2]]) / sum(area), .groups = "drop") |>
  left_join(aoa |> select(comp, faixa_zona, area_aoa = area_dentro), by = c("comp", "faixa_zona")) |>
  mutate(ordem = match(faixa_zona, c(ZB, ZPel))) |> arrange(comp, ordem) |> select(-ordem)
write.csv(fx, file.path(raiz, "saida/box_lacunas_faixas.csv"), row.names = FALSE)
prop <- tab |> group_by(comp, faixa_zona, ev) |> summarise(a = sum(area), .groups = "drop") |> group_by(comp, faixa_zona) |> mutate(a = a / sum(a)) |> ungroup()

# ---------------- gap table: where, how large, what closes it ----------------
f <- function(cp, z) fx[fx$comp == cp & fx$faixa_zona == z, ]
pct <- function(x) sprintf("%.0f%%", 100 * x)
g <- list(
  c("Bentos profundo (> 2.000 m)", sprintf("%d de %d classes em stand-by; %s da área em stand-by; %s da área dentro da AOA", f("bentos", ZB[3])$classes_standby, f("bentos", ZB[3])$classes, pct(f("bentos", ZB[3])$area_standby), pct(f("bentos", ZB[3])$area_aoa)),
    "Amostras de fauna de fundo abissal e batial com profundidade, aparelho (trenó epibentônico, ROV, box corer) e lista completa por evento", "Publicar no OBIS como evento (eventID, profundidade mínima e máxima, samplingProtocol); 10 eventos numa célula de 25 km formam uma unidade"),
  c("Bentos do talude (200–2.000 m)", sprintf("%d de %d classes em stand-by; %s da área em stand-by; %s dentro da AOA", f("bentos", ZB[2])$classes_standby, f("bentos", ZB[2])$classes, pct(f("bentos", ZB[2])$area_standby), pct(f("bentos", ZB[2])$area_aoa)),
    "Arrastos de levantamento com ausências e esforço; grupos ainda fora do método (Gastropoda, Amphipoda, Hydrozoa); a fauna de teste do Gate 12 já resolveu 17 classes", "Conjuntos com ausências explícitas e esforço; cada grupo novo entra como fauna de teste sem refazer o mapa"),
  c("Pelágico mesopelágico (200–1.000 m)", sprintf("%d de %d classes em stand-by; %s do espaço-tempo em stand-by; %s dentro da AOA", f("pelagico", ZPel[2])$classes_standby, f("pelagico", ZPel[2])$classes, pct(f("pelagico", ZPel[2])$area_standby), pct(f("pelagico", ZPel[2])$area_aoa)),
    "Redes estratificadas (MOCNESS, Multinet) com cada lance dentro de uma zona, repetidas nas estações do ano; DNA ambiental com resolução de gênero", "Profundidade mínima e máxima de cada lance; hoje 43% dos lances abaixo de 200 m cruzam o limite e não servem para testar zonas"),
  c("Pelágico batipelágico (≥ 1.000 m)", sprintf("%d de %d classes em stand-by; %s do espaço-tempo em stand-by; %s dentro da AOA", f("pelagico", ZPel[3])$classes_standby, f("pelagico", ZPel[3])$classes, pct(f("pelagico", ZPel[3])$area_standby), pct(f("pelagico", ZPel[3])$area_aoa)),
    "Qualquer amostra de fauna abaixo de 1.000 m com profundidade do lance; bases regionais (DEEPEND, BATHYPELAGIC, MAFIA) ainda fora do OBIS", "Depositar no OBIS ou indicar a base; licença aberta (CC0 ou CC-BY)"),
  c("Pelágico epipelágico (< 200 m)", sprintf("%d de %d classes em stand-by; %s do espaço-tempo em stand-by; %s dentro da AOA", f("pelagico", ZPel[1])$classes_standby, f("pelagico", ZPel[1])$classes, pct(f("pelagico", ZPel[1])$area_standby), pct(f("pelagico", ZPel[1])$area_aoa)),
    "Séries repetidas por trimestre fora do Atlântico Norte (as classes sazonais pedem as quatro estações)", "Séries temporais no OBIS com data completa (o mês de date_start)"),
  c("Limites verticais de 200 e 1.000 m", "Só o micronécton apoia 200 m (ω² 0,147); nenhuma fauna apoia 1.000 m; procariotos sem corte preferencial (ω² ≤ 0,019)",
    "Perfis verticais completos na mesma estação (várias profundidades, mesmo aparelho) e DNA ambiental identificado até gênero; Syndiniales ficou de fora por falta de gênero", "Amostra com profundidade de coleta e identificação até gênero"),
  c("Limiares ambientais", "Limiares do GLORYS12 confirmados globalmente (24), regionais (16), não confirmados (4) e não testáveis (16)",
    "Limiares publicados com variável, valor e região", "Entram na tabela de limiares publicados e são testados pelo mesmo caminho (etapa 3)"))
gap <- setNames(as.data.frame(do.call(rbind, g)), c("lacuna", "tamanho", "dado_que_resolve", "como_contribuir"))
write.csv(gap, file.path(raiz, "saida/box_lacunas_tabela.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ---------------- figure ----------------
hatch <- function(r, sbv) { v <- values(r)[, 1]; rc <- rowColFromCell(r, seq_len(ncell(r))); k0 <- which(!is.na(v)); k <- k0[sbv[v[k0]] & (rc[k0, 1] + rc[k0, 2]) %% 7 < 2]
  h <- rast(r); hv <- rep(NA, ncell(r)); hv[k] <- 1; h[] <- hv; h }
mapa <- function(r, f, titulo) { ix <- match(f$ev, EV); v <- values(r)[, 1]; e <- rast(r); e[] <- ix[match(v, f$valor)]
  plot(e, col = COR, breaks = 0:6 + 0.5, legend = FALSE, axes = FALSE, box = FALSE, maxcell = ncell(r), main = titulo, cex.main = 1.1)
  plot(hatch(r, setNames(f$sb, f$valor)[as.character(seq_len(max(f$valor)))]), col = "#3a3a3a", breaks = c(0.5, 1.5), legend = FALSE, axes = FALSE, box = FALSE, add = TRUE, maxcell = ncell(r)) }
png(file.path(raiz, "saida/box_lacunas_figura.png"), 2400, 2000, res = 170, bg = "white")
layout(matrix(c(1, 1, 1, 0, 2, 0, 3, 3, 3, 4, 4, 4), 4, byrow = TRUE), widths = c(1, 3, 1), heights = c(1000, 640, 820, 150))
par(mar = c(0.3, 0.5, 2, 0.5))
mapa(ub, fb, "(a) Bentos, nível base: grau de evidência de cada classe")
mapa(up[["m07_z500"]], fp, "(b) Pelágico, julho, 500 m (mesopelágico)")
par(mar = c(4, 12, 2.5, 7)); grp <- c(paste("Bentos", ZB), paste("Pelágico", c("epi", "meso", "bati")))
M <- sapply(seq_along(grp), function(i) { cp <- if (i <= 3) "bentos" else "pelagico"; z <- c(ZB, ZPel)[i]; x <- prop |> filter(comp == cp, faixa_zona == z); setNames(x$a, x$ev)[EV] }); M[is.na(M)] <- 0
bp <- barplot(M, horiz = TRUE, col = COR, names.arg = grp, las = 1, xlim = c(0, 1), border = NA, xlab = "fração da área (bentos) ou do espaço-tempo (pelágico)", cex.names = 0.95,
  main = "(c) Área por grau de evidência e área de aplicabilidade (AOA)", cex.main = 1.1)
hb <- 0.5 * diff(bp[1:2]) / 1.2   # half bar height (barplot default space 0.2)
for (j in seq_len(ncol(M))) { x1 <- cumsum(M[, j]); x0 <- x1 - M[, j]
  for (i in 5:6) if (M[i, j] > 0) rect(x0[i], bp[j] - hb, x1[i], bp[j] + hb, density = 18, angle = 45, col = "#3a3a3a", border = NA) }
ao <- c(fx$area_aoa[fx$comp == "bentos"][match(ZB, fx$faixa_zona[fx$comp == "bentos"])], fx$area_aoa[fx$comp == "pelagico"][match(ZPel, fx$faixa_zona[fx$comp == "pelagico"])])
points(ao, bp, pch = 23, bg = "white", cex = 1.6); text(1.01, bp, sprintf("AOA %.0f%%", 100 * ao), xpd = TRUE, adj = 0, cex = 0.9)
par(mar = c(0, 1, 0, 1)); plot.new()
legend("center", legend = c(ROT, "losango = fração dentro da AOA"), fill = c(COR, NA), border = NA,
  pch = c(rep(NA, 6), 23), pt.bg = "white", ncol = 4, cex = 0.9, bty = "n", title = "achurado = stand-by (lacuna de conhecimento, não ausência de habitat)")
dev.off()
options(width = 200); print(as.data.frame(fx |> mutate(across(where(is.double), ~ round(.x, 3))))); cat(format(Sys.time()), "FIM\n")
