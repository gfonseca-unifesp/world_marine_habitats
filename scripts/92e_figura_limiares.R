# 92e_figura_limiares.R  (GUIA_LIMIARES_PUBLICADOS.md, step L5: new panel of Figure 4)
# One panel per variable, one row per basin: robust global thresholds of the model (grey bands and points in the row "modelo"),
# local thresholds of the basin GF (ticks, 88c), published thresholds (line = interval; circle at the midpoint, filled when confirmed
# by enrichment; colour = Gate 10 use). Colours are Okabe-Ito (colourblind-safe); identity is also in the legend.
# Output: saida/figura_limiares.png
suppressPackageStartupMessages(library(dplyr))
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
USO <- c(validacao = "#0072B2", apoio = "#D55E00", relato = "#7A7A7A"); NOME_USO <- c(validacao = "validação", apoio = "apoio", relato = "relato")
INK <- "#222222"; MUTED <- "#666666"; BANDA <- "#DDDDDD"; GRADE <- "#EEEEEE"
VARS <- list(prof = list(cp = "bentos", t = "Profundidade do fundo (m), bentos", log = ""), tfundo_med = list(cp = "bentos", t = "Temperatura de fundo (°C), bentos", log = ""),
  z = list(cp = "pelagico", t = "Profundidade z (m, escala log), pelágico", log = "x"), sst = list(cp = "pelagico", t = "SST (°C), pelágico", log = ""),
  t_z = list(cp = "pelagico", t = "Temperatura na profundidade z (°C), pelágico", log = ""))
en <- read.csv(file.path(raiz, "saida/limiares_enriquecimento.csv"), encoding = "UTF-8") |> group_by(id, compartimento, variavel_modelo) |>
  summarise(limiar_min = first(limiar_min), limiar_max = first(limiar_max), bacia = first(bacia), uso = first(uso), confirmado = any(confirmado_enriq), .groups = "drop") |>
  filter(!is.na(limiar_min)) |> mutate(bacia = ifelse(is.na(bacia) | bacia == "", "global", bacia))
rob <- read.csv(file.path(raiz, "saida/etapa3_n1_limiares_robustos.csv"))
bac <- read.csv(file.path(raiz, "camadas/bacias_legenda.csv"))
loc <- read.csv(file.path(raiz, "saida/etapa3_n1_transf_limiares_locais.csv")) |> mutate(bacia_nome = bac$bacia[match(bacia, bac$codigo)])
painel <- function(v, o) { p <- en |> filter(variavel_modelo == v, compartimento == o$cp); r <- rob |> filter(preditor == v, grepl(paste0("\\(", o$cp, "\\)"), lista))
  l <- loc |> filter(preditor == v, comp == o$cp, bacia_nome %in% p$bacia)
  bs <- sort(setdiff(unique(p$bacia), "global")); linhas <- c("modelo", if ("global" %in% p$bacia) "global", bs); y <- rev(seq_along(linhas)); names(y) <- linhas
  xs <- c(p$limiar_min, p$limiar_max, r$de, r$ate); if (o$log == "x") xs <- xs[xs > 0]; xl <- range(xs, na.rm = TRUE)
  par(mar = c(3, 8.5, 2.2, 1), mgp = c(2, 0.5, 0), tcl = -0.2)
  plot(NA, xlim = xl, ylim = c(0.4, length(linhas) + 0.6), log = o$log, axes = FALSE, xlab = "", ylab = "")
  rect(pmax(r$de, xl[1]), 0.4, r$ate, length(linhas) + 0.6, col = BANDA, border = NA)
  abline(h = y, col = GRADE); axis(1, col = MUTED, col.axis = MUTED, cex.axis = 0.8); axis(2, y, linhas, las = 1, tick = FALSE, col.axis = INK, cex.axis = 0.8)
  title(main = o$t, adj = 0, cex.main = 0.95, font.main = 1, col.main = INK)
  points(r$limiar, rep(y["modelo"], nrow(r)), pch = 21, bg = INK, col = "white", cex = 1.3, lwd = 1)
  if (nrow(l)) segments(l$limiar, y[l$bacia_nome] - 0.34, l$limiar, y[l$bacia_nome] - 0.16, col = INK, lwd = 1)
  for (b in unique(p$bacia)) { q <- p[p$bacia == b, ] |> arrange(limiar_min); n <- nrow(q); dy <- if (n == 1) 0.1 else seq(0.32, -0.1, length.out = n)   # stacked in the row
    for (k in seq_len(n)) { yy <- y[b] + dy[k]; cor <- USO[q$uso[k]]; a <- q$limiar_min[k]; d <- q$limiar_max[k]
      if (d > a) segments(a, yy, d, yy, col = cor, lwd = 2, lend = 1)
      m <- if (o$log == "x") sqrt(a * d) else (a + d) / 2
      points(m, yy, pch = 21, bg = if (q$confirmado[k]) cor else "white", col = cor, cex = 1.1, lwd = 1.5) } } }
png(file.path(raiz, "saida/figura_limiares.png"), 2200, 1400, res = 170, bg = "white"); layout(matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow = TRUE), heights = c(1.25, 1.1, 1))
for (v in names(VARS)) painel(v, VARS[[v]])
par(mar = c(0, 1, 0, 1)); plot.new()
legend("left", bty = "n", cex = 0.95, text.col = INK, y.intersp = 1.25,
  legend = c(paste("publicado,", NOME_USO), "publicado, confirmado por enriquecimento (preenchido)", "modelo global: limiar robusto (faixa cinza = intervalo)", "modelo local da bacia (traço)"),
  pch = c(21, 21, 21, 21, 21, NA), pt.bg = c("white", "white", "white", MUTED, INK, NA), col = c(USO, MUTED, INK, INK), lty = c(NA, NA, NA, NA, NA, 1), pt.cex = 1.2)
dev.off(); cat("figura gravada\n")
