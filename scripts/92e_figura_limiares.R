# 92e_figura_limiares.R  (GUIA_LIMIARES_PUBLICADOS.md, step L5: new panel of Figure 4)
# One panel per variable, one row per study (basin in the row label): robust global thresholds of the model (grey bands and points in
# the row "modelo"), local thresholds of the basin GF (ticks, 88c), published thresholds (line = interval; circle at the midpoint,
# filled when confirmed by enrichment; colour = Gate 10 use; label = threshold id, the key to the table of saida/figura_limiares_tabela.csv).
# Colours are Okabe-Ito (colourblind-safe); identity is also in the legend.
# Outputs: saida/figura_limiares.png, saida/figura_limiares_tabela.csv (id, panel and row of each plotted threshold)
suppressPackageStartupMessages(library(dplyr))
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
USO <- c(validacao = "#0072B2", apoio = "#D55E00", relato = "#7A7A7A"); NOME_USO <- c(validacao = "validação", apoio = "apoio", relato = "relato")
INK <- "#222222"; MUTED <- "#666666"; BANDA <- "#DDDDDD"; GRADE <- "#EEEEEE"
ABREV <- c("Atlântico Norte" = "Atl. N.", "Atlântico Sul" = "Atl. S.", "Pacífico Norte" = "Pac. N.", "Pacífico Sul" = "Pac. S.", "Índico" = "Índico", "Ártico" = "Ártico", "Austral" = "Austral", "global" = "global")
VARS <- list(prof = list(cp = "bentos", t = "Profundidade do fundo (m), bentos", log = ""), z = list(cp = "pelagico", t = "Profundidade z (m, escala log), pelágico", log = "x"),
  tfundo_med = list(cp = "bentos", t = "Temperatura de fundo (°C), bentos", log = ""), sst = list(cp = "pelagico", t = "SST (°C), pelágico", log = ""),
  t_z = list(cp = "pelagico", t = "Temperatura na profundidade z (°C), pelágico", log = ""))
en <- read.csv(file.path(raiz, "saida/limiares_enriquecimento.csv"), encoding = "UTF-8") |> group_by(id, compartimento, variavel_modelo) |>
  summarise(limiar_min = first(limiar_min), limiar_max = first(limiar_max), bacia = first(bacia), uso = first(uso), confirmado = any(confirmado_enriq), .groups = "drop") |>
  filter(!is.na(limiar_min)) |> mutate(bacia = ifelse(is.na(bacia) | bacia == "", "global", bacia)) |>
  left_join(read.csv(file.path(raiz, "saida/limiares_publicados_padronizados.csv"), encoding = "UTF-8") |> select(id, referencia), by = "id")
rob <- read.csv(file.path(raiz, "saida/etapa3_n1_limiares_robustos.csv"))
bac <- read.csv(file.path(raiz, "camadas/bacias_legenda.csv"))
loc <- read.csv(file.path(raiz, "saida/etapa3_n1_transf_limiares_locais.csv")) |> mutate(bacia_nome = bac$bacia[match(bacia, bac$codigo)])
linhas_de <- function(v, o) en |> filter(variavel_modelo == v, compartimento == o$cp) |> distinct(bacia, referencia) |> arrange(bacia != "global", bacia, referencia)
tabela <- list()
painel <- function(v, o) { p <- en |> filter(variavel_modelo == v, compartimento == o$cp); r <- rob |> filter(preditor == v, grepl(paste0("\\(", o$cp, "\\)"), lista))
  li <- linhas_de(v, o); rotulos <- c("modelo", sprintf("%s (%s)", li$referencia, ABREV[li$bacia])); n <- length(rotulos); y <- rev(seq_len(n))
  xs <- c(p$limiar_min, p$limiar_max, r$de, r$ate); if (o$log == "x") xs <- xs[xs > 0]; xl <- range(xs, na.rm = TRUE); esc <- if (o$log == "x") log10 else identity
  par(mar = c(2.6, 13, 2, 1), mgp = c(2, 0.5, 0), tcl = -0.2)
  plot(NA, xlim = xl, ylim = c(0.4, n + 0.6), log = o$log, axes = FALSE, xlab = "", ylab = "")
  rect(pmax(r$de, xl[1]), 0.4, r$ate, n + 0.6, col = BANDA, border = NA)
  abline(h = y, col = GRADE); axis(1, col = MUTED, col.axis = MUTED, cex.axis = 0.8); axis(2, y, rotulos, las = 1, tick = FALSE, col.axis = INK, cex.axis = 0.72)
  title(main = o$t, adj = 0, cex.main = 0.95, font.main = 1, col.main = INK)
  points(r$limiar, rep(y[1], nrow(r)), pch = 21, bg = INK, col = "white", cex = 1.3, lwd = 1)
  for (k in seq_len(nrow(li))) { yy <- y[k + 1]; q <- p |> filter(referencia == li$referencia[k], bacia == li$bacia[k]) |> arrange(limiar_min)
    l <- loc |> filter(preditor == v, comp == o$cp, bacia_nome == li$bacia[k]); if (nrow(l)) segments(l$limiar, yy - 0.42, l$limiar, yy - 0.28, col = MUTED, lwd = 0.8)
    meio <- if (o$log == "x") sqrt(q$limiar_min * q$limiar_max) else (q$limiar_min + q$limiar_max) / 2
    desl <- 0.26; for (j in seq_len(nrow(q))) { cor <- USO[q$uso[j]]; a <- q$limiar_min[j]; d <- q$limiar_max[j]
      if (d > a) segments(a, yy, d, yy, col = cor, lwd = 2, lend = 1)
      points(meio[j], yy, pch = 21, bg = if (q$confirmado[j]) cor else "white", col = cor, cex = 1.05, lwd = 1.5)
      perto <- j > 1 && abs(esc(meio[j]) - esc(meio[j - 1])) < 0.07 * diff(esc(xl)); desl <- if (perto && desl == 0.26) -0.27 else 0.26   # alternate labels of close thresholds
      text(meio[j], yy + desl, q$id[j], cex = 0.5, col = INK, xpd = NA) } }
  tabela[[v]] <<- p |> left_join(li |> mutate(linha = row_number()), by = c("bacia", "referencia")) |> transmute(id, painel = o$t, ordem_painel = match(v, names(VARS)), linha) }
alt <- sapply(names(VARS), function(v) nrow(linhas_de(v, VARS[[v]])) + 1)
png(file.path(raiz, "saida/figura_limiares.png"), 2500, 1900, res = 170, bg = "white")
layout(matrix(1:6, 3, 2, byrow = TRUE), heights = c(max(alt[1:2]) + 1.5, max(alt[3:4]) + 2.2, alt[5] + 2.2))
for (v in names(VARS)) painel(v, VARS[[v]])
par(mar = c(0, 1, 0, 1)); plot.new()
legend("left", bty = "n", cex = 0.85, text.col = INK, y.intersp = 1.15,
  legend = c(paste("publicado,", NOME_USO), "preenchido = confirmado por enriquecimento", "modelo global: limiar robusto (faixa cinza = intervalo)", "traço = limiar do modelo local da bacia", "LP000 = Id na tabela seguinte"),
  pch = c(21, 21, 21, 21, 21, NA, NA), pt.bg = c("white", "white", "white", MUTED, INK, NA, NA), col = c(USO, MUTED, INK, MUTED, NA), lty = c(NA, NA, NA, NA, NA, 1, NA), pt.cex = 1.1)
dev.off()
write.csv(bind_rows(tabela) |> arrange(ordem_painel, linha), file.path(raiz, "saida/figura_limiares_tabela.csv"), row.names = FALSE, fileEncoding = "UTF-8"); cat("figura gravada\n")
