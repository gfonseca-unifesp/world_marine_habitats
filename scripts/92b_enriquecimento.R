# 92b_enriquecimento.R  (GUIA_LIMIARES_PUBLICADOS.md, step L2b: tolerance-free test of published thresholds)
# Does the published interval concentrate more faunal turnover on our Gradient Forest curve than random intervals of the same
# width? For each published threshold x project group (saida/limiares_comparacao.csv) with a value and a model predictor:
#   gain = rise of the cumulative-importance curve of the group inside [limiar_min, limiar_max]; single values get a window of
#   width WMIN x sampled range centred on the value (PROVISIONAL until Gate 10);
#   null = 999 windows of the same width placed at random inside the range sampled by our units (the curve's range);
#   enrichment ratio = gain / mean(null gains); empirical p = (1 + #null >= gain) / 1000. Effect size (ratio) is read first.
# Repeated with the local GF curve of the study basin (88c, saida/etapa3_n1_transf_curvas_locais.csv), when it exists.
# Output: saida/limiares_enriquecimento.csv
suppressPackageStartupMessages({library(dplyr); library(arrow)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); WMIN <- 0.05; NNULL <- 999   # WMIN PROVISIONAL (Gate 10)
cmp <- read.csv(file.path(raiz, "saida/limiares_comparacao.csv"), encoding = "UTF-8")
cmp <- cmp |> filter(!situacao %in% c("variavel_ausente", "sem_valor", "sem_quebra"))
bac_nome <- read.csv(file.path(raiz, "camadas/bacias_legenda.csv"))
gc_ds <- open_dataset(file.path(raiz, "saida/etapa3_n1_gf_curvas.csv"), format = "csv")             # 855 MB: filter before collect
fl <- file.path(raiz, "saida/etapa3_n1_transf_curvas_locais.csv"); loc <- if (file.exists(fl)) read.csv(fl) |> mutate(bacia_nome = bac_nome$bacia[match(bacia, bac_nome$codigo)]) else NULL
enriq <- function(x, y, a, b) {
  o <- order(x); x <- x[o]; y <- y[o]; f <- function(v) approx(x, y, v, rule = 2, ties = max)$y
  lo <- min(x); hi <- max(x); w <- b - a; if (w <= 0) { w <- WMIN * (hi - lo); a <- a - w / 2; b <- b + w / 2 }
  if (b < lo || a > hi || w >= hi - lo) return(c(razao = NA, p = NA, ganho = NA))
  g <- f(b) - f(a); s <- runif(NNULL, lo, hi - w); nul <- f(s + w) - f(s)
  c(razao = g / mean(nul), p = (1 + sum(nul >= g)) / (NNULL + 1), ganho = g / (max(y) - min(y)))
}
set.seed(2026)
cur_cache <- list()
curva <- function(cp, g, v) { k <- paste(cp, g, v); if (is.null(cur_cache[[k]])) cur_cache[[k]] <<- gc_ds |> filter(comp == cp, grupo == g, preditor == v) |> select(x, y) |> collect(); cur_cache[[k]] }
res <- bind_rows(lapply(seq_len(nrow(cmp)), function(i) { r <- cmp[i, ]
  cv <- curva(r$compartimento, r$grupo_projeto, r$variavel_modelo)
  e <- if (nrow(cv) >= 5) enriq(cv$x, cv$y, r$limiar_min, r$limiar_max) else c(razao = NA, p = NA, ganho = NA)
  eb <- c(razao = NA, p = NA, ganho = NA)
  if (!is.null(loc) && !r$bacia %in% c("global", "")) { cl <- loc |> filter(comp == r$compartimento, grupo == r$grupo_projeto, preditor == r$variavel_modelo, bacia_nome == r$bacia)
    if (nrow(cl) >= 5) eb <- enriq(cl$x, cl$y, r$limiar_min, r$limiar_max) }
  data.frame(id = r$id, compartimento = r$compartimento, grupo_projeto = r$grupo_projeto, comunidade = r$comunidade, variavel_modelo = r$variavel_modelo,
    limiar_min = r$limiar_min, limiar_max = r$limiar_max, situacao = r$situacao, bacia = r$bacia,
    razao = e[["razao"]], p = e[["p"]], ganho = e[["ganho"]], razao_bacia = eb[["razao"]], p_bacia = eb[["p"]]) }))
write.csv(res, file.path(raiz, "saida/limiares_enriquecimento.csv"), row.names = FALSE, fileEncoding = "UTF-8")
options(width = 200)
print(res |> group_by(id, variavel_modelo) |> summarise(grupos = n(), razao_med = median(razao, na.rm = TRUE), enriquecidos = sum(razao >= 1 & p <= 0.05, na.rm = TRUE),
  razao_bacia_med = median(razao_bacia, na.rm = TRUE), .groups = "drop") |> mutate(across(where(is.double), ~ round(.x, 2))) |> as.data.frame())
cat(format(Sys.time()), "FIM\n")
