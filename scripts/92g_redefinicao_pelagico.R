# 92g_redefinicao_pelagico.R  (GUIA_LIMIARES_PUBLICADOS.md, step L3 for the pelagic; Gate 10 item 5: after option B; PI decision of 06/11/2026)
# 1 L2b again with the fauna that did not exist in 92b: enrichment of each published pelagic threshold on the Gradient Forest curves of the
#   option-B groups and of the level-2 groups (rule of 92b: gain inside the published interval against 999 random windows of the same width;
#   single values get a window of WMIN x sampled range). Thresholds on z, sst and t_z with a value.
# 2 Candidates, rule of 92c: use validacao or apoio, variable z (the only pelagic variable with such thresholds), not detected by the model and
#   ratio >= 1 in at least one group, on the original or on the new curves. All published z values are tested anyway, as a report: those
#   at a zone limit of the map (200 and 1,000 m in 91b) need no partition, and their fauna test says whether that limit has support.
# 3 Test per value and context, inside the study box (global study = whole ocean): context = HP (93b), plus "todos" as a report (no HP
#   control, so z and region are confounded). Fauna units above vs. below the value, pares_classes() per group, one row per ecological
#   group (level 2 and option B share records: option B stays). Consensus and decision as in 92c, with equivalence also needing median
#   omega2 < 0.01 (Gate 11 item 6). A partition needs >= 10 cells on each side of the context in some month x depth layer.
# Outputs: saida/limiares_pelagico_enriquecimento.csv, saida/limiares_pelagico_testes.csv, saida/limiares_pelagico_redefinicao.csv
suppressPackageStartupMessages({library(dplyr); library(terra); library(arrow); library(parallel)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); source(file.path(raiz, "scripts/funcoes_grade.R"))
WMIN <- 0.05; NNULL <- 999; AREA_MIN <- 10; ZONAS <- c(200, 1000); SEM_SINAL <- c("Euphausiacea", "Hyperiidea"); ZP <- c(0, 50, 100, 200, 500, 1000, 2000)
g25 <- ease2_grade("25km")
en <- read.csv(file.path(raiz, "saida/limiares_enriquecimento.csv"), encoding = "UTF-8") |> filter(compartimento == "pelagico")
lp <- read.csv(file.path(raiz, "saida/limiares_publicados_padronizados.csv"), encoding = "UTF-8") |> filter(compartimento == "pelagico") |> select(id, referencia, regiao, caixa_wkt)

# ---------------- 1 enrichment on the new curves ----------------
cB <- read.csv(file.path(raiz, "saida/etapa3_n1_gf_curvas_b.csv"), encoding = "UTF-8") |> filter(comp == "pelagico") |> mutate(fonte = "B") |> select(-comp)
cN <- read.csv(file.path(raiz, "saida/pelagico_n2_gf_curvas.csv"), encoding = "UTF-8") |> mutate(fonte = "n2")
cv <- bind_rows(cB, cN)
enriq <- function(x, y, a, b) { o <- order(x); x <- x[o]; y <- y[o]; f <- function(v) approx(x, y, v, rule = 2, ties = max)$y   # as in 92b
  lo <- min(x); hi <- max(x); w <- b - a; if (w <= 0) { w <- WMIN * (hi - lo); a <- a - w / 2; b <- b + w / 2 }
  if (b < lo || a > hi || w >= hi - lo) return(c(razao = NA, p = NA))
  g <- f(b) - f(a); s <- runif(NNULL, lo, hi - w); nul <- f(s + w) - f(s); c(razao = g / mean(nul), p = (1 + sum(nul >= g)) / (NNULL + 1)) }
set.seed(2027)
lim <- en |> distinct(id, variavel_modelo, limiar_min, limiar_max, uso, situacao) |> group_by(id) |> slice(1) |> ungroup()
ee <- bind_rows(lapply(seq_len(nrow(lim)), function(i) { l <- lim[i, ]; d <- cv |> filter(preditor == l$variavel_modelo)
  bind_rows(lapply(split(d, paste(d$fonte, d$grupo)), function(k) { if (nrow(k) < 5) return(NULL); e <- enriq(k$x, k$y, l$limiar_min, l$limiar_max)
    data.frame(id = l$id, variavel = l$variavel_modelo, limiar_min = l$limiar_min, limiar_max = l$limiar_max, fonte = k$fonte[1], grupo = k$grupo[1], razao = e[["razao"]], p = e[["p"]]) })) }))
write.csv(ee, file.path(raiz, "saida/limiares_pelagico_enriquecimento.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ---------------- 2 candidates ----------------
tab <- en |> filter(variavel_modelo == "z", uso %in% c("validacao", "apoio")) |> group_by(id) |>
  summarise(limiar_min = first(limiar_min), limiar_max = first(limiar_max), uso = first(uso), situacao = first(situacao), razao_orig_max = suppressWarnings(max(razao, na.rm = TRUE)), .groups = "drop") |>
  left_join(ee |> group_by(id) |> summarise(razao_novos_max = suppressWarnings(max(razao, na.rm = TRUE))), by = "id") |> left_join(lp, by = "id") |>
  mutate(valor = (limiar_min + limiar_max) / 2, limite_de_zona = valor %in% ZONAS, candidato = situacao == "nao_detectado" & pmax(razao_orig_max, razao_novos_max, na.rm = TRUE) >= 1)
print(as.data.frame(tab |> select(id, referencia, valor, uso, situacao, razao_orig_max, razao_novos_max, limite_de_zona, candidato)))

# ---------------- 3 fauna tests ----------------
curvas_p <- open_dataset(file.path(raiz, "saida/etapa3_n1_gf_curvas.csv"), format = "csv") |> filter(comp == "pelagico") |> collect() |> as.data.frame()
transf <- function(E) { vs <- intersect(unique(curvas_p$preditor), names(E))
  sapply(vs, function(v) { s <- 0; for (g in unique(curvas_p$grupo[curvas_p$preditor == v])) { c1 <- curvas_p[curvas_p$grupo == g & curvas_p$preditor == v, ]
    s <- s + approx(c1$x, c1$y, E[[v]], rule = 2, ties = mean)$y }; s }) }
mais_proximo <- function(T, C) { d <- sapply(seq_len(nrow(C)), function(k) colSums((t(T) - C[k, ])^2)); max.col(-d, ties.method = "first") }
mod <- readRDS(file.path(raiz, "intermediario/etapa3_n1/habitats_novos_modelo.rds")); CP <- mod$kmP$centers; rm(mod); gc(verbose = FALSE)
ordena <- function(h) h[order(as.integer(sub("^HP", "", h)))]
fus <- read.csv(file.path(raiz, "saida/pelagico_b_fusoes.csv")); nf <- tapply(fus$hp_n1, fus$final, function(x) ordena(unique(x))[1]); fus$hp <- as.character(nf[as.character(fus$final)])
cel_de <- function(id) { p <- do.call(rbind, strsplit(id, "_")); cellFromRowCol(g25, as.integer(p[, 2]), as.integer(p[, 3])) }
arqs <- c(list.files(file.path(raiz, "intermediario/etapa3_n1"), pattern = "^pelagico_.*[.]rds$", full.names = TRUE), list.files(file.path(raiz, "intermediario/etapa3_n1b"), pattern = "^pelagico_.*[.]rds$", full.names = TRUE),
          list.files(file.path(raiz, "intermediario/pelagico_n2"), pattern = "[.]rds$", full.names = TRUE))
arqs <- arqs[!grepl(paste(SEM_SINAL, collapse = "|"), arqs)]
sem_ac <- function(x) chartr("áéíóúâêôãõç", "aeiouaeoaoc", x)
U <- lapply(arqs, function(a) { r <- readRDS(a); ll <- project(xyFromCell(g25, cel_de(r$meta$cell_id)), from = crs(g25), to = "EPSG:4326")
  list(arq = a, grupo = r$grupo, fam = sem_ac(r$grupo), B = grepl("etapa3_n1b", a), z = r$E$z, lon = ll[, 1], lat = ll[, 2], hp = fus$hp[mais_proximo(transf(r$E), CP)]) })
na_caixa <- function(lon, lat, wkt) { if (is.na(wkt) || wkt == "") return(rep(TRUE, length(lon))); e <- ext(vect(wkt, crs = "EPSG:4326")); lon >= e$xmin & lon <= e$xmax & lat >= e$ymin & lat <= e$ymax }
jobs <- bind_rows(lapply(seq_len(nrow(tab)), function(k) { x <- tab[k, ]; cx <- unique(unlist(lapply(U, function(u) { s <- na_caixa(u$lon, u$lat, x$caixa_wkt); u$hp[s] })))
  data.frame(k = k, contexto = c("todos", ordena(cx))) }))
testa <- function(j, U, tab) { suppressPackageStartupMessages({library(vegan); library(dplyr); library(terra)}); source(file.path(raiz, "scripts/funcoes_permanova.R")); x <- tab[j$k, ]
  bind_rows(lapply(U, function(u) { ok <- na_caixa(u$lon, u$lat, x$caixa_wkt); if (j$contexto != "todos") ok <- ok & u$hp == j$contexto
    lado <- ifelse(u$z >= x$valor, "fundo", "raso")[ok]; if (sum(ok) < 10 || min(table(factor(lado, c("raso", "fundo")))) < 5) return(NULL)
    r <- readRDS(u$arq); Dm <- as.matrix(sqrt(betadiver(r$X[ok, , drop = FALSE], "sim"))); set.seed(7500 + j$k)
    p <- pares_classes(Dm, lado, r$meta$bloco[ok], un_min = 5, n_classe = 10, B = 100)
    if (!is.null(p)) cbind(id = x$id, valor = x$valor, contexto = j$contexto, grupo = u$grupo, fam = u$fam, opcao_B = u$B, unidades_raso = sum(lado == "raso"), unidades_fundo = sum(lado == "fundo"), p) })) }
cl <- makeCluster(min(14, nrow(jobs))); clusterExport(cl, c("raiz", "na_caixa")); tes <- bind_rows(parLapplyLB(cl, split(jobs, seq_len(nrow(jobs))), testa, U = U, tab = tab)); stopCluster(cl)
cat(format(Sys.time()), "linhas de teste:", nrow(tes), "\n")
if (nrow(tes)) { tes <- tes |> arrange(desc(opcao_B)) |> distinct(id, contexto, fam, .keep_all = TRUE); write.csv(tes, file.path(raiz, "saida/limiares_pelagico_testes.csv"), row.names = FALSE, fileEncoding = "UTF-8") }
ct <- if (nrow(tes)) tes |> filter(!is.na(separa)) |> group_by(id, contexto) |> summarise(grupos = n(), distintos = sum(separa), grupos_nomes = paste(grupo, collapse = "; "), omega2_med = median(omega2),
    unidades_raso = sum(unidades_raso), unidades_fundo = sum(unidades_fundo), .groups = "drop") else data.frame(id = character(), contexto = character())

# ---------------- minimum area and decision ----------------
hpt <- rast(file.path(raiz, "saida/habitats_novos_pelagico_mes_z.tif")); hs0 <- ordena(unique(fus$hp_n1)); hp_val <- fus$hp[match(hs0, fus$hp_n1)]
xy <- project(xyFromCell(g25, seq_len(ncell(g25))), from = crs(g25), to = "EPSG:4326")
area <- function(v, wkt, h) { s <- na_caixa(xy[, 1], xy[, 2], wkt); m <- function(zs) max(sapply(1:12, function(mm) max(c(0, sapply(zs, function(z) sum(hp_val[values(hpt[[sprintf("m%02d_z%d", mm, z)]])[s, 1]] == h, na.rm = TRUE))))))
  c(m(ZP[ZP < v]), m(ZP[ZP >= v])) }
red <- jobs |> mutate(id = tab$id[k]) |> left_join(ct, by = c("id", "contexto")) |> left_join(tab |> select(id, referencia, regiao, valor, uso, situacao, limite_de_zona, candidato, caixa_wkt), by = "id") |>
  mutate(grupos = coalesce(grupos, 0L), distintos = coalesce(distintos, 0L), frac = ifelse(grupos > 0, distintos / grupos, NA),
    consenso = case_when(grupos == 0 ~ "sem teste", grupos >= 2 & frac <= 1/3 & omega2_med < 0.01 ~ "equivalente", frac >= 2/3 ~ "distinto", grupos == 1 & distintos == 0 ~ "fraco", TRUE ~ "depende do grupo"))
red$cel_raso <- red$cel_fundo <- NA_real_
for (i in which(red$candidato & !red$limite_de_zona & red$contexto != "todos")) { a <- area(red$valor[i], red$caixa_wkt[i], red$contexto[i]); red$cel_raso[i] <- a[1]; red$cel_fundo[i] <- a[2] }
red <- red |> mutate(decisao = case_when(contexto == "todos" ~ "relato", limite_de_zona ~ "já é limite de zona", !candidato ~ "não é candidato", pmin(cel_raso, cel_fundo) < AREA_MIN ~ "manter (área mínima)",
    consenso == "equivalente" ~ "manter (conflito com a literatura)", TRUE ~ "partir")) |> select(-caixa_wkt, -k)
write.csv(red, file.path(raiz, "saida/limiares_pelagico_redefinicao.csv"), row.names = FALSE, fileEncoding = "UTF-8")
options(width = 220); print(ee |> filter(!is.na(razao), razao >= 2, p <= 0.05) |> mutate(across(where(is.double), ~ round(.x, 3))) |> as.data.frame())
print(red |> filter(grupos > 0) |> mutate(across(where(is.double), ~ round(.x, 3))) |> select(id, valor, contexto, grupos, distintos, omega2_med, unidades_raso, unidades_fundo, consenso, decisao) |> as.data.frame())
print(table(red$decisao, red$consenso)); cat(format(Sys.time()), "FIM\n")
