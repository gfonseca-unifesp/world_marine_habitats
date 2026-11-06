# 92c_redefinicao.R  (GUIA_LIMIARES_PUBLICADOS.md, step L3; Gate 10 item 3: depth and bottom temperature in the benthos only)
# Published thresholds that our model does not detect are tested as new class boundaries with our own fauna.
#  candidates: benthic rows of saida/limiares_enriquecimento.csv with variable prof or tfundo_med, use validacao or apoio (relato
#    types such as richness do not split habitats), situacao nao_detectado and enrichment ratio >= 1 in at least one group;
#    split value = the published value (midpoint of the published interval);
#  per candidate and per level-III context inside the study box (depth band x HB, or depth band x protected feature, as in
#    91_mapa_unido_bentos.R): fauna units shallower (colder) vs. deeper (warmer) than the value, pares_classes() per taxonomic
#    group (omega2 >= 0.01 and power >= 0.8), consensus as in 91;
#  decision: distinct -> split, evidence "distinta pela fauna"; equivalent (not distinct in >= 2/3 of >= 2 groups) -> no split,
#    conflict with the literature recorded; untestable, single group without separating, or "depende do grupo" -> split, evidence
#    "separada pela literatura" (Gate 9 principle: separate when any evidence separates); either side of the box with < 10 cells
#    of that context -> no split (minimum area).
# 91_mapa_unido_bentos.R applies the splits (criterion 9) to its base classes.
# Output: saida/limiares_candidatos_redefinicao.csv
suppressPackageStartupMessages({library(dplyr); library(terra); library(arrow); library(parallel)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); source(file.path(raiz, "scripts/funcoes_grade.R")); source(file.path(raiz, "scripts/funcoes_permanova.R"))
dir_t <- file.path(raiz, "intermediario/etapa3_n1"); PROT <- c("C", "D"); AREA_MIN <- 10; VARS_L3 <- c("prof", "tfundo_med")
faixa <- function(prof) findInterval(prof, c(200, 2000), left.open = TRUE) + 1
g25 <- ease2_grade("25km"); M <- lapply(1:12, function(m) rast(file.path(raiz, sprintf("camadas/glorys/mes/glorys_%02d_25km.tif", m))))

# ---------------- candidates ----------------
en <- read.csv(file.path(raiz, "saida/limiares_enriquecimento.csv"), encoding = "UTF-8")
lp <- read.csv(file.path(raiz, "saida/limiares_publicados_padronizados.csv"), encoding = "UTF-8") |> select(id, referencia, regiao, caixa_wkt)
cands <- en |> filter(compartimento == "bentos", variavel_modelo %in% VARS_L3, uso %in% c("validacao", "apoio")) |> group_by(id) |>
  summarise(variavel = first(variavel_modelo), limiar_min = first(limiar_min), limiar_max = first(limiar_max), bacia = first(bacia), uso = first(uso),
    grupos_nao_detectado = sum(situacao == "nao_detectado"), razao_max_nao_detectado = suppressWarnings(max(razao[situacao == "nao_detectado"], na.rm = TRUE)), .groups = "drop") |>
  filter(grupos_nao_detectado > 0, is.finite(razao_max_nao_detectado), razao_max_nao_detectado >= 1) |>
  left_join(lp, by = "id") |> filter(!is.na(caixa_wkt), caixa_wkt != "") |> mutate(valor = (limiar_min + limiar_max) / 2)
cat("candidatos L3:", nrow(cands), "\n"); print(cands |> select(id, referencia, variavel, valor, uso, razao_max_nao_detectado) |> as.data.frame())

# ---------------- cells: level-III context, value, box ----------------
e25 <- read_parquet(file.path(raiz, "camadas/ebu/ebu_celula_25km.parquet"))
p <- do.call(rbind, strsplit(e25$cell_id, "_")); cel_b <- cellFromRowCol(g25, as.integer(p[, 2]), as.integer(p[, 3]))
lh <- read.csv(file.path(raiz, "saida/habitats_novos_legenda.csv")) |> filter(arquivo == "bentos_25km")
tf <- sapply(1:12, function(k) values(M[[k]][["tfundo"]])[cel_b, 1])
C <- data.frame(cel = cel_b, prof = values(M[[1]][["prof_glorys"]])[cel_b, 1], tfundo_med = rowMeans(tf),
  hb = lh$habitat[values(rast(file.path(raiz, "saida/habitats_novos_bentos_25km.tif")))[cel_b, 1]], grupo_ebu = e25$ebu_grupo); rm(tf)
C <- C[complete.cases(C), ]; C$faixa <- faixa(C$prof)
C$nivel3 <- ifelse(C$grupo_ebu %in% PROT, paste(C$faixa, C$grupo_ebu, sep = "|"), paste(C$faixa, C$hb, sep = "|"))
ll <- project(xyFromCell(g25, C$cel), from = crs(g25), to = "EPSG:4326"); C$lon <- ll[, 1]; C$lat <- ll[, 2]
na_caixa <- function(lon, lat, wkt) { e <- ext(vect(wkt, crs = "EPSG:4326")); lon >= e$xmin & lon <= e$xmax & lat >= e$ymin & lat <= e$ymax }

# ---------------- fauna units ----------------
arqs <- list.files(dir_t, pattern = "^bentos_.*\\.rds$", full.names = TRUE)
U <- lapply(arqs, function(a) { r <- readRDS(a); p <- do.call(rbind, strsplit(rownames(r$X), "_"))
  i <- match(cellFromRowCol(g25, as.integer(p[, 2]), as.integer(p[, 3])), C$cel); list(arq = a, grupo = r$grupo, i = i) })

# ---------------- tests ----------------
jobs <- bind_rows(lapply(seq_len(nrow(cands)), function(k) { cd <- cands[k, ]; dentro <- na_caixa(C$lon, C$lat, cd$caixa_wkt)
  lado <- C[[cd$variavel]] >= cd$valor
  C[dentro, ] |> mutate(lado = lado[dentro]) |> group_by(nivel3) |> summarise(celulas_raso = sum(!lado), celulas_fundo = sum(lado), .groups = "drop") |>
    filter(celulas_raso > 0, celulas_fundo > 0) |> mutate(k = k) }))
cat("contextos que atravessam o valor:", nrow(jobs), "\n")
testa <- function(j, U, C, cands) {
  suppressPackageStartupMessages({library(vegan); library(dplyr); library(terra)}); source(file.path(raiz, "scripts/funcoes_permanova.R"))
  cd <- cands[j$k, ]; e <- ext(vect(cd$caixa_wkt, crs = "EPSG:4326"))
  bind_rows(lapply(U, function(u) { ok <- !is.na(u$i); ii <- u$i; ok[ok] <- C$nivel3[ii[ok]] == j$nivel3 & C$lon[ii[ok]] >= e$xmin & C$lon[ii[ok]] <= e$xmax &
      C$lat[ii[ok]] >= e$ymin & C$lat[ii[ok]] <= e$ymax
    if (sum(ok) < 10) return(NULL)
    lado <- ifelse(C[[cd$variavel]][ii[ok]] >= cd$valor, "fundo", "raso"); if (min(table(factor(lado, c("raso", "fundo")))) < 5) return(NULL)
    r <- readRDS(u$arq); Dm <- as.matrix(sqrt(betadiver(r$X[ok, , drop = FALSE], "sim"))); set.seed(7300 + j$k)
    x <- pares_classes(Dm, lado, r$meta$bloco[ok], un_min = 5, n_classe = 10, B = 100)
    if (!is.null(x)) cbind(id = cd$id, nivel3 = j$nivel3, grupo = u$grupo, unidades_raso = sum(lado == "raso"), unidades_fundo = sum(lado == "fundo"), x) }))
}
cl <- makeCluster(min(14, max(1, nrow(jobs)))); clusterExport(cl, "raiz")
tes <- bind_rows(parLapplyLB(cl, split(jobs, seq_len(nrow(jobs))), testa, U = U, C = C, cands = cands)); stopCluster(cl)
if (nrow(tes)) tes <- tes |> filter(!is.na(separa))
cat(format(Sys.time()), "testes:", nrow(tes), "\n")
cons <- if (nrow(tes)) tes |> group_by(id, nivel3) |> summarise(grupos = n(), distintos = sum(separa), grupos_nomes = paste(grupo, collapse = "; "),
  omega2_med = median(omega2), .groups = "drop") else data.frame(id = character(), nivel3 = character(), grupos = integer(), distintos = integer(), grupos_nomes = character(), omega2_med = numeric())
red <- jobs |> mutate(id = cands$id[k]) |> left_join(cons, by = c("id", "nivel3")) |> mutate(grupos = coalesce(grupos, 0L), distintos = coalesce(distintos, 0L),
    frac = ifelse(grupos > 0, distintos / grupos, NA),
    consenso = case_when(grupos == 0 ~ "sem teste", grupos >= 2 & frac <= 1/3 ~ "equivalente", frac >= 2/3 & grupos >= 1 & distintos >= 1 ~ "distinto",
      grupos == 1 & distintos == 0 ~ "fraco", TRUE ~ "depende do grupo"),
    decisao = case_when(pmin(celulas_raso, celulas_fundo) < AREA_MIN ~ "manter (área mínima)", consenso == "equivalente" ~ "manter (conflito com a literatura)", TRUE ~ "partir"),
    evidencia = case_when(decisao != "partir" ~ NA_character_, consenso == "distinto" ~ "distinta pela fauna", TRUE ~ "separada pela literatura")) |>
  left_join(cands |> select(id, referencia, regiao, variavel, valor, limiar_min, limiar_max, caixa_wkt, uso, razao_max_nao_detectado), by = "id") |>
  select(id, referencia, regiao, variavel, valor, limiar_min, limiar_max, uso, razao_max_nao_detectado, nivel3, celulas_raso, celulas_fundo,
    grupos, distintos, grupos_nomes, omega2_med, consenso, decisao, evidencia, caixa_wkt)
write.csv(red, file.path(raiz, "saida/limiares_candidatos_redefinicao.csv"), row.names = FALSE, fileEncoding = "UTF-8")
if (nrow(tes)) write.csv(tes, file.path(raiz, "saida/limiares_redefinicao_testes.csv"), row.names = FALSE, fileEncoding = "UTF-8")
options(width = 220); print(red |> select(id, valor, nivel3, celulas_raso, celulas_fundo, grupos, distintos, consenso, decisao) |> as.data.frame())
print(table(red$decisao, red$evidencia, useNA = "ifany")); cat(format(Sys.time()), "FIM\n")
