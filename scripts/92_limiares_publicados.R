# 92_limiares_publicados.R  (published thresholds x global-model thresholds; GUIA_LIMIARES_PUBLICADOS.md, steps L1 and L2)
# Input: limiares_publicados_v0.csv (33 thresholds from 23 studies, exploratory search of 30/10/2026; values from abstracts).
# L1 standardise: project group (grupo_projeto; "comunidade" when the study mixes groups), predictor name in the model
#    (variavel_modelo; NA when the variable is not in the model), approximate lat/lon box per study region (caixa_wkt,
#    caixa_aproximada = sim), basin of the box centre from 86b_bacias.R (bacia_conferida), independence from our data
#    (dados_no_obis: parcial = an OBIS dataset whose title matches the study region has records of the compartment in the box;
#    nao = no records of the compartment in the box; indeterminado = otherwise; always to be confirmed by hand), and for the
#    pelagic whether the published depth coincides with a usual sampling-stratum limit (limite_igual_estrato).
# L2 compare each published threshold x project group with the robust thresholds of the model (89_ligacao_habitats.R) on the
#    same predictor (any group of the compartment when grupo_projeto = comunidade). Tolerance TOL (relative, PROVISIONAL until
#    Gate 10): confirmado = published interval widened by TOL overlaps the robust interval [de, ate]; proximo = overlaps with
#    2 x TOL; fora_do_alcance = published value outside the range sampled by our units; nao_detectado = otherwise;
#    variavel_ausente; tipo_incompativel (riqueza, nicho_termico, habitat; reported only, PROVISIONAL); sem_valor;
#    sem_quebra (published negative result). Repeated with the local thresholds of the study basin (88c, situacao_bacia).
# Outputs: saida/limiares_publicados_padronizados.csv, saida/limiares_comparacao.csv
suppressPackageStartupMessages({library(dplyr); library(terra); library(arrow); library(tidyr)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); source(file.path(raiz, "scripts/funcoes_grade.R"))
TOL <- 0.10; ESTRATOS <- c(50, 100, 200, 300, 500, 1000)                       # PROVISIONAL (Gate 10)
TIPO_INCOMP <- c("riqueza", "nicho_termico", "habitat"); SEM_SINAL <- c("Euphausiacea", "Hyperiidea")
lp <- read.csv(file.path(raiz, "limiares_publicados_v0.csv"), encoding = "UTF-8", stringsAsFactors = FALSE)

# ---------------- L1 standardise ----------------
grupo_de <- function(cp, g) {
  if (cp == "bentos") return(if (grepl("^Bivalvia$", g)) "Bivalvia" else if (grepl("Asteroidea|Ophiuroidea|Echinoidea|Holothuroidea", g)) "Echinodermata"
    else if (grepl("^Polychaeta$|^Annelida$", g)) "Polychaeta" else if (grepl("^Porifera$", g)) "Porifera" else if (grepl("^Anthozoa$", g)) "Anthozoa" else "comunidade")
  if (grepl("Siphonophorae|salp|Thaliacea", g, ignore.case = TRUE)) "Siphonophorae e Thaliacea" else if (grepl("Euphausi", g)) "Euphausiacea"
  else if (grepl("Chaetognatha", g)) "Chaetognatha" else if (grepl("anfípodes|Hyperiid", g, ignore.case = TRUE)) "Hyperiidea"
  else if (grepl("^copépodes|Calanoida", g, ignore.case = TRUE)) "Calanoida" else "comunidade"
}
var_de <- function(cp, v) switch(v, depth = if (cp == "bentos") "prof" else "prof", z = "z", bottom_temperature = "tfundo_med",
  sst = if (cp == "bentos") "sst_med" else "sst", temperature = "t_z", water_mass_temperature = "t_z", NA_character_)
# approximate boxes (lon_min, lon_max, lat_min, lat_max; lon_min > lon_max crosses the date line); NA = global or basin-wide
CAIXAS <- list(
  "Porcupine" = c(-16, -10, 48, 52), "Ártico Central" = c(-20, 180, 78, 90), "Kara, Laptev" = c(55, 180, 70, 82), "Havaí" = c(-161, -154, 18, 23),
  "Kyushu-Palau" = c(133, 137, 15, 30), "Clarion-Clipperton" = c(-160, -115, 5, 20), "Fossas do Pacífico NW" = c(125, 150, 10, 40),
  "Mar do Leste" = c(128, 138, 35, 43), "Grande Barreira" = c(142, 154, -25, -10), "Pacífico tropical oeste" = c(125, 160, 0, 20),
  "Prydz" = c(70, 80, -70, -66), "Mar do Sul da China" = c(105, 121, 5, 23), "Setores Pacífico Leste e Índico do Oceano Austral" = c(20, -70, -70, -50),
  "Sul do Golfo do México" = c(-97, -86, 18, 24), "Zona de Transição do Pacífico Norte" = c(150, -130, 28, 45), "subártico" = c(140, -130, 40, 60),
  "Mar de Ross" = c(160, -170, -78, -70), "Cabo Verde" = c(-26, -16, 12, 20))
CHAVES <- c("Porcupine" = "Porcupine", "Ártico Central" = "Fram|Beaufort|Central Arctic|Arctic Ocean", "Kara, Laptev" = "Kara|Laptev|Siberian",
  "Havaí" = "Hawai", "Kyushu-Palau" = "Kyushu|Palau", "Clarion-Clipperton" = "Clarion|Clipperton|CCZ|ABYSSLINE", "Fossas do Pacífico NW" = "Mariana|Izu|Ogasawara|Ryukyu|Japan Trench|hadal",
  "Mar do Leste" = "East Sea|Korea", "Grande Barreira" = "Great Barrier|GBR", "Pacífico tropical oeste" = "Western Pacific|West Pacific",
  "Prydz" = "Prydz", "Mar do Sul da China" = "South China Sea", "Setores Pacífico Leste e Índico do Oceano Austral" = "Southern Ocean|Antarctic",
  "Sul do Golfo do México" = "Gulf of Mexico|Campeche", "Zona de Transição do Pacífico Norte" = "Transition Zone|North Pacific", "subártico" = "subarctic|Subarctic|North Pacific",
  "Mar de Ross" = "Ross Sea", "Cabo Verde" = "Cape Verde|Cabo Verde")
chave_de <- function(reg) { k <- names(CAIXAS)[vapply(names(CAIXAS), function(x) grepl(x, reg, fixed = TRUE), TRUE)]; if (length(k)) k[1] else NA_character_ }
wkt <- function(b) { poli <- function(x1, x2) sprintf("POLYGON((%s %s, %s %s, %s %s, %s %s, %s %s))", x1, b[3], x2, b[3], x2, b[4], x1, b[4], x1, b[3])
  if (b[1] <= b[2]) poli(b[1], b[2]) else sprintf("MULTIPOLYGON((%s),(%s))", sub("^POLYGON", "", poli(b[1], 180)), sub("^POLYGON", "", poli(-180, b[2]))) }
na_caixa <- function(lon, lat, b) lat >= b[3] & lat <= b[4] & (if (b[1] <= b[2]) lon >= b[1] & lon <= b[2] else lon >= b[1] | lon <= b[2])
bac <- rast(file.path(raiz, "camadas/bacias_25km.tif")); bac_nome <- read.csv(file.path(raiz, "camadas/bacias_legenda.csv"))
bacia_centro <- function(b) { lon <- if (b[1] <= b[2]) mean(b[1:2]) else { m <- (b[1] + b[2] + 360) / 2; if (m > 180) m - 360 else m }
  xy <- project(cbind(lon, mean(b[3:4])), from = "EPSG:4326", to = "EPSG:6933"); v <- terra::extract(bac, xy)[, 1]
  if (is.na(v)) NA_character_ else bac_nome$bacia[match(v, bac_nome$codigo)] }

lp <- lp |> mutate(grupo_projeto = mapply(grupo_de, compartimento, grupo_publicado), variavel_modelo = mapply(var_de, compartimento, variavel_padrao),
  variavel_no_modelo = ifelse(is.na(variavel_modelo), "nao", "sim"), chave = vapply(regiao, chave_de, ""),
  limiar_min = ifelse(is.na(limiar_min) & !is.na(limiar_max), limiar_max, limiar_min),          # "up to X" read as X
  limite_igual_estrato = ifelse(compartimento == "pelagico" & variavel_padrao == "z" & !is.na(limiar_min) & limiar_min == limiar_max & limiar_min %in% ESTRATOS, "sim", "nao"))
lp$caixa_wkt <- vapply(lp$chave, function(k) if (is.na(k)) NA_character_ else wkt(CAIXAS[[k]]), "")
lp$caixa_aproximada <- ifelse(is.na(lp$chave), NA, "sim")
lp$bacia_conferida <- vapply(lp$chave, function(k) if (is.na(k)) NA_character_ else bacia_centro(CAIXAS[[k]]), "")
# independence: records of the compartment in the box (curated level-1 tables) and OBIS datasets whose title matches the region
cj <- read.csv(file.path(raiz, "bruto/obis_conjuntos.csv"), encoding = "UTF-8") |> select(dataset_id, titulo)
indep <- bind_rows(lapply(unique(na.omit(lp$chave)), function(k) { b <- CAIXAS[[k]]
  bind_rows(lapply(c("bentos", "pelagico"), function(cp) {
    d <- open_dataset(file.path(raiz, sprintf("saida/tabela_certa_%s_n1.parquet", cp))) |> filter(lat >= b[3], lat <= b[4]) |> select(dataset, lon, lat) |> collect()
    d <- d[na_caixa(d$lon, d$lat, b), ]; tab <- d |> count(dataset, name = "registros") |> left_join(cj, by = c("dataset" = "dataset_id"))
    m <- tab |> filter(grepl(CHAVES[[k]], titulo, ignore.case = TRUE))
    data.frame(chave = k, compartimento = cp, registros_caixa = nrow(d), conjuntos_caixa = nrow(tab),
      conjuntos_compativeis = if (nrow(m)) paste(sprintf("%s (%d)", substr(m$titulo, 1, 60), m$registros), collapse = "; ") else NA_character_) })) }))
lp <- lp |> left_join(indep, by = c("chave", "compartimento")) |>
  mutate(dados_no_obis = case_when(is.na(chave) ~ "indeterminado", registros_caixa == 0 ~ "nao", !is.na(conjuntos_compativeis) ~ "parcial", TRUE ~ "indeterminado"))
write.csv(lp, file.path(raiz, "saida/limiares_publicados_padronizados.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat(format(Sys.time()), "L1:", nrow(lp), "linhas;", sum(lp$variavel_no_modelo == "sim"), "com variavel no modelo;", sum(!is.na(lp$chave)), "com caixa\n")

# ---------------- L2 compare with the model ----------------
rob <- read.csv(file.path(raiz, "saida/etapa3_n1_limiares_robustos.csv"))
rob_g <- rob |> mutate(id_rob = row_number()) |> separate_rows(lista, sep = ", ") |>
  mutate(grupo = sub(" \\((bentos|pelagico)\\)$", "", lista), comp = sub("^.*\\((bentos|pelagico)\\)$", "\\1", lista)) |> select(-lista)
loc <- read.csv(file.path(raiz, "saida/etapa3_n1_transf_limiares_locais.csv")) |> mutate(bacia_nome = bac_nome$bacia[match(bacia, bac_nome$codigo)])
arqs <- list.files(file.path(raiz, "intermediario/etapa3_n1"), pattern = "^(bentos|pelagico)_.*\\.rds$", full.names = TRUE)
faixa_amostrada <- bind_rows(lapply(arqs, function(a) { r <- readRDS(a); vs <- names(r$E)[!grepl("^bac_", names(r$E))]
  data.frame(comp = r$comp, grupo = r$grupo, preditor = vs, xmin = sapply(r$E[vs], min), xmax = sapply(r$E[vs], max)) }))
grupos_comp <- faixa_amostrada |> distinct(comp, grupo)
compara <- function(a, b, de, ate) { if (!length(de)) return("nao_detectado")
  if (any(a * (1 - TOL) <= ate & b * (1 + TOL) >= de)) "confirmado" else if (any(a * (1 - 2 * TOL) <= ate & b * (1 + 2 * TOL) >= de)) "proximo" else "nao_detectado" }
cmp <- bind_rows(lapply(seq_len(nrow(lp)), function(i) { x <- lp[i, ]
  gs <- if (x$grupo_projeto == "comunidade") setdiff(grupos_comp$grupo[grupos_comp$comp == x$compartimento], SEM_SINAL) else x$grupo_projeto
  bind_rows(lapply(gs, function(g) {
    base <- data.frame(id = x$id, compartimento = x$compartimento, grupo_projeto = g, comunidade = x$grupo_projeto == "comunidade", variavel_modelo = x$variavel_modelo,
      limiar_min = x$limiar_min, limiar_max = x$limiar_max, tipo_limiar = x$tipo_limiar, bacia = x$bacia, escala = x$escala, limite_igual_estrato = x$limite_igual_estrato)
    if (x$tipo_limiar == "sem_quebra") return(cbind(base, situacao = "sem_quebra", situacao_bacia = NA, limiar_modelo = NA, limiar_bacia = NA))
    if (is.na(x$variavel_modelo)) return(cbind(base, situacao = "variavel_ausente", situacao_bacia = NA, limiar_modelo = NA, limiar_bacia = NA))
    if (is.na(x$limiar_min)) return(cbind(base, situacao = "sem_valor", situacao_bacia = NA, limiar_modelo = NA, limiar_bacia = NA))
    if (x$tipo_limiar %in% TIPO_INCOMP) tipo_inc <- TRUE else tipo_inc <- FALSE
    fa <- faixa_amostrada |> filter(comp == x$compartimento, grupo == g, preditor == x$variavel_modelo)
    if (!nrow(fa)) return(cbind(base, situacao = "variavel_ausente", situacao_bacia = NA, limiar_modelo = NA, limiar_bacia = NA))
    r <- rob_g |> filter(comp == x$compartimento, grupo == g, preditor == x$variavel_modelo)
    s <- if (x$limiar_max < fa$xmin || x$limiar_min > fa$xmax) "fora_do_alcance" else compara(x$limiar_min, x$limiar_max, r$de, r$ate)
    if (tipo_inc) s <- "tipo_incompativel"
    lmod <- if (nrow(r)) r$limiar[which.min(abs(r$limiar - mean(c(x$limiar_min, x$limiar_max))))] else NA
    l <- loc |> filter(comp == x$compartimento, grupo == g, preditor == x$variavel_modelo, bacia_nome == x$bacia)
    sb <- if (x$bacia %in% c("global", "") || is.na(x$bacia)) NA else if (!nrow(l)) "sem_modelo_local" else compara(x$limiar_min, x$limiar_max, l$limiar, l$limiar)
    lb <- if (nrow(l)) l$limiar[which.min(abs(l$limiar - mean(c(x$limiar_min, x$limiar_max))))] else NA
    cbind(base, situacao = s, situacao_bacia = sb, limiar_modelo = lmod, limiar_bacia = lb) })) }))
write.csv(cmp, file.path(raiz, "saida/limiares_comparacao.csv"), row.names = FALSE, fileEncoding = "UTF-8")
options(width = 220)
print(lp |> select(id, compartimento, grupo_projeto, variavel_modelo, bacia, bacia_conferida, dados_no_obis, limite_igual_estrato), right = FALSE)
print(table(cmp$situacao)); print(table(cmp$situacao_bacia, useNA = "ifany"))
res <- cmp |> group_by(id) |> summarise(situacao = { s <- unique(situacao); ord <- c("confirmado", "proximo", "nao_detectado", "fora_do_alcance", "tipo_incompativel", "sem_valor", "variavel_ausente", "sem_quebra"); ord[min(match(s, ord))] }, .groups = "drop")
print(table(res$situacao)); cat(format(Sys.time()), "FIM\n")
