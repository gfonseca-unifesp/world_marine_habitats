# 90_relatorio.R
# Relatorio HTML de validacao visual das etapas ja executadas. Pode ser reexecutado a cada etapa nova.
# Saida: saida/relatorio_validacao.html (arquivo unico, imagens embutidas).
suppressPackageStartupMessages({library(terra); library(arrow); library(dplyr); library(jsonlite)})
raiz <- "C:/Users/model/OBIS_Tot"
source(file.path(raiz, "scripts/funcoes_grade.R"))
dir.create(file.path(raiz, "saida"), showWarnings = FALSE)
D <- list(gerado = format(Sys.time(), "%d/%m/%Y %H:%M"))
png64 <- function(arq) paste0("data:image/png;base64,", base64enc::base64encode(arq))
tmp <- tempfile(fileext = ".png")

# ---- Fase 0: manifesto e download -------------------------------------------------------------
man <- read.csv(file.path(raiz, "bruto/manifesto_s3.csv"))
man$local <- file.path(raiz, "bruto", man$key)
man$tam_local <- ifelse(file.exists(man$local), file.size(man$local), NA)
man$completo <- !is.na(man$tam_local) & man$tam_local == man$size
man$parte <- (seq_len(nrow(man)) - 1) %% 3 + 1
D$download <- list(
  arquivos = nrow(man), gb_total = round(sum(man$size) / 1024^3, 1),
  completos = sum(man$completo), gb_completos = round(sum(man$size[man$completo]) / 1024^3, 1),
  multipart = sum(grepl("-", man$etag)), listado = readLines(file.path(raiz, "bruto/manifesto_data.txt"))[1],
  partes = lapply(1:3, function(p) list(parte = p, feitos = sum(man$completo & man$parte == p),
                                        total = sum(man$parte == p),
                                        gb = round(sum(man$size[man$completo & man$parte == p]) / 1024^3, 1),
                                        gb_total = round(sum(man$size[man$parte == p]) / 1024^3, 1))))
br <- c(0, 1e5, 1e6, 1e7, 1e8, 1e9, Inf)
rot <- c("< 100 kB", "100 kB a 1 MB", "1 a 10 MB", "10 a 100 MB", "100 MB a 1 GB", "> 1 GB")
fx <- cut(man$size, br, rot, right = FALSE)
D$tamanhos <- data.frame(faixa = rot, arquivos = as.integer(table(fx)),
                         gb = round(as.numeric(tapply(man$size, fx, sum)) / 1024^3, 2))
D$tamanhos$gb[is.na(D$tamanhos$gb)] <- 0

# ---- Fase 0.3: preenchimento de campos numa amostra local ---------------------------------------
locais <- man$local[man$completo]
if (length(locais) >= 50) {
  set.seed(1); am <- sample(locais, min(300, length(locais)))
  ds <- open_dataset(am, format = "parquet", unify_schemas = FALSE)
  x <- ds |>
    mutate(lon = interpreted$decimalLongitude, ano = interpreted$date_year,
           ds = interpreted$date_start, de = interpreted$date_end, pmin = interpreted$minimumDepthInMeters, pmax = interpreted$maximumDepthInMeters,
           bat = interpreted$bathymetry, ev = interpreted$eventID, fam = interpreted$family,
           gen = interpreted$genus, apar = interpreted$samplingProtocol, est = interpreted$lifeStage,
           inc = interpreted$coordinateUncertaintyInMeters, bor = interpreted$basisOfRecord,
           lic = interpreted$license) |>
    select(dataset_id, lon, ano, ds, de, pmin, pmax, bat, ev, fam, gen, apar, est, inc, bor, lic, dropped, absence) |>
    collect()
  # Mes (ou dia) conhecido = inicio e fim da coleta no mesmo mes (ou dia); date_start/date_end em ms desde 1970
  di <- as.POSIXct(as.numeric(x$ds) / 1000, origin = "1970-01-01", tz = "UTC")
  dfim <- as.POSIXct(as.numeric(x$de) / 1000, origin = "1970-01-01", tz = "UTC")
  x$mes <- ifelse(!is.na(di) & format(di, "%Y-%m") == format(dfim, "%Y-%m"), 1, NA)
  x$dia <- ifelse(!is.na(di) & format(di, "%Y-%m-%d") == format(dfim, "%Y-%m-%d"), 1, NA)
  campo <- c(lon = "Coordenadas", ano = "Ano", mes = "Mês (início e fim no mesmo mês)", dia = "Dia (início e fim no mesmo dia)", pmin = "Profundidade mínima",
             pmax = "Profundidade máxima", bat = "Batimetria (OBIS)", ev = "eventID", fam = "Família",
             gen = "Gênero", apar = "Aparelho (samplingProtocol)", est = "Estágio de vida (lifeStage)",
             inc = "Incerteza de coordenada", lic = "Licença no registro (a do conjunto vem dos metadados)")
  cheio <- sapply(names(campo), function(k) { v <- x[[k]]; mean(!is.na(v) & as.character(v) != "") })
  D$preenchimento <- list(n_registros = nrow(x), n_conjuntos = length(unique(x$dataset_id)),
                          n_arquivos = length(am),
                          tabela = data.frame(campo = unname(campo), pct = round(100 * unname(cheio), 1)),
                          dropped = round(100 * mean(x$dropped, na.rm = TRUE), 1),
                          ausencia = round(100 * mean(x$absence, na.rm = TRUE), 1),
                          desde1993 = round(100 * mean(x$ano >= 1993, na.rm = TRUE), 1),
                          basis = as.data.frame(sort(table(x$bor), decreasing = TRUE)[1:6]))
  names(D$preenchimento$basis) <- c("tipo", "registros")
}

# ---- Fase 1.1: grade ---------------------------------------------------------------------------
D$grade <- read.csv(file.path(raiz, "camadas/grades/grades_parametros.csv"))

# ---- Fase 1.2: EBU ------------------------------------------------------------------------------
cls <- read.csv(file.path(raiz, "camadas/ebu/ebu_classes.csv"))
D$ebu_classes <- as.data.frame(table(cls$nivel)); names(D$ebu_classes) <- c("nivel", "classes")
res_ord <- c("6.25km", "12.5km", "25km", "50km", "100km")
ebu_res <- lapply(res_ord, function(r) read_parquet(file.path(raiz, sprintf("camadas/ebu/ebu_celula_%s.parquet", r))))
names(ebu_res) <- res_ord
D$ebu_pureza <- do.call(rbind, lapply(res_ord, function(r) {
  e <- ebu_res[[r]]
  data.frame(res = r, celulas = nrow(e),
             med_grupo = median(e$pureza_grupo, na.rm = TRUE), med_regiao = median(e$pureza_regiao, na.rm = TRUE),
             med_unidade = median(e$pureza_unidade, na.rm = TRUE),
             p80_grupo = round(100 * mean(e$pureza_grupo >= 0.8, na.rm = TRUE), 1),
             p80_regiao = round(100 * mean(e$pureza_regiao >= 0.8, na.rm = TRUE), 1),
             p80_unidade = round(100 * mean(e$pureza_unidade >= 0.8, na.rm = TRUE), 1))
}))
e25 <- ebu_res[["25km"]]
D$ebu_grupos <- e25 |> count(ebu_grupo, name = "celulas") |> arrange(ebu_grupo) |> as.data.frame()
geo <- c(A = "Plataforma", B = "Talude", C = "Montes submarinos e guyots", D = "Dorsais e vales de rift",
         E = "Abissal e hadal", F = "Platôs")
D$ebu_grupos$nome <- geo[D$ebu_grupos$ebu_grupo]
reg <- e25 |> count(ebu_regiao, name = "celulas") |> arrange(desc(celulas))
D$ebu_regioes <- list(n = nrow(reg), max = reg$celulas[1], min = min(reg$celulas),
                      mediana = median(reg$celulas), menos_de_10 = sum(reg$celulas < 10))
uni <- e25 |> count(ebu_unidade, name = "celulas")
D$ebu_unidades <- list(n = nrow(uni), menos_de_5 = sum(uni$celulas < 5), mediana = median(uni$celulas))

# Mapas: um painel por grupo geomorfologico (presenca), e pureza da unidade a 25 km
g25 <- ease2_grade("25km")
cel <- function(ids) { p <- do.call(rbind, strsplit(ids, "_")); cellFromRowCol(g25, as.integer(p[, 2]), as.integer(p[, 3])) }
e25$cel <- cel(e25$cell_id)
png(tmp, 1800, 1000, res = 150, bg = "white"); par(mfrow = c(2, 3), mar = c(1, 1, 2.5, 1))
for (gq in names(geo)) {
  r <- rast(g25); v <- rep(NA, ncell(r)); v[e25$cel[e25$ebu_grupo == gq]] <- 1; r[] <- v
  base <- rast(g25); vb <- rep(NA, ncell(base)); vb[e25$cel] <- 1; base[] <- vb
  plot(base, col = "#e1e0d9", legend = FALSE, axes = FALSE, box = FALSE, mar = NA)
  plot(r, col = "#2a78d6", legend = FALSE, axes = FALSE, add = TRUE)
  title(sprintf("%s · %s", gq, geo[gq]), cex.main = 2.2, font.main = 1, col.main = "#0b0b0b", line = 0.5)
}
dev.off(); D$mapa_ebu_grupos <- png64(tmp)

seq_blue <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
r <- rast(g25); v <- rep(NA, ncell(r)); v[e25$cel] <- e25$pureza_unidade; r[] <- v
png(tmp, 1800, 820, res = 150, bg = "white"); par(mar = c(1, 1, 1, 1))
plot(r, col = seq_blue, breaks = c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 1.0001), axes = FALSE, box = FALSE,
     plg = list(title = "pureza", cex = 0.9))
dev.off(); D$mapa_ebu_pureza <- png64(tmp)

# ---- Fase 1.3: EMU ------------------------------------------------------------------------------
niv <- read.csv(file.path(raiz, "camadas/emu/emu_niveis.csv"))
emu_cls <- read.csv(file.path(raiz, "camadas/emu/emu_classes.csv"))
s <- rast(file.path(raiz, "camadas/emu/emu_25km.tif")); p <- rast(file.path(raiz, "camadas/emu/emu_pureza_25km.tif"))
D$emu_prof <- do.call(rbind, lapply(seq_len(nlyr(s)), function(b) {
  v <- values(s[[b]])[, 1]; pv <- values(p[[b]])[, 1][!is.na(v)]
  data.frame(prof = niv$profundidade_m[b], celulas = sum(!is.na(v)), emus = length(unique(na.omit(v))),
             pureza_baixa = if (length(pv)) round(100 * mean(pv < 80), 2) else NA)
}))
D$emu_prof <- D$emu_prof[order(D$emu_prof$prof), ]
b0 <- which(niv$profundidade_m == 0)
sup <- values(s[[b0]])[, 1]
tab <- as.data.frame(table(sup)); names(tab) <- c("Cluster37", "celulas_superficie")
tab$Cluster37 <- as.integer(as.character(tab$Cluster37))
# celulas em que a EMU aparece em qualquer profundidade
qualquer <- sapply(emu_cls$Cluster37, function(k) sum(any(s == k)[] %in% TRUE))
emu_tab <- data.frame(Cluster37 = emu_cls$Cluster37, celulas_qualquer = qualquer) |>
  left_join(tab, by = "Cluster37") |> left_join(emu_cls, by = "Cluster37") |>
  filter(Cluster37 > 0) |> mutate(celulas_superficie = coalesce(celulas_superficie, 0L))
D$emu_tabela <- emu_tab[, c("Cluster37", "rotulo", "prof_media_m", "PRCNT_H20", "celulas_qualquer",
                            "celulas_superficie", "EMU_Name")]
# Mapa: numero de EMUs empilhadas por coluna d'agua
nemu <- app(s, function(x) length(unique(x[!is.na(x)])))
nemu[nemu == 0] <- NA
png(tmp, 1800, 820, res = 150, bg = "white"); par(mar = c(1, 1, 1, 1))
plot(nemu, col = seq_blue, breaks = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 20), axes = FALSE, box = FALSE,
     plg = list(title = "EMUs na coluna", legend = c("1", "2", "3", "4", "5", "6", "7+"), cex = 0.9))
dev.off(); D$mapa_emu_coluna <- png64(tmp)
D$emu_coluna <- as.data.frame(table(cut(values(nemu)[, 1], c(0, 1, 2, 3, 4, 5, 6, 50), c("1", "2", "3", "4", "5", "6", "7+"))))
names(D$emu_coluna) <- c("emus", "celulas")

# ---- Coerencia entre as mascaras de fundo (EBU) e de agua (EMU) a 25 km ------------------------
agua <- !is.na(sup); fundo <- seq_len(ncell(g25)) %in% e25$cel
D$mascaras <- data.frame(categoria = c("EBU e EMU", "Só EBU (fundo sem água na EMU)", "Só EMU (água sem fundo na EBU)"),
                         celulas = c(sum(agua & fundo), sum(fundo & !agua), sum(agua & !fundo)))
r <- rast(g25); v <- rep(NA, ncell(r)); v[agua & fundo] <- 1; v[fundo & !agua] <- 2; v[agua & !fundo] <- 3; r[] <- v
png(tmp, 1800, 820, res = 150, bg = "white"); par(mar = c(1, 1, 1, 1))
plot(r, col = c("#e1e0d9", "#eb6834", "#2a78d6"), type = "classes", axes = FALSE, box = FALSE,
     levels = c("EBU e EMU", "Só EBU", "Só EMU"), plg = list(cex = 0.9))
dev.off(); D$mapa_mascaras <- png64(tmp)

# ---- 0.1 verificacao e 0.2 ambiente -------------------------------------------------------------
vf <- file.path(raiz, "bruto/verificacao.csv")
if (file.exists(vf)) { v <- read.csv(vf); D$verificacao <- as.data.frame(table(status = v$status)) }
af <- file.path(raiz, "saida/ambiente.csv"); if (file.exists(af)) D$ambiente <- read.csv(af)

# ---- 1.4 GLORYS12: amplitude sazonal ------------------------------------------------------------
ag <- file.path(raiz, "camadas/glorys/amplitude_100km.parquet")
if (file.exists(ag)) {
  a <- read_parquet(ag)
  qs <- function(v) round(quantile(v, c(0.1, 0.5, 0.9), na.rm = TRUE), 2)
  vv <- c(sst_amp = "Temperatura de superfície (°C)", sss_amp = "Salinidade de superfície",
          mld_amp = "Camada de mistura (m)", tfundo_amp = "Temperatura de fundo (°C)")
  D$glorys <- do.call(rbind, lapply(names(vv), function(k) { q <- qs(a[[k]])
    data.frame(variavel = vv[[k]], p10 = q[1], mediana = q[2], p90 = q[3], celulas = sum(!is.na(a[[k]]))) }))
  g100 <- ease2_grade("100km")
  mapa_amp <- function(k, titulo) {
    p <- do.call(rbind, strsplit(a$cell_id, "_"))
    r <- rast(g100); vals <- rep(NA, ncell(r)); vals[cellFromRowCol(g100, as.integer(p[, 2]), as.integer(p[, 3]))] <- a[[k]]; r[] <- vals
    br <- unique(quantile(a[[k]], seq(0, 1, length.out = 8), na.rm = TRUE)); br[length(br)] <- br[length(br)] * 1.0001
    png(tmp, 1800, 820, res = 150, bg = "white"); par(mar = c(1, 1, 2, 1))
    plot(r, col = seq_blue[seq_len(length(br) - 1)], breaks = br, axes = FALSE, box = FALSE,
         plg = list(title = titulo, cex = 0.9)); dev.off(); png64(tmp)
  }
  D$mapa_amp_sst <- mapa_amp("sst_amp", "°C"); D$mapa_amp_tf <- mapa_amp("tfundo_amp", "°C")
}

# ---- 2.1 a 2.4 ----------------------------------------------------------------------------------
tf <- file.path(raiz, "saida/nivel1_totais.csv")
if (file.exists(tf)) D$n1_totais <- read.csv(tf)
pf <- file.path(raiz, "saida/perfil_nivel1.csv")
if (file.exists(pf)) {
  D$perfil <- read.csv(pf)
  cj <- read.csv(file.path(raiz, "intermediario/conjunto_tipo_escopo.csv"))
  D$perfil$maior_titulo <- cj$titulo[match(D$perfil$maior_conjunto, cj$dataset_id)]
  fl <- read.csv(file.path(raiz, "saida/perfil_flags.csv"))
  D$flags <- fl |> group_by(flag) |> summarise(registros = sum(registros), grupos = n()) |> arrange(desc(registros)) |> as.data.frame()
  bs <- read.csv(file.path(raiz, "saida/perfil_basis.csv"))
  D$basis <- bs |> group_by(basisOfRecord) |> summarise(registros = sum(Freq)) |> arrange(desc(registros)) |> as.data.frame()
  D$cj_tipos <- as.data.frame(table(tipo = cj$tipo)); D$cj_lic <- as.data.frame(table(licenca = cj$licenca))
  D$cj_lic$registros <- as.numeric(tapply(cj$registros, cj$licenca, sum)[as.character(D$cj_lic$licenca)])
  D$cj_incertos <- cj[cj$incerto %in% TRUE, c("dataset_id", "titulo", "tipo", "registros")] |> arrange(desc(registros)) |> head(30)
  D$cj_cpr <- cj[cj$cpr %in% TRUE, c("dataset_id", "titulo", "registros")]
  D$cj_n <- nrow(cj)
  tx <- read.csv(file.path(raiz, "saida/taxa.csv"))
  D$taxa <- tx |> group_by(grupo) |> summarise(taxons = n_distinct(aphiaid), familias = n_distinct(family[family != ""]),
                                              registros = sum(registros)) |> as.data.frame()
  D$taxa_dup <- sum(table(tx$aphiaid) > 1)
}

# ---- 3 Sondagem e janelas ------------------------------------------------------------------------
sf <- file.path(raiz, "saida/sondagem_resolucao.csv")
if (file.exists(sf)) {
  D$sond <- read.csv(sf)
  D$autoc <- read.csv(file.path(raiz, "saida/sondagem_autocorrelacao.csv")) |> distinct(grupo, alcance_km, patamar)
  D$curva <- read.csv(file.path(raiz, "saida/sondagem_curva_esforco.csv")) |> filter(nivel %in% c("grupo", "emu"))
  jp <- read.csv(file.path(raiz, "saida/janelas_pontuadas.csv")); js <- read.csv(file.path(raiz, "saida/janelas_selecionadas.csv"))
  D$jan_n <- list(com_dados = nrow(jp), passam = sum(jp$passa_filtro))
  D$jan_sel <- js[, c("ordem", "janela", "lon", "lat", "faixa_lat", "analoga", "celulas_serie", "pares_fronteira",
                      "conjuntos", "pct_maior_conjunto", "ebu_grupos", "emus", "pontuacao", "motivo")]
  D$jan_top <- head(jp[jp$passa_filtro, c("janela", "lon", "lat", "analoga", "celulas_serie", "pares_fronteira",
                                         "conjuntos", "pct_maior_conjunto", "pontuacao")], 20)
  D$jan_ana <- jp[!is.na(jp$analoga), c("analoga", "passa_filtro", "celulas_serie", "pares_fronteira", "conjuntos", "pontuacao")]
  # mapa: pontuacao das janelas que passam, selecionadas em destaque
  e25 <- ext(ease2_grade("25km")); lw <- 20 * EASE2_BASE_M        # janelas de 20 x 20 celulas (a ultima linha/coluna fica parcial)
  g500 <- rast(ncols = 70, nrows = 30, xmin = e25$xmin, xmax = e25$xmin + 70 * lw, ymin = e25$ymax - 30 * lw,
               ymax = e25$ymax, crs = "EPSG:6933")
  rcj <- function(j) { p <- do.call(rbind, strsplit(j, "_")); cellFromRowCol(g500, as.integer(p[, 1]), as.integer(p[, 2])) }
  r <- rast(g500); v <- rep(NA, ncell(r)); v[rcj(jp$janela)] <- 0; ok2 <- jp$passa_filtro
  v[rcj(jp$janela[ok2])] <- jp$pontuacao[ok2]; r[] <- v
  pol <- as.polygons(rast(g500, vals = ifelse(seq_len(ncell(g500)) %in% rcj(js$janela), 1, NA)), dissolve = FALSE)
  png(tmp, 1800, 820, res = 150, bg = "white"); par(mar = c(1, 1, 1, 1))
  plot(r, col = c("#e1e0d9", seq_blue[2:7]), breaks = c(-0.5, 0.5, 5, 6, 7, 8, 9, 12), axes = FALSE, box = FALSE,
       plg = list(title = "pontuação", legend = c("com dados, não passa", "< 5", "5–6", "6–7", "7–8", "8–9", "≥ 9"), cex = 0.8))
  plot(pol, add = TRUE, border = "#eb6834", lwd = 2.5)
  ll <- project(cbind(js$lon, js$lat), from = "EPSG:4326", to = "EPSG:6933")
  text(ll[, 1], ll[, 2] + 4e5, js$ordem, col = "#0b0b0b", cex = 1, font = 2)
  dev.off(); D$mapa_janelas <- png64(tmp)
}

# ---- 4 Curadoria e 5 Tratamento ------------------------------------------------------------------
cf <- file.path(raiz, "saida/cascata_filtros_n1.csv")
if (file.exists(cf)) {
  cas <- read.csv(cf)
  D$cascata <- cas
  D$cascata_tot <- cas |> group_by(filtro) |> summarise(entram = sum(entram), saem = sum(saem), restam = sum(restam)) |>
    mutate(ordem = as.numeric(sub(" .*", "", filtro))) |> arrange(ordem) |> as.data.frame()
  q <- open_dataset(file.path(raiz, "intermediario/quarentena/n1"))
  D$quarentena <- q |> count(filtro, motivo) |> collect() |> arrange(desc(n)) |> as.data.frame()
  rt <- open_dataset(file.path(raiz, "intermediario/n1_curado")) |> group_by(grupo, compartimento) |>
    summarise(retidos = n(), so_esforco = sum(so_esforco), sem_prof = sum(sem_profundidade), relidas = sum(data_relida),
              A = sum(nivel_certeza == "A"), B = sum(nivel_certeza == "B"), C = sum(nivel_certeza == "C")) |> collect()
  ex <- read.csv(file.path(raiz, "saida/nivel1_totais.csv"))
  D$retidos <- rt |> left_join(ex |> select(grupo, extraidos = registros), by = "grupo") |> arrange(desc(retidos)) |> as.data.frame()
}
tf5 <- file.path(raiz, "saida/tratamento_n1_resumo.csv")
if (file.exists(tf5)) {
  D$trat <- read.csv(tf5)
  es <- read_parquet(file.path(raiz, "saida/esforco_nivel1.parquet"))
  D$esf <- es |> group_by(grupo, compartimento, cell_25) |> summarise(eventos = sum(eventos), .groups = "drop") |>
    group_by(grupo, compartimento) |> summarise(celulas = n(), celulas_10ev = sum(eventos >= 10), eventos = sum(eventos), .groups = "drop") |>
    as.data.frame()
  D$conj_grupo_n <- nrow(read.csv(file.path(raiz, "saida/conjunto_grupo.csv")))
}

# ---- 6 Tabelas de certeza ----------------------------------------------------------------------
vf6 <- file.path(raiz, "saida/validacao_n1.csv")
if (file.exists(vf6)) {
  D$valid6 <- read.csv(vf6)
  tb <- read_parquet(file.path(raiz, "saida/tabela_certa_bentos_n1.parquet"), col_select = c("grupo", "nivel_certeza", "cell_id"))
  tp <- read_parquet(file.path(raiz, "saida/tabela_certa_pelagico_n1.parquet"), col_select = c("grupo", "nivel_certeza", "cell_id"))
  D$tab6 <- bind_rows(tb |> mutate(tabela = "bentos"), tp |> mutate(tabela = "pelágico")) |> group_by(tabela, grupo) |>
    summarise(registros = n(), celulas = n_distinct(cell_id), A = sum(nivel_certeza == "A"), B = sum(nivel_certeza == "B"),
              C = sum(nivel_certeza == "C"), .groups = "drop") |> as.data.frame()
  rm(tb, tp)
  ts <- read.csv(file.path(raiz, "saida/testabilidade_n1.csv"))
  D$test6 <- ts |> distinct(nivel, grupo, fatia, classes_aptas, contraste_viavel) |> group_by(nivel, grupo) |>
    summarise(fatias = n(), fatias_viaveis = sum(contraste_viavel == "sim"), classes_aptas_max = max(classes_aptas), .groups = "drop") |>
    as.data.frame()
  D$inc6 <- length(list.files(file.path(raiz, "saida/incidencia_n1")))
}

# ---- 7 Nivel 2 -----------------------------------------------------------------------------------
vf7 <- file.path(raiz, "saida/validacao_n2.csv")
if (file.exists(vf7)) {
  D$valid7 <- read.csv(vf7)
  D$n2_tot <- read.csv(file.path(raiz, "saida/nivel2_totais.csv")) |> group_by(regiao) |>
    summarise(extraidos = sum(registros), grupos = n(), conjuntos = max(conjuntos)) |> as.data.frame()
  c2 <- read.csv(file.path(raiz, "saida/cascata_filtros_n2.csv"))
  D$n2_ret <- c2 |> filter(filtro == "7.3 compartimento") |> group_by(regiao) |> summarise(retidos = sum(restam)) |> as.data.frame()
  D$n2_comp <- c2 |> filter(filtro == "7.3 compartimento") |> group_by(regiao) |> summarise(entram = sum(entram), sem_decisao = sum(saem)) |>
    mutate(pct = round(100 * sem_decisao / entram, 1)) |> as.data.frame()
  t2 <- read.csv(file.path(raiz, "saida/testabilidade_n2.csv"))
  D$n2_test <- t2 |> filter(fatia %in% c("estatico", paste0("M", sprintf("%02d", 1:12)))) |>
    group_by(regiao, nivel) |> summarise(combinacoes = n_distinct(paste(grupo_ecologico, fatia)),
      viaveis = n_distinct(paste(grupo_ecologico, fatia)[contraste_viavel == "sim"]), .groups = "drop") |> as.data.frame()
}

af <- file.path(raiz, "entrada/aparelho_conjunto.csv")
if (file.exists(af)) D$apar_cj <- read.csv(af) |> count(fonte, aparelho_conjunto) |> as.data.frame()

# ---- 8 Etapa 1 no nivel 1 ----------------------------------------------------------------------
pf8 <- file.path(raiz, "saida/etapa1_n1_permanova.csv")
if (file.exists(pf8)) {
  pm <- read.csv(pf8)
  D$e1_perm <- pm |> filter(is.na(indice) | indice == "sim") |> select(any_of(c("comp", "grupo", "nivel", "fatia", "sem_mlpa", "unidades", "classes",
    "taxons", "R2_global", "unidades_mistas", "blocos_mistos", "classes_mistas", "R2_dentro", "R2_dentro_min", "R2_dentro_max", "p_dentro_max",
    "R2_bloco", "disp_p", "omega2_global", "omega2_dentro", "omega2_dentro_min", "omega2_dentro_max", "n_med", "omega2_eq",
    "omega2_eq_q05", "omega2_eq_q95", "poder_eq", "p_eq_med", "contraste"))) |> as.data.frame()
  D$e1_sor <- pm |> filter(indice == "sor") |> select(comp, grupo, nivel, fatia, sem_mlpa, R2_global, R2_dentro, p_dentro_max) |> as.data.frame()
  iv <- read.csv(file.path(raiz, "saida/etapa1_n1_indval.csv"))
  D$e1_iv <- iv |> group_by(comp, grupo, nivel, fatia, sem_mlpa, classe) |>
    summarise(taxons = n(), indicadores = sum(p_holm < 0.05), indicadores_bh = sum(if ("p_bh" %in% names(iv)) p_bh < 0.05 else NA),
              melhores = paste(head(taxon[order(-stat)][(if ("p_bh" %in% names(iv)) p_bh else p_holm)[order(-stat)] < 0.05], 3), collapse = ", "),
              .groups = "drop") |> as.data.frame()
}

# ---- 8 Habitats avaliados e distintos (pares), niveis 1 e 2 ------------------------------------
# status do habitat num contraste: sem par testavel | distinto de todos | distinto de parte | nao distinto
status_hab <- function(hb, pr, chave) {
  h <- hb |> mutate(habitat = as.character(habitat))
  if (!is.null(pr) && nrow(pr)) {
    pr <- pr |> mutate(classe_a = as.character(classe_a), classe_b = as.character(classe_b))
    lado <- bind_rows(pr |> mutate(habitat = classe_a, outro = classe_b), pr |> mutate(habitat = classe_b, outro = classe_a)) |>
      group_by(across(all_of(c(chave, "habitat")))) |>
      summarise(pares = sum(!is.na(separa)), distinto_de = sum(separa %in% TRUE),
                distinto_lista = paste(sort(outro[separa %in% TRUE]), collapse = ", "),
                igual_lista = paste(sort(outro[separa %in% FALSE]), collapse = ", "), .groups = "drop")
    h <- h |> left_join(lado, by = c(chave, "habitat"))
  } else { h$pares <- NA_integer_; h$distinto_de <- NA_integer_ }
  h |> mutate(pares = coalesce(pares, 0L), distinto_de = coalesce(distinto_de, 0L),
              status = ifelse(pares == 0, "sem par testável", ifelse(distinto_de == pares, "distinto de todos",
                       ifelse(distinto_de > 0, "distinto de parte", "não distinto"))))
}
emu_cl <- read.csv(file.path(raiz, "camadas/emu/emu_classes.csv"))
D$emu_nomes <- setNames(as.list(paste0(emu_cl$EMU_Name, " (", emu_cl$prof_media_m, " m)")), emu_cl$Cluster37)
k1 <- c("comp", "grupo", "nivel", "fatia", "sem_mlpa")
if (file.exists(file.path(raiz, "saida/etapa1_n1_pares.csv"))) {
  h1 <- read.csv(file.path(raiz, "saida/etapa1_n1_habitats.csv")); p1 <- read.csv(file.path(raiz, "saida/etapa1_n1_pares.csv"))
  D$hab1 <- status_hab(h1, p1, k1) |> filter(!sem_mlpa) |> as.data.frame()
  D$par1 <- p1 |> filter(!sem_mlpa, !is.na(separa)) |> select(grupo, nivel, fatia, classe_a, classe_b, unidades_a, unidades_b, blocos, omega2, poder_eq, separa) |> as.data.frame()
}
k2 <- c("regiao", "nivel", "grupo_ecologico", "fatia", "res_escolhida")
if (file.exists(file.path(raiz, "saida/etapa1_n2_b100_pares.csv"))) {
  h2 <- read.csv(file.path(raiz, "saida/etapa1_n2_b100_habitats.csv")); p2p <- read.csv(file.path(raiz, "saida/etapa1_n2_b100_pares.csv"))
  D$hab2 <- status_hab(h2, p2p, k2) |> as.data.frame()
  D$par2 <- p2p |> filter(!is.na(separa)) |> select(regiao, nivel, grupo_ecologico, fatia, classe_a, classe_b, unidades_a, unidades_b, blocos, omega2, poder_eq, separa) |> as.data.frame()
}

# ---- 8 Etapa 1 no nivel 2 ----------------------------------------------------------------------
mag <- function(w) as.character(cut(w, c(-Inf, 0.01, 0.06, 0.14, Inf), c("desprezível", "pequeno", "médio", "grande"), right = FALSE))
f2 <- list.files(file.path(raiz, "saida"), pattern = "^etapa1_n2_b[0-9]+_permanova\\.csv$", full.names = TRUE)
if (length(f2)) {
  p2 <- bind_rows(lapply(f2, read.csv)) |> filter(is.na(indice) | indice == "sim")
  t2 <- p2 |> filter(contraste == "testado") |> mutate(separa = omega2_dentro >= 0.01 & poder_eq >= 0.8)
  D$e2_res <- t2 |> group_by(bloco_km, nivel) |> summarise(contrastes = n(), omega2_med = median(omega2_dentro),
    omega2_q25 = quantile(omega2_dentro, 0.25), omega2_q75 = quantile(omega2_dentro, 0.75), separa = sum(separa),
    poder_med = median(poder_eq), unidades_med = median(unidades_mistas), .groups = "drop") |>
    left_join(p2 |> group_by(bloco_km, nivel) |> summarise(total = n(), sem_teste = sum(contraste != "testado"), .groups = "drop"), by = c("bloco_km", "nivel")) |> as.data.frame()
  D$e2_reg <- t2 |> filter(bloco_km == 100) |> group_by(regiao, nivel) |> summarise(contrastes = n(), omega2_med = median(omega2_dentro), separa = sum(separa), .groups = "drop") |> as.data.frame()
  D$e2_tab <- t2 |> filter(bloco_km == 100) |> mutate(efeito = mag(omega2_dentro)) |>
    select(regiao, nivel, grupo_ecologico, fatia, res_escolhida, unidades_mistas, blocos_mistos, classes_mistas, R2_dentro, omega2_dentro, efeito,
           omega2_eq, poder_eq, R2_bloco, disp_p, separa) |> arrange(regiao, nivel, grupo_ecologico, fatia) |> as.data.frame()
}

# ---- 8.2 Etapa 2: diagnostico -------------------------------------------------------------------
le_s <- function(f) { a <- file.path(raiz, "saida", f); if (file.exists(a)) read.csv(a) else NULL }
D$e2_gen <- le_s("etapa2_genero_comparacao.csv")
for (nv in 1:2) {
  cg <- le_s(sprintf("etapa2_n%d_classif.csv", nv)); if (is.null(cg)) next
  D[[sprintf("e2c%d", nv)]] <- cg
  pc <- le_s(sprintf("etapa2_n%d_classif_classe.csv", nv))
  D[[sprintf("e2cc%d", nv)]] <- pc |> select(any_of(c("comp", "grupo", "regiao", "nivel", "grupo_ecologico", "fatia", "classe", "unidades", "sensibilidade",
    "precisao", "confunde_com", "dispersao", "k", "silhueta", "coerencia", "coerencia_esperada", "tamanhos", "descontinuidade")))
  D[[sprintf("e2eq%d", nv)]] <- le_s(sprintf("etapa2_n%d_equivalencias.csv", nv))
  fu <- le_s(sprintf("etapa2_n%d_fusoes.csv", nv)); if (!is.null(fu) && nrow(fu)) D[[sprintf("e2fu%d", nv)]] <- fu
  # mapa do acerto da classificacao (bloco deixado de fora)
  un <- read_parquet(file.path(raiz, sprintf("saida/etapa2_n%d_unidades.parquet", nv))) |> filter(!is.na(prevista))
  if (nv == 1) {
    g100 <- ease2_grade("100km"); xy <- project(cbind(un$lon, un$lat), from = "EPSG:4326", to = "EPSG:6933")
    un$c100 <- cellFromXY(g100, xy)
    base <- rast(g100); vb <- rep(NA, ncell(base)); rc <- rowColFromCell(g25, which(!is.na(sup)))
    vb[cellFromRowCol(g100, (rc[, 1] - 1) %/% 4 + 1, (rc[, 2] - 1) %/% 4 + 1)] <- 1; base[] <- vb
    mapas <- list()
    for (nvl in c("EBU grupo", "EBU região", "EMU")) {
      a <- un |> filter(nivel == nvl) |> group_by(c100) |> summarise(acerto = mean(correta)); r <- rast(g100); v <- rep(NA, ncell(r)); v[a$c100] <- a$acerto; r[] <- v
      png(tmp, 1800, 820, res = 150, bg = "white"); par(mar = c(1, 1, 2, 1))
      plot(base, col = "#e1e0d9", legend = FALSE, axes = FALSE, box = FALSE, main = paste0(nvl, ": acerto (bloco deixado de fora)"))
      plot(r, col = c("#d03b3b", "#eb6834", "#fab219", "#9fd39a", "#0ca30c"), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0001), add = TRUE,
           plg = list(title = "acerto", cex = 0.85, legend = c("0–0,2", "0,2–0,4", "0,4–0,6", "0,6–0,8", "0,8–1")))
      dev.off(); mapas[[nvl]] <- png64(tmp)
    }
    D$e2_mapas1 <- mapas
  } else {
    regs <- sort(unique(un$regiao)); png(tmp, 2000, 1000, res = 150, bg = "white"); par(mfrow = c(2, 4), mar = c(2, 2, 2, 1))
    for (rg in regs) { a <- un |> filter(regiao == rg) |> group_by(cell, lon, lat) |> summarise(acerto = mean(correta), .groups = "drop")
      cor <- c("#d03b3b", "#eb6834", "#fab219", "#9fd39a", "#0ca30c")[cut(a$acerto, c(-0.01, 0.2, 0.4, 0.6, 0.8, 1), labels = FALSE)]
      plot(a$lon, a$lat, pch = 15, cex = 0.6, col = cor, asp = 1 / cos(mean(a$lat) * pi / 180), main = rg, xlab = "", ylab = "", cex.main = 0.9, cex.axis = 0.7) }
    dev.off(); D$e2_mapa2 <- png64(tmp)
  }
}

# ---- Sintese bottom-up (85_sintese.R) -----------------------------------------------------------
sf8 <- file.path(raiz, "saida/sintese_habitats_n1.csv")
if (file.exists(sf8)) {
  D$sint <- read.csv(sf8)
  mp <- c("EBU grupo" = "ebu_grupo_25km", "EBU região" = "ebu_regiao_25km", "EMU superfície" = "emu_superficie_25km", "EMU 1.000 m" = "emu_1000m_25km")
  D$sint_mapas <- lapply(mp, function(a) { f <- file.path(raiz, sprintf("saida/mapa_habitats_%s.png", a)); if (file.exists(f)) png64(f) })
  # reconciliacao: da contagem de habitats (status) as classes de cada mapa
  D$sint_rec <- bind_rows(lapply(names(mp), function(k) {
    nv <- sub(" (superfície|1.000 m)$", "", k); s <- D$sint |> filter(nivel == nv); a <- s |> filter(!grepl("^stand-by", status))
    r <- rast(file.path(raiz, sprintf("saida/habitats_final_%s.tif", mp[[k]]))); l <- read.csv(file.path(raiz, sprintf("saida/habitats_final_%s_legenda.csv", mp[[k]])))
    v <- unique(values(r)[, 1]); v <- v[!is.na(v)]
    hab_mapa <- if (nv == "EMU") { e <- rast(file.path(raiz, "camadas/emu/emu_25km.tif")); camada <- if (grepl("1.000", k)) e[["prof_1000m"]] else e[[nlyr(e)]]
      as.character(sort(unique(na.omit(values(camada)[, 1])))) } else s$habitat
    data.frame(mapa = k, habitats_camada = nrow(s), standby = nrow(s) - nrow(a), testados = nrow(a),
      sozinhos = sum(is.na(a$fusao)), fundidos = sum(!is.na(a$fusao)), uniões = length(unique(na.omit(a$fusao))),
      classes_finais = length(unique(a$classe_final)), habitats_no_mapa = length(hab_mapa),
      classes_no_mapa = sum(v %in% l$valor[!l$standby]), standby_no_mapa = sum(s$habitat[grepl("^stand-by", s$status)] %in% hab_mapa))
  }))
  D$sint_dest <- D$sint |> mutate(destino_tipo = ifelse(grepl("^fundido", destino), "fundido", destino)) |> count(nivel, status, destino_tipo) |> as.data.frame()
  D$sint_leg <- bind_rows(lapply(names(mp)[1:3], function(k) { f <- file.path(raiz, sprintf("saida/habitats_final_%s_legenda.csv", mp[[k]]))
    if (file.exists(f)) cbind(camada = sub(" superfície", "", k), read.csv(f)) }))
}

# ---- 8.3 Etapa 3: Gradient Forest e GDM (nivel 1) ----------------------------------------------
D$e3_un <- le_s("etapa3_n1_unidades.csv")
# mapas das observacoes que entram na etapa 3 (unidades de 87_preditores.R): numero de grupos taxonomicos com unidade na celula,
# agregado a 100 km (maximo) para ficar visivel; oceano do GLORYS em cinza claro
dir_e3 <- file.path(raiz, "intermediario/etapa3_n1"); D$e3_obs <- list(); D$e3_obs_n <- list()
oce <- !is.na(rast(file.path(raiz, "camadas/glorys/mes/glorys_01_25km.tif"))[["prof_glorys"]])
for (cp in c("bentos", "pelagico")) {
  fs <- list.files(dir_e3, pattern = paste0("^", cp, "_.*\\.rds$"), full.names = TRUE); if (!length(fs)) next
  ob <- bind_rows(lapply(fs, function(a) { r <- readRDS(a); data.frame(grupo = r$grupo, cel = cellFromXY(g25, cbind(r$meta$x, r$meta$y))) }))
  ng <- ob |> distinct(grupo, cel) |> count(cel); r <- rast(g25); v <- rep(NA_real_, ncell(r)); v[ng$cel] <- ng$n; r[] <- v
  r100 <- aggregate(r, 4, max, na.rm = TRUE); o100 <- aggregate(oce, 4, max, na.rm = TRUE); ngr <- length(fs)
  cores <- hcl.colors(ngr + 1, "YlOrRd", rev = TRUE)[-1]   # sem o amarelo mais claro, que some no cinza
  png(tmp, 2000, 1080, res = 150, bg = "white"); layout(matrix(1:2), heights = c(1000, 80)); par(mar = c(0.5, 1, 2, 1))
  plot(o100, col = c("white", "#e4e4e4"), breaks = c(-0.5, 0.5, 1.5), legend = FALSE, axes = FALSE, box = FALSE,
       main = sprintf("Observações na análise, %s: %s unidades em %s células de 25 km", if (cp == "bentos") "bentos" else "pelágico",
                      format(nrow(ob), big.mark = "."), format(length(unique(ob$cel)), big.mark = ".")))
  plot(r100, col = cores, breaks = c(0.5, seq_len(ngr) + 0.5), legend = FALSE, axes = FALSE, box = FALSE, add = TRUE)
  par(mar = c(0, 1, 0, 1)); plot.new()
  legend("top", legend = seq_len(ngr), fill = cores, horiz = TRUE, bty = "n", cex = 1.1, title = "grupos taxonômicos com dados na célula (100 km)")
  dev.off(); D$e3_obs[[cp]] <- png64(tmp)
  D$e3_obs_n[[cp]] <- list(unidades = nrow(ob), celulas = length(unique(ob$cel)), grupos = ngr)
}
gi <- le_s("etapa3_n1_gf_importancia.csv")
if (!is.null(gi)) {
  D$e3_imp <- gi |> mutate(preditor = ifelse(grepl("^bac_", preditor), "bacia", preditor)) |> group_by(comp, grupo, preditor) |>
    summarise(importancia = sum(importancia), .groups = "drop") |> group_by(comp, grupo) |> mutate(parcela = importancia / sum(importancia)) |> ungroup() |> as.data.frame()
  D$e3_r2 <- gi |> distinct(comp, grupo, unidades, taxons, especies_ajustadas, r2_medio, r2_total) |> as.data.frame()
  D$e3_lim <- le_s("etapa3_n1_gf_limiares.csv"); D$e3_rob <- le_s("etapa3_n1_limiares_robustos.csv")
  gd <- le_s("etapa3_n1_gdm.csv")
  if (!is.null(gd)) D$e3_gdm <- gd |> group_by(comp, grupo, preditor) |> summarise(importancia = mean(importancia), deviance = mean(deviance_explicada), .groups = "drop") |>
    group_by(comp, grupo) |> mutate(parcela = importancia / sum(importancia)) |> ungroup() |> as.data.frame()
  D$e3_auc <- le_s("etapa3_n1_auc.csv"); D$e3_pg <- le_s("etapa3_n1_pares_gf.csv")
  ha <- le_s("etapa3_n1_habitats_ambiente.csv"); if (!is.null(ha)) D$e3_ha <- ha |> mutate(across(where(is.numeric), ~ round(.x, 2)))
  # curvas de importancia acumulada: um painel por preditor, uma linha por grupo
  cv <- le_s("etapa3_n1_gf_curvas.csv"); cores_g <- c("#2a78d6", "#eb6834", "#1baf7a", "#d03b3b", "#8a5cd1", "#fab219", "#5c5c5c")
  D$e3_curvas <- list()
  for (cp in c("bentos", "pelagico")) {
    c1 <- cv |> filter(comp == cp); if (!nrow(c1)) next
    vs <- gi |> filter(comp == cp) |> group_by(preditor) |> summarise(i = sum(importancia)) |> arrange(desc(i)) |> pull(preditor)
    vs <- vs[vs %in% c1$preditor]   # as bac_ entram na importancia mas nao tem curva
    gs <- sort(unique(c1$grupo)); nc <- 4; nr <- ceiling(length(vs) / nc)
    png(tmp, 2000, 420 * nr + 120, res = 150, bg = "white"); par(mfrow = c(nr, nc), mar = c(3.2, 3, 2, 0.5), mgp = c(1.9, 0.6, 0), oma = c(3, 0, 0, 0))
    for (v in vs) { c2 <- c1 |> filter(preditor == v)
      plot(NA, xlim = range(c2$x), ylim = c(0, max(c2$y)), xlab = v, ylab = "importância acumulada", main = v, cex.main = 1, cex.lab = 0.8, cex.axis = 0.75)
      for (j in seq_along(gs)) { c3 <- c2[c2$grupo == gs[j], ]; if (nrow(c3)) lines(c3$x, c3$y, col = cores_g[j], lwd = 2) }
      l <- D$e3_rob; if (!is.null(l)) abline(v = l$limiar[l$preditor == v], lty = 2, col = "#898781") }
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot.new()
    legend("bottom", legend = gs, col = cores_g[seq_along(gs)], lwd = 3, horiz = TRUE, bty = "n", cex = 0.85)
    dev.off(); D$e3_curvas[[cp]] <- png64(tmp)
  }
}

# ---- 3.5 Habitats novos -------------------------------------------------------------------------
hv <- le_s("etapa3_n1_habitats_novos_validacao.csv")
if (!is.null(hv)) {
  D$hn_val <- hv; D$hn_conc <- le_s("etapa3_n1_habitats_novos_concordancia.csv"); D$hn_fid <- le_s("etapa3_n1_habitats_novos_fidelidade.csv")
  D$hn_perfil <- le_s("etapa3_n1_habitats_novos_perfil.csv") |> mutate(across(where(is.numeric), ~ round(.x, 2)))
  D$hn_regras <- le_s("etapa3_n1_habitats_novos_regras.csv"); D$hn_fus <- le_s("etapa3_n1_habitats_novos_fusoes.csv")
  D$hn_pares <- le_s("etapa3_n1_habitats_novos_pares.csv") |> count(comp, consenso) |> as.data.frame()
  D$hn_leg <- le_s("habitats_novos_legenda.csv"); D$hn_bac <- if (file.exists(file.path(raiz, "saida/mapa_bacias.png"))) png64(file.path(raiz, "saida/mapa_bacias.png"))
  mp <- c("Bentos (25 km)" = "bentos_25km", "Pelágico, superfície, janeiro" = "pelagico_0m_jan", "Pelágico, superfície, julho" = "pelagico_0m_jul",
          "Pelágico, 200 m, janeiro" = "pelagico_200m_jan", "Pelágico, 200 m, julho" = "pelagico_200m_jul")
  D$hn_mapas <- lapply(mp, function(a) { f <- file.path(raiz, sprintf("saida/mapa_habitats_novos_%s.png", a)); if (file.exists(f)) png64(f) })
  s1 <- le_s("etapa3_n1_gf_importancia_comsfundo.csv")   # Gate 7 (com a salinidade de fundo); o principal ficou sem ela (Gate 8)
  if (!is.null(s1)) { s1 <- s1 |> filter(comp == "bentos")
    D$hn_sens <- s1 |> distinct(grupo, r2_com = r2_medio) |> left_join(gi |> filter(comp == "bentos") |> distinct(grupo, r2_sem = r2_medio), by = "grupo") |>
      left_join(s1 |> group_by(grupo) |> summarise(parcela_sfundo = sum(importancia[preditor == "sfundo_med"]) / sum(importancia)), by = "grupo") |> as.data.frame() }
  D$tr_aj <- le_s("etapa3_n1_transf_ajustes.csv"); D$tr_rob <- le_s("etapa3_n1_transf_robustos.csv")
  D$saz_var <- le_s("etapa3_n1_sazonal_variacao.csv"); D$saz_est <- le_s("etapa3_n1_sazonal_estabilidade.csv"); D$saz_fauna <- le_s("etapa3_n1_sazonal_fauna.csv")
  sl <- le_s("etapa3_n1_sazonal_limiares.csv"); if (!is.null(sl)) D$saz_lim <- sl |> group_by(variavel_mensal, faixa) |>
    summarise(limiares = n(), frac_max = max(frac_cruza), frac_med = mean(frac_cruza), .groups = "drop") |> as.data.frame()
  f <- file.path(raiz, "saida/mapa_bentos_meses_fora.png"); if (file.exists(f)) D$saz_mapa <- png64(f)
  D$e3_bac <- gi |> mutate(tipo = ifelse(grepl("^bac_", preditor), "bacia", "ambiente")) |> group_by(comp, grupo, tipo) |> summarise(i = sum(importancia), .groups = "drop") |>
    group_by(comp, grupo) |> mutate(parcela = i / sum(i)) |> filter(tipo == "bacia") |> select(comp, grupo, parcela_bacia = parcela) |> as.data.frame()
}

# ---- Status das tarefas ------------------------------------------------------------------------
ok <- function(f) if (file.exists(file.path(raiz, f))) "concluída" else "pendente"
D$tarefas <- data.frame(
  tarefa = c("0.1 Cópia local do OBIS", "0.2 Ambiente", "0.3 Esquema e colunas", "1.1 Grades", "1.2 EBU", "1.3 EMU",
             "1.4 Amplitude sazonal (GLORYS12)", "2.1 Extração dos táxons focais", "2.2 Perfil", "2.3 Tabela de táxons",
             "2.4 Tipo e escopo dos conjuntos", "3 Sondagem", "Gate 1", "4 Curadoria do nível 1", "5 Tratamento do nível 1", "6 Tabelas de certeza do nível 1", "Gate 2", "7 Nível 2 (8 regiões)", "Gate 3", "Gate 3 item 4: aparelho do conjunto", "8.1 Etapa 1 no nível 1", "8.1 Etapa 1 no nível 2", "8.2 Etapa 2 (diagnóstico) e mapa final", "8.3 Etapa 3 (Gradient Forest + GDM)"),
  status = c(if (file.exists(vf)) "concluída (verificada)" else "em andamento", ok("saida/ambiente.csv"),
             "parcial (amostra; eMoF pendente)", "concluída", "concluída", "concluída",
             ok("camadas/glorys/amplitude_100km.parquet"), ok("saida/nivel1_totais.csv"), ok("saida/perfil_nivel1.csv"),
             ok("saida/taxa.csv"), ok("intermediario/conjunto_tipo_escopo.csv"),
             if (file.exists(file.path(raiz, "saida/janelas_selecionadas.csv"))) "concluída" else "pendente",
             "concluída (decisões em LEIA-ME_curadoria.md)", ok("saida/cascata_filtros_n1.csv"), ok("saida/tratamento_n1_resumo.csv"), ok("saida/validacao_n1.csv"),
             if (file.exists(vf6)) "concluída (decisões em LEIA-ME_curadoria.md)" else "pendente",
             ok("saida/validacao_n2.csv"), if (file.exists(vf7)) "concluída (decisões em LEIA-ME_curadoria.md)" else "pendente",
             ok("entrada/aparelho_conjunto.csv"), if (file.exists(pf8)) "concluída (Gate 4 em LEIA-ME_curadoria.md)" else "pendente",
             if (length(f2)) "concluída (Gate 5 em LEIA-ME_curadoria.md)" else "pendente",
             if (file.exists(file.path(raiz, "saida/etapa2_n1_classif.csv"))) "concluída (Gate 6 em LEIA-ME_curadoria.md)" else "em andamento",
             if (file.exists(file.path(raiz, "saida/etapa3_n1_auc.csv"))) "aguardando decisões" else "em andamento"),
  saida = c("bruto/occurrence, bruto/verificacao.csv", "saida/ambiente.csv", "amostra local (seção 0.3)",
            "camadas/grades/grades_parametros.csv", "camadas/ebu/ebu_celula_<res>.parquet",
            "camadas/emu/emu_<res>.tif, emu_pureza_<res>.tif", "camadas/glorys/amplitude_<res>.parquet",
            "intermediario/nivel1_extracao/", "saida/perfil_*.csv", "saida/taxa.csv",
            "intermediario/conjunto_tipo_escopo.csv, bruto/obis_conjuntos.csv",
            "saida/sondagem_*.csv, saida/janelas_*.csv", "entrada/regioes_nivel2.csv",
            "intermediario/n1_curado/, quarentena/n1/, saida/cascata_filtros_n1.csv",
            "intermediario/n1_tratado/, saida/esforco_nivel1.parquet, saida/conjunto_grupo.csv",
            "saida/tabela_certa_*_n1.parquet, validacao_n1.csv, testabilidade_n1.csv, incidencia_n1/", "-",
            "saida/tabela_certa_*_n2.parquet, validacao_n2.csv, testabilidade_n2.csv, resolucao_n2.csv", "-",
            "entrada/aparelho_conjunto.csv", "saida/etapa1_n1_permanova.csv, saida/etapa1_n1_indval.csv",
            "saida/etapa1_n2_b100|b50_permanova.csv e _indval.csv",
            "saida/etapa2_*.csv, saida/sintese_habitats_n1.csv, saida/habitats_final_*",
            "camadas/glorys/mes/, saida/etapa3_n1_*.csv"))

# ---- Gate 9: mapa unido do bentos (91_mapa_unido_bentos.R) ---------------------------------------
fmu <- file.path(raiz, "saida/mapa_unido_bentos_classes.csv")
if (file.exists(fmu) && "valor" %in% names(read.csv(fmu, nrows = 1))) {
  mu <- read.csv(fmu)
  D$mu_classes <- mu |> select(any_of(c("codigo", "valor", "celulas", "faixa", "nivel3_rotulo", "regioes", "hb", "evidencia", "evidencia_nivel3",
    "protegida", "sazonalidade", "standby", "cor", "candidatos", "classe_mae", "referencia")), starts_with("tol_")) |> as.data.frame()
  D$mu_niveis <- le_s("mapa_unido_bentos_niveis.csv"); D$mu_aval <- le_s("mapa_unido_bentos_avaliacao.csv")
  D$mu_repres <- le_s("mapa_unido_bentos_representacao.csv"); D$mu_hier <- le_s("mapa_unido_bentos_hierarquia.csv")
  rb <- rast(file.path(raiz, "saida/mapa_unido_bentos_25km.tif"))[["base"]]; vb <- values(rb)[, 1]; vb[is.na(vb)] <- 0L
  rl <- rle(as.integer(vb)); D$mu_grade <- list(ncol = ncol(rb), nrow = nrow(rb), comp = rl$lengths, val = rl$values)   # base classes, run-length coded
  f3 <- file.path(raiz, "saida/mapa_unido_bentos_nivel3.png"); if (file.exists(f3)) D$mu_mapa3 <- png64(f3)
}

# ---- Pelagico no nivel 2 (93_pelagico_n2.R) ------------------------------------------------------
D$pn2_un <- le_s("pelagico_n2_unidades.csv")
pi2 <- le_s("pelagico_n2_gf_importancia.csv")
if (!is.null(pi2)) D$pn2_imp <- pi2 |> group_by(grupo) |> arrange(desc(importancia), .by_group = TRUE) |> mutate(parcela = importancia / sum(importancia)) |>
  summarise(r2_medio = first(r2_medio), especies_ajustadas = first(especies_ajustadas),
    preditores = paste(sprintf("%s (%.0f%%)", head(preditor, 6), 100 * head(parcela, 6)), collapse = ", "), .groups = "drop") |> as.data.frame()
D$pn2_vs <- le_s("pelagico_n2_limiares_vs_n1.csv") |> select(any_of(c("grupo", "preditor", "limiar", "limiar_n1", "coincide_n1"))) |> as.data.frame()
pf <- le_s("pelagico_n2_fusoes.csv")
if (!is.null(pf)) { nf <- tapply(pf$hp_n1, pf$final_n1n2, function(x) x[order(as.integer(sub("HP", "", x)))][1]); D$pn2_fus <- pf |> mutate(hp = as.character(nf[as.character(final_n1n2)])) }

# ---- Mapa unido do pelagico (91b_mapa_unido_pelagico.R) ------------------------------------------
fmp <- file.path(raiz, "saida/mapa_unido_pelagico_classes.csv")
if (file.exists(fmp)) {
  D$mp_classes <- read.csv(fmp) |> select(any_of(c("codigo", "zona", "nivel3_rotulo", "hp", "emu", "evidencia", "evidencia_nivel3", "regime", "persistencia",
    "celulas_max", "candidatos", "cor", "standby"))) |> mutate(persistencia = round(persistencia, 2)) |> as.data.frame()
  D$mp_niveis <- le_s("mapa_unido_pelagico_niveis.csv"); D$mp_aval <- le_s("mapa_unido_pelagico_avaliacao.csv")
  mm <- c("0m_jan", "0m_jul", "200m_jan", "200m_jul", "nivel3_0m_jan", "nivel3_500m_jan")
  D$mp_mapas <- setNames(lapply(mm, function(a) { f <- file.path(raiz, sprintf("saida/mapa_unido_pelagico_%s.png", a)); if (file.exists(f)) png64(f) }), mm)
}

# ---- Gate 10: limiares publicados (92, 92b, 92c, 92d) --------------------------------------------
le_u <- function(a) { f <- file.path(raiz, "saida", a); if (file.exists(f)) read.csv(f, encoding = "UTF-8") }
en <- le_u("limiares_enriquecimento.csv"); lpp <- le_u("limiares_publicados_padronizados.csv")
if (!is.null(en)) {
  por_id <- en |> group_by(id) |> summarise(compartimento = first(compartimento), variavel_modelo = first(variavel_modelo), bacia = first(bacia), uso = first(uso), peso = first(peso),
    intervalo = ifelse(first(limiar_min) == first(limiar_max), format(first(limiar_min)), paste(first(limiar_min), "–", first(limiar_max))), grupos = n(),
    razao_max = suppressWarnings(max(razao, na.rm = TRUE)), grupo_max = if (all(is.na(razao))) NA_character_ else grupo_projeto[which.max(razao)],
    confirmado = any(confirmado_enriq), .groups = "drop") |> mutate(razao_max = ifelse(is.finite(razao_max), razao_max, NA)) |> as.data.frame()
  num <- function(x) vapply(x, function(v) if (is.na(v)) "" else prettyNum(v, big.mark = ".", decimal.mark = ","), "")
  ftb <- le_u("figura_limiares_tabela.csv"); lid <- por_id |> select(id, uso, peso, razao_max, grupo_max, confirmado)
  D$lim_tab <- lpp |> group_by(id) |> slice(1) |> ungroup() |> transmute(id, referencia, url = sub("\\?.*$", "", url), grupo_publicado, regiao, tipo_limiar,
      intervalo = ifelse(is.na(limiar_min), "sem valor no resumo", paste0(ifelse(limiar_min == limiar_max, num(limiar_min), paste(num(limiar_min), "–", num(limiar_max))), " ", unidade))) |>
    left_join(lid, by = "id") |> left_join(ftb, by = "id") |> arrange(is.na(ordem_painel), ordem_painel, linha, id) |> mutate(uso = recode(uso, validacao = "validação")) |> as.data.frame()
  D$lim_uso <- por_id |> group_by(uso) |> summarise(limiares = n(), confirmados = sum(confirmado), ids = paste(id[confirmado], collapse = ", "), .groups = "drop") |>
    arrange(match(uso, c("validacao", "apoio", "relato"))) |> mutate(uso = recode(uso, validacao = "validação")) |> as.data.frame()
}
f <- file.path(raiz, "saida/figura_limiares.png"); if (file.exists(f)) D$lim_fig <- png64(f)
D$lim_diag <- le_s("limiares_diagnostico_z.csv")
red <- le_u("limiares_candidatos_redefinicao.csv"); if (!is.null(red)) D$lim_red <- red |> select(-caixa_wkt) |> as.data.frame()
af <- le_u("limiares_afericao.csv"); afg <- le_u("limiares_afericao_grupos.csv")
if (!is.null(af)) {
  rz <- if (!is.null(afg)) afg |> group_by(id) |> filter(any(!is.na(razao)), any(!is.na(razao_glorys))) |>
    summarise(razao = gsub(".", ",", sprintf("%.2f → %.2f", max(razao, na.rm = TRUE), max(razao_glorys, na.rm = TRUE)), fixed = TRUE), .groups = "drop") else NULL
  iv <- function(a, b) ifelse(is.na(a) | is.na(b), "–", ifelse(round(a) == round(b), format(round(a)), paste(round(a), "–", round(b))))
  D$lim_af <- af |> filter(!is.na(eventos)) |> mutate(intervalo = iv(limiar_min, limiar_max), intervalo_glorys = iv(limiar_min_glorys, limiar_max_glorys), desloc_rel = round(desloc_rel, 2)) |>
    left_join(rz, by = "id") |> select(id, referencia, intervalo, eventos, desloc_med, desloc_rel, batim_obis_desloc_med, intervalo_glorys, any_of("razao")) |>
    mutate(across(c(desloc_med, batim_obis_desloc_med), round)) |> as.data.frame()
  if (!is.null(afg)) D$lim_afg <- as.data.frame.matrix(table(afg$situacao, afg$situacao_glorys)) |> tibble::rownames_to_column("situacao")
}
D$lim_var <- le_u("limiares_variaveis_candidatas.csv")

# ---- Comparacao dos mapas (94_comparacao_mapas.R) -------------------------------------------------
fc <- list.files(file.path(raiz, "saida/comparacao"), pattern = "\\.png$", full.names = TRUE)
if (length(fc)) D$cmp <- setNames(lapply(fc, png64), sub("\\.png$", "", basename(fc)))


# ---- Secoes novas: opcao B, Gate 12, limites verticais, etapa 4, Box de lacunas (13/11/2026) -----
# Each section = list(id, titulo, sub, blocos); each block = table (rows with labelled columns) or image; rendered by a generic JS block.
lu <- function(a) { f <- file.path(raiz, "saida", a); if (file.exists(f)) read.csv(f, encoding = "UTF-8", check.names = FALSE) }
r3 <- function(d) d |> mutate(across(where(is.double), ~ round(.x, 3))) |> as.data.frame(check.names = FALSE)
tb <- function(titulo, rows, nota = NULL) if (!is.null(rows) && nrow(rows)) list(tipo = "tab", titulo = titulo, nota = nota, rows = r3(rows))
im <- function(titulo, arq, alt) { f <- file.path(raiz, "saida", arq); if (file.exists(f)) list(tipo = "img", titulo = titulo, alt = alt, img = png64(f)) }
med <- function(x) round(median(x, na.rm = TRUE), 3)
S <- list()

# opcao B (Gate 11) and test of the approved HP
b_cmp <- lu("opcaoB_comparacao.csv"); b_par <- lu("pelagico_b_pares.csv"); b_fus <- lu("pelagico_b_fusoes.csv"); b_tot <- lu("nivel1b_totais.csv")
l3p <- lu("limiares_pelagico_redefinicao.csv"); l3t <- lu("limiares_pelagico_testes.csv")
S$b <- list(id = "fB", titulo = "Opção B (Gate 11) · grupos ecológicos com filtro de aparelho e teste dos HP aprovados",
  sub = "Micronécton e peixes pelágicos pequenos extraídos do OBIS com o filtro de aparelho do Gate 11. A versão B dos habitats pelágicos fica como sensibilidade; a fauna B testa os HP aprovados (93b) e entra no mapa unido do pelágico (91b).",
  blocos = Filter(Negate(is.null), list(
    tb("Extração da opção B", b_tot |> transmute(Grupo = grupo, Registros = registros, `Ausências` = ausencias, `Desde 1993` = desde_1993, Conjuntos = conjuntos)),
    tb("HP aprovados × versão B, por grupo", b_cmp |> transmute(Grupo = grupo, `Grupo do nível 1` = ifelse(nivel_1_original, "sim", "não"), `Versão` = versao, Classes = classes,
      Unidades = unidades, `ω²` = omega2, `ω² esforço igual` = omega2_eq, `Poder` = poder_eq)),
    tb("Pares de HP aprovados testados pela fauna B", b_par |> transmute(Grupo = grupo, `HP a` = classe_a, `HP b` = classe_b, `Unidades a` = unidades_a, `Unidades b` = unidades_b,
      Blocos = blocos, `ω²` = omega2, `ω² esforço igual` = omega2_eq, Poder = poder_eq, Separa = ifelse(separa, "sim", "não"))),
    tb("HP depois da fusão com o consenso n1 + n2 + B", b_fus |> transmute(`HP do nível 1` = hp_n1, `HP final` = final, `Stand-by` = ifelse(standby, "sim", "não"), `Testado pela fauna B` = ifelse(testado_b, "sim", "não")),
      sprintf("%d HP finais; %d em stand-by.", length(unique(b_fus$final)), length(unique(b_fus$final[b_fus$standby])))),
    tb("L3 do pelágico: limiares publicados testados como fronteira", l3p |> transmute(Id = id, Estudo = referencia, `Região` = regiao, `Valor (m)` = valor, `Limite de zona` = ifelse(limite_de_zona, "sim", "não"),
      Grupos = grupos, Distintos = distintos, `ω² mediano` = omega2_med, Consenso = consenso, `Decisão` = decisao)),
    tb("L3 do pelágico: testes por grupo", l3t |> transmute(Id = id, `Valor (m)` = valor, Grupo = grupo, `Unidades rasas` = unidades_raso, `Unidades fundas` = unidades_fundo,
      `ω²` = omega2, `ω² esforço igual` = omega2_eq, Poder = poder_eq, Separa = ifelse(separa, "sim", "não"))))))

# Gate 12: test fauna
c_un <- lu("etapa3_n1c_unidades.csv"); c_av <- lu("teste_fauna_c_avaliacao.csv"); c_pa <- lu("teste_fauna_c_pares.csv"); c_cl <- lu("teste_fauna_c_classes.csv")
ref <- lu("mapa_unido_bentos_avaliacao.csv")
if (!is.null(c_av)) {
  rf <- if (!is.null(ref)) ref |> group_by(classificacao) |> summarise(ref = med(omega2_eq)) else data.frame(classificacao = character(), ref = numeric())
  S$c <- list(id = "fG12", titulo = "Gate 12 · Fauna de teste independente dos mapas aprovados",
    sub = "Grupos que não entraram em nenhum mapa: peixes demersais e Decapoda bentônicos (bentos); Cyclopoida pelágicos, Appendicularia e Pteropoda (pelágico, com o filtro de aparelho). Nada aprovado é refeito; a evidência nova vai para o material suplementar (S20) e para a versão seguinte do mapa.",
    blocos = Filter(Negate(is.null), list(
      tb("Unidades da fauna de teste (87, sufixo c)", c_un |> transmute(Compartimento = comp, Grupo = grupo, Unidades = unidades, `Sem ambiente` = sem_ambiente, `Táxons` = taxons)),
      tb("Aceitação dos mapas pela fauna de teste", c_av |> left_join(rf, by = "classificacao") |> transmute(Compartimento = comp, Grupo = grupo, `Classificação` = classificacao, Classes = classes,
        Unidades = unidades, `ω²` = omega2, `ω² esforço igual` = omega2_eq, Poder = poder_eq, `Referência (mediana dos grupos de construção)` = if (is.null(ref)) NA else ifelse(comp == "bentos", ref, NA)),
        "Esforço igual = 10 unidades por classe. A referência é a mediana do ω² com esforço igual dos grupos que construíram o mapa do bentos."),
      tb("Pares de classes testados pela fauna de teste", c_pa |> group_by(Compartimento = comp, Grupo = grupo, `Classificação` = classificacao) |> summarise(Pares = n(), Separam = sum(separa), .groups = "drop")),
      tb("O que a fauna de teste acrescenta às classes do nível base", c_cl |> group_by(Compartimento = comp, `Evidência antes` = evidencia_antes) |>
        summarise(Classes = n(), `Testadas agora` = sum(testada_c), `Distintas agora` = sum(distinta_c), .groups = "drop")))))
}

# vertical limits
v_un <- lu("limites_verticais_unidades.csv"); v_pa <- lu("limites_verticais_pares.csv"); v_pe <- lu("limites_verticais_perfil.csv")
if (!is.null(v_pa)) {
  pw <- v_pe |> mutate(omega2 = round(omega2, 3)) |> select(grupo, corte, omega2) |> tidyr::pivot_wider(names_from = corte, values_from = omega2) |> rename(Grupo = grupo)
  names(pw)[-1] <- paste(names(pw)[-1], "m")
  S$v <- list(id = "fLV", titulo = "Gate 12 · Limites verticais de 200 e 1.000 m",
    sub = "DNA ambiental de procariotos com a amostra como unidade (≥ 10 táxons, ao menos 50 m acima do fundo GLORYS, até 300 amostras por conjunto × zona); Syndiniales ficou de fora porque nenhuma amostra chega a 10 táxons. Grupos curados em célula de 100 km × zona, no ano todo (sensibilidade da regra da unidade). Separa = ω² ≥ 0,01 e poder ≥ 0,8.",
    blocos = Filter(Negate(is.null), list(
      tb("Unidades por zona", v_un |> transmute(Grupo = grupo, Fonte = fonte, Epi = epi, Meso = meso, Bati = bati)),
      tb("Pares de zonas adjacentes", v_pa |> transmute(Grupo = grupo, Fonte = fonte, `Zona a` = classe_a, `Zona b` = classe_b, `Unidades a` = unidades_a, `Unidades b` = unidades_b,
        Blocos = blocos, `ω²` = omega2, `ω² esforço igual` = omega2_eq, Poder = poder_eq, Separa = ifelse(is.na(separa), "sem teste", ifelse(separa, "sim", "não")))),
      tb("Perfil de ω² por corte de profundidade", pw, "Para cada corte, ω² de acima × abaixo entre as unidades a menos de um fator 2 do corte. Vazio = menos de 5 unidades de um lado."))))
}

# stage 4
e_vb <- lu("etapa4_validacao_blocos.csv"); e_va <- lu("etapa4_validacao_anos.csv"); e_ao <- lu("etapa4_aoa_resumo.csv")
if (!is.null(e_vb)) {
  rs <- function(d) { if (!"classes" %in% names(d)) d$classes <- NA_real_; d |> group_by(Compartimento = comp, `Classificação` = classificacao) |> summarise(Grupos = n_distinct(grupo), `Classes (mediana)` = median(classes, na.rm = TRUE),
    `Acerto (mediana)` = med(acerto), `Acaso (mediana)` = med(acaso), `Kappa (mediana)` = med(kappa), `Kappa mín.` = round(min(kappa, na.rm = TRUE), 3), `Kappa máx.` = round(max(kappa, na.rm = TRUE), 3), .groups = "drop") }
  S$e <- list(id = "fE4", titulo = "Etapa 4 · Mapas mensais, validação e área de aplicabilidade (material suplementar S19)",
    sub = "Validação: cada unidade recebe a classe do mapa e é classificada pela fauna do vizinho mais parecido no treino, deixando de fora o bloco de 2.500 km da unidade, ou o período (corte na mediana do ano). O kappa cai com o número de classes: compare classificações com número parecido de classes. Área de aplicabilidade (AOA, Meyer e Pebesma): distância no espaço do Gradient Forest ao treino mais próximo, com limiar do terceiro quartil + 1,5 IQR da validação por blocos.",
    blocos = Filter(Negate(is.null), list(
      im("Mapas mensais do pelágico", "etapa4_mensal_pelagico_0m.png", "Mapa unido do pelágico, 12 meses, 0 m"),
      im(NULL, "etapa4_mensal_pelagico_200m.png", "Mapa unido do pelágico, 12 meses, 200 m"),
      tb("Validação deixando blocos de fora (resumo)", rs(e_vb)),
      tb("Validação deixando períodos de fora (resumo)", rs(e_va)),
      list(tipo = "det", titulo = "Validação por grupo (blocos de fora)", rows = r3(e_vb |> select(Compartimento = comp, Grupo = grupo, Fonte = fonte, `Classificação` = classificacao, Classes = classes,
        `Classificáveis` = classificaveis, `Não classificáveis` = nao_classificaveis, Acerto = acerto, `Acerto equilibrado` = acerto_equilibrado, Acaso = acaso, Kappa = kappa))),
      tb("Área de aplicabilidade por faixa, zona e grau de evidência", e_ao |> transmute(Compartimento = comp, `Faixa, zona ou evidência` = faixa_zona, Classes = classes,
        `Área dentro da AOA` = area_dentro, `Classes com a maior parte fora` = classes_maioria_fora)),
      im("Área de aplicabilidade", "etapa4_aoa_bentos.png", "Área de aplicabilidade do bentos (25 km)"),
      im(NULL, "etapa4_aoa_pelagico_500m_jan.png", "Área de aplicabilidade do pelágico, janeiro, 500 m"))))
}

# gaps box
x_fx <- lu("box_lacunas_faixas.csv"); x_tb <- lu("box_lacunas_tabela.csv")
if (!is.null(x_fx)) S$x <- list(id = "fBox", titulo = "Figura 8 e Box 1 · Lacunas e como contribuir",
  sub = "Grau de evidência de cada classe dos mapas unidos, área em stand-by (achurado: lacuna de conhecimento, não ausência de habitat) e área de aplicabilidade. As classes do bentos que a fauna de teste do Gate 12 mostrou distintas têm cor própria, sem mudar o mapa aprovado. Área: células de 25 km no bentos; célula × mês × profundidade no pelágico.",
  blocos = Filter(Negate(is.null), list(
    im(NULL, "box_lacunas_figura.png", "Figura 8: grau de evidência no bentos e no pelágico a 500 m, e área por grau de evidência com a área de aplicabilidade"),
    tb("Lacunas por faixa ou zona", x_fx |> transmute(Compartimento = comp, `Faixa ou zona` = faixa_zona, Classes = classes, `Em stand-by` = classes_standby,
      `Área em stand-by` = area_standby, `Área com evidência da fauna` = area_fauna, `Área separada pela literatura` = area_literatura, `Área dentro da AOA` = area_aoa)),
    tb("Lacunas, dado que resolve e como contribuir", x_tb |> transmute(Lacuna = lacuna, Tamanho = tamanho, `Dado que resolve` = dado_que_resolve, `Como contribuir` = como_contribuir)))))
D$novas <- unname(S)

# ---- HTML ---------------------------------------------------------------------------------------
modelo <- readLines(file.path(raiz, "scripts/90_relatorio_modelo.html"), encoding = "UTF-8")
html <- sub("__DADOS__", toJSON(D, auto_unbox = TRUE, dataframe = "rows", na = "null", digits = 4),
            paste(modelo, collapse = "\n"), fixed = TRUE)
writeLines(html, file.path(raiz, "saida/relatorio_validacao.html"), useBytes = TRUE)
cat("relatorio gravado:", file.path(raiz, "saida/relatorio_validacao.html"), "\n")
