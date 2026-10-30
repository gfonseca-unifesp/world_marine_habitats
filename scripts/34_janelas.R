# 34_janelas.R  (task 3.4)
# Level 2 candidate windows: 500 x 500 km (20 x 20 cells of 25 km), aligned to the 100 km grid.
# Mandatory filter: >= 1 cell with a time series and >= 10 boundary pairs (25 km).
# Scoring (percentiles among the windows that pass the filter):
#   3 x dynamic potential (GLORYS12 seasonal amplitude: surface, mixed layer and bottom)
#   2 x class diversity (EBU groups and regions + EMU with data)
#   2 x independence (number of datasets; penalizes dominance of the largest dataset)
#   1 x time series, 1 x boundary pairs, +1 benthos and pelagic together, +1 analogous region
# Complementary selection: Foz do Amazonas fixed; then, at each step, the window with the highest
#   score x (1 + 0.5 x novelties), where novelties = EBU geomorphological groups, EMU classes and
#   latitude bands not yet covered. Windows adjacent to one already chosen are excluded. Up to 8 regions.
# Output: saida/janelas_pontuadas.csv, saida/janelas_selecionadas.csv
suppressPackageStartupMessages({library(arrow); library(dplyr); library(terra)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
source(file.path(raiz, "scripts/funcoes_grade.R"))
g25 <- ease2_grade("25km")
W <- 20                                                     # 25 km cells per side of the window

cel <- read_parquet(file.path(raiz, "intermediario/sondagem_celulas.parquet")) |> filter(res == "25km")
cel <- cel |> mutate(wr = (row - 1) %/% W + 1, wc = (col - 1) %/% W + 1, janela = paste(wr, wc, sep = "_"))
base <- cel |> filter(nivel %in% c("grupo", "emu"))         # one level per compartment, to avoid double counting

# time series and compartments per window
j_serie <- base |> group_by(janela) |> summarise(celulas_serie = n_distinct(cell_id[serie]),
                                                  celulas_dados = n_distinct(cell_id),
                                                  bentos = any(compartimento == "bentos"), pelagico = any(compartimento == "pelagico"))
# boundary pairs within the window (EBU region in benthos; EMU in pelagic), by group
pares <- cel |> filter(nivel %in% c("regiao", "emu")) |> group_by(grupo) |> group_modify(function(d, k) {
  key <- paste(d$row, d$col); cl <- setNames(d$classe, key)
  v <- rbind(data.frame(a = key, b = paste(d$row, d$col + 1)), data.frame(a = key, b = paste(d$row + 1, d$col)))
  v <- v[v$b %in% key & cl[v$a] != cl[v$b], ]
  rc <- do.call(rbind, strsplit(v$a, " ")); data.frame(janela = paste((as.integer(rc[, 1]) - 1) %/% W + 1, (as.integer(rc[, 2]) - 1) %/% W + 1, sep = "_"))
}) |> ungroup() |> count(janela, name = "pares_fronteira")
# class diversity
div <- cel |> filter(nivel %in% c("grupo", "regiao", "emu")) |> group_by(janela) |>
  summarise(ebu_grupos = n_distinct(classe[nivel == "grupo"]), ebu_regioes = n_distinct(classe[nivel == "regiao"]),
            emus = n_distinct(classe[nivel == "emu"]),
            lista_ebu_grupos = paste(sort(unique(classe[nivel == "grupo"])), collapse = ""),
            lista_emus = paste(sort(unique(classe[nivel == "emu"])), collapse = ","))
# independence: distinct datasets per window (presences with valid coordinate)
ds <- open_dataset(file.path(raiz, "intermediario/nivel1_extracao"))
p <- ds |> filter(!(absence %in% TRUE), !is.na(lon), !is.na(lat), !is.na(year)) |> select(dataset_id, lon, lat) |> collect()
p <- p[abs(p$lat) <= 90 & abs(p$lon) <= 180, ]
c25 <- cellFromXY(g25, project(cbind(p$lon, p$lat), from = "EPSG:4326", to = "EPSG:6933"))
rc <- rowColFromCell(g25, c25); p$janela <- paste((rc[, 1] - 1) %/% W + 1, (rc[, 2] - 1) %/% W + 1, sep = "_")
p <- p[!is.na(c25), ]
ind <- p |> count(janela, dataset_id) |> group_by(janela) |>
  summarise(conjuntos = n(), pct_maior_conjunto = round(100 * max(n) / sum(n), 1), registros = sum(n))
rm(p); gc(verbose = FALSE)
# dynamic potential (GLORYS12, 100 km cells -> window = 5 x 5 cells of 100 km)
amp <- read_parquet(file.path(raiz, "camadas/glorys/amplitude_100km.parquet"))
a_rc <- do.call(rbind, strsplit(amp$cell_id, "_"))
amp$janela <- paste((as.integer(a_rc[, 2]) - 1) %/% 5 + 1, (as.integer(a_rc[, 3]) - 1) %/% 5 + 1, sep = "_")
pr <- function(v) rank(v, na.last = "keep") / sum(!is.na(v))
amp <- amp |> mutate(d_sup = (pr(sst_amp) + pr(sss_amp) + pr(mld_amp)) / 3, d_fundo = pr(tfundo_amp))
din <- amp |> group_by(janela) |> summarise(dinamico_sup = mean(d_sup, na.rm = TRUE), dinamico_fundo = mean(d_fundo, na.rm = TRUE),
                                            sst_amp = mean(sst_amp, na.rm = TRUE), tfundo_amp = mean(tfundo_amp, na.rm = TRUE))
# window centre, latitude band and analogous regions
jan <- j_serie |> left_join(pares, by = "janela") |> left_join(div, by = "janela") |> left_join(ind, by = "janela") |>
  left_join(din, by = "janela") |> mutate(pares_fronteira = coalesce(pares_fronteira, 0L))
w_rc <- do.call(rbind, strsplit(jan$janela, "_")); wr <- as.integer(w_rc[, 1]); wc <- as.integer(w_rc[, 2])
cx <- xFromCol(g25, (wc - 1) * W + W / 2); cy <- yFromRow(g25, (wr - 1) * W + W / 2)
ll <- project(cbind(cx, cy), from = "EPSG:6933", to = "EPSG:4326"); jan$lon <- round(ll[, 1], 2); jan$lat <- round(ll[, 2], 2)
jan$faixa_lat <- cut(abs(jan$lat), c(0, 23.5, 60, 90), c("tropical", "temperada", "polar"), include.lowest = TRUE)
analogas <- data.frame(nome = c("Foz do Amazonas", "Pluma do Congo", "Pluma do Orinoco", "Pluma do Mississippi",
  "Golfo de Bengala (Ganges-Brahmaputra)", "Rio da Prata", "Pluma do Yangtzé", "Ressurgência de Benguela",
  "Ressurgência de Humboldt", "Ressurgência da Califórnia", "Ressurgência das Canárias"),
  tipo = c(rep("pluma", 7), rep("ressurgência", 4)),
  lon = c(-49, 12.3, -60.5, -89.2, 90.5, -56.5, 122.3, 14, -76.5, -122.8, -16.5),
  lat = c(1.5, -6.1, 8.8, 29, 21, -35.5, 31.3, -23.5, -14, 37.5, 25.5))
a_cel <- cellFromXY(g25, project(as.matrix(analogas[, c("lon", "lat")]), from = "EPSG:4326", to = "EPSG:6933"))
a_rc2 <- rowColFromCell(g25, a_cel)
analogas$janela <- paste((a_rc2[, 1] - 1) %/% W + 1, (a_rc2[, 2] - 1) %/% W + 1, sep = "_")
jan$analoga <- analogas$nome[match(jan$janela, analogas$janela)]

# filter and scoring
cand <- jan |> mutate(passa_filtro = celulas_serie >= 1 & pares_fronteira >= 10)
f <- cand |> filter(passa_filtro)
pr2 <- function(v) rank(v, na.last = "keep", ties.method = "average") / sum(!is.na(v))
f <- f |> mutate(din = pmax(coalesce(dinamico_sup, 0), coalesce(dinamico_fundo, 0)),
                 s_din = pr2(din), s_div = pr2(ebu_grupos + ebu_regioes + emus),
                 s_ind = pr2(conjuntos * (1 - pct_maior_conjunto / 100)),
                 s_ser = pr2(celulas_serie), s_par = pr2(pares_fronteira),
                 pontuacao = round(3 * s_din + 2 * s_div + 2 * s_ind + s_ser + s_par + (bentos & pelagico) + !is.na(analoga), 3))
cand <- cand |> left_join(f |> select(janela, pontuacao, s_din, s_div, s_ind, s_ser, s_par), by = "janela")
write.csv(cand |> arrange(desc(pontuacao)), file.path(raiz, "saida/janelas_pontuadas.csv"), row.names = FALSE)

# complementary selection
ebu_set <- function(s) unique(strsplit(paste(s, collapse = ""), "")[[1]])
emu_set <- function(s) unique(unlist(strsplit(paste(s, collapse = ","), ",")))
sel <- cand |> filter(janela == analogas$janela[1]) |> mutate(ordem = 1, motivo = "inclusão fixa (Foz do Amazonas)")
resto <- f |> filter(!janela %in% sel$janela)
while (nrow(sel) < 8 && nrow(resto) > 0) {
  cob_ebu <- ebu_set(sel$lista_ebu_grupos); cob_emu <- emu_set(sel$lista_emus); cob_lat <- unique(as.character(sel$faixa_lat))
  nov <- sapply(seq_len(nrow(resto)), function(i) length(setdiff(ebu_set(resto$lista_ebu_grupos[i]), cob_ebu)) +
                  length(setdiff(emu_set(resto$lista_emus[i]), cob_emu)) + !(as.character(resto$faixa_lat[i]) %in% cob_lat))
  adj <- resto$pontuacao * (1 + 0.5 * nov)
  # do not choose a window adjacent (including diagonally) to one already chosen
  s_rc <- do.call(rbind, strsplit(sel$janela, "_")); r_rc <- do.call(rbind, strsplit(resto$janela, "_"))
  encosta <- sapply(seq_len(nrow(resto)), function(i) any(abs(as.integer(r_rc[i, 1]) - as.integer(s_rc[, 1])) <= 1 &
                                                          abs(as.integer(r_rc[i, 2]) - as.integer(s_rc[, 2])) <= 1))
  adj[encosta] <- -Inf
  if (all(!is.finite(adj))) break
  i <- which.max(adj)
  sel <- bind_rows(sel, resto[i, ] |> mutate(ordem = nrow(sel) + 1,
                   motivo = sprintf("pontuação %.2f; %d novidade(s)", resto$pontuacao[i], nov[i])))
  resto <- resto[-i, ]
}
write.csv(sel, file.path(raiz, "saida/janelas_selecionadas.csv"), row.names = FALSE)
cat("janelas com dados:", nrow(cand), "| passam no filtro:", sum(cand$passa_filtro), "\n")
print(sel |> select(ordem, janela, lon, lat, faixa_lat, analoga, celulas_serie, pares_fronteira, conjuntos, pontuacao, motivo) |> as.data.frame())
cat("\nAnálogas e o filtro:\n")
print(cand |> filter(!is.na(analoga)) |> select(analoga, passa_filtro, celulas_serie, pares_fronteira, pontuacao) |> as.data.frame())
