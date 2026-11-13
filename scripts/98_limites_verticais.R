# 98_limites_verticais.R  (Gate 12, 13/11/2026, item 3: are the 200 m and 1,000 m limits of the pelagic map seen by the fauna below 200 m?)
# A Environmental DNA (OBIS inventory 25/25b): prokaryotes (Bacteria and Archaea phyla of 25b) and Syndiniales. Each sample (dataset x event)
#   is a unit, because one sample holds hundreds of taxa; taxon = genus, else family. Records from 1993 on with a sampling depth; water
#   samples only: depth at least 50 m above the GLORYS bottom of the 25 km cell (drops sediment samples). Samples with >= 10 taxa; at most
#   300 samples per dataset x zone (seed 1), so that one cruise does not dominate.
# B Animal groups already curated (level 1, option B and the test fauna c): units = 100 km cell x vertical zone, all year (sensitivity of
#   the unit rule; the 25 km cell x quarter rule leaves almost nothing below 200 m), >= 10 events, 10 events drawn (seed 1001).
# Tests, per group, on presence-absence beta_sim with 2,500 km blocks: pairs of adjacent zones (epi x meso, meso x bati; pares_classes:
#   omega2 >= 0.01 and power >= 0.8), and a depth profile: for cuts at 50-2,000 m, omega2 of above x below among units within a factor 2
#   of the cut (where turnover concentrates). Outputs: saida/limites_verticais_{pares,perfil,unidades}.csv
suppressPackageStartupMessages({library(arrow); library(dplyr); library(terra); library(tidyr); library(parallel)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); source(file.path(raiz, "scripts/funcoes_grade.R")); set_cpu_count(8)
g25 <- ease2_grade("25km"); pg <- values(rast(file.path(raiz, "camadas/glorys/mes/glorys_01_25km.tif"))[["prof_glorys"]])[, 1]
zona <- function(z) c("epi", "meso", "bati")[findInterval(z, c(200, 1000)) + 1]; CORTES <- c(50, 100, 150, 200, 300, 500, 750, 1000, 1500, 2000)
MICRO <- c("Proteobacteria", "Bacteroidetes", "Cyanobacteria", "Euryarchaeota", "Thaumarchaeota", "Verrucomicrobia", "Actinobacteria", "Planctomycetes", "Chloroflexi")

# ---------------- A environmental DNA ----------------
fa <- file.path(raiz, "intermediario/limites_verticais_edna.parquet")
if (!file.exists(fa)) {
  arqs <- list.files(file.path(raiz, "bruto/occurrence"), full.names = TRUE); lotes <- split(arqs, ceiling(seq_along(arqs) / 150)); partes <- list()
  for (k in seq_along(lotes)) { partes[[k]] <- open_dataset(lotes[[k]], format = "parquet", unify_schemas = FALSE) |> filter(!dropped, is.na(absence) | !absence) |>
      mutate(phylum = interpreted$phylum, order = interpreted$order, family = interpreted$family, genus = interpreted$genus, year = interpreted$date_year,
             lon = interpreted$decimalLongitude, lat = interpreted$decimalLatitude, ev = interpreted$eventID, ed = interpreted$eventDate,
             z = coalesce(interpreted$maximumDepthInMeters, interpreted$depth, interpreted$minimumDepthInMeters)) |>
      filter(phylum %in% MICRO | order == "Syndiniales", year >= 1993L, !is.na(z), !is.na(lon), !is.na(lat)) |>
      mutate(grupo = ifelse(order %in% "Syndiniales", "Syndiniales (DNA)", "Procariotos (DNA)"), taxon = coalesce(genus, family),
             amostra = coalesce(ev, paste(ed, lon, lat, z))) |> filter(!is.na(taxon)) |>
      distinct(grupo, dataset_id, amostra, lon, lat, z, taxon) |> collect()
    if (k %% 10 == 0) cat(format(Sys.time()), "lote", k, "de", length(lotes), "\n") }
  write_parquet(bind_rows(partes), fa); rm(partes); gc(verbose = FALSE) }
a <- read_parquet(fa)
cel <- cellFromXY(g25, project(cbind(a$lon, a$lat), from = "EPSG:4326", to = crs(g25))); rc <- rowColFromCell(g25, cel)
a <- a |> mutate(fundo = pg[cel], bloco = paste((rc[, 1] - 1) %/% 100 + 1, (rc[, 2] - 1) %/% 100 + 1, sep = "_")) |> filter(!is.na(fundo), z <= fundo - 50)
am <- a |> group_by(grupo, dataset_id, amostra) |> summarise(z = median(z), bloco = first(bloco), taxa = n_distinct(taxon), .groups = "drop") |> filter(taxa >= 10) |>
  mutate(zona = zona(z)) |> group_by(grupo, dataset_id, zona) |> group_modify(~ { set.seed(1); .x[sample(nrow(.x), min(nrow(.x), 300)), ] }) |> ungroup()
UA <- lapply(split(am, am$grupo), function(s) { d <- a |> semi_join(s, by = c("grupo", "dataset_id", "amostra")) |> distinct(dataset_id, amostra, taxon) |>
    mutate(u = paste(dataset_id, amostra), v = 1L) |> select(u, taxon, v) |> pivot_wider(names_from = taxon, values_from = v, values_fill = 0L)
  X <- as.matrix(d[, -1]); rownames(X) <- d$u; s$u <- paste(s$dataset_id, s$amostra); s <- s[match(rownames(X), s$u), ]
  list(grupo = s$grupo[1], fonte = "DNA ambiental, amostra", X = X[, colSums(X) >= 5, drop = FALSE], z = s$z, bloco = s$bloco) })
cat(format(Sys.time()), "DNA:", paste(names(UA), sapply(UA, function(u) nrow(u$X)), collapse = "; "), "amostras\n")

# ---------------- B curated animal groups, 100 km x zone, all year ----------------
tabs <- c("saida/tabela_certa_pelagico_n1.parquet", "saida/tabela_certa_pelagico_n1b.parquet", "saida/tabela_certa_pelagico_n1c.parquet")
p <- bind_rows(lapply(tabs[file.exists(file.path(raiz, tabs))], function(f) open_dataset(file.path(raiz, f)) |> filter(year >= 1993, nivel_certeza %in% c("A", "B")) |>
  select(grupo, event_id, cell_100, bloco_espacial, taxon_unidade, depth_min_m, depth_max_m) |> collect())) |>
  filter(!is.na(taxon_unidade), taxon_unidade != "", !grupo %in% c("Euphausiacea", "Hyperiidea")) |>
  mutate(z = coalesce(depth_max_m, depth_min_m)) |> filter(!is.na(z), is.na(depth_min_m) | zona(depth_min_m) == zona(z)) |>   # hauls within one zone
  mutate(zona = zona(z), unidade = paste(cell_100, zona))
UB <- lapply(split(p, p$grupo), function(d) { un <- d |> group_by(unidade) |> summarise(ev = n_distinct(event_id), z = median(z), bloco = first(bloco_espacial), .groups = "drop") |> filter(ev >= 10)
  if (sum(zona(un$z) != "epi") < 5) return(NULL)
  set.seed(1001); s <- d |> filter(unidade %in% un$unidade) |> distinct(unidade, event_id) |> group_by(unidade) |> slice_sample(n = 10) |> ungroup()
  m <- d |> semi_join(s, by = c("unidade", "event_id")) |> distinct(unidade, taxon_unidade) |> mutate(v = 1L) |> pivot_wider(names_from = taxon_unidade, values_from = v, values_fill = 0L)
  X <- as.matrix(m[, -1]); rownames(X) <- m$unidade; un <- un[match(rownames(X), un$unidade), ]
  list(grupo = d$grupo[1], fonte = "100 km x zona, ano todo", X = X[, colSums(X) > 0, drop = FALSE], z = un$z, bloco = un$bloco) })
U <- c(UA, Filter(Negate(is.null), UB))
uni <- bind_rows(lapply(U, function(u) data.frame(grupo = u$grupo, fonte = u$fonte, zona = zona(u$z)) |> count(grupo, fonte, zona) |> pivot_wider(names_from = zona, values_from = n, values_fill = 0)))
print(as.data.frame(uni))

# ---------------- tests ----------------
testa <- function(u) { suppressPackageStartupMessages({library(vegan); library(dplyr)}); source(file.path(raiz, "scripts/funcoes_permanova.R"))
  Dm <- as.matrix(sqrt(betadiver(u$X, "sim"))); Dm[is.na(Dm)] <- 1; zz <- c("epi", "meso", "bati")[findInterval(u$z, c(200, 1000)) + 1]; set.seed(7700)
  pr <- pares_classes(Dm, zz, u$bloco, un_min = 5, n_classe = 10, B = 100, pares = rbind(c("epi", "meso"), c("meso", "bati")))
  perf <- bind_rows(lapply(c(50, 100, 150, 200, 300, 500, 750, 1000, 1500, 2000), function(cc) { k <- u$z >= cc / 2 & u$z < cc * 2; lado <- ifelse(u$z[k] >= cc, "abaixo", "acima")
    if (min(table(factor(lado, c("acima", "abaixo")))) < 5) return(NULL); set.seed(cc)
    x <- tryCatch(permanova_rapida(as.dist(Dm[k, k]), lado, u$bloco[k], 199), error = function(e) NULL); if (is.null(x)) return(NULL)
    data.frame(corte = cc, acima = sum(lado == "acima"), abaixo = sum(lado == "abaixo"), omega2 = x$omega2) }))
  list(pares = if (!is.null(pr)) cbind(grupo = u$grupo, fonte = u$fonte, pr), perfil = if (nrow(perf)) cbind(grupo = u$grupo, fonte = u$fonte, perf)) }
cl <- makeCluster(min(10, length(U))); clusterExport(cl, "raiz"); R <- parLapplyLB(cl, U, testa); stopCluster(cl)
pares <- bind_rows(lapply(R, `[[`, "pares")); perfil <- bind_rows(lapply(R, `[[`, "perfil"))
write.csv(pares, file.path(raiz, "saida/limites_verticais_pares.csv"), row.names = FALSE); write.csv(perfil, file.path(raiz, "saida/limites_verticais_perfil.csv"), row.names = FALSE)
write.csv(uni, file.path(raiz, "saida/limites_verticais_unidades.csv"), row.names = FALSE)
options(width = 200); print(pares |> mutate(across(where(is.double), ~ round(.x, 3))) |> as.data.frame())
print(perfil |> mutate(omega2 = round(omega2, 3)) |> select(grupo, corte, omega2) |> pivot_wider(names_from = corte, values_from = omega2) |> as.data.frame())
cat(format(Sys.time()), "FIM\n")
