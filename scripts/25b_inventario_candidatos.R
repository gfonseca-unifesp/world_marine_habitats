# 25b_inventario_candidatos.R  (PI request of 06/11/2026: potential of the groups outside level 1, before any curation)
# Candidates chosen from saida/inventario_grupos_familia.csv. For each, without quality or gear filters (upper bound; for the new pelagic groups
# of option B the gear filter removed 93-96% of the records): presence records from 1993 on with month and position, distinct events per
# 0.25 degree cell x quarter x depth zone (< 200, 200-1,000, >= 1,000 m, no depth), and records per dataset.
# A potential unit = cell x quarter x zone with >= 10 events (the level-1 rule, on a rough 0.25 degree grid).
# Outputs: saida/inventario_candidatos.csv, saida/inventario_candidatos_conjuntos.csv
suppressPackageStartupMessages({library(arrow); library(dplyr)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); set_cpu_count(8)
tc <- read.csv(file.path(raiz, "entrada/taxon_compartimento.csv"), encoding = "UTF-8")
dem <- unique(tc$family[tc$compartimento == "bentos" & tc$class %in% c("Teleostei", "Actinopteri", "Elasmobranchii", "Holocephali")])
if (!length(dem)) dem <- unique(tc$family[tc$compartimento == "bentos" & tc$peixe %in% TRUE])
DEC_P <- c("Sergestidae", "Sergiidae", "Oplophoridae", "Acanthephyridae", "Pasiphaeidae", "Benthesicymidae", "Luciferidae", "Gnathophausiidae", "Lophogastridae")
MICRO <- c("Proteobacteria", "Bacteroidetes", "Cyanobacteria", "Euryarchaeota", "Thaumarchaeota", "Verrucomicrobia", "Actinobacteria", "Planctomycetes", "Chloroflexi")
cat("familias de peixes demersais:", length(dem), "\n")
# files are opened in chunks below
t0 <- Sys.time()
consulta <- function(ds) ds |> filter(!dropped, is.na(absence) | !absence) |>
  mutate(phylum = interpreted$phylum, subphylum = interpreted$subphylum, class = interpreted$class, order = interpreted$order, suborder = interpreted$suborder, family = interpreted$family,
         year = interpreted$date_year, ds_ms = interpreted$date_start, lon = interpreted$decimalLongitude, lat = interpreted$decimalLatitude, ev = interpreted$eventID,
         z = coalesce(interpreted$maximumDepthInMeters, interpreted$depth, interpreted$minimumDepthInMeters)) |>
  filter(year >= 1993L, !is.na(lon), !is.na(lat), !is.na(ds_ms)) |>                 # month from date_start: the month text field is empty in many datasets
  mutate(mes = month(cast(ds_ms, timestamp("ms")))) |>
  mutate(grupo = case_when(
    class == "Appendicularia" ~ "P Appendicularia", order == "Halocyprida" ~ "P Ostracoda pelágicos (Halocyprida)", order == "Cyclopoida" ~ "P Copepoda Cyclopoida",
    order %in% c("Pteropoda", "Thecosomata", "Gymnosomata") ~ "P Pteropoda", order %in% c("Trachymedusae", "Narcomedusae") | class == "Scyphozoa" | phylum == "Ctenophora" ~ "P Medusas e ctenóforos",
    family %in% DEC_P ~ "P Decápodes e lofogástridos pelágicos", phylum == "Radiozoa" ~ "P Radiozoa", phylum %in% MICRO ~ "P Procariotos", order == "Syndiniales" ~ "P Syndiniales",
    class %in% c("Bacillariophyceae", "Dinophyceae", "Coccolithophyceae") ~ "P Fitoplâncton", phylum == "Foraminifera" ~ "M Foraminifera",
    class == "Gastropoda" ~ "B Gastropoda", order == "Decapoda" ~ "B Decapoda", order == "Amphipoda" ~ "B Amphipoda (sem Hyperiidea)", class == "Thecostraca" ~ "B Cirripedia",
    family %in% dem ~ "B Peixes demersais", phylum == "Nematoda" | order == "Harpacticoida" ~ "B Meiofauna (Nematoda e Harpacticoida)", order %in% c("Leptothecata", "Anthoathecata") ~ "B Hydrozoa bentônicos",
    phylum %in% c("Brachiopoda", "Nemertea") | class == "Sipuncula" | class == "Pycnogonida" ~ "B Brachiopoda, Nemertea, Sipuncula, Pycnogonida", class %in% c("Phaeophyceae", "Florideophyceae", "Ulvophyceae") ~ "B Macroalgas")) |>
  filter(!is.na(grupo), !(grupo == "B Amphipoda (sem Hyperiidea)" & suborder == "Hyperiidea" & !is.na(suborder))) |>
  mutate(zona = case_when(is.na(z) ~ "sem_z", z < 200 ~ "z0_200", z < 1000 ~ "z200_1000", TRUE ~ "z1000"), celx = floor(lon * 4), cely = floor(lat * 4), tri = floor((mes - 1) / 3) + 1,
         chave = coalesce(ev, paste(dataset_id, interpreted$eventDate, lon, lat))) |>
  group_by(grupo, dataset_id, zona, celx, cely, tri, chave) |> summarise(registros = n(), .groups = "drop")
# a single scan of all files ran out of memory (06/11 and 13/11): the files are read in chunks, each chunk reduced to distinct
# events x cell x quarter x zone and written to disk; events are then counted per dataset and summed (an event belongs to one dataset).
tmp <- file.path(raiz, "intermediario/inventario_candidatos_eventos"); if (dir.exists(tmp)) unlink(tmp, recursive = TRUE); dir.create(tmp)
arqs <- list.files(file.path(raiz, "bruto/occurrence"), full.names = TRUE); lotes <- split(arqs, ceiling(seq_along(arqs) / 150))
for (k in seq_along(lotes)) { x <- consulta(open_dataset(lotes[[k]], format = "parquet", unify_schemas = FALSE)) |> collect()
  if (nrow(x)) write_parquet(x, file.path(tmp, sprintf("lote_%03d.parquet", k))); rm(x); gc(verbose = FALSE)
  if (k %% 10 == 0) cat(format(Sys.time()), "lote", k, "de", length(lotes), "\n") }
cat(format(Sys.time()), "eventos gravados em", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
dd <- open_dataset(tmp); gs <- dd |> distinct(grupo) |> collect() |> pull(grupo)
a0 <- bind_rows(lapply(gs, function(g) dd |> filter(grupo == g) |> group_by(grupo, dataset_id, zona, celx, cely, tri) |>
  summarise(eventos = n_distinct(chave), registros = sum(registros), .groups = "drop") |> collect()))
cat(format(Sys.time()), "agregado:", nrow(a0), "linhas em", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
write_parquet(a0, file.path(raiz, "intermediario/inventario_candidatos_celulas.parquet"))
ag <- a0 |> group_by(grupo, zona, celx, cely, tri) |> summarise(eventos = sum(eventos), registros = sum(registros), .groups = "drop")
cj <- a0 |> group_by(grupo, dataset_id, zona) |> summarise(registros = sum(registros), .groups = "drop")
res <- ag |> group_by(grupo, zona) |> summarise(registros = sum(registros), unidades_10ev = sum(eventos >= 10), .groups = "drop") |>
  tidyr::pivot_wider(names_from = zona, values_from = c(registros, unidades_10ev), id_cols = grupo, values_fill = 0) |> arrange(grupo)
write.csv(res, file.path(raiz, "saida/inventario_candidatos.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(cj, file.path(raiz, "saida/inventario_candidatos_conjuntos.csv"), row.names = FALSE, fileEncoding = "UTF-8")
options(width = 200); print(as.data.frame(res)); cat(format(Sys.time()), "FIM\n")
