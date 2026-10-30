# 71_extrair_nivel2.R  (task 7.1)
# Extracts from the local OBIS copy, in the 8 level 2 regions (entrada/regioes_nivel2.csv; 500 x 500 km windows),
# the focal taxa and the extended list (E5.2):
#   Decapoda; Actinopteri and Elasmobranchii (includes ichthyoplankton); Gastropoda; Foraminifera; Amphipoda;
#   Harpacticoida; Ostracoda; Mysida; Cephalopoda; Hydrozoa.
# Excludes the highly migratory species of UNCLOS Annex I (E5.4), flagged with reason "unclos_anexo1"
# in the excluded-records file (nothing is deleted).
# Output: intermediario/nivel2_extracao/ (by region), intermediario/nivel2_unclos/, saida/nivel2_totais.csv
suppressPackageStartupMessages({library(arrow); library(dplyr); library(terra)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
source(file.path(raiz, "scripts/funcoes_grade.R"))
reg <- read.csv(file.path(raiz, "entrada/regioes_nivel2.csv"))
g25 <- ease2_grade("25km"); W <- 20

# limits of each window in degrees (EASE-Grid 2.0 is cylindrical: the corners give the exact rectangle)
p <- do.call(rbind, strsplit(reg$janela, "_")); wr <- as.integer(p[, 1]); wc <- as.integer(p[, 2])
x0 <- xFromCol(g25, (wc - 1) * W + 1) - res(g25)[1] / 2; x1 <- xFromCol(g25, wc * W) + res(g25)[1] / 2
y1 <- yFromRow(g25, (wr - 1) * W + 1) + res(g25)[2] / 2; y0 <- yFromRow(g25, wr * W) - res(g25)[2] / 2
c0 <- project(cbind(x0, y0), from = "EPSG:6933", to = "EPSG:4326"); c1 <- project(cbind(x1, y1), from = "EPSG:6933", to = "EPSG:4326")
reg$lon_min <- c0[, 1]; reg$lat_min <- c0[, 2]; reg$lon_max <- c1[, 1]; reg$lat_max <- c1[, 2]
write.csv(reg, file.path(raiz, "entrada/regioes_nivel2.csv"), row.names = FALSE)

# level 2 taxa (focal + extended), by the WoRMS taxonomic field
tx <- list(phylum = c("Echinodermata", "Porifera", "Chaetognatha", "Bryozoa", "Foraminifera"),
           subphylum = c("Anthozoa"),
           class = c("Polychaeta", "Bivalvia", "Ascidiacea", "Thaliacea", "Actinopteri", "Elasmobranchii", "Gastropoda",
                     "Ostracoda", "Cephalopoda", "Hydrozoa"),
           order = c("Cumacea", "Isopoda", "Tanaidacea", "Calanoida", "Euphausiacea", "Siphonophorae", "Decapoda",
                     "Amphipoda", "Harpacticoida", "Mysida"))
f_ph <- tx$phylum; f_sp <- tx$subphylum; f_cl <- tx$class; f_or <- tx$order
# UNCLOS, Annex I (highly migratory species): genera and families listed
unclos_gen <- c("Thunnus", "Katsuwonus", "Euthynnus", "Auxis", "Tetrapturus", "Makaira", "Kajikia", "Istiompax", "Istiophorus",
                "Xiphias", "Scomberesox", "Cololabis", "Coryphaena", "Cetorhinus", "Rhincodon")
unclos_fam <- c("Bramidae", "Istiophoridae", "Xiphiidae", "Alopiidae", "Carcharhinidae", "Sphyrnidae", "Lamnidae")
unclos_esp <- c("Hexanchus griseus")

ds <- open_dataset(list.files(file.path(raiz, "bruto/occurrence"), full.names = TRUE), format = "parquet", unify_schemas = FALSE)
t0 <- Sys.time()
# filter for the 8 boxes built as an expression (arrow does not accept user functions inside filter)
caixa <- rlang::parse_expr(paste(sprintf("(lon >= %.6f & lon <= %.6f & lat >= %.6f & lat <= %.6f)",
                                         reg$lon_min, reg$lon_max, reg$lat_min, reg$lat_max), collapse = " | "))
x <- ds |>
  mutate(lon = interpreted$decimalLongitude, lat = interpreted$decimalLatitude,
         phylum = interpreted$phylum, subphylum = interpreted$subphylum, superclass = interpreted$superclass,
         class = interpreted$class, order = interpreted$order, suborder = interpreted$suborder) |>
  filter(!is.na(lon), !is.na(lat), !!caixa,
         phylum %in% f_ph | subphylum %in% f_sp | class %in% f_cl | order %in% f_or | superclass %in% "Actinopteri") |>
  mutate(record_id = `_id`, eventDate = interpreted$eventDate, date_start = interpreted$date_start, date_end = interpreted$date_end,
         year = interpreted$date_year, year_txt = interpreted$year, month_txt = interpreted$month, day_txt = interpreted$day,
         verbatimEventDate = interpreted$verbatimEventDate, depth_min = interpreted$minimumDepthInMeters,
         depth_max = interpreted$maximumDepthInMeters, depth = interpreted$depth, bathymetry = interpreted$bathymetry,
         shoredistance = interpreted$shoredistance, eventID = interpreted$eventID, parentEventID = interpreted$parentEventID,
         occurrenceID = interpreted$occurrenceID, aphiaid = interpreted$aphiaid, scientificName = interpreted$scientificName,
         taxonRank = interpreted$taxonRank, family = interpreted$family, genus = interpreted$genus, species = interpreted$species,
         samplingProtocol = interpreted$samplingProtocol, lifeStage = interpreted$lifeStage,
         coordUncertainty = interpreted$coordinateUncertaintyInMeters, basisOfRecord = interpreted$basisOfRecord,
         occurrenceStatus = interpreted$occurrenceStatus, identificationQualifier = interpreted$identificationQualifier,
         verbatimLon = interpreted$verbatimLongitude, verbatimLat = interpreted$verbatimLatitude) |>
  select(record_id, dataset_id, lon, lat, verbatimLon, verbatimLat, eventDate, date_start, date_end, year, year_txt, month_txt,
         day_txt, verbatimEventDate, depth_min, depth_max, depth, bathymetry, shoredistance, eventID, parentEventID, occurrenceID,
         aphiaid, scientificName, taxonRank, phylum, subphylum, superclass, class, order, suborder, family, genus, species, samplingProtocol,
         lifeStage, coordUncertainty, basisOfRecord, occurrenceStatus, identificationQualifier, flags, dropped, absence)
# stage 1: the box subset goes straight to disk (without passing through R memory)
tmp_cx <- file.path(raiz, "intermediario/nivel2_caixas"); unlink(tmp_cx, recursive = TRUE)
write_dataset(x, tmp_cx, format = "parquet", max_rows_per_file = 2e6)
cat(format(Sys.time()), "- recorte das caixas gravado em", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")

# stage 2: one box at a time -> exact window (25 km cell), region, broad group and UNCLOS
for (d in c("nivel2_extracao", "nivel2_unclos")) unlink(file.path(raiz, "intermediario", d), recursive = TRUE)
cx <- open_dataset(tmp_cx); tot <- list(); n_un <- 0
for (k in seq_len(nrow(reg))) {
  y <- cx |> filter(lon >= reg$lon_min[k], lon <= reg$lon_max[k], lat >= reg$lat_min[k], lat <= reg$lat_max[k]) |> collect()
  c25 <- cellFromXY(g25, project(cbind(y$lon, y$lat), from = "EPSG:4326", to = "EPSG:6933")); rc <- rowColFromCell(g25, c25)
  y$janela <- paste((rc[, 1] - 1) %/% W + 1, (rc[, 2] - 1) %/% W + 1, sep = "_")
  y <- y[!is.na(c25) & y$janela == reg$janela[k], ]
  y$regiao <- reg$nome[k]
  y$grupo_amplo <- case_when(
    y$phylum %in% c("Echinodermata", "Porifera", "Chaetognatha", "Bryozoa", "Foraminifera") ~ y$phylum,
    y$subphylum %in% "Anthozoa" ~ "Anthozoa",
    y$class %in% c("Polychaeta", "Bivalvia", "Ascidiacea", "Thaliacea", "Gastropoda", "Ostracoda", "Cephalopoda", "Hydrozoa") ~ y$class,
    y$superclass %in% "Actinopteri" | y$class %in% "Elasmobranchii" ~ "Peixes",   # in OBIS, bony fishes: superclass Actinopteri, class Teleostei
    y$order %in% "Amphipoda" & y$suborder %in% "Hyperiidea" ~ "Hyperiidea",
    TRUE ~ y$order)
  un <- (y$genus %in% unclos_gen) | (y$family %in% unclos_fam) | (y$species %in% unclos_esp)
  n_un <- n_un + sum(un)
  write_dataset(y[!un, ], file.path(raiz, "intermediario/nivel2_extracao"), partitioning = "regiao", format = "parquet")
  if (any(un)) write_dataset(y[un, ] |> mutate(motivo = "unclos_anexo1"), file.path(raiz, "intermediario/nivel2_unclos"),
                             partitioning = "regiao", format = "parquet")
  tot[[k]] <- y[!un, ] |> group_by(regiao, grupo_amplo) |> summarise(registros = n(), ausencias = sum(absence, na.rm = TRUE),
          conjuntos = n_distinct(dataset_id), aparelho = round(100 * mean(!is.na(samplingProtocol) & samplingProtocol != ""), 1),
          estagio = round(100 * mean(!is.na(lifeStage) & lifeStage != ""), 1),
          profundidade = round(100 * mean(!is.na(depth_min) | !is.na(depth_max) | !is.na(depth)), 1), .groups = "drop")
  cat(format(Sys.time()), "-", reg$nome[k], ":", sum(!un), "registros\n")
  rm(y); gc(verbose = FALSE)
}
unlink(tmp_cx, recursive = TRUE)
tot <- bind_rows(tot)
write.csv(tot, file.path(raiz, "saida/nivel2_totais.csv"), row.names = FALSE)
cat("registros nas 8 regioes:", sum(tot$registros), "| UNCLOS excluidos:", n_un, "\n")
print(tot |> group_by(regiao) |> summarise(registros = sum(registros), grupos = n()) |> as.data.frame())
