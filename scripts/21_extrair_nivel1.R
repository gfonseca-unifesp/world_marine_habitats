# 21_extrair_nivel1.R  (task 2.1)
# Extracts from the local OBIS copy all records of the level 1 focal taxa, in a single pass,
# including failed records (dropped) and absences. No quality filter is applied here.
# Output: intermediario/nivel1_extracao/ (parquet partitioned by grupo) and saida/nivel1_totais.csv
suppressPackageStartupMessages({library(arrow); library(dplyr)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
dir.create(file.path(raiz, "intermediario"), showWarnings = FALSE); dir.create(file.path(raiz, "saida"), showWarnings = FALSE)

# Focal groups: taxonomic level used and compartment (strategy E2; current WoRMS classification)
focais <- tribble(
  ~grupo,                      ~compartimento, ~campo,      ~taxon,
  "Polychaeta",                "bentos",       "class",     "Polychaeta",
  "Bivalvia",                  "bentos",       "class",     "Bivalvia",
  "Echinodermata",             "bentos",       "phylum",    "Echinodermata",
  "Anthozoa",                  "bentos",       "subphylum", "Anthozoa",
  "Porifera",                  "bentos",       "phylum",    "Porifera",
  "Peracarida bentonicos",     "bentos",       "order",     "Cumacea",
  "Peracarida bentonicos",     "bentos",       "order",     "Isopoda",
  "Peracarida bentonicos",     "bentos",       "order",     "Tanaidacea",
  "Bryozoa e Ascidiacea",      "bentos",       "phylum",    "Bryozoa",
  "Bryozoa e Ascidiacea",      "bentos",       "class",     "Ascidiacea",
  "Calanoida",                 "pelagico",     "order",     "Calanoida",
  "Euphausiacea",              "pelagico",     "order",     "Euphausiacea",
  "Chaetognatha",              "pelagico",     "phylum",    "Chaetognatha",
  "Hyperiidea",                "pelagico",     "suborder",  "Hyperiidea",
  "Siphonophorae e Thaliacea", "pelagico",     "order",     "Siphonophorae",
  "Siphonophorae e Thaliacea", "pelagico",     "class",     "Thaliacea")
write.csv(focais, file.path(raiz, "saida/nivel1_definicao_grupos.csv"), row.names = FALSE)
tx <- function(cp) focais$taxon[focais$campo == cp]
f_phylum <- tx("phylum"); f_subphylum <- tx("subphylum"); f_class <- tx("class")
f_order <- tx("order"); f_suborder <- tx("suborder")

arqs <- list.files(file.path(raiz, "bruto/occurrence"), full.names = TRUE)
ds <- open_dataset(arqs, format = "parquet", unify_schemas = FALSE)
t0 <- Sys.time()
q <- ds |>
  mutate(phylum = interpreted$phylum, subphylum = interpreted$subphylum, class = interpreted$class,
         order = interpreted$order, suborder = interpreted$suborder) |>
  filter(phylum %in% f_phylum | subphylum %in% f_subphylum | class %in% f_class |
         order %in% f_order | suborder %in% f_suborder) |>
  mutate(grupo = case_when(
           phylum == "Echinodermata" ~ "Echinodermata", phylum == "Porifera" ~ "Porifera",
           phylum == "Chaetognatha" ~ "Chaetognatha", subphylum == "Anthozoa" ~ "Anthozoa",
           class == "Polychaeta" ~ "Polychaeta", class == "Bivalvia" ~ "Bivalvia",
           phylum == "Bryozoa" | class == "Ascidiacea" ~ "Bryozoa e Ascidiacea",
           order %in% c("Cumacea", "Isopoda", "Tanaidacea") ~ "Peracarida bentonicos",
           order == "Calanoida" ~ "Calanoida", order == "Euphausiacea" ~ "Euphausiacea",
           suborder == "Hyperiidea" ~ "Hyperiidea",
           order == "Siphonophorae" | class == "Thaliacea" ~ "Siphonophorae e Thaliacea"),
         record_id = `_id`, lon = interpreted$decimalLongitude, lat = interpreted$decimalLatitude,
         eventDate = interpreted$eventDate, date_start = interpreted$date_start, date_end = interpreted$date_end,
         year = interpreted$date_year, depth_min = interpreted$minimumDepthInMeters,
         depth_max = interpreted$maximumDepthInMeters, depth = interpreted$depth,
         bathymetry = interpreted$bathymetry, shoredistance = interpreted$shoredistance,
         eventID = interpreted$eventID, parentEventID = interpreted$parentEventID,
         occurrenceID = interpreted$occurrenceID, aphiaid = interpreted$aphiaid,
         scientificName = interpreted$scientificName, taxonRank = interpreted$taxonRank,
         family = interpreted$family, genus = interpreted$genus, species = interpreted$species,
         samplingProtocol = interpreted$samplingProtocol, lifeStage = interpreted$lifeStage,
         coordUncertainty = interpreted$coordinateUncertaintyInMeters, basisOfRecord = interpreted$basisOfRecord,
         occurrenceStatus = interpreted$occurrenceStatus, identificationQualifier = interpreted$identificationQualifier,
         license = interpreted$license, marine = interpreted$marine, brackish = interpreted$brackish,
         verbatimLon = interpreted$verbatimLongitude, verbatimLat = interpreted$verbatimLatitude,
         year_txt = interpreted$year, month_txt = interpreted$month, day_txt = interpreted$day,
         verbatimEventDate = interpreted$verbatimEventDate) |>
  select(grupo, record_id, dataset_id, lon, lat, verbatimLon, verbatimLat, eventDate, date_start, date_end, year,
         depth_min, depth_max, depth, bathymetry, shoredistance, eventID, parentEventID, occurrenceID,
         aphiaid, scientificName, taxonRank, phylum, subphylum, class, order, suborder, family, genus, species,
         samplingProtocol, lifeStage, coordUncertainty, basisOfRecord, occurrenceStatus,
         identificationQualifier, license, marine, brackish, year_txt, month_txt, day_txt, verbatimEventDate,
         flags, dropped, absence)
dest <- file.path(raiz, "intermediario/nivel1_extracao")
if (dir.exists(dest)) unlink(dest, recursive = TRUE)
write_dataset(q, dest, format = "parquet", partitioning = "grupo", max_rows_per_file = 5e6)
cat(format(Sys.time()), "- extracao gravada em", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")

tot <- open_dataset(dest) |> group_by(grupo) |>
  summarise(registros = n(), dropped = sum(dropped, na.rm = TRUE), ausencias = sum(absence, na.rm = TRUE),
            conjuntos = n_distinct(dataset_id)) |> collect() |> arrange(grupo)
tot$compartimento <- focais$compartimento[match(tot$grupo, focais$grupo)]
write.csv(tot, file.path(raiz, "saida/nivel1_totais.csv"), row.names = FALSE)
print(as.data.frame(tot)); cat("TOTAL:", sum(tot$registros), "\n")
