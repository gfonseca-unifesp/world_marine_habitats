# 21b_extrair_nivel1_b.R  (option B of 30/10/2026: widen the level-1 pelagic with the groups that gained signal at level 2)
# Extracts from the local OBIS copy, in a single pass and without quality filters, the families that entrada/taxon_compartimento.csv
# (approved at Gate 3) assigns to the pelagic groups "peixes pelágicos pequenos" and "micronécton" (Euphausiacea excluded, already a
# level-1 group). Gelatinous plankton adds no new group: its families are Thaliacea, already in level 1, plus one Hydrozoa family.
# The family list comes from the level-2 windows, so pelagic families absent from those windows are not extracted.
# Written apart from the approved level-1 extraction, which stays untouched.
# Outputs: intermediario/nivel1b_extracao/ (parquet partitioned by grupo), saida/nivel1b_totais.csv
suppressPackageStartupMessages({library(arrow); library(dplyr)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
tc <- read.csv(file.path(raiz, "entrada/taxon_compartimento.csv"), encoding = "UTF-8")
f_peixes <- unique(tc$family[tc$compartimento == "pelagico" & tc$grupo_ecologico == "peixes pelágicos pequenos"])
f_micro <- unique(tc$family[tc$compartimento == "pelagico" & tc$grupo_ecologico == "micronécton" & tc$grupo_amplo != "Euphausiacea"])
cat("familias: peixes", length(f_peixes), "| micronecton", length(f_micro), "\n")
ds <- open_dataset(list.files(file.path(raiz, "bruto/occurrence"), full.names = TRUE), format = "parquet", unify_schemas = FALSE)
t0 <- Sys.time()
q <- ds |>
  mutate(phylum = interpreted$phylum, subphylum = interpreted$subphylum, class = interpreted$class, order = interpreted$order, suborder = interpreted$suborder, family = interpreted$family) |>
  filter(family %in% c(f_peixes, f_micro)) |>
  mutate(grupo = ifelse(family %in% f_micro, "micronecton", "peixes pelagicos pequenos"),
         record_id = `_id`, lon = interpreted$decimalLongitude, lat = interpreted$decimalLatitude,
         eventDate = interpreted$eventDate, date_start = interpreted$date_start, date_end = interpreted$date_end,
         year = interpreted$date_year, depth_min = interpreted$minimumDepthInMeters,
         depth_max = interpreted$maximumDepthInMeters, depth = interpreted$depth,
         bathymetry = interpreted$bathymetry, shoredistance = interpreted$shoredistance,
         eventID = interpreted$eventID, parentEventID = interpreted$parentEventID,
         occurrenceID = interpreted$occurrenceID, aphiaid = interpreted$aphiaid,
         scientificName = interpreted$scientificName, taxonRank = interpreted$taxonRank,
         genus = interpreted$genus, species = interpreted$species,
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
dest <- file.path(raiz, "intermediario/nivel1b_extracao")
if (dir.exists(dest)) unlink(dest, recursive = TRUE)
write_dataset(q, dest, format = "parquet", partitioning = "grupo", max_rows_per_file = 5e6)
cat(format(Sys.time()), "- extracao gravada em", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
tot <- open_dataset(dest) |> group_by(grupo) |>
  summarise(registros = n(), dropped = sum(dropped, na.rm = TRUE), ausencias = sum(absence, na.rm = TRUE), conjuntos = n_distinct(dataset_id), desde_1993 = sum(year >= 1993, na.rm = TRUE)) |> collect()
cj <- read.csv(file.path(raiz, "intermediario/conjunto_tipo_escopo.csv")); novos <- open_dataset(dest) |> distinct(dataset_id) |> collect()
tot$conjuntos_todos <- nrow(novos); tot$conjuntos_fora_da_tabela <- sum(!novos$dataset_id %in% cj$dataset_id)
write.csv(tot, file.path(raiz, "saida/nivel1b_totais.csv"), row.names = FALSE); print(as.data.frame(tot)); cat(format(Sys.time()), "FIM\n")
