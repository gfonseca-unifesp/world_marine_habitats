# 21c_extrair_nivel1_c.R  (Gate 12, 13/11/2026: test fauna for the approved maps, chosen from the OBIS inventory of 25/25b)
# The new groups test the approved habitats and united maps; they do not rebuild them. Extracted in a single pass, without quality
# filters, apart from the approved extractions:
#   benthos:  "peixes demersais" = families that entrada/taxon_compartimento.csv (Gate 3) puts in the ecological group "peixes demersais";
#             "Decapoda bentonicos" = Decapoda families with compartment "bentos" in that table ("ambiguo" families left out);
#   pelagic:  "Cyclopoida pelagicos" = Oithonidae, Corycaeidae, Oncaeidae, Sapphirinidae, Lubbockiidae (the order also holds parasitic and
#             benthic families); "Appendicularia" (class); "Pteropoda" (orders Pteropoda, Thecosomata, Gymnosomata).
# The pelagic groups get the pelagic-gear filter of Gate 11 (filter 4.0, 40_curadoria_n1.R); the benthic ones none, as the level-1 benthos.
# Outputs: intermediario/nivel1c_extracao/ (parquet partitioned by grupo), saida/nivel1c_totais.csv, saida/nivel1c_definicao_grupos.csv
suppressPackageStartupMessages({library(arrow); library(dplyr)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
tc <- read.csv(file.path(raiz, "entrada/taxon_compartimento.csv"), encoding = "UTF-8")
f_dem <- unique(tc$family[tc$grupo_ecologico == "peixes demersais"])
f_dec <- unique(tc$family[tc$order == "Decapoda" & tc$compartimento == "bentos"])
f_cyc <- c("Oithonidae", "Corycaeidae", "Oncaeidae", "Sapphirinidae", "Lubbockiidae"); o_pte <- c("Pteropoda", "Thecosomata", "Gymnosomata")
write.csv(data.frame(grupo = c("peixes demersais", "Decapoda bentonicos", "Cyclopoida pelagicos", "Appendicularia", "Pteropoda"),
  compartimento = c("bentos", "bentos", "pelagico", "pelagico", "pelagico")), file.path(raiz, "saida/nivel1c_definicao_grupos.csv"), row.names = FALSE)
cat("familias: peixes demersais", length(f_dem), "| Decapoda bentonicos", length(f_dec), "\n")
ds <- open_dataset(list.files(file.path(raiz, "bruto/occurrence"), full.names = TRUE), format = "parquet", unify_schemas = FALSE)
t0 <- Sys.time()
q <- ds |>
  mutate(phylum = interpreted$phylum, subphylum = interpreted$subphylum, class = interpreted$class, order = interpreted$order, suborder = interpreted$suborder, family = interpreted$family) |>
  filter(family %in% c(f_dem, f_dec, f_cyc) | class == "Appendicularia" | order %in% o_pte) |>
  mutate(grupo = case_when(class == "Appendicularia" ~ "Appendicularia", order %in% o_pte ~ "Pteropoda", family %in% f_cyc ~ "Cyclopoida pelagicos",
                           family %in% f_dec ~ "Decapoda bentonicos", TRUE ~ "peixes demersais"),
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
dest <- file.path(raiz, "intermediario/nivel1c_extracao")
if (dir.exists(dest)) unlink(dest, recursive = TRUE)
write_dataset(q, dest, format = "parquet", partitioning = "grupo", max_rows_per_file = 5e6)
cat(format(Sys.time()), "- extracao gravada em", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
tot <- open_dataset(dest) |> group_by(grupo) |>
  summarise(registros = n(), dropped = sum(dropped, na.rm = TRUE), ausencias = sum(absence, na.rm = TRUE), conjuntos = n_distinct(dataset_id), desde_1993 = sum(year >= 1993, na.rm = TRUE)) |> collect()
cj <- read.csv(file.path(raiz, "intermediario/conjunto_tipo_escopo.csv")); novos <- open_dataset(dest) |> distinct(dataset_id) |> collect()
tot$conjuntos_todos <- nrow(novos); tot$conjuntos_fora_da_tabela <- sum(!novos$dataset_id %in% cj$dataset_id)
write.csv(tot, file.path(raiz, "saida/nivel1c_totais.csv"), row.names = FALSE); print(as.data.frame(tot)); cat(format(Sys.time()), "FIM\n")
