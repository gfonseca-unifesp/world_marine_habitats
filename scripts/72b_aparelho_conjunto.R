# 72b_aparelho_conjunto.R  (Gate 3, item 4: sampling gear inferred from the dataset)
# For each level 2 dataset, searches the title and, if the title is not decisive, the abstract (OBIS API, cached in
# bruto/api_datasets/) for specific sampling-gear terms. The dataset is assigned "bentos" or "pelagico" only when the terms
# point in a single direction. Datasets whose type says nothing about the compartment are left without inference
# (telemetry, tagging, eDNA, aerial or bird and mammal observation, citizen science, museums, compilations).
# Used in 73_curadoria_n2.R: step (3b) of the compartment, only for records not decided by
# record gear, depth or taxon attribute (flags fonte = "aparelho do conjunto").
# Output: entrada/aparelho_conjunto.csv
suppressPackageStartupMessages({library(arrow); library(dplyr); library(jsonlite)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
cjs <- open_dataset(file.path(raiz, "intermediario/nivel2_extracao")) |> distinct(dataset_id) |> collect() |> pull()
meta <- function(id, campo) { a <- file.path(raiz, "bruto/api_datasets", paste0(id, ".json")); if (!file.exists(a)) return("")
  r <- fromJSON(a, simplifyVector = FALSE)$results; if (!length(r) || is.null(r[[1]][[campo]])) "" else paste(unlist(r[[1]][[campo]]), collapse = " ") }
d <- data.frame(dataset_id = cjs, titulo = vapply(cjs, meta, "", campo = "title"), resumo = vapply(cjs, meta, "", campo = "abstract"))

# specific sampling gears (without the generic words "benthic"/"pelagic", common in abstracts of any type)
GP <- "plankton net|bongo|wp-?2|wp-?3|multi ?net|mocness|neuston|\\bcpr\\b|continuous plankton|midwater trawl|mid-water trawl|pelagic trawl|ring net|vertical (haul|tow)|oblique tow|tucker|ikmt|isaacs|juday|hensen|manta net|ichthyoplankton|cone net|purse seine|calvet|calcofi|plankton survey|plankton tow|zooplankton"
GB <- "grab|van ?veen|day grab|smith.?mcintyre|box ?-?cor|corer|\\bcores?\\b|dredge|beam trawl|bottom trawl|otter trawl|groundfish|sledge|epibenthic|benthic survey|zoobenthos|macrobenthos|meiofauna|\\brov\\b|divers?\\b|dive imagery|scuba|quadrats?|swath|transects?|suction|airlift|sigsbee|agassiz|\\btraps?\\b|\\bpots?\\b|research vessel survey|survey trawl|4 seam|kelp forest|intertidal|rocky shore|tide ?pool|near-seafloor"
NAO_T <- "telemetry|tracking|tagging|\\btags?\\b|acoustic|edna|environmental dna|metabarcod|\\bcoi\\b|18s|aerial|citizen|inaturalist|naturewatch|museum|collection|fishnet|fishbase|database|atlas|contaminant|seabird|marine mammal|bycatch|compiled|snapshot|recorder|information system|port sampling|expedition report|oceanic 19|wildlife"
NAO_R <- "telemetry|acoustic receiver|tagged|tracking|seabird|marine mammal|cetacean|bycatch|environmental dna|edna"
tt <- tolower(d$titulo); ab <- tolower(d$resumo)
pt <- grepl(GP, tt); bt <- grepl(GB, tt); pa <- grepl(GP, ab); ba <- grepl(GB, ab)
d$aparelho_conjunto <- ifelse(pt != bt, ifelse(pt, "pelagico", "bentos"), ifelse(!pt & !bt & pa != ba, ifelse(pa, "pelagico", "bentos"), NA))
d$fonte <- ifelse(pt != bt, "titulo", ifelse(!is.na(d$aparelho_conjunto), "resumo", ifelse(pt | bt | (pa & ba), "termos nos dois sentidos", "sem termos")))
nt <- grepl(NAO_T, tt) | grepl(NAO_R, ab); d$aparelho_conjunto[nt] <- NA; d$fonte[nt] <- "tipo de conjunto sem aparelho"
d$resumo <- NULL
write.csv(d, file.path(raiz, "entrada/aparelho_conjunto.csv"), row.names = FALSE)
print(table(d$fonte, d$aparelho_conjunto, useNA = "ifany"))
