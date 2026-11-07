# 24b_conjuntos_b.R  (option B; Gate 11, 06/11/2026: metadata, type and pelagic-gear filter of the datasets of 21b)
# Same metadata and type rules as 24_conjuntos.R, on the option-B extraction (API answers cached in bruto/api_datasets). Taxon table as 23_taxa.R.
# Gate 11 filters, record by record, applied in 40_curadoria_n1.R when WMH_N1B=1 (filter 4.0):
#  1 gear: keep records whose sampling text (samplingProtocol, or else dataset title and abstract) names a pelagic gear (midwater or
#    pelagic trawl, plankton or ichthyoplankton net, opening-closing net, MOCNESS, Bongo, neuston, IKMT, multinet); drop bottom trawl,
#    visual census, seine, gillnet and longline, and anything with no pelagic-gear evidence;
#  2 the Falkland Islands fisheries survey is dropped;
#  3 micronekton = mesopelagic fishes and oceanic squids: Loliginidae and Sepiolidae (coastal or neritic) are dropped;
#  4 taxonomic unit: genus for the new groups (40_curadoria_n1.R, UNID_GENERO).
# Outputs: intermediario/conjunto_tipo_escopo_b.csv, saida/taxa_b.csv, saida/nivel1b_filtro_aparelho.csv
suppressPackageStartupMessages({library(arrow); library(dplyr); library(jsonlite)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
ds <- open_dataset(file.path(raiz, "intermediario/nivel1b_extracao"))
por_cj <- ds |> mutate(tem_ev = !is.na(eventID) & eventID != "", tem_prof = !is.na(depth_min) | !is.na(depth_max)) |> group_by(dataset_id) |>
  summarise(registros = n(), ausencias = sum(absence, na.rm = TRUE), pct_eventID = mean(tem_ev), pct_prof = mean(tem_prof), grupos = n_distinct(grupo),
            ano_min = min(year, na.rm = TRUE), ano_max = max(year, na.rm = TRUE)) |> collect()
cache <- file.path(raiz, "bruto/api_datasets")
pega <- function(id) { arq <- file.path(cache, paste0(id, ".json"))
  if (!file.exists(arq)) for (t in 1:3) { ok <- tryCatch({ curl::curl_download(paste0("https://api.obis.org/v3/dataset/", id), arq, quiet = TRUE); TRUE }, error = function(e) FALSE); if (ok) break; Sys.sleep(5 * t) }
  if (!file.exists(arq)) return(NULL)
  r <- tryCatch(fromJSON(arq, simplifyVector = FALSE)$results, error = function(e) NULL); if (!length(r)) return(NULL); r <- r[[1]]
  s <- function(v) if (is.null(v)) NA_character_ else paste(unlist(v), collapse = "; ")
  data.frame(dataset_id = id, titulo = s(r$title), licenca_texto = s(r$intellectualrights), palavras = s(sapply(r$keywords, function(k) k$keyword)), resumo = substr(s(r$abstract), 1, 3000)) }
meta <- bind_rows(lapply(por_cj$dataset_id, pega)); cat(format(Sys.time()), "metadados:", nrow(meta), "de", nrow(por_cj), "\n")
lic <- function(t) { t <- tolower(ifelse(is.na(t), "", t)); case_when(t == "" ~ "ausente", grepl("cc0|public domain|publicdomain|zero", t) ~ "CC0",
  grepl("non.?commercial|by.?nc|-nc|nc 4|nc-", t) ~ "CC-BY-NC", grepl("attribution|cc.?by|by 4|by/4", t) ~ "CC-BY", TRUE ~ "outra") }
PEL <- "midwater|mid-water|pelagic trawl|plankton|ichthyoplankton|bongo|neuston|mocness|multinet|multi-net|manta|ikmt|isaacs|opening.?closing|oblique|calcofi|ring net|wp-?2|tucker|rmt|cobb|hydrobios"
FUNDO <- "bottom trawl|demersal|groundfish|ground fish|otter trawl|beam trawl|benthic|visual|census|diver|scuba|video|rov\\b|seine|gillnet|gill net|longline|long-line|hook|trap|kelp|reef fish|fishery|fisheries observer|landings"
tb <- por_cj |> left_join(meta, by = "dataset_id") |> mutate(licenca = lic(licenca_texto), t = tolower(paste(titulo, palavras, substr(resumo, 1, 1500))),
  sinal_telemetria = grepl("telemetr|satellite tag|satellite track|argos|gps track|tracking|biologging|acoustic tag", t),
  sinal_lista = grepl("checklist|check-list|species list|inventory|list of species|catalogue of|catalog of", t),
  tipo = case_when(sinal_telemetria ~ "c_telemetria", sinal_lista & pct_eventID < 0.5 ~ "d_lista_especies", pct_eventID >= 0.5 ~ "a_eventos", TRUE ~ "b_ocorrencias"),
  cj_pelagico = grepl(PEL, t), cj_fundo = grepl(FUNDO, t), malvinas = grepl("falkland", tolower(titulo)), prof_padrao_m = NA_real_, cpr = grepl("continuous plankton recorder", tolower(titulo))) |> select(-t)
write.csv(tb, file.path(raiz, "intermediario/conjunto_tipo_escopo_b.csv"), row.names = FALSE)
taxa <- ds |> filter(!is.na(aphiaid)) |> group_by(aphiaid, scientificName, taxonRank, phylum, subphylum, class, order, suborder, family, genus, species, grupo) |> summarise(registros = n(), .groups = "drop") |> collect()
write.csv(taxa, file.path(raiz, "saida/taxa_b.csv"), row.names = FALSE)
print(table(tb$tipo)); print(table(tb$licenca)); cat("conjuntos com sinal pelagico:", sum(tb$cj_pelagico), "| de fundo:", sum(tb$cj_fundo), "| Malvinas:", paste(substr(tb$titulo[tb$malvinas], 1, 70), collapse = "; "), "\n")
cat(format(Sys.time()), "FIM\n")
