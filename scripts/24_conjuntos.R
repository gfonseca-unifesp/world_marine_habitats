# 24_conjuntos.R  (task 2.4)
# Metadata, type and scope of each dataset with records at level 1.
#  - Metadata from the OBIS API (/v3/dataset/<id>): title, institution, license, citation, DOI, abstract, keywords.
#  - Standardized license (CC0, CC-BY, CC-BY-NC, outra, ausente) from the intellectualrights text.
#  - Type: (c) telemetry, (d) species list, (a) event-based with effort, (b) occurrences only.
#    Text and data rules; ambiguous cases are flagged for the user's decision at Gate 1.
#  - Documented default depth (e.g. Continuous Plankton Recorder ~7 m), used in rule D4.3.
# Output: bruto/obis_conjuntos.csv (raw metadata) and intermediario/conjunto_tipo_escopo.csv
suppressPackageStartupMessages({library(arrow); library(dplyr); library(jsonlite)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
ds <- open_dataset(file.path(raiz, "intermediario/nivel1_extracao"))
cheio <- function(x) !is.na(x) & as.character(x) != ""

# Scope and indicators per dataset, from the records themselves
por_cj <- ds |>
  mutate(tem_ev = !is.na(eventID) & eventID != "", tem_prof = !is.na(depth_min) | !is.na(depth_max)) |>
  group_by(dataset_id) |>
  summarise(registros = n(), ausencias = sum(absence, na.rm = TRUE), pct_eventID = mean(tem_ev),
            pct_prof = mean(tem_prof), grupos = n_distinct(grupo), ano_min = min(year, na.rm = TRUE),
            ano_max = max(year, na.rm = TRUE)) |> collect()
esc <- ds |> distinct(dataset_id, grupo) |> collect() |> group_by(dataset_id) |>
  summarise(grupos_focais = paste(sort(grupo), collapse = "; "))
bas <- ds |> count(dataset_id, basisOfRecord) |> collect() |> group_by(dataset_id) |>
  slice_max(n, n = 1, with_ties = FALSE) |> select(dataset_id, basis_principal = basisOfRecord)
por_cj <- por_cj |> left_join(esc, by = "dataset_id") |> left_join(bas, by = "dataset_id")
cat("conjuntos com registros focais:", nrow(por_cj), "\n")

# Metadata from the API (cached on disk, so it can be rerun)
cache <- file.path(raiz, "bruto/api_datasets"); dir.create(cache, showWarnings = FALSE)
pega <- function(id) {
  arq <- file.path(cache, paste0(id, ".json"))
  if (!file.exists(arq)) for (t in 1:3) {
    ok <- tryCatch({ curl::curl_download(paste0("https://api.obis.org/v3/dataset/", id), arq, quiet = TRUE); TRUE },
                   error = function(e) FALSE)
    if (ok) break; Sys.sleep(5 * t)
  }
  if (!file.exists(arq)) return(NULL)
  r <- fromJSON(arq, simplifyVector = FALSE)$results
  if (!length(r)) return(NULL); r <- r[[1]]
  s <- function(v) if (is.null(v)) NA_character_ else paste(unlist(v), collapse = "; ")
  data.frame(dataset_id = id, titulo = s(r$title), instituicao = s(sapply(r$institutes, function(i) i$name)),
             licenca_texto = s(r$intellectualrights), citacao = s(r$citation), doi = s(r$doi),
             url = s(r$url), publicado = s(r$published), palavras = s(sapply(r$keywords, function(k) k$keyword)),
             resumo = substr(s(r$abstract), 1, 3000))
}
meta <- bind_rows(lapply(por_cj$dataset_id, pega))
write.csv(meta, file.path(raiz, "bruto/obis_conjuntos.csv"), row.names = FALSE)

lic <- function(t) {
  t <- tolower(ifelse(is.na(t), "", t))
  case_when(t == "" ~ "ausente",
            grepl("cc0|public domain|publicdomain|zero", t) ~ "CC0",
            grepl("non.?commercial|by.?nc|-nc|nc 4|nc-", t) ~ "CC-BY-NC",
            grepl("attribution|cc.?by|by 4|by/4", t) ~ "CC-BY",
            TRUE ~ "outra")
}
txt <- function(...) tolower(paste(..., sep = " "))
out <- por_cj |> left_join(meta, by = "dataset_id") |>
  mutate(licenca = lic(licenca_texto),
         t = txt(titulo, palavras, substr(resumo, 1, 1500)),
         sinal_telemetria = grepl("telemetr|satellite tag|satellite track|argos|gps track|tracking|biologging|acoustic tag", t),
         sinal_lista = grepl("checklist|check-list|species list|inventory|list of species|catalogue of|catalog of", t),
         tipo = case_when(sinal_telemetria ~ "c_telemetria",
                          sinal_lista & pct_eventID < 0.5 ~ "d_lista_especies",
                          pct_eventID >= 0.5 ~ "a_eventos",
                          TRUE ~ "b_ocorrencias"),
         incerto = (sinal_telemetria & grupos > 0) | (sinal_lista & pct_eventID >= 0.5) |
                   is.na(titulo),
         cpr = grepl("continuous plankton recorder|continuous zooplankton recorder|\\(cpr\\)|cpr survey", tolower(titulo)),
         prof_padrao_m = ifelse(cpr, 7, NA_real_),
         fonte_prof_padrao = ifelse(cpr, "CPR: amostragem a ~7 m (documentacao do programa)", NA_character_)) |>
  select(-t)
write.csv(out, file.path(raiz, "intermediario/conjunto_tipo_escopo.csv"), row.names = FALSE)
cat("\nTipos:\n"); print(table(out$tipo)); cat("\nLicencas (conjuntos):\n"); print(table(out$licenca))
cat("\nIncertos para o Gate 1:", sum(out$incerto), "| CPR:", sum(out$cpr), "| sem metadados:", sum(is.na(out$titulo)), "\n")
