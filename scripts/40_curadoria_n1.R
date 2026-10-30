# 40_curadoria_n1.R  (Phase 4, filters 4.1 to 4.10, with the Gate 1 decisions)
# Each filter is a function that receives the records still without a reason and returns the quarantine reason (or NA).
# The filters run in cascade, in the order of the guide; a record leaves at the first filter it fails.
# Outputs (per focal group):
#   intermediario/n1_curado/grupo=<g>/      retained records, with cell_id 25/50/100 km, so_esforco, sem_profundidade,
#                                           known month/day, decimals and nivel_certeza (A/B/C)
#   intermediario/quarentena/n1/grupo=<g>/  discarded records, with the motivo column and the filter
#   intermediario/ausencias/n1/grupo=<g>/   declared absences
#   saida/cascata_filtros_n1.csv            entering, leaving and remaining per filter and group
suppressPackageStartupMessages({library(arrow); library(dplyr); library(terra)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
source(file.path(raiz, "scripts/funcoes_grade.R"))
ds  <- open_dataset(file.path(raiz, "intermediario/nivel1_extracao"))
cj  <- read.csv(file.path(raiz, "intermediario/conjunto_tipo_escopo.csv"))
grp <- read.csv(file.path(raiz, "saida/nivel1_definicao_grupos.csv")) |> distinct(grupo, compartimento)
for (d in c("n1_curado", "quarentena/n1", "ausencias/n1")) unlink(file.path(raiz, "intermediario", d), recursive = TRUE)

# ---- parameters (Gate 1 and D decisions) ----
# NO_ACCEPTED_NAME was removed from the list (stop of 09/10/2026): it flagged names accepted in WoRMS (e.g., Metridia lucens);
# replaced by the classification check in filter 4.5.
FLAGS_EXCLUI <- c("NO_COORD", "ZERO_COORD", "LAT_OUT_OF_RANGE", "LON_OUT_OF_RANGE",
                  "MIN_DEPTH_EXCEEDS_MAX", "DATE_IN_FUTURE", "DATE_BEFORE_MIN", "NOT_MARINE")
DEC_MIN_N1     <- 1                           # 1 decimal place: maximum error 0.05 degree (~5.6 km) < 1/4 of the 25 km cell
TERRA_M        <- -1000                       # shoredistance below this = on land (replaces ON_LAND)
CONJ_CORROMPIDO <- c("d29d7aa5-08e0-4c1f-9b1c-3b46b3b93563")
BASIS_OK       <- c("humanobservation", "materialsample", "preservedspecimen", "occurrence", "machineobservation",
                    "livingspecimen", "materialentity", "materialcitation")
INC_MAX_M      <- EASE2_BASE_M / 4            # 1/4 of the side of the 25 km cell (D11)
UNID_GENERO    <- c("Euphausiacea", "Chaetognatha")   # Gate 1, item 2
lic_ok <- with(cj, licenca %in% c("CC0", "CC-BY", "CC-BY-NC") |
  (licenca == "outra" & grepl("unrestricted|open.government|waived|public domain", tolower(licenca_texto))))
CONJ_LIC_OK <- cj$dataset_id[lic_ok]

# ---- single classification per AphiaID (the most frequent) ----
tx <- read.csv(file.path(raiz, "saida/taxa.csv")) |> group_by(aphiaid) |> slice_max(registros, n = 1, with_ties = FALSE) |>
  ungroup() |> select(aphiaid, taxonRank_u = taxonRank, phylum_u = phylum, class_u = class, order_u = order,
                      family_u = family, genus_u = genus, species_u = species)

# ---- parsing of dates that OBIS did not interpret (stop of 09/10/2026) ----
# E.g., "2007-09/ 2007-10" (space after the slash) ends up without date_year in OBIS. Here: remove spaces, split start and end,
# accept YYYY, YYYY-MM and YYYY-MM-DD (with or without time). The year is valid when start and end fall in the same year.
le_data <- function(txt) {
  t <- gsub("\\s+", "", txt); p <- strsplit(ifelse(is.na(t), "", t), "/", fixed = TRUE)
  ini <- vapply(p, function(v) if (length(v)) v[1] else "", ""); fim <- vapply(p, function(v) if (length(v) > 1 && nzchar(v[2])) v[2] else if (length(v)) v[1] else "", "")
  limite <- function(s, final) {
    s <- substr(s, 1, 10)
    d <- rep(as.Date(NA), length(s))
    a <- grepl("^\\d{4}-\\d{2}-\\d{2}$", s); d[a] <- as.Date(s[a], "%Y-%m-%d")
    b <- grepl("^\\d{4}-\\d{2}$", s); d[b] <- as.Date(paste0(s[b], "-01"), "%Y-%m-%d")
    if (final) d[b] <- as.Date(format(d[b] + 31, "%Y-%m-01"), "%Y-%m-%d") - 1
    c <- grepl("^\\d{4}$", s); d[c] <- as.Date(paste0(s[c], if (final) "-12-31" else "-01-01"), "%Y-%m-%d")
    d
  }
  di <- limite(ini, FALSE); df <- limite(fim, TRUE)
  # non-ISO formats (e.g., "7/12/2016", "02.08.2007"): year only, when there is a single 4-digit year in the text
  outro <- is.na(di) & grepl("(18|19|20)\\d{2}", t)
  if (any(outro)) {
    anos <- regmatches(t[outro], gregexpr("(18|19|20)\\d{2}", t[outro]))
    a1 <- vapply(anos, function(a) if (length(unique(a)) == 1) a[1] else NA_character_, "")
    di[outro] <- as.Date(paste0(a1, "-01-01"), "%Y-%m-%d"); df[outro] <- as.Date(paste0(a1, "-12-31"), "%Y-%m-%d")
  }
  ok <- !is.na(di) & !is.na(df) & df >= di
  data.frame(ano = ifelse(ok & format(di, "%Y") == format(df, "%Y"), as.integer(format(di, "%Y")), NA_integer_),
             ini_ms = ifelse(ok, as.numeric(as.POSIXct(di, tz = "UTC")) * 1000, NA_real_),
             fim_ms = ifelse(ok, as.numeric(as.POSIXct(df, tz = "UTC")) * 1000, NA_real_))
}

# ---- ocean masks at 25 km (seafloor in EBU; water in EMU at any level) ----
g25 <- ease2_grade("25km")
fundo <- read_parquet(file.path(raiz, "camadas/ebu/ebu_celula_25km.parquet"))$cell_id
emu25 <- rast(file.path(raiz, "camadas/emu/emu_25km.tif"))
agua_cel <- which(!is.na(values(app(emu25, function(v) any(!is.na(v)))))[, 1] & values(app(emu25, function(v) any(!is.na(v))))[, 1] == 1)
rcA <- rowColFromCell(g25, agua_cel); agua <- paste("25km", rcA[, 1], rcA[, 2], sep = "_")

cheio <- function(x) !is.na(x) & as.character(x) != ""
ms <- function(x) as.POSIXct(as.numeric(x) / 1000, origin = "1970-01-01", tz = "UTC")
decimais <- function(verb, num) {
  s <- ifelse(cheio(verb) & grepl("^-?[0-9]+([.,][0-9]+)?$", trimws(verb)), trimws(verb),
              sub("0+$", "", formatC(abs(num), format = "f", digits = 6)))
  ifelse(grepl("[.,]", s), nchar(sub("^[^.,]*[.,]", "", s)), 0L)
}
tem_flag <- function(fl, nomes) vapply(fl, function(f) any(f %in% nomes), TRUE)
primeira_flag <- function(fl, nomes) vapply(fl, function(f) { k <- f[f %in% nomes]; if (length(k)) k[1] else NA_character_ }, "")

# ---- filters: each returns a vector of reasons (NA = passes) ----
f41_geografia <- function(x, comp) {
  inval <- is.na(x$lon) | is.na(x$lat) | abs(x$lat) > 90 | abs(x$lon) > 180 | (x$lon == 0 & x$lat == 0) |
           tem_flag(x$flags, c("NO_COORD", "ZERO_COORD", "LAT_OUT_OF_RANGE", "LON_OUT_OF_RANGE"))
  m <- ifelse(inval, "coord_invalida", NA)
  ok <- is.na(m)
  cel <- rep(NA_character_, nrow(x))
  xy <- project(cbind(x$lon[ok], x$lat[ok]), from = "EPSG:4326", to = "EPSG:6933")
  c25 <- cellFromXY(g25, xy); rc <- rowColFromCell(g25, c25)
  cel[ok] <- ifelse(is.na(c25), NA, paste("25km", rc[, 1], rc[, 2], sep = "_"))
  mascara <- if (comp == "bentos") fundo else agua
  m[ok & !(cel %in% mascara)] <- "fora_oceano"
  attr(m, "cel") <- cel; m
}
f42_conjuntos <- function(x) ifelse(x$dataset_id %in% CONJ_CORROMPIDO, "conjunto_corrompido", NA)
f43_sinalizadores <- function(x) {
  m <- primeira_flag(x$flags, FLAGS_EXCLUI); m <- ifelse(is.na(m), NA, paste0("flag_", tolower(m)))
  m[is.na(m) & !is.na(x$shoredistance) & x$shoredistance < TERRA_M] <- "em_terra"; m
}
f44_status <- function(x) {
  aus <- x$absence %in% TRUE | tolower(x$occurrenceStatus) %in% "absent"
  idq <- grepl("cf\\.?|aff\\.?|\\?", x$identificationQualifier) | grepl(" cf\\.? | aff\\.? |\\?", x$scientificName)
  b <- gsub("[^a-z]", "", tolower(x$basisOfRecord))
  ifelse(aus, "ausencia", ifelse(idq, "id_incerta", ifelse(b == "fossilspecimen", "fossil",
    ifelse(b == "nomenclaturalchecklist", "basis_recusado", ifelse(cheio(b) & !(b %in% BASIS_OK), "basis_invalido", NA)))))
}
f45_taxon <- function(x) ifelse(is.na(x$aphiaid), "sem_aphia",
  ifelse(!cheio(x$phylum_u) & !cheio(x$class_u) & !cheio(x$order_u) & !cheio(x$family_u) & !cheio(x$genus_u),
         "sem_classificacao", NA))
f46_licenca <- function(x) ifelse(x$dataset_id %in% CONJ_LIC_OK, NA, "licenca_recusada")
f47_data <- function(x) ifelse(is.na(x$year), "sem_ano", NA)
f48_incerteza <- function(x) {
  inf <- !is.na(x$coordUncertainty)
  ifelse(inf & x$coordUncertainty > INC_MAX_M, "incerteza_alta", ifelse(!inf & x$decimais < DEC_MIN_N1, "precisao_baixa", NA))
}
f49_profundidade <- function(x) {
  p <- pmax(x$depth_max, x$depth_min, x$depth, na.rm = TRUE)
  tol <- pmax(20, 0.1 * pmax(x$bathymetry, 0, na.rm = TRUE))
  ruim <- (!is.na(x$depth_min) & !is.na(x$depth_max) & x$depth_min > x$depth_max) |
          (!is.na(p) & pmin(x$depth_min, x$depth_max, x$depth, na.rm = TRUE) < 0) |
          (!is.na(p) & !is.na(x$bathymetry) & x$bathymetry > 0 & p > x$bathymetry + tol)
  ifelse(ruim, "profundidade_incoerente", NA)
}
f410_duplicatas <- function(x) {
  score <- cheio(x$eventID) + !is.na(x$depth_min) + !is.na(x$depth_max) + !is.na(x$coordUncertainty) +
           cheio(x$samplingProtocol) + cheio(x$family) + cheio(x$genus) + cheio(x$lifeStage)
  o <- order(-score, x$record_id)
  m <- rep(NA_character_, nrow(x))
  k1 <- ifelse(cheio(x$occurrenceID), paste(x$dataset_id, x$occurrenceID), NA)
  d1 <- rep(FALSE, nrow(x)); d1[o] <- !is.na(k1[o]) & duplicated(k1[o], incomparables = NA)
  m[d1] <- "dup_interna"
  # (b) between datasets: the same occurrence (AphiaID, coordinate to 4 decimals, date, depth to 1 m) published
  #     in more than one dataset. The dataset with the most metadata for that key stays; the records of the others leave.
  #     Replicates within the same dataset (several samples at the same station and date) are NOT duplicates.
  pz <- round(pmax(x$depth_max, x$depth_min, x$depth, na.rm = TRUE))
  k2 <- paste(x$aphiaid, round(x$lon, 4), round(x$lat, 4), x$date_start, pz)
  vivo <- which(is.na(m) & !is.na(x$date_start))
  t2 <- data.frame(i = vivo, k = k2[vivo], cj = x$dataset_id[vivo], s = score[vivo])
  melhor <- t2 |> group_by(k, cj) |> summarise(s = mean(s), .groups = "drop") |> group_by(k) |>
    filter(n() > 1) |> arrange(desc(s), cj) |> summarise(fica = first(cj), .groups = "drop")
  t2 <- inner_join(t2, melhor, by = "k")
  m[t2$i[t2$cj != t2$fica]] <- "dup_entre_conjuntos"; m
}
FILTROS <- list("4.1 geografia" = "f41", "4.2 conjuntos" = "f42", "4.3 sinalizadores" = "f43", "4.4 status" = "f44",
                "4.5 taxon" = "f45", "4.6 licenca" = "f46", "4.7 data" = "f47", "4.8 incerteza" = "f48",
                "4.9 profundidade" = "f49", "4.10 duplicatas" = "f410")

cascata <- list()
for (g in grp$grupo) {
  comp <- grp$compartimento[grp$grupo == g]
  x <- ds |> filter(grupo == g) |> collect() |> left_join(tx, by = "aphiaid")
  x$decimais <- pmin(decimais(x$verbatimLon, x$lon), decimais(x$verbatimLat, x$lat))
  # dates that OBIS did not read: fill in year, start and end from the original eventDate
  x$data_relida <- FALSE
  # (1) separate year, month and day fields (e.g., CPR, without eventDate)
  sem <- which(is.na(x$year) & grepl("^\\d{4}$", trimws(x$year_txt)))
  if (length(sem)) {
    a <- as.integer(trimws(x$year_txt[sem])); m <- suppressWarnings(as.integer(x$month_txt[sem]))
    dd <- suppressWarnings(as.integer(x$day_txt[sem])); m[!m %in% 1:12] <- NA
    ini <- as.Date(sprintf("%04d-%02d-%02d", a, ifelse(is.na(m), 1L, m), ifelse(is.na(m) | is.na(dd), 1L, dd)), "%Y-%m-%d")
    fim_ano <- as.Date(sprintf("%04d-12-31", a), "%Y-%m-%d")
    fim_mes <- as.Date(format(as.Date(sprintf("%04d-%02d-01", a, ifelse(is.na(m), 1L, m)), "%Y-%m-%d") + 31, "%Y-%m-01"), "%Y-%m-%d") - 1
    fim <- ifelse(is.na(m), fim_ano, ifelse(is.na(dd), fim_mes, ini))
    okd <- !is.na(ini) & a >= 1700 & a <= 2026
    x$year[sem[okd]] <- a[okd]
    x$date_start[sem[okd]] <- bit64::as.integer64(as.numeric(as.POSIXct(ini[okd], tz = "UTC")) * 1000)
    x$date_end[sem[okd]] <- bit64::as.integer64(as.numeric(as.POSIXct(as.Date(fim[okd], origin = "1970-01-01"), tz = "UTC")) * 1000)
    x$data_relida[sem[okd]] <- TRUE
  }
  # (2) eventDate (or verbatimEventDate) as text
  txt <- ifelse(cheio(x$eventDate), x$eventDate, x$verbatimEventDate)
  sem <- which(is.na(x$year) & cheio(txt))
  if (length(sem)) {
    d <- le_data(txt[sem]); okd <- !is.na(d$ano)
    x$year[sem[okd]] <- d$ano[okd]
    x$date_start[sem[okd]] <- bit64::as.integer64(d$ini_ms[okd]); x$date_end[sem[okd]] <- bit64::as.integer64(d$fim_ms[okd])
    x$data_relida[sem[okd]] <- TRUE
  }
  n <- nrow(x); motivo <- rep(NA_character_, n); filtro <- rep(NA_character_, n)
  for (nm in names(FILTROS)) {
    vivo <- which(is.na(motivo)); entra <- length(vivo)
    m <- switch(FILTROS[[nm]],
      f41 = f41_geografia(x[vivo, ], comp), f42 = f42_conjuntos(x[vivo, ]), f43 = f43_sinalizadores(x[vivo, ]),
      f44 = f44_status(x[vivo, ]), f45 = f45_taxon(x[vivo, ]), f46 = f46_licenca(x[vivo, ]),
      f47 = f47_data(x[vivo, ]), f48 = f48_incerteza(x[vivo, ]), f49 = f49_profundidade(x[vivo, ]),
      f410 = f410_duplicatas(x[vivo, ]))
    if (FILTROS[[nm]] == "f41") { x$cell_25 <- NA_character_; x$cell_25[vivo] <- attr(m, "cel") }
    sai <- !is.na(m); motivo[vivo[sai]] <- m[sai]; filtro[vivo[sai]] <- nm
    por_motivo <- table(m[sai])
    cascata[[length(cascata) + 1]] <- data.frame(grupo = g, filtro = nm, entram = entra, saem = sum(sai),
      restam = entra - sum(sai), pct_saem = round(100 * sum(sai) / max(entra, 1), 2),
      motivos = paste(names(por_motivo), por_motivo, sep = "=", collapse = "; "))
  }
  x$motivo <- motivo; x$filtro <- filtro
  # retained records: attributes for the following phases
  r <- x[is.na(motivo), ]
  rc <- do.call(rbind, strsplit(r$cell_25, "_")); r25 <- as.integer(rc[, 2]); c25 <- as.integer(rc[, 3])
  r$cell_50 <- paste("50km", (r25 - 1) %/% 2 + 1, (c25 - 1) %/% 2 + 1, sep = "_")
  r$cell_100 <- paste("100km", (r25 - 1) %/% 4 + 1, (c25 - 1) %/% 4 + 1, sep = "_")
  r$bloco_2500 <- paste((r25 - 1) %/% 100 + 1, (c25 - 1) %/% 100 + 1, sep = "_")
  r$unidade_tax <- if (g %in% UNID_GENERO) "genero" else "familia"
  r$taxon_unidade <- if (g %in% UNID_GENERO) r$genus_u else r$family_u
  r$so_esforco <- !cheio(r$taxon_unidade)
  r$sem_profundidade <- is.na(r$depth_min) & is.na(r$depth_max) & is.na(r$depth)
  r$basis_citacao <- gsub("[^a-z]", "", tolower(r$basisOfRecord)) == "materialcitation"
  di <- ms(r$date_start); dfim <- ms(r$date_end)
  r$mes <- ifelse(!is.na(di) & format(di, "%Y-%m") == format(dfim, "%Y-%m"), as.integer(format(di, "%m")), NA)
  r$dia_conhecido <- !is.na(di) & format(di, "%Y-%m-%d") == format(dfim, "%Y-%m-%d")
  pos_inf <- !is.na(r$coordUncertainty)
  r$nivel_certeza <- ifelse(pos_inf & r$dia_conhecido, "A", ifelse(!is.na(r$mes), "B", "C"))
  r$evento <- ifelse(cheio(r$eventID), paste(r$dataset_id, r$eventID), paste(r$dataset_id, r$date_start, round(r$lon, 4), round(r$lat, 4)))
  r$compartimento <- comp
  write_dataset(r |> mutate(grupo = g), file.path(raiz, "intermediario/n1_curado"), partitioning = "grupo", format = "parquet")
  q <- x[!is.na(motivo) & motivo != "ausencia", ]
  if (nrow(q)) write_dataset(q |> mutate(grupo = g), file.path(raiz, "intermediario/quarentena/n1"), partitioning = "grupo", format = "parquet")
  a <- x[motivo %in% "ausencia", ]
  if (nrow(a)) write_dataset(a |> mutate(grupo = g), file.path(raiz, "intermediario/ausencias/n1"), partitioning = "grupo", format = "parquet")
  stopifnot(nrow(r) + nrow(q) + nrow(a) == n)              # the cascade balances
  cat(format(Sys.time()), "-", g, ": extraidos", n, "| retidos", nrow(r), "| quarentena", nrow(q), "| ausencias", nrow(a), "\n")
  rm(x, r, q, a); gc(verbose = FALSE)
}
cas <- bind_rows(cascata)
write.csv(cas, file.path(raiz, "saida/cascata_filtros_n1.csv"), row.names = FALSE)
alerta <- cas |> filter(pct_saem > 30, !filtro %in% c("4.3 sinalizadores"))
if (nrow(alerta)) { cat("\n*** PARADA: filtros que tiraram mais de 30% dos registros restantes ***\n"); print(alerta) }
