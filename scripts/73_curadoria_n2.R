# 73_curadoria_n2.R  (task 7.3: level 2 curation and compartment per record)
# Same filters as Phase 4 (40_curadoria_n1.R), with the decisions of Gate 1 and of the 09/10/2026 stop, and:
#  - 4.8 precision: position error = reported uncertainty or, without it, half a unit of the last decimal place (degrees -> km).
#        Each record keeps the finest resolution the error allows (error <= 1/4 of the side: 6.25 / 12.5 / 25 km);
#        above 1/4 of 25 km it goes to quarantine. Level 2 rule D11: 2 decimal places only count for 12.5 and
#        6.25 km if the year is >= 2000; before that, the record is limited to 25 km.
#  - compartment per record (E5.3): focal taxa -> level 1 compartment; others, in order:
#        (1) gear (samplingProtocol or eMoF); (2) depth relative to the seabed (tolerance D12);
#        (3) taxon attribute in the approved table (entrada/taxon_compartimento.csv), only for adult or unstaged;
#        (3b) gear inferred from the dataset title/abstract (entrada/aparelho_conjunto.csv; Gate 3, item 4),
#             also for larvae, column fonte_compartimento = "aparelho do conjunto";
#        (4) larva without (1), (2) and (3b) -> quarantine "estagio_incerto"; others undecided -> "compartimento_incerto".
#        Families "fora" (Molidae) -> quarantine "taxon_fora".
#  - ecological group: from the table; egg or larva in the pelagic -> "meroplancton e ictioplancton".
# Outputs per region: intermediario/n2_curado/, intermediario/quarentena/n2/, intermediario/ausencias/n2/,
#                    saida/cascata_filtros_n2.csv
suppressPackageStartupMessages({library(arrow); library(dplyr); library(terra); library(jsonlite)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
source(file.path(raiz, "scripts/funcoes_grade.R"))
ds  <- open_dataset(file.path(raiz, "intermediario/nivel2_extracao"))
emof <- read_parquet(file.path(raiz, "intermediario/nivel2_aparelho_emof.parquet"))
tab <- read.csv(file.path(raiz, "entrada/taxon_compartimento.csv"))
APAR_CONJ <- read.csv(file.path(raiz, "entrada/aparelho_conjunto.csv"))          # from 72b_aparelho_conjunto.R
for (d in c("n2_curado", "quarentena/n2", "ausencias/n2")) unlink(file.path(raiz, "intermediario", d), recursive = TRUE)

# ---- dataset licenses (OBIS API, cached; includes the new level 2 datasets) ----
cache <- file.path(raiz, "bruto/api_datasets")
lic_de <- function(id) {
  arq <- file.path(cache, paste0(id, ".json"))
  if (!file.exists(arq)) tryCatch(curl::curl_download(paste0("https://api.obis.org/v3/dataset/", id), arq, quiet = TRUE), error = function(e) NULL)
  if (!file.exists(arq)) return(NA_character_)
  r <- fromJSON(arq, simplifyVector = FALSE)$results; if (!length(r)) return(NA_character_)
  t <- r[[1]]$intellectualrights; if (is.null(t)) NA_character_ else paste(unlist(t), collapse = " ")
}
cjs <- ds |> distinct(dataset_id) |> collect() |> pull()
lic_txt <- setNames(vapply(cjs, lic_de, ""), cjs)
lt <- tolower(ifelse(is.na(lic_txt), "", lic_txt))
lic_ok <- (grepl("cc0|public domain|publicdomain|zero|attribution|cc.?by|by 4|by/4|non.?commercial|by.?nc", lt) |
           grepl("unrestricted|open.government|waived", lt)) & !grepl("restricted", gsub("unrestricted", "", lt)) & lt != ""
CONJ_LIC_OK <- cjs[lic_ok]
cat("conjuntos:", length(cjs), "| licenca aprovada:", length(CONJ_LIC_OK), "\n")

# ---- parameters (same as level 1, unless noted) ----
FLAGS_EXCLUI <- c("NO_COORD", "ZERO_COORD", "LAT_OUT_OF_RANGE", "LON_OUT_OF_RANGE",
                  "MIN_DEPTH_EXCEEDS_MAX", "DATE_IN_FUTURE", "DATE_BEFORE_MIN", "NOT_MARINE")
TERRA_M <- -1000; CONJ_CORROMPIDO <- c("d29d7aa5-08e0-4c1f-9b1c-3b46b3b93563")
BASIS_OK <- c("humanobservation", "materialsample", "preservedspecimen", "occurrence", "machineobservation",
              "livingspecimen", "materialentity", "materialcitation")
LADO <- c("6.25km" = EASE2_BASE_M / 4, "12.5km" = EASE2_BASE_M / 2, "25km" = EASE2_BASE_M)
FOCAL <- c(Polychaeta = "bentos", Bivalvia = "bentos", Echinodermata = "bentos", Anthozoa = "bentos", Porifera = "bentos",
           Cumacea = "bentos", Isopoda = "bentos", Tanaidacea = "bentos", Bryozoa = "bentos", Ascidiacea = "bentos",
           Calanoida = "pelagico", Euphausiacea = "pelagico", Chaetognatha = "pelagico", Hyperiidea = "pelagico",
           Siphonophorae = "pelagico", Thaliacea = "pelagico")
APAR_PEL <- "plankton|bongo|wp-?2|wp-?3|multi ?net|mocness|neuston|\\bcpr\\b|continuous plankton|midwater|mid-water|pelagic trawl|ring net|vertical (haul|tow)|oblique|tucker|ikmt|isaacs|juday|hensen|manta|ichthyoplankton|cone net|purse seine|calvet|calcofi"
APAR_BEN <- "grab|van ?veen|day grab|smith|box ?-?cor|corer|\\bcore\\b|dredge|beam trawl|bottom trawl|otter trawl|sledge|epibenthic|benthic|\\brov\\b|diver|scuba|quadrat|transect|suction|airlift|sigsbee|agassiz|\\btrap|\\bpot\\b|survey trawl|4 seam|bridle|rock|intertidal"
LARVA <- "larv|egg|ovo|zoea|nauplii|naupl|megalop|juvenile stage|postlarv|veliger|ichthyo"

cheio <- function(x) !is.na(x) & as.character(x) != ""
ms <- function(x) as.POSIXct(as.numeric(x) / 1000, origin = "1970-01-01", tz = "UTC")
decimais <- function(verb, num) {
  s <- ifelse(cheio(verb) & grepl("^-?[0-9]+([.,][0-9]+)?$", trimws(verb)), trimws(verb), sub("0+$", "", formatC(abs(num), format = "f", digits = 6)))
  ifelse(grepl("[.,]", s), nchar(sub("^[^.,]*[.,]", "", s)), 0L)
}
tem_flag <- function(fl, nomes) vapply(fl, function(f) any(f %in% nomes), TRUE)
primeira_flag <- function(fl, nomes) vapply(fl, function(f) { k <- f[f %in% nomes]; if (length(k)) k[1] else NA_character_ }, "")
le_data <- function(txt) {                                          # same as 40_curadoria_n1.R
  t <- gsub("\\s+", "", txt); p <- strsplit(ifelse(is.na(t), "", t), "/", fixed = TRUE)
  ini <- vapply(p, function(v) if (length(v)) v[1] else "", ""); fim <- vapply(p, function(v) if (length(v) > 1 && nzchar(v[2])) v[2] else if (length(v)) v[1] else "", "")
  limite <- function(s, final) {
    s <- substr(s, 1, 10); d <- rep(as.Date(NA), length(s))
    a <- grepl("^\\d{4}-\\d{2}-\\d{2}$", s); d[a] <- as.Date(s[a], "%Y-%m-%d")
    b <- grepl("^\\d{4}-\\d{2}$", s); d[b] <- as.Date(paste0(s[b], "-01"), "%Y-%m-%d")
    if (final) d[b] <- as.Date(format(d[b] + 31, "%Y-%m-01"), "%Y-%m-%d") - 1
    c <- grepl("^\\d{4}$", s); d[c] <- as.Date(paste0(s[c], if (final) "-12-31" else "-01-01"), "%Y-%m-%d"); d }
  di <- limite(ini, FALSE); df <- limite(fim, TRUE)
  outro <- is.na(di) & grepl("(18|19|20)\\d{2}", t)
  if (any(outro)) { an <- regmatches(t[outro], gregexpr("(18|19|20)\\d{2}", t[outro]))
    a1 <- vapply(an, function(a) if (length(unique(a)) == 1) a[1] else NA_character_, "")
    di[outro] <- as.Date(paste0(a1, "-01-01"), "%Y-%m-%d"); df[outro] <- as.Date(paste0(a1, "-12-31"), "%Y-%m-%d") }
  ok <- !is.na(di) & !is.na(df) & df >= di
  data.frame(ano = ifelse(ok & format(di, "%Y") == format(df, "%Y"), as.integer(format(di, "%Y")), NA_integer_),
             ini_ms = ifelse(ok, as.numeric(as.POSIXct(di, tz = "UTC")) * 1000, NA_real_), fim_ms = ifelse(ok, as.numeric(as.POSIXct(df, tz = "UTC")) * 1000, NA_real_))
}
completa_datas <- function(x) {                                    # same as 40_curadoria_n1.R
  x$data_relida <- FALSE
  sem <- which(is.na(x$year) & grepl("^\\d{4}$", trimws(x$year_txt)))
  if (length(sem)) {
    a <- as.integer(trimws(x$year_txt[sem])); m <- suppressWarnings(as.integer(x$month_txt[sem])); dd <- suppressWarnings(as.integer(x$day_txt[sem])); m[!m %in% 1:12] <- NA
    ini <- as.Date(sprintf("%04d-%02d-%02d", a, ifelse(is.na(m), 1L, m), ifelse(is.na(m) | is.na(dd), 1L, dd)), "%Y-%m-%d")
    fa <- as.Date(sprintf("%04d-12-31", a), "%Y-%m-%d")
    fm <- as.Date(format(as.Date(sprintf("%04d-%02d-01", a, ifelse(is.na(m), 1L, m)), "%Y-%m-%d") + 31, "%Y-%m-01"), "%Y-%m-%d") - 1
    fim <- ifelse(is.na(m), fa, ifelse(is.na(dd), fm, ini)); okd <- !is.na(ini) & a >= 1700 & a <= 2026
    x$year[sem[okd]] <- a[okd]
    x$date_start[sem[okd]] <- bit64::as.integer64(as.numeric(as.POSIXct(ini[okd], tz = "UTC")) * 1000)
    x$date_end[sem[okd]] <- bit64::as.integer64(as.numeric(as.POSIXct(as.Date(fim[okd], origin = "1970-01-01"), tz = "UTC")) * 1000)
    x$data_relida[sem[okd]] <- TRUE
  }
  txt <- ifelse(cheio(x$eventDate), x$eventDate, x$verbatimEventDate); sem <- which(is.na(x$year) & cheio(txt))
  if (length(sem)) { d <- le_data(txt[sem]); okd <- !is.na(d$ano)
    x$year[sem[okd]] <- d$ano[okd]; x$date_start[sem[okd]] <- bit64::as.integer64(d$ini_ms[okd])
    x$date_end[sem[okd]] <- bit64::as.integer64(d$fim_ms[okd]); x$data_relida[sem[okd]] <- TRUE }
  x
}

# ---- ocean mask at 25 km (seabed in the EBU or water in the EMU) ----
g25 <- ease2_grade("25km"); g6 <- ease2_grade("6.25km")
fundo <- read_parquet(file.path(raiz, "camadas/ebu/ebu_celula_25km.parquet"))$cell_id
emu25 <- rast(file.path(raiz, "camadas/emu/emu_25km.tif")); ag <- which(values(app(emu25, function(v) any(!is.na(v))))[, 1] == 1)
rcA <- rowColFromCell(g25, ag); oceano <- union(fundo, paste("25km", rcA[, 1], rcA[, 2], sep = "_"))

# ---- filters (reason or NA) ----
f41 <- function(x) {
  inval <- is.na(x$lon) | is.na(x$lat) | abs(x$lat) > 90 | abs(x$lon) > 180 | (x$lon == 0 & x$lat == 0) |
           tem_flag(x$flags, c("NO_COORD", "ZERO_COORD", "LAT_OUT_OF_RANGE", "LON_OUT_OF_RANGE"))
  m <- ifelse(inval, "coord_invalida", NA); ok <- is.na(m)
  c6 <- rep(NA_integer_, nrow(x)); c6[ok] <- cellFromXY(g6, project(cbind(x$lon[ok], x$lat[ok]), from = "EPSG:4326", to = "EPSG:6933"))
  rc <- rowColFromCell(g6, c6); r6 <- rc[, 1]; k6 <- rc[, 2]
  c25 <- ifelse(is.na(c6), NA, paste("25km", (r6 - 1) %/% 4 + 1, (k6 - 1) %/% 4 + 1, sep = "_"))
  m[ok & !(c25 %in% oceano)] <- "fora_oceano"
  attr(m, "r6") <- r6; attr(m, "c6") <- k6; m
}
f42 <- function(x) ifelse(x$dataset_id %in% CONJ_CORROMPIDO, "conjunto_corrompido", NA)
f43 <- function(x) { m <- primeira_flag(x$flags, FLAGS_EXCLUI); m <- ifelse(is.na(m), NA, paste0("flag_", tolower(m)))
  m[is.na(m) & !is.na(x$shoredistance) & x$shoredistance < TERRA_M] <- "em_terra"; m }
f44 <- function(x) {
  aus <- x$absence %in% TRUE | tolower(x$occurrenceStatus) %in% "absent"
  idq <- grepl("cf\\.?|aff\\.?|\\?", x$identificationQualifier) | grepl(" cf\\.? | aff\\.? |\\?", x$scientificName)
  b <- gsub("[^a-z]", "", tolower(x$basisOfRecord))
  ifelse(aus, "ausencia", ifelse(idq, "id_incerta", ifelse(b == "fossilspecimen", "fossil", ifelse(b == "nomenclaturalchecklist",
    "basis_recusado", ifelse(cheio(b) & !(b %in% BASIS_OK), "basis_invalido", NA)))))
}
f45 <- function(x) ifelse(is.na(x$aphiaid), "sem_aphia", ifelse(!cheio(x$phylum) & !cheio(x$class) & !cheio(x$order) &
  !cheio(x$family) & !cheio(x$genus), "sem_classificacao", NA))
f46 <- function(x) ifelse(x$dataset_id %in% CONJ_LIC_OK, NA, "licenca_recusada")
f47 <- function(x) ifelse(is.na(x$year), "sem_ano", NA)
f48 <- function(x) {
  erro <- ifelse(!is.na(x$coordUncertainty), x$coordUncertainty, 0.5 * 10^(-x$decimais) * 111320)
  rmax <- ifelse(erro <= LADO["6.25km"] / 4, "6.25km", ifelse(erro <= LADO["12.5km"] / 4, "12.5km", ifelse(erro <= LADO["25km"] / 4, "25km", NA)))
  # level 2 D11: without uncertainty, 2 decimal places before 2000 are limited to 25 km
  lim <- is.na(x$coordUncertainty) & x$decimais == 2 & x$year < 2000 & rmax %in% c("6.25km", "12.5km"); rmax[lim] <- "25km"
  m <- ifelse(is.na(rmax), ifelse(!is.na(x$coordUncertainty), "incerteza_alta", "precisao_baixa"), NA)
  attr(m, "rmax") <- rmax; m
}
f49 <- function(x) {
  p <- pmax(x$depth_max, x$depth_min, x$depth, na.rm = TRUE); tol <- pmax(20, 0.1 * pmax(x$bathymetry, 0, na.rm = TRUE))
  ruim <- (!is.na(x$depth_min) & !is.na(x$depth_max) & x$depth_min > x$depth_max) |
          (!is.na(p) & pmin(x$depth_min, x$depth_max, x$depth, na.rm = TRUE) < 0) |
          (!is.na(p) & !is.na(x$bathymetry) & x$bathymetry > 0 & p > x$bathymetry + tol)
  ifelse(ruim, "profundidade_incoerente", NA)
}
f410 <- function(x) {                                               # same as 40_curadoria_n1.R (corrected version)
  score <- cheio(x$eventID) + !is.na(x$depth_min) + !is.na(x$depth_max) + !is.na(x$coordUncertainty) + cheio(x$samplingProtocol) +
           cheio(x$family) + cheio(x$genus) + cheio(x$lifeStage)
  o <- order(-score, x$record_id); m <- rep(NA_character_, nrow(x))
  k1 <- ifelse(cheio(x$occurrenceID), paste(x$dataset_id, x$occurrenceID), NA)
  d1 <- rep(FALSE, nrow(x)); d1[o] <- !is.na(k1[o]) & duplicated(k1[o], incomparables = NA); m[d1] <- "dup_interna"
  pz <- round(pmax(x$depth_max, x$depth_min, x$depth, na.rm = TRUE)); k2 <- paste(x$aphiaid, round(x$lon, 4), round(x$lat, 4), x$date_start, pz)
  vivo <- which(is.na(m) & !is.na(x$date_start)); t2 <- data.frame(i = vivo, k = k2[vivo], cj = x$dataset_id[vivo], s = score[vivo])
  melhor <- t2 |> group_by(k, cj) |> summarise(s = mean(s), .groups = "drop") |> group_by(k) |> filter(n() > 1) |>
    arrange(desc(s), cj) |> summarise(fica = first(cj), .groups = "drop")
  t2 <- inner_join(t2, melhor, by = "k"); m[t2$i[t2$cj != t2$fica]] <- "dup_entre_conjuntos"; m
}
f411 <- function(x) {                                               # compartment per record (E5.3)
  comp <- rep(NA_character_, nrow(x)); fonte <- rep(NA_character_, nrow(x))
  foc <- x$grupo_amplo %in% names(FOCAL); comp[foc] <- FOCAL[x$grupo_amplo[foc]]; fonte[foc] <- "taxon focal"
  apar <- tolower(paste(x$samplingProtocol, x$aparelho_emof))
  pel <- grepl(APAR_PEL, apar); ben <- grepl(APAR_BEN, apar)
  s <- is.na(comp) & (pel != ben); comp[s] <- ifelse(pel[s], "pelagico", "bentos"); fonte[s] <- "aparelho"
  p <- pmax(x$depth_max, x$depth_min, x$depth, na.rm = TRUE); tol <- pmax(20, 0.1 * pmax(x$bathymetry, 0, na.rm = TRUE))
  temp <- !is.na(p) & !is.na(x$bathymetry) & x$bathymetry > 0
  s <- is.na(comp) & temp; comp[s] <- ifelse(p[s] >= x$bathymetry[s] - tol[s], "bentos", "pelagico"); fonte[s] <- "profundidade"
  larva <- grepl(LARVA, tolower(x$lifeStage))
  tc <- tab$compartimento[match(x$family, tab$family)]
  s <- is.na(comp) & !larva & tc %in% c("bentos", "pelagico"); comp[s] <- tc[s]; fonte[s] <- "atributo do taxon"
  ac <- APAR_CONJ$aparelho_conjunto[match(x$dataset_id, APAR_CONJ$dataset_id)]   # (3b) Gate 3, item 4
  s <- is.na(comp) & ac %in% c("bentos", "pelagico") & !(tc %in% "fora"); comp[s] <- ac[s]; fonte[s] <- "aparelho do conjunto"
  m <-ifelse(tc %in% "fora", "taxon_fora", ifelse(is.na(comp) & larva, "estagio_incerto", ifelse(is.na(comp), "compartimento_incerto", NA)))
  attr(m, "comp") <- comp; attr(m, "fonte") <- fonte; attr(m, "larva") <- larva; m
}
FILTROS <- c("4.1 geografia" = "f41", "4.2 conjuntos" = "f42", "4.3 sinalizadores" = "f43", "4.4 status" = "f44", "4.5 taxon" = "f45",
             "4.6 licenca" = "f46", "4.7 data" = "f47", "4.8 precisao" = "f48", "4.9 profundidade" = "f49", "4.10 duplicatas" = "f410",
             "7.3 compartimento" = "f411")

grupo_eco <- function(x, comp, larva) {
  g <- tab$grupo_ecologico[match(x$family, tab$family)]
  peixe <- x$grupo_amplo == "Peixes"
  def <- is.na(g) | g %in% c("a definir pelo registro", "fora") | (comp == "bentos" & !grepl("infauna|epifauna|meiofauna|demersais", g)) |
         (comp == "pelagico" & grepl("infauna|epifauna|meiofauna|demersais", g))
  g[def & comp == "bentos"] <- ifelse(peixe[def & comp == "bentos"], "peixes demersais",
    ifelse(x$grupo_amplo[def & comp == "bentos"] %in% c("Harpacticoida", "Ostracoda", "Foraminifera"), "meiofauna",
    ifelse(x$grupo_amplo[def & comp == "bentos"] %in% c("Hydrozoa", "Porifera", "Anthozoa", "Bryozoa", "Ascidiacea"), "epifauna séssil", "epifauna móvel e megafauna")))
  g[def & comp == "pelagico"] <- ifelse(peixe[def & comp == "pelagico"], "peixes pelágicos pequenos",
    ifelse(x$grupo_amplo[def & comp == "pelagico"] %in% c("Hydrozoa", "Siphonophorae", "Thaliacea"), "gelatinosos", "zooplâncton holoplanctônico"))
  g[comp == "pelagico" & larva] <- "meroplâncton e ictioplâncton"
  g
}

cascata <- list()
partes <- ds |> distinct(regiao, grupo_amplo) |> collect() |> arrange(regiao, grupo_amplo)
for (i in seq_len(nrow(partes))) {
  rg <- partes$regiao[i]; ga <- partes$grupo_amplo[i]
  x <- ds |> filter(regiao == rg, grupo_amplo == ga) |> collect()
  x$regiao <- rg
  x <- x |> left_join(emof, by = "record_id")
  x <- completa_datas(x)
  x$decimais <- pmin(decimais(x$verbatimLon, x$lon), decimais(x$verbatimLat, x$lat))
  n <- nrow(x); motivo <- rep(NA_character_, n); filtro <- rep(NA_character_, n)
  x$r6 <- NA_integer_; x$c6 <- NA_integer_; x$res_max <- NA_character_; x$compartimento <- NA_character_
  x$fonte_compartimento <- NA_character_; x$larva <- FALSE
  for (nm in names(FILTROS)) {
    vivo <- which(is.na(motivo)); entra <- length(vivo); if (!entra) next
    m <- get(FILTROS[[nm]])(x[vivo, ])
    if (FILTROS[[nm]] == "f41") { x$r6[vivo] <- attr(m, "r6"); x$c6[vivo] <- attr(m, "c6") }
    if (FILTROS[[nm]] == "f48") x$res_max[vivo] <- attr(m, "rmax")
    if (FILTROS[[nm]] == "f411") { x$compartimento[vivo] <- attr(m, "comp"); x$fonte_compartimento[vivo] <- attr(m, "fonte"); x$larva[vivo] <- attr(m, "larva") }
    sai <- !is.na(m); motivo[vivo[sai]] <- m[sai]; filtro[vivo[sai]] <- nm
    tb <- table(m[sai])
    cascata[[length(cascata) + 1]] <- data.frame(regiao = rg, grupo_amplo = ga, filtro = nm, entram = entra, saem = sum(sai),
      restam = entra - sum(sai), pct_saem = round(100 * sum(sai) / entra, 2), motivos = paste(names(tb), tb, sep = "=", collapse = "; "))
  }
  x$motivo <- motivo; x$filtro <- filtro
  r <- x[is.na(motivo), ]
  if (nrow(r)) {
    r$cell_6.25 <- paste("6.25km", r$r6, r$c6, sep = "_")
    r$cell_12.5 <- paste("12.5km", (r$r6 - 1) %/% 2 + 1, (r$c6 - 1) %/% 2 + 1, sep = "_")
    r$cell_25 <- paste("25km", (r$r6 - 1) %/% 4 + 1, (r$c6 - 1) %/% 4 + 1, sep = "_")
    r$grupo_ecologico <- grupo_eco(r, r$compartimento, r$larva)
    di <- ms(r$date_start); dfim <- ms(r$date_end)
    r$mes <- ifelse(!is.na(di) & format(di, "%Y-%m") == format(dfim, "%Y-%m"), as.integer(format(di, "%m")), NA)
    r$dia_conhecido <- !is.na(di) & format(di, "%Y-%m-%d") == format(dfim, "%Y-%m-%d")
    r$nivel_certeza <- ifelse(!is.na(r$coordUncertainty) & r$dia_conhecido, "A", ifelse(!is.na(r$mes), "B", "C"))
    r$evento <- ifelse(cheio(r$eventID), paste(r$dataset_id, r$eventID), paste(r$dataset_id, r$date_start, round(r$lon, 4), round(r$lat, 4)))
    r$unidade_tax <- ifelse(ga %in% c("Euphausiacea", "Chaetognatha"), "genero", "familia")
    r$taxon_unidade <- ifelse(r$unidade_tax == "genero", r$genus, r$family)
    r$so_esforco <- !cheio(r$taxon_unidade)
    r$sem_profundidade <- is.na(r$depth_min) & is.na(r$depth_max) & is.na(r$depth)
    write_dataset(r, file.path(raiz, "intermediario/n2_curado"), partitioning = "regiao", format = "parquet", basename_template = sprintf("p%03d-{i}.parquet", i))
  }
  q <- x[!is.na(motivo) & motivo != "ausencia", ]; if (nrow(q)) write_dataset(q, file.path(raiz, "intermediario/quarentena/n2"), partitioning = "regiao", format = "parquet", basename_template = sprintf("p%03d-{i}.parquet", i))
  a <- x[motivo %in% "ausencia", ]; if (nrow(a)) write_dataset(a, file.path(raiz, "intermediario/ausencias/n2"), partitioning = "regiao", format = "parquet", basename_template = sprintf("p%03d-{i}.parquet", i))
  stopifnot(nrow(r) + nrow(q) + nrow(a) == n)   # each part with its own name (p<i>-...), so as not to overwrite the previous parts of the same region
  if (n > 2e5) cat(format(Sys.time()), "-", rg, "/", ga, ": extraidos", n, "| retidos", nrow(r), "\n")
  rm(x, r, q, a); gc(verbose = FALSE)
}
cas <- bind_rows(cascata)
write.csv(cas, file.path(raiz, "saida/cascata_filtros_n2.csv"), row.names = FALSE)
res <- cas |> group_by(regiao, filtro) |> summarise(entram = sum(entram), saem = sum(saem), .groups = "drop") |>
  mutate(pct = round(100 * saem / entram, 1))
alerta <- res |> filter(pct > 30, !filtro %in% c("4.3 sinalizadores", "4.4 status"))
if (nrow(alerta)) { cat("\n*** PARADA: filtros acima de 30% (por regiao) ***\n"); print(as.data.frame(alerta)) }
print(cas |> filter(filtro == "7.3 compartimento") |> group_by(regiao) |> summarise(restam = sum(restam)) |> as.data.frame())
