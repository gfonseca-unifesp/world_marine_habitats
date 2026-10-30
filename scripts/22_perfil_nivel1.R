# 22_perfil_nivel1.R  (task 2.2)
# Profile of the level 1 extraction, by focal group: field completeness, dates, coordinate precision,
# flags, basisOfRecord, licenses on the record and dataset dominance.
# Output: saida/perfil_nivel1.csv (one row per grupo), saida/perfil_flags.csv, saida/perfil_basis.csv,
#        saida/perfil_decimais.csv, saida/perfil_conjuntos_top.csv
suppressPackageStartupMessages({library(arrow); library(dplyr)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
ds <- open_dataset(file.path(raiz, "intermediario/nivel1_extracao"))
cheio <- function(x) !is.na(x) & as.character(x) != ""
ms <- function(x) as.POSIXct(as.numeric(x) / 1000, origin = "1970-01-01", tz = "UTC")
# decimal places of the coordinate: from the original text when available; otherwise from the number rounded to 6 places
decimais <- function(verb, num) {
  s <- ifelse(cheio(verb) & grepl("^-?[0-9]+([.,][0-9]+)?$", trimws(verb)), trimws(verb),
              sub("0+$", "", formatC(abs(num), format = "f", digits = 6)))
  ifelse(grepl("[.,]", s), nchar(sub("^[^.,]*[.,]", "", s)), 0L)
}

grupos <- ds |> distinct(grupo) |> collect() |> pull(grupo) |> sort()
perfil <- list(); flg <- list(); bas <- list(); dec <- list(); top <- list()
for (g in grupos) {
  x <- ds |> filter(grupo == g) |>
    select(dataset_id, lon, lat, verbatimLon, verbatimLat, date_start, date_end, year, depth_min, depth_max,
           depth, bathymetry, eventID, family, genus, samplingProtocol, lifeStage, coordUncertainty,
           basisOfRecord, license, flags, dropped, absence) |> collect()
  di <- ms(x$date_start); df <- ms(x$date_end)
  mes <- !is.na(di) & format(di, "%Y-%m") == format(df, "%Y-%m")
  dia <- !is.na(di) & format(di, "%Y-%m-%d") == format(df, "%Y-%m-%d")
  d <- pmin(decimais(x$verbatimLon, x$lon), decimais(x$verbatimLat, x$lat))
  pres <- !x$absence %in% TRUE
  cj <- sort(table(x$dataset_id), decreasing = TRUE)
  perfil[[g]] <- data.frame(
    grupo = g, registros = nrow(x), presencas = sum(pres), conjuntos = length(cj),
    pct_dropped = 100 * mean(x$dropped %in% TRUE), pct_ausencia = 100 * mean(!pres),
    pct_ano = 100 * mean(!is.na(x$year)), pct_desde1993 = 100 * mean(x$year >= 1993, na.rm = TRUE),
    pct_mes = 100 * mean(mes), pct_dia = 100 * mean(dia),
    pct_prof = 100 * mean(!is.na(x$depth_min) | !is.na(x$depth_max) | !is.na(x$depth)),
    pct_prof_minmax = 100 * mean(!is.na(x$depth_min) & !is.na(x$depth_max)),
    pct_batimetria = 100 * mean(!is.na(x$bathymetry)), pct_incerteza = 100 * mean(!is.na(x$coordUncertainty)),
    pct_eventID = 100 * mean(cheio(x$eventID)), pct_familia = 100 * mean(cheio(x$family)),
    pct_genero = 100 * mean(cheio(x$genus)), pct_aparelho = 100 * mean(cheio(x$samplingProtocol)),
    pct_estagio = 100 * mean(cheio(x$lifeStage)), pct_licenca_registro = 100 * mean(cheio(x$license)),
    pct_dec_ge2 = 100 * mean(d >= 2), pct_dec_ge3 = 100 * mean(d >= 3),
    pct_maior_conjunto = 100 * as.numeric(cj[1]) / nrow(x), maior_conjunto = names(cj)[1])
  fl <- unlist(x$flags); fl <- fl[!is.na(fl)]
  flg[[g]] <- if (length(fl)) data.frame(grupo = g, flag = names(table(fl)), registros = as.integer(table(fl))) else NULL
  bas[[g]] <- as.data.frame(table(basisOfRecord = ifelse(cheio(x$basisOfRecord), x$basisOfRecord, "(vazio)"))) |> mutate(grupo = g)
  dec[[g]] <- as.data.frame(table(decimais = pmin(d, 6))) |> mutate(grupo = g)
  top[[g]] <- data.frame(grupo = g, dataset_id = names(cj)[1:min(5, length(cj))],
                         registros = as.integer(cj[1:min(5, length(cj))]))
  cat(format(Sys.time()), "-", g, "\n"); rm(x); gc(verbose = FALSE)
}
num <- function(df) { df[] <- lapply(df, function(v) if (is.numeric(v)) round(v, 1) else v); df }
write.csv(num(bind_rows(perfil)), file.path(raiz, "saida/perfil_nivel1.csv"), row.names = FALSE)
write.csv(bind_rows(flg), file.path(raiz, "saida/perfil_flags.csv"), row.names = FALSE)
write.csv(bind_rows(bas), file.path(raiz, "saida/perfil_basis.csv"), row.names = FALSE)
write.csv(bind_rows(dec), file.path(raiz, "saida/perfil_decimais.csv"), row.names = FALSE)
write.csv(bind_rows(top), file.path(raiz, "saida/perfil_conjuntos_top.csv"), row.names = FALSE)
print(num(bind_rows(perfil))[, c("grupo", "registros", "pct_mes", "pct_prof", "pct_aparelho", "pct_dec_ge2", "pct_maior_conjunto")])
