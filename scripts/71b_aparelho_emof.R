# 71b_aparelho_emof.R  (tasks 0.3 and 7.3)
# Sampling gear declared in the measurements extension (eMoF: measurementType with "sampling instrument/device/
# gear/protocol", "gear", "net type", "instrument"), only for level 2 records WITHOUT samplingProtocol.
# Reads each dataset in batches of up to 20 thousand records, to keep memory low (the previous version reached 56 GB).
# The gear usually belongs to the event: within a batch, one record per eventID is enough.
# Output: intermediario/nivel2_aparelho_emof.parquet (record_id, aparelho_emof)
suppressPackageStartupMessages({library(arrow); library(dplyr)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
x <- open_dataset(file.path(raiz, "intermediario/nivel2_extracao")) |>
  filter(is.na(samplingProtocol) | samplingProtocol == "") |> select(record_id, dataset_id, eventID) |> collect()
EMOF <- "http://rs.iobis.org/obis/terms/ExtendedMeasurementOrFact"
apar_emof <- function(t) {
  if (is.null(t) || !NROW(t) || is.null(t$source)) return(NA_character_)
  tp <- tolower(t$source$measurementType); vl <- t$source$measurementValue
  k <- !is.na(tp) & grepl("sampling (instrument|device|gear|protocol)|\\bgear\\b|net type|instrument", tp)
  if (any(k)) paste(unique(trimws(vl[k])), collapse = " | ") else NA_character_
}
# one representative record per event (or the record itself, without eventID)
x$chave <- ifelse(!is.na(x$eventID) & x$eventID != "", paste(x$dataset_id, x$eventID), x$record_id)
rep_ <- x[!duplicated(x$chave), c("record_id", "dataset_id", "chave")]
cat("registros sem samplingProtocol:", nrow(x), "| representantes (eventos):", nrow(rep_), "\n")
out <- list(); LOTE <- 20000
for (cj in unique(rep_$dataset_id)) {
  f <- file.path(raiz, "bruto/occurrence", paste0(cj, ".parquet"))
  ids <- rep_$record_id[rep_$dataset_id == cj]
  for (i in seq(1, length(ids), by = LOTE)) {
    lote <- ids[i:min(i + LOTE - 1, length(ids))]
    e <- tryCatch(open_dataset(f) |> filter(`_id` %in% lote) |> select(`_id`, extensions) |> collect(), error = function(err) NULL)
    if (is.null(e) || !nrow(e) || !is.data.frame(e$extensions) || !(EMOF %in% names(e$extensions))) { rm(e); next }
    a <- vapply(e$extensions[[EMOF]], apar_emof, "")
    if (any(!is.na(a))) out[[length(out) + 1]] <- data.frame(record_id = e$`_id`[!is.na(a)], aparelho_emof = a[!is.na(a)])
    rm(e); gc(verbose = FALSE)
  }
}
r <- bind_rows(out)
r$chave <- rep_$chave[match(r$record_id, rep_$record_id)]
res <- x |> select(record_id, chave) |> inner_join(r |> select(chave, aparelho_emof), by = "chave") |> select(record_id, aparelho_emof)
write_parquet(res, file.path(raiz, "intermediario/nivel2_aparelho_emof.parquet"))
cat("registros com aparelho na eMoF:", nrow(res), "\n")
print(head(sort(table(res$aparelho_emof), decreasing = TRUE), 20))
