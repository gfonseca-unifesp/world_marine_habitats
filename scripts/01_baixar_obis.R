# 01_baixar_obis.R
# Local copy of the full OBIS database (AWS Open Data, s3://obis-open-data/occurrence/).
# Input:  bruto/manifesto_s3.csv (list of files with size and ETag, generated before the download).
# Output: bruto/occurrence/*.parquet and bruto/log_download_parte<k>.csv.
# Can be rerun: files already complete (same size as in the manifest) are skipped.

# Usage: Rscript 01_baixar_obis.R [parte] [n_partes]  (e.g. 1 3, 2 3, 3 3 in parallel; with no arguments downloads everything)
args  <- commandArgs(trailingOnly = TRUE)
parte <- if (length(args) >= 2) as.integer(args[1]) else 1L
n_par <- if (length(args) >= 2) as.integer(args[2]) else 1L

raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
man  <- read.csv(file.path(raiz, "bruto/manifesto_s3.csv"), stringsAsFactors = FALSE)
man  <- man[(seq_len(nrow(man)) - 1L) %% n_par == parte - 1L, ]
man$destino <- file.path(raiz, "bruto", man$key)
man$url     <- paste0("https://obis-open-data.s3.amazonaws.com/", man$key)
completos <- function() file.exists(man$destino) & file.size(man$destino) == man$size

# Small batches (sized to the connection pool) so that no request waits in the queue;
# timeout only for a stalled connection (under 1 kB/s for 5 min).
lote <- 6
for (passada in 1:5) {
  falta <- man[!completos(), ]
  cat(format(Sys.time()), sprintf("- passada %d | arquivos: %d | completos: %d | a baixar: %d (%.1f GB)\n",
      passada, nrow(man), nrow(man) - nrow(falta), nrow(falta), sum(falta$size) / 1024^3))
  if (nrow(falta) == 0) break
  n_falhas <- 0
  for (i in seq(1, nrow(falta), by = lote)) {
    j <- min(i + lote - 1, nrow(falta))
    r <- tryCatch(curl::multi_download(falta$url[i:j], falta$destino[i:j], resume = TRUE,
                    progress = FALSE, multi_timeout = Inf, connecttimeout = 60,
                    low_speed_limit = 1000, low_speed_time = 300),
                  error = function(e) data.frame(success = rep(FALSE, j - i + 1)))
    n_falhas <- n_falhas + sum(!(r$success %in% TRUE))
    if (i %% 300 == 1) cat(format(Sys.time()), sprintf("-   %d de %d\n", j, nrow(falta)))
  }
  cat(format(Sys.time()), sprintf("- fim da passada %d: %d falhas\n", passada, n_falhas))
}

ok <- completos()
log <- data.frame(key = man$key, size_manifesto = man$size,
                  size_local = ifelse(file.exists(man$destino), file.size(man$destino), NA),
                  completo = ok)
write.csv(log, file.path(raiz, sprintf("bruto/log_download_parte%d.csv", parte)), row.names = FALSE)
cat(format(Sys.time()), "- FIM. completos:", sum(ok), "de", nrow(man),
    sprintf("(%.2f GB)\n", sum(man$size[ok]) / 1024^3))
