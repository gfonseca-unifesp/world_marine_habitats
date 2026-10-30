# 11_ebu_baixar.R  (task 1.2, part 1)
# Downloads the published EBU (Harris et al. 2026; Esri, CC BY 4.0) from the feature service
# EBU_Publication_WFL1, layer 29 (geomorphic_features_with_sescapes_and_EBU), in pages of 2000 features.
# Geometry in EPSG:4326, generalized to 0.002 degree (~200 m), sufficient for cells of 3 km or larger.
# Output: camadas/ebu/paginas/*.geojson and camadas/ebu/ebu.gpkg. Can be rerun (pages already done are skipped).
suppressPackageStartupMessages({library(sf); library(jsonlite)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
dir_pg <- file.path(raiz, "camadas/ebu/paginas"); dir.create(dir_pg, recursive = TRUE, showWarnings = FALSE)
url <- "https://services5.arcgis.com/7weheFjxuNkGGiZi/arcgis/rest/services/EBU_Publication_WFL1/FeatureServer/29/query"
campos <- "OBJECTID,geomorphic_group,Geomorphic,group_cluster,ebu_type,ebu_group,display_flag,MAJORITY_PERCENT,geodesic_area_km2"

ids <- fromJSON(paste0(url, "?where=1%3D1&returnIdsOnly=true&f=json"))$objectIds
ids <- sort(ids); cat(format(Sys.time()), "- feicoes:", length(ids), "\n")
writeLines(c(format(Sys.time(), tz = "UTC", usetz = TRUE), url, paste("feicoes:", length(ids))),
           file.path(raiz, "camadas/ebu/ebu_download_info.txt"))

lotes <- split(ids, ceiling(seq_along(ids) / 2000))
for (k in seq_along(lotes)) {
  arq <- file.path(dir_pg, sprintf("p%04d.geojson", k))
  if (file.exists(arq) && file.size(arq) > 1000) next
  q <- list(where = sprintf("OBJECTID>=%d AND OBJECTID<=%d", min(lotes[[k]]), max(lotes[[k]])),
            outFields = campos, outSR = 4326, maxAllowableOffset = 0.002, geometryPrecision = 5,
            returnGeometry = "true", f = "geojson")
  for (tent in 1:5) {
    ok <- tryCatch({
      curl::curl_download(paste0(url, "?", paste(names(q), sapply(q, function(v) URLencode(as.character(v), TRUE)),
                                                  sep = "=", collapse = "&")), arq, quiet = TRUE)
      TRUE }, error = function(e) FALSE)
    if (ok && file.size(arq) > 1000) break
    Sys.sleep(10 * tent)
  }
  if (k %% 25 == 0) cat(format(Sys.time()), "- pagina", k, "de", length(lotes), "\n")
}

cat(format(Sys.time()), "- juntando paginas (GDAL, modo append)\n")
pg <- list.files(dir_pg, full.names = TRUE)
gpkg <- file.path(raiz, "camadas/ebu/ebu.gpkg"); if (file.exists(gpkg)) file.remove(gpkg)
for (i in seq_along(pg)) {
  sf::gdal_utils("vectortranslate", pg[i], gpkg,
                 options = c(if (i > 1) "-append", "-nln", "ebu", "-nlt", "MULTIPOLYGON", "-f", "GPKG"))
}
ids_lidos <- st_read(gpkg, query = "SELECT OBJECTID FROM ebu", quiet = TRUE)$OBJECTID
cat("feicoes lidas:", length(ids_lidos), "de", length(ids), "| duplicadas:", sum(duplicated(ids_lidos)), "\n")
cat(format(Sys.time()), "- FIM\n")
