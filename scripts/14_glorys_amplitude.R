# 14_glorys_amplitude.R  (task 1.4)
# Seasonal amplitude of the environment, used only to score the level 2 windows (E4, "dynamic potential").
# Source: GLORYS12V1, monthly climatology 1993-2016 (GLOBAL_MULTIYEAR_PHY_001_030,
#         cmems_mod_glo_phy_my_0.083deg-climatology_P1M-m_202311), native files (one netCDF per month),
#         downloaded from the Copernicus Marine S3 to camadas/glorys/nativo/.
#         (The subsetting API was discarded: ~30 h for the global fields; the native files take ~30 min.)
# Variables: temperature and salinity at the first level (0.49 m), mixed layer depth (mlotst)
#            and bottom temperature (bottomT). Bottom salinity and oxygen are left out.
# Output: camadas/glorys/clim_<var>_25km.tif (12 bands) and
#         camadas/glorys/amplitude_<res>.parquet (25 and 100 km): annual mean and amplitude (monthly max - min).
suppressPackageStartupMessages({library(terra); library(arrow)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
source(file.path(raiz, "scripts/funcoes_grade.R"))
dir_g <- file.path(raiz, "camadas/glorys"); terraOptions(memfrac = 0.6, progress = 0)
arqs <- sprintf("%s/nativo/mercatorglorys12v1_gl12_mean_1993_2016_%02d.nc", dir_g, 1:12)
stopifnot(all(file.exists(arqs)))

vars <- c(sst = "thetao", sss = "so", mld = "mlotst", tfundo = "bottomT")
g25 <- ease2_grade("25km")
amp25 <- NULL
for (nm in names(vars)) {
  meses <- lapply(arqs, function(a) { r <- rast(a, subds = vars[[nm]]); r[[1]] })   # 1st layer = surface (3D) or the only layer (2D)
  r <- rast(meses); names(r) <- month.abb
  if (is.na(crs(r)) || crs(r) == "") crs(r) <- "EPSG:4326"
  rp <- project(r, g25, method = "average")
  writeRaster(rp, file.path(dir_g, sprintf("clim_%s_25km.tif", nm)), overwrite = TRUE)
  d <- data.frame(cel = 1:ncell(g25), media = values(mean(rp))[, 1], amplitude = values(max(rp) - min(rp))[, 1])
  names(d)[2:3] <- paste0(nm, c("_media", "_amp"))
  amp25 <- if (is.null(amp25)) d else merge(amp25, d, by = "cel")
  cat(format(Sys.time()), "-", nm, ": camadas lidas", nlyr(r), "| unidade:", units(r)[1],
      "| media global da amplitude:", round(mean(d[[3]], na.rm = TRUE), 2), "\n")
}
amp25 <- amp25[rowSums(!is.na(amp25[, -1])) > 0, ]
rc <- rowColFromCell(g25, amp25$cel); amp25$cell_id <- paste("25km", rc[, 1], rc[, 2], sep = "_")
amp25$cell_id_100 <- paste("100km", (rc[, 1] - 1) %/% 4 + 1, (rc[, 2] - 1) %/% 4 + 1, sep = "_")
write_parquet(amp25[, -1], file.path(dir_g, "amplitude_25km.parquet"))
col <- grep("_(media|amp)$", names(amp25), value = TRUE)
amp100 <- aggregate(amp25[, col], by = list(cell_id = amp25$cell_id_100), FUN = mean, na.rm = TRUE)
write_parquet(amp100, file.path(dir_g, "amplitude_100km.parquet"))
cat(format(Sys.time()), "- FIM. celulas 25 km:", nrow(amp25), "| 100 km:", nrow(amp100), "\n")
print(summary(amp100[, grep("_amp$", names(amp100))]))
