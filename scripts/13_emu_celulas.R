# 13_emu_celulas.R  (task 1.3)
# Assigns the EMU (Sayre et al. 2017, v1.0, WOA2013; Esri netCDF emu_v5_20200914.nc, CC BY 4.0)
# to the cells of the EASE-Grid 2.0 grid, at each of the 102 depth levels.
# Method: for each level, reprojects the 1/4 degree grid to the fine 6.25 km grid (nearest neighbor)
# and aggregates to 12.5 / 25 / 50 / 100 km with the dominant EMU (mode) and the purity
# (fraction of subcells with water at that level that have the dominant EMU).
# Output: camadas/emu/emu_<res>.tif and camadas/emu/emu_pureza_<res>.tif (one band per level),
#         camadas/emu/emu_niveis.csv (depth of each band).
suppressPackageStartupMessages({library(ncdf4); library(terra)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
source(file.path(raiz, "scripts/funcoes_grade.R"))
terraOptions(memfrac = 0.6, progress = 0)

nc <- nc_open(file.path(raiz, "camadas/emu/emu_v5_20200914.nc"))
lon <- ncvar_get(nc, "longitude"); lat <- ncvar_get(nc, "latitude"); prof <- ncvar_get(nc, "depth")
cat("lon:", range(lon), "| lat:", range(lat), "| prof:", range(prof), "| n niveis:", length(prof), "\n")
write.csv(data.frame(banda = seq_along(prof), profundidade_m = abs(prof)),
          file.path(raiz, "camadas/emu/emu_niveis.csv"), row.names = FALSE)
dx <- abs(diff(lon[1:2])); dy <- abs(diff(lat[1:2]))
modelo <- rast(ncols = length(lon), nrows = length(lat), xmin = min(lon) - dx / 2, xmax = max(lon) + dx / 2,
               ymin = min(lat) - dy / 2, ymax = max(lat) + dy / 2, crs = "EPSG:4326")

g6 <- ease2_grade("6.25km")
fat <- c("12.5km" = 2, "25km" = 4, "50km" = 8, "100km" = 16)
saida <- setNames(vector("list", length(fat)), names(fat)); pur <- saida
for (k in seq_along(prof)) {
  a <- ncvar_get(nc, "cluster37", start = c(1, 1, k), count = c(-1, -1, 1))   # [lon, lat]
  a[a <= 0] <- NA   # 0 = unclassified (no name in the EMU table); negative = no data
  m <- t(a); if (lat[1] < lat[length(lat)]) m <- m[nrow(m):1, ]                # rows from north to south
  r <- setValues(modelo, as.vector(t(m)))
  r6 <- project(r, g6, method = "near")
  for (res in names(fat)) {
    f <- fat[[res]]
    mo <- aggregate(r6, f, fun = "modal", na.rm = TRUE)
    nf <- aggregate(!is.na(r6), f, fun = "sum")
    ig <- aggregate(r6 == disagg(mo, f), f, fun = "sum", na.rm = TRUE)
    saida[[res]] <- c(saida[[res]], mo); pur[[res]] <- c(pur[[res]], round(100 * ig / nf))
  }
  if (k %% 10 == 0) cat(format(Sys.time()), "- nivel", k, "de", length(prof), "\n")
}
nc_close(nc)
for (res in names(fat)) {
  s <- rast(saida[[res]]); p <- rast(pur[[res]])
  names(s) <- names(p) <- paste0("prof_", abs(prof), "m")
  writeRaster(s, file.path(raiz, sprintf("camadas/emu/emu_%s.tif", res)), overwrite = TRUE, datatype = "INT1U")
  writeRaster(p, file.path(raiz, sprintf("camadas/emu/emu_pureza_%s.tif", res)), overwrite = TRUE, datatype = "INT1U")
  b0 <- which.min(abs(prof))                                   # surface band (0 m)
  sup <- values(s[[b0]])[, 1]; ps <- values(p[[b0]])[, 1]
  cat(res, ": celulas com agua na superficie", sum(!is.na(sup)),
      "| mediana da pureza na superficie (%):", median(ps[!is.na(sup)]), "\n")
}
cat(format(Sys.time()), "- FIM\n")
