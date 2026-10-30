# 86_glorys_campos.R  (Phase 8, stage 3.1: monthly environmental fields at 25 km, for stages 3 and 4)
# Source: GLORYS12V1, monthly climatology 1993-2016, native files (1/12 degree, 50 levels) in camadas/glorys/nativo/.
# Per month, resampled to the 25 km EASE-Grid 2.0 grid (mean of the native cells):
#   thetao_<z>, so_<z>  temperature and salinity at 14 standard depths (native level closest to
#                       0, 10, 25, 50, 100, 200, 300, 500, 1000, 1500, 2000, 3000, 4000, 5000 m);
#   sfundo              bottom salinity (last valid level of the column);
#   vfundo              bottom current speed (last valid level) (m/s);
#   vsup                surface current speed (m/s);
#   gelo                sea ice concentration (0-1);
#   mld, tfundo         mixed layer depth (m) and bottom temperature (already in clim_*_25km.tif; repeated here for convenience);
#   prof_glorys         depth of the last valid level (m), month 1 only.
# Output: camadas/glorys/mes/glorys_<mm>_25km.tif. Parallel: 6 months at a time (memory).
suppressPackageStartupMessages({library(terra); library(parallel)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
dir_m <- file.path(raiz, "camadas/glorys/mes"); dir.create(dir_m, showWarnings = FALSE)
Z_PADRAO <- c(0, 10, 25, 50, 100, 200, 300, 500, 1000, 1500, 2000, 3000, 4000, 5000)

um_mes <- function(m) {
  suppressPackageStartupMessages(library(terra)); source(file.path(raiz, "scripts/funcoes_grade.R"))
  terraOptions(memfrac = 0.12, progress = 0, tempdir = file.path(tempdir(), paste0("m", m)))
  saida <- file.path(dir_m, sprintf("glorys_%02d_25km.tif", m)); if (file.exists(saida)) return(saida)
  a <- file.path(raiz, sprintf("camadas/glorys/nativo/mercatorglorys12v1_gl12_mean_1993_2016_%02d.nc", m))
  lr <- function(v) { r <- rast(a, subds = v); if (is.na(crs(r)) || crs(r) == "") crs(r) <- "EPSG:4326"; r }
  th <- lr("thetao"); so <- lr("so"); uo <- lr("uo"); vo <- lr("vo")
  z <- terra::depth(th)
  iz <- vapply(Z_PADRAO, function(v) which.min(abs(z - v)), 1L)
  nval <- sum(!is.na(th))                                   # valid levels per column (contiguous from the surface)
  nval <- classify(nval, cbind(0, NA))
  fundo <- function(x) selectRange(x, nval)
  vf <- sqrt(fundo(uo)^2 + fundo(vo)^2); vs <- sqrt(uo[[1]]^2 + vo[[1]]^2)
  camadas <- c(th[[iz]], so[[iz]], fundo(so), vf, vs, lr("siconc"), lr("mlotst"), lr("bottomT"))
  nomes <- c(paste0("thetao_", Z_PADRAO), paste0("so_", Z_PADRAO), "sfundo", "vfundo", "vsup", "gelo", "mld", "tfundo")
  if (m == 1) { camadas <- c(camadas, classify(nval, cbind(1:length(z), z))); nomes <- c(nomes, "prof_glorys") }
  g25 <- ease2_grade("25km")
  out <- project(camadas, g25, method = "average"); names(out) <- nomes
  writeRaster(out, saida, overwrite = TRUE); saida
}

cl <- makeCluster(6); clusterExport(cl, c("raiz", "dir_m", "Z_PADRAO", "um_mes"))
t0 <- Sys.time(); r <- parLapply(cl, 1:12, function(m) tryCatch(um_mes(m), error = function(e) paste("ERRO mes", m, conditionMessage(e))))
stopCluster(cl); print(unlist(r)); cat("tempo:", format(Sys.time() - t0), "\n")
x <- rast(file.path(dir_m, "glorys_01_25km.tif")); print(x); print(global(x[[c("thetao_0", "thetao_1000", "sfundo", "vfundo", "gelo", "prof_glorys")]], "range", na.rm = TRUE))
