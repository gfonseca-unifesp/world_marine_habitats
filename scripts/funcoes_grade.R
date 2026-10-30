# funcoes_grade.R
# Definition of the nested equal-area EASE-Grid 2.0 Global grid (EPSG:6933), 25 km family.
# Official NSIDC parameters: 25 km cell = 25025.26 m; 1388 columns x 584 rows.
# Usage: source(file.path(raiz, "scripts/funcoes_grade.R"))

EASE2_BASE_M <- 25025.26
EASE2_NCOL25 <- 1388L
EASE2_NROW25 <- 584L
EASE2_FATOR  <- c("100km" = 1/4, "50km" = 1/2, "25km" = 1, "12.5km" = 2, "6.25km" = 4, "3.125km" = 8)

# Template raster (no values) for a given resolution
ease2_grade <- function(res = "25km") {
  f <- EASE2_FATOR[[res]]
  xmax <- EASE2_NCOL25 * EASE2_BASE_M / 2; ymax <- EASE2_NROW25 * EASE2_BASE_M / 2
  terra::rast(ncols = EASE2_NCOL25 * f, nrows = EASE2_NROW25 * f, xmin = -xmax, xmax = xmax,
              ymin = -ymax, ymax = ymax, crs = "EPSG:6933")
}

# cell_id "<res>_<linha>_<coluna>" from lon/lat (degrees). NA outside the grid.
ease2_cell_id <- function(lon, lat, res = "25km") {
  g <- ease2_grade(res)
  xy <- terra::project(cbind(lon, lat), from = "EPSG:4326", to = "EPSG:6933")
  cel <- terra::cellFromXY(g, xy)
  rc <- terra::rowColFromCell(g, cel)
  ifelse(is.na(cel), NA_character_, paste(res, rc[, 1], rc[, 2], sep = "_"))
}
