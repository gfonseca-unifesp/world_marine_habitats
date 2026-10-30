# 12_ebu_celulas.R  (task 1.2, part 2)
# Assigns the EBU to the cells of the EASE-Grid 2.0 grid, at the three hierarchical levels:
#   geomorphic group (geomorphic_group, 6), benthic region (group_cluster, 57) and unit (ebu_type).
# Method: rasterizes the polygons on the fine 3.125 km grid (cell center) and aggregates to
# 6.25 / 12.5 / 25 / 50 / 100 km with the dominant class (mode) and the purity
# (fraction of seafloor subcells that have the dominant class). frac_fundo = fraction of the cell with EBU.
# Output: camadas/ebu/ebu_classes.csv (codes) and camadas/ebu/ebu_celula_<res>.parquet
suppressPackageStartupMessages({library(sf); library(terra); library(arrow)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
source(file.path(raiz, "scripts/funcoes_grade.R"))
terraOptions(memfrac = 0.6, progress = 0)

ebu <- st_read(file.path(raiz, "camadas/ebu/ebu.gpkg"), quiet = TRUE)
ebu <- ebu[!st_is_empty(ebu), ]
niveis <- c(grupo = "geomorphic_group", regiao = "group_cluster", unidade = "ebu_type")
classes <- do.call(rbind, lapply(names(niveis), function(n) {
  v <- sort(unique(ebu[[niveis[[n]]]])); data.frame(nivel = n, codigo = seq_along(v), classe = v) }))
write.csv(classes, file.path(raiz, "camadas/ebu/ebu_classes.csv"), row.names = FALSE)
print(table(classes$nivel))

v <- vect(st_transform(ebu, 6933))
g_fina <- ease2_grade("3.125km")
for (n in names(niveis)) {
  cod <- match(v[[niveis[[n]]]][, 1], classes$classe[classes$nivel == n])
  v$cod <- cod
  r <- rasterize(v, g_fina, field = "cod", fun = "min", touches = FALSE)
  writeRaster(r, file.path(raiz, sprintf("camadas/ebu/ebu_%s_3.125km.tif", n)), overwrite = TRUE, datatype = "INT2U")
  cat(format(Sys.time()), "- rasterizado:", n, "\n")
}

fat <- c("6.25km" = 2, "12.5km" = 4, "25km" = 8, "50km" = 16, "100km" = 32)
for (res in names(fat)) {
  f <- fat[[res]]; out <- NULL
  for (n in names(niveis)) {
    r  <- rast(file.path(raiz, sprintf("camadas/ebu/ebu_%s_3.125km.tif", n)))
    m  <- aggregate(r, f, fun = "modal", na.rm = TRUE)
    nf <- aggregate(!is.na(r), f, fun = "sum")
    ig <- aggregate(r == disagg(m, f), f, fun = "sum", na.rm = TRUE)
    pur <- ig / nf
    d <- data.frame(cel = 1:ncell(m), cod = values(m)[, 1], pureza = values(pur)[, 1],
                    frac_fundo = values(nf)[, 1] / f^2)
    d <- d[!is.na(d$cod), ]
    d$classe <- classes$classe[classes$nivel == n][d$cod]
    names(d)[names(d) %in% c("cod", "pureza", "classe")] <- paste0(c("cod_", "pureza_", "ebu_"), n)
    out <- if (is.null(out)) d else merge(out, d[, c("cel", paste0(c("ebu_", "pureza_"), n))], by = "cel", all = TRUE)
  }
  rc <- rowColFromCell(ease2_grade(res), out$cel)
  out$cell_id <- paste(res, rc[, 1], rc[, 2], sep = "_")
  out <- out[, c("cell_id", "frac_fundo", "ebu_grupo", "pureza_grupo", "ebu_regiao", "pureza_regiao",
                 "ebu_unidade", "pureza_unidade")]
  write_parquet(out, file.path(raiz, sprintf("camadas/ebu/ebu_celula_%s.parquet", res)))
  cat(format(Sys.time()), "-", res, ": celulas com fundo", nrow(out),
      "| mediana pureza grupo/regiao/unidade:",
      round(median(out$pureza_grupo, na.rm = TRUE), 2), round(median(out$pureza_regiao, na.rm = TRUE), 2),
      round(median(out$pureza_unidade, na.rm = TRUE), 2), "\n")
}
