# 10_grades.R  (task 1.1)
# Documents and tests the nested EASE-Grid 2.0 grid defined in funcoes_grade.R.
# Output: camadas/grades/grades_parametros.csv. The ocean mask comes from the EBU and EMU layers (1.2 and 1.3).
suppressPackageStartupMessages(library(terra))
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
source(file.path(raiz, "scripts/funcoes_grade.R"))
dir.create(file.path(raiz, "camadas/grades"), recursive = TRUE, showWarnings = FALSE)

par <- do.call(rbind, lapply(names(EASE2_FATOR), function(nm) {
  r <- ease2_grade(nm)
  data.frame(resolucao = nm, celula_m = res(r)[1], ncol = ncol(r), nrow = nrow(r),
             n_celulas = ncell(r), area_km2 = round((res(r)[1] / 1000)^2, 2))
}))
ymax <- EASE2_NROW25 * EASE2_BASE_M / 2
lim <- project(cbind(0, ymax - 1), from = "EPSG:6933", to = "EPSG:4326")[2]
par$lat_limite_graus <- round(lim, 2)
write.csv(par, file.path(raiz, "camadas/grades/grades_parametros.csv"), row.names = FALSE)
print(par)

# Nesting test: the center of each fine cell falls in the expected coarse cell
set.seed(1)
g6 <- ease2_grade("6.25km"); g25 <- ease2_grade("25km"); g100 <- ease2_grade("100km")
am <- sample(ncell(g6), 2e5); xy <- xyFromCell(g6, am); rc6 <- rowColFromCell(g6, am)
rc25 <- rowColFromCell(g25, cellFromXY(g25, xy)); rc100 <- rowColFromCell(g100, cellFromXY(g100, xy))
ok25  <- all((rc6 - 1) %/% 4 + 1 == rc25)
ok100 <- all((rc6 - 1) %/% 16 + 1 == rc100)
cat("\naninhamento 6,25 -> 25 km:", ok25, "| 6,25 -> 100 km:", ok100, "\n")
# Test of the cell_id function
cat("Foz do Amazonas (-48.5, 2.5):", ease2_cell_id(-48.5, 2.5, "25km"), "\n")
