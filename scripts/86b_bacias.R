# 86b_bacias.R  (Phase 8, stage 3, Gate 7 item 2: ocean basin as a biogeographic covariate)
# With no official basin map downloaded, basins are delimited by rules and by the connectivity of the GLORYS12 ocean mask
# at 25 km (documented approximation):
#   Arctic: lat >= 66.5; Southern: lat <= -45; Mediterranean and Black: lon -6 to 42, lat 30 to 47;
#   Indian: lon 20 to 100 (lat -45 to 31) and lon 100 to 147 south of 10 S;
#   the remainder is split into connected blocks (explicit barrier at the Isthmus of Panama, because the 25 km mean opens water there;
#   the Pacific has one seed on each side of the dateline): Atlantic and Pacific (seeds 30 W 0 and 150 W 0);
#   Atlantic and Pacific split at the equator (N and S). Leftover isolated patches (closed seas): "outro".
# Output: camadas/bacias_25km.tif (codes) and camadas/bacias_legenda.csv
suppressPackageStartupMessages(library(terra))
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); source(file.path(raiz, "scripts/funcoes_grade.R"))
g <- ease2_grade("25km")
agua <- !is.na(rast(file.path(raiz, "camadas/glorys/mes/glorys_01_25km.tif"))[["thetao_0"]])
ll <- project(xyFromCell(g, 1:ncell(g)), from = "EPSG:6933", to = "EPSG:4326"); lon <- ll[, 1]; lat <- ll[, 2]
w <- values(agua)[, 1] %in% TRUE
b <- rep(NA_integer_, ncell(g))
b[w & lat >= 66.5] <- 1L                                                   # Arctic
b[w & lat <= -45] <- 2L                                                    # Southern
b[w & is.na(b) & lon >= -6 & lon <= 42 & lat >= 30 & lat <= 47] <- 3L     # Mediterranean and Black
b[w & is.na(b) & ((lon >= 20 & lon < 100 & lat < 31) | (lon >= 100 & lon <= 147 & lat < -10))] <- 4L   # Indian
istmo <- lon >= -84 & lon <= -76.5 & lat >= 7 & lat <= 10.5               # at 25 km the mean opens water across Panama
r <- rast(g); r[] <- ifelse(w & is.na(b) & !istmo, 1L, NA); p <- patches(r, directions = 4)
semente <- function(lo, la) p[cellFromXY(g, project(cbind(lo, la), from = "EPSG:4326", to = "EPSG:6933"))][1, 1]
pa <- semente(-30, 0); pp <- c(semente(-150, 0), semente(160, 0)); pv <- values(p)[, 1]   # Pacific: two sides of the dateline
stopifnot(!pa %in% pp)
b[!is.na(pv) & pv == pa] <- ifelse(lat[!is.na(pv) & pv == pa] >= 0, 5L, 6L)   # Atlantic N / S
b[!is.na(pv) & pv %in% pp] <- ifelse(lat[!is.na(pv) & pv %in% pp] >= 0, 7L, 8L)   # Pacific N / S
b[w & is.na(b) & istmo] <- ifelse(lon[w & is.na(b) & istmo] < -80 & lat[w & is.na(b) & istmo] < 9, 7L, 5L)   # coasts of the isthmus
b[w & is.na(b)] <- 9L
leg <- data.frame(codigo = 1:9, bacia = c("Ártico", "Austral", "Mediterrâneo e Negro", "Índico", "Atlântico Norte", "Atlântico Sul",
                                          "Pacífico Norte", "Pacífico Sul", "outro (mares fechados)"))
out <- rast(g); out[] <- b; names(out) <- "bacia"
writeRaster(out, file.path(raiz, "camadas/bacias_25km.tif"), overwrite = TRUE, datatype = "INT1U")
write.csv(leg, file.path(raiz, "camadas/bacias_legenda.csv"), row.names = FALSE)
leg$celulas <- as.integer(table(factor(b, levels = 1:9))); print(leg)
png(file.path(raiz, "saida/mapa_bacias.png"), 2000, 1000, res = 150, bg = "white"); par(mar = c(1, 1, 2, 1))
plot(out, type = "classes", col = c("#8ecae6", "#bdbdbd", "#e9c46a", "#f4a261", "#2a78d6", "#7fb3e8", "#1baf7a", "#8fd6b8", "#d03b3b"),
     levels = leg$bacia, axes = FALSE, box = FALSE, main = "Bacias oceânicas (aproximação por regras)", plg = list(cex = 0.75)); dev.off()
