# 00_ambiente.R  (task 0.2)
# Records the versions of R and of the packages used. Output: saida/ambiente.csv
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
pk <- c("arrow", "dplyr", "terra", "sf", "ncdf4", "jsonlite", "curl", "CopernicusMarine", "stars",
        "base64enc", "blockCV", "CAST", "worrms", "rfishbase")
v <- sapply(pk, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "nao instalado")
amb <- data.frame(item = c("R", "sistema", "GDAL (terra)", pk),
                  versao = c(R.version.string, paste(Sys.info()[c("sysname", "release")], collapse = " "),
                             terra::gdal(), v))
write.csv(amb, file.path(raiz, "saida/ambiente.csv"), row.names = FALSE)
print(amb)
