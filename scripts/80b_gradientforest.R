# 80b_gradientforest.R  (installs extendedForest and gradientForest from the official R-Forge mirror on GitHub)
# R-Forge was down on 17/10/2026; github.com/r-forge/gradientforest is the mirror of the SVN repository.
# Requires Rtools45 (C:/rtools45), installed on 19/10/2026 with the user's authorization.
Sys.setenv(PATH = paste("C:/rtools45/usr/bin", "C:/rtools45/x86_64-w64-mingw32.static.posix/bin", Sys.getenv("PATH"), sep = ";"))
stopifnot(pkgbuild::has_build_tools(debug = FALSE))
zip <- file.path(tempdir(), "gf.zip")
curl::curl_download("https://codeload.github.com/r-forge/gradientforest/zip/refs/heads/master", zip)
dir <- file.path(tempdir(), "gf"); unzip(zip, exdir = dir)
pk <- list.dirs(dir, recursive = TRUE); pk <- pk[basename(pk) %in% c("extendedForest", "gradientForest") & file.exists(file.path(pk, "DESCRIPTION"))]
print(pk)
for (p in c("extendedForest", "gradientForest")) install.packages(pk[basename(pk) == p][1], repos = NULL, type = "source")
suppressPackageStartupMessages({library(gradientForest); library(vegan)})
cat("gradientForest", as.character(packageVersion("gradientForest")), "| extendedForest", as.character(packageVersion("extendedForest")), "\n")
data(dune); data(dune.env); X <- data.frame(A1 = dune.env$A1, Moisture = as.numeric(dune.env$Moisture))
gf <- gradientForest(cbind(X, dune), predictor.vars = names(X), response.vars = names(dune), ntree = 200, trace = FALSE)
cat("teste ok: especies com R2 > 0 =", length(gf$result), "| importancia:", round(importance(gf), 3), "\n")
