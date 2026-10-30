# 80_pacotes.R  (Phase 8: installation and testing of the analysis packages)
# vegan (PERMANOVA, adonis2), indicspecies (IndVal.g), gdm (CRAN) and gradientForest + extendedForest (R-Forge).
lib <- .libPaths()[1]
cran <- c("vegan", "indicspecies", "gdm", "permute", "future.apply")
falta <- cran[!cran %in% rownames(installed.packages())]
if (length(falta)) install.packages(falta, repos = "https://cloud.r-project.org", Ncpus = 8)
if (!all(c("extendedForest", "gradientForest") %in% rownames(installed.packages()))) {
  r <- try(install.packages(c("extendedForest", "gradientForest"), repos = c("https://R-Forge.R-project.org", "https://cloud.r-project.org"), Ncpus = 8))
  if (!"gradientForest" %in% rownames(installed.packages()))   # no binary on R-Forge: source (requires Rtools)
    install.packages(c("extendedForest", "gradientForest"), repos = "https://R-Forge.R-project.org", type = "source")
}
for (p in c(cran, "extendedForest", "gradientForest")) cat(sprintf("%-15s %s\n", p,
  if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "NAO INSTALADO"))

# minimal test of each one
suppressPackageStartupMessages({library(vegan); library(indicspecies); library(gdm); library(gradientForest)})
data(dune); data(dune.env)
print(adonis2(dune ~ Management, data = dune.env, permutations = 199, parallel = 4))
iv <- multipatt(dune, dune.env$Management, func = "IndVal.g", control = how(nperm = 199)); cat("IndVal.g ok:", nrow(iv$sign), "\n")
X <- data.frame(dune.env[, c("A1", "Moisture")]); X$Moisture <- as.numeric(X$Moisture)
gf <- gradientForest(cbind(X, dune), predictor.vars = names(X), response.vars = names(dune), ntree = 100, trace = FALSE)
cat("gradientForest ok: especies com R2>0 =", length(gf$result), "\n")
