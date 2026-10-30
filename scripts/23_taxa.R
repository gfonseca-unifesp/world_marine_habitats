# 23_taxa.R  (task 2.3)
# Table of level 1 taxa from the WoRMS classification that accompanies OBIS.
# Output: saida/taxa.csv (one AphiaID per row, with the number of records)
suppressPackageStartupMessages({library(arrow); library(dplyr)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
ds <- open_dataset(file.path(raiz, "intermediario/nivel1_extracao"))
taxa <- ds |> filter(!is.na(aphiaid)) |>
  group_by(aphiaid, scientificName, taxonRank, phylum, subphylum, class, order, suborder, family, genus, species, grupo) |>
  summarise(registros = n(), .groups = "drop") |> collect() |> arrange(grupo, family, genus, species)
dup <- taxa |> count(aphiaid) |> filter(n > 1)
write.csv(taxa, file.path(raiz, "saida/taxa.csv"), row.names = FALSE)
sem <- ds |> filter(is.na(aphiaid)) |> count() |> collect()
cat("taxons (AphiaID):", n_distinct(taxa$aphiaid), "| AphiaID com mais de uma classificacao:", nrow(dup),
    "| registros sem AphiaID:", sem$n, "\n")
print(taxa |> group_by(grupo) |> summarise(taxons = n_distinct(aphiaid), familias = n_distinct(family, na.rm = TRUE),
                                          pct_reg_familia = round(100 * sum(registros[!is.na(family)]) / sum(registros), 1)) |>
        as.data.frame())
