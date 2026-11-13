# 25_inventario_grupos.R  (PI request of 06/11/2026: revisit level-1 curation; are there taxonomic groups that could enter and did not?)
# One pass over the local OBIS copy: presence records (not dropped, not absences) per phylum x subphylum x class x order x suborder x family,
# with the numbers that decide a candidate group: records from 1993 on, with eventID, with sampling depth >= 200 m and >= 1,000 m, datasets.
# Compartment and ecological group come from entrada/taxon_compartimento.csv (Gate 3; families of the level-2 windows only).
# Output: saida/inventario_grupos_familia.csv (nothing else changes; a group only enters after the PI decides)
suppressPackageStartupMessages({library(arrow); library(dplyr)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); set_cpu_count(6)
ds <- open_dataset(list.files(file.path(raiz, "bruto/occurrence"), full.names = TRUE), format = "parquet", unify_schemas = FALSE)
t0 <- Sys.time()
r <- ds |> filter(!dropped, is.na(absence) | !absence) |>
  mutate(phylum = interpreted$phylum, subphylum = interpreted$subphylum, class = interpreted$class, order = interpreted$order, suborder = interpreted$suborder, family = interpreted$family,
         year = interpreted$date_year, ev = interpreted$eventID, z = coalesce(interpreted$maximumDepthInMeters, interpreted$depth, interpreted$minimumDepthInMeters)) |>
  mutate(rec = year >= 1993L, com_ev = !is.na(ev) & ev != "") |>
  group_by(phylum, subphylum, class, order, suborder, family) |>
  summarise(registros = n(), reg_1993 = sum(rec, na.rm = TRUE), reg_1993_ev = sum(rec & com_ev, na.rm = TRUE), reg_1993_z = sum(rec & !is.na(z), na.rm = TRUE),
            z200 = sum(rec & z >= 200, na.rm = TRUE), z1000 = sum(rec & z >= 1000, na.rm = TRUE), conjuntos = n_distinct(dataset_id), .groups = "drop") |> collect()
cat(format(Sys.time()), "passada em", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min;", nrow(r), "linhas;", sum(r$registros), "registros\n")
write.csv(r, file.path(raiz, "saida/inventario_grupos_familia.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat(format(Sys.time()), "FIM\n")
