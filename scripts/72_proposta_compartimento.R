# 72_proposta_compartimento.R  (task 7.2; STOP G5: the table is only applied after user approval)
# Compartment and ecological group proposal by family, for level 2.
# Sources (never from memory for the taxon attribute):
#   invertebrates -> WoRMS functional group (wm_attr_data, "Functional group", adult or unstaged only), for up to the 5 most
#                    frequent species of each family in the data of the 8 regions (raw cache in bruto/api_worms/);
#   fishes        -> FishBase DemersPelag (rfishbase::species), for all fish species in the data.
# The ecological group (E5.5) follows a convention by higher taxon, flagged in the fonte_grupo column, for review.
# Families with no attribute or with a mixed attribute are left "ambiguo": the compartment comes from the record (gear, depth).
# Output: entrada/taxon_compartimento_proposta.csv and saida/proposta_compartimento_resumo.csv
suppressPackageStartupMessages({library(arrow); library(dplyr); library(jsonlite)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
ds <- open_dataset(file.path(raiz, "intermediario/nivel2_extracao"))
tx <- ds |> filter(!(absence %in% TRUE), !is.na(aphiaid)) |>
  count(grupo_amplo, phylum, class, order, suborder, family, genus, species, aphiaid, taxonRank) |> collect()
fam <- tx |> filter(!is.na(family), family != "") |>
  group_by(grupo_amplo, phylum, class, order, family) |> summarise(registros = sum(n), especies = n_distinct(species, na.rm = TRUE), .groups = "drop")
cat("familias nas 8 regioes:", nrow(fam), "\n")

# ---- WoRMS: functional group of the most frequent species of each family (invertebrates) ----
cache <- file.path(raiz, "bruto/api_worms"); dir.create(cache, showWarnings = FALSE)
fg_de <- function(id) {
  arq <- file.path(cache, paste0(id, ".rds"))
  if (!file.exists(arq)) {
    r <- tryCatch(worrms::wm_attr_data(id = id, include_inherited = TRUE), error = function(e) NULL)
    saveRDS(r, arq)
  }
  r <- readRDS(arq); if (is.null(r) || !nrow(r)) return(NA_character_)
  fg <- r[r$measurementType == "Functional group", ]
  if (!nrow(fg)) return(NA_character_)
  # life stage reported as a child of the attribute; the adult (or unstaged) counts, as in E5.3
  estagio <- vapply(fg$children, function(ch) {
    if (is.null(ch) || !NROW(ch)) return(NA_character_)
    s <- ch$measurementValue[ch$measurementType == "Life stage"]; if (length(s)) paste(tolower(s), collapse = ",") else NA_character_ }, "")
  adulto <- is.na(estagio) | grepl("adult", estagio)
  v <- fg$measurementValue[adulto]
  if (!length(v)) return(NA_character_)
  paste(sort(unique(tolower(v))), collapse = "; ")
}
peixe <- tx$grupo_amplo == "Peixes"
amostra <- tx |> filter(!peixe, taxonRank == "Species", !is.na(family)) |> group_by(family) |> slice_max(n, n = 5, with_ties = FALSE) |> ungroup()
cat("consultas ao WoRMS:", nrow(amostra), "\n")
amostra$fg <- vapply(amostra$aphiaid, fg_de, "")
classif <- function(v) {
  v <- tolower(v); b <- grepl("benthos", v); p <- grepl("plankton", v); n <- grepl("nekton", v)
  ifelse(is.na(v), NA, ifelse(b & !p & !n, "bentos", ifelse((p | n) & !b, "pelagico", "misto")))
}
amostra$comp <- classif(amostra$fg)
fg_fam <- amostra |> group_by(family) |>
  summarise(n_consultadas = n(), n_com_atributo = sum(!is.na(comp)), bentos = sum(comp == "bentos", na.rm = TRUE),
            pelagico = sum(comp == "pelagico", na.rm = TRUE), misto = sum(comp == "misto", na.rm = TRUE),
            valores_worms = paste(unique(na.omit(fg)), collapse = " / "), .groups = "drop")

# ---- FishBase: DemersPelag of all fish species ----
esp_peixe <- unique(tx$species[peixe & !is.na(tx$species)])
fb <- tryCatch(rfishbase::species(esp_peixe, fields = c("Species", "DemersPelag")), error = function(e) NULL)
fb_fam <- NULL
if (!is.null(fb)) {
  fb$family <- tx$family[match(fb$Species, tx$species)]
  fb$comp <- case_when(fb$DemersPelag %in% c("demersal", "bathydemersal", "reef-associated") ~ "bentos",
                       fb$DemersPelag %in% c("pelagic-neritic", "pelagic-oceanic", "bathypelagic", "pelagic") ~ "pelagico",
                       fb$DemersPelag %in% "benthopelagic" ~ "misto", TRUE ~ NA_character_)
  fb_fam <- fb |> group_by(family) |> summarise(n_fb = n(), fb_bentos = sum(comp == "bentos", na.rm = TRUE),
    fb_pelagico = sum(comp == "pelagico", na.rm = TRUE), fb_misto = sum(comp == "misto", na.rm = TRUE),
    batipelagicos = sum(DemersPelag == "bathypelagic", na.rm = TRUE),
    valores_fishbase = paste(names(sort(table(DemersPelag), decreasing = TRUE)), collapse = " / "), .groups = "drop")
}

# ---- proposal ----
focal_b <- c("Polychaeta", "Bivalvia", "Echinodermata", "Anthozoa", "Porifera", "Cumacea", "Isopoda", "Tanaidacea", "Bryozoa", "Ascidiacea")
focal_p <- c("Calanoida", "Euphausiacea", "Chaetognatha", "Hyperiidea", "Siphonophorae", "Thaliacea")
prop <- fam |> left_join(fg_fam, by = "family")
if (!is.null(fb_fam)) prop <- prop |> left_join(fb_fam, by = "family")
maioria <- function(b, p, m, tot) ifelse(is.na(tot) | tot == 0, "ambiguo",
                                  ifelse(b / tot >= 0.7, "bentos", ifelse(p / tot >= 0.7, "pelagico", "ambiguo")))
prop <- prop |> mutate(
  peixe = grupo_amplo == "Peixes",
  compartimento = case_when(
    grupo_amplo %in% focal_b ~ "bentos", grupo_amplo %in% focal_p ~ "pelagico",
    order == "Globigerinida" ~ "pelagico",
    peixe ~ maioria(coalesce(fb_bentos, 0L), coalesce(fb_pelagico, 0L), coalesce(fb_misto, 0L), n_fb),
    TRUE ~ maioria(bentos, pelagico, misto, n_com_atributo)),
  fonte_compartimento = case_when(
    grupo_amplo %in% c(focal_b, focal_p) ~ "táxon focal do nível 1 (E2)",
    order == "Globigerinida" ~ "foraminíferos planctônicos (E5.2)",
    peixe & !is.na(n_fb) ~ sprintf("FishBase: %s de %s espécies demersais, %s pelágicas", fb_bentos, n_fb, fb_pelagico),
    !peixe & coalesce(n_com_atributo, 0L) > 0 ~ sprintf("WoRMS: %s de %s espécies bentos, %s plâncton/nécton", bentos, n_com_atributo, pelagico),
    TRUE ~ "sem atributo: compartimento pelo registro (aparelho, profundidade)"),
  grupo_ecologico = case_when(
    compartimento == "bentos" & peixe ~ "peixes demersais",
    compartimento == "bentos" & grupo_amplo %in% c("Polychaeta", "Bivalvia", "Cumacea", "Tanaidacea") ~ "infauna",
    compartimento == "bentos" & grupo_amplo %in% c("Porifera", "Anthozoa", "Bryozoa", "Ascidiacea", "Hydrozoa") ~ "epifauna séssil",
    compartimento == "bentos" & grupo_amplo %in% c("Harpacticoida", "Ostracoda", "Foraminifera") ~ "meiofauna",
    compartimento == "bentos" ~ "epifauna móvel e megafauna",
    compartimento == "pelagico" & peixe & coalesce(batipelagicos, 0L) > 0 & family %in% c("Myctophidae", "Gonostomatidae", "Sternoptychidae", "Phosichthyidae") ~ "micronécton",
    compartimento == "pelagico" & peixe ~ "peixes pelágicos pequenos",
    compartimento == "pelagico" & grupo_amplo %in% c("Euphausiacea", "Cephalopoda") ~ "micronécton",
    compartimento == "pelagico" & (grupo_amplo %in% c("Siphonophorae", "Thaliacea", "Hydrozoa")) ~ "gelatinosos",
    compartimento == "pelagico" ~ "zooplâncton holoplanctônico",
    TRUE ~ "a definir pelo registro"),
  fonte_grupo = "convenção por táxon superior (E5.5); meroplâncton e ictioplâncton definidos pelo estágio de vida no registro")
write.csv(prop |> arrange(grupo_amplo, desc(registros)), file.path(raiz, "entrada/taxon_compartimento_proposta.csv"), row.names = FALSE)
res <- prop |> group_by(grupo_amplo, compartimento) |> summarise(familias = n(), registros = sum(registros), .groups = "drop")
write.csv(res, file.path(raiz, "saida/proposta_compartimento_resumo.csv"), row.names = FALSE)
print(as.data.frame(res))
print(prop |> group_by(grupo_ecologico) |> summarise(familias = n(), registros = sum(registros)) |> as.data.frame())
