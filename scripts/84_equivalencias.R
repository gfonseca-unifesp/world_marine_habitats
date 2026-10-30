# 84_equivalencias.R  (Phase 8, stage 2b: classes with no difference between them, and 2a: family x genus)
# Habitat pair, in a taxonomic group: "distinto" if it separates at family level OR at genus level (guide v3: genus when family does not
# show a difference); "nao distinto" if it does not separate at family level and (does not separate at genus level or genus is not testable).
# Consensus across groups (level 1: taxonomic groups; level 2: ecological groups within each region):
#   equivalente = tested in >= 2 groups and not distinct in >= 2/3 of them; distinto = distinct in >= 2/3; otherwise "depende do grupo".
# Merge proposal: complete hierarchical clustering of the habitats with distance = fraction of groups in which the pair is distinct
#   (never-tested pair = 1, does not merge), cut at 1/3: only habitats in which ALL pairs of the group are equivalent are merged.
# Outputs: saida/etapa2_genero_comparacao.csv, saida/etapa2_n<k>_pares_final.csv, saida/etapa2_n<k>_equivalencias.csv,
#         saida/etapa2_n<k>_fusoes.csv
suppressPackageStartupMessages({library(dplyr)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); le <- function(f) { a <- file.path(raiz, "saida", f); if (file.exists(a)) read.csv(a) else NULL }

# ---- 2a: family x genus, per contrast ----
cmp <- function(fam, gen, chave) {
  if (is.null(gen)) return(NULL)
  f <- fam |> filter(indice == "sim", contraste == "testado") |> select(all_of(chave), omega2_fam = omega2_dentro, poder_fam = poder_eq, taxons_fam = taxons, unid_fam = unidades_mistas)
  g <- gen |> filter(indice == "sim", contraste == "testado") |> select(all_of(chave), omega2_gen = omega2_dentro, poder_gen = poder_eq, taxons_gen = taxons, unid_gen = unidades_mistas)
  full_join(f, g, by = chave) |> mutate(separa_fam = omega2_fam >= 0.01 & poder_fam >= 0.8, separa_gen = omega2_gen >= 0.01 & poder_gen >= 0.8,
    resultado = case_when(is.na(separa_fam) ~ "só gênero testável", separa_fam ~ "separa em família", separa_gen %in% TRUE ~ "separa só em gênero",
                          is.na(separa_gen) ~ "não separa (gênero não testável)", TRUE ~ "não separa em família nem em gênero"))
}
k1 <- c("comp", "grupo", "nivel", "fatia"); k2 <- c("regiao", "nivel", "grupo_ecologico", "fatia", "res_escolhida")
f1 <- le("etapa1_n1_permanova.csv"); f1 <- f1[!f1$sem_mlpa, ]; g1 <- le("etapa1_n1_genero_permanova.csv")
c1 <- cmp(f1, g1, k1); c2 <- cmp(le("etapa1_n2_b100_permanova.csv"), le("etapa1_n2_b100_genero_permanova.csv"), k2)
comp_gen <- bind_rows(if (!is.null(c1)) mutate(c1, nivel_analise = 1L), if (!is.null(c2)) mutate(c2, nivel_analise = 2L))
write.csv(comp_gen, file.path(raiz, "saida/etapa2_genero_comparacao.csv"), row.names = FALSE)
print(comp_gen |> count(nivel_analise, resultado))

# ---- final pairs (family + genus) ----
par_final <- function(pf, pg, chave) {
  pf <- pf |> mutate(classe_a = as.character(classe_a), classe_b = as.character(classe_b))
  x <- pf |> select(all_of(chave), classe_a, classe_b, blocos, omega2_fam = omega2, poder_fam = poder_eq, separa_fam = separa)
  if (!is.null(pg)) x <- x |> left_join(pg |> mutate(classe_a = as.character(classe_a), classe_b = as.character(classe_b)) |>
    select(all_of(chave), classe_a, classe_b, omega2_gen = omega2, poder_gen = poder_eq, separa_gen = separa), by = c(chave, "classe_a", "classe_b"))
  else x$separa_gen <- NA
  x |> mutate(separa = ifelse(is.na(separa_fam) & is.na(separa_gen), NA, (separa_fam %in% TRUE) | (separa_gen %in% TRUE)),
              decidido_por = case_when(is.na(separa) ~ "não testável", separa_fam %in% TRUE ~ "família", separa_gen %in% TRUE ~ "gênero",
                                       TRUE ~ "nenhum"))
}
consenso <- function(pf, grupo_col, estrato) {
  pf |> filter(!is.na(separa)) |> group_by(across(all_of(c(estrato, "classe_a", "classe_b")))) |>
    summarise(grupos_testados = n(), grupos_distintos = sum(separa), grupos_lista = paste(sort(unique(.data[[grupo_col]])), collapse = ", "),
              distinto_em = paste(sort(unique(.data[[grupo_col]][separa])), collapse = ", "), .groups = "drop") |>
    mutate(frac_distinto = grupos_distintos / grupos_testados,
           consenso = ifelse(grupos_testados >= 2 & frac_distinto <= 1/3, "equivalente", ifelse(frac_distinto >= 2/3, "distinto", "depende do grupo")))
}
fusoes <- function(eq, estrato) {
  bind_rows(lapply(split(eq, eq[estrato], drop = TRUE), function(e) {
    h <- sort(unique(c(e$classe_a, e$classe_b))); if (length(h) < 2) return(NULL)
    M <- matrix(1, length(h), length(h), dimnames = list(h, h)); diag(M) <- 0
    M[cbind(e$classe_a, e$classe_b)] <- ifelse(e$grupos_testados >= 2, e$frac_distinto, 1)
    M[cbind(e$classe_b, e$classe_a)] <- M[cbind(e$classe_a, e$classe_b)]
    cl <- cutree(hclust(as.dist(M), method = "complete"), h = 1/3)
    g <- split(names(cl), cl); g <- g[lengths(g) > 1]; if (!length(g)) return(NULL)
    data.frame(e[1, estrato, drop = FALSE], fusao = seq_along(g), habitats = vapply(g, paste, "", collapse = " + "), n_habitats = lengths(g), row.names = NULL)
  }))
}
for (nv in 1:2) {
  if (nv == 1) { pf <- le("etapa1_n1_pares.csv"); pf <- pf[!pf$sem_mlpa, ]; pg <- le("etapa1_n1_genero_pares.csv"); ch <- k1; gcol <- "grupo"; est <- c("nivel", "fatia_tipo") }
  else { pf <- le("etapa1_n2_b100_pares.csv"); pg <- le("etapa1_n2_b100_genero_pares.csv"); ch <- k2; gcol <- "grupo_ecologico"; est <- c("regiao", "nivel", "fatia_tipo") }
  if (is.null(pf)) next
  pfin <- par_final(pf, pg, ch) |> mutate(fatia_tipo = ifelse(fatia == "estatico", "estatico", "sazonal"))
  write.csv(pfin, file.path(raiz, sprintf("saida/etapa2_n%d_pares_final.csv", nv)), row.names = FALSE)
  eq <- consenso(pfin, gcol, est); write.csv(eq, file.path(raiz, sprintf("saida/etapa2_n%d_equivalencias.csv", nv)), row.names = FALSE)
  fu <- fusoes(eq, est); write.csv(if (is.null(fu)) data.frame() else fu, file.path(raiz, sprintf("saida/etapa2_n%d_fusoes.csv", nv)), row.names = FALSE)
  cat("\nnivel", nv, "\n"); print(pfin |> count(nivel, decidido_por)); print(eq |> count(across(all_of(est)), consenso)); if (!is.null(fu)) print(fu)
}
