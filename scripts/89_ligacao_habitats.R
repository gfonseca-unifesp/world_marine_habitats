# 89_ligacao_habitats.R  (Phase 8, stage 3.4: does the environment explain the distinct habitats and the merges?)
# For each group and class level (bentos: EBU group and region; pelagico: EMU):
#  - environmental profile of each habitat = median of the predictors over the habitat's units (purity >= 0.8 in benthos);
#  - Gradient Forest transformation: each predictor becomes its cumulative importance (the group's cumimp curve), so that
#    the distance between two habitats in this space is the change in composition that the environment predicts between them;
#  - stage 2 pairs (84_equivalencias.R, family or genus) of the same group: distinct x not distinct; AUC of the GF
#    distance for separating the two (0.5 = chance; 1 = perfect separation) and the predictor contributing most in each pair.
# Robust thresholds: GF thresholds that recur in >= 2 groups (same predictor; complete linkage: all within
#   10% of the predictor's global range, no chaining).
# Outputs: saida/etapa3_n1_habitats_ambiente.csv, saida/etapa3_n1_pares_gf.csv, saida/etapa3_n1_auc.csv,
#         saida/etapa3_n1_limiares_robustos.csv
suppressPackageStartupMessages({library(dplyr)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); dir_t <- file.path(raiz, "intermediario/etapa3_n1")
curvas <- read.csv(file.path(raiz, "saida/etapa3_n1_gf_curvas.csv"))
pares <- read.csv(file.path(raiz, "saida/etapa2_n1_pares_final.csv")) |> filter(fatia == "estatico", !is.na(separa)) |>
  mutate(classe_a = as.character(classe_a), classe_b = as.character(classe_b))
auc <- function(s, y) { if (length(unique(y)) < 2) return(NA_real_); r <- rank(s); n1 <- sum(y); n0 <- sum(!y); (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0) }

hab_amb <- list(); pg <- list()
for (a in list.files(dir_t, pattern = "^(bentos|pelagico)_.*\\.rds$", full.names = TRUE)) {
  r <- readRDS(a); cv <- curvas |> filter(comp == r$comp, grupo == r$grupo); if (!nrow(cv)) next
  niveis <- if (r$comp == "bentos") list("EBU grupo" = c("ebu_grupo", "pureza_grupo"), "EBU região" = c("ebu_regiao", "pureza_regiao")) else list("EMU" = c("emu", "pureza_emu"))
  for (nv in names(niveis)) {
    cl <- as.character(r$meta[[niveis[[nv]][1]]]); ok <- !is.na(cl) & r$meta[[niveis[[nv]][2]]] >= 0.8
    if (sum(ok) < 10) next
    perf <- aggregate(r$E[ok, , drop = FALSE], list(habitat = cl[ok]), median); nn <- table(cl[ok]); perf$unidades <- as.integer(nn[perf$habitat])
    hab_amb[[length(hab_amb) + 1]] <- cbind(comp = r$comp, grupo = r$grupo, nivel = nv, perf)
    # GF transformation
    vs <- intersect(unique(cv$preditor), names(r$E))
    tr <- sapply(vs, function(v) { c1 <- cv[cv$preditor == v, ]; approx(c1$x, c1$y, perf[[v]], rule = 2, ties = mean)$y })
    if (is.null(dim(tr))) tr <- matrix(tr, nrow = nrow(perf), dimnames = list(NULL, vs)); rownames(tr) <- perf$habitat
    p <- pares |> filter(nivel == nv, grupo == r$grupo, classe_a %in% perf$habitat, classe_b %in% perf$habitat)
    if (!nrow(p)) next
    dd <- tr[p$classe_a, , drop = FALSE] - tr[p$classe_b, , drop = FALSE]
    p$dist_gf <- sqrt(rowSums(dd^2)); p$preditor_principal <- colnames(dd)[apply(abs(dd), 1, which.max)]
    p$parcela_principal <- apply(dd^2, 1, max) / pmax(rowSums(dd^2), 1e-12)
    pg[[length(pg) + 1]] <- cbind(comp = r$comp, p |> select(nivel, grupo, classe_a, classe_b, separa, decidido_por, dist_gf, preditor_principal, parcela_principal))
  }
}
ha <- bind_rows(hab_amb); write.csv(ha, file.path(raiz, "saida/etapa3_n1_habitats_ambiente.csv"), row.names = FALSE)
pg <- bind_rows(pg); write.csv(pg, file.path(raiz, "saida/etapa3_n1_pares_gf.csv"), row.names = FALSE)
au <- bind_rows(pg |> group_by(comp, nivel, grupo) |> summarise(pares = n(), distintos = sum(separa), auc = auc(dist_gf, separa),
                  dist_distintos = median(dist_gf[separa]), dist_nao = median(dist_gf[!separa]), .groups = "drop"),
                pg |> group_by(comp, nivel) |> summarise(grupo = "todos", pares = n(), distintos = sum(separa), auc = auc(dist_gf, separa),
                  dist_distintos = median(dist_gf[separa]), dist_nao = median(dist_gf[!separa]), .groups = "drop"))
write.csv(au, file.path(raiz, "saida/etapa3_n1_auc.csv"), row.names = FALSE)

# robust thresholds
lim <- read.csv(file.path(raiz, "saida/etapa3_n1_gf_limiares.csv"))
amp <- curvas |> group_by(preditor) |> summarise(amp = diff(range(x)))
rob <- bind_rows(lapply(split(lim, lim$preditor), function(l) {
  tol <- 0.1 * amp$amp[amp$preditor == l$preditor[1]]; l <- l[order(l$limiar), ]
  g <- if (nrow(l) > 1) cutree(hclust(dist(l$limiar), method = "complete"), h = tol) else 1   # complete linkage: no chaining of distant thresholds
  bind_rows(lapply(split(l, g), function(k) data.frame(preditor = k$preditor[1], limiar = median(k$limiar), de = min(k$limiar), ate = max(k$limiar),
    grupos = length(unique(paste(k$comp, k$grupo))), lista = paste(unique(paste0(k$grupo, " (", k$comp, ")")), collapse = ", "), intensidade_media = mean(k$intensidade))))
})) |> filter(grupos >= 2) |> arrange(preditor, limiar)
write.csv(rob, file.path(raiz, "saida/etapa3_n1_limiares_robustos.csv"), row.names = FALSE)
options(width = 200); print(au |> mutate(across(where(is.numeric), ~ round(.x, 3))) |> as.data.frame()); print(rob)
print(pg |> filter(separa) |> count(nivel, preditor_principal, sort = TRUE) |> head(15))
