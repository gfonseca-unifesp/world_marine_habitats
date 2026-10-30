# 81_etapa1_n1.R  (Phase 8, stage 1 at level 1: do the EBU and EMU classes separate the fauna?)
# Sampling unit: 25 km cell (benthos) or cell x EMU (pelagic) with purity >= 0.8 and >= 10 events.
# Effort standardization: in each unit, draw 10 events (rarefaction by events); the unit's incidence is the
#   union of the taxa from those 10 events. R independent draws (fixed seeds).
# Classes with < 5 units are excluded (rule 6.3); a contrast requires >= 2 classes.
# Dissimilarity: Simpson (beta_sim, taxon turnover only, little sensitive to richness differences); Sorensen as support.
# PERMANOVA (vegan::adonis2 statistic, computed by funcoes_permanova.R) on the square root of the dissimilarity, in two forms:
#   (a) global, D ~ class, free permutation: descriptive (global R2; p is optimistic because of autocorrelation);
#   (b) within the 2,500 km blocks (G3), main test: D ~ block + class, permutations restricted to the blocks,
#       only with blocks that have >= 2 classes (in a single-class block the permutation changes nothing and p goes to 1
#       by construction; seen in the 16/10/2026 trial). betadisper (dispersion) with the same permutations.
# Effect size (takes priority over p, user request on 19/10/2026; large n makes almost everything significant):
#   partial omega2 (R2 corrected for n and number of classes) and equal-effort test: 10 units per class, 200 repetitions,
#   fraction with p < 0.05 = power at the same effort ("effort-weighted p"), comparable across groups.
# IndVal.g (indicspecies::multipatt, one class per taxon) on the units of (b), same permutations; p adjusted
#   by Holm and by Benjamini-Hochberg within each contrast; 9,999 permutations (with 999, the smallest p = 0.001 never passes
#   Holm with ~180 taxa). Rerunning with a new NPERM_IV redoes only the IndVal (same draw 1).
# Slices: static benthos; static pelagic (all seasons and "sem_mes") for the 5 groups and T1-T4 only for
#   Calanoida, Chaetognatha and Euphausiacea (Gate 2). G8 sensitivity: benthos without the MLPA kelp dataset.
# Parallel: PSOCK cluster, one contrast per process; result of each contrast in intermediario/etapa1_n1/<id>.rds
#   (rerunning skips what already exists). Outputs: saida/etapa1_n1_permanova.csv, saida/etapa1_n1_indval.csv
suppressPackageStartupMessages({library(arrow); library(dplyr); library(parallel)})
source(file.path(Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"), "scripts/funcoes_permanova.R"))   # omega2_parcial in the consolidation
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
NUCLEOS <- 14; R_SORTEIOS <- 5; NPERM <- 999; NPERM_IV <- 9999; N_EQ <- 10; B_EQ <- 200; EQ_VERSAO <- "v1-10porclasse"; PR_VERSAO <- "v1"; EV_MIN <- 10; UN_MIN <- 5; PUR_MIN <- 0.8
# unidade=genero (stage 2a, guide v3: "if there is no difference at family level, at genus level, in datasets with consistent identification
#   to genus"): datasets with >= 80% of the group's records identified to genus; Euphausiacea and Chaetognatha are already genus.
UNIDADE <- if ("unidade=genero" %in% commandArgs(TRUE)) "genero" else "familia"; SUF <- if (UNIDADE == "genero") "_genero" else ""; GEN_MIN <- 0.8
dir_t <- file.path(raiz, paste0("intermediario/etapa1_n1", SUF)); dir.create(dir_t, showWarnings = FALSE)
cj <- read.csv(file.path(raiz, "intermediario/conjunto_tipo_escopo.csv"))
MLPA <- cj$dataset_id[grepl("MLPA", cj$titulo)]; MLPA <- MLPA[which.max(cj$registros[cj$dataset_id %in% MLPA])]

# ---- list of contrasts ----
gb <- c("Anthozoa", "Bivalvia", "Bryozoa e Ascidiacea", "Echinodermata", "Peracarida bentonicos", "Polychaeta", "Porifera")
gp <- c("Calanoida", "Chaetognatha", "Euphausiacea", "Hyperiidea", "Siphonophorae e Thaliacea")
tarefas <- bind_rows(
  expand.grid(comp = "bentos", grupo = gb, nivel = c("EBU grupo", "EBU região"), fatia = "estatico", sem_mlpa = c(FALSE, TRUE), stringsAsFactors = FALSE),
  expand.grid(comp = "pelagico", grupo = gp, nivel = "EMU", fatia = "estatico", sem_mlpa = FALSE, stringsAsFactors = FALSE),
  expand.grid(comp = "pelagico", grupo = gp[1:3], nivel = "EMU", fatia = paste0("T", 1:4), sem_mlpa = FALSE, stringsAsFactors = FALSE))
tarefas$id <- gsub("[^A-Za-z0-9]+", "_", paste(tarefas$comp, tarefas$grupo, tarefas$nivel, tarefas$fatia, ifelse(tarefas$sem_mlpa, "semMLPA", "todos")))
if (UNIDADE == "genero") tarefas <- tarefas[!tarefas$sem_mlpa & !tarefas$grupo %in% c("Euphausiacea", "Chaetognatha"), ]
if ("teste" %in% commandArgs(TRUE)) {   # quick trial: 2 contrasts, 99 permutations, 1 draw, temporary folder
  tarefas <- tarefas[tarefas$id %in% c("pelagico_Euphausiacea_EMU_T1_todos", "bentos_Bivalvia_EBU_grupo_estatico_todos"), ]
  NPERM <- 99; R_SORTEIOS <- 1; NUCLEOS <- 2; dir_t <- file.path(tempdir(), "etapa1"); dir.create(dir_t)
}
cat("contrastes:", nrow(tarefas), "| ja prontos:", sum(file.exists(file.path(dir_t, paste0(tarefas$id, ".rds")))), "\n")

roda <- function(k) {
  tk <- tarefas[k, ]; saida <- file.path(dir_t, paste0(tk$id, ".rds")); so_iv <- FALSE; refaz_iv <- TRUE; refaz_eq <- TRUE; refaz_pr <- TRUE
  if (file.exists(saida)) {   # contrast done: redo only what changed version (IndVal and/or equal-effort test), in draw 1
    old <- readRDS(saida); refaz_iv <- !identical(old$nperm_iv, NPERM_IV); refaz_eq <- !identical(old$eq_versao, EQ_VERSAO)
    refaz_pr <- !identical(old$pr_versao, PR_VERSAO)
    if (is.null(old$perm) || !(refaz_iv || refaz_eq || refaz_pr)) return(tk$id); so_iv <- TRUE }
  suppressPackageStartupMessages({library(arrow); library(dplyr); library(vegan); library(indicspecies); library(permute)})
  source(file.path(raiz, "scripts/funcoes_permanova.R"))   # same statistic as adonis2, ~65x faster with blocks
  arq <- file.path(raiz, sprintf("saida/tabela_certa_%s_n1.parquet", tk$comp)); g <- tk$grupo
  d <- open_dataset(arq) |> filter(grupo == g)
  if (tk$comp == "bentos") {
    col <- if (tk$nivel == "EBU grupo") c("ebu_grupo", "pureza_grupo") else c("ebu_regiao", "pureza_regiao")
    d <- d |> select(cell_id, event_id, taxon_unidade, genus, bloco_espacial, dataset, classe = all_of(col[1]), pureza = all_of(col[2]))
    if (tk$sem_mlpa) d <- d |> filter(dataset != MLPA)
    d <- d |> collect() |> mutate(unidade = cell_id)
  } else {
    d <- d |> select(cell_id, event_id, taxon_unidade, genus, dataset, bloco_espacial, fatia, emu, emu_rotulo, pureza = pureza_emu)
    if (tk$fatia != "estatico") { f <- tk$fatia; d <- d |> filter(fatia == f) }
    d <- d |> collect() |> mutate(classe = as.character(emu), unidade = paste(cell_id, emu))
  }
  if (UNIDADE == "genero") {   # only datasets with consistent identification to genus
    cg <- d |> group_by(dataset) |> summarise(fg = mean(!is.na(genus) & genus != "")) |> filter(fg >= GEN_MIN)
    d <- d |> filter(dataset %in% cg$dataset) |> mutate(taxon_unidade = genus)
  }
  d <- d |> filter(!is.na(classe), pureza >= PUR_MIN, !is.na(taxon_unidade), taxon_unidade != "")
  # eligible units: >= 10 events; eligible classes: >= 5 units
  un <- d |> distinct(unidade, classe, bloco_espacial, event_id) |> count(unidade, classe, bloco_espacial, name = "ev") |> filter(ev >= EV_MIN)
  ca <- un |> count(classe) |> filter(n >= UN_MIN); un <- un |> filter(classe %in% ca$classe)
  res <- list(tarefa = tk, unidades = nrow(un), classes = nrow(ca), n_por_classe = ca, perm = NULL, indval = NULL)
  if (nrow(ca) < 2) { saveRDS(res, saida); return(tk$id) }
  ev <- d |> filter(unidade %in% un$unidade) |> distinct(unidade, event_id, taxon_unidade)
  evs <- ev |> distinct(unidade, event_id)
  meta <- un |> arrange(unidade); grp <- factor(meta$classe); blk <- factor(meta$bloco_espacial)
  # within-block test: only blocks with >= 2 classes; classes with >= UN_MIN units in those blocks
  bm <- meta |> group_by(bloco_espacial) |> summarise(k = n_distinct(classe)) |> filter(k >= 2)
  cm <- meta |> filter(bloco_espacial %in% bm$bloco_espacial) |> count(classe) |> filter(n >= UN_MIN)
  im <- which(meta$bloco_espacial %in% bm$bloco_espacial & meta$classe %in% cm$classe)
  im <- im[meta$bloco_espacial[im] %in% (meta[im, ] |> group_by(bloco_espacial) |> filter(n_distinct(classe) >= 2) |> pull(bloco_espacial))]
  grp_m <- droplevels(grp[im]); blk_m <- droplevels(blk[im]); testa_m <- nlevels(grp_m) >= 2
  res$unidades_mistas <- length(im); res$blocos_mistos <- nlevels(blk_m); res$classes_mistas <- nlevels(grp_m)
  ctrl_m <- how(blocks = blk_m, nperm = NPERM)   # without "permute::": multipatt re-evaluates the call
  perm <- list()
  for (r in if (so_iv) 1 else seq_len(R_SORTEIOS)) {
    set.seed(1000 + r)
    sor <- evs |> group_by(unidade) |> slice_sample(n = EV_MIN) |> ungroup()
    m <- ev |> semi_join(sor, by = c("unidade", "event_id")) |> distinct(unidade, taxon_unidade) |> mutate(v = 1L) |>
      tidyr::pivot_wider(names_from = taxon_unidade, values_from = v, values_fill = 0L) |> arrange(unidade)
    stopifnot(identical(m$unidade, meta$unidade)); X <- as.matrix(m[, -1]); X <- X[, colSums(X) > 0, drop = FALSE]
    for (ind in if (so_iv) character(0) else if (r == 1) c("sim", "sor") else "sim") {
      D <- sqrt(if (ind == "sim") betadiver(X, "sim") else vegdist(X, "bray", binary = TRUE))   # square root: avoids negative SS
      # (a) global, free permutation: descriptive (p optimistic because of spatial autocorrelation)
      set.seed(2000 + r); a <- permanova_rapida(D, grp, NULL, NPERM)
      lin <- data.frame(sorteio = r, indice = ind, unidades = nrow(X), taxons = ncol(X), classes = nlevels(grp),
        R2_global = a$R2, F_global = a$F, p_livre = a$p,
        unidades_mistas = length(im), blocos_mistos = nlevels(blk_m), classes_mistas = nlevels(grp_m),
        R2_dentro = NA, F_dentro = NA, p_dentro = NA, R2_bloco = NA, disp_F = NA, disp_p = NA)
      # (b) within blocks: class after block, permutations restricted to the blocks (main test, G3)
      if (testa_m) {
        Dm <- as.dist(as.matrix(D)[im, im])
        set.seed(3000 + r); b <- permanova_rapida(Dm, grp_m, blk_m, controle = ctrl_m)
        lin$R2_dentro <- b$R2; lin$F_dentro <- b$F; lin$p_dentro <- b$p; lin$R2_bloco <- b$R2_bloco
        if (r == 1) { set.seed(4000); pb <- permutest(betadisper(Dm, grp_m), permutations = ctrl_m)$tab   # dispersion: only in the 1st draw
          lin$disp_F <- pb$F[1]; lin$disp_p <- pb$`Pr(>F)`[1] }
        rm(Dm)
      }
      perm[[length(perm) + 1]] <- lin
    }
    if (r == 1 && testa_m && refaz_eq) {   # equal-effort test (effort-weighted p), draw 1, sqrt(beta_sim)
      Mi <- as.matrix(sqrt(betadiver(X[im, , drop = FALSE], "sim"))); set.seed(6000)
      res$eq <- permanova_esforco_igual(Mi, grp_m, blk_m, n_classe = N_EQ, B = B_EQ, nperm = 199); res$eq_versao <- EQ_VERSAO; rm(Mi)
    } else if (so_iv) { res$eq <- old$eq; res$eq_versao <- old$eq_versao }
    if (r == 1 && testa_m && !refaz_iv) { res$indval <- old$indval; res$nperm_iv <- old$nperm_iv }
    if (r == 1) { res$habitats_dados <- levels(grp); res$habitats_testados <- if (testa_m) levels(grp_m) else character(0) }
    if (r == 1 && testa_m && refaz_pr) {   # habitat pairs (which are distinct from each other), draw 1, sqrt(beta_sim)
      Mi <- as.matrix(sqrt(betadiver(X[im, , drop = FALSE], "sim"))); set.seed(7000)
      res$pares <- pares_classes(Mi, grp_m, blk_m, un_min = UN_MIN, n_classe = N_EQ, B = 100); res$pr_versao <- PR_VERSAO; rm(Mi)
    } else if (r == 1 && testa_m && so_iv) { res$pares <- old$pares; res$pr_versao <- old$pr_versao }
    if (r == 1 && testa_m && refaz_iv) {   # IndVal.g on the units of the mixed blocks, restricted permutations
      Xm <- X[im, colSums(X[im, , drop = FALSE]) > 0, drop = FALSE]
      ctrl_iv <- how(blocks = blk_m, nperm = NPERM_IV)
      set.seed(5000); iv <- multipatt(as.data.frame(Xm), grp_m, func = "IndVal.g", duleg = TRUE, control = ctrl_iv)
      s <- iv$sign; s$taxon <- rownames(s); s$classe <- levels(grp_m)[s$index]
      s$p_holm <- p.adjust(s$p.value, "holm"); s$p_bh <- p.adjust(s$p.value, "BH"); s$freq <- colMeans(Xm)[s$taxon]
      s$freq_classe <- vapply(seq_len(nrow(s)), function(j) mean(Xm[grp_m == s$classe[j], s$taxon[j]]), 0)
      res$indval <- s |> select(taxon, classe, stat, p.value, p_holm, p_bh, freq, freq_classe); res$nperm_iv <- NPERM_IV
    }
  }
  res$perm <- if (so_iv) old$perm else bind_rows(perm); saveRDS(res, saida); tk$id
}

ordem <- order(-(tarefas$comp == "bentos" & tarefas$fatia == "estatico"))  # largest first
cl <- makeCluster(NUCLEOS); clusterExport(cl, c("roda", "tarefas","raiz", "dir_t", "MLPA", "R_SORTEIOS", "NPERM", "NPERM_IV", "N_EQ", "B_EQ", "EQ_VERSAO", "PR_VERSAO", "UNIDADE", "GEN_MIN", "EV_MIN", "UN_MIN", "PUR_MIN"))
t0 <- Sys.time()
feitos <- parLapplyLB(cl, ordem, function(k) tryCatch(roda(k), error = function(e) paste("ERRO", tarefas$id[k], conditionMessage(e))))
stopCluster(cl)
print(unlist(feitos)[grepl("^ERRO", unlist(feitos))]); cat("tempo:", format(Sys.time() - t0), "\n")

# ---- consolidation ----
rs <- lapply(file.path(dir_t, paste0(tarefas$id, ".rds")), function(a) if (file.exists(a)) readRDS(a))
pm <- bind_rows(lapply(rs, function(r) { if (is.null(r)) return(NULL); b <- r$tarefa[, c("comp", "grupo", "nivel", "fatia", "sem_mlpa")]
  if (is.null(r$perm)) return(cbind(b, unidades = r$unidades, classes = r$classes, contraste = "não viável"))
  s <- r$perm |> mutate(omega2_global = omega2_parcial(F_global, classes - 1, unidades),
                         omega2_dentro = omega2_parcial(F_dentro, classes_mistas - 1, unidades_mistas)) |> group_by(indice) |> summarise(across(c(unidades, taxons, classes, unidades_mistas, blocos_mistos, classes_mistas), first),
    R2_dentro_min = min(R2_dentro), R2_dentro_max = max(R2_dentro),   # min and max before the median (summarise is sequential)
    omega2_dentro_min = min(omega2_dentro), omega2_dentro_max = max(omega2_dentro), omega2_global = median(omega2_global), omega2_dentro = median(omega2_dentro),
    R2_global = median(R2_global), p_livre_max = max(p_livre), R2_dentro = median(R2_dentro), F_dentro = median(F_dentro), p_dentro_max = max(p_dentro), R2_bloco = median(R2_bloco), disp_F = suppressWarnings(median(disp_F, na.rm = TRUE)),
    disp_p = suppressWarnings(max(disp_p, na.rm = TRUE)), .groups = "drop") |> mutate(disp_p = ifelse(is.finite(disp_p), disp_p, NA))
  eq <- if (!is.null(r$eq)) r$eq[rep(1, nrow(s)), ] else NULL; if (!is.null(eq)) eq[s$indice != "sim", ] <- NA   # equal effort only with sqrt(beta_sim)
  if (!is.null(eq)) s <- cbind(s, eq)
  cbind(b, s, contraste = ifelse(s$classes_mistas < 2, "só global (sem blocos mistos)", ifelse(is.na(s$p_dentro_max), "sem variação dentro dos blocos", "testado"))) }))
if ("teste" %in% commandArgs(TRUE)) { print(as.data.frame(pm)); for (r in rs) print(head(r$indval[order(r$indval$p.value), ], 6)); quit(save = "no") }
write.csv(pm, file.path(raiz, paste0("saida/etapa1_n1", SUF, "_permanova.csv")), row.names = FALSE)
hb <- bind_rows(lapply(rs, function(r) if (!is.null(r$habitats_dados)) data.frame(r$tarefa[rep(1, length(r$habitats_dados)), c("comp", "grupo", "nivel", "fatia", "sem_mlpa")],
  habitat = r$habitats_dados, testado = r$habitats_dados %in% r$habitats_testados,
  unidades = r$n_por_classe$n[match(r$habitats_dados, r$n_por_classe$classe)])))
write.csv(hb, file.path(raiz, paste0("saida/etapa1_n1", SUF, "_habitats.csv")), row.names = FALSE)
pr <- bind_rows(lapply(rs, function(r) if (!is.null(r$pares)) cbind(r$tarefa[rep(1, nrow(r$pares)), c("comp", "grupo", "nivel", "fatia", "sem_mlpa")], r$pares)))
write.csv(pr, file.path(raiz, paste0("saida/etapa1_n1", SUF, "_pares.csv")), row.names = FALSE)
iv <- bind_rows(lapply(rs, function(r) if (!is.null(r$indval)) cbind(r$tarefa[, c("comp", "grupo", "nivel", "fatia", "sem_mlpa")], r$indval)))
write.csv(iv, file.path(raiz, paste0("saida/etapa1_n1", SUF, "_indval.csv")), row.names = FALSE)
print(pm |> filter(is.na(indice) | indice == "sim") |> select(comp, grupo, nivel, fatia, sem_mlpa, unidades, classes, R2_global,
  unidades_mistas, classes_mistas, R2_dentro, omega2_dentro, p_dentro_max, any_of(c("n_med", "omega2_eq", "poder_eq", "p_eq_med")), disp_p) |> as.data.frame())
