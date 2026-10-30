# 82_etapa1_n2.R  (Phase 8, stage 1 at level 2: do the EBU/EMU classes separate the fauna within the 8 regions?)
# Same design as 81_etapa1_n1.R, with the level 2 differences:
#  - contrasts = rows of saida/resolucao_n2.csv (region x level x ecological group x slice, at the resolution chosen by
#    rule 3.2); slice "sem_mes" is excluded. Only records whose position allows the resolution (res_max) are used.
#  - unit: cell at the contrast's resolution (benthos) or cell x EMU (pelagic, only status "ok"), purity >= 0.8,
#    >= 10 events; class with >= 5 units.
#  - spatial block: square of BLOCO_KM (argument "bloco=100" or "bloco=50") on the 6.25 km grid; the 2,500 km block
#    of level 1 does not fit in a 500 km window (autocorrelation range > 1,000 km in the survey).
#  - reading by effect size (Gate 4): partial omega2 within blocks and equal-effort test (10 units per
#    class, 200 repetitions); IndVal.g with 9,999 permutations, BH and Holm.
# Outputs: intermediario/etapa1_n2_b<km>/<id>.rds; saida/etapa1_n2_b<km>_permanova.csv, saida/etapa1_n2_b<km>_indval.csv
suppressPackageStartupMessages({library(arrow); library(dplyr); library(parallel)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
source(file.path(raiz, "scripts/funcoes_permanova.R"))
arg <- grep("^bloco=", commandArgs(TRUE), value = TRUE); BLOCO_KM <- if (length(arg)) as.numeric(sub("bloco=", "", arg)) else 100
NUCLEOS <- 14; R_SORTEIOS <- 5; NPERM <- 999; NPERM_IV <- 9999; N_EQ <- 10; B_EQ <- 200; EV_MIN <- 10; UN_MIN <- 5; PUR_MIN <- 0.8; PR_VERSAO <- "v1"
ORD <- c("25km" = 1, "12.5km" = 2, "6.25km" = 3)
# unidade=genero (stage 2a): datasets with >= 80% of the ecological group's records identified to genus
UNIDADE <- if ("unidade=genero" %in% commandArgs(TRUE)) "genero" else "familia"; SUF <- if (UNIDADE == "genero") "_genero" else ""; GEN_MIN <- 0.8
dir_t <- file.path(raiz, sprintf("intermediario/etapa1_n2_b%d%s", BLOCO_KM, SUF)); dir.create(dir_t, showWarnings = FALSE)

tarefas <- read.csv(file.path(raiz, "saida/resolucao_n2.csv")) |> filter(fatia != "sem_mes") |>
  mutate(comp = ifelse(nivel == "emu", "pelagico", "bentos"),
         id = gsub("[^A-Za-z0-9]+", "_", iconv(paste(regiao, nivel, grupo_ecologico, fatia, res_escolhida), to = "ASCII//TRANSLIT")))
stopifnot(!anyDuplicated(tarefas$id))
cat("bloco:", BLOCO_KM, "km | contrastes:", nrow(tarefas), "| ja prontos:", sum(file.exists(file.path(dir_t, paste0(tarefas$id, ".rds")))), "\n")

roda <- function(k) {
  tk <- tarefas[k, ]; saida <- file.path(dir_t, paste0(tk$id, ".rds")); so_pr <- FALSE
  if (file.exists(saida)) { old <- readRDS(saida); if (is.null(old$perm) || identical(old$pr_versao, PR_VERSAO)) return(tk$id); so_pr <- TRUE }   # only the pairs
  suppressPackageStartupMessages({library(arrow); library(dplyr); library(vegan); library(indicspecies); library(permute)})
  source(file.path(raiz, "scripts/funcoes_permanova.R"))
  s <- sub("km", "", tk$res_escolhida); rg <- tk$regiao; ge <- tk$grupo_ecologico; ft <- tk$fatia
  okres <- names(ORD)[ORD >= ORD[tk$res_escolhida]]
  d <- open_dataset(file.path(raiz, sprintf("saida/tabela_certa_%s_n2.parquet", tk$comp))) |>
    filter(regiao == rg, grupo_ecologico == ge, res_max %in% okres)
  if (tk$comp == "bentos") {
    nv <- sub("ebu_", "", tk$nivel); fcol <- if (ft == "estatico") "fatia_estatica" else "fatia_trimestre"
    d <- d |> select(cell = all_of(paste0("cell_", s)), r6, c6, evento, taxon_unidade, genus, dataset_id, fatia = all_of(fcol),
                     classe = all_of(paste0("ebu_", nv, "_", s)), pureza = all_of(paste0("pureza_", nv, "_", s))) |>
      filter(fatia == ft) |> collect() |> mutate(unidade = cell)
  } else {
    d <- d |> select(cell = all_of(paste0("cell_", s)), r6, c6, evento, taxon_unidade, genus, dataset_id, fatia = fatia_mes,
                     classe = all_of(paste0("emu_", s)), pureza = all_of(paste0("pureza_emu_", s)), situacao = all_of(paste0("situacao_emu_", s))) |>
      filter(fatia == ft, situacao == "ok") |> collect() |> mutate(classe = as.character(classe), unidade = paste(cell, classe))
  }
  if (UNIDADE == "genero") {
    cg <- d |> group_by(dataset_id) |> summarise(fg = mean(!is.na(genus) & genus != "")) |> filter(fg >= GEN_MIN)
    d <- d |> filter(dataset_id %in% cg$dataset_id) |> mutate(taxon_unidade = genus)
  }
  lado <- BLOCO_KM / 6.25
  d <- d |> filter(!is.na(classe), pureza >= PUR_MIN, !is.na(taxon_unidade), taxon_unidade != "") |>
    mutate(classe = as.character(classe), bloco_espacial = paste((r6 - 1) %/% lado, (c6 - 1) %/% lado, sep = "_"))
  # the unit's block is that of most of its records (12.5/25 km cells can cross the edge of a 50 km block)
  un <- d |> group_by(unidade, classe) |> summarise(ev = n_distinct(evento), bloco_espacial = names(which.max(table(bloco_espacial))), .groups = "drop") |>
    filter(ev >= EV_MIN)
  ca <- un |> count(classe) |> filter(n >= UN_MIN); un <- un |> filter(classe %in% ca$classe)
  res <- list(tarefa = tk, unidades = nrow(un), classes = nrow(ca), n_por_classe = ca, perm = NULL, indval = NULL)
  if (nrow(ca) < 2) { saveRDS(res, saida); return(tk$id) }
  ev <- d |> filter(unidade %in% un$unidade) |> distinct(unidade, evento, taxon_unidade); evs <- ev |> distinct(unidade, evento)
  meta <- un |> arrange(unidade); grp <- factor(meta$classe); blk <- factor(meta$bloco_espacial)
  bm <- meta |> group_by(bloco_espacial) |> summarise(k = n_distinct(classe)) |> filter(k >= 2)
  cm <- meta |> filter(bloco_espacial %in% bm$bloco_espacial) |> count(classe) |> filter(n >= UN_MIN)
  im <- which(meta$bloco_espacial %in% bm$bloco_espacial & meta$classe %in% cm$classe)
  im <- im[meta$bloco_espacial[im] %in% (meta[im, ] |> group_by(bloco_espacial) |> filter(n_distinct(classe) >= 2) |> pull(bloco_espacial))]
  grp_m <- droplevels(grp[im]); blk_m <- droplevels(blk[im]); testa_m <- nlevels(grp_m) >= 2
  ctrl_m <- how(blocks = blk_m, nperm = NPERM)
  perm <- list()
  for (r in if (so_pr) 1 else seq_len(R_SORTEIOS)) {
    set.seed(1000 + r)
    sor <- evs |> group_by(unidade) |> slice_sample(n = EV_MIN) |> ungroup()
    m <- ev |> semi_join(sor, by = c("unidade", "evento")) |> distinct(unidade, taxon_unidade) |> mutate(v = 1L) |>
      tidyr::pivot_wider(names_from = taxon_unidade, values_from = v, values_fill = 0L) |> arrange(unidade)
    stopifnot(identical(m$unidade, meta$unidade)); X <- as.matrix(m[, -1]); X <- X[, colSums(X) > 0, drop = FALSE]
    for (ind in if (so_pr) character(0) else if (r == 1) c("sim", "sor") else "sim") {
      D <- sqrt(if (ind == "sim") betadiver(X, "sim") else vegdist(X, "bray", binary = TRUE))
      set.seed(2000 + r); a <- permanova_rapida(D, grp, NULL, NPERM)
      lin <- data.frame(sorteio = r, indice = ind, unidades = nrow(X), taxons = ncol(X), classes = nlevels(grp),
        R2_global = a$R2, F_global = a$F, p_livre = a$p, unidades_mistas = length(im), blocos_mistos = nlevels(blk_m),
        classes_mistas = nlevels(grp_m), R2_dentro = NA, F_dentro = NA, p_dentro = NA, R2_bloco = NA, disp_F = NA, disp_p = NA)
      if (testa_m) {
        Dm <- as.dist(as.matrix(D)[im, im])
        set.seed(3000 + r); b <- permanova_rapida(Dm, grp_m, blk_m, controle = ctrl_m)
        lin$R2_dentro <- b$R2; lin$F_dentro <- b$F; lin$p_dentro <- b$p; lin$R2_bloco <- b$R2_bloco
        if (r == 1 && length(im) > 3) { set.seed(4000); pb <- permutest(betadisper(Dm, grp_m), permutations = ctrl_m)$tab
          lin$disp_F <- pb$F[1]; lin$disp_p <- pb$`Pr(>F)`[1] }
      }
      perm[[length(perm) + 1]] <- lin
    }
    if (r == 1) { res$habitats_dados <- levels(grp); res$habitats_testados <- if (testa_m) levels(grp_m) else character(0) }
    if (r == 1 && testa_m) {   # habitat pairs (which are distinct from each other)
      Mp <- as.matrix(sqrt(betadiver(X[im, , drop = FALSE], "sim"))); set.seed(7000)
      res$pares <- pares_classes(Mp, grp_m, blk_m, un_min = UN_MIN, n_classe = N_EQ, B = 100); res$pr_versao <- PR_VERSAO; rm(Mp)
    }
    if (so_pr) { res$eq <- old$eq; res$indval <- old$indval; res$perm <- old$perm; saveRDS(res, saida); return(tk$id) }
    if (r == 1 && testa_m) {
      Mi <- as.matrix(sqrt(betadiver(X[im, , drop = FALSE], "sim"))); set.seed(6000)
      res$eq <- permanova_esforco_igual(Mi, grp_m, blk_m, n_classe = N_EQ, B = B_EQ, nperm = 199)
      Xm <- X[im, colSums(X[im, , drop = FALSE]) > 0, drop = FALSE]
      if (ncol(Xm) >= 2) {
        ctrl_iv <- how(blocks = blk_m, nperm = NPERM_IV)
        set.seed(5000); iv <- multipatt(as.data.frame(Xm), grp_m, func = "IndVal.g", duleg = TRUE, control = ctrl_iv)
        sg <- iv$sign; sg$taxon <- rownames(sg); sg$classe <- levels(grp_m)[sg$index]
        sg$p_holm <- p.adjust(sg$p.value, "holm"); sg$p_bh <- p.adjust(sg$p.value, "BH"); sg$freq <- colMeans(Xm)[sg$taxon]
        sg$freq_classe <- vapply(seq_len(nrow(sg)), function(j) mean(Xm[grp_m == sg$classe[j], sg$taxon[j]]), 0)
        res$indval <- sg |> select(taxon, classe, stat, p.value, p_holm, p_bh, freq, freq_classe)
      }
    }
  }
  res$perm <- bind_rows(perm); saveRDS(res, saida); tk$id
}

cl <- makeCluster(NUCLEOS)
clusterExport(cl, c("roda", "tarefas", "raiz", "dir_t", "BLOCO_KM", "R_SORTEIOS", "NPERM", "NPERM_IV", "N_EQ", "B_EQ", "EV_MIN", "UN_MIN", "PUR_MIN", "ORD", "PR_VERSAO", "UNIDADE", "GEN_MIN"))
t0 <- Sys.time()
feitos <- unlist(parLapplyLB(cl, seq_len(nrow(tarefas)), function(k) tryCatch(roda(k), error = function(e) paste("ERRO", tarefas$id[k], conditionMessage(e)))))
stopCluster(cl)
print(feitos[grepl("^ERRO", feitos)]); cat("tempo:", format(Sys.time() - t0), "\n")

# ---- consolidation ----
chave <- c("regiao", "nivel", "grupo_ecologico", "fatia", "res_escolhida")
rs <- lapply(file.path(dir_t, paste0(tarefas$id, ".rds")), function(a) if (file.exists(a)) readRDS(a))
pm <- bind_rows(lapply(rs, function(r) { if (is.null(r)) return(NULL); b <- r$tarefa[, chave]
  if (is.null(r$perm)) return(cbind(b, unidades = r$unidades, classes = r$classes, contraste = "não viável"))
  s <- r$perm |> mutate(omega2_global = omega2_parcial(F_global, classes - 1, unidades),
                        omega2_dentro = omega2_parcial(F_dentro, classes_mistas - 1, unidades_mistas)) |> group_by(indice) |>
    summarise(across(c(unidades, taxons, classes, unidades_mistas, blocos_mistos, classes_mistas), first),
      R2_dentro_min = min(R2_dentro), R2_dentro_max = max(R2_dentro), omega2_dentro_min = min(omega2_dentro), omega2_dentro_max = max(omega2_dentro),
      R2_global = median(R2_global), omega2_global = median(omega2_global), p_livre_max = max(p_livre), R2_dentro = median(R2_dentro),
      omega2_dentro = median(omega2_dentro), F_dentro = median(F_dentro), p_dentro_max = max(p_dentro), R2_bloco = median(R2_bloco),
      disp_p = suppressWarnings(max(disp_p, na.rm = TRUE)), .groups = "drop") |> mutate(disp_p = ifelse(is.finite(disp_p), disp_p, NA))
  if (!is.null(r$eq)) { eq <- r$eq[rep(1, nrow(s)), ]; eq[s$indice != "sim", ] <- NA; s <- cbind(s, eq) }
  cbind(b, s, contraste = ifelse(s$classes_mistas < 2, "sem blocos mistos", ifelse(is.na(s$p_dentro_max), "sem variação dentro dos blocos", "testado"))) }))
pm$bloco_km <- BLOCO_KM
write.csv(pm, file.path(raiz, sprintf("saida/etapa1_n2_b%d%s_permanova.csv", BLOCO_KM, SUF)), row.names = FALSE)
hb <- bind_rows(lapply(rs, function(r) if (!is.null(r$habitats_dados)) data.frame(r$tarefa[rep(1, length(r$habitats_dados)), chave],
  habitat = r$habitats_dados, testado = r$habitats_dados %in% r$habitats_testados, unidades = r$n_por_classe$n[match(r$habitats_dados, r$n_por_classe$classe)])))
hb$bloco_km <- BLOCO_KM; write.csv(hb, file.path(raiz, sprintf("saida/etapa1_n2_b%d%s_habitats.csv", BLOCO_KM, SUF)), row.names = FALSE)
pr <- bind_rows(lapply(rs, function(r) if (!is.null(r$pares) && nrow(r$pares)) cbind(r$tarefa[rep(1, nrow(r$pares)), chave], r$pares)))
pr$bloco_km <- BLOCO_KM; write.csv(pr, file.path(raiz, sprintf("saida/etapa1_n2_b%d%s_pares.csv", BLOCO_KM, SUF)), row.names = FALSE)
iv <- bind_rows(lapply(rs, function(r) if (!is.null(r$indval)) cbind(r$tarefa[, chave], r$indval)))
write.csv(iv, file.path(raiz, sprintf("saida/etapa1_n2_b%d%s_indval.csv", BLOCO_KM, SUF)), row.names = FALSE)
print(pm |> filter(is.na(indice) | indice == "sim") |> count(nivel, contraste))
