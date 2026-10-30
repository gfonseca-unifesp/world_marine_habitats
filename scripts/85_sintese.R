# 85_sintese.R  (Phase 8, end of stage 2: conclusion of the bottom-up method)
# Final list of habitats (level 1, global, static) and reclassified global map.
# Status of each habitat, from the consensus of pairs across taxonomic groups (84_equivalencias.R; family or genus):
#   "distingue-se por si so"  = distinct on its own: all tested pairs are distinct;
#   "distingue-se em parte"   = partly distinct: at least one distinct pair and others equivalent or "depende do grupo";
#   "nao se distingue"        = not distinct: no distinct pair and at least one equivalent (enters a merge when all pairs of the
#                               set are equivalent);
#   "inconclusivo"            = inconclusive: no distinct nor equivalent pair, only "depende do grupo" (stays alone on the map);
#   "stand-by: sem par"       = has data, but shares no block with another habitat (not yet comparable);
#   "stand-by: sem dados"     = does not have enough units in any group.
# Map: merged habitats get the same color (new class Un; "U" for union, because F1-F8 are already EBU regions); the others keep their identity; stand-by with its own color and hatched (Gate 8).
# Classes: EBU group (6) and EBU region (57) in the benthos; EMU (37) in the pelagic, mapped at the surface and at 1,000 m.
# Outputs: saida/sintese_habitats_n1.csv, saida/habitats_final_<nivel>_25km.tif (+ legend .csv), saida/mapa_habitats_<nivel>.png
suppressPackageStartupMessages({library(dplyr); library(arrow); library(terra)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); source(file.path(raiz, "scripts/funcoes_grade.R"))
eq <- read.csv(file.path(raiz, "saida/etapa2_n1_equivalencias.csv")) |> filter(fatia_tipo == "estatico") |>
  mutate(classe_a = as.character(classe_a), classe_b = as.character(classe_b))
fu <- read.csv(file.path(raiz, "saida/etapa2_n1_fusoes.csv")); if (nrow(fu)) fu <- fu |> filter(fatia_tipo == "estatico")
hb <- read.csv(file.path(raiz, "saida/etapa1_n1_habitats.csv")) |> filter(!sem_mlpa, fatia == "estatico") |> mutate(habitat = as.character(habitat))
e25 <- read_parquet(file.path(raiz, "camadas/ebu/ebu_celula_25km.parquet")); emu <- rast(file.path(raiz, "camadas/emu/emu_25km.tif"))
EBUN <- c(A = "plataforma", B = "talude", C = "montes submarinos", D = "dorsais e riftes", E = "abissal e hadal", F = "platôs")
emu_cl <- read.csv(file.path(raiz, "camadas/emu/emu_classes.csv"))
universo <- list("EBU grupo" = sort(unique(e25$ebu_grupo)), "EBU região" = sort(unique(e25$ebu_regiao)),
                 "EMU" = as.character(sort(setdiff(unique(as.vector(values(emu))), c(NA, 0)))))

# regional (level 2): shelf x slope per region
p2 <- tryCatch(read.csv(file.path(raiz, "saida/etapa2_n2_pares_final.csv")), error = function(e) NULL)
reg2 <- if (!is.null(p2)) p2 |> filter(nivel == "ebu_grupo", !is.na(separa)) |> group_by(regiao) |> summarise(d = any(separa)) else NULL

# Gate 6, decision 1(a): shelf (A) x slope (B) distinct by the local evidence of level 2; the unions of the EBU group
# are recomputed with this constraint (same rule as 84_equivalencias.R: complete clustering, cut at 1/3).
ab <- eq$nivel == "EBU grupo" & ((eq$classe_a == "A" & eq$classe_b == "B") | (eq$classe_a == "B" & eq$classe_b == "A"))
eq$consenso[ab] <- "distinto"; eq$frac_distinto[ab] <- 1; eq$distinto_em[ab] <- "nível 2 (Gate 6)"
e <- eq |> filter(nivel == "EBU grupo"); h <- sort(unique(c(e$classe_a, e$classe_b)))
M <- matrix(1, length(h), length(h), dimnames = list(h, h)); diag(M) <- 0
M[cbind(e$classe_a, e$classe_b)] <- ifelse(e$grupos_testados >= 2, e$frac_distinto, 1); M[cbind(e$classe_b, e$classe_a)] <- M[cbind(e$classe_a, e$classe_b)]
cl <- cutree(hclust(as.dist(M), method = "complete"), h = 1/3); g <- split(names(cl), cl); g <- g[lengths(g) > 1]
fu <- bind_rows(fu |> filter(nivel != "EBU grupo"),
                if (length(g)) data.frame(nivel = "EBU grupo", fatia_tipo = "estatico", fusao = seq_along(g),
                                          habitats = vapply(g, paste, "", collapse = " + "), n_habitats = lengths(g)))

sint <- bind_rows(lapply(names(universo), function(nv) {
  e <- eq |> filter(nivel == nv); h <- universo[[nv]]
  lado <- bind_rows(e |> transmute(habitat = classe_a, outro = classe_b, consenso, grupos_testados, distinto_em),
                    e |> transmute(habitat = classe_b, outro = classe_a, consenso, grupos_testados, distinto_em))
  r <- lado |> group_by(habitat) |> summarise(pares = n(), distinto_de = paste(sort(outro[consenso == "distinto"]), collapse = ", "),
    equivalente_a = paste(sort(outro[consenso == "equivalente"]), collapse = ", "), depende_de = paste(sort(outro[consenso == "depende do grupo"]), collapse = ", "),
    n_distinto = sum(consenso == "distinto"), n_equiv = sum(consenso == "equivalente"), n_dep = sum(consenso == "depende do grupo"), .groups = "drop")
  dados <- hb |> filter(nivel == nv) |> group_by(habitat) |> summarise(grupos_com_dados = paste(sort(unique(grupo)), collapse = ", "), unidades = sum(unidades, na.rm = TRUE))
  x <- data.frame(nivel = nv, habitat = h) |> left_join(r, by = "habitat") |> left_join(dados, by = "habitat") |>
    mutate(pares = coalesce(pares, 0L),
           status = case_when(pares > 0 & n_distinto == pares ~ "distingue-se por si só",
                              pares > 0 & n_distinto > 0 ~ "distingue-se em parte",
                              pares > 0 & n_equiv > 0 ~ "não se distingue",
                              pares > 0 ~ "inconclusivo",
                              !is.na(grupos_com_dados) ~ "stand-by: sem par",
                              TRUE ~ "stand-by: sem dados"))
  f <- if (nrow(fu)) fu |> filter(nivel == nv) else fu
  x$fusao <- NA_character_
  if (nrow(f)) for (j in seq_len(nrow(f))) { m <- strsplit(f$habitats[j], " \\+ ")[[1]]; x$fusao[x$habitat %in% m] <- paste0("U", j, ": ", f$habitats[j]) }
  x$classe_final <- ifelse(is.na(x$fusao), x$habitat, sub(":.*$", "", x$fusao))
  # destination on the map and equivalences the merge does not honor (equivalence is not transitive: the merge requires all pairs equivalent)
  x$destino <- ifelse(grepl("^stand-by", x$status), "stand-by (achurado)", ifelse(is.na(x$fusao), "sozinho", paste("fundido em", sub(":.*$", "", x$fusao))))
  x$equiv_fora_fusao <- vapply(seq_len(nrow(x)), function(i) { eqv <- if (is.na(x$equivalente_a[i]) || x$equivalente_a[i] == "") character(0) else strsplit(x$equivalente_a[i], ", ")[[1]]
    mesmos <- x$habitat[!is.na(x$fusao) & x$fusao %in% x$fusao[i]]; paste(setdiff(eqv, mesmos), collapse = ", ") }, "")
  x$rotulo <- if (nv == "EBU grupo") paste(x$habitat, EBUN[x$habitat]) else if (nv == "EMU") paste("EMU", x$habitat, "-", emu_cl$EMU_Name[match(x$habitat, emu_cl$Cluster37)]) else x$habitat
  if (nv == "EBU grupo" && !is.null(reg2)) x$nivel2 <- ifelse(x$habitat %in% c("A", "B"),
    sprintf("plataforma × talude distintos em %d de %d regiões", sum(reg2$d), nrow(reg2)), "não comparável nas janelas")
  x
}))
write.csv(sint, file.path(raiz, "saida/sintese_habitats_n1.csv"), row.names = FALSE)
print(sint |> count(nivel, status)); print(sint |> filter(!is.na(fusao)) |> distinct(nivel, fusao))

# ---- maps and GeoTIFF ----
g25 <- ease2_grade("25km"); source(file.path(raiz, "scripts/funcoes_mapa.R"))
rc_de <- function(id) { p <- do.call(rbind, strsplit(id, "_")); cellFromRowCol(g25, as.integer(p[, 2]), as.integer(p[, 3])) }
mapa <- function(nv, valores_hab, titulo, arq) {   # Gate 8: each stand-by with its own color, hatched (funcoes_mapa.R)
  s <- sint |> filter(nivel == nv); ativo <- !grepl("^stand-by", s$status)
  ca <- sort(unique(s$classe_final[ativo])); csb <- sort(unique(s$habitat[!ativo])); cls <- c(ca, csb); sb <- cls %in% csb
  cod <- match(ifelse(ativo, s$classe_final, s$habitat), cls); names(cod) <- s$habitat
  r <- rast(g25); r[] <- unname(cod[as.character(valores_hab)])
  cores <- unname(PAL_HAB[seq_along(cls)])
  leg <- data.frame(valor = seq_along(cls), classe_final = cls,
    habitats = vapply(seq_along(cls), function(j) if (sb[j]) cls[j] else paste(s$habitat[ativo & s$classe_final == cls[j]], collapse = " + "), ""),
    status = vapply(seq_along(cls), function(j) paste(unique(s$status[if (sb[j]) s$habitat == cls[j] else ativo & s$classe_final == cls[j]]), collapse = " / "), ""),
    cor = cores, standby = sb)
  writeRaster(r, file.path(raiz, sprintf("saida/habitats_final_%s.tif", arq)), overwrite = TRUE, datatype = "INT2U")
  write.csv(leg, file.path(raiz, sprintf("saida/habitats_final_%s_legenda.csv", arq)), row.names = FALSE)
  lab <- ifelse(sb, paste(cls, "(stand-by)"), ifelse(leg$classe_final %in% s$habitat, leg$classe_final, paste0(leg$classe_final, " (", leg$habitats, ")")))
  mapa_classes(r, cores, sb, lab, titulo, file.path(raiz, sprintf("saida/mapa_habitats_%s.png", arq)), ncol = 7)
  leg
}
c25 <- rc_de(e25$cell_id)
v <- rep(NA, ncell(g25)); v[c25] <- e25$ebu_grupo; l1 <- mapa("EBU grupo", v, "Habitats bentônicos finais: grupos da EBU (nível 1)", "ebu_grupo_25km")
v <- rep(NA, ncell(g25)); v[c25] <- e25$ebu_regiao; l2 <- mapa("EBU região", v, "Habitats bentônicos finais: regiões da EBU (nível 1)", "ebu_regiao_25km")
sup <- values(emu[[nlyr(emu)]])[, 1]; l3 <- mapa("EMU", sup, "Habitats pelágicos finais: EMU na superfície (nível 1)", "emu_superficie_25km")
b1000 <- which(names(emu) == "prof_1000m"); m1000 <- values(emu[[b1000]])[, 1]; l4 <- mapa("EMU", m1000, "Habitats pelágicos finais: EMU a 1.000 m (nível 1)", "emu_1000m_25km")
cat("mapas gravados\n")
