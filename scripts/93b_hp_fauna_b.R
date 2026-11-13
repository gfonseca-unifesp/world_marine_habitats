# 93b_hp_fauna_b.R  (PI decision of 06/11/2026: the option-B fauna tests the approved pelagic habitats, instead of rebuilding them)
# The k-means groups of 89b stay as approved (K0 = 30, level-1 GF transformation). The option-B units (etapa3_n1b: micronekton and small
# pelagic fishes, genus level, pelagic gear only) get the nearest level-1 centroid, as the level-2 units in 93_pelagico_n2.R; pairs per
# group (pares_classes, omega2 >= 0.01 and power >= 0.8); consensus and merge rule of 89b/93 on level 1 + level 2 + option B.
# Double counting: the level-2 and option-B versions of a group share OBIS records, so one pair x ecological group counts once,
# keeping the option-B row (global, gear-filtered) when both exist.
# Outputs: saida/pelagico_b_{pares,consenso,fusoes}.csv (fusoes: k0, hp_n1, final, standby, testado_b; read by 91b_mapa_unido_pelagico.R)
suppressPackageStartupMessages({library(dplyr); library(arrow); library(parallel)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); K0 <- 30
curvas_p <- open_dataset(file.path(raiz, "saida/etapa3_n1_gf_curvas.csv"), format = "csv") |> filter(comp == "pelagico") |> collect() |> as.data.frame()
transf <- function(E) { vs <- intersect(unique(curvas_p$preditor), names(E))   # as in 89b (level-1 pelagic curves)
  sapply(vs, function(v) { s <- 0; for (g in unique(curvas_p$grupo[curvas_p$preditor == v])) { c1 <- curvas_p[curvas_p$grupo == g & curvas_p$preditor == v, ]
    s <- s + approx(c1$x, c1$y, E[[v]], rule = 2, ties = mean)$y }; s }) }
mais_proximo <- function(T, C) { d <- sapply(seq_len(nrow(C)), function(k) colSums((t(T) - C[k, ])^2)); max.col(-d, ties.method = "first") }
mod <- readRDS(file.path(raiz, "intermediario/etapa3_n1/habitats_novos_modelo.rds")); CP <- mod$kmP$centers; FP <- mod$FP; rm(mod); gc(verbose = FALSE)
arqs <- list.files(file.path(raiz, "intermediario/etapa3_n1b"), pattern = "^pelagico_.*[.]rds$", full.names = TRUE)
U <- lapply(arqs, function(a) { r <- readRDS(a); list(arq = a, grupo = r$grupo, k0 = mais_proximo(transf(r$E), CP)) })
testa <- function(u) { suppressPackageStartupMessages({library(vegan); library(dplyr)}); source(file.path(raiz, "scripts/funcoes_permanova.R")); r <- readRDS(u$arq)
  Dm <- as.matrix(sqrt(betadiver(r$X, "sim"))); set.seed(7110); x <- pares_classes(Dm, as.character(u$k0), r$meta$bloco, un_min = 5, n_classe = 10, B = 100)
  if (!is.null(x)) cbind(grupo = u$grupo, x) }
cl <- makeCluster(length(U)); clusterExport(cl, "raiz"); pb <- bind_rows(parLapply(cl, U, testa)); stopCluster(cl)
pb <- pb |> filter(!is.na(separa)); write.csv(pb, file.path(raiz, "saida/pelagico_b_pares.csv"), row.names = FALSE)
sem_ac <- function(x) chartr("áéíóúâêôãõç", "aeiouaeoaoc", x)
ord <- function(d) d |> mutate(a = pmin(as.integer(classe_a), as.integer(classe_b)), b = pmax(as.integer(classe_a), as.integer(classe_b)))
pn <- ord(read.csv(file.path(raiz, "saida/pelagico_n2_pares.csv"), encoding = "UTF-8")) |> transmute(a, b, fam = sem_ac(grupo), separa, fonte = "n2")
pB <- ord(pb) |> transmute(a, b, fam = sem_ac(grupo), separa, fonte = "B")
p2 <- bind_rows(pB, pn) |> distinct(a, b, fam, .keep_all = TRUE)                                # option-B row first: it wins
cat("linhas de nivel 2 retiradas por repeticao com a opcao B:", nrow(pn) + nrow(pB) - nrow(p2), "\n")
n2b <- p2 |> group_by(classe_a = a, classe_b = b) |> summarise(grupos_n2b = n(), distintos_n2b = sum(separa), .groups = "drop")
n1 <- read.csv(file.path(raiz, "saida/etapa3_n1_habitats_novos_pares.csv")) |> filter(comp == "pelagico") |>
  transmute(a = pmin(classe_a, classe_b), b = pmax(classe_a, classe_b), grupos_n1 = grupos, distintos_n1 = distintos) |> rename(classe_a = a, classe_b = b)
cons <- full_join(n1, n2b, by = c("classe_a", "classe_b")) |> mutate(across(starts_with(c("grupos", "distintos")), ~ coalesce(.x, 0L))) |>
  mutate(grupos = grupos_n1 + grupos_n2b, distintos = distintos_n1 + distintos_n2b, frac = distintos / grupos,
         consenso = ifelse(grupos >= 2 & frac <= 1/3, "equivalente", ifelse(frac >= 2/3, "distinto", "depende do grupo")))
write.csv(cons, file.path(raiz, "saida/pelagico_b_consenso.csv"), row.names = FALSE)
Mx <- matrix(1, K0, K0); diag(Mx) <- 0; Mx[cbind(cons$classe_a, cons$classe_b)] <- ifelse(cons$grupos >= 2, cons$frac, 1); Mx[cbind(cons$classe_b, cons$classe_a)] <- Mx[cbind(cons$classe_a, cons$classe_b)]
g <- cutree(hclust(as.dist(Mx), method = "complete"), h = 1/3); fin <- match(g, unique(g)); te <- (1:K0) %in% c(cons$classe_a, cons$classe_b)
fus <- data.frame(k0 = 1:K0, hp_n1 = FP$habitat, final = fin, standby = !tapply(te, fin, any)[as.character(fin)], testado_b = (1:K0) %in% c(pB$a, pB$b))
ant <- read.csv(file.path(raiz, "saida/pelagico_n2_fusoes.csv"))
write.csv(fus, file.path(raiz, "saida/pelagico_b_fusoes.csv"), row.names = FALSE)
cat("pares opcao B:", nrow(pb), "por grupo:", paste(names(table(pb$grupo)), table(pb$grupo), collapse = "; "), "| separam:", sum(pb$separa), "\n")
cat("HP antes:", length(unique(ant$final_n1n2)), "stand-by", sum(!tapply(!ant$standby_n1n2, ant$final_n1n2, any)), "| agora:", length(unique(fus$final)), "stand-by", sum(!tapply(!fus$standby, fus$final, any)), "\n")
print(table(cons$consenso)); print(fus |> mutate(antes = ant$final_n1n2, standby_antes = ant$standby_n1n2) |> filter(standby != standby_antes | ave(antes, final, FUN = function(x) length(unique(x))) > 1))
cat(format(Sys.time()), "FIM\n")
