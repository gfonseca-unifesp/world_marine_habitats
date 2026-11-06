# 92f_tabela_validacao.R  (GUIA_LIMIARES_PUBLICADOS.md, rule 2: no published value enters the article unchecked against the full text)
# Worksheet for the PI to validate the threshold base one by one: one row per published threshold, sorted by study, with the DOI, the
# values extracted from the abstracts and yellow columns to fill in. DOIs are in entrada/limiares_referencias.csv and were checked on
# Crossref (title, journal and year kept in intermediario/limiares_crossref.rds).
# Output: limiares_validacao_v0.xlsx (project root, beside limiares_publicados_v0.csv)
suppressPackageStartupMessages({library(dplyr); library(openxlsx)})
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot"); arq <- file.path(raiz, "limiares_validacao_v0.xlsx")
if (file.exists(arq)) stop("limiares_validacao_v0.xlsx already exists and may hold the PI's notes: rename it before regenerating")
lp <- read.csv(file.path(raiz, "saida/limiares_publicados_padronizados.csv"), encoding = "UTF-8")
rf <- read.csv(file.path(raiz, "entrada/limiares_referencias.csv"), encoding = "UTF-8", check.names = FALSE)
cr <- readRDS(file.path(raiz, "intermediario/limiares_crossref.rds"))
en <- read.csv(file.path(raiz, "saida/limiares_enriquecimento.csv"), encoding = "UTF-8") |> group_by(id) |>
  summarise(uso = first(uso), peso = first(peso), confirmado = any(confirmado_enriq), .groups = "drop")
USO <- c(validacao = "validação", apoio = "apoio", relato = "relato")
t0 <- lp |> left_join(rf, by = "referencia") |> left_join(cr, by = "doi") |> left_join(en, by = "id") |> arrange(referencia, id) |>
  transmute(Id = id, Estudo = referencia, `Título` = gsub("&amp;", "&", titulo), Revista = gsub("&amp;", "&", revista), `Ano (Crossref)` = ano_doi, DOI = doi, Compartimento = ifelse(compartimento == "bentos", "bentos", "pelágico"),
    `Grupo no estudo` = grupo_publicado, `Região` = regiao, `Variável no estudo` = variavel_publicada, `Variável no modelo` = ifelse(is.na(variavel_modelo), "fora do modelo", variavel_modelo),
    `Limiar mín` = limiar_min, `Limiar máx` = limiar_max, Unidade = unidade, `Tipo de limiar` = tipo_limiar, `Método` = metodo, `Medição` = medicao, `Nota da extração` = nota,
    `Origem do valor` = ifelse(status_extracao == "resumo_sem_valor", "resumo sem valor: extrair do texto", "resumo"), `Dados no OBIS (triagem automática)` = dados_no_obis,
    `Uso atual` = ifelse(is.na(uso), "–", USO[uso]), `Peso atual` = peso, `Confirmado por enriquecimento` = ifelse(is.na(confirmado), "–", ifelse(confirmado, "sim", "não")))
PI <- c("Valor confere?", "Limiar mín corrigido", "Limiar máx corrigido", "Dados no OBIS (conferido)", "Tipo corrigido", "Observações")
t1 <- t0; for (p in PI) t1[[p]] <- NA_character_
n <- nrow(t1); nc <- ncol(t1); ipi <- which(names(t1) %in% PI); tipos <- sort(unique(lp$tipo_limiar))
wb <- createWorkbook(); modifyBaseFont(wb, fontName = "Arial", fontSize = 10)
addWorksheet(wb, "Validação"); addWorksheet(wb, "Leia-me"); addWorksheet(wb, "Listas", visible = FALSE)
writeData(wb, "Listas", data.frame(confere = c("sim", "não", "corrigir", "sem acesso ao texto")[1:4]), startCol = 1); writeData(wb, "Listas", data.frame(obis = c("sim", "parcial", "não", "indeterminado")), startCol = 2); writeData(wb, "Listas", data.frame(tipo = tipos), startCol = 3)
writeData(wb, "Validação", t1, headerStyle = createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold", fgFill = "#E7E6E6", wrapText = TRUE, valign = "top", border = "bottom"))
x <- paste0("https://doi.org/", t0$DOI); names(x) <- t0$DOI; class(x) <- "hyperlink"; writeData(wb, "Validação", x, startCol = which(names(t1) == "DOI"), startRow = 2)
addStyle(wb, "Validação", createStyle(wrapText = TRUE, valign = "top"), rows = 2:(n + 1), cols = 1:nc, gridExpand = TRUE, stack = TRUE)
addStyle(wb, "Validação", createStyle(fgFill = "#FFF2CC", border = "TopBottomLeftRight", borderColour = "#BFBFBF"), rows = 2:(n + 1), cols = ipi, gridExpand = TRUE, stack = TRUE)
addStyle(wb, "Validação", createStyle(fgFill = "#FFD966", textDecoration = "bold", wrapText = TRUE, valign = "top", border = "bottom"), rows = 1, cols = ipi, stack = FALSE)
for (j in seq_len(n)[-1]) if (t1$Estudo[j] != t1$Estudo[j - 1]) addStyle(wb, "Validação", createStyle(border = "top", borderColour = "#7F7F7F"), rows = j + 1, cols = 1:nc, stack = TRUE)
larg <- c(Id = 7, Estudo = 20, `Título` = 42, Revista = 20, `Ano (Crossref)` = 8, DOI = 26, Compartimento = 11, `Grupo no estudo` = 20, `Região` = 24, `Variável no estudo` = 16, `Variável no modelo` = 11,
  `Limiar mín` = 8, `Limiar máx` = 8, Unidade = 7, `Tipo de limiar` = 15, `Método` = 24, `Medição` = 9, `Nota da extração` = 36, `Origem do valor` = 12, `Dados no OBIS (triagem automática)` = 13, `Uso atual` = 9, `Peso atual` = 6, `Confirmado por enriquecimento` = 12,
  `Valor confere?` = 11, `Limiar mín corrigido` = 10, `Limiar máx corrigido` = 10, `Dados no OBIS (conferido)` = 13, `Tipo corrigido` = 15, `Observações` = 60)
setColWidths(wb, "Validação", cols = 1:nc, widths = larg[names(t1)]); setRowHeights(wb, "Validação", rows = 1, heights = 54)
freezePane(wb, "Validação", firstActiveRow = 2, firstActiveCol = 3); addFilter(wb, "Validação", rows = 1, cols = 1:nc)
dv <- function(col, ref) dataValidation(wb, "Validação", cols = which(names(t1) == col), rows = 2:(n + 1), type = "list", value = ref, showErrorMsg = FALSE)
dv("Valor confere?", "'Listas'!$A$2:$A$5"); dv("Dados no OBIS (conferido)", "'Listas'!$B$2:$B$5"); dv("Tipo corrigido", sprintf("'Listas'!$C$2:$C$%d", length(tipos) + 1))
tx <- c("Validação dos limiares publicados (base v0)", "",
  "Como preencher: confira cada linha no texto completo e use só as colunas amarelas, no fim da aba Validação (role para a direita ou filtre por Estudo). Não altere as demais colunas nem os Ids, porque a tabela volta para a análise pelo Id.",
  "Valor confere?: sim, não, corrigir (preencha os limiares corrigidos) ou sem acesso ao texto.",
  "Limiar mín e máx corrigidos: só quando o valor mudar, na unidade da linha. Para um valor único, repita o número nas duas colunas.",
  "Dados no OBIS (conferido): sim, parcial, não ou indeterminado. A resposta define se o limiar fica como apoio ou passa para validação (Gate 10, item 5).",
  "Tipo corrigido: só quando o tipo estiver errado. Somente troca_composicao, agrupamento e limite_distribuicao validam (Gate 10, item 2).",
  "Observações: texto livre com os pontos importantes do artigo.", "",
  "Exemplo de preenchimento (fictício):",
  "Id LP999 | Valor confere? = corrigir | Limiar mín corrigido = 1200 | Limiar máx corrigido = 1500 | Dados no OBIS (conferido) = não | Tipo corrigido = (vazio) | Observações = Tabela 2 do artigo; o resumo arredonda o intervalo; amostras só de verão.", "",
  "Sobre os dados da planilha:",
  "Todos os valores vêm dos resumos, por isso nada aqui foi conferido no texto completo.",
  "Os DOIs foram localizados pelo autor e pelo tema e conferidos no Crossref (título, revista e ano). Confirme na leitura que cada DOI é o estudo certo, em especial Carter et al. 2025, Matthews et al. 2023 e Tao et al. 2025, que foram localizados só pelo tema ou têm artigo próximo do mesmo grupo.",
  "Ano (Crossref) difere do ano usado na base em Olabarria (2004 × 2005), Swanborn (2024 × 2025), Ellis (2011 × 2012) e Yang (2017 × 2016), em geral por publicação on-line anterior ao volume.",
  "Limiares que não entram na figura do relatório (variável fora do modelo ou sem valor) também estão aqui. LP026, LP028 e LP033 estão sem valor: LP026 e LP033 precisam de extração do texto.")
writeData(wb, "Leia-me", data.frame(tx), colNames = FALSE); setColWidths(wb, "Leia-me", cols = 1, widths = 120)
addStyle(wb, "Leia-me", createStyle(wrapText = TRUE, valign = "top"), rows = seq_along(tx), cols = 1, stack = TRUE)
addStyle(wb, "Leia-me", createStyle(textDecoration = "bold", fontSize = 12), rows = 1, cols = 1, stack = TRUE); addStyle(wb, "Leia-me", createStyle(textDecoration = "bold"), rows = c(10, 13), cols = 1, stack = TRUE)
saveWorkbook(wb, arq); cat("gravada:", arq, "| linhas:", n, "| estudos:", n_distinct(t1$Estudo), "| sem DOI:", sum(is.na(t0$DOI)), "| sem titulo:", sum(is.na(t0$`Título`)), "\n")
