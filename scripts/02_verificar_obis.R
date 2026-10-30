# 02_verificar_obis.R  (task 0.1, acceptance)
# Checks the local copy against the manifest: size of every file, and MD5 = ETag for files
# uploaded to S3 as a single part (ETag without hyphen). Multipart files (ETag with hyphen) are checked by size only.
# Output: bruto/verificacao.csv. After verification, the folder bruto/occurrence becomes read-only.
raiz <- Sys.getenv("WMH_ROOT", "C:/Users/model/OBIS_Tot")
man <- read.csv(file.path(raiz, "bruto/manifesto_s3.csv"), stringsAsFactors = FALSE)
man$local <- file.path(raiz, "bruto", man$key)
man$existe <- file.exists(man$local)
man$tamanho_ok <- man$existe & file.size(man$local) == man$size
man$multipart <- grepl("-", man$etag)
man$md5_local <- NA_character_
simples <- which(man$existe & !man$multipart)
man$md5_local[simples] <- unname(tools::md5sum(man$local[simples]))
man$md5_ok <- ifelse(man$multipart, NA, man$md5_local == man$etag)
man$status <- ifelse(!man$existe, "ausente", ifelse(!man$tamanho_ok, "tamanho_diferente",
                ifelse(man$multipart, "ok_tamanho", ifelse(man$md5_ok, "ok_md5", "md5_diferente"))))
write.csv(man[, c("key", "size", "etag", "md5_local", "status")], file.path(raiz, "bruto/verificacao.csv"), row.names = FALSE)
print(table(man$status))
if (all(man$status %in% c("ok_md5", "ok_tamanho"))) {
  Sys.chmod(man$local, mode = "0444")   # read-only
  cat("todos conferidos; arquivos marcados como somente leitura\n")
}
