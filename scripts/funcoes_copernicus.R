# funcoes_copernicus.R
# Loads the Copernicus Marine credentials from the local file login_cop.txt
# (line 1 = user, line 2 = password) into the R session options, without printing the values.
# The credentials file must never be copied to the Drive or included in deliverables.
cop_carregar_login <- function(arq = Sys.getenv("WMH_COPERNICUS_LOGIN", "C:/Users/model/OBIS_Tot/login_cop.txt")) {
  l <- trimws(readLines(arq, warn = FALSE)); l <- l[nzchar(l)]
  stopifnot(length(l) >= 2)
  options(CopernicusMarine_uid = l[1], CopernicusMarine_pwd = l[2])
  invisible(TRUE)
}
