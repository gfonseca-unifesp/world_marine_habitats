# funcoes_mapa.R  (habitat class maps; used by 85_sintese.R and 89b_habitats_novos.R)
# Gate 8 (30/10/2026): stand-by is no longer a single gray; each stand-by habitat has its own color, with diagonal hatching
#   (dark stripes every 7 cells), so that potential habitats can be identified on the map.
PAL_HAB <- c(palette.colors(36, "Polychrome 36")[-c(1, 2)], palette.colors(26, "Alphabet"))   # without white/dark gray
COR_ACH <- "#000000A6"
mapa_classes <- function(r, cores, hach, rotulos, titulo, arq, ncol = 8) {
  # r: SpatRaster with values 1..length(cores); hach: logical per class (TRUE = stand-by, hatched)
  v <- values(r)[, 1]; rc <- rowColFromCell(r, seq_len(ncell(r)))
  hv <- rep(NA_integer_, ncell(r)); k <- which(!is.na(v)); k <- k[hach[v[k]] & (rc[k, 1] + rc[k, 2]) %% 7 < 2]; hv[k] <- 1L
  h <- rast(r); h[] <- hv
  nl <- ceiling(length(rotulos) / ncol)
  png(arq, 2000, 1000 + 45 * nl, res = 150, bg = "white")
  layout(matrix(1:2), heights = c(1000, 45 * nl + 40)); par(mar = c(0.5, 1, 2, 1))
  plot(r, col = cores, breaks = c(0.5, seq_along(cores) + 0.5), legend = FALSE, axes = FALSE, box = FALSE, main = titulo, maxcell = ncell(r))
  if (length(k)) plot(h, col = COR_ACH, breaks = c(0.5, 1.5), legend = FALSE, axes = FALSE, box = FALSE, add = TRUE, maxcell = ncell(r))
  par(mar = c(0, 1, 0, 1)); plot.new()
  legend("top", legend = rotulos, fill = cores, ncol = ncol, cex = 0.7, bty = "n", x.intersp = 0.5)
  if (any(hach)) legend("top", legend = rotulos, fill = ifelse(hach, COR_ACH, "transparent"), border = ifelse(hach, "black", "transparent"), density = ifelse(hach, 40, NA), angle = 45,
                        ncol = ncol, cex = 0.7, bty = "n", x.intersp = 0.5, text.col = NA)
  dev.off()
}
