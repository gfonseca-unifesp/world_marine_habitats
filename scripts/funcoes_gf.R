# funcoes_gf.R  (Phase 8, stage 3: shared Gradient Forest functions; used by 88_gf_gdm.R and 88c_transferencia.R)
# limiares_curva: thresholds = peaks of the rate of change of the cumulative importance curve (5-point moving average on a grid of
#   200 points) with density >= 2x the gradient mean; at most 3 per predictor. ganho = share of the importance in the segment.
limiares_curva <- function(x, y, n = 200, max_pk = 3) {
  xs <- seq(min(x), max(x), length.out = n); ys <- approx(x, y, xs, ties = mean)$y
  dy <- diff(ys) / diff(xs); dy <- stats::filter(dy, rep(1 / 5, 5), sides = 2); dy[is.na(dy)] <- 0
  media <- (max(ys) - min(ys)) / (max(xs) - min(xs)); if (!is.finite(media) || media <= 0) return(NULL)
  pk <- which(diff(sign(diff(dy))) == -2) + 1; pk <- pk[dy[pk] >= 2 * media]
  pk <- head(pk[order(-dy[pk])], max_pk)
  if (!length(pk)) return(NULL)
  data.frame(limiar = (xs[pk] + xs[pk + 1]) / 2, intensidade = as.numeric(dy[pk]) / media,
    ganho = vapply(pk, function(k) { j <- max(1, k - 5):min(length(ys), k + 6); diff(range(ys[j])) / max(ys) }, 0))
}
