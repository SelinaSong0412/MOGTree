combine.mat <- function(a, b) {
  nrow1 <- nrow(a)
  ncol1 <- ncol(a)
  nrow2 <- nrow(b)
  ncol2 <- ncol(b)

  combine <- matrix(NA, max(nrow1, nrow2), ncol1 + ncol2)
  combine[1:nrow1, 1:ncol1] <- a
  combine[1:nrow2, (ncol1 + 1):(ncol1 + ncol2)] <- b
  combine
}
