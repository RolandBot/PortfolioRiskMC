ES99 <- function(L) {
  es99 <- apply(L, 2, function(l) {
    mean(l[tail99(l)])
  })
  es99
}

tail99 <- function(l, sorted = TRUE) {
  idx <- tail(order(l), ceiling(0.01*length(l)))
  if (sorted) {
    idx <- sort(idx)
  }
  idx
}
