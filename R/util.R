colQuants <- function (x, q) apply(x, 2, \(col) quantile(col, probs = q))
rowQuants <- function (x, q) apply(x, 1, \(row) quantile(row, probs = q))
