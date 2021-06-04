# To cut data into batches
chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
