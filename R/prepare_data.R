prepare_data <- function(x, y, log = F, center = F, scale.z = F){

  n = dim(x)[1]

  x <- as.matrix(x)
  y <- as.matrix(y)

  # Prepare for log_e(1+x) transform and library normalisation
  cx <- Matrix::colSums(x)
  cy <- Matrix::colSums(y)

  if (log == T){
    x <- log(1+x)
    y <- log(1+y)
  }
  if (center == T){

    x <- scale(x, center = TRUE, scale = FALSE)
    y <- scale(y, center = TRUE, scale = FALSE)

  }
  if (scale.z == T){

    x.size <- norm(x, type = "F") / (ncol(x) * nrow(x))
    x <- x / x.size

    y.size <- norm(y, type = "F") / (ncol(y) * nrow(y))
    y <- y / y.size

  }

  return(list(x=x,y=y))
}
