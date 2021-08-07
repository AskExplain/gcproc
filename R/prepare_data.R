#' @export
prepare_data <- function(x, log = F, center = F, scale.z = F){

  n = dim(x)[1]

  x <- as.matrix(x)

  if (log == T){
    x <- log(1+x)
  }
  if (center == T){

    x <- scale(x, center = TRUE, scale = FALSE)

  }
  if (scale.z == T){

    x.size <- norm(x, type = "F") / (ncol(x) * nrow(x))
    x <- x / x.size

  }

  return(x)
}
