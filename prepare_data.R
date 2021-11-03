#' @export
prepare_data <- function(x, log = F, center = F, scale.norm = F){

  n = dim(x)[1]

  x <- as.matrix(x)

  # Prepare for log_e(1+x) transform and library normalisation
  cx <- Matrix::colSums(x)

  if (log == T){
    x <- log(1+x)
  }
  if (center == T){

    x <- scale(x, center = TRUE, scale = FALSE)

  }
  if (scale.norm == T){

    x.size <- norm(x, type = "F") / (ncol(x) * nrow(x))
    x <- x / x.size


  }

  return(x)
}
