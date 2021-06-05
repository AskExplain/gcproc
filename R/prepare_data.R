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
  if (center == T | scale.z == T){

    x <- scale(x,center = center,scale = scale.z)
    y <- scale(y,center = center,scale = scale.z)

  }

  return(list(x=x,y=y))
}
