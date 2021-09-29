#' Extract configuration parameters of gcproc
#'
#' @param i_dim Dimension reduction for samples (assumed to be along rows)
#' @param j_dim Dimension reduction for features (assumed to be along columns)
#' @param max_iter Maximum iteration of gcproc
#' @param tol Tolerance threshold for convergence (metric: Root Mean Squared Error)
#' @param verbose Print statements?
#' @param init Initialisation method for the model ("random","eigen-quick","eigen-dense","svd-quick","svd-dense")
#' @return  Configuration parameters for gcproc
#' @export
extract_config <- function(verbose=T){
  config <- list(
    i_dim = 100,
    j_dim = 100,
    k_dim = 100,
    min_iter=2,
    max_iter=350,
    tol=1,
    verbose=T,
    init="random")

  if (verbose == T){
    print(config)
  }

  return(config)
}

#' Extract anchor framework to put into gcproc
#'
#' Anchors allow the transfer of learned parameters from a pre-trained model.
#' NOTE: This is an empty framework that the user must fill in.
#'
#' @param code Anchor the code
#' @return  Anchor framework for gcproc
#' @export
extract_anchors_framework <- function(verbose=T){
  anchors <- list(
    code = NULL
  )

  if (verbose == T){
    print(anchors)
  }

  return(anchors)
}

#' Extract pivot framework to put into gcproc.
#'
#' Pivots allow initialisation of parameters as input.
#' NOTE: This is an empty framework that the user must fill in.
#'
#' @param code code
#' @return  Pivot framework for gcproc
#' @export
extract_pivots_framework <- function(verbose=T){
  pivots <- list(
    code = NULL
  )

  if (verbose == T){
    print(pivots)
  }

  return(pivots)
}



#' Extract recovery framework to put into gcproc
#'
#' Can recover data points by imputing or predicting missing values
#'
#' @param x Design matrix of x where 1 is to predict the test set, 0 is to be modelled as the train set
#' @param y Design matrix of y where 1 is to predict the test set, 0 is to be modelled as the train set
#' @param fn Allows user to specify a prediction function
#' @param param Parameters to be put into prediction function
#' @return  Prediction framework for gcproc
#' @export
extract_recovery_framework <- function(verbose=T){
  recover <- list(
    method = c("matrix.projection","knn"),
    x = NULL,
    y = NULL,
    fn = NULL,
    param = NULL,
    design.list = NULL,
    predict.list = NULL
  )

  if (verbose == T){
    print(recover)
  }

  return(recover)
}






#' Extract fixed framework to put into gcproc
#'
#' Fix data to improve modelling capacity for similar axes
#' @export
extract_fixed_framework <- function(verbose=T){
  fixed <- list(alpha=NULL,
                beta=NULL)

  if (verbose == T){
    print(fixed)
  }

  return(fixed)
}

transform.data <- function(x,method="scale"){

  if (method == "scale"){
    center = T
    scale = T

    x <- as.matrix(x)
    nc <- ncol(x)
    if (is.logical(center)) {
      if (center) {
        center <- colMeans(x, na.rm=TRUE)
        x <- sweep(x, 2L, center, check.margin=FALSE)
      }
    }
    else if (is.numeric(center) && (length(center) == nc))
      x <- sweep(x, 2L, center, check.margin=FALSE)
    else
      stop("length of 'center' must equal the number of columns of 'x'")
    if (is.logical(scale)) {
      if (scale) {
        f <- function(v) {
          v <- v[!is.na(v)]
          sqrt(sum(v^2) / max(1, length(v) - 1L))
        }
        scale <- apply(x, 2L, f)
        scale <- sapply(scale,function(scale){if(scale==0|is.na(scale)){1}else{scale}})
        x <- sweep(x, 2L, scale, "/", check.margin=FALSE)
      }
    }
    else if (is.numeric(scale) && length(scale) == nc)
      x <- sweep(x, 2L, scale, "/", check.margin=FALSE)
    else
      stop("length of 'scale' must equal the number of columns of 'x'")
    if(is.numeric(center)) attr(x, "scaled:center") <- center
    if(is.numeric(scale)) attr(x, "scaled:scale") <- scale
    x
  }

  return(x)
}
