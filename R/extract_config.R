#' Extract configuration parameters of gcproc
#'
#' @param i_dim Dimension reduction for samples (assumed to be along rows)
#' @param j_dim Dimension reduction for features (assumed to be along columns)
#' @param max_iter Maximum iteration of gcproc
#' @param tol Tolerance threshold for convergence (metric: Root Mean Squared Error)
#' @param verbose Print statements?
#' @param init Initialisation method for the model ("svd","eigen-sparse","eigen-dense")
#' @return  Configuration parameters for gcproc
#' @export
extract_config <- function(verbose=T){
  config <- list(
    i_dim = 30,
    j_dim = NULL,
    max_iter=150,
    tol=1e-5,
    verbose=T,
    init="svd-quick")

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
#' @param anchor_y.sample Equivalent to alpha.K
#' @param anchor_y.feature Equivalent to v.beta
#' @param anchor_x.sample Equivalent to alpha.L
#' @param anchor_x.feature Equivalent to u.beta
#' @return  Anchor framework for gcproc
#' @export
extract_anchors_framework <- function(verbose=T){
  anchors <- list(
    anchor_y.sample = NULL,
    anchor_y.feature = NULL,
    anchor_x.sample = NULL,
    anchor_x.feature = NULL
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
#' @param pivot_y.sample Equivalent to alpha.K
#' @param pivot_y.feature Equivalent to v.beta
#' @param pivot_x.sample Equivalent to alpha.L
#' @param pivot_x.feature Equivalent to u.beta
#' @return  Pivot framework for gcproc
#' @export
extract_pivots_framework <- function(verbose=T){
  pivots <- list(
    pivot_y.sample = NULL,
    pivot_y.feature = NULL,
    pivot_x.sample = NULL,
    pivot_x.feature = NULL
  )

  if (verbose == T){
    print(pivots)
  }

  return(pivots)
}



#' Extract prediction framework to put into gcproc
#'
#' Can either impute or predict by replacing missing values
#'
#' @param x Design matrix of x where 1 is to predict the test set, 0 is to be modelled as the train set
#' @param y Design matrix of y where 1 is to predict the test set, 0 is to be modelled as the train set
#' @return  Prediction framework for gcproc
#' @export
extract_prediction_framework <- function(verbose=T){
  prediction <- list(
    x = NULL,
    y = NULL
  )

  if (verbose == T){
    print(prediction)
  }

  return(prediction)
}







