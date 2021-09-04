#' Extract configuration parameters of gcproc
#'
#' @param i_dim Dimension reduction for samples (assumed to be along rows)
#' @param j_dim Dimension reduction for features (assumed to be along columns)
#' @param max_iter Maximum iteration of gcproc
#' @param tol Tolerance threshold for convergence (metric: Root Mean Squared Error)
#' @param verbose Print statements?
#' @param init Initialisation method for the model ("random","eigen-sparse","eigen-dense","svd-quick","svd-dense)
#' @return  Configuration parameters for gcproc
#' @export
extract_config <- function(verbose=T){
  config <- list(
    all.i_dim = c(30,25),
    all.j_dim = c(30,25),
    layers = 2,
    max_iter=150,
    tol=1e-2,
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
    method = c("matrix.projection"),
    x = NULL,
    y = NULL,
    fn = NULL,
    param = NULL
  )

  if (verbose == T){
    print(recover)
  }

  return(recover)
}







