#' Generalised Canonical Procrustes
#'
#' Extract configuration parameters of gcproc
#'
#' @param k_dim Dimension reduction for samples (assumed to be along rows)
#' @param j_dim Dimension reduction for features (assumed to be along columns)
#' @param eta Step size of variational inference method. Please ensure this is appropriate (default = 5e-3)
#' @param max_iter Maximum iteration of gcproc
#' @param min_iter Minimum iteration of gcproc
#' @param tol Tolerance threshold for convergence (metric: Root Mean Squared Error)
#' @param log log with pseudo-count tranformation of data
#' @param center Center data by setting mean location to zero (akin to Procrustes)
#' @param scale.z Scale data by normalising (akin to Procrustes)
#' @param batches Number of total batches per run. Please ensure this is appropriate for your dataset
#' @param cores Number of cores the batches are split over each run
#' @param verbose Print statements?
#' @param init Initialisation method for the model ("svd","eigen-sparse","eigen-dense")
#' @return  Configuration parameters for gcproc
#' @export
extract_config <- function(method=NULL,verbose=T){
  if (method=="prediction"){
    config <- list(
      k_dim = 70,
      j_dim = 70,
      eta=1e-2,
      max_iter=15,
      min_iter = 5,
      tol=1e-3,
      log=F,
      center=F,
      scale.z=F,
      batches=2,
      batch_size=700,
      cores=2,
      verbose=T,
      init="svd-quick")
  }
  if (method=="clustering"){
    config <- list(
      k_dim = 70,
      j_dim = 70,
      eta=1e-2,
      max_iter=150,
      min_iter = 15,
      tol=1e-3,
      log=F,
      center=F,
      scale.z=F,
      batches=64,
      batch_size=64,
      cores=2,
      verbose=T,
      init="svd-quick")
  }


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
#' @param anchor_y.sample Equivalent to alpha.L.K
#' @param anchor_y.feature Equivalent to v.beta
#' @param anchor_x.sample Equivalent to alpha.L.J
#' @param anchor_x.feature Equivalent to u.beta
#' @param anchor_y.cov.sample Equivalent to y.gamma
#' @param anchor_y.cov.feature Equivalent to y.delta
#' @param anchor_x.cov.sample Equivalent to x.gamma
#' @param anchor_x.cov.feature Equivalent to x.delta
#' @return  Anchor framework for gcproc
#' @export
extract_anchors_framework <- function(verbose=T){
  anchors <- list(
    anchor_y.sample = NULL,
    anchor_y.feature = NULL,
    anchor_x.sample = NULL,
    anchor_x.feature = NULL,
    anchor_y.cov.sample = NULL,
    anchor_y.cov.feature = NULL,
    anchor_x.cov.sample = NULL,
    anchor_x.cov.feature = NULL
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
#' @param pivot_y.sample Equivalent to alpha.L.K
#' @param pivot_y.feature Equivalent to v.beta
#' @param pivot_x.sample Equivalent to alpha.L.J
#' @param pivot_x.feature Equivalent to u.beta
#' @param pivot_y.cov.sample Equivalent to y.gamma
#' @param pivot_y.cov.feature Equivalent to y.delta
#' @param pivot_x.cov.sample Equivalent to x.gamma
#' @param pivot_x.cov.feature Equivalent to x.delta
#' @return  Pivot framework for gcproc
#' @export
extract_pivots_framework <- function(verbose=T){
  pivots <- list(
    pivot_y.sample = NULL,
    pivot_y.feature = NULL,
    pivot_x.sample = NULL,
    pivot_x.feature = NULL,
    pivot_y.cov.sample = NULL,
    pivot_y.cov.feature = NULL,
    pivot_x.cov.sample = NULL,
    pivot_x.cov.feature = NULL
  )

  if (verbose == T){
    print(pivots)
  }

  return(pivots)
}


#' Extract covariates framework to put into gcproc
#'
#' Covariates can be used to further improve batch modelling.
#' NOTE: This is an empty framework that the user must fill in.
#'
#' @param covariates_y.sample Covariates with dimension equivalent to (p1 by samples of y)
#' @param covariates_y.feature Covariates with dimension equivalent to (features of y by p2)
#' @param covariates_x.sample Covariates with dimension equivalent to (p3 by samples of x)
#' @param covariates_x.feature Covariates with dimension equivalent to (features of x by p4)
#' @return  Covariate framework for gcproc
#' @export
extract_covariates_framework <- function(verbose=T){
  covariates <- list(
    covariates_y.sample = NULL,
    covariates_y.feature = NULL,
    covariates_x.sample = NULL,
    covariates_x.feature = NULL
  )

  if (verbose == T){
    print(covariates)
  }

  return(covariates)
}











#' Extract prediction framework to put into gcproc
#'
#' Predictions can made be on data not currently calculated
#' NOTE: This is an empty framework that the user must fill in.
#'
#' @param covariates_y.sample Covariates with dimension equivalent to (p1 by samples of y)
#' @param covariates_y.feature Covariates with dimension equivalent to (features of y by p2)
#' @param covariates_x.sample Covariates with dimension equivalent to (p3 by samples of x)
#' @param covariates_x.feature Covariates with dimension equivalent to (features of x by p4)
#' @return  Covariate framework for gcproc
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









