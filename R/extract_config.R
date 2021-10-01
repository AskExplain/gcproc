#' Extract configuration parameters of gcproc
#'
#' @param i_dim Dimension reduction for samples (assumed to be along rows)
#' @param j_dim Dimension reduction for features (assumed to be along columns)
#' @param min_iter Minimum iteration of gcproc
#' @param max_iter Maximum iteration of gcproc
#' @param tol Tolerance threshold for convergence (metric: Root Mean Squared Error)
#' @param n_cores Number of CPU cores to use for prediction part only
#' @param verbose Print statements?
#' @param init Initialisation method for the model ("random","eigen-quick","eigen-dense","svd-quick","svd-dense")
#' @return  Configuration parameters for gcproc
#' @export
extract_config <- function(verbose=T){
  config <- list(
    i_dim = 30,
    j_dim = 30,
    min_iter=2,
    max_iter=350,
    tol=1,
    batch = 3,
    n_cores = 2,
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
#' @param task Allows user to specify either a regression or classification task
#' @param design.list A list of design structures where each element is given a 1 to indicate the test set, 0 indicates the train set.
#' @param labels For classification, these are the pre-defined labels
#' @param predict.list This will be filled in by gpcroc with the predictions and return a prediction for indicated design matrices only. Leave as NULL to begin.
#' @return  Prediction framework for gcproc
#' @export
extract_recovery_framework <- function(verbose=T){
  recover <- list(
    task = c("regression"),    # c("classification")
    method = c("knn.reg"),     # c("label.projection)
    design.list = NULL,
    labels = NULL,
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
#' @param alpha Fixing the alpha parameters. A vector of integers, where identical integers indicate same the data axis. Axes that are not shared are given NA.
#' @param beta Fixing the beta parameters. A vector of integers, where identical integers indicate same the data axis. Axes that are not shared are given NA.
#' @export
extract_fixed_framework <- function(verbose=T){
  fixed <- list(alpha=NULL,
                beta=NULL)

  if (verbose == T){
    print(fixed)
  }

  return(fixed)
}



#' @export
chunk <- function(x,n){
  if (n==1){
    list(x)
  }
  else{
    split(x, cut(seq_along(x), n, labels = FALSE))
  }
}
