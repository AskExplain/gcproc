#' Extract configuration parameters of gcproc
#'
#' @param i_dim Dimension reduction for samples (assumed to be along rows)
#' @param j_dim Dimension reduction for features (assumed to be along columns)
#' @param min_iter Minimum iteration of gcproc
#' @param max_iter Maximum iteration of gcproc
#' @param tol Tolerance threshold for convergence (metric: Root Mean Squared Error)
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
    verbose=T,
    init="random")

  if (verbose == T){
    print(config)
  }

  return(config)
}

#' Extract anchor framework to put into gcproc
#'
#' Transfers learned parameters from a pre-trained model.
#' NOTE: This is an empty framework that the user must fill in.
#'
#' @param code Anchor the code
#' @return  Anchor framework for gcproc
#' @export
extract_transfer_framework <- function(verbose=T){
  transfer <- list(
    code = NULL,
    pivot = NULL
  )

  if (verbose == T){
    print(transfer)
  }

  return(transfer)
}


#' Extract recovery framework to put into gcproc
#'
#' Can recover data points by imputing or predicting missing values
#'
#' @param task Allows user to specify either a regression or classification task
#' @param method The algorithm for the task (Options are regression: "knn.reg","matrix.projection", -- provide your own --   ;   classification: "label.projection")
#' @param design.list A list of design structures where each element is given a 1 to indicate the test set, 0 indicates the train set.
#' @param labels For classification, these are the pre-defined labels
#' @return  Prediction framework for gcproc
#' @export
extract_recovery_framework <- function(verbose=T){
  recover <- list(
    task = c("regression"),    # c("classification")
    method = c("knn.reg"),     # c("label.projection)
    design.list = NULL,
    labels = NULL
  )

  if (verbose == T){
    print(recover)
  }

  return(recover)
}






#' Extract join framework to put into gcproc
#'
#' Join data to improve modelling capacity for similar axes
#' @param alpha Joining the alpha parameters. A vector of integers, where identical integers indicate same the data axis to be joined. Axes that should not be shared are given NA.
#' @param beta Joining the beta parameters. A vector of integers, where identical integers indicate same the data axis to be joined. Axes that should not be shared are given NA.
#' @export
extract_join_framework <- function(verbose=T){
  join <- list(alpha=NULL,
                beta=NULL)

  if (verbose == T){
    print(join)
  }

  return(join)
}
