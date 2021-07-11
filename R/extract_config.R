#' Generalised Canonical Procrustes
#'
#' Runs gcproc with multiple initial seeds
#'
#' @param k_dim Dimension reduction for samples
#' @param j_dim Dimension reduction for features
#' @param eta Step size of variational inference method. Please ensure this is appropriate (default = 5e-3)
#' @param max_iter Maximum iteration of gcproc
#' @param min_iter Minimum iteration of gcproc
#' @param tol Tolerance threshold
#' @param log log with pseudo-count tranformation of data
#' @param center Center data by setting mean location to zero (akin to Procrustes)
#' @param scale.z Scale data by normalising (akin to Procrustes)
#' @param batches Number of total batches per run. Please ensure this is appropriate for your dataset
#' @param cores Number of cores the batches are split over each run
#' @param verbose Print statements?
#' @param init Initialisation method for the model ("svd","eigen-sparse","eigen-dense")
#' @return  The best seeded gcproc model with the highest negative log-likelihood is returned.
#' @export
extract_config <- function(verbose=T){
  config <- list(
       k_dim = 70,
       j_dim = 70,
       eta=5e-3,
       max_iter=1500,
       min_iter = 15,
       tol=1e-3,
       log=F,
       center=T,
       scale.z=T,
       batches=16,
       cores=2,
       verbose=T,
       init="svd")
  
  if (verbose == T){
    print(config)
  }
  
  return(config)
}