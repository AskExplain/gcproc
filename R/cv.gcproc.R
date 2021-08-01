#' Generalised Canonical Procrustes
#'
#' Runs gcproc with multiple initial seeds
#'
#' @param x Reference dataset of cell by gene matrix
#' @param y Experimental dataset of cell by gene matrix
#' @param config Configuration parameters (please read gcproc code for more details)
#' @param initial_starts Multiple fixed seeds for initial seeding
#' @param anchors Fixing and anchoring the main model parameters to transfer prior information
#' @param pivots Initialisation of the main model parameters to speed up learning process
#' @return  The best seeded gcproc model with the highest negative log-likelihood is returned.
#' @export
cv.gcproc <- function(x,
                      y,
                      config = list(k_dim = 70,
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
                                             init="svd-quick"),
                      anchors = NULL,
                      pivots = NULL,
                      initial_starts = 3
){

  param.dim <- data.frame(c(1:initial_starts),config$k_dim,config$j_dim)


  internal_config <- config
  internal_config$max_iter <- 15

  main_llik <- c()
  for (row in 1:nrow(param.dim)){
    seed <- param.dim[row,1]
    set.seed(seed)
    k <- param.dim[row,2]
    j <- param.dim[row,3]

    if (config$verbose==T){
      print("Beginning tuning dimension to: ")
      print(paste("seed: ",seed, "   k_dim: ",k,"   j_dim: ",j,sep=""))
    }

    gcproc.model <- try(gcproc(y = y,
                               x = x,
                               config = internal_config,
                               seed = seed,
                               anchors = anchors,
                               pivots = pivots
                               ),silent = F)
    if (!is.character(gcproc.model)){
      cost <- tail(gcproc.model$convergence.parameters$llik.vec,1)
      main_llik <- rbind(main_llik,c(seed,k,j,cost))
    } else {
      main_llik <- rbind(main_llik,c(-Inf))
    }
  }

  main_dim <- main_llik[which(main_llik[,4]==max(main_llik[,4])),]
  main_seed <- main_dim[1]
  config$k_dim <- main_dim[2]
  config$j_dim <- main_dim[3]

  if (config$verbose==T){
    print("Running final gcproc at optimal dimension of: ")
    print(paste("seed:", main_seed, "   k_dim: ",config$k_dim,"   j_dim: ",config$j_dim,sep=""))
  }

  set.seed(main_seed)
  final.gcproc.model <- gcproc(x = x,
                               y = y,
                               config = config,
                               pivots = pivots,
                               anchors = anchors,
                               seed = main_seed)

  final.gcproc.model$cv.llik <- main_llik
  final.gcproc.model$meta.parameters$initial_starts <- initial_starts

  return(final.gcproc.model)
}
