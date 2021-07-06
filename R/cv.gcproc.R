cv.gcproc <- function(x,
                      y,
                      config = NULL,
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
