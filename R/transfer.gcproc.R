transfer.gcproc <- function(y,
                            x,
                            config = NULL,
                            initial_starts = 3,
                            anchors = NULL,
                            pivots = NULL
){

  y <- as.matrix(y)
  x <- as.matrix(x)

  if (!is.null(anchors)){
    config$anchors = anchors
  } else {
    print("Anchors must be provided with format:")
    print(
    "anchors <- list( anchor_y.sample = NULL,
                      anchor_y.feature = NULL,
                      anchor_x.sample = NULL,
                      anchor_x.feature = NULL  )"
          )
  }

  internal_config <- config
  internal_config$max_iter <- 15

  main_llik <- c()
  for (seed in c(1:initial_starts)){
    print(paste("seed: ",seed))
    set.seed(seed)
    final.gcproc.model <- try(gcproc(x = x,
                                     y = y,
                                     config = internal_config,
                                     seed = seed,
                                     anchors = anchors,
                                     pivots = pivots),silent = F)

    if (!is.character(final.gcproc.model)){
      main_llik <- rbind(main_llik,c(seed,tail(final.gcproc.model$convergence.parameters$llik.vec,1)))
    } else {
      main_llik <- rbind(main_llik,c(-Inf))
    }

  }

  main_dim <- main_llik[which(main_llik[,2]==max(main_llik[,2])),]
  main_seed <- main_dim[1]


  if (config$verbose){
    print("Running final gcproc at optimal dimension of: ")
    print(paste("seed:", main_seed,sep=""))
  }

  set.seed(main_seed)
  final.gcproc.model <- gcproc(x = x,
                               y = y,
                               config = config,
                               seed = main_seed,
                               anchors = anchors,
                               pivots = pivots)

  final.initial_starts <- initial_starts
  return(final.gcproc.model)

}
