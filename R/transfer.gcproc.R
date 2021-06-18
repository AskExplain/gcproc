transfer.gcproc <- function(gcproc.model,y,x,anchors=NULL,initial.param=NULL){

  y <- as.matrix(y)
  x <- as.matrix(x)

  k_dim <- gcproc.model$meta.parameters$k_dim
  j_dim <- gcproc.model$meta.parameters$j_dim

  if (is.null(anchors)){
    initial.anchors <- initialise.gcproc(x=x,y=y,init=init,k_dim=k_dim,j_dim=j_dim)

    anchor_y.sample = NULL
    anchor_y.feature = NULL
    anchor_x.sample = NULL
    anchor_x.feature = NULL

    if (is.null(anchors)){
      if (dim(y)[1]==dim(gcproc.model$transformed.data$y)[1]){
        anchor_y.sample = gcproc.model$main.parameters$alpha.L.K
      }
      if (dim(y)[2]==dim(gcproc.model$transformed.data$y)[2]){
        anchor_y.feature = gcproc.model$main.parameters$v.beta
      }
      if (dim(x)[1]==dim(gcproc.model$transformed.data$x)[1]){
        anchor_x.sample = gcproc.model$main.parameters$alpha.L.J
      }
      if (dim(x)[2]==dim(gcproc.model$transformed.data$x)[2]){
        anchor_x.feature = gcproc.model$main.parameters$u.beta
      }

      anchors$anchor_y.sample <- anchor_y.sample
      anchors$anchor_y.feature <- anchor_y.feature
      anchors$anchor_x.sample <- anchor_x.sample
      anchors$anchor_x.feature <- anchor_x.feature

    }
  }

  if (is.null(initial.param)){
    initial.param <- initialise.gcproc(x=x,y=y,init="svd",k_dim=70,j_dim=70)
  }

  main_llik <- c()
  for (seed in c(1:gcproc.model$meta.parameters$seeds)){
    print(paste("seed: ",seed))
    set.seed(seed)
    final.gcproc.model <- try(gcproc(x = x,
                                     y = y,
                                     k_dim = k_dim,
                                     j_dim = j_dim,
                                     eta = gcproc.model$meta.parameters$eta,
                                     max_iter = 15,
                                     min_iter = gcproc.model$meta.parameters$min_iter,
                                     tol = gcproc.model$meta.parameters$tol,
                                     batches = gcproc.model$meta.parameters$batches,
                                     cores = gcproc.model$meta.parameters$cores,
                                     verbose = T,
                                     init = "svd",
                                     log = gcproc.model$meta.parameters$log,
                                     center = gcproc.model$meta.parameters$center,
                                     scale.z = gcproc.model$meta.parameters$scale.z,
                                     anchors = anchors,
                                     initial.param = initial.param,
                                     seed = seed),silent = F)

    if (!is.character(final.gcproc.model)){
      main_llik <- rbind(main_llik,c(seed,tail(final.gcproc.model$convergence.parameters$llik.vec,1)))
    } else {
      main_llik <- rbind(main_llik,c(-Inf))
    }

  }

  main_dim <- main_llik[which(main_llik[,2]==max(main_llik[,2])),]
  main_seed <- main_dim[1]

  set.seed(main_seed)
  final.gcproc.model <- gcproc(x = x,
                               y = y,
                               k_dim = gcproc.model$meta.parameters$k_dim,
                               j_dim = gcproc.model$meta.parameters$j_dim,
                               eta = gcproc.model$meta.parameters$eta,
                               max_iter = gcproc.model$meta.parameters$max_iter,
                               min_iter = gcproc.model$meta.parameters$min_iter,
                               tol = gcproc.model$meta.parameters$tol,
                               batches = gcproc.model$meta.parameters$batches,
                               cores = gcproc.model$meta.parameters$cores,
                               verbose = gcproc.model$meta.parameters$verbose,
                               init = "svd",
                               log = gcproc.model$meta.parameters$log,
                               center = gcproc.model$meta.parameters$center,
                               scale.z = gcproc.model$meta.parameters$scale.z,
                               anchors = anchors,
                               initial.param = initial.param,
                               seed = main_seed)



  final.gcproc.model$meta.parameters$seeds <- gcproc.model$meta.parameters$seeds
  return(final.gcproc.model)

}
