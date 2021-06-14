cv.gcproc <- function(x,
                      y,
                      seeds = 3,
                      eta=1e-2,
                      max_iter= 1500,
                      min_iter = 15,
                      tol=1e-3,
                      log=F,
                      center=F,
                      scale.z=F,
                      batches=16,
                      cores=2,
                      verbose=F,
                      init="eigen"){

  main_llik <- c()
  for (seed in c(1:seeds)){
    for (j in seq(2,min(min(dim(x)[2],dim(y)[2])-1,30),length.out=3)){
      for (k in seq(2,min(min(dim(x)[1],dim(y)[1])-1,30),length.out=3)){
        set.seed(seed)
        j <- floor(j)
        k <- floor(k)

        if (verbose==T){
          print("Beginning tuning dimension to: ")
          print(paste("seed: ",seed, "   k_dim: ",k,"   j_dim: ",j,sep=""))
        }

        gcproc.model <- try(gcproc(y = y,
                                   x = x,
                                   k_dim = k,
                                   j_dim = j,
                                   eta = eta,
                                   max_iter = 105,
                                   min_iter = 15,
                                   tol = tol,
                                   log = log,
                                   center = center,
                                   scale.z = scale.z,
                                   batches = batches,
                                   cores = cores,
                                   verbose = F,
                                   init=init),silent = F)
        if (!is.character(gcproc.model)){
          main_llik <- rbind(main_llik,c(seed,j,k,tail(gcproc.model$convergence.parameters$llik.vec,1)))
        } else {
          main_llik <- rbind(main_llik,c(-Inf))
        }


      }
    }
  }

  main_dim <- main_llik[which(main_llik[,4]==max(main_llik[,4])),]
  main_seed <- main_dim[1]
  main_j_dim <- main_dim[2]
  main_k_dim <- main_dim[3]

  if (verbose==T){
    print("Running final gcproc at optimal dimension of: ")
    print(paste("seed:", main_seed, "   k_dim: ",main_k_dim,"   j_dim: ",main_j_dim,sep=""))
  }

  set.seed(main_seed)
  final.gcproc.model <- gcproc(x = x,
                         y = y,
                         k_dim = main_k_dim,
                         j_dim = main_j_dim,
                         eta = eta,
                         max_iter = max_iter,
                         min_iter = min_iter,
                         tol = tol,
                         log = log,
                         center = center,
                         scale.z = scale.z,
                         batches = batches,
                         cores = cores,
                         verbose = verbose,
                         init=init)


  final.gcproc.model$cv.llik <- main_llik
  final.gcproc.model$meta.parameters$seeds <- seeds

  return(final.gcproc.model)
}
