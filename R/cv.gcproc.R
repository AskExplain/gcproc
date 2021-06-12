cv.gcproc <- function(x,
                      y,
                      eta=1e-2,
                      max_iter= 1500,
                      min_iter = 15,
                      tol=1e-2,
                      log=F,
                      center=F,
                      scale.z=F,
                      batches=16,
                      cores=2,
                      verbose=F,
                      init="eigen"){

  main_llik <- c()
  for (j in seq(2,min(min(dim(x)[2],dim(y)[2]),30),length.out=7)){
    for (k in seq(2,min(min(dim(x)[1],dim(y)[1]),30),length.out=7)){

      j <- floor(j)
      k <- floor(k)

      if (verbose==T){
        print("Beginning tuning dimension to: ")
        print(paste("k_dim: ",k,"   j_dim: ",j,sep=""))
      }

      gcproc.model <- try(gcproc(y = y,
                                 x = x,
                                 k_dim = k,
                                 j_dim = j,
                                 eta = eta,
                                 max_iter = 55,
                                 min_iter = 55,
                                 tol = tol,
                                 log = log,
                                 center = center,
                                 scale.z = scale.z,
                                 batches = batches,
                                 cores = cores,
                                 verbose = T,
                                 init="random"),silent = F)
      if (!is.character(gcproc.model)){

        main_score.1 <- gcproc.model$main.parameters$alpha.L.J%*%rowSums(gcproc.model$transformed.data$x)
        main_score.2 <- gcproc.model$main.parameters$alpha.L.K%*%rowSums(gcproc.model$transformed.data$y)
        main_score.3 <- colSums(gcproc.model$transformed.data$x)%*%gcproc.model$main.parameters$u.beta
        main_score.4 <- colSums(gcproc.model$transformed.data$y)%*%gcproc.model$main.parameters$v.beta

        main_llik <- rbind(main_llik,c(j,k,sum(sqrt((main_score.1-main_score.2)^2))+sum(sqrt((main_score.3-main_score.4)^2))))
      } else {
        main_llik <- rbind(main_llik,c(Inf))
      }


    }
  }

  print(main_llik)
  main_dim <- main_llik[which(main_llik[,3]==min(main_llik[,3])),-3]
  main_j_dim <- main_dim[1]
  main_k_dim <- main_dim[2]

  if (verbose==T){
    print("Running final gcproc at optimal dimension of: ")
    print(paste("k_dim: ",main_k_dim,"   j_dim: ",main_j_dim,sep=""))
  }

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

  return(final.gcproc.model)
}
