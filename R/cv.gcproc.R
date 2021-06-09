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
  for (l in seq(2,30,length.out=5)){
    for (k in seq(2,30,length.out=5)){


      if (verbose==T){
        print("Beginning tuning dimension to: ")
        print(paste("k_dim: ",k,"   l_dim: ",l,sep=""))
      }

      gcproc.model <- try(gcproc(y = y,
                                 x = x,
                                 k_dim = k,
                                 l_dim = l,
                                 eta = eta,
                                 max_iter = 35,
                                 min_iter = 35,
                                 tol = tol,
                                 log = F,
                                 center = T,
                                 scale.z = F,
                                 batches = batches,
                                 cores = cores,
                                 verbose = F,
                                 init="random"),silent = F)
      if (!is.character(gcproc.model)){

        main_score.1 <- gcproc.model$main.parameters$alpha.L.J%*%rowSums(gcproc.model$transformed.data$x)
        main_score.2 <- gcproc.model$main.parameters$alpha.L.K%*%rowSums(gcproc.model$transformed.data$y)
        main_score.3 <- colSums(gcproc.model$transformed.data$x)%*%gcproc.model$main.parameters$u.beta
        main_score.4 <- colSums(gcproc.model$transformed.data$y)%*%gcproc.model$main.parameters$v.beta

        main_llik <- rbind(main_llik,c(l,k,sum(sqrt((main_score.1-main_score.2)^2))+sum(sqrt((main_score.3-main_score.4)^2))))
      } else {
        break
      }


    }
  }

  print(main_llik)
  main_dim <- main_llik[which(main_llik[,3]==min(main_llik[,3])),-3]
  main_l_dim <- main_dim[1]
  main_k_dim <- main_dim[2]

  if (verbose==T){
    print("Running final gcproc at optimal dimension of: ")
    print(paste("k_dim: ",main_k_dim,"   l_dim: ",main_l_dim,sep=""))
  }

  final.gcproc.model <- gcproc(x = x,
                         y = y,
                         k_dim = main_k_dim,
                         l_dim = main_l_dim,
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

  return(list(gcproc.model = final.gcproc.model,
              cv.llik = main_llik)
  )
}
