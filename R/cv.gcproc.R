cv.gcproc <- function(x,
                      y,
                      k_dim = 30,
                      l_dim = 30,
                      eta=1e-2,
                      max_iter= 1500,
                      min_iter = 15,
                      tol=1e-3,
                      log=F,
                      center=F,
                      scale.z=F,
                      batches=2,
                      cores=2,
                      verbose=F,
                      init="eigen"){

  min_dim <- min(min(dim(x)),min(dim(y)))-2

  if (verbose==T){
    print("Running initial maximum dimensions")
  }
  gcproc.model <- gcproc(x = x,
         y = y,
         k_dim = min_dim,
         l_dim = min_dim,
         eta = eta,
         max_iter = 1,
         min_iter = 1,
         tol = tol,
         log = log,
         center = center,
         scale.z = scale.z,
         batches = batches,
         cores = cores,
         verbose = F,
         init="eigen")


  if (verbose==T){
    print("Beginning tuning dimensions")
  }

  main_llik <- c()
  for (i in c(seq((min_dim-2),2,length.out=min(30,(min_dim-2)/2)))){

    i <- round(i,0)

    if (verbose==T){
      print("Beginning tuning dimension to: ")
      print(i)
    }

    transformed.y <- gcproc.model$main.parameters$alpha.L.K%*%gcproc.model$transformed.data$y%*%gcproc.model$main.parameters$v.beta
    transformed.x <- gcproc.model$main.parameters$alpha.L.J%*%gcproc.model$transformed.data$x%*%gcproc.model$main.parameters$u.beta

    gcproc.model <- try(gcproc(y = transformed.y,
                           x = transformed.x,
                           k_dim = i,
                           l_dim = i,
                           eta = eta,
                           max_iter = 1,
                           min_iter = 1,
                           tol = tol,
                           log = F,
                           center = F,
                           scale.z = F,
                           batches = batches,
                           cores = cores,
                           verbose = F,
                           init="eigen"),silent = F)
    if (!is.character(gcproc.model)){
      main_llik <- rbind(main_llik,c(i,tail(gcproc.model$convergence.parameters$llik.vec,1)))
    } else {
      break
    }

  }

  main_dim <- main_llik[which(main_llik[,2]==max(main_llik[,2])),1]

  if (verbose==T){
    print("Running final gcproc at optimal dimension of: ")
    print(main_dim)
  }

  final.gcproc.model <- gcproc(x = x,
                         y = y,
                         k_dim = main_dim,
                         l_dim = main_dim,
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
