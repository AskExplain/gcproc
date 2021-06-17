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
                      init="svd"){

  j_dim <- 70
  k_dim <- 70
  param.dim <- data.frame(c(1:seeds),k_dim,j_dim)


  if (init == "svd"){

    set.seed(1)
    u.beta.svd <- irlba::irlba(
      Matrix::crossprod(
        x = Matrix::Matrix(x,sparse=T),
        y = Matrix::Matrix(x,sparse=T)
      ),k_dim,tol=1e-10,maxit = 10000)
    v.beta.svd <- irlba::irlba(
      Matrix::crossprod(
        x = Matrix::Matrix(y,sparse=T),
        y = Matrix::Matrix(y,sparse=T)
      ),k_dim,tol=1e-10,maxit = 10000)

    u.beta.star.beta <- u.beta.svd$v
    v.beta.star.beta <- v.beta.svd$v

    alpha.L.J.svd <- irlba::irlba(
      Matrix::crossprod(
      x = Matrix::t(x),
      y = Matrix::t(x)
    ),k_dim,tol=1e-10,maxit = 10000)
    alpha.L.K.svd <- irlba::irlba(
      Matrix::crossprod(
      x = Matrix::t(y),
      y = Matrix::t(y)
    ),k_dim,tol=1e-10,maxit = 10000)

    alpha.L.J.star.alpha.L.J = t(alpha.L.J.svd$u)
    alpha.L.K.star.alpha.L.K = t(alpha.L.K.svd$u)

  }
  if (init == "random"){
    set.seed(1)

    u.beta.star.beta <- matrix(rnorm(j_dim*dim(x)[2]),ncol = j_dim,nrow=dim(x)[2])
    v.beta.star.beta <- matrix(rnorm(j_dim*dim(y)[2]),ncol = j_dim,nrow=dim(y)[2])

    alpha.L.J.star.alpha.L.J <- matrix(rnorm(k_dim*dim(x)[1]),nrow = k_dim,ncol=dim(x)[1])
    alpha.L.K.star.alpha.L.K <- matrix(rnorm(k_dim*dim(y)[1]),nrow = k_dim,ncol=dim(y)[1])

  }


  anchors <- list(  anchor_y.sample = alpha.L.K.star.alpha.L.K,
                    anchor_y.feature = v.beta.star.beta,
                    anchor_x.sample = alpha.L.J.star.alpha.L.J,
                    anchor_x.feature = u.beta.star.beta  )


  main_llik <- c()
  for (row in 1:nrow(param.dim)){
    seed <- param.dim[row,1]
    set.seed(seed)
    k <- param.dim[row,2]
    j <- param.dim[row,3]

    if (verbose==T){
      print("Beginning tuning dimension to: ")
      print(paste("seed: ",seed, "   k_dim: ",k,"   j_dim: ",j,sep=""))
    }

    gcproc.model <- try(gcproc(y = y,
                               x = x,
                               k_dim = k,
                               j_dim = j,
                               eta = eta,
                               max_iter = 15,
                               min_iter = 15,
                               tol = tol,
                               log = log,
                               center = center,
                               scale.z = scale.z,
                               batches = batches,
                               cores = cores,
                               verbose = verbose,
                               init= "anchors",
                               seed=seed,
                               anchors = anchors),silent = F)
    if (!is.character(gcproc.model)){
      cost <- tail(gcproc.model$convergence.parameters$llik.vec,1)
      main_llik <- rbind(main_llik,c(seed,k,j,cost))
    } else {
      main_llik <- rbind(main_llik,c(-Inf))
    }
  }

  main_dim <- main_llik[which(main_llik[,4]==max(main_llik[,4])),]
  main_seed <- main_dim[1]
  main_k_dim <- main_dim[2]
  main_j_dim <- main_dim[3]

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
                               init="anchors",
                               anchor=anchors,
                               seed = main_seed)

  final.gcproc.model$cv.llik <- main_llik
  final.gcproc.model$meta.parameters$seeds <- seeds

  return(final.gcproc.model)
}
