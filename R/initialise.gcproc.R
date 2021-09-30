#' @export
initialise.gcproc <- function(data_list,
                              config,
                              anchors=NULL){


  if (config$verbose){
    print("Initialising data")
  }

  main.parameters <- list()
  for (i in 1:length(data_list)){

    initial.param <- list()
    # Initialise parameters
    initial.param <-initialise.parameters(x = as.matrix(data_list[[i]]), i_dim=config$i_dim,j_dim=config$j_dim,init=config$init,verbose=config$verbose)


    # Check pivoting parameters
    initial.param$pivot_x.sample <- initial.param$pivot_x.sample
    initial.param$pivot_x.feature <- initial.param$pivot_x.feature

    # Check anchoring parameters
    alpha <- initial.param$pivot_x.sample
    beta <- initial.param$pivot_x.feature

    # Find intercept in endecoded space
    X_encode <- (alpha%*%as.matrix(data_list[[i]])%*%(beta))
    X_code <- (MASS::ginv((alpha)%*%t(alpha))%*%(X_encode)%*%MASS::ginv(t(beta)%*%(beta)))

    decode <- X_code

    main.parameters[[i]] = list(
      alpha = alpha,
      beta = beta
    )

    code = list(
      encode = X_encode,
      decode = decode
    )





  }


  return(
    list(
      main.parameters = main.parameters,
      code = code
    )
  )

}


#' @export
initialise.parameters <- function(x,i_dim=70,j_dim=70,init="svd-quick",verbose=F){

  x <- Matrix::Matrix(x,sparse=T)

  set.seed(1)

  if (init=="random"){
    param.beta <- matrix(rnorm(dim(x)[2]*j_dim),nrow=dim(x)[2],ncol=j_dim)
    param.alpha = matrix(rnorm(dim(x)[1]*i_dim),nrow=i_dim,ncol=dim(x)[1])
  } else {
    cov_x <- corpcor::cov.shrink(x,verbose = F)
    cov_tx <- corpcor::cov.shrink(Matrix::t(x),verbose = F)
  }

  if (init=="svd-quick"){
    param.beta.svd <- irlba::irlba(
      cov_x,j_dim,maxit = 10000,verbose = F)
    rm(cov_x)

    param.beta <- param.beta.svd$v


    param.alpha.J.svd <- irlba::irlba(
      cov_tx,i_dim,maxit = 10000,verbose = F)
    rm(cov_tx)

    param.alpha = t(param.alpha.J.svd$u)

  }
  if (init=="svd-dense"){
    param.beta.svd <- svd(
      cov_x,j_dim)
    rm(cov_x)

    param.beta <- param.beta.svd$v[,c(1:j_dim)]

    param.alpha.J.svd <- svd(
      cov_tx,i_dim)
    rm(cov_tx)

    param.alpha = t(param.alpha.J.svd$u[,c(1:i_dim)])

  }
  if (init=="eigen-quick"){
    param.beta.svd <- RSpectra::eigs(
      cov_x,j_dim)
    rm(cov_x)

    param.beta <- param.beta.svd$vectors

    param.alpha.J.svd <- RSpectra::eigs(
      cov_tx,i_dim)
    rm(cov_tx)

    param.alpha = t(param.alpha.J.svd$vectors)

  }
  if (init=="eigen-dense"){

    param.beta.svd <- eigen(
      cov_x,j_dim)
    rm(cov_x)

    param.beta <- param.beta.svd$vectors[,c(1:j_dim)]

    param.alpha.J.svd <- eigen(
      cov_tx)
    rm(cov_tx)

    param.alpha = t(param.alpha.J.svd$vectors[,c(1:i_dim)])

  }


  pivots <- list(
                   pivot_x.sample = param.alpha,
                   pivot_x.feature = param.beta  )
  return(pivots)

}
