#' @export
initialise.gcproc <- function(data_list,
                              config,
                              anchors=NULL,
                              pivots=NULL){


  if (config$verbose){
    print("Initialising data")
  }

  left_decode <- matrix(rnorm(config$i_dim*config$k_dim),nrow=config$k_dim,ncol=config$i_dim)
  right_decode <- matrix(rnorm(config$k_dim*config$j_dim),nrow=config$j_dim,ncol=config$k_dim)

  main.parameters <- list()
  for (i in 1:length(data_list)){

    initial.param <- list()
    # Initialise parameters
    if (any(do.call('c',lapply(pivots,function(piv){is.null(piv)})))){
      initial.param <-initialise.parameters(x = data_list[[i]], i_dim=config$i_dim,j_dim=config$j_dim,init=config$init,verbose=config$verbose)
    }

    # Check pivoting parameters
    initial.param$pivot_x.sample <- if (is.null(pivots$alpha)){initial.param$pivot_x.sample}else{pivots$alpha}
    initial.param$pivot_x.feature <- if (is.null(pivots$beta)){initial.param$pivot_x.feature}else{pivots$beta}

    # Check anchoring parameters
    alpha <- if (is.null( anchors$alpha)){initial.param$pivot_x.sample}else{ anchors$alpha}
    beta <- if (is.null( anchors$beta)){initial.param$pivot_x.feature}else{ anchors$beta}

    # Find intercept in endecoded space
    X_encode <- (alpha%*%(data_list[[i]])%*%(beta))
    X_code <- (MASS::ginv((alpha)%*%t(alpha))%*%(X_encode)%*%MASS::ginv(t(beta)%*%(beta)))

    decode <- X_code

    right_decode <- t(MASS::ginv(left_decode%*%t(left_decode))%*%left_decode%*%(MASS::ginv((alpha)%*%t(alpha))%*%(X_encode)%*%MASS::ginv(t(beta)%*%(beta))))
    left_decode <- t((MASS::ginv((alpha)%*%t(alpha))%*%(X_encode)%*%MASS::ginv(t(beta)%*%(beta)))%*%right_decode%*%MASS::ginv(t(right_decode)%*%right_decode))

    code <- MASS::ginv((left_decode)%*%t(left_decode))%*%(left_decode)%*%MASS::ginv((alpha)%*%t(alpha))%*%X_encode%*%MASS::ginv(t(beta)%*%(beta))%*%(right_decode)%*%MASS::ginv(t(right_decode)%*%(right_decode))

    main.parameters[[i]] = list(
      right_decode = right_decode,
      left_decode = left_decode,
      alpha = alpha,
      beta = beta
    )

    code = list(
      encode = X_encode,
      decode = decode,
      code = code
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
