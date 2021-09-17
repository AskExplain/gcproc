#' @export
initialise.gcproc <- function(x,
                              y,
                              config,
                              anchors=NULL,
                              pivots=NULL){


  if (config$verbose){
    print("Initialising data")
  }

  initial.param <- list()
  # Initialise parameters
  if (any(do.call('c',lapply(pivots,function(piv){is.null(piv)})))){
    initial.param <-initialise.parameters(x=x,y=y,i_dim=config$i_dim,j_dim=config$j_dim,init=config$init,verbose=config$verbose)
  }

  # Check pivoting parameters
  initial.param$pivot_y.sample <- if (is.null(pivots$alpha.K)){initial.param$pivot_y.sample}else{pivots$alpha.K}
  initial.param$pivot_x.sample <- if (is.null(pivots$alpha.L)){initial.param$pivot_x.sample}else{pivots$alpha.L}
  initial.param$pivot_y.feature <- if (is.null(pivots$v.beta)){initial.param$pivot_y.feature}else{pivots$v.beta}
  initial.param$pivot_x.feature <- if (is.null(pivots$u.beta)){initial.param$pivot_x.feature}else{pivots$u.beta}

  # Check anchoring parameters
  alpha.K <- if (is.null( anchors$alpha.K)){initial.param$pivot_y.sample}else{ anchors$alpha.K}
  alpha.L <- if (is.null( anchors$alpha.L)){initial.param$pivot_x.sample}else{ anchors$alpha.L}
  v.beta <- if (is.null( anchors$v.beta)){initial.param$pivot_y.feature}else{ anchors$v.beta}
  u.beta <- if (is.null( anchors$u.beta)){initial.param$pivot_x.feature}else{ anchors$u.beta}

  # Find intercept in endecoded space
  y_encode.final <- alpha.K%*%y%*%v.beta
  x_encode.final <- alpha.L%*%x%*%u.beta

  Y_encode <- (alpha.K%*%(y)%*%(v.beta))
  X_encode <- (alpha.L%*%(x)%*%(u.beta))

  Y_code <- (MASS::ginv((alpha.K)%*%t(alpha.K))%*%Y_encode%*%MASS::ginv(t(v.beta)%*%(v.beta)))

  X_code <- (MASS::ginv((alpha.L)%*%t(alpha.L))%*%(X_encode)%*%MASS::ginv(t(u.beta)%*%(u.beta)))

  main_code <- Y_code

  main.parameters = list(
    alpha.L = alpha.L,
    alpha.K = alpha.K,
    u.beta = u.beta,
    v.beta = v.beta
  )

  code = list(
    main_code = main_code,
    X_code = X_code,
    Y_code = Y_code,
    X_encode = X_encode,
    Y_encode = Y_encode
  )



  return(
    list(
      main.parameters = main.parameters,
      code = code
    )
  )

}


#' @export
initialise.parameters <- function(x,y,i_dim=70,j_dim=70,init="svd-quick",verbose=F){

  x <- Matrix::Matrix(x,sparse=T)
  y <- Matrix::Matrix(y,sparse=T)

  set.seed(1)

  if (init=="random"){
    u.beta <- matrix(rnorm(dim(x)[2]*j_dim),nrow=dim(x)[2],ncol=j_dim)
    v.beta <- matrix(rnorm(dim(y)[2]*j_dim),nrow=dim(y)[2],ncol=j_dim)

    alpha.L = matrix(rnorm(dim(x)[1]*i_dim),nrow=i_dim,ncol=dim(x)[1])
    alpha.K = matrix(rnorm(dim(y)[1]*i_dim),nrow=i_dim,ncol=dim(y)[1])

  } else {

    cov_x <- corpcor::cov.shrink(x,verbose = F)
    cov_y <- corpcor::cov.shrink(y,verbose = F)
    cov_tx <- corpcor::cov.shrink(Matrix::t(x),verbose = F)
    cov_ty <- corpcor::cov.shrink(Matrix::t(y),verbose = F)

  }

  if (init=="svd-quick"){
    u.beta.svd <- irlba::irlba(
      cov_x,j_dim,maxit = 10000,verbose = F)
    rm(cov_x)

    v.beta.svd <- irlba::irlba(
      cov_y,j_dim,maxit = 10000,verbose = F)
    rm(cov_y)

    u.beta <- u.beta.svd$v
    v.beta <- v.beta.svd$v


    alpha.L.J.svd <- irlba::irlba(
      cov_tx,i_dim,maxit = 10000,verbose = F)
    rm(cov_tx)


    alpha.L.K.svd <- irlba::irlba(
      cov_ty,i_dim,maxit = 10000,verbose = F)
    rm(cov_ty)

    alpha.L = t(alpha.L.J.svd$u)
    alpha.K = t(alpha.L.K.svd$u)

  }
  if (init=="svd-dense"){
    u.beta.svd <- svd(
      cov_x,j_dim)
    rm(cov_x)

    v.beta.svd <- svd(
      cov_y,j_dim)
    rm(cov_y)

    u.beta <- u.beta.svd$v[,c(1:j_dim)]
    v.beta <- v.beta.svd$v[,c(1:j_dim)]


    alpha.L.J.svd <- svd(
      cov_tx,i_dim)
    rm(cov_tx)


    alpha.L.K.svd <- svd(
      cov_ty,i_dim)
    rm(cov_ty)

    alpha.L = t(alpha.L.J.svd$u[,c(1:i_dim)])
    alpha.K = t(alpha.L.K.svd$u[,c(1:i_dim)])

  }
  if (init=="eigen-quick"){
    u.beta.svd <- RSpectra::eigs(
      cov_x,j_dim)
    rm(cov_x)

    v.beta.svd <- RSpectra::eigs(
      cov_y,j_dim)
    rm(cov_y)

    u.beta <- u.beta.svd$vectors
    v.beta <- v.beta.svd$vectors

    alpha.L.J.svd <- RSpectra::eigs(
      cov_tx,i_dim)
    rm(cov_tx)

    alpha.L.K.svd <- RSpectra::eigs(
      cov_ty,i_dim)
    rm(cov_ty)

    alpha.L = t(alpha.L.J.svd$vectors)
    alpha.K = t(alpha.L.K.svd$vectors)

  }
  if (init=="eigen-dense"){

    u.beta.svd <- eigen(
      cov_x,j_dim)
    rm(cov_x)


    v.beta.svd <- eigen(
      cov_y,j_dim)
    rm(cov_y)

    u.beta <- u.beta.svd$vectors[,c(1:j_dim)]
    v.beta <- v.beta.svd$vectors[,c(1:j_dim)]


    alpha.L.J.svd <- eigen(
      cov_tx)
    rm(cov_tx)


    alpha.L.K.svd <- eigen(
      cov_ty)
    rm(cov_ty)

    alpha.L = t(alpha.L.J.svd$vectors[,c(1:i_dim)])
    alpha.K = t(alpha.L.K.svd$vectors[,c(1:i_dim)])

  }


  pivots <- list(  pivot_y.sample = alpha.K,
                   pivot_y.feature = v.beta,
                   pivot_x.sample = alpha.L,
                   pivot_x.feature = u.beta  )
  return(pivots)

}
