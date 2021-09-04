#' @export
initialise.gcproc <- function(x,
                              y,
                              fixed,
                              reference,
                              config,
                              anchors=NULL,
                              pivots=NULL){

  # Initialise parameters
  if (any(do.call('c',lapply(pivots,function(piv){is.null(piv)})))){
    initial.param <-initialise.parameters(x=x,y=y,i_dim=config$i_dim,j_dim=config$j_dim,init=config$init,verbose=config$verbose)
  }

  # Check pivoting parameters
  initial.param$pivot_y.sample <- if (is.null(pivots$pivot_y.sample)){initial.param$pivot_y.sample}else{pivots$pivot_y.sample}
  initial.param$pivot_x.sample <- if (is.null(pivots$pivot_x.sample)){initial.param$pivot_x.sample}else{pivots$pivot_x.sample}
  initial.param$pivot_y.feature <- if (is.null(pivots$pivot_y.feature)){initial.param$pivot_y.feature}else{pivots$pivot_y.feature}
  initial.param$pivot_x.feature <- if (is.null(pivots$pivot_x.feature)){initial.param$pivot_x.feature}else{pivots$pivot_x.feature}

  # Check anchoring parameters
  alpha.K <- if (is.null( anchors$anchor_y.sample)){initial.param$pivot_y.sample}else{ anchors$anchor_y.sample}
  alpha.L <- if (is.null( anchors$anchor_x.sample)){initial.param$pivot_x.sample}else{ anchors$anchor_x.sample}
  v.beta <- if (is.null( anchors$anchor_y.feature)){initial.param$pivot_y.feature}else{ anchors$anchor_y.feature}
  u.beta <- if (is.null( anchors$anchor_x.feature)){initial.param$pivot_x.feature}else{ anchors$anchor_x.feature}

  # Find intercept in endecoded space
  y_encode.final <- alpha.K%*%y%*%v.beta
  x_encode.final <- alpha.L%*%x%*%u.beta

  Y_encode <- (alpha.K%*%(y)%*%(v.beta))
  X_encode <- (alpha.L%*%(x)%*%(u.beta))

  Y_code <- (MASS::ginv((alpha.K)%*%t(alpha.K))%*%Y_encode%*%MASS::ginv(t(v.beta)%*%(v.beta)))

  X_code <- (MASS::ginv((alpha.L)%*%t(alpha.L))%*%(X_encode)%*%MASS::ginv(t(u.beta)%*%(u.beta)))


  Y_decoded <- t(alpha.K)%*%Y_code%*%t(v.beta)
  X_decoded <- t(alpha.L)%*%X_code%*%t(u.beta)


  if (reference == "y"){
    main_code <- Y_code

  }
  if (reference == "x"){
    main_code <- X_code
  }

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
    Y_encode = Y_encode,
    X_decoded = X_decoded,
    Y_decoded = Y_decoded
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

  }
  if (init=="svd-quick"){
    cov_x <- Matrix::crossprod(x,x)
    u.beta.svd <- irlba::irlba(
      cov_x,j_dim,maxit = 10000,verbose = F)
    rm(cov_x)

    cov_y <- Matrix::crossprod(y,y)
    v.beta.svd <- irlba::irlba(
      cov_y,j_dim,maxit = 10000,verbose = F)
    rm(cov_y)

    u.beta <- u.beta.svd$v
    v.beta <- v.beta.svd$v

    cov_tx <- Matrix::crossprod(Matrix::t(x),Matrix::t(x))
    alpha.L.J.svd <- irlba::irlba(
      cov_tx,i_dim,maxit = 10000,verbose = F)
    rm(cov_tx)

    cov_ty <- Matrix::crossprod(Matrix::t(y),Matrix::t(y))
    alpha.L.K.svd <- irlba::irlba(
      cov_ty,i_dim,maxit = 10000,verbose = F)
    rm(cov_ty)

    alpha.L = t(alpha.L.J.svd$u)
    alpha.K = t(alpha.L.K.svd$u)

  }
  if (init=="svd-dense"){
    cov_x <- Matrix::crossprod(x,x)
    u.beta.svd <- svd(
      cov_x,j_dim)
    rm(cov_x)

    cov_y <- Matrix::crossprod(y,y)
    v.beta.svd <- svd(
      cov_y,j_dim)
    rm(cov_y)

    u.beta <- u.beta.svd$v[,c(1:j_dim)]
    v.beta <- v.beta.svd$v[,c(1:j_dim)]

    cov_tx <- Matrix::crossprod(Matrix::t(x),Matrix::t(x))
    alpha.L.J.svd <- svd(
      cov_tx,i_dim)
    rm(cov_tx)

    cov_ty <- Matrix::crossprod(Matrix::t(y),Matrix::t(y))
    alpha.L.K.svd <- svd(
      cov_ty,i_dim)
    rm(cov_ty)

    alpha.L = t(alpha.L.J.svd$u[,c(1:i_dim)])
    alpha.K = t(alpha.L.K.svd$u[,c(1:i_dim)])

  }
  if (init=="eigen-quick"){
    cov_x <- Matrix::crossprod(x,x)
    u.beta.svd <- RSpectra::eigs(
      cov_x,j_dim)
    rm(cov_x)

    cov_y <- Matrix::crossprod(y,y)
    v.beta.svd <- RSpectra::eigs(
      cov_y,j_dim)
    rm(cov_y)

    u.beta <- u.beta.svd$vectors
    v.beta <- v.beta.svd$vectors

    cov_tx <- Matrix::crossprod(Matrix::t(x),Matrix::t(x))
    alpha.L.J.svd <- RSpectra::eigs(
      cov_tx,i_dim)
    rm(cov_tx)

    cov_ty <- Matrix::crossprod(Matrix::t(y),Matrix::t(y))
    alpha.L.K.svd <- RSpectra::eigs(
      cov_ty,i_dim)
    rm(cov_ty)

    alpha.L = t(alpha.L.J.svd$vectors)
    alpha.K = t(alpha.L.K.svd$vectors)

  }
  if (init=="eigen-dense"){
    cov_x <- Matrix::crossprod(x,x)
    u.beta.svd <- eigen(
      cov_x,j_dim)
    rm(cov_x)

    cov_y <- Matrix::crossprod(y,y)
    v.beta.svd <- eigen(
      cov_y,j_dim)
    rm(cov_y)

    u.beta <- u.beta.svd$vectors[,c(1:j_dim)]
    v.beta <- v.beta.svd$vectors[,c(1:j_dim)]

    cov_tx <- Matrix::crossprod(Matrix::t(x),Matrix::t(x))
    alpha.L.J.svd <- eigen(
      cov_tx)
    rm(cov_tx)

    cov_ty <- Matrix::crossprod(Matrix::t(y),Matrix::t(y))
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
