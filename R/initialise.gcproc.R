#' @export
initialise.gcproc <- function(x,y,k_dim=70,j_dim=70,init="svd-quick",verbose=F){

  x <- Matrix::Matrix(x,sparse=T)
  y <- Matrix::Matrix(y,sparse=T)

  set.seed(1)


  if (init=="random"){
    u.beta.star.beta <- matrix(rnorm(dim(x)[2]*j_dim),nrow=dim(x)[2],ncol=j_dim)
    v.beta.star.beta <- matrix(rnorm(dim(y)[2]*j_dim),nrow=dim(y)[2],ncol=j_dim)

    alpha.L.J.star.alpha.L.J = matrix(rnorm(dim(x)[1]*k_dim),nrow=k_dim,ncol=dim(x)[1])
    alpha.L.K.star.alpha.L.K = matrix(rnorm(dim(y)[1]*k_dim),nrow=k_dim,ncol=dim(y)[1])

  }
  if (init=="svd-quick"){
    cov_x <- Matrix::crossprod(x,x)
    u.beta.svd <- irlba::irlba(
      cov_x,j_dim,tol=1e-10,maxit = 10000)
    rm(cov_x)

    cov_y <- Matrix::crossprod(y,y)
    v.beta.svd <- irlba::irlba(
      cov_y,j_dim,tol=1e-10,maxit = 10000)
    rm(cov_y)

    u.beta.star.beta <- u.beta.svd$v
    v.beta.star.beta <- v.beta.svd$v

    cov_tx <- Matrix::crossprod(Matrix::t(x),Matrix::t(x))
    alpha.L.J.svd <- irlba::irlba(
      cov_tx,k_dim,tol=1e-10,maxit = 10000)
    rm(cov_tx)

    cov_ty <- Matrix::crossprod(Matrix::t(y),Matrix::t(y))
    alpha.L.K.svd <- irlba::irlba(
      cov_ty,k_dim,tol=1e-10,maxit = 10000)
    rm(cov_ty)

    alpha.L.J.star.alpha.L.J = t(alpha.L.J.svd$u)
    alpha.L.K.star.alpha.L.K = t(alpha.L.K.svd$u)

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

    u.beta.star.beta <- u.beta.svd$v[,c(1:j_dim)]
    v.beta.star.beta <- v.beta.svd$v[,c(1:j_dim)]

    cov_tx <- Matrix::crossprod(Matrix::t(x),Matrix::t(x))
    alpha.L.J.svd <- svd(
      cov_tx,k_dim)
    rm(cov_tx)

    cov_ty <- Matrix::crossprod(Matrix::t(y),Matrix::t(y))
    alpha.L.K.svd <- svd(
      cov_ty,k_dim)
    rm(cov_ty)

    alpha.L.J.star.alpha.L.J = t(alpha.L.J.svd$u[,c(1:k_dim)])
    alpha.L.K.star.alpha.L.K = t(alpha.L.K.svd$u[,c(1:k_dim)])

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

    u.beta.star.beta <- u.beta.svd$vectors
    v.beta.star.beta <- v.beta.svd$vectors

    cov_tx <- Matrix::crossprod(Matrix::t(x),Matrix::t(x))
    alpha.L.J.svd <- RSpectra::eigs(
      cov_tx,k_dim)
    rm(cov_tx)

    cov_ty <- Matrix::crossprod(Matrix::t(y),Matrix::t(y))
    alpha.L.K.svd <- RSpectra::eigs(
      cov_ty,k_dim)
    rm(cov_ty)

    alpha.L.J.star.alpha.L.J = t(alpha.L.J.svd$vectors)
    alpha.L.K.star.alpha.L.K = t(alpha.L.K.svd$vectors)

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

    u.beta.star.beta <- u.beta.svd$vectors[,c(1:j_dim)]
    v.beta.star.beta <- v.beta.svd$vectors[,c(1:j_dim)]

    cov_tx <- Matrix::crossprod(Matrix::t(x),Matrix::t(x))
    alpha.L.J.svd <- eigen(
      cov_tx)
    rm(cov_tx)

    cov_ty <- Matrix::crossprod(Matrix::t(y),Matrix::t(y))
    alpha.L.K.svd <- eigen(
      cov_ty)
    rm(cov_ty)

    alpha.L.J.star.alpha.L.J = t(alpha.L.J.svd$vectors[,c(1:k_dim)])
    alpha.L.K.star.alpha.L.K = t(alpha.L.K.svd$vectors[,c(1:k_dim)])

  }


  pivots <- list(  pivot_y.sample = alpha.L.K.star.alpha.L.K,
                    pivot_y.feature = v.beta.star.beta,
                    pivot_x.sample = alpha.L.J.star.alpha.L.J,
                    pivot_x.feature = u.beta.star.beta  )
  return(pivots)

}
