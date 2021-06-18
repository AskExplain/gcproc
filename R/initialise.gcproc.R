initialise.gcproc <- function(x,y,k_dim=70,j_dim=70,init="svd"){

  x <- Matrix::Matrix(x,sparse=T)
  y <- Matrix::Matrix(y,sparse=T)

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


  anchors <- list(  anchor_y.sample = alpha.L.K.star.alpha.L.K,
                    anchor_y.feature = v.beta.star.beta,
                    anchor_x.sample = alpha.L.J.star.alpha.L.J,
                    anchor_x.feature = u.beta.star.beta  )
  return(anchors)

}
