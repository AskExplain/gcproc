initialise.gcproc <- function(x,y,k_dim=70,j_dim=70,init="svd",verbose=F){

  x <- Matrix::Matrix(x,sparse=T)
  y <- Matrix::Matrix(y,sparse=T)

  set.seed(1)

  cov_x <- corpcor::cov.shrink((x),lambda = 0.99,lambda.var = 0.99,verbose = verbose)
  u.beta.svd <- irlba::irlba(
    cov_x,j_dim,tol=1e-10,maxit = 10000)
  rm(cov_x)

  cov_y <- corpcor::cov.shrink((y),lambda = 0.99,lambda.var = 0.99,verbose = verbose)
  v.beta.svd <- irlba::irlba(
    cov_y,j_dim,tol=1e-10,maxit = 10000)
  rm(cov_y)

  u.beta.star.beta <- u.beta.svd$v
  v.beta.star.beta <- v.beta.svd$v

  cov_tx <- corpcor::cov.shrink(Matrix::t(x),lambda = 0.99,lambda.var = 0.99,verbose = verbose)
  alpha.L.J.svd <- irlba::irlba(
    cov_tx,k_dim,tol=1e-10,maxit = 10000)
  rm(cov_tx)

  cov_ty <- corpcor::cov.shrink(Matrix::t(y),lambda = 0.99,lambda.var = 0.99,verbose = verbose)
  alpha.L.K.svd <- irlba::irlba(
    cov_ty,k_dim,tol=1e-10,maxit = 10000)
  rm(cov_ty)

  alpha.L.J.star.alpha.L.J = t(alpha.L.J.svd$u)
  alpha.L.K.star.alpha.L.K = t(alpha.L.K.svd$u)

  anchors <- list(  anchor_y.sample = alpha.L.K.star.alpha.L.K,
                    anchor_y.feature = v.beta.star.beta,
                    anchor_x.sample = alpha.L.J.star.alpha.L.J,
                    anchor_x.feature = u.beta.star.beta  )
  return(anchors)

}
