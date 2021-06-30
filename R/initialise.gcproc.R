initialise.gcproc <- function(x,y,k_dim=70,j_dim=70,init="svd"){

  x <- Matrix::Matrix(x,sparse=T)
  y <- Matrix::Matrix(y,sparse=T)

  if (init == "svd"){
    set.seed(1)
    u.beta.svd <- irlba::irlba(
      Matrix::crossprod(
        x = Matrix::Matrix(x,sparse=T),
        y = Matrix::Matrix(x,sparse=T)
      ),j_dim,tol=1e-10,maxit = 10000)
    v.beta.svd <- irlba::irlba(
      Matrix::crossprod(
        x = Matrix::Matrix(y,sparse=T),
        y = Matrix::Matrix(y,sparse=T)
      ),j_dim,tol=1e-10,maxit = 10000)

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

  if (init == "eigen-sparse"){
    set.seed(1)
    u.beta.svd <- RSpectra::eigs_sym(
      Matrix::crossprod(
        x = Matrix::Matrix(x,sparse=T),
        y = Matrix::Matrix(x,sparse=T)
      ),j_dim,tol=1e-10,maxit = 10000)
    v.beta.svd <- RSpectra::eigs_sym(
      Matrix::crossprod(
        x = Matrix::Matrix(y,sparse=T),
        y = Matrix::Matrix(y,sparse=T)
      ),j_dim,tol=1e-10,maxit = 10000)

    u.beta.star.beta <- u.beta.svd$vectors[,c(1:j_dim)]
    v.beta.star.beta <- v.beta.svd$vectors[,c(1:j_dim)]

    alpha.L.J.svd <- RSpectra::eigs_sym(
      Matrix::crossprod(
        x = Matrix::t(x),
        y = Matrix::t(x)
      ),k_dim,tol=1e-10,maxit = 10000)
    alpha.L.K.svd <- RSpectra::eigs_sym(
      Matrix::crossprod(
        x = Matrix::t(y),
        y = Matrix::t(y)
      ),k_dim,tol=1e-10,maxit = 10000)

    alpha.L.J.star.alpha.L.J = t(alpha.L.J.svd$vectors[,c(1:k_dim)])
    alpha.L.K.star.alpha.L.K = t(alpha.L.K.svd$vectors[,c(1:k_dim)])
  }

  if (init == "eigen-dense"){
    set.seed(1)
    u.beta.svd <- eigen(
      Matrix::crossprod(
        x = Matrix::Matrix(x,sparse=T),
        y = Matrix::Matrix(x,sparse=T)
      ))
    v.beta.svd <- eigen(
      Matrix::crossprod(
        x = Matrix::Matrix(y,sparse=T),
        y = Matrix::Matrix(y,sparse=T)
      ))

    u.beta.star.beta <- u.beta.svd$vectors[,c(1:j_dim)]
    v.beta.star.beta <- v.beta.svd$vectors[,c(1:j_dim)]

    alpha.L.J.svd <- eigen(
      Matrix::crossprod(
        x = Matrix::t(x),
        y = Matrix::t(x)
      ))
    alpha.L.K.svd <- eigen(
      Matrix::crossprod(
        x = Matrix::t(y),
        y = Matrix::t(y)
      ))

    alpha.L.J.star.alpha.L.J = t(alpha.L.J.svd$vectors[,c(1:k_dim)])
    alpha.L.K.star.alpha.L.K = t(alpha.L.K.svd$vectors[,c(1:k_dim)])
  }


  if (init == "random"){

    u.beta.star.beta <- matrix(rnorm(dim(x)[2]*j_dim),nrow=dim(x)[2],ncol=j_dim)
    v.beta.star.beta <- matrix(rnorm(dim(y)[2]*j_dim),nrow=dim(y)[2],ncol=j_dim)

    alpha.L.J.star.alpha.L.J = matrix(rnorm(dim(x)[1]*k_dim),ncol=dim(x)[1],nrow=k_dim)
    alpha.L.K.star.alpha.L.K = matrix(rnorm(dim(y)[1]*k_dim),ncol=dim(y)[1],nrow=k_dim)

  }



  anchors <- list(  anchor_y.sample = alpha.L.K.star.alpha.L.K,
                    anchor_y.feature = v.beta.star.beta,
                    anchor_x.sample = alpha.L.J.star.alpha.L.J,
                    anchor_x.feature = u.beta.star.beta  )
  return(anchors)

}
