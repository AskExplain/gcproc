#' @export
initialise.gcproc <- function(data_list,
                              config,
                              covariate,
                              transfer,
                              join,
                              pivots){

  if (config$verbose){
    print(paste("Initialising data with : ",config$init,sep=""))
  }

  alpha.list <- list()
  beta.list <- list()

  main.code <- list(code=list(),encode=list())
  main.parameters <- list(alpha = list(), beta = list())
  for (i in 1:length(data_list)){

    initial.param <-initialise.parameters(x = as.matrix(data_list[[i]]),transfer = transfer, i_dim=config$i_dim,j_dim=config$j_dim,init=config$init,verbose=config$verbose)

    alpha <- initial.param$pivot_x.sample
    beta <- initial.param$pivot_x.feature

    encode.d <- (alpha%*%as.matrix(data_list[[i]])%*%(beta))
    code.d <- (pinv(t(alpha))%*%(encode.d)%*%pinv((beta)))

    if (is.null(transfer$code)){

      main.code = list(
        encode = encode.d,
        code = code.d
      )

    } else {

      main.code <- transfer$code

    }

    alpha.list <- c(alpha.list,list(alpha))
    beta.list <- c(beta.list,list(beta))

  }

  main.parameters$alpha <- alpha.list
  main.parameters$beta <- beta.list


  return(
    list(
      main.parameters = main.parameters,
      main.code = main.code
      )
  )

}


#' @export
initialise.parameters <- function(x,transfer,i_dim,j_dim,init="svd",verbose=F){

  x <- Matrix::Matrix(x,sparse=T)

  set.seed(1)

  if (init=="random"){
    param.beta <- if(is.null(transfer$beta)){array(rnorm(config$j_dim),dim=c(dim(x)[2],config$j_dim))}else{transfer$beta}
    param.alpha = if(is.null(transfer$alpha)){array(rnorm(config$i_dim),dim=c(config$i_dim,dim(x)[1]))}else{transfer$alpha}
  } else {
    cov_x <- corpcor::cov.shrink(x,verbose = F)
    cov_tx <- corpcor::cov.shrink(Matrix::t(x),verbose = F)
  }

  if (init=="svd"){
    param.beta.svd <- irlba::irlba(
      cov_x,j_dim,verbose = F)
    rm(cov_x)

    param.beta <- if(is.null(transfer$beta)){param.beta.svd$v}else{transfer$beta}


    param.alpha.J.svd <- irlba::irlba(
      cov_tx,i_dim,verbose = F)
    rm(cov_tx)

    param.alpha = if(is.null(transfer$alpha)){t(param.alpha.J.svd$u)}else{transfer$alpha}

  }

  pivots <- list(
    pivot_x.sample = as.matrix(param.alpha),
    pivot_x.feature = as.matrix(param.beta)  )
  return(pivots)

}
