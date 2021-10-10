#' @export
initialise.gcproc <- function(data_list,
                              config,
                              transfer
                              ){


  if (config$verbose){
    print(paste("Initialising data with : ",config$init,sep=""))
  }

  main.parameters <- list()
  for (i in 1:length(data_list)){

    initial.param <-initialise.parameters(x = as.matrix(data_list[[i]]), i_dim=config$i_dim,j_dim=config$j_dim,init=config$init,verbose=config$verbose)

    # Check anchoring parameters
    alpha <- initial.param$pivot_x.sample

    beta <- initial.param$pivot_x.feature

    if (is.null(transfer$code)){
      # Find intercept in endecoded space
      encode <- (alpha%*%as.matrix(data_list[[i]])%*%(beta))
      code <- (pinv(t(alpha))%*%(encode)%*%pinv((beta)))

      code = list(
        encode = encode,
        code = code
        )

    } else {

      code <- transfer$code

    }

    main.parameters[[i]] = list(
      alpha = alpha,
      beta = beta
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
initialise.parameters <- function(x,transfer=NULL,i_dim=70,j_dim=70,init="svd-quick",verbose=F){

  x <- Matrix::Matrix(x,sparse=T)

  set.seed(1)

  if (init=="random"){
    param.beta <- if(is.null(transfer$beta)){matrix(rnorm(j_dim),nrow=dim(x)[2],ncol=j_dim)}else{transfer$beta}
    param.alpha = if(is.null(transfer$alpha)){matrix(rnorm(i_dim),nrow=i_dim,ncol=dim(x)[1])}else{transfer$alpha}
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
