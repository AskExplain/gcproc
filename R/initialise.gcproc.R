#' @export
initialise.gcproc <- function(data_list,
                              config,
                              covariate,
                              transfer,
                              join,
                              pivots){

  index <- list()

  index$code_indicator <- do.call('c',lapply(c(1:dim(covariate$factor)[2]),function(X){c(unique(covariate$factor[,X]))}))



  if (config$verbose){
    print(paste("Initialising data with : ",config$init,sep=""))
  }

  main.code <- list(code=list(),encode=list())
  main.index <- list()
  main.parameters <- list(alpha = list(), beta = list())
  for (i in 1:length(data_list)){

    main.index[[i]] <- rbind(data.frame(factor = index$code_indicator,update = as.integer(index$code_indicator %in% c(covariate$factor[i,]))))


    initial.param <-initialise.parameters(x = as.matrix(data_list[[i]]),transfer = transfer, i_dim=config$i_dim,j_dim=config$j_dim,init=config$init,verbose=config$verbose)

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

    main.code$encode <- code$encode
    for (code.id in which(index$code_indicator %in% c(covariate$factor[i,]))){
      main.code$code[[code.id]] <- code$code*rnorm(prod(dim(code$code)))
    }

    alpha[-pivots$alpha,] <- 0
    beta[,-pivots$beta] <- 0

    main.parameters$alpha[[join$alpha[i]]] <- alpha
    main.parameters$beta[[join$beta[i]]] <- beta

  }


  names(main.code$code) <- index$code_indicator
  names(main.parameters$alpha) <- unique(join$alpha)
  names(main.parameters$beta) <- unique(join$beta)


  return(
    list(
      main.parameters = main.parameters,
      main.code = main.code,
      main.index = main.index
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
