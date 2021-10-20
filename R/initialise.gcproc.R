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



  main.code <- list(code=list(),encode=list())
  main.index <- list()
  main.parameters <- list(alpha = list(), beta = list())
  main.proportion <- list()

  index <- list()

  for (i in 1:length(data_list)){

    if (covariate$fix){
      main.proportion[[i]] <- cbind(covariate$factor[[i]])
    } else {
      main.proportion[[i]] <- cbind(covariate$factor[[i]]+array(prod(dim(rnorm(covariate$factor[[i]]))),dim=c(dim(covariate$factor[[i]]))))
    }
    # main.proportion[[i]] <- main.proportion[[i]] / rowSums(main.proportion[[i]])
    index$code_indicator <- colnames(covariate$factor[[i]])


    initial.param <-initialise.parameters(x = as.matrix(data_list[[i]]),transfer = transfer, i_dim=config$i_dim,j_dim=config$j_dim,init=config$init,verbose=config$verbose)

    alpha <- initial.param$pivot_x.sample
    beta <- initial.param$pivot_x.feature
    main.code$encode <- (alpha%*%as.matrix(data_list[[i]])%*%(beta))/length(index$code_indicator)

    alpha.list <- list()
    beta.list <- list()
    code.list <- list()

    for (j in c(1:length(index$code_indicator))){
      alpha.list <- c(alpha.list,list(alpha*rnorm(prod(dim(alpha)))))
      beta.list <- c(beta.list,list(beta*rnorm(prod(dim(beta)))))

      code.list <- c(code.list,list((pinv(t(alpha))%*%(main.code$encode)%*%pinv((beta)))*rnorm(prod(dim(main.code$encode)))))
    }


    main.parameters$alpha[[i]] <- alpha.list
    main.parameters$beta[[i]] <- beta.list

    if (is.null(transfer$code)){

      main.code$code <- code.list

    } else {

      main.code <- transfer$code

    }

    names(main.code$code) <- index$code_indicator

  }



  return(
    list(
      main.parameters = main.parameters,
      main.code = main.code,
      main.index = main.index,
      main.proportion = main.proportion
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
