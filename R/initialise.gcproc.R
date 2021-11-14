#' @export
initialise.gcproc <- function(data_list,
                              config,
                              covariate,
                              transfer,
                              join
                              ){
  
  if (config$verbose){
    print(paste("Initialising data with : ",config$init,sep=""))
  }
  
  main.code <- list(code=list(),encode=list())
  main.parameters <- list(alpha = list(), beta = list())
  
  for (i in 1:length(data_list)){
    
    if ( !is.null(transfer$main.parameters$alpha[[1]]) & !is.null(transfer$main.parameters$beta[[1]]) ){
      transfer.param <- list(main.parameters = list(
        alpha = transfer$main.parameters$alpha[[join$alpha[i]]],
        beta = transfer$main.parameters$beta[[join$beta[i]]]
        )
      )
    } else {
      transfer.param <- transfer$main.parameters
    }
    
    initial.param <-initialise.parameters(
      x = as.matrix(data_list[[i]]), 
      config = config, 
      transfer = transfer.param
    )
    
    # Check anchoring parameters
    alpha <- initial.param$pivot_x.sample
    beta <- initial.param$pivot_x.feature
    
    main.parameters$alpha[[join$alpha[i]]] <- alpha
    main.parameters$beta[[join$beta[i]]] <- beta
    
    if (is.null(transfer$main.code)){
      
      encode <- (alpha%*%as.matrix(data_list[[i]])%*%(beta))
      code <- (pinv(t(alpha))%*%(encode)%*%pinv((beta)))
      
      main.code = list(
        encode = encode,
        code = code
      )
      
    } else {    
      main.code <- transfer$main.code
    }
    
  }
  
  for (iter in 1:1){
    for (i in 1:length(data_list)){
      
      internal.parameters <- list(alpha=main.parameters$alpha[[join$alpha[i]]],
                                  beta=main.parameters$beta[[join$beta[i]]])
      
      return_update <- update_set(x = as.matrix(data_list[[i]]),
                                  main.parameters = internal.parameters,
                                  main.code = main.code, 
                                  method = "svd",
                                  pivots = list(alpha = c(1:config$i_dim) ,beta = c(1:config$j_dim))
      )
      
      main.parameters$alpha[[join$alpha[i]]] <- return_update$main.parameters$alpha
      main.parameters$beta[[join$beta[i]]] <- return_update$main.parameters$beta
      
      main.code <- return_update$main.code
      
    }
  }
  
  return(
    list(
      main.parameters = main.parameters,
      main.code = main.code
    )
  )
  
}


#' @export
initialise.parameters <- function(x,config,transfer){
  
  if (is.null(transfer$main.parameters$beta) | is.null(transfer$main.parameters$alpha) ){
    x <- Matrix::Matrix(x,sparse=T)
    
    set.seed(config$seed)
  }
  
  param.beta <-   if (!is.null(transfer$main.parameters$beta)){
    transfer$main.parameters$beta
  } else if (config$init=="random"){
    array(rnorm(dim(x)[2]*config$j_dim),dim=c(dim(x)[2],config$j_dim))
  } 
  
  param.alpha <- if (!is.null(transfer$main.parameters$alpha)){
    transfer$main.parameters$alpha
  } else if (config$init=="random") {
    array(rnorm(config$i_dim*dim(x)[1]),dim=c(config$i_dim,dim(x)[1]))
  } 
  
  pivots <- list(
    pivot_x.sample = as.matrix(param.alpha),
    pivot_x.feature = as.matrix(param.beta)  
    )
  
  
  return(pivots)
  
}
