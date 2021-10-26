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
  main.parameters <- list(alpha = list(), beta = list())
  
  for (i in 1:length(data_list)){
    
    initial.param <-initialise.parameters(x = as.matrix(data_list[[i]]),transfer = transfer, i_dim=config$i_dim,j_dim=config$j_dim,init=config$init,verbose=config$verbose)
    
    # Check anchoring parameters
    alpha <- initial.param$pivot_x.sample
    beta <- initial.param$pivot_x.feature
    
    main.parameters$alpha[[join$alpha[i]]] <- alpha
    main.parameters$beta[[join$beta[i]]] <- beta
    
    if (is.null(transfer$code)){
      
      # Find intercept in endecoded space
      encode <- (alpha%*%as.matrix(data_list[[i]])%*%(beta))
      code <- (pinv(t(alpha))%*%(encode)%*%pinv((beta)))
      
      main.code = list(
        encode = encode,
        code = code
      )
      
    } else {    
      main.code <- transfer$code
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
initialise.parameters <- function(x,transfer,i_dim,j_dim,init="random",verbose=F){
  
  x <- Matrix::Matrix(x,sparse=T)
  
  set.seed(1)
  
  if (init=="random"){
    param.beta <- array(rnorm(dim(x)[2]*config$j_dim),dim=c(dim(x)[2],config$j_dim))
    param.alpha = array(rnorm(config$i_dim*dim(x)[1]),dim=c(config$i_dim,dim(x)[1]))
  } 
  
  if (init=="irlba"){
    param.beta.svd <- irlba::irlba(
      x,j_dim)
    param.beta <- param.beta.svd$v
    
    param.alpha.svd <- irlba::irlba(
      x,i_dim)
    param.alpha = t(param.alpha.svd$u)
  }
  if (init=="svdr"){
    param.beta.svd <- irlba::svdr(
      x,j_dim)
    param.beta <- param.beta.svd$v
    
    param.alpha.svd <- irlba::svdr(
      x,i_dim)
    param.alpha = t(param.alpha.svd$u)
  }

  pivots <- list(
    pivot_x.sample = as.matrix(param.alpha),
    pivot_x.feature = as.matrix(param.beta)  )
  return(pivots)
  
}
