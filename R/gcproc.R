#' Generalised Canonical Procrustes
#'
#' A method that uses a likelihood model to align two datasets via an encoding in a lower dimensional space
#'
#' @param x Reference dataset of sample by feature matrix (required)
#' @param y Experimental dataset of sample by feature matrix (required)
#' @param config Configuration parameters. Modifiable list accessed from extract_config()
#' @param seed Fixed seed (default = 1)
#' @param anchors Transferring pre-trained model parameters (not required)
#' @param pivots Initialisation of model parameters (not required)
#' @param covariates Allows for covariate adjustment (not required)
#' @return  Main parameters contains the learned model parameters.
#' @return alpha.L.K:     To reduce rows of y to k_dim (the reduced form is t(KY))
#' @return alpha.L.J:     To reduce rows of x to to k_dim (the reduced form is t(JX))
#' @return v.beta:        To reduce columns of y to to j_dim (the reduced form is (Yv))
#' @return u.beta:        To reduce columns of x to to j_dim (the reduced form is (Yu))
#' @return y.gamma:       Covariate parameters for the sample covariates of y
#' @return y.delta:       Covariate parameters for the feature covariates of y
#' @return x.gamma:       Covariate parameters for the sample covariates of x
#' @return x.delta:       Covariate parameters for the feature covariates of y
#' @return intercept.x:   Intercept for the encoding of x (where KYv = JXu + intercept)
#' @export
gcproc <- function(x,
                   y,
                   config = list(k_dim = 70,
                                 j_dim = 70,
                                 eta=1e-1,
                                 max_iter=15,
                                 min_iter = 5,
                                 tol=1e-3,
                                 log=F,
                                 center=F,
                                 scale.z=F,
                                 batches=2,
                                 cores=2,
                                 verbose=T,
                                 init="svd-quick"),
                   seed = 1,
                   anchors = NULL,
                   pivots = NULL,
                   covariates = NULL,
                   predict = NULL

){

  if (is.null(config$batch_size)){
    config$batch_size <- 1000
    config$batches <- ceiling(max(c(dim(x),dim(y)))/config$batch_size)
  }


  b.a <- config$eta
  a.b <- c(1-b.a)

  prepare_data = TRUE
  initialise = TRUE
  variational_gradient_descent_updates = TRUE
  run_covariates = TRUE



  anchor_y.sample = NULL
  anchor_y.feature = NULL
  anchor_x.sample = NULL
  anchor_x.feature = NULL
  anchor_y.cov.sample = NULL
  anchor_y.cov.feature = NULL
  anchor_x.cov.sample = NULL
  anchor_x.cov.feature = NULL

  if (!is.null(anchors)){
    check_anchors = TRUE

    anchors$anchor_y.sample = anchors$anchor_y.sample
    anchors$anchor_y.feature = anchors$anchor_y.feature
    anchors$anchor_x.sample = anchors$anchor_x.sample
    anchors$anchor_x.feature = anchors$anchor_x.feature
    anchors$anchor_y.cov.sample = anchors$anchor_y.cov.sample
    anchors$anchor_y.cov.feature = anchors$anchor_y.cov.feature
    anchors$anchor_x.cov.sample = anchors$anchor_x.cov.sample
    anchors$anchor_x.cov.feature = anchors$anchor_x.cov.feature

  }

  if (config$verbose){
    print(paste("Using gcproc with following dimension reductions:   Sample dimension (config$k_dim): ",config$k_dim, "   Feature dimension (config$j_dim): ", config$j_dim,sep=""))
  }

  n.x <- dim(x)[1]
  p.x <- dim(x)[2]
  n.y <- dim(y)[1]
  p.y <- dim(y)[2]

  if (prepare_data == T){
    X.x <- prepare_data(x=x,
                        log=config$log,
                        center=config$center,
                        scale=config$scale.z)

    Y.y <- prepare_data(x=y,
                        log=config$log,
                        center=config$center,
                        scale=config$scale.z)


    if (config$verbose){
      print("Prepared transformations on data are (True/False):")
      print(paste("Natural Log (with pseudocount): ", config$log, sep = ""))
      print(paste("Center: ", config$center, sep = ""))
      print(paste("Standardise by size of norm: ", config$scale.z, sep = ""))
    }

  }


  if (initialise==T){
    if (config$verbose){
      print("Initialising data")
    }

    # Prepare convergence checking parameters
    count = 1
    llik.vec <- rep(0,10)
    score.vec <- rep(0,10)


    # Prepare priors
    a0.beta = 10e-2
    b0.beta = 10e-4
    c0.beta = 10e-2
    d0.beta = 10e-4

    a0.alpha.L.J = 10e-2
    b0.alpha.L.J = 10e-4
    c0.alpha.L.J = 10e-2
    d0.alpha.L.J = 10e-4

    a0.alpha.L.K = 10e-2
    b0.alpha.L.K = 10e-4
    c0.alpha.L.K = 10e-2
    d0.alpha.L.K = 10e-4


    if (run_covariates==T){

      covariates_y.sample = matrix(0,nrow=3,ncol=n.y)
      covariates_y.feature = matrix(0,nrow=p.y,ncol=3)
      covariates_x.sample = matrix(0,nrow=3,ncol=n.x)
      covariates_x.feature = matrix(0,nrow=p.x,ncol=3)


      if (is.null(covariates)){
        covariates_list <- list(
          covariates_y.sample = covariates_y.sample,
          covariates_y.feature = covariates_y.feature,
          covariates_x.sample = covariates_x.sample,
          covariates_x.feature = covariates_x.feature
        )
      } else {
        covariates_list <- list(
          covariates_y.sample = if(is.null(covariates[[1]])){covariates_y.sample}else{covariates[[1]]},
          covariates_y.feature = if(is.null(covariates[[2]])){covariates_y.feature}else{covariates[[2]]},
          covariates_x.sample = if(is.null(covariates[[3]])){covariates_x.sample}else{covariates[[3]]},
          covariates_x.feature = if(is.null(covariates[[4]])){covariates_x.feature}else{covariates[[4]]}
        )
      }

    }

    # Initialise parameters
    if (is.null(pivots)){
      initial.param <-initialise.gcproc(x=X.x,y=Y.y,k_dim=config$k_dim,j_dim=config$j_dim,init=config$init,verbose=config$verbose)
    }

    # Check pivoting parameters
    initial.param$pivot_y.sample <- if (is.null(pivots$pivot_y.sample)){initial.param$pivot_y.sample}else{pivots$pivot_y.sample}
    initial.param$pivot_x.sample <- if (is.null(pivots$pivot_x.sample)){initial.param$pivot_x.sample}else{pivots$pivot_x.sample}
    initial.param$pivot_y.feature <- if (is.null(pivots$pivot_y.feature)){initial.param$pivot_y.feature}else{pivots$pivot_y.feature}
    initial.param$pivot_x.feature <- if (is.null(pivots$pivot_x.feature)){initial.param$pivot_x.feature}else{pivots$pivot_x.feature}
    initial.param$y.gamma <- if (is.null(pivots$y.gamma)){(MASS::ginv((covariates_list$covariates_y.sample)%*%t(covariates_list$covariates_y.sample))%*%(covariates_list$covariates_y.sample)%*%(((Y.y))))}else{pivots$y.gamma}
    initial.param$y.delta <- if (is.null(pivots$y.delta)){(MASS::ginv(t(covariates_list$covariates_y.feature)%*%(covariates_list$covariates_y.feature))%*%t(covariates_list$covariates_y.feature)%*%((t(Y.y))))}else{pivots$y.delta}
    initial.param$x.gamma <- if (is.null(pivots$x.gamma)){(MASS::ginv((covariates_list$covariates_x.sample)%*%t(covariates_list$covariates_x.sample))%*%(covariates_list$covariates_x.sample)%*%(((X.x))))}else{pivots$x.gamma}
    initial.param$x.delta <- if (is.null(pivots$x.delta)){(MASS::ginv(t(covariates_list$covariates_x.feature)%*%(covariates_list$covariates_x.feature))%*%t(covariates_list$covariates_x.feature)%*%((t(X.x))))}else{pivots$x.delta}

    # Check anchoring parameters
    alpha.L.K.star.alpha.L.K.final <- if (is.null( anchors$anchor_y.sample)){initial.param$pivot_y.sample}else{ anchors$anchor_y.sample}
    alpha.L.J.star.alpha.L.J.final <- if (is.null( anchors$anchor_x.sample)){initial.param$pivot_x.sample}else{ anchors$anchor_x.sample}
    v.beta.star.beta.final <- if (is.null( anchors$anchor_y.feature)){initial.param$pivot_y.feature}else{ anchors$anchor_y.feature}
    u.beta.star.beta.final <- if (is.null( anchors$anchor_x.feature)){initial.param$pivot_x.feature}else{ anchors$anchor_x.feature}
    y.gamma.final <- if (is.null( anchors$anchor_y.cov.sample)){initial.param$y.gamma}else{ anchors$anchor_y.cov.sample}
    y.delta.final <- if (is.null( anchors$anchor_y.cov.feature)){initial.param$y.delta}else{ anchors$anchor_y.cov.feature}
    x.gamma.final <- if (is.null( anchors$anchor_x.cov.sample)){initial.param$x.gamma}else{ anchors$anchor_x.cov.sample}
    x.delta.final <- if (is.null( anchors$anchor_x.cov.feature)){initial.param$x.delta}else{ anchors$anchor_x.cov.feature}

    # #Initialise inverse covariance of parameters
    v.V.star.inv.beta.final = t(alpha.L.K.star.alpha.L.K.final%*%Y.y)%*%(alpha.L.K.star.alpha.L.K.final%*%Y.y)
    u.V.star.inv.beta.final = t(alpha.L.J.star.alpha.L.J.final%*%X.x)%*%(alpha.L.J.star.alpha.L.J.final%*%X.x)
    V.star.inv.alpha.L.K.final = (Y.y%*%v.beta.star.beta.final)%*%t(Y.y%*%v.beta.star.beta.final)
    V.star.inv.alpha.L.J.final = (X.x%*%u.beta.star.beta.final)%*%t(X.x%*%u.beta.star.beta.final)

    # Find intercept in endecoded space
    y_encode.final <- alpha.L.K.star.alpha.L.K.final%*%Y.y%*%v.beta.star.beta.final
    x_encode.final <- alpha.L.J.star.alpha.L.J.final%*%X.x%*%u.beta.star.beta.final

    Y_encode <- (alpha.L.K.star.alpha.L.K.final%*%(Y.y)%*%(v.beta.star.beta.final))
    X_encode <- (alpha.L.J.star.alpha.L.J.final%*%(X.x)%*%(u.beta.star.beta.final))

    Y_code <- (MASS::ginv((alpha.L.K.star.alpha.L.K.final)%*%t(alpha.L.K.star.alpha.L.K.final))%*%Y_encode%*%MASS::ginv(t(v.beta.star.beta.final)%*%(v.beta.star.beta.final)))
    X_code <- (MASS::ginv((alpha.L.J.star.alpha.L.J.final)%*%t(alpha.L.J.star.alpha.L.J.final))%*%(X_encode)%*%MASS::ginv(t(u.beta.star.beta.final)%*%(u.beta.star.beta.final)))

    Y_decoded <- t(alpha.L.K.star.alpha.L.K.final)%*%Y_code%*%t(v.beta.star.beta.final)
    X_decoded <- t(alpha.L.J.star.alpha.L.J.final)%*%X_code%*%t(u.beta.star.beta.final)


    main.parameters.copy = main.parameters = list(
      alpha.L.J.star.alpha.L.J.final = alpha.L.J.star.alpha.L.J.final,
      alpha.L.K.star.alpha.L.K.final = alpha.L.K.star.alpha.L.K.final,
      u.beta.star.beta.final = u.beta.star.beta.final,
      v.beta.star.beta.final = v.beta.star.beta.final,
      y.gamma.final = y.gamma.final,
      y.delta.final = y.delta.final,
      x.gamma.final = x.gamma.final,
      x.delta.final = x.delta.final
    )

    predict =  predict

    code = list(
      X_code = X_code,
      Y_code = Y_code,
      X_encode = X_encode,
      Y_encode = Y_encode,
      X_decoded = X_decoded,
      Y_decoded = Y_decoded
    )

  }

  if (config$verbose){
    print(paste("Beginning gcproc learning with:   Cores: ",config$cores, "   Batches: ", config$batches, "   Learning Rate (config$eta): ", config$eta, "   Tolerance Threshold: ", config$tol, "   Minimum Number of iterations: ",config$min_iter, "   Maximum number of iterations: ", config$max_iter, "   Verbose: ", config$verbose,sep=""))
  }

  while (T){

    vi.step <- variational_inference_update(
      x = X.x,
      y = Y.y,
      seed = seed,
      count = count,
      config = config,
      parameters = main.parameters,
      code = code,
      predict = predict,
      covariates = covariates_list,
      anchors = anchors)

    X.x <- vi.step$x
    Y.y <- vi.step$y
    parameters = vi.step$main.parameters
    code = vi.step$code
    predict = vi.step$predict




    # vi.step <- variational_inference_update(
    #   x = (1-a.b*predict$x)*X.x+(b.a*predict$x)*Y.y,
    #   y = Y.y,
    #   seed = seed,
    #   count = count,
    #   config = config,
    #   parameters = main.parameters,
    #   code = code,
    #   predict = predict,
    #   covariates = covariates_list,
    #   anchors = anchors)
    #
    # X.x <- vi.step$x
    # Y.y <- vi.step$y
    # parameters = vi.step$main.parameters
    # code = vi.step$code
    # predict = vi.step$predict



    matrix.residuals <- code$Y_code - code$X_code
    if (config$j_dim==1){
      llik.vec <- c(llik.vec, mean(mclust::dmvnorm(matrix.residuals,sigma = sd(matrix.residuals),log = T)))
    } else {
      llik.vec <- c(llik.vec, mean(mclust::dmvnorm(matrix.residuals,sigma = diag(diag(t(matrix.residuals)%*%(matrix.residuals)/dim(matrix.residuals)[2])),log = T)))
    }
    score.vec <- c(score.vec, (mean(abs(matrix.residuals))))
    MSE <- mean(tail(score.vec,5))
    prev.MSE <- mean(tail(score.vec,10)[1:5])


    if (config$verbose == T){
      print(paste("Iteration: ",count," with Tolerance of: ",abs(prev.MSE - MSE)," and Log-Lik of: ",tail(llik.vec,1),sep=""))
    }

    # Check convergence
    if (count>config$min_iter){
      if ((count>config$max_iter) | abs(prev.MSE - MSE) < config$tol ){
        break
      }
    }
    count = count + 1

  }

  return(list(

    main.parameters = main.parameters,

    predict =  predict,

    code = code,

    meta.parameters = config,

    convergence.parameters = list(
      iterations = count,
      score.vec = score.vec,
      llik.vec = llik.vec
    )

  ))


}




#' @export
variational_inference_update <- function(x,
                                         y,
                                         seed,
                                         count,
                                         config,
                                         parameters,
                                         code,
                                         predict,
                                         covariates_data,
                                         anchors){


  update_batched_parameters = TRUE
  epoch_run = TRUE
  check_anchors = !is.null(anchors)


  b.a <- max(1e-3,config$eta*(1/(count*10)))
  a.b <- 1-b.a

  set.seed(seed+count)

  x.g.sample <- lapply(c(1:config$batches),function(X){set.seed(X+seed+count);sample(c(1:dim(x)[1]),size = if(dim(x)[1]<config$batch_size){dim(x)[1]}else{config$batch_size})})
  y.g.sample <- lapply(c(1:config$batches),function(X){set.seed(X+seed+count);sample(c(1:dim(y)[1]),size = if(dim(y)[1]<config$batch_size){dim(y)[1]}else{config$batch_size})})
  x.v.sample <- lapply(c(1:config$batches),function(X){set.seed(X+seed+count);sample(c(1:dim(x)[2]),size = if(dim(x)[2]<config$batch_size){dim(x)[2]}else{config$batch_size})})
  y.v.sample <- lapply(c(1:config$batches),function(X){set.seed(X+seed+count);sample(c(1:dim(y)[2]),size = if(dim(y)[2]<config$batch_size){dim(y)[2]}else{config$batch_size})})

  x.c.g.sample <- lapply(c(1:config$batches),function(X){set.seed(X+seed+count);sample(c(1:dim(covariates_data$covariates_x.sample)[1]),size = if(dim(covariates_data$covariates_x.sample)[1]<config$batch_size){dim(covariates_data$covariates_x.sample)[1]}else{config$batch_size})})
  y.c.g.sample <- lapply(c(1:config$batches),function(X){set.seed(X+seed+count);sample(c(1:dim(covariates_data$covariates_y.sample)[1]),size = if(dim(covariates_data$covariates_y.sample)[1]<config$batch_size){dim(covariates_data$covariates_y.sample)[1]}else{config$batch_size})})
  x.c.v.sample <- lapply(c(1:config$batches),function(X){set.seed(X+seed+count);sample(c(1:dim(covariates_data$covariates_x.feature)[2]),size = if(dim(covariates_data$covariates_x.feature)[2]<config$batch_size){dim(covariates_data$covariates_x.feature)[2]}else{config$batch_size})})
  y.c.v.sample <- lapply(c(1:config$batches),function(X){set.seed(X+seed+count);sample(c(1:dim(covariates_data$covariates_y.feature)[2]),size = if(dim(covariates_data$covariates_y.feature)[2]<config$batch_size){dim(covariates_data$covariates_y.feature)[2]}else{config$batch_size})})

  internal_list <- list(
    intercept.x = 0,
    y_encode = 0,
    x_encode = 0,
    alpha.L.J.star.alpha.L.J = 0,
    alpha.L.K.star.alpha.L.K = 0,
    V.star.inv.alpha.L.J = 0,
    V.star.inv.alpha.L.K = 0,
    v.beta.star.beta = 0,
    u.beta.star.beta = 0,
    u.V.star.inv.beta = 0,
    v.V.star.inv.beta = 0,
    y.gamma = 0,
    y.delta = 0,
    x.gamma = 0,
    x.delta = 0,
    y.g.ids = 0,
    x.g.ids = 0,
    x.v.ids = 0,
    y.v.ids = 0
  )


  to_return <- parallel::mclapply(c(1:config$batches),function(i){


    x.g.ids <- x.g.sample[[i]]
    y.g.ids <- y.g.sample[[i]]

    x.v.ids <- x.v.sample[[i]]
    y.v.ids <- y.v.sample[[i]]



    x.c.g.ids <- x.c.g.sample[[i]]
    y.c.g.ids <- y.c.g.sample[[i]]

    x.c.v.ids <- x.c.v.sample[[i]]
    y.c.v.ids <- y.c.v.sample[[i]]


    alpha.L.J.star.alpha.L.J <- parameters$alpha.L.J.star.alpha.L.J.final[,x.g.ids]
    alpha.L.K.star.alpha.L.K <- parameters$alpha.L.K.star.alpha.L.K.final[,y.g.ids]

    u.beta.star.beta <- parameters$u.beta.star.beta.final[x.v.ids,]
    v.beta.star.beta <- parameters$v.beta.star.beta.final[y.v.ids,]

    V.star.inv.alpha.L.J <- parameters$V.star.inv.alpha.L.J.final[x.g.ids,x.g.ids]
    V.star.inv.alpha.L.K <- parameters$V.star.inv.alpha.L.K.final[y.g.ids,y.g.ids]

    v.V.star.inv.beta <- parameters$v.V.star.inv.beta.final[y.v.ids,y.v.ids]
    u.V.star.inv.beta <- parameters$u.V.star.inv.beta.final[x.v.ids,x.v.ids]

    y <- (y)[y.g.ids,y.v.ids]
    x <- (x)[x.g.ids,x.v.ids]

    cov.y.s <- covariates_data$covariates_y.sample[y.c.g.ids,y.g.ids]
    cov.y.f <- covariates_data$covariates_y.feature[y.v.ids,y.c.v.ids]
    cov.x.s <- covariates_data$covariates_x.sample[x.c.g.ids,x.g.ids]
    cov.x.f <- covariates_data$covariates_x.feature[x.v.ids,x.c.v.ids]

    cov.y.s.p <- (t(covariates_data$covariates_y.sample)[y.g.ids,y.c.g.ids]%*%parameters$y.gamma.final[y.c.g.ids,y.v.ids])
    cov.x.s.p <- (t(covariates_data$covariates_x.sample)[x.g.ids,x.c.g.ids]%*%parameters$x.gamma.final[x.c.g.ids,x.v.ids])

    cov.y.f.p <- t(covariates_data$covariates_y.feature[y.v.ids,y.c.v.ids]%*%parameters$y.delta.final[y.c.v.ids,y.g.ids])
    cov.x.f.p <- t(covariates_data$covariates_x.feature[x.v.ids,x.c.v.ids]%*%parameters$x.delta.final[x.c.v.ids,x.g.ids])

    y_encode <- code$Y_encode
    x_encode <- code$X_encode

    V.star.inv.alpha.L.J = (( x - cov.x.f.p - cov.x.s.p)%*%u.beta.star.beta)%*%t(( x - cov.x.f.p - cov.x.s.p)%*%u.beta.star.beta)
    alpha.L.J.star.alpha.L.J = if (is.null(anchors$anchor_x.sample)){(y_encode)%*%t(( x - cov.x.f.p - cov.x.s.p)%*%u.beta.star.beta)%*%MASS::ginv(V.star.inv.alpha.L.J)}else{anchors$anchor_x.sample[,x.g.ids]}

    V.star.inv.alpha.L.K = ((y - cov.y.f.p - cov.y.s.p)%*%v.beta.star.beta)%*%t((y - cov.y.f.p - cov.y.s.p)%*%v.beta.star.beta)
    alpha.L.K.star.alpha.L.K = if (is.null(anchors$anchor_y.sample)){((x_encode ))%*%t((y - cov.y.f.p - cov.y.s.p)%*%v.beta.star.beta)%*%MASS::ginv(V.star.inv.alpha.L.K)}else{anchors$anchor_y.sample[,y.g.ids]}

    u.V.star.inv.beta = t(alpha.L.J.star.alpha.L.J%*%( x - cov.x.f.p - cov.x.s.p))%*%(alpha.L.J.star.alpha.L.J%*%( x - cov.x.f.p - cov.x.s.p))
    u.beta.star.beta = if (is.null(anchors$anchor_x.feature)){MASS::ginv(u.V.star.inv.beta)%*%t(alpha.L.J.star.alpha.L.J%*%( x - cov.x.f.p - cov.x.s.p))%*%(y_encode)}else{anchors$anchor_x.feature[x.v.ids,]}

    v.V.star.inv.beta = t(alpha.L.K.star.alpha.L.K%*%(y - cov.y.f.p - cov.y.s.p))%*%(alpha.L.K.star.alpha.L.K%*%(y - cov.y.f.p - cov.y.s.p))
    v.beta.star.beta = if (is.null(anchors$anchor_y.feature)){MASS::ginv(v.V.star.inv.beta)%*%t(alpha.L.K.star.alpha.L.K%*%(y - cov.y.f.p - cov.y.s.p))%*%((x_encode ))}else{anchors$anchor_y.feature[y.v.ids,]}

    y.gamma <- if(is.null(anchors$anchor_y.cov.sample)){MASS::ginv((cov.y.s)%*%t(cov.y.s))%*%(cov.y.s)%*%(((MASS::ginv(t(alpha.L.K.star.alpha.L.K)%*%(alpha.L.K.star.alpha.L.K))%*%t(alpha.L.K.star.alpha.L.K)%*%((alpha.L.K.star.alpha.L.K %*% (y) %*% v.beta.star.beta - (alpha.L.K.star.alpha.L.K%*%cov.y.f.p%*%v.beta.star.beta)  - (x_encode )))%*%t(v.beta.star.beta)%*%MASS::ginv((v.beta.star.beta)%*%t(v.beta.star.beta)))))}else{anchors$anchor_y.cov.sample}
    y.delta <- if(is.null(anchors$anchor_y.cov.feature)){MASS::ginv(t(cov.y.f)%*%(cov.y.f))%*%t(cov.y.f)%*%((t((MASS::ginv(t(alpha.L.K.star.alpha.L.K)%*%(alpha.L.K.star.alpha.L.K))%*%t(alpha.L.K.star.alpha.L.K)%*%(alpha.L.K.star.alpha.L.K %*% (y) %*% v.beta.star.beta - (alpha.L.K.star.alpha.L.K%*%cov.y.s.p%*%v.beta.star.beta) - (x_encode )))%*%t(v.beta.star.beta)%*%MASS::ginv((v.beta.star.beta)%*%t(v.beta.star.beta)))))}else{anchors$anchor_y.cov.feature}

    x.gamma <- if(is.null(anchors$anchor_x.cov.sample)){MASS::ginv((cov.x.s)%*%t(cov.x.s))%*%(cov.x.s)%*%(((MASS::ginv(t(alpha.L.J.star.alpha.L.J)%*%(alpha.L.J.star.alpha.L.J))%*%t(alpha.L.J.star.alpha.L.J)%*%((alpha.L.J.star.alpha.L.J%*%(x)%*%u.beta.star.beta - (alpha.L.J.star.alpha.L.J%*%cov.x.f.p%*%u.beta.star.beta)) - y_encode)%*%t(u.beta.star.beta)%*%MASS::ginv((u.beta.star.beta)%*%t(u.beta.star.beta)))))}else{anchors$anchor_x.cov.sample}
    x.delta <- if(is.null(anchors$anchor_x.cov.feature)){MASS::ginv(t(cov.x.f)%*%(cov.x.f))%*%t(cov.x.f)%*%((t((MASS::ginv(t(alpha.L.J.star.alpha.L.J)%*%(alpha.L.J.star.alpha.L.J))%*%t(alpha.L.J.star.alpha.L.J)%*%(alpha.L.J.star.alpha.L.J%*%(x)%*%u.beta.star.beta - (alpha.L.J.star.alpha.L.J%*%cov.x.s.p%*%u.beta.star.beta) - y_encode))%*%t(u.beta.star.beta)%*%MASS::ginv((u.beta.star.beta)%*%t(u.beta.star.beta)))))}else{anchors$anchor_x.cov.feature}

    internal_list$alpha.L.J.star.alpha.L.J = alpha.L.J.star.alpha.L.J
    internal_list$alpha.L.K.star.alpha.L.K = alpha.L.K.star.alpha.L.K
    internal_list$V.star.inv.alpha.L.J = V.star.inv.alpha.L.J
    internal_list$V.star.inv.alpha.L.K = V.star.inv.alpha.L.K
    internal_list$v.beta.star.beta = v.beta.star.beta
    internal_list$u.beta.star.beta = u.beta.star.beta
    internal_list$u.V.star.inv.beta = u.V.star.inv.beta
    internal_list$v.V.star.inv.beta = v.V.star.inv.beta
    internal_list$y.gamma = y.gamma
    internal_list$y.delta = y.delta
    internal_list$x.gamma = x.gamma
    internal_list$x.delta = x.delta
    internal_list$y.g.ids = y.g.ids
    internal_list$x.g.ids = x.g.ids
    internal_list$x.v.ids = x.v.ids
    internal_list$y.v.ids = y.v.ids
    internal_list$y.c.g.ids = y.c.g.ids
    internal_list$x.c.g.ids = x.c.g.ids
    internal_list$x.c.v.ids = x.c.v.ids
    internal_list$y.c.v.ids = y.c.v.ids

    return(internal_list)
  },mc.silent = config$verbose,mc.cores = config$cores)


  if (update_batched_parameters==T){

    for (i in 1:config$batches){

      if (is.character(to_return[[i]])){
        to_return[[i]] <- internal_list
      }

      alpha.L.J.star.alpha.L.J <- to_return[[i]]$alpha.L.J.star.alpha.L.J
      alpha.L.K.star.alpha.L.K <- to_return[[i]]$alpha.L.K.star.alpha.L.K
      V.star.inv.alpha.L.J <- to_return[[i]]$V.star.inv.alpha.L.J
      V.star.inv.alpha.L.K <- to_return[[i]]$V.star.inv.alpha.L.K

      v.beta.star.beta <- to_return[[i]]$v.beta.star.beta
      u.beta.star.beta <- to_return[[i]]$u.beta.star.beta
      u.V.star.inv.beta <- to_return[[i]]$u.V.star.inv.beta
      v.V.star.inv.beta <- to_return[[i]]$v.V.star.inv.beta

      y.gamma <- to_return[[i]]$y.gamma
      y.delta <- to_return[[i]]$y.delta
      x.gamma <- to_return[[i]]$x.gamma
      x.delta <- to_return[[i]]$x.delta

      x.g.ids <- to_return[[i]]$x.g.ids
      y.g.ids <- to_return[[i]]$y.g.ids

      x.v.ids <- to_return[[i]]$x.v.ids
      y.v.ids <- to_return[[i]]$y.v.ids


      x.c.g.ids <- to_return[[i]]$x.c.g.ids
      y.c.g.ids <- to_return[[i]]$y.c.g.ids

      x.c.v.ids <- to_return[[i]]$x.c.v.ids
      y.c.v.ids <- to_return[[i]]$y.c.v.ids



      parameters$alpha.L.J.star.alpha.L.J.final[,x.g.ids] <- a.b*parameters$alpha.L.J.star.alpha.L.J.final[,x.g.ids] + b.a*alpha.L.J.star.alpha.L.J
      parameters$V.star.inv.alpha.L.J.final[x.g.ids,x.g.ids] <- a.b*parameters$V.star.inv.alpha.L.J.final[x.g.ids,x.g.ids] + b.a*V.star.inv.alpha.L.J

      parameters$alpha.L.K.star.alpha.L.K.final[,y.g.ids] <- a.b*parameters$alpha.L.K.star.alpha.L.K.final[,y.g.ids] + b.a*alpha.L.K.star.alpha.L.K
      parameters$V.star.inv.alpha.L.K.final[y.g.ids,y.g.ids] <- a.b*parameters$V.star.inv.alpha.L.K.final[y.g.ids,y.g.ids] + b.a*V.star.inv.alpha.L.K

      parameters$u.beta.star.beta.final[x.v.ids,] <- a.b*parameters$u.beta.star.beta.final[x.v.ids,] + b.a*u.beta.star.beta
      parameters$u.V.star.inv.beta.final[x.v.ids,x.v.ids] <- a.b*parameters$u.V.star.inv.beta.final[x.v.ids,x.v.ids] + b.a*u.V.star.inv.beta

      parameters$v.beta.star.beta.final[y.v.ids,] <- a.b*parameters$v.beta.star.beta.final[y.v.ids,] + b.a*v.beta.star.beta
      parameters$v.V.star.inv.beta.final[y.v.ids,y.v.ids] <- a.b*parameters$v.V.star.inv.beta.final[y.v.ids,y.v.ids] + b.a*v.V.star.inv.beta

      parameters$y.gamma.final[y.c.g.ids,y.v.ids] <- a.b*parameters$y.gamma.final[y.c.g.ids,y.v.ids] + b.a*y.gamma
      parameters$y.delta.final[y.c.v.ids,y.g.ids] <- a.b*parameters$y.delta.final[y.c.v.ids,y.g.ids] + b.a*y.delta
      parameters$x.gamma.final[x.c.g.ids,x.v.ids] <- a.b*parameters$x.gamma.final[x.c.g.ids,x.v.ids] + b.a*x.gamma
      parameters$x.delta.final[x.c.v.ids,x.g.ids] <- a.b*parameters$x.delta.final[x.c.v.ids,x.g.ids] + b.a*x.delta


    }

  }

  if (check_anchors==T){

    parameters$alpha.L.J.star.alpha.L.J.final <- if (is.null(anchors$anchor_x.sample)){parameters$alpha.L.J.star.alpha.L.J.final}else{anchors$anchor_x.sample}
    parameters$alpha.L.K.star.alpha.L.K.final <- if (is.null(anchors$anchor_y.sample)){parameters$alpha.L.K.star.alpha.L.K.final}else{anchors$anchor_y.sample}

    parameters$u.beta.star.beta.final <- if (is.null(anchors$anchor_x.feature)){parameters$u.beta.star.beta.final}else{anchors$anchor_x.feature}
    parameters$v.beta.star.beta.final <- if (is.null(anchors$anchor_y.feature)){parameters$v.beta.star.beta.final}else{anchors$anchor_y.feature}

  }


  if (epoch_run==T){

    cov.y.s.p <- (t(covariates_data$covariates_y.sample)%*%parameters$y.gamma.final)
    cov.x.s.p <- (t(covariates_data$covariates_x.sample)%*%parameters$x.gamma.final)

    cov.y.f.p <- t(covariates_data$covariates_y.feature%*%parameters$y.delta.final)
    cov.x.f.p <- t(covariates_data$covariates_x.feature%*%parameters$x.delta.final)


    y_encode <- (parameters$alpha.L.K.star.alpha.L.K.final%*%( y )%*%(parameters$v.beta.star.beta.final))
    x_encode <- (parameters$alpha.L.J.star.alpha.L.J.final%*%( x )%*%(parameters$u.beta.star.beta.final))

    cov.y_encode <- (parameters$alpha.L.K.star.alpha.L.K.final%*%( (cov.y.s.p + cov.y.f.p) )%*%(parameters$v.beta.star.beta.final))
    cov.x_encode <- (parameters$alpha.L.J.star.alpha.L.J.final%*%( (cov.x.s.p + cov.x.f.p) )%*%(parameters$u.beta.star.beta.final))

    code$Y_encode <- y_encode - cov.y_encode
    code$X_encode <- x_encode - cov.x_encode

    code$Y_code <- (MASS::ginv((parameters$alpha.L.K.star.alpha.L.K.final)%*%t(parameters$alpha.L.K.star.alpha.L.K.final))%*%(y_encode)%*%MASS::ginv(t(parameters$v.beta.star.beta.final)%*%(parameters$v.beta.star.beta.final)))
    code$X_code <- (MASS::ginv((parameters$alpha.L.J.star.alpha.L.J.final)%*%t(parameters$alpha.L.J.star.alpha.L.J.final))%*%(x_encode)%*%MASS::ginv(t(parameters$u.beta.star.beta.final)%*%(parameters$u.beta.star.beta.final)))

    code$Y_decoded <- t(parameters$alpha.L.K.star.alpha.L.K.final)%*%(code$Y_code)%*%t(parameters$v.beta.star.beta.final)
    code$X_decoded <- (t(parameters$alpha.L.J.star.alpha.L.J.final)%*%(code$X_code)%*%t(parameters$u.beta.star.beta.final))

    if (!is.null(predict$y)){
      predict$predict.y <- code$Y_decoded
      y[which(predict$y==1,arr.ind = T)] <- predict$predict.y[which(predict$y==1,arr.ind = T)]
    }
    if (!is.null(predict$x)){
      predict$predict.x <- code$X_decoded
      x[which(predict$x==1,arr.ind = T)] <- predict$predict.x[which(predict$x==1,arr.ind = T)]
    }

  }



  main.parameters = parameters

  predict =  predict

  code = code

  final_vi.update <- list(
    x = x,
    y = y,
    main.parameters = main.parameters,
    predict = predict,
    code = code,
    iteration = count
  )

  return(final_vi.update)

}
