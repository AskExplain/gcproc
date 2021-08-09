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
#' @return intercept.encode:   Intercept for the encoding of x (where KYv = JXu + intercept)
#' @export
gcproc <- function(x,
                   y,
                   config = list(k_dim = 70,
                                 j_dim = 70,
                                 eta=5e-3,
                                 max_iter=1500,
                                 min_iter = 15,
                                 tol=1e-3,
                                 log=F,
                                 center=T,
                                 scale.z=T,
                                 batches=16,
                                 cores=2,
                                 verbose=T,
                                 init="svd-quick"),
                   seed = 1,
                   anchors = NULL,
                   pivots = NULL,
                   covariates = NULL,
                   predict = NULL

){

  x <- as.matrix(x)
  y <- as.matrix(y)


  # Step sizes
  b.a <- config$eta
  a.b <- c(1-b.a)


  config$batches <- ceiling(max(c(dim(x),dim(y)))/128)

  prepare_data = TRUE
  initialise = TRUE
  variational_gradient_descent_updates = TRUE
  update_batched_parameters = TRUE
  epoch_run = TRUE

  run_covariates = TRUE
  check_anchors = FALSE



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

    anchor_y.sample = anchors$anchor_y.sample
    anchor_y.feature = anchors$anchor_y.feature
    anchor_x.sample = anchors$anchor_x.sample
    anchor_x.feature = anchors$anchor_x.feature
    anchor_y.cov.sample = anchors$anchor_y.cov.sample
    anchor_y.cov.feature = anchors$anchor_y.cov.feature
    anchor_x.cov.sample = anchors$anchor_x.cov.sample
    anchor_x.cov.feature = anchors$anchor_x.cov.feature

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
      print(paste("Standardise by standard deviation: ", config$scale.z, sep = ""))
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
    best.predict_MSE <- sum(abs(y)) + sum(abs(x))

    if (run_covariates==T){

      covariates_y.sample = matrix(0,nrow=3,ncol=dim(y)[1])
      covariates_y.feature = matrix(0,nrow=dim(y)[2],ncol=3)
      covariates_x.sample = matrix(0,nrow=3,ncol=dim(x)[1])
      covariates_x.feature = matrix(0,nrow=dim(x)[2],ncol=3)


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
    alpha.L.K.star.alpha.L.K.final <- if (is.null(anchor_y.sample)){initial.param$pivot_y.sample}else{anchor_y.sample}
    alpha.L.J.star.alpha.L.J.final <- if (is.null(anchor_x.sample)){initial.param$pivot_x.sample}else{anchor_x.sample}
    v.beta.star.beta.final <- if (is.null(anchor_y.feature)){initial.param$pivot_y.feature}else{anchor_y.feature}
    u.beta.star.beta.final <- if (is.null(anchor_x.feature)){initial.param$pivot_x.feature}else{anchor_x.feature}
    y.gamma.final <- if (is.null(anchor_y.cov.sample)){initial.param$y.gamma}else{anchor_y.cov.sample}
    y.delta.final <- if (is.null(anchor_y.cov.feature)){initial.param$y.delta}else{anchor_y.cov.feature}
    x.gamma.final <- if (is.null(anchor_x.cov.sample)){initial.param$x.gamma}else{anchor_x.cov.sample}
    x.delta.final <- if (is.null(anchor_x.cov.feature)){initial.param$x.delta}else{anchor_x.cov.feature}

    # #Initialise inverse covariance of parameters
    v.V.star.inv.beta.final = t(alpha.L.K.star.alpha.L.K.final%*%Y.y)%*%(alpha.L.K.star.alpha.L.K.final%*%Y.y)
    u.V.star.inv.beta.final = t(alpha.L.J.star.alpha.L.J.final%*%X.x)%*%(alpha.L.J.star.alpha.L.J.final%*%X.x)
    V.star.inv.alpha.L.K.final = (Y.y%*%v.beta.star.beta.final)%*%t(Y.y%*%v.beta.star.beta.final)
    V.star.inv.alpha.L.J.final = (X.x%*%u.beta.star.beta.final)%*%t(X.x%*%u.beta.star.beta.final)

    # Find intercept in encoded space
    Y_encode <- alpha.L.K.star.alpha.L.K.final%*%Y.y%*%v.beta.star.beta.final
    X_encode <- alpha.L.J.star.alpha.L.J.final%*%X.x%*%u.beta.star.beta.final
    intercept.encode <- Y_encode - X_encode

    Y_code <- (MASS::ginv((alpha.L.K.star.alpha.L.K.final)%*%t(alpha.L.K.star.alpha.L.K.final))%*%alpha.L.K.star.alpha.L.K.final%*%(Y.y)%*%(v.beta.star.beta.final)%*%MASS::ginv(t(v.beta.star.beta.final)%*%(v.beta.star.beta.final)))
    X_code <- (MASS::ginv((alpha.L.J.star.alpha.L.J.final)%*%t(alpha.L.J.star.alpha.L.J.final))%*%(alpha.L.J.star.alpha.L.J.final%*%(X.x)%*%(u.beta.star.beta.final))%*%MASS::ginv(t(u.beta.star.beta.final)%*%(u.beta.star.beta.final)))
    intercept.code <- Y_code - X_code
  }

  if (config$verbose){
    print(paste("Beginning gcproc learning with:   Cores: ",config$cores, "   Batches: ", config$batches, "   Learning Rate (config$eta): ", config$eta, "   Tolerance Threshold: ", config$tol, "   Minimum Number of iterations: ",config$min_iter, "   Maximum number of iterations: ", config$max_iter, "   Verbose: ", config$verbose,sep=""))
  }


  while (T){


    internal_list <- list(
      intercept.encode = 0,
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



    if (variational_gradient_descent_updates == T){
      set.seed(seed+count)

      x.g.sample <- lapply(c(1:config$batches),function(X){set.seed(X+seed+count);sample(c(1:dim(X.x)[1]),size = if(dim(X.x)[1]<128){dim(X.x)[1]}else{128})})
      y.g.sample <- lapply(c(1:config$batches),function(X){set.seed(X+seed+count);sample(c(1:dim(Y.y)[1]),size = if(dim(Y.y)[1]<128){dim(Y.y)[1]}else{128})})
      x.v.sample <- lapply(c(1:config$batches),function(X){set.seed(X+seed+count);sample(c(1:dim(X.x)[2]),size = if(dim(X.x)[2]<128){dim(X.x)[2]}else{128})})
      y.v.sample <- lapply(c(1:config$batches),function(X){set.seed(X+seed+count);sample(c(1:dim(Y.y)[2]),size = if(dim(Y.y)[2]<128){dim(Y.y)[2]}else{128})})

      to_return <- parallel::mclapply(c(1:config$batches),function(i){

        x.g.ids <- x.g.sample[[i]]
        y.g.ids <- y.g.sample[[i]]

        x.v.ids <- x.v.sample[[i]]
        y.v.ids <- y.v.sample[[i]]

        alpha.L.J.star.alpha.L.J <- alpha.L.J.star.alpha.L.J.final[,x.g.ids]
        alpha.L.K.star.alpha.L.K <- alpha.L.K.star.alpha.L.K.final[,y.g.ids]

        u.beta.star.beta <- u.beta.star.beta.final[x.v.ids,]
        v.beta.star.beta <- v.beta.star.beta.final[y.v.ids,]

        V.star.inv.alpha.L.J <- V.star.inv.alpha.L.J.final[x.g.ids,x.g.ids]
        V.star.inv.alpha.L.K <- V.star.inv.alpha.L.K.final[y.g.ids,y.g.ids]

        v.V.star.inv.beta <- v.V.star.inv.beta.final[y.v.ids,y.v.ids]
        u.V.star.inv.beta <- u.V.star.inv.beta.final[x.v.ids,x.v.ids]

        cov.y.s <- covariates_list$covariates_y.sample[,y.g.ids]
        cov.y.f <- covariates_list$covariates_y.feature[y.v.ids,]
        cov.x.s <- covariates_list$covariates_x.sample[,x.g.ids]
        cov.x.f <- covariates_list$covariates_x.feature[x.v.ids,]

        cov.y.s.p <- (t(covariates_list$covariates_y.sample)[y.g.ids,]%*%y.gamma.final[,y.v.ids])
        cov.x.s.p <- (t(covariates_list$covariates_x.sample)[x.g.ids,]%*%x.gamma.final[,x.v.ids])

        cov.y.f.p <- t(covariates_list$covariates_y.feature[y.v.ids,]%*%y.delta.final[,y.g.ids])
        cov.x.f.p <- t(covariates_list$covariates_x.feature[x.v.ids,]%*%x.delta.final[,x.g.ids])

        y <- (Y.y)[y.g.ids,y.v.ids]
        x <- (X.x)[x.g.ids,x.v.ids]

        y_encode <- Y_encode
        x_encode <- X_encode

        V.star.inv.alpha.L.J = ((x + cov.x.f.p + cov.x.s.p)%*%u.beta.star.beta)%*%t((x + cov.x.f.p + cov.x.s.p)%*%u.beta.star.beta)
        alpha.L.J.star.alpha.L.J = if (is.null(anchor_x.sample)){(y_encode)%*%t((x + cov.x.f.p + cov.x.s.p)%*%u.beta.star.beta)%*%MASS::ginv(V.star.inv.alpha.L.J)}else{anchor_x.sample[,x.g.ids]}

        V.star.inv.alpha.L.K =  ((y + cov.y.f.p + cov.y.s.p)%*%v.beta.star.beta)%*%t((y + cov.y.f.p + cov.y.s.p)%*%v.beta.star.beta)
        alpha.L.K.star.alpha.L.K = if (is.null(anchor_y.sample)){((x_encode + intercept.encode))%*%t((y + cov.y.f.p + cov.y.s.p)%*%v.beta.star.beta)%*%MASS::ginv(V.star.inv.alpha.L.K)}else{anchor_y.sample[,y.g.ids]}


        u.V.star.inv.beta = t(alpha.L.J.star.alpha.L.J%*%(x + cov.x.s.p + cov.x.f.p))%*%(alpha.L.J.star.alpha.L.J%*%(x + cov.x.s.p + cov.x.f.p))
        u.beta.star.beta = if (is.null(anchor_x.feature)){MASS::ginv(u.V.star.inv.beta)%*%t(alpha.L.J.star.alpha.L.J%*%(x + cov.x.s.p + cov.x.f.p))%*%(y_encode)}else{anchor_x.feature[x.v.ids,]}


        v.V.star.inv.beta = t(alpha.L.K.star.alpha.L.K%*%(y + cov.y.s.p + cov.y.f.p))%*%(alpha.L.K.star.alpha.L.K%*%(y + cov.y.s.p + cov.y.f.p))
        v.beta.star.beta = if (is.null(anchor_y.feature)){MASS::ginv(v.V.star.inv.beta)%*%t(alpha.L.K.star.alpha.L.K%*%(y + cov.y.s.p + cov.y.f.p))%*%((x_encode + intercept.encode))}else{anchor_y.feature[y.v.ids,]}


        y.gamma <- if(is.null(anchor_y.cov.sample)){MASS::ginv((cov.y.s)%*%t(cov.y.s))%*%(cov.y.s)%*%(((MASS::ginv(t(alpha.L.K.star.alpha.L.K)%*%(alpha.L.K.star.alpha.L.K))%*%t(alpha.L.K.star.alpha.L.K)%*%((alpha.L.K.star.alpha.L.K %*% (y) %*% v.beta.star.beta + (alpha.L.K.star.alpha.L.K%*%cov.y.f.p%*%v.beta.star.beta)  - (x_encode + intercept.encode)))%*%t(v.beta.star.beta)%*%MASS::ginv((v.beta.star.beta)%*%t(v.beta.star.beta)))))}else{anchor_y.cov.sample}
        y.delta <- if(is.null(anchor_y.cov.feature)){MASS::ginv(t(cov.y.f)%*%(cov.y.f))%*%t(cov.y.f)%*%((t((MASS::ginv(t(alpha.L.K.star.alpha.L.K)%*%(alpha.L.K.star.alpha.L.K))%*%t(alpha.L.K.star.alpha.L.K)%*%(alpha.L.K.star.alpha.L.K %*% (y) %*% v.beta.star.beta + (alpha.L.K.star.alpha.L.K%*%cov.y.s.p%*%v.beta.star.beta) - (x_encode + intercept.encode)))%*%t(v.beta.star.beta)%*%MASS::ginv((v.beta.star.beta)%*%t(v.beta.star.beta)))))}else{anchor_y.cov.feature}

        x.gamma <- if(is.null(anchor_x.cov.sample)){MASS::ginv((cov.x.s)%*%t(cov.x.s))%*%(cov.x.s)%*%(((MASS::ginv(t(alpha.L.J.star.alpha.L.J)%*%(alpha.L.J.star.alpha.L.J))%*%t(alpha.L.J.star.alpha.L.J)%*%((alpha.L.J.star.alpha.L.J%*%(x)%*%u.beta.star.beta + (alpha.L.J.star.alpha.L.J%*%cov.x.f.p%*%u.beta.star.beta)) - y_encode)%*%t(u.beta.star.beta)%*%MASS::ginv((u.beta.star.beta)%*%t(u.beta.star.beta)))))}else{anchor_x.cov.sample}
        x.delta <- if(is.null(anchor_x.cov.feature)){MASS::ginv(t(cov.x.f)%*%(cov.x.f))%*%t(cov.x.f)%*%((t((MASS::ginv(t(alpha.L.J.star.alpha.L.J)%*%(alpha.L.J.star.alpha.L.J))%*%t(alpha.L.J.star.alpha.L.J)%*%(alpha.L.J.star.alpha.L.J%*%(x)%*%u.beta.star.beta + (alpha.L.J.star.alpha.L.J%*%cov.x.s.p%*%u.beta.star.beta) - y_encode))%*%t(u.beta.star.beta)%*%MASS::ginv((u.beta.star.beta)%*%t(u.beta.star.beta)))))}else{anchor_x.cov.feature}


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

        return(internal_list)
      },mc.silent = config$verbose,mc.cores = config$cores)


    }



    if (update_batched_parameters==T){

      y.gamma.internal <- 0
      y.delta.internal <- 0
      x.gamma.internal <- 0
      x.delta.internal <- 0

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

        alpha.L.J.star.alpha.L.J.final[,x.g.ids] <- a.b*alpha.L.J.star.alpha.L.J.final[,x.g.ids] + b.a*alpha.L.J.star.alpha.L.J
        V.star.inv.alpha.L.J.final[x.g.ids,x.g.ids] <- a.b*V.star.inv.alpha.L.J.final[x.g.ids,x.g.ids] + b.a*V.star.inv.alpha.L.J

        alpha.L.K.star.alpha.L.K.final[,y.g.ids] <- a.b*alpha.L.K.star.alpha.L.K.final[,y.g.ids] + b.a*alpha.L.K.star.alpha.L.K
        V.star.inv.alpha.L.K.final[y.g.ids,y.g.ids] <- a.b*V.star.inv.alpha.L.K.final[y.g.ids,y.g.ids] + b.a*V.star.inv.alpha.L.K

        u.beta.star.beta.final[x.v.ids,] <- a.b*u.beta.star.beta.final[x.v.ids,] + b.a*u.beta.star.beta
        u.V.star.inv.beta.final[x.v.ids,x.v.ids] <- a.b*u.V.star.inv.beta.final[x.v.ids,x.v.ids] + b.a*u.V.star.inv.beta

        v.beta.star.beta.final[y.v.ids,] <- a.b*v.beta.star.beta.final[y.v.ids,] + b.a*v.beta.star.beta
        v.V.star.inv.beta.final[y.v.ids,y.v.ids] <- a.b*v.V.star.inv.beta.final[y.v.ids,y.v.ids] + b.a*v.V.star.inv.beta

        y.gamma.final[,y.v.ids] <- a.b*y.gamma.final[,y.v.ids] + b.a*y.gamma
        y.delta.final[,y.g.ids] <- a.b*y.delta.final[,y.g.ids] + b.a*y.delta
        x.gamma.final[,x.v.ids] <- a.b*x.gamma.final[,x.v.ids] + b.a*x.gamma
        x.delta.final[,x.g.ids] <- a.b*x.delta.final[,x.g.ids] + b.a*x.delta

      }

    }


    if (epoch_run==T){

      cov.y.s.p <- (t(covariates_list$covariates_y.sample)%*%y.gamma.final)
      cov.x.s.p <- (t(covariates_list$covariates_x.sample)%*%x.gamma.final)

      cov.y.f.p <- t(covariates_list$covariates_y.feature%*%y.delta.final)
      cov.x.f.p <- t(covariates_list$covariates_x.feature%*%x.delta.final)

      Y_encode <- (alpha.L.K.star.alpha.L.K.final%*%(y + cov.y.s.p + cov.y.f.p)%*%(v.beta.star.beta.final))
      Y_code <- (MASS::ginv((alpha.L.K.star.alpha.L.K.final)%*%t(alpha.L.K.star.alpha.L.K.final))%*%Y_encode%*%MASS::ginv(t(v.beta.star.beta.final)%*%(v.beta.star.beta.final)))

      X_encode <- (alpha.L.J.star.alpha.L.J.final%*%(x + cov.x.s.p + cov.x.f.p)%*%(u.beta.star.beta.final))
      X_code <- (MASS::ginv((alpha.L.J.star.alpha.L.J.final)%*%t(alpha.L.J.star.alpha.L.J.final))%*%(X_encode+intercept.encode)%*%MASS::ginv(t(u.beta.star.beta.final)%*%(u.beta.star.beta.final)))

      # intercept.code = Y_code - X_code
      # intercept.encode = Y_encode - X_encode

      intercept.code <- intercept.encode <- 0
    }



    if (check_anchors==T){

      alpha.L.J.star.alpha.L.J.final <- if (is.null(anchor_x.sample)){alpha.L.J.star.alpha.L.J.final}else{anchor_x.sample}
      alpha.L.K.star.alpha.L.K.final <- if (is.null(anchor_y.sample)){alpha.L.K.star.alpha.L.K.final}else{anchor_y.sample}

      u.beta.star.beta.final <- if (is.null(anchor_x.feature)){u.beta.star.beta.final}else{anchor_x.feature}
      v.beta.star.beta.final <- if (is.null(anchor_y.feature)){v.beta.star.beta.final}else{anchor_y.feature}

    }

    if (!is.null(predict) & count > 100 ){

      predictions.x <- (t(alpha.L.J.star.alpha.L.J.final)%*%((X_code+intercept.code))%*%t(u.beta.star.beta.final))
      predictions.y <- (t(alpha.L.K.star.alpha.L.K.final)%*%((Y_code))%*%t(v.beta.star.beta.final))

      main.prediction <- (predictions.x+predictions.y)/2

      predict_MSE <- mean(abs((((x+y)/2 - main.prediction)*(predict$x==0))))


      if (best.predict_MSE*Inf > predict_MSE){
        best.predict_MSE <- predict_MSE

        if (!is.null(predict$x)){
          x[which(predict$x==1,arr.ind = T)] <- main.prediction[which(predict$x==1,arr.ind = T)]
        }
        predict$predictions.x <- main.prediction

        if (!is.null(predict$y)){
          y[which(predict$y==1,arr.ind = T)] <- main.prediction[which(predict$y==1,arr.ind = T)]
        }
        predict$predictions.y <- main.prediction
      }


    }


    matrix.residuals <- Y_code - X_code
    llik.vec <- c(llik.vec, mean(mclust::dmvnorm(matrix.residuals,sigma = diag(diag(t(matrix.residuals)%*%(matrix.residuals)/dim(matrix.residuals)[2])),log = T)))
    score.vec <- c(score.vec, (mean(abs(matrix.residuals))))
    MSE <- mean(tail(score.vec,2))
    prev_MSE <- mean(tail(score.vec,4)[1:2])

    llik <- mean(tail(llik.vec,2))
    prev_llik <- mean(tail(llik.vec,4)[1:2])

    if (config$verbose == T){
      print(paste("Iteration: ",count," with Tolerance of: ",abs(llik-prev_llik)," and Log-Lik of: ",tail(llik.vec,1),sep=""))
    }

    # Check convergence
    if (count>config$min_iter){
      if ((count>config$max_iter) | abs(llik - prev_llik) < config$tol ){
        break
      }
    }

    # Step sizes
    b.a <- max(1e-3,config$eta * (1/(count)))
    a.b <- c(1-b.a)

    count = count + 1


    debug=T
    if (debug==T & count > 100+1){

      predictions.x <- (predict$predictions.x)
      predictions.y <- (predict$predictions.y)

      par(mfcol=c(3,2))
      plot(predictions.y[-train_id,id_col]~data[-train_id,id_col],main=summary(lm(predictions.y[-train_id,id_col]~data[-train_id,id_col]))$r.squared)
      abline(lm(predictions.y[-train_id,id_col]~data[-train_id,id_col]),col="red")

      plot(predictions.x[-train_id,id_col]~data[-train_id,id_col],main=summary(lm(predictions.x[-train_id,id_col]~data[-train_id,id_col]))$r.squared)
      abline(lm(predictions.x[-train_id,id_col]~data[-train_id,id_col]),col="red")

      plot(glmnet_predict~data[-train_id,id_col],main=summary(lm(glmnet_predict~data[-train_id,id_col]))$r.squared)
      abline(lm(glmnet_predict~data[-train_id,id_col]),col="red")

      mean(abs(glmnet_predict-data[-train_id,id_col]))
      mean(abs(predictions.x[-train_id,id_col]-data[-train_id,id_col]))

      plot(predictions.y[,id_col]~data[,id_col],main=summary(lm(predictions.y[,id_col]~data[,id_col]))$r.squared)
      abline(lm(predictions.y[,id_col]~data[,id_col]),col="red")

      plot(predictions.x[,id_col]~data[,id_col],main=summary(lm(predictions.x[,id_col]~data[,id_col]))$r.squared)
      abline(lm(predictions.x[,id_col]~data[,id_col]),col="red")

    }
  }



  return(list(

    main.parameters = list(
      alpha.L.K = alpha.L.K.star.alpha.L.K.final,
      v.beta = v.beta.star.beta.final,
      alpha.L.J = alpha.L.J.star.alpha.L.J.final,
      u.beta = u.beta.star.beta.final,
      y.gamma = y.gamma.final,
      y.delta = y.delta.final,
      x.gamma = x.gamma.final,
      x.delta = x.delta.final
    ),

    code = list(
      intercept.code = intercept.code,
      intercept.encode = intercept.encode,
      X_code = X_code,
      Y_code = Y_code,
      X_encode = X_encode,
      Y_encode = Y_encode
    ),

    meta.parameters = config,

    convergence.parameters = list(
      iterations = count,
      score.vec = score.vec,
      llik.vec = llik.vec
    )

  ))


}

