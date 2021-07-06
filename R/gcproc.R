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
                                 batches=64,
                                 cores=8,
                                 verbose=T,
                                 init="svd"),
                   seed = 1,
                   anchors = NULL,
                   pivots = NULL
                   ){

  prepare_data = TRUE
  initialise = TRUE
  variational_gradient_descent_updates = TRUE
  update_batched_parameters = TRUE

  check_anchors = FALSE

  anchor_y.sample = NULL
  anchor_y.feature = NULL
  anchor_x.sample = NULL
  anchor_x.feature = NULL

  if (!is.null(anchors)){
    check_anchors = TRUE

    anchor_y.sample = anchors$anchor_y.sample
    anchor_y.feature = anchors$anchor_y.feature
    anchor_x.sample = anchors$anchor_x.sample
    anchor_x.feature = anchors$anchor_x.feature

  }

  if (config$verbose){
    print(paste("Using gcproc to dimensionally reduce both samples and features with following dimensions:   Sample dimension (config$k_dim): ",config$k_dim, "   Feature dimension (config$j_dim): ", config$j_dim,sep=""))
  }

  n.x <- dim(x)[1]
  p.x <- dim(x)[2]
  n.y <- dim(y)[1]
  p.y <- dim(y)[2]

  if (prepare_data == T){
    transformed.data <- prepare_data(x=x,
                                     y=y,
                                     log=config$log,
                                     center=config$center,
                                     scale=config$scale.z)


    if (config$verbose){
      print("Prepared transformations on data are (True/False):")
      print(paste("Natural Log (with pseudocount): ", config$log, sep = ""))
      print(paste("Center: ", config$center, sep = ""))
      print(paste("Standardise by standard deviation: ", config$scale.z, sep = ""))
    }
    x <- transformed.data$x
    y <- transformed.data$y
  }


  if (initialise==T){
    if (config$verbose){
      print("Initialising data")
    }

    if (is.null(pivots)){
      initial.param <-initialise.gcproc(x=x,y=y,k_dim=config$k_dim,j_dim=config$j_dim,init=config$init,verbose=config$verbose)
    } else {
      initial.param <- pivots
    }

    alpha.L.K.star.alpha.L.K <- if (is.null(anchor_y.sample)){initial.param$anchor_y.sample}else{anchor_y.sample}
    alpha.L.J.star.alpha.L.J <- if (is.null(anchor_x.sample)){initial.param$anchor_x.sample}else{anchor_x.sample}
    v.beta.star.beta <- if (is.null(anchor_y.feature)){initial.param$anchor_y.feature}else{anchor_y.feature}
    u.beta.star.beta <- if (is.null(anchor_x.feature)){initial.param$anchor_x.feature}else{anchor_x.feature}


    a0.beta = 10e-2
    b0.beta = 10e-4
    c0.beta = 10e-2
    d0.beta = 10e-4

    # #Initialise u.beta
    v.V.star.inv.beta = t(alpha.L.K.star.alpha.L.K%*%y)%*%(alpha.L.K.star.alpha.L.K%*%y)

    #Initialise u.beta
    u.V.star.inv.beta = t(alpha.L.J.star.alpha.L.J%*%x)%*%(alpha.L.J.star.alpha.L.J%*%x)


    V.star.inv.alpha.L.K = (y%*%v.beta.star.beta)%*%t(y%*%v.beta.star.beta)
    V.star.inv.alpha.L.J = (x%*%u.beta.star.beta)%*%t(x%*%u.beta.star.beta)


    # Initialise alpha.L.J
    a0.alpha.L.J = 10e-2
    b0.alpha.L.J = 10e-4
    c0.alpha.L.J = 10e-2
    d0.alpha.L.J = 10e-4

    # Initialise alpha.L.K
    a0.alpha.L.K = 10e-2
    b0.alpha.L.K = 10e-4
    c0.alpha.L.K = 10e-2
    d0.alpha.L.K = 10e-4

    count = 1
    llik.vec <- c()
    score.vec <- 10e8

    X.x <- x
    Y.y <- y

    alpha.L.J.star.alpha.L.J.final <- if (is.null(anchor_x.sample)){alpha.L.J.star.alpha.L.J}else{anchor_x.sample}
    alpha.L.K.star.alpha.L.K.final <- if (is.null(anchor_y.sample)){alpha.L.K.star.alpha.L.K}else{anchor_y.sample}

    u.beta.star.beta.final <- if (is.null(anchor_x.feature)){u.beta.star.beta}else{anchor_x.feature}
    v.beta.star.beta.final <- if (is.null(anchor_y.feature)){v.beta.star.beta}else{anchor_y.feature}

    V.star.inv.alpha.L.J.final <- V.star.inv.alpha.L.J
    V.star.inv.alpha.L.K.final <- V.star.inv.alpha.L.K

    u.V.star.inv.beta.final <- u.V.star.inv.beta
    v.V.star.inv.beta.final <- v.V.star.inv.beta

  }


  if (config$verbose){
    print(paste("Beginning gcproc learning with:   Cores: ",config$cores, "   Batches: ", config$batches, "   Learning Rate (config$eta): ", config$eta, "   Tolerance Threshold: ", config$tol, "   Minimum Number of iterations: ",config$min_iter, "   Maximum number of iterations: ", config$max_iter, "   Verbose: ", config$verbose,sep=""))
  }

  MSE <- Inf
  bar_count <- 0
  bar_c = F

  while (T){

    if (variational_gradient_descent_updates == T){
      set.seed(seed+count)

      x.g.sample <- lapply(c(1:config$batches),function(X){sample(c(1:dim(X.x)[1]),size = dim(X.x)[1]/(config$batches))})
      y.g.sample <- lapply(c(1:config$batches),function(X){sample(c(1:dim(Y.y)[1]),size = dim(Y.y)[1]/(config$batches))})
      x.v.sample <- lapply(c(1:config$batches),function(X){sample(c(1:dim(X.x)[2]),size = dim(X.x)[2]/(config$batches))})
      y.v.sample <- lapply(c(1:config$batches),function(X){sample(c(1:dim(Y.y)[2]),size = dim(Y.y)[2]/(config$batches))})

      to_return <- parallel::mclapply(c(1:config$batches),function(i){


        x.g.ids <- x.g.sample[[i]]
        y.g.ids <- y.g.sample[[i]]

        x.v.ids <- x.v.sample[[i]]
        y.v.ids <- y.v.sample[[i]]


        x <- X.x[x.g.ids,x.v.ids]
        y <- Y.y[y.g.ids,y.v.ids]

        alpha.L.J.star.alpha.L.J <- alpha.L.J.star.alpha.L.J.final[,x.g.ids]
        alpha.L.K.star.alpha.L.K <- alpha.L.K.star.alpha.L.K.final[,y.g.ids]

        u.beta.star.beta <- u.beta.star.beta.final[x.v.ids,]
        v.beta.star.beta <- v.beta.star.beta.final[y.v.ids,]

        V.star.inv.alpha.L.J <- V.star.inv.alpha.L.J.final[x.g.ids,x.g.ids]
        V.star.inv.alpha.L.K <- V.star.inv.alpha.L.K.final[y.g.ids,y.g.ids]

        v.V.star.inv.beta <- v.V.star.inv.beta.final[y.v.ids,y.v.ids]
        u.V.star.inv.beta <- u.V.star.inv.beta.final[x.v.ids,x.v.ids]

        a.star.alpha.L.J = a0.alpha.L.J + min(dim(x)[1],length(x.g.ids))/2
        b.star.alpha.L.J = b0.alpha.L.J + (1/2)*(((alpha.L.K.star.alpha.L.K%*%y%*%v.beta.star.beta)%*%t(alpha.L.K.star.alpha.L.K%*%y%*%v.beta.star.beta)) - (alpha.L.J.star.alpha.L.J)%*%V.star.inv.alpha.L.J%*%t(alpha.L.J.star.alpha.L.J))
        c.star.alpha.L.J = c0.alpha.L.J + 1/2

        E.alpha.L.J.D.alpha.L.J <- a.star.alpha.L.J*((t(alpha.L.J.star.alpha.L.J)%*%(MASS::ginv(b.star.alpha.L.J))%*%(alpha.L.J.star.alpha.L.J)))+(MASS::ginv(V.star.inv.alpha.L.J))
        E.alpha.L.J.D.alpha.L.J <- n.x/min(dim(x)[1],length(x.g.ids))*E.alpha.L.J.D.alpha.L.J

        d.star.alpha.L.J = d0.alpha.L.J + (1/2)*(E.alpha.L.J.D.alpha.L.J)

        E.diag.alpha.alpha.L.J = c.star.alpha.L.J*MASS::ginv(d.star.alpha.L.J)
        E.diag.alpha.alpha.L.J <- n.x/min(dim(x)[1],length(x.g.ids))*E.diag.alpha.alpha.L.J

        a.star.alpha.L.K = a0.alpha.L.K + min(dim(y)[1],length(y.g.ids))/2
        b.star.alpha.L.K = b0.alpha.L.K + (1/2)*(((alpha.L.J.star.alpha.L.J%*%x%*%u.beta.star.beta)%*%t(alpha.L.J.star.alpha.L.J%*%x%*%u.beta.star.beta)) - (alpha.L.K.star.alpha.L.K)%*%V.star.inv.alpha.L.K%*%t(alpha.L.K.star.alpha.L.K))
        c.star.alpha.L.K = c0.alpha.L.K + 1/2

        E.alpha.L.K.D.alpha.L.K <- a.star.alpha.L.K*((t(alpha.L.K.star.alpha.L.K)%*%(MASS::ginv(b.star.alpha.L.K))%*%(alpha.L.K.star.alpha.L.K)))+(MASS::ginv(V.star.inv.alpha.L.K))
        E.alpha.L.K.D.alpha.L.K <- n.y/min(dim(y)[1],length(y.g.ids))*E.alpha.L.K.D.alpha.L.K

        d.star.alpha.L.K = d0.alpha.L.K + (1/2)*(E.alpha.L.K.D.alpha.L.K)

        E.diag.alpha.alpha.L.K = c.star.alpha.L.K*MASS::ginv(d.star.alpha.L.K)
        E.diag.alpha.alpha.L.K <- n.y/min(dim(y)[1],length(y.g.ids))*E.diag.alpha.alpha.L.K

        a.star.beta = a0.beta + min(dim(x)[2],length(x.v.ids))/2
        b.star.beta = b0.beta + (1/2)*((t(alpha.L.J.star.alpha.L.J%*%x%*%u.beta.star.beta)%*%(alpha.L.J.star.alpha.L.J%*%x%*%u.beta.star.beta)) - t(v.beta.star.beta)%*%t(alpha.L.K.star.alpha.L.K%*%y)%*%(alpha.L.K.star.alpha.L.K%*%y)%*%v.beta.star.beta)
        c.star.beta = c0.beta + 1/2

        u.E.beta.D.beta <- a.star.beta*(((u.beta.star.beta)%*%MASS::ginv(b.star.beta)%*%t(u.beta.star.beta)))+(MASS::ginv(u.V.star.inv.beta))
        u.E.beta.D.beta <- p.x/min(dim(x)[2],length(x.v.ids))*u.E.beta.D.beta

        d.star.beta = d0.beta + (1/2)*(u.E.beta.D.beta)

        u.E.diag.alpha.beta = c.star.beta*MASS::ginv(d.star.beta)
        u.E.diag.alpha.beta <- p.x/min(dim(x)[2],length(x.v.ids))*u.E.diag.alpha.beta

        a.star.beta = a0.beta + min(dim(y)[2],length(y.v.ids))/2
        b.star.beta = b0.beta + (1/2)*((t(alpha.L.K.star.alpha.L.K%*%y%*%v.beta.star.beta)%*%(alpha.L.K.star.alpha.L.K%*%y%*%v.beta.star.beta)) - t(u.beta.star.beta)%*%t(alpha.L.J.star.alpha.L.J%*%x)%*%(alpha.L.J.star.alpha.L.J%*%x)%*%u.beta.star.beta)
        c.star.beta = c0.beta + 1/2

        v.E.beta.D.beta <- a.star.beta*(((v.beta.star.beta)%*%(MASS::ginv(b.star.beta))%*%t(v.beta.star.beta)))+MASS::ginv(v.V.star.inv.beta)
        v.E.beta.D.beta <- p.y/min(dim(y)[2],length(y.v.ids))*v.E.beta.D.beta

        d.star.beta = d0.beta + (1/2)*(v.E.beta.D.beta)

        v.E.diag.alpha.beta = c.star.beta*MASS::ginv(d.star.beta)
        v.E.diag.alpha.beta <- p.y/min(dim(y)[2],length(y.v.ids))*v.E.diag.alpha.beta


        V.star.inv.alpha.L.J = E.diag.alpha.alpha.L.J + (x%*%u.beta.star.beta)%*%t(x%*%u.beta.star.beta)
        alpha.L.J.star.alpha.L.J = if (is.null(anchor_x.sample)){(alpha.L.K.star.alpha.L.K%*%y%*%v.beta.star.beta)%*%t(x%*%u.beta.star.beta)%*%MASS::ginv(V.star.inv.alpha.L.J)}else{anchor_x.sample[,x.g.ids]}

        V.star.inv.alpha.L.K = E.diag.alpha.alpha.L.K + (y%*%v.beta.star.beta)%*%t(y%*%v.beta.star.beta)
        alpha.L.K.star.alpha.L.K = if (is.null(anchor_y.sample)){( alpha.L.J.star.alpha.L.J%*%x%*%u.beta.star.beta)%*%t(y%*%v.beta.star.beta)%*%MASS::ginv(V.star.inv.alpha.L.K)}else{anchor_y.sample[,y.g.ids]}

        u.V.star.inv.beta = u.E.diag.alpha.beta + t(alpha.L.J.star.alpha.L.J%*%x)%*%(alpha.L.J.star.alpha.L.J%*%x)
        u.beta.star.beta = if (is.null(anchor_x.feature)){MASS::ginv(u.V.star.inv.beta)%*%t((alpha.L.J.star.alpha.L.J%*%x))%*%(alpha.L.K.star.alpha.L.K%*%y%*%v.beta.star.beta)}else{anchor_x.feature[x.v.ids,]}

        v.V.star.inv.beta = v.E.diag.alpha.beta + t(alpha.L.K.star.alpha.L.K%*%y)%*%(alpha.L.K.star.alpha.L.K%*%y)
        v.beta.star.beta = if (is.null(anchor_y.feature)){MASS::ginv(v.V.star.inv.beta)%*%t((alpha.L.K.star.alpha.L.K%*%y))%*%( alpha.L.J.star.alpha.L.J%*%x%*%u.beta.star.beta)}else{anchor_y.feature[y.v.ids,]}



        return(list(
          alpha.L.J.star.alpha.L.J = alpha.L.J.star.alpha.L.J,
          alpha.L.K.star.alpha.L.K = alpha.L.K.star.alpha.L.K,
          V.star.inv.alpha.L.J = V.star.inv.alpha.L.J,
          V.star.inv.alpha.L.K = V.star.inv.alpha.L.K,
          v.beta.star.beta = v.beta.star.beta,
          u.beta.star.beta = u.beta.star.beta,
          u.V.star.inv.beta = u.V.star.inv.beta,
          v.V.star.inv.beta = v.V.star.inv.beta,
          y.g.ids = y.g.ids,
          x.g.ids = x.g.ids,
          x.v.ids = x.v.ids,
          y.v.ids = y.v.ids
        ))
      },mc.silent = config$verbose,mc.cores = config$cores)


    }
    # print("up2")


    if (update_batched_parameters==T){


      for (i in 1:config$batches){
        # print("up3")

        alpha.L.J.star.alpha.L.J <- to_return[[i]]$alpha.L.J.star.alpha.L.J
        alpha.L.K.star.alpha.L.K <- to_return[[i]]$alpha.L.K.star.alpha.L.K
        V.star.inv.alpha.L.J <- to_return[[i]]$V.star.inv.alpha.L.J
        V.star.inv.alpha.L.K <- to_return[[i]]$V.star.inv.alpha.L.K

        v.beta.star.beta <- to_return[[i]]$v.beta.star.beta
        u.beta.star.beta <- to_return[[i]]$u.beta.star.beta
        u.V.star.inv.beta <- to_return[[i]]$u.V.star.inv.beta
        v.V.star.inv.beta <- to_return[[i]]$v.V.star.inv.beta

        x.g.ids <- to_return[[i]]$x.g.ids
        y.g.ids <- to_return[[i]]$y.g.ids

        x.v.ids <- to_return[[i]]$x.v.ids
        y.v.ids <- to_return[[i]]$y.v.ids

        b.a <- config$eta
        a.b <- c(1-b.a)

        alpha.L.J.star.alpha.L.J.final[,x.g.ids] <- a.b*alpha.L.J.star.alpha.L.J.final[,x.g.ids] + b.a*alpha.L.J.star.alpha.L.J
        V.star.inv.alpha.L.J.final[x.g.ids,x.g.ids] <- a.b*V.star.inv.alpha.L.J.final[x.g.ids,x.g.ids] + b.a*V.star.inv.alpha.L.J

        alpha.L.K.star.alpha.L.K.final[,y.g.ids] <- a.b*alpha.L.K.star.alpha.L.K.final[,y.g.ids] + b.a*alpha.L.K.star.alpha.L.K
        V.star.inv.alpha.L.K.final[y.g.ids,y.g.ids] <- a.b*V.star.inv.alpha.L.K.final[y.g.ids,y.g.ids] + b.a*V.star.inv.alpha.L.K

        u.beta.star.beta.final[x.v.ids,] <- a.b*u.beta.star.beta.final[x.v.ids,] + b.a*u.beta.star.beta
        u.V.star.inv.beta.final[x.v.ids,x.v.ids] <- a.b*u.V.star.inv.beta.final[x.v.ids,x.v.ids] + b.a*u.V.star.inv.beta

        v.beta.star.beta.final[y.v.ids,] <- a.b*v.beta.star.beta.final[y.v.ids,] + b.a*v.beta.star.beta
        v.V.star.inv.beta.final[y.v.ids,y.v.ids] <- a.b*v.V.star.inv.beta.final[y.v.ids,y.v.ids] + b.a*v.V.star.inv.beta



      }

    }

    if (check_anchors==T){

      alpha.L.J.star.alpha.L.J.final <- if (is.null(anchor_x.sample)){alpha.L.J.star.alpha.L.J.final}else{anchor_x.sample}
      alpha.L.K.star.alpha.L.K.final <- if (is.null(anchor_y.sample)){alpha.L.K.star.alpha.L.K.final}else{anchor_y.sample}

      u.beta.star.beta.final <- if (is.null(anchor_x.feature)){u.beta.star.beta.final}else{anchor_x.feature}
      v.beta.star.beta.final <- if (is.null(anchor_y.feature)){v.beta.star.beta.final}else{anchor_y.feature}

    }


    matrix.residuals <- alpha.L.K.star.alpha.L.K.final%*%Y.y%*%v.beta.star.beta.final - (alpha.L.J.star.alpha.L.J.final%*%X.x%*%u.beta.star.beta.final)
    intercept <- matrix.residuals

    matrix.residuals <- matrix.residuals


    llik.vec <- c(llik.vec, mean(mclust::dmvnorm(matrix.residuals,sigma = diag(diag(t(matrix.residuals)%*%matrix.residuals/dim(matrix.residuals)[2])),log = T)))
    score.vec <- c(score.vec, (mean(abs(c(matrix.residuals)))))
    MSE <- tail(score.vec,1)

    if (runif(1)<1/(count+1)) {
      bar_count = 0
    } else {
      bar_count <- bar_count + 1
    }
    # print(c(count,bar_count))


    # Check convergence
    if (count>config$min_iter){
      if ((count>config$max_iter) | MSE<config$tol){
        break
      }
    }

    count = count + 1

    if (config$verbose == T){
      print(paste("Iteration: ",count," with MSE of: ",MSE," and Log-Lik of: ",tail(llik.vec,1),sep=""))
    }


  }


  return(list(

    transformed.data = transformed.data,

    residuals = matrix.residuals,

    main.parameters = list(
      alpha.L.J = alpha.L.J.star.alpha.L.J.final,
      alpha.L.K = alpha.L.K.star.alpha.L.K.final,
      u.beta = u.beta.star.beta.final,
      v.beta = v.beta.star.beta.final,
      intercept = intercept
    ),

    meta.parameters = config,

    convergence.parameters = list(
      iterations = count,
      score.vec = score.vec,
      llik.vec = llik.vec
    )

  ))


}





