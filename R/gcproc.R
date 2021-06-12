gcproc <- function(x,
                   y,
                   k_dim = 30,
                   j_dim = 30,
                   eta=1e-2,
                   max_iter=1500,
                   min_iter = 15,
                   tol=1e-3,
                   log=F,
                   center=F,
                   scale.z=F,
                   batches=2,
                   cores=2,
                   verbose=F,
                   init="eigen",
                   anchors = NULL){

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

  if (verbose){
    print(paste("Using gcproc to dimensionally reduce both samples and features with following dimensions:   Sample dimension (k_dim): ",k_dim, "   Feature dimension (j_dim): ", j_dim,sep=""))
  }

  n.x <- dim(x)[1]
  p.x <- dim(x)[2]
  n.y <- dim(y)[1]
  p.y <- dim(y)[2]

  if (prepare_data == T){
    transformed.data <- prepare_data(x=x,
                                     y=y,
                                     log=log,
                                     center=center,
                                     scale=scale.z)


    if (verbose){
      print("Prepared transformations on data are (True/False):")
      print(paste("Natural Log (with pseudocount): ", log, sep = ""))
      print(paste("Center: ", center, sep = ""))
      print(paste("Standardise by standard deviation: ", scale.z, sep = ""))
    }
    x <- transformed.data$x
    y <- transformed.data$y
  }


  if (initialise==T){
    if (verbose){
      print("Initialising data")
    }

    if (init == "eigen"){

      u.beta.eigen <- RSpectra::eigs_sym(t(x)%*%x,j_dim)
      v.beta.eigen <- RSpectra::eigs_sym(t(y)%*%y,j_dim)

      u.beta.star.beta <- u.beta.eigen$vectors
      v.beta.star.beta <- v.beta.eigen$vectors

    }
    if (init == "svd"){

      u.beta.svd <- irlba::irlba(t(x)%*%x,j_dim)
      v.beta.svd <- irlba::irlba(t(y)%*%y,j_dim)

      u.beta.star.beta <- u.beta.svd$v
      v.beta.star.beta <- v.beta.svd$v

    }
    if (init == "random"){

      u.beta.star.beta <- matrix(rnorm(j_dim*dim(x)[2]),ncol = j_dim,nrow=dim(x)[2])
      v.beta.star.beta <- matrix(rnorm(j_dim*dim(y)[2]),ncol = j_dim,nrow=dim(y)[2])

    }

    a0.beta = 10e-2
    b0.beta = 10e-4
    c0.beta = 10e-2
    d0.beta = 10e-4

    # #Initialise u.beta
    v.V.star.inv.beta = t(y)%*%(y)

    #Initialise u.beta
    u.V.star.inv.beta = t(x)%*%(x)

    if (init == "eigen"){

      alpha.L.J.eigen <- RSpectra::eigs_sym(x%*%t(x),k_dim)
      alpha.L.K.eigen <- RSpectra::eigs_sym(y%*%t(y),k_dim)

      alpha.L.J.star.alpha.L.J = t(alpha.L.J.eigen$vectors)
      alpha.L.K.star.alpha.L.K = t(alpha.L.K.eigen$vectors)

    }
    if (init == "svd"){

      alpha.L.J.svd <- irlba::irlba(x%*%t(x),k_dim)
      alpha.L.K.svd <- irlba::irlba(y%*%t(y),k_dim)

      alpha.L.J.star.alpha.L.J = t(alpha.L.J.svd$u)
      alpha.L.K.star.alpha.L.K = t(alpha.L.K.svd$u)

    }
    if (init == "random"){

      alpha.L.J.star.alpha.L.J <- matrix(rnorm(k_dim*dim(x)[1]),nrow = k_dim,ncol=dim(x)[1])
      alpha.L.K.star.alpha.L.K <- matrix(rnorm(k_dim*dim(y)[1]),nrow = k_dim,ncol=dim(y)[1])

    }

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


    count = 0
    llik.vec <- c()
    prev.internal.score.vec <- rep(Inf,10)
    score.vec <- 0

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


  if (verbose){
    print(paste("Beginning gcproc learning with:   Cores: ",cores, "   Batches: ", batches, "   Learning Rate (eta): ", eta, "   Tolerance Threshold: ", tol, "   Minimum Number of iterations: ",min_iter, "   Maximum number of iterations: ", max_iter, "   Verbose: ", verbose,sep=""))
  }

  while (T){

    if (variational_gradient_descent_updates == T){

      x.g.sample <- if(dim(X.x)[1]>batches){chunk(sample(c(1:dim(X.x)[1])),batches)}else{chunk(sample(c(1:dim(X.x)[1])),1)}
      y.g.sample <- if(dim(Y.y)[1]>batches){chunk(sample(c(1:dim(Y.y)[1])),batches)}else{chunk(sample(c(1:dim(Y.y)[1])),1)}

      x.v.sample <- if(dim(X.x)[2]>batches){chunk(sample(c(1:dim(X.x)[2])),batches)}else{chunk(sample(c(1:dim(X.x)[2])),1)}
      y.v.sample <- if(dim(Y.y)[2]>batches){chunk(sample(c(1:dim(Y.y)[2])),batches)}else{chunk(sample(c(1:dim(Y.y)[2])),1)}


      to_return <- pbmcapply::pbmclapply(c(1:batches),function(i){
        set.seed(i*count)

        x.g.ids <- x.g.sample[[min(length(x.g.sample),i)]]
        y.g.ids <- y.g.sample[[min(length(y.g.sample),i)]]

        x.v.ids <- x.v.sample[[min(length(x.v.sample),i)]]
        y.v.ids <- y.v.sample[[min(length(y.v.sample),i)]]

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




        a.star.beta = a0.beta + dim(y)[2]/2
        b.star.beta = b0.beta + (1/2)*((t(alpha.L.K.star.alpha.L.K%*%y%*%v.beta.star.beta)%*%(alpha.L.K.star.alpha.L.K%*%y%*%v.beta.star.beta)) - t(u.beta.star.beta)%*%t(alpha.L.J.star.alpha.L.J%*%x)%*%(alpha.L.J.star.alpha.L.J%*%x)%*%u.beta.star.beta)
        c.star.beta = c0.beta + 1/2

        v.E.beta.D.beta <- a.star.beta*(((v.beta.star.beta)%*%(MASS::ginv(b.star.beta))%*%t(v.beta.star.beta)))+MASS::ginv(v.V.star.inv.beta)
        v.E.beta.D.beta <- p.y/min(dim(y)[2],length(y.v.ids))*v.E.beta.D.beta

        d.star.beta = d0.beta + (1/2)*(v.E.beta.D.beta)

        v.E.diag.alpha.beta = c.star.beta*MASS::ginv(d.star.beta)
        v.E.diag.alpha.beta <- p.y/min(dim(y)[2],length(y.v.ids))*v.E.diag.alpha.beta



        a.star.beta = a0.beta + dim(x)[2]/2
        b.star.beta = b0.beta + (1/2)*((t(alpha.L.J.star.alpha.L.J%*%x%*%u.beta.star.beta)%*%(alpha.L.J.star.alpha.L.J%*%x%*%u.beta.star.beta)) - t(v.beta.star.beta)%*%t(alpha.L.K.star.alpha.L.K%*%y)%*%(alpha.L.K.star.alpha.L.K%*%y)%*%v.beta.star.beta)
        c.star.beta = c0.beta + 1/2

        u.E.beta.D.beta <- a.star.beta*(((u.beta.star.beta)%*%MASS::ginv(b.star.beta)%*%t(u.beta.star.beta)))+(MASS::ginv(u.V.star.inv.beta))
        u.E.beta.D.beta <- p.x/min(dim(x)[2],length(x.v.ids))*u.E.beta.D.beta

        d.star.beta = d0.beta + (1/2)*(u.E.beta.D.beta)

        u.E.diag.alpha.beta = c.star.beta*MASS::ginv(d.star.beta)
        u.E.diag.alpha.beta <- p.x/min(dim(x)[2],length(x.v.ids))*u.E.diag.alpha.beta




        a.star.alpha.L.K = a0.alpha.L.K + dim(y)[1]/2
        b.star.alpha.L.K = b0.alpha.L.K + (1/2)*(((alpha.L.J.star.alpha.L.J%*%x%*%u.beta.star.beta)%*%t(alpha.L.J.star.alpha.L.J%*%x%*%u.beta.star.beta)) - (alpha.L.K.star.alpha.L.K)%*%V.star.inv.alpha.L.K%*%t(alpha.L.K.star.alpha.L.K))
        c.star.alpha.L.K = c0.alpha.L.K + 1/2

        E.alpha.L.K.D.alpha.L.K <- a.star.alpha.L.K*((t(alpha.L.K.star.alpha.L.K)%*%(MASS::ginv(b.star.alpha.L.K))%*%(alpha.L.K.star.alpha.L.K)))+(MASS::ginv(V.star.inv.alpha.L.K))
        E.alpha.L.K.D.alpha.L.K <- n.y/min(dim(y)[1],length(y.g.ids))*E.alpha.L.K.D.alpha.L.K

        d.star.alpha.L.K = d0.alpha.L.K + (1/2)*(E.alpha.L.K.D.alpha.L.K)

        E.diag.alpha.alpha.L.K = c.star.alpha.L.K*MASS::ginv(d.star.alpha.L.K)
        E.diag.alpha.alpha.L.K <- n.y/min(dim(y)[1],length(y.g.ids))*E.diag.alpha.alpha.L.K


        a.star.alpha.L.J = a0.alpha.L.J + dim(x)[1]/2
        b.star.alpha.L.J = b0.alpha.L.J + (1/2)*(((alpha.L.K.star.alpha.L.K%*%y%*%v.beta.star.beta)%*%t(alpha.L.K.star.alpha.L.K%*%y%*%v.beta.star.beta)) - (alpha.L.J.star.alpha.L.J)%*%V.star.inv.alpha.L.J%*%t(alpha.L.J.star.alpha.L.J))
        c.star.alpha.L.J = c0.alpha.L.J + 1/2

        E.alpha.L.J.D.alpha.L.J <- a.star.alpha.L.J*((t(alpha.L.J.star.alpha.L.J)%*%(MASS::ginv(b.star.alpha.L.J))%*%(alpha.L.J.star.alpha.L.J)))+(MASS::ginv(V.star.inv.alpha.L.J))
        E.alpha.L.J.D.alpha.L.J <- n.x/min(dim(x)[1],length(x.g.ids))*E.alpha.L.J.D.alpha.L.J

        d.star.alpha.L.J = d0.alpha.L.J + (1/2)*(E.alpha.L.J.D.alpha.L.J)

        E.diag.alpha.alpha.L.J = c.star.alpha.L.J*MASS::ginv(d.star.alpha.L.J)
        E.diag.alpha.alpha.L.J <- n.x/min(dim(x)[1],length(x.g.ids))*E.diag.alpha.alpha.L.J




        u.V.star.inv.beta = u.E.diag.alpha.beta + t(alpha.L.J.star.alpha.L.J%*%x)%*%(alpha.L.J.star.alpha.L.J%*%x)
        u.beta.star.beta = if (is.null(anchor_x.feature)){MASS::ginv(u.V.star.inv.beta)%*%t((alpha.L.J.star.alpha.L.J%*%x))%*%(alpha.L.K.star.alpha.L.K%*%y%*%v.beta.star.beta)}else{anchor_x.feature[x.v.ids,]}

        v.V.star.inv.beta = v.E.diag.alpha.beta + t(alpha.L.K.star.alpha.L.K%*%y)%*%(alpha.L.K.star.alpha.L.K%*%y)
        v.beta.star.beta = if (is.null(anchor_y.feature)){MASS::ginv(v.V.star.inv.beta)%*%t((alpha.L.K.star.alpha.L.K%*%y))%*%(alpha.L.J.star.alpha.L.J%*%x%*%u.beta.star.beta)}else{anchor_y.feature[y.v.ids,]}


        V.star.inv.alpha.L.J = E.diag.alpha.alpha.L.J + (x%*%u.beta.star.beta)%*%t(x%*%u.beta.star.beta)
        alpha.L.J.star.alpha.L.J = if (is.null(anchor_x.sample)){(alpha.L.K.star.alpha.L.K%*%y%*%v.beta.star.beta)%*%t(x%*%u.beta.star.beta)%*%MASS::ginv(V.star.inv.alpha.L.J)}else{anchor_x.sample[,x.g.ids]}

        V.star.inv.alpha.L.K = E.diag.alpha.alpha.L.K + (y%*%v.beta.star.beta)%*%t(y%*%v.beta.star.beta)
        alpha.L.K.star.alpha.L.K = if (is.null(anchor_y.sample)){(alpha.L.J.star.alpha.L.J%*%x%*%u.beta.star.beta)%*%t(y%*%v.beta.star.beta)%*%MASS::ginv(V.star.inv.alpha.L.K)}else{anchor_y.sample[,y.g.ids]}


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
      },mc.cores = cores)


    }


    if (update_batched_parameters==T){
      for (i in 1:batches){
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

        a.b <- c(1-eta)
        b.a <- c(eta)

        alpha.L.J.star.alpha.L.J.final[,x.g.ids] <- a.b*alpha.L.J.star.alpha.L.J.final[,x.g.ids] + b.a*alpha.L.J.star.alpha.L.J
        V.star.inv.alpha.L.J.final[x.g.ids,x.g.ids] <- a.b*V.star.inv.alpha.L.J.final[x.g.ids,x.g.ids] + b.a*V.star.inv.alpha.L.J

        alpha.L.K.star.alpha.L.K.final[,y.g.ids] <- a.b*alpha.L.K.star.alpha.L.K.final[,y.g.ids] + b.a*alpha.L.K.star.alpha.L.K
        V.star.inv.alpha.L.K.final[y.g.ids,y.g.ids] <- a.b*V.star.inv.alpha.L.K.final[y.g.ids,y.g.ids] + b.a*V.star.inv.alpha.L.K

        u.beta.star.beta.final[x.v.ids,] <- a.b*u.beta.star.beta.final[x.v.ids,] + b.a*u.beta.star.beta
        v.beta.star.beta.final[y.v.ids,] <- a.b*v.beta.star.beta.final[y.v.ids,] + b.a*v.beta.star.beta

        u.V.star.inv.beta.final[x.v.ids,x.v.ids] <- a.b*u.V.star.inv.beta.final[x.v.ids,x.v.ids] + b.a*u.V.star.inv.beta
        v.V.star.inv.beta.final[y.v.ids,y.v.ids] <- a.b*v.V.star.inv.beta.final[y.v.ids,y.v.ids] + b.a*v.V.star.inv.beta
      }

    }

    if (check_anchors==T){

      alpha.L.J.star.alpha.L.J.final <- if (is.null(anchor_x.sample)){alpha.L.J.star.alpha.L.J.final}else{anchor_x.sample}
      alpha.L.K.star.alpha.L.K.final <- if (is.null(anchor_y.sample)){alpha.L.K.star.alpha.L.K.final}else{anchor_y.sample}

      u.beta.star.beta.final <- if (is.null(anchor_x.feature)){u.beta.star.beta.final}else{anchor_x.feature}
      v.beta.star.beta.final <- if (is.null(anchor_y.feature)){v.beta.star.beta.final}else{anchor_y.feature}

    }

    matrix.residuals <- alpha.L.K.star.alpha.L.K.final%*%Y.y%*%v.beta.star.beta.final - alpha.L.J.star.alpha.L.J.final%*%X.x%*%u.beta.star.beta.final

    llik.vec <- c(llik.vec, mean(mclust::dmvnorm(matrix.residuals,sigma = diag(diag(t(matrix.residuals)%*%matrix.residuals/dim(matrix.residuals)[2])),log = T)))
    score.vec <- c(score.vec, (mean(abs(c(matrix.residuals)))))

    # Check convergence
    if (count>min_iter){
      if ((count>max_iter) | tail(score.vec,1)<tol){
        break
      }
    }

    count = count + 1

    if (verbose == T){
      print(paste("Iteration: ",count," with MSE of: ",tail(score.vec,1),sep=""))
    }

  }

  return(list(

    transformed.data = transformed.data,

    residuals = residuals,

    main.parameters = list(
      alpha.L.J = alpha.L.J.star.alpha.L.J.final,
      alpha.L.K = alpha.L.K.star.alpha.L.K.final,
      u.beta = u.beta.star.beta.final,
      v.beta = v.beta.star.beta.final
    ),

    meta.parameters = list(
      k_dim = k_dim,
      j_dim = j_dim,
      eta= eta,
      max_iter = max_iter,
      min_iter = min_iter,
      tol = tol,
      log = log,
      center = center,
      scale.z = scale.z,
      batches = batches,
      cores = cores,
      verbose = verbose,
      init = init
    ),

    convergence.parameters = list(
      iterations = count,
      score.vec = score.vec,
      llik.vec = llik.vec
    )

  ))


}





