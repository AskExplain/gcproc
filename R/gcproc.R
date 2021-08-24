#' Generalised Canonical Procrustes
#'
#' A method that uses a likelihood model to align two datasets via an encoding in a lower dimensional space. The coding of the datasets simultaneously can be used to reconstruct the data to predict missing test points. The parameters can be used to reduce either the feature or the sample dimensions into a smaller subspace for further embedding or projection. To run as default, only requires x and y as inputs.
#'
#' @param x Reference dataset of sample by feature matrix (required)
#' @param y Experimental dataset of sample by feature matrix (required)
#' @param reference Choice of shared space representing x or y (required)
#' @param fixed Enforcing feature (i_dim) and, or sample parameters (j_dim) to be identical for x and y (required)
#' @param config Configuration parameters (required, default provided)
#' @param anchors Transferring pre-trained model parameters (not required)
#' @param pivots Initialisation of model parameters (not required)
#' @param predict A binary matrix representing data points to be predicted. Here a matrix representing the dimensions of x and y are provided, where 0 represents the training set, and 1 represents the test set (not required)
#'
#' @return Main parameters contains the learned model parameters. The alpha and beta matrix multiply x and y by, (K)(Y)(v) and (J)(X)(u). By multiplying with the parameter, the dimension of the samples and features can be dimensionally reduced for further visualisation analysis such as embedding or projection.
#'
#' @return Code contains the learned shared encoding space. The encoded space refers to the full dimension reduction of both samples and features after matrix multiplication by parameters K and v for y, as well as, J and u for x. The decode is an estimation of the full matrix dataset, where the code is used and matrix multiplied as t(K)(Y_code)t(v), and t(J)(X_code)t(u) to calculate the decoded estimation.
#'
#' @return Prediction contains the predictions for the test dataset as indicated by a 1 in the binary prediction matrices. In addition to the x and y variables within this list that correspond to the input binary design matrix, the predict.y and predict.x datasets represent the decoded estimation, and predicts any missing test data.
#'
#' @export
gcproc <- function(x,
                   y,
                   reference = "y",
                   fixed = list(i_dim = F, j_dim = F),
                   config = gcproc::extract_config(verbose = F),
                   anchors = gcproc::extract_anchors_framework(verbose = F),
                   pivots = gcproc::extract_pivots_framework(verbose = F),
                   predict = gcproc::extract_prediction_framework(verbose = F)
){


  if (is.null(config$j_dim)){
    eig.val_x <- RSpectra::eigs(crossprod((x),(x)),k = min(ncol(x)-2,30))$values
    eig.val_y <- RSpectra::eigs(crossprod((y),(y)),k = min(ncol(y)-2,30))$values

    main_dim.x <- sum(cumsum(eig.val_x/sum(eig.val_x))<(1-1e-10))
    main_dim.y <- sum(cumsum(eig.val_y/sum(eig.val_y))<(1-1e-10))

    main_dim <- min(main_dim.x,main_dim.y)

    config$j_dim <- min(5,main_dim)
  }


  prepare_data = TRUE
  initialise = TRUE
  variational_gradient_descent_updates = TRUE
  run_covariates = TRUE


  anchor_y.sample = NULL
  anchor_y.feature = NULL
  anchor_x.sample = NULL
  anchor_x.feature = NULL

  if (!is.null(anchors)){
    check_anchors = TRUE

    anchors$anchor_y.sample = anchors$anchor_y.sample
    anchors$anchor_y.feature = anchors$anchor_y.feature
    anchors$anchor_x.sample = anchors$anchor_x.sample
    anchors$anchor_x.feature = anchors$anchor_x.feature

  }

  n.x <- dim(x)[1]
  p.x <- dim(x)[2]
  n.y <- dim(y)[1]
  p.y <- dim(y)[2]

  if (initialise==T){
    if (config$verbose){
      print("Initialising data")
    }

    # Prepare convergence checking parameters
    count = 1
    llik.vec <- c()
    score.vec <- c()

    # Initialise parameters
    if (any(do.call('c',lapply(pivots,function(piv){is.null(piv)})))){
      initial.param <-initialise.gcproc(x=x,y=y,i_dim=config$i_dim,j_dim=config$j_dim,init=config$init,verbose=config$verbose)
    }

    # Check pivoting parameters
    initial.param$pivot_y.sample <- if (is.null(pivots$pivot_y.sample)){initial.param$pivot_y.sample}else{pivots$pivot_y.sample}
    initial.param$pivot_x.sample <- if (is.null(pivots$pivot_x.sample)){initial.param$pivot_x.sample}else{pivots$pivot_x.sample}
    initial.param$pivot_y.feature <- if (is.null(pivots$pivot_y.feature)){initial.param$pivot_y.feature}else{pivots$pivot_y.feature}
    initial.param$pivot_x.feature <- if (is.null(pivots$pivot_x.feature)){initial.param$pivot_x.feature}else{pivots$pivot_x.feature}

    # Check anchoring parameters
    alpha.K <- if (is.null( anchors$anchor_y.sample)){initial.param$pivot_y.sample}else{ anchors$anchor_y.sample}
    alpha.L <- if (is.null( anchors$anchor_x.sample)){initial.param$pivot_x.sample}else{ anchors$anchor_x.sample}
    v.beta <- if (is.null( anchors$anchor_y.feature)){initial.param$pivot_y.feature}else{ anchors$anchor_y.feature}
    u.beta <- if (is.null( anchors$anchor_x.feature)){initial.param$pivot_x.feature}else{ anchors$anchor_x.feature}

    # #Initialise inverse covariance of parameters
    v.V.star.inv.beta.final = t(alpha.K%*%y)%*%(alpha.K%*%y)
    u.V.star.inv.beta.final = t(alpha.L%*%x)%*%(alpha.L%*%x)
    V.star.inv.alpha.K.final = (y%*%v.beta)%*%t(y%*%v.beta)
    V.star.inv.alpha.L.final = (x%*%u.beta)%*%t(x%*%u.beta)

    # Find intercept in endecoded space
    y_encode.final <- alpha.K%*%y%*%v.beta
    x_encode.final <- alpha.L%*%x%*%u.beta

    Y_encode <- (alpha.K%*%(y)%*%(v.beta))
    X_encode <- (alpha.L%*%(x)%*%(u.beta))

    Y_code <- (MASS::ginv((alpha.K)%*%t(alpha.K))%*%Y_encode%*%MASS::ginv(t(v.beta)%*%(v.beta)))

    X_code <- (MASS::ginv((alpha.L)%*%t(alpha.L))%*%(X_encode)%*%MASS::ginv(t(u.beta)%*%(u.beta)))


    Y_decoded <- t(alpha.K)%*%Y_code%*%t(v.beta)
    X_decoded <- t(alpha.L)%*%X_code%*%t(u.beta)


    if (reference == "y"){
      main_code <- Y_code

    }
    if (reference == "x"){
      main_code <- X_code
    }




    main.parameters.copy = main.parameters = list(
      alpha.L = alpha.L,
      alpha.K = alpha.K,
      u.beta = u.beta,
      v.beta = v.beta
    )

    predict =  predict

    code = list(
      main_code = main_code,
      X_code = X_code,
      Y_code = Y_code,
      X_encode = X_encode,
      Y_encode = Y_encode,
      X_decoded = X_decoded,
      Y_decoded = Y_decoded
    )

  }


  if (config$verbose){
    print(paste("Beginning gcproc learning with:    Sample dimension reduction (config$i_dim): ",config$i_dim, "    Feature dimension reduction (config$j_dim): ", config$j_dim,"    Tolerance Threshold: ", config$tol, "   Maximum number of iterations: ", config$max_iter, "   Verbose: ", config$verbose,sep=""))
  }


  while (T){


    if (!is.null(predict$y)){
      Y.y <- scale(y)


      y[which((rowSums(predict$y)>0)==T),] <- do.call('rbind',parallel::mclapply(c(which((rowSums(predict$y)>0)==T)),function(id_row){

        test_id <- as.logical(predict$y[id_row,])
        train_id <- as.logical(1 - predict$y[id_row,])

        sparse.y <- y[id_row,]
        if (any(test_id)){
        covariate_predictors <- rbind(1,(main.parameters$alpha.K[,-id_row]%*%Y.y[-id_row,train_id]))
        test_predictors <- rbind(1,(main.parameters$alpha.K[,-id_row]%*%Y.y[-id_row,test_id]))

        sparse.y[test_id] <- (y[id_row,train_id])%*%t(covariate_predictors)%*%MASS::ginv(covariate_predictors%*%t(covariate_predictors))%*%test_predictors
        }
        return(sparse.y)
      }))

      y <- as.matrix(y)

      y[,which((colSums(predict$y)>0)==T)]  <- do.call('cbind',parallel::mclapply(c(which((colSums(predict$y)>0)==T)),function(id_col){

        test_id <- as.logical(predict$y[,id_col])
        train_id <- as.logical(1 - predict$y[,id_col])
        sparse.y <- y[,id_col]

        if (any(test_id)){

          covariate_predictors <- cbind(1,Y.y[train_id,-id_col]%*%main.parameters$v.beta[-id_col,])
          test_predictors <- cbind(1,Y.y[test_id,-id_col]%*%main.parameters$v.beta[-id_col,])

          sparse.y[test_id] <- ((test_predictors)%*%(MASS::ginv(t(covariate_predictors)%*%(covariate_predictors))%*%t(covariate_predictors)%*%((as.matrix(y[train_id,id_col])))))

        }

        return(sparse.y)
      }))

      y <- as.matrix(y)

    }

    if (!is.null(predict$x)){
      X.x <- scale(x)


      x[which((rowSums(predict$x)>0)==T),]  <- do.call('rbind',parallel::mclapply(c(which((rowSums(predict$x)>0)==T)),function(id_row){

        test_id <- as.logical(predict$x[id_row,])
        train_id <- as.logical(1 - predict$x[id_row,])
        sparse.x <- x[id_row,]

        if (any(test_id)){
          covariate_predictors <- rbind(1,(main.parameters$alpha.L[,-id_row]%*%X.x[-id_row,train_id]))
          test_predictors <- rbind(1,(main.parameters$alpha.L[,-id_row]%*%X.x[-id_row,test_id]))

          sparse.x[test_id] <- (x[id_row,train_id])%*%t(covariate_predictors)%*%MASS::ginv(covariate_predictors%*%t(covariate_predictors))%*%test_predictors
        }
        return(sparse.x)
      }))

      x <- as.matrix(x)

      x[,which((colSums(predict$x)>0)==T)]  <- do.call('cbind',parallel::mclapply(c(which((colSums(predict$x)>0)==T)),function(id_col){

        test_id <- as.logical(predict$x[,id_col])
        train_id <- as.logical(1 - predict$x[,id_col])
        sparse.x <- x[,id_col]

        if (any(test_id)){
        covariate_predictors <- cbind(1,X.x[train_id,-id_col]%*%main.parameters$u.beta[-id_col,])
        test_predictors <- cbind(1,X.x[test_id,-id_col]%*%main.parameters$u.beta[-id_col,])


        sparse.x[test_id] <- ((test_predictors)%*%(MASS::ginv(t(covariate_predictors)%*%(covariate_predictors))%*%t(covariate_predictors)%*%((as.matrix(x[train_id,id_col])))))
        }
        return(sparse.x)
      }))

      x <- as.matrix(x)

    }

    code$Y_encode <- (main.parameters$alpha.K%*%( y )%*%(main.parameters$v.beta))
    code$X_encode <- (main.parameters$alpha.L%*%( x )%*%(main.parameters$u.beta))

    code$Y_code <- (MASS::ginv((main.parameters$alpha.K)%*%t(main.parameters$alpha.K))%*%(code$Y_encode)%*%MASS::ginv(t(main.parameters$v.beta)%*%(main.parameters$v.beta)))
    code$X_code <- (MASS::ginv((main.parameters$alpha.L)%*%t(main.parameters$alpha.L))%*%(code$X_encode)%*%MASS::ginv(t(main.parameters$u.beta)%*%(main.parameters$u.beta)))

    if (reference == "y"){
      code$main_code <- (code$Y_code)
    }
    if (reference == "x"){
      code$main_code <- (code$X_code)
    }


    main.parameters$alpha.K <- if(is.null(anchors$anchor_y.sample)){t(y%*%t((code$main_code)%*%t(main.parameters$v.beta))%*%MASS::ginv(((code$main_code)%*%t(main.parameters$v.beta))%*%t((code$main_code)%*%t(main.parameters$v.beta))))}else{anchors$anchor_y.sample}
    main.parameters$alpha.L <- if(is.null(anchors$anchor_x.sample)){t(x%*%t((code$main_code)%*%t(main.parameters$u.beta))%*%MASS::ginv(((code$main_code)%*%t(main.parameters$u.beta))%*%t((code$main_code)%*%t(main.parameters$u.beta))))}else{anchors$anchor_x.sample}

    main.parameters$v.beta <- if(is.null(anchors$anchor_y.feature)){t(MASS::ginv(t((t(main.parameters$alpha.K)%*%(code$main_code)))%*%((t(main.parameters$alpha.K)%*%(code$main_code))))%*%t(t(main.parameters$alpha.K)%*%(code$main_code))%*%y)}else{anchors$anchor_y.feature}
    main.parameters$u.beta <- if(is.null(anchors$anchor_x.feature)){t(MASS::ginv(t((t(main.parameters$alpha.L)%*%(code$main_code)))%*%((t(main.parameters$alpha.L)%*%(code$main_code))))%*%t(t(main.parameters$alpha.L)%*%(code$main_code))%*%x)}else{anchors$anchor_x.feature}



    if (reference=="x"){
      if (fixed$i_dim == T){
        main.parameters$alpha.K <- main.parameters$alpha.L
      }
      if (fixed$j_dim == T){
        main.parameters$v.beta <- main.parameters$u.beta
      }
    }


    if (reference=="y"){
      if (fixed$i_dim == T){
        main.parameters$alpha.L <- main.parameters$alpha.K
      }
      if (fixed$j_dim == T){
        main.parameters$u.beta <- main.parameters$v.beta
      }
    }

    # Check convergence

    score_lag <- 5
    accept_score <- 3
    matrix.residuals <- (code$Y_encode - code$X_encode)
    llik.vec <- c(llik.vec, mean(mclust::dmvnorm((matrix.residuals),sigma = diag(diag(t(matrix.residuals)%*%(matrix.residuals)/dim(matrix.residuals)[2])),log = T)))
    score.vec <- c(score.vec, (mean(abs(matrix.residuals))))
    MSE <- mean(tail(score.vec,accept_score))
    prev.MSE <- mean(tail(score.vec,score_lag)[1:accept_score])

    if ( count > ( score_lag ) ){
      if (config$verbose == T){
        print(paste("Iteration: ",count," with Tolerance of: ", abs(prev.MSE - MSE)," and Log-Lik of: ",tail(llik.vec,1),sep=""))
      }
      if ((count > config$max_iter ) | abs(prev.MSE - MSE) < config$tol){
        break
      }
    }


    count = count + 1

  }




  code$K.y_dim_red <- t(main.parameters$alpha.K%*%y)
  code$y.v_dim_red <- y%*%main.parameters$v.beta
  code$L.x_dim_red <- t(main.parameters$alpha.L%*%x)
  code$x.u_dim_red <- x%*%main.parameters$u.beta

  predict$predict.y <- Matrix::Matrix(y*predict$y,sparse = T)
  predict$predict.x <- Matrix::Matrix(x*predict$x,sparse = T)

  return(list(

    main.parameters = main.parameters,

    predict =  predict,

    code = code,

    meta.parameters = list(
      config = config,
      fixed = fixed,
      reference = reference
    ),

    convergence.parameters = list(
      iterations = count,
      score.vec = score.vec,
      llik.vec = llik.vec
    )


  ))


}
