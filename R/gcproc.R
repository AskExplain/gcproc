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
#' @param recover Important information for prediction or imputation (not required)
#'
#' @return Main parameters contains the learned model parameters. The alpha and beta matrix multiply x and y by, (K)(Y)(v) and (L)(X)(u). By multiplying with the parameter, the dimension of the samples and features can be dimensionally reduced for further visualisation analysis such as embedding or projection.
#'
#' @return Code contains the learned shared encoding space. The encoded space refers to the full dimension reduction of both samples and features after matrix multiplication by parameters K and v for y, as well as, L and u for x. The decode is an estimation of the full matrix dataset, where the code is used and matrix multiplied as t(K)(Y_code)t(v), and t(L)(X_code)t(u) to calculate the decoded estimation.
#'
#' @return Recover contains the predictions for the test dataset as indicated by a 1 in the binary prediction matrices. In addition to the x and y variables within this list that correspond to the input binary design matrix, the predict.y and predict.x datasets represent the decoded estimation, and predicts any missing test data.
#'
#' @export
gcproc <- function(x,
                   y,
                   reference = "y",
                   fixed = list(i_dim = F, j_dim = F),
                   config = gcproc::extract_config(verbose = F),
                   anchors = gcproc::extract_anchors_framework(verbose = F),
                   pivots = gcproc::extract_pivots_framework(verbose = F),
                   recover = gcproc::extract_recovery_framework(verbose = F)
){



  # if (is.null(config$i_dim)){
  #   eig.val_x <- RSpectra::eigs(crossprod(t(x),t(x)),k = min(nrow(x)-2,30))$values
  #   eig.val_y <- RSpectra::eigs(crossprod(t(y),t(y)),k = min(nrow(y)-2,30))$values
  #
  #   main_dim <- sum(cumsum(((eig.val_x)/sum((eig.val_x))+(eig.val_y)/sum((eig.val_y)))/2)<(1-1e-2))
  #
  #   config$i_dim <- max(5,main_dim)
  # }
  #
  #
  # if (is.null(config$j_dim)){
  #   eig.val_x <- RSpectra::eigs(crossprod((x),(x)),k = min(ncol(x)-2,30))$values
  #   eig.val_y <- RSpectra::eigs(crossprod((y),(y)),k = min(ncol(y)-2,30))$values
  #
  #   main_dim <- sum(cumsum(((eig.val_x)/sum((eig.val_x))+(eig.val_y)/sum((eig.val_y)))/2)<(1-1e-2))
  #
  #   config$j_dim <- max(5,main_dim)
  # }

  config$init.i_dim <- config$i_dim
  config$init.j_dim <- config$j_dim


  prepare_data = TRUE
  initialise = TRUE


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

  if (initialise==T){

    # Prepare convergence checking parameters
    count = 1
    llik.vec <- c()
    score.vec <- c()
    score_lag <- 5 # How many previous scores kept track of
    accept_score <- 4 # How many scores used to calculate previous and current "mean score"

    recover$predict.y <- recover$y
    recover$predict.x <- recover$x

    initial.parameters <- update.parameters(  x=x,
                                               y=y,
                                               config=config,
                                               reference=reference,
                                               fixed=fixed,
                                               anchors=anchors,
                                               pivots=pivots)

    config <- initial.parameters$config
    all.parameters <- initial.parameters$parameters$all.parameters
    main.parameters <- initial.parameters$parameters$main.parameters
    code <- initial.parameters$code

  }

  if (config$verbose){
    print(paste("Beginning gcproc learning with:    Sample dimension reduction (config$i_dim): ",config$i_dim, "    Feature dimension reduction (config$j_dim): ", config$j_dim,"    Tolerance Threshold: ", config$tol, "   Maximum number of iterations: ", config$max_iter, "   Verbose: ", config$verbose,sep=""))
  }


  while (T){


    if (!is.null(recover$y) | !is.null(recover$x) ){
      recover <- recover_points(
        x = x,
        y = y,
        main.parameters = main.parameters,
        config = config,
        recover = recover,
        fixed = fixed)

      if (!is.null(recover$x)){
        x <- recover$predict.x
      }
      if (!is.null(recover$y)){
        y <- recover$predict.y
      }
    }


    new.update.parameters <- update.parameters(  x=x,
                                              y=y,
                                              main.all.parameters = all.parameters,
                                              code = code,
                                              config=config,
                                              reference=reference,
                                              fixed=fixed,
                                              anchors=anchors,
                                              pivots=pivots)

    config <- new.update.parameters$config
    all.parameters <- new.update.parameters$parameters$all.parameters
    main.parameters <- new.update.parameters$parameters$main.parameters
    code <- new.update.parameters$code

    # Check convergence
    matrix.residuals <- (code[[config$layers]]$Y_code - code[[config$layers]]$X_code)
    llik.vec <- c(llik.vec, mean(mclust::dmvnorm((matrix.residuals),sigma = diag(diag(t(matrix.residuals)%*%(matrix.residuals)/dim(matrix.residuals)[2])),log = T)))
    score.vec <- c(score.vec, (mean(abs(matrix.residuals))))
    MSE <- mean(tail(score.vec,accept_score))
    prev.MSE <- mean(tail(score.vec,score_lag)[1:accept_score])

    if ( count > ( score_lag ) ){
      if (config$verbose == T){
        print(paste("Iteration: ",count," with Tolerance of: ", abs(prev.MSE - MSE)," and Log-Lik of: ",tail(llik.vec,1),sep=""))
      }
      if (count > 5){
        if ((count > config$max_iter ) | abs(prev.MSE - MSE) < config$tol){
          break
        }
      }
    }


    count = count + 1

  }

  dimension_reduction <- list()
  dimension_reduction$K.y_dim_red <- t(main.parameters$alpha.K%*%y)
  dimension_reduction$y.v_dim_red <- y%*%main.parameters$v.beta
  dimension_reduction$L.x_dim_red <- t(main.parameters$alpha.L%*%x)
  dimension_reduction$x.u_dim_red <- x%*%main.parameters$u.beta

  # recover$predict.x <- Matrix::Matrix(recover$predict.x*recover$x,sparse=T)
  # recover$predict.y <- Matrix::Matrix(recover$predict.y*recover$y,sparse=T)

  return(list(

    main.parameters = main.parameters,

    recover =  recover,

    dimension_reduction = dimension_reduction,

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



update.parameters <- function(x,
                              y,
                              main.all.parameters = NULL,
                              code = NULL,
                              config,
                              reference,
                              fixed,
                              anchors,
                              pivots){

  side.config <- config

  if (is.null(main.all.parameters)){
    all.parameters <- list()
  } else {
    all.parameters <- main.all.parameters
  }


  if (is.null(code)){
    all.code <- list()
  } else {
    all.code <- code
  }


  for (layer in c(1:config$layers)){

    side.config$i_dim <- side.config$all.i_dim[layer]
    side.config$j_dim <- side.config$all.j_dim[layer]

    if (is.null(code) & is.null(main.all.parameters)){
      initialise.model <- initialise.gcproc(x = x,
                                            y = y,
                                            fixed = fixed,
                                            reference = reference,
                                            config = side.config,
                                            anchors = anchors,
                                            pivots = pivots)

      all.parameters[[layer]] <- initialise.model$main.parameters
      all.code[[layer]] <- initialise.model$code
    }

    all.code[[layer]]$Y_encode <- (all.parameters[[layer]]$alpha.K%*%( y )%*%(all.parameters[[layer]]$v.beta))
    all.code[[layer]]$X_encode <- (all.parameters[[layer]]$alpha.L%*%( x )%*%(all.parameters[[layer]]$u.beta))

    all.code[[layer]]$Y_code <- (MASS::ginv((all.parameters[[layer]]$alpha.K)%*%t(all.parameters[[layer]]$alpha.K))%*%(all.code[[layer]]$Y_encode)%*%MASS::ginv(t(all.parameters[[layer]]$v.beta)%*%(all.parameters[[layer]]$v.beta)))
    all.code[[layer]]$X_code <- (MASS::ginv((all.parameters[[layer]]$alpha.L)%*%t(all.parameters[[layer]]$alpha.L))%*%(all.code[[layer]]$X_encode)%*%MASS::ginv(t(all.parameters[[layer]]$u.beta)%*%(all.parameters[[layer]]$u.beta)))

    if (reference == "y"){
      all.code[[layer]]$main_code <- all.code[[layer]]$Y_code
    }
    if (reference == "x"){
      all.code[[layer]]$main_code <- all.code[[layer]]$X_code
    }

    if (reference=="x"){
      all.parameters[[layer]]$alpha.L <- if(is.null(anchors$anchor_x.sample)){t(x%*%t((all.code[[layer]]$main_code)%*%t(all.parameters[[layer]]$u.beta))%*%MASS::ginv(((all.code[[layer]]$main_code)%*%t(all.parameters[[layer]]$u.beta))%*%t((all.code[[layer]]$main_code)%*%t(all.parameters[[layer]]$u.beta))))}else{anchors$anchor_x.sample}
      all.parameters[[layer]]$u.beta <- if(is.null(anchors$anchor_x.feature)){t(MASS::ginv(t((t(all.parameters[[layer]]$alpha.L)%*%(all.code[[layer]]$main_code)))%*%((t(all.parameters[[layer]]$alpha.L)%*%(all.code[[layer]]$main_code))))%*%t(t(all.parameters[[layer]]$alpha.L)%*%(all.code[[layer]]$main_code))%*%x)}else{anchors$anchor_x.feature}

      if (fixed$i_dim == T){
        all.parameters[[layer]]$alpha.K <- all.parameters[[layer]]$alpha.L
      } else {
        all.parameters[[layer]]$alpha.K <- if(is.null(anchors$anchor_y.sample)){t(y%*%t((all.code[[layer]]$main_code)%*%t(all.parameters[[layer]]$v.beta))%*%MASS::ginv(((all.code[[layer]]$main_code)%*%t(all.parameters[[layer]]$v.beta))%*%t((all.code[[layer]]$main_code)%*%t(all.parameters[[layer]]$v.beta))))}else{anchors$anchor_y.sample}
      }
      if (fixed$j_dim == T){
        all.parameters[[layer]]$v.beta <- all.parameters[[layer]]$u.beta
      } else {
        all.parameters[[layer]]$v.beta <- if(is.null(anchors$anchor_y.feature)){t(MASS::ginv(t((t(all.parameters[[layer]]$alpha.K)%*%(all.code[[layer]]$main_code)))%*%((t(all.parameters[[layer]]$alpha.K)%*%(all.code[[layer]]$main_code))))%*%t(t(all.parameters[[layer]]$alpha.K)%*%(all.code[[layer]]$main_code))%*%y)}else{anchors$anchor_y.feature}
      }
    }


    if (reference=="y"){
      all.parameters[[layer]]$alpha.K <- if(is.null(anchors$anchor_y.sample)){t(y%*%t((all.code[[layer]]$main_code)%*%t(all.parameters[[layer]]$v.beta))%*%MASS::ginv(((all.code[[layer]]$main_code)%*%t(all.parameters[[layer]]$v.beta))%*%t((all.code[[layer]]$main_code)%*%t(all.parameters[[layer]]$v.beta))))}else{anchors$anchor_y.sample}
      all.parameters[[layer]]$v.beta <- if(is.null(anchors$anchor_y.feature)){t(MASS::ginv(t((t(all.parameters[[layer]]$alpha.K)%*%(all.code[[layer]]$main_code)))%*%((t(all.parameters[[layer]]$alpha.K)%*%(all.code[[layer]]$main_code))))%*%t(t(all.parameters[[layer]]$alpha.K)%*%(all.code[[layer]]$main_code))%*%y)}else{anchors$anchor_y.feature}

      if (fixed$i_dim == T){
        all.parameters[[layer]]$alpha.L <- all.parameters[[layer]]$alpha.K
      } else {
        all.parameters[[layer]]$alpha.L <- if(is.null(anchors$anchor_x.sample)){t(x%*%t((all.code[[layer]]$main_code)%*%t(all.parameters[[layer]]$u.beta))%*%MASS::ginv(((all.code[[layer]]$main_code)%*%t(all.parameters[[layer]]$u.beta))%*%t((all.code[[layer]]$main_code)%*%t(all.parameters[[layer]]$u.beta))))}else{anchors$anchor_x.sample}
      }
      if (fixed$j_dim == T){
        all.parameters[[layer]]$u.beta <- all.parameters[[layer]]$v.beta
      } else {
        all.parameters[[layer]]$u.beta <- if(is.null(anchors$anchor_x.feature)){t(MASS::ginv(t((t(all.parameters[[layer]]$alpha.L)%*%(all.code[[layer]]$main_code)))%*%((t(all.parameters[[layer]]$alpha.L)%*%(all.code[[layer]]$main_code))))%*%t(t(all.parameters[[layer]]$alpha.L)%*%(all.code[[layer]]$main_code))%*%x)}else{anchors$anchor_x.feature}
      }
    }

    x <- all.code[[layer]]$X_code
    y <- all.code[[layer]]$Y_code

  }


  main.parameters <- list()

  for (layer in c(1:config$layers)){

    if (layer == 1){
      main.parameters$alpha.K <- all.parameters[[layer]]$alpha.K
      main.parameters$alpha.L <- all.parameters[[layer]]$alpha.L
      main.parameters$v.beta <- all.parameters[[layer]]$v.beta
      main.parameters$u.beta <- all.parameters[[layer]]$u.beta
    }
    if (layer > 1){
      main.parameters$alpha.K <- all.parameters[[layer]]$alpha.K%*%main.parameters$alpha.K
      main.parameters$alpha.L <- all.parameters[[layer]]$alpha.L%*%main.parameters$alpha.L
      main.parameters$v.beta <- main.parameters$v.beta%*%all.parameters[[layer]]$v.beta
      main.parameters$u.beta <- main.parameters$u.beta%*%all.parameters[[layer]]$u.beta
    }


  }


  return(list(

    config = config,

    code = all.code,

    parameters = list (

      main.parameters = main.parameters,

      all.parameters = all.parameters

    )

  ))


}

