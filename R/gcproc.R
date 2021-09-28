#' Generalised Canonical Procrustes
#'
#' A method that uses a likelihood model to align two datasets via an encoding in a lower dimensional space. The coding of the datasets simultaneously can be used to reconstruct the data to predict missing test points. The parameters can be used to reduce either the feature or the sample dimensions into a smaller subspace for further embedding or projection. To run as default, only requires x and y as inputs.
#'
#' @param x Reference dataset of sample by feature matrix (required)
#' @param y Experimental dataset of sample by feature matrix (required)
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
gcproc <- function(data_list,
                   config = gcproc::extract_config(verbose = F),
                   anchors = gcproc::extract_anchors_framework(verbose = F),
                   pivots = gcproc::extract_pivots_framework(verbose = F),
                   recover = gcproc::extract_recovery_framework(verbose = F),
                   fixed = gcproc::extract_fixed_framework(verbose=F)
                   ){
  runtime.start <- Sys.time()


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
    score.vec <- c()
    score_lag <- 2 # How many previous scores kept track of
    accept_score <- 1 # How many scores used to calculate previous and current "mean score"

    recover$predict.list <- recover$design.list

    initialise.model <- initialise.gcproc(data_list = data_list,
                                          config = config,
                                          anchors = anchors,
                                          pivots = pivots)

    main.parameters <- initialise.model$main.parameters
    code <- initialise.model$code
  }

  if (config$verbose){
    print(paste("Beginning gcproc learning with:    Sample dimension reduction (config$i_dim): ",config$i_dim, "    Feature dimension reduction (config$j_dim): ", config$j_dim,"    Tolerance Threshold: ", config$tol, "   Maximum number of iterations: ", config$max_iter, "   Verbose: ", config$verbose,sep=""))
  }

  while (T){

    prev_code <- code
    for (i in 1:length(data_list)){
      return_update <- update_set(x = data_list[[i]],
                                  main.parameters = main.parameters[[i]],
                                  code = code,
                                  anchors = anchors
                                )

      main.parameters[[i]] <- return_update$main.parameters
      code <- return_update$code

      if (i %in% fixed$alpha | i %in% fixed$beta){

        if (!is.null(fixed$alpha)){

          for (unq_a in unique(fixed$alpha)){
            a_id <- which(fixed$alpha == unq_a)
            main.alpha <- main.parameters[[a_id[1]]]$alpha

            for (a in a_id[-1]){
              main.parameters[[a]]$alpha <- main.alpha
            }
          }
        }

        if (!is.null(fixed$beta)){
          for (unq_b in unique(fixed$beta)){
            b_id <- which(fixed$beta == unq_b)
            main.beta <- main.parameters[[b_id[1]]]$beta

            for (b in b_id[-1]){
              main.parameters[[b]]$beta <- main.beta
            }
          }
        }

      }


      if (any(do.call('c',lapply(recover$design.list,function(X){!is.null(X)})))){

        recover <- recover_points(
          data_list,
          code = code,
          main.parameters = main.parameters,
          config = config,
          recover = recover
        )

        for (i in 1:length(data_list)){
          if (!is.null(recover$predict.list[[i]])){
            data_list[[i]] <- recover$predict.list[[i]]
          }
        }

      }


    }




    matrix.residuals <- code$encode - prev_code$encode

    total.mse <- mean(abs(matrix.residuals))

    # Check convergence
    score.vec <- c(score.vec, total.mse)
    MSE <- mean(tail(score.vec,accept_score))
    prev.MSE <- mean(tail(score.vec,score_lag)[1:accept_score])

    if ( count > ( score_lag ) ){
      if (config$verbose == T){
        print(paste("Iteration: ",count," with Tolerance of: ", abs(prev.MSE - MSE),sep=""))
      }
    }

    if (count > config$min_iter){
      if ((count > config$max_iter ) | abs(prev.MSE - MSE) < config$tol){
        break
      }
    }


    count = count + 1



  }




#   dimension_reduction <- list()
#   dimension_reduction$K.y_dim_red <- t(main.parameters$alpha.K%*%y)
#   dimension_reduction$y.v_dim_red <- y%*%main.parameters$v.beta
#   dimension_reduction$L.x_dim_red <- t(main.parameters$alpha.L%*%x)
#   dimension_reduction$x.u_dim_red <- x%*%main.parameters$u.beta

  # code$Y_decoded <- t(main.parameters$alpha.K)%*%(code$main_code)%*%t(main.parameters$v.beta)
  # code$X_decoded <- t(main.parameters$alpha.L)%*%(code$main_code)%*%t(main.parameters$u.beta)

  # code$main_code <- t(main.parameters$s_code)%*%main.parameters$s_code

  runtime.end <- Sys.time()

  return(list(

    main.parameters = main.parameters,

    recover =  recover,

    # dimension_reduction = dimension_reduction,

    code = code,

    meta.parameters = list(
      config = config,
      fixed = fixed,
      runtime = list(
        runtime.start = runtime.start,
        runtime.end = runtime.end,
        runtime_total = runtime.end - runtime.start
      )
    ),

    convergence.parameters = list(
      iterations = count,
      score.vec = score.vec
    )


  ))


}





update_set <- function(x,
                       main.parameters,
                       code,
                       anchors
                       ){

  main.parameters$alpha <- if(is.null(anchors$anchor_x.sample)){t(x%*%t((code$main_code)%*%t(main.parameters$beta))%*%MASS::ginv(((code$main_code)%*%t(main.parameters$beta))%*%t((code$main_code)%*%t(main.parameters$beta))))}else{anchors$anchor_x.sample}
  main.parameters$beta <- if(is.null(anchors$anchor_x.feature)){t(MASS::ginv(t((t(main.parameters$alpha)%*%(code$main_code)))%*%((t(main.parameters$alpha)%*%(code$main_code))))%*%t(t(main.parameters$alpha)%*%(code$main_code))%*%x)}else{anchors$anchor_x.feature}

  code$encode <- (main.parameters$alpha%*%( x )%*%(main.parameters$beta))
  main.parameters$r.s_code <- t(MASS::ginv(t(code$code)%*%(code$code))%*%t(code$code)%*%MASS::ginv(main.parameters$l.s_code%*%t(main.parameters$l.s_code))%*%main.parameters$l.s_code%*%(MASS::ginv((main.parameters$alpha)%*%t(main.parameters$alpha))%*%(code$encode)%*%MASS::ginv(t(main.parameters$beta)%*%(main.parameters$beta))))
  main.parameters$l.s_code <- t((MASS::ginv((main.parameters$alpha)%*%t(main.parameters$alpha))%*%(code$encode)%*%MASS::ginv(t(main.parameters$beta)%*%(main.parameters$beta)))%*%main.parameters$r.s_code%*%MASS::ginv(t(main.parameters$r.s_code)%*%(main.parameters$r.s_code))%*%t(code$code)%*%MASS::ginv((code$code)%*%t(code$code)))

  code$main_code <- MASS::ginv((main.parameters$alpha)%*%t(main.parameters$alpha))%*%(code$encode)%*%MASS::ginv(t(main.parameters$beta)%*%(main.parameters$beta))
  code$code <- MASS::ginv((main.parameters$l.s_code)%*%t(main.parameters$l.s_code))%*%(main.parameters$l.s_code)%*%code$main_code%*%(main.parameters$r.s_code)%*%MASS::ginv(t(main.parameters$r.s_code)%*%(main.parameters$r.s_code))


  return(list(main.parameters = main.parameters,
              code = code
              ))

}


