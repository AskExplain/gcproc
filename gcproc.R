#' Generalised Canonical Procrustes
#'
#' A method that uses a likelihood model to align multiple datasets via an encoding in a lower dimensional space. The parameters can be used to reduce either the feature or the sample dimensions into a smaller subspace for further embedding or prediction. To run as default, only a data list is required - please review the config parameters at gcproc::extract_config(T)  .
#'
#' @param data_list List of data matrices of varying dimensionality. Attempts to find similarities among all datasets with a core structure.
#' @param config Configuration parameters (required, default provided)
#' @param transfer Transferring pre-trained model parameters (not required)
#' @param recover Important information used for prediction or imputation (not required)
#' @param join Join parameters that share the same axes to be similar (not required)
#'
#' @return Main parameters contains the learned model parameters. The alpha and beta matrix multiply example datasets x and y by, (K)(Y)(v) and (L)(X)(u). By multiplying with the parameter, the dimension of the samples and features can be dimensionally reduced for further visualisation analysis such as embedding or projection.
#'
#' @return Code contains the learned shared encoding space. The encoded space refers to the full dimension reduction of both samples and features after matrix multiplication by parameters K and v for y, as well as, L and u for x. The decode is an estimation of the full matrix dataset, where the code is used and matrix multiplied as t(K)(Y_code)t(v), and t(L)(X_code)t(u) to calculate the decoded estimation.
#'
#' @return Recover contains the list of predictions for the test dataset as indicated by a 1 in the binary prediction matrices. The prediction occurs in the shared lower dimensional space where all data sets in the list are projected to using a common latent code.
#'
#' @export
gcproc <- function(data_list,
                   config = gcproc::extract_config(verbose = F),
                   transfer = gcproc::extract_transfer_framework(verbose = F),
                   recover = gcproc::extract_recovery_framework(verbose = F),
                   join = gcproc::extract_join_framework(verbose=F)
){

  runtime.start <- Sys.time()

  set.seed(config$seed)

  initialise = TRUE

  # Prepare convergence checking parameters
  count = 1
  score.vec <- c()
  score_lag <- 2 # How many previous scores kept track of
  accept_score <- 1 # How many scores used to calculate previous and current "mean score"

  recover$predict.list <- lapply(c(1:length(data_list)),function(X){NULL})

  initialise.model <- initialise.gcproc(data_list = data_list,
                                        config = config,
                                        transfer = transfer)

  main.parameters <- initialise.model$main.parameters
  code <- initialise.model$code

  if (config$verbose){
    print(paste("Beginning gcproc learning with:    Sample dimension reduction (config$i_dim): ",config$i_dim, "    Feature dimension reduction (config$j_dim): ", config$j_dim,"    Tolerance Threshold: ", config$tol, "   Maximum number of iterations: ", config$max_iter, "   Verbose: ", config$verbose, sep=""))
  }

  while (T){

    prev_code <- code

    for (i in 1:length(data_list)){


      return_update <- update_set(x = as.matrix(data_list[[i]]),
                                  main.parameters = main.parameters[[i]],
                                  code = code,
                                  transfer = transfer
      )

      main.parameters[[i]] <- return_update$main.parameters
      code <- return_update$code

      if (!is.null(join$alpha)){

        a_id <- which(join$alpha == join$alpha[i])
        main.alpha <- main.parameters[[i]]$alpha

        for (a in a_id){
          main.parameters[[a]]$alpha <- main.alpha
        }
      }

      if (!is.null(join$beta)){

        b_id <- which(join$beta == unique(join$beta)[i])
        main.beta <- main.parameters[[b_id[1]]]$beta

        for (b in b_id){
          main.parameters[[b]]$beta <- main.beta
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
    } else {
      if (config$verbose){
        print(paste("Iteration: ",count," ... initialising ... ",sep=""))
      }
    }

    if (count > config$min_iter){
      if ((count > config$max_iter ) | abs(prev.MSE - MSE) < config$tol){
        break
      }
    }

    count = count + 1

  }

  if (any(do.call('c',lapply(recover$design.list,function(X){!is.null(X)})))){

    recover_data <- recover_points(
      data_list,
      code = code,
      main.parameters = main.parameters,
      config = config,
      recover = recover
    )

    recover <- recover_data$recover
    data_list <- recover_data$data_list
  }





  if (config$verbose){
    print("Learning has converged for gcproc, beginning prediction (if requested) and dimension reduction")
  }




  dimension_reduction <- lapply(c(1:length(data_list)),function(Y){

    x <- as.matrix(data_list[[Y]])

    feature_x.dim_reduce.encode <- t(main.parameters[[Y]]$alpha%*%x)
    sample_x.dim_reduce.encode <- x%*%main.parameters[[Y]]$beta

    return(list(
      feature_x.dim_reduce.encode = feature_x.dim_reduce.encode,
      sample_x.dim_reduce.encode = sample_x.dim_reduce.encode
    ))
  })

  runtime.end <- Sys.time()


  if (config$verbose){
    print(paste("Done! Total runtime of   ", runtime.end - runtime.start ,sep=""))
  }

  return(list(

    main.parameters = main.parameters,

    code = code,

    recover =  recover,

    dimension_reduction = dimension_reduction,

    meta.parameters = list(
      config = config,
      join = join,
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
                       transfer = NULL
){

  main.parameters$alpha <- t(x%*%t((code$code)%*%t(main.parameters$beta))%*%MASS::ginv(((code$code)%*%t(main.parameters$beta))%*%t((code$code)%*%t(main.parameters$beta))))
  main.parameters$beta <- t(MASS::ginv(t((t(main.parameters$alpha)%*%(code$code)))%*%((t(main.parameters$alpha)%*%(code$code))))%*%t(t(main.parameters$alpha)%*%(code$code))%*%x)

  code$encode <- (main.parameters$alpha%*%( x )%*%(main.parameters$beta))
  code$code <- MASS::ginv((main.parameters$alpha)%*%t(main.parameters$alpha))%*%code$encode%*%MASS::ginv(t(main.parameters$beta)%*%(main.parameters$beta))

  return(list(main.parameters = main.parameters,
              code = code
  ))

}
