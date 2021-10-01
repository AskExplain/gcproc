#' Generalised Canonical Procrustes
#'
#' A method that uses a likelihood model to align two datasets via an encoding in a lower dimensional space. The coding of the datasets simultaneously can be used to reconstruct the data to predict missing test points. The parameters can be used to reduce either the feature or the sample dimensions into a smaller subspace for further embedding or projection. To run as default, only data list is required.
#'
#' @param data_list List of data matrices of varying dimensionality. Attempts to find similarities among all datasets with a core structure.
#' @param config Configuration parameters (required, default provided)
#' @param anchors Transferring pre-trained model parameters (not required)
#' @param pivots Initialisation of model parameters (not required)
#' @param recover Important information for prediction or imputation (not required)
#' @param fixed Constrain parameters that share the same axes to be similar (not required)
#'
#' @return Main parameters contains the learned model parameters. The alpha and beta matrix multiply example datasets x and y by, (K)(Y)(v) and (L)(X)(u). By multiplying with the parameter, the dimension of the samples and features can be dimensionally reduced for further visualisation analysis such as embedding or projection.
#'
#' @return Code contains the learned shared encoding space. The encoded space refers to the full dimension reduction of both samples and features after matrix multiplication by parameters K and v for y, as well as, L and u for x. The decode is an estimation of the full matrix dataset, where the code is used and matrix multiplied as t(K)(Y_code)t(v), and t(L)(X_code)t(u) to calculate the decoded estimation.
#'
#' @return Recover contains the list of predictions for the test dataset as indicated by a 1 in the binary prediction matrices.
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


  initialise = TRUE

  if (initialise==T){

    # Prepare convergence checking parameters
    count = 1
    score.vec <- c()
    score_lag <- 2 # How many previous scores kept track of
    accept_score <- 1 # How many scores used to calculate previous and current "mean score"

    recover$predict.list <- lapply(c(1:length(data_list)),function(X){NULL})

    initialise.model <- initialise.gcproc(data_list = data_list,
                                          config = config,
                                          anchors = anchors)

    main.parameters <- initialise.model$main.parameters
    code <- if(is.null(pivots$code)){initialise.model$code}else{pivots$code}
  }

  if (config$verbose){
    print(paste("Beginning gcproc learning with:    Sample dimension reduction (config$i_dim): ",config$i_dim, "    Feature dimension reduction (config$j_dim): ", config$j_dim,"    Tolerance Threshold: ", config$tol, "   Maximum number of iterations: ", config$max_iter, "   Verbose: ", config$verbose,sep=""))
  }

  while (T){

    prev_code <- code
    for (i in 1:length(data_list)){

      return_update <- update_set(x = as.matrix(data_list[[i]]),
                                  main.parameters = main.parameters[[i]],
                                  code = code
                                  )

      main.parameters[[i]] <- return_update$main.parameters
      code <- if(is.null(anchors$code)){return_update$code}else{anchors$code}


      if (i %in% fixed$alpha | i %in% fixed$beta){

        if (!is.null(fixed$alpha)){

          a_id <- which(fixed$alpha == fixed$alpha[i])
          main.alpha <- main.parameters[[i]]$alpha

          for (a in a_id){
            main.parameters[[a]]$alpha <- main.alpha
          }
        }

        if (!is.null(fixed$beta)){

          b_id <- which(fixed$beta == unique(fixed$beta)[i])
          main.beta <- main.parameters[[b_id[1]]]$beta

          for (b in b_id){
            main.parameters[[b]]$beta <- main.beta
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
    } else {
      print(paste("Iteration: ",count," ... initialising ... ",sep=""))
    }

    if (count > config$min_iter){
      if ((count > config$max_iter ) | abs(prev.MSE - MSE) < config$tol){
        break
      }
    }

    count = count + 1

  }



  if (any(do.call('c',lapply(recover$design.list,function(X){!is.null(X)})))){

    recover <- recover_points(
      data_list,
      code = code,
      main.parameters = main.parameters,
      config = config,
      recover = recover
    )

  }



  dimension_reduction <- lapply(c(1:length(data_list)),function(X){
    feature_x.dim_reduce.encode <- t(main.parameters[[X]]$alpha%*%data_list[[X]])
    sample_x.dim_reduce.encode <- data_list[[X]]%*%main.parameters[[X]]$beta

    feature_x.dim_reduce.code <- t(MASS::ginv((main.parameters[[X]]$alpha)%*%t(main.parameters[[X]]$alpha))%*%main.parameters[[X]]$alpha%*%data_list[[X]])
    sample_x.dim_reduce.code <- data_list[[X]]%*%main.parameters[[X]]$beta%*%MASS::ginv(t(main.parameters[[X]]$beta)%*%(main.parameters[[X]]$beta))

    return(list(
      feature_x.dim_reduce.encode = feature_x.dim_reduce.encode,
      sample_x.dim_reduce.encode = sample_x.dim_reduce.encode,
      feature_x.dim_reduce.code = feature_x.dim_reduce.code,
      sample_x.dim_reduce.code = sample_x.dim_reduce.code
    ))
  })

  runtime.end <- Sys.time()

  return(list(

    main.parameters = main.parameters,

    code = code,

    recover =  recover,

    dimension_reduction = dimension_reduction,

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


  alpha.batch <- chunk(c(1:config$i_dim),config$batch)
  beta.batch <- chunk(c(1:config$j_dim),config$batch)



  main.parameters$alpha <- as.matrix(Reduce('+',parallel::mclapply(mc.cores = config$n_cores, alpha.batch,function(a){

    internal_alpha <- main.parameters$alpha
    internal_alpha[-a,] <- 0
    internal_alpha <- Matrix::Matrix(internal_alpha,sparse=T)

    main_decode_beta <- ((code$decode)%*%t(main.parameters$beta))[a,]
    internal_alpha[a,] <- t(x%*%t(main_decode_beta)%*%MASS::ginv((main_decode_beta)%*%t(main_decode_beta)))

    return(internal_alpha)
  })))




  main.parameters$beta <- as.matrix(Reduce('+',parallel::mclapply(mc.cores = config$n_cores, beta.batch,function(b){

    internal_beta <- main.parameters$beta
    internal_beta[,-b] <- 0
    internal_beta <- Matrix::Matrix(internal_beta,sparse=T)


    main_decode_beta <- (t(main.parameters$alpha)%*%(code$decode))[,b]
    internal_beta[,b] <- t(MASS::ginv(t((main_decode_beta))%*%((main_decode_beta)))%*%t(main_decode_beta)%*%x)

    return(internal_beta)
  })))





  all_batch.id <- cbind(rep(c(1:length(alpha.batch)),length(beta.batch)),
                        rep(c(1:length(beta.batch)),length(alpha.batch)))

  internal.code <- parallel::mclapply(mc.cores = config$n_cores, c(1:dim(all_batch.id)[1]),function(c){

    a <- alpha.batch[[all_batch.id[c,1]]]
    b <- beta.batch[[all_batch.id[c,2]]]

    internal_encode <- code$encode
    internal_encode <- Matrix::Matrix(0,nrow = config$i_dim, ncol = config$j_dim, sparse=T)

    internal_decode <- code$decode
    internal_decode <- Matrix::Matrix(0,nrow = config$i_dim, ncol = config$j_dim, sparse=T)

    internal_encode[a,b] <- (main.parameters$alpha[a,]%*%( x )%*%(main.parameters$beta[,b]))
    internal_decode[a,b] <- MASS::ginv((main.parameters$alpha[a,])%*%t(main.parameters$alpha[a,]))%*%internal_encode[a,b]%*%MASS::ginv(t(main.parameters$beta[,b])%*%(main.parameters$beta[,b]))

    return(list(encode=internal_encode,
                decode=internal_decode))
  })

  code$encode <- as.matrix(Reduce('+',lapply(internal.code,function(X){X$encode})))
  code$decode <- as.matrix(Reduce('+',lapply(internal.code,function(X){X$decode})))

  return(list(main.parameters = main.parameters,
              code = code
              ))

}


