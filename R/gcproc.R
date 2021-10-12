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
                   covariate = gcproc::extract_covariate_framework(verbose=F),
                   transfer = gcproc::extract_transfer_framework(verbose = F),
                   recover = gcproc::extract_recovery_framework(verbose = F),
                   join = gcproc::extract_join_framework(verbose=F)
){

  runtime.start <- Sys.time()

  set.seed(config$config$seed)

  initial_setup <- TRUE

  if (is.null(covariate$factor)){

    covariate$factor <- data.frame(rep("ALL",length(data_list)))

  } else {

    covariate$factor <- data.frame(cbind("ALL",covariate$factor))

  }



  # Prepare convergence checking parameters
  count = 0
  update = 0
  score.vec <- c()

  score_lag <- 2 # How many previous scores kept track of
  accept_score <- 1 # How many scores used to calculate previous and current "mean score"

  internal_anchor.update <- 1
  internal_bootstrap <- 1
  main_sample_seed <- 1

  convergence.parameters <- list(count = count, score.vec = score.vec)


  recover$predict.list <- lapply(c(1:length(data_list)),function(X){NULL})


  initialise.model <- initialise.gcproc(data_list = data_list,
                                        config = config,
                                        covariate = covariate,
                                        transfer = transfer,
                                        join = join)

  main.parameters <- initialise.model$main.parameters
  main.code <- initialise.model$main.code
  main.index <- initialise.model$main.index




  if (config$verbose){
    print(paste("Beginning gcproc learning with:    Sample dimension reduction (config$i_dim): ",config$i_dim, "    Feature dimension reduction (config$j_dim): ", config$j_dim,"    Tolerance Threshold: ", config$tol, "   Maximum number of iterations: ", config$max_iter, "   Verbose: ", config$verbose, sep=""))
  }


  while (T){

    prev_code = main.code

    for (j in 1:length(data_list)){

      index <- main.index[[j]]


      internal.param <- list(
        alpha = main.parameters$alpha[[join$alpha[j]]],
        beta = main.parameters$beta[[join$beta[j]]]
      )

      gcproc.update <- update_set(x = as.matrix(data_list[[j]]),
                                  main.parameters = internal.param,
                                  main.code = main.code,
                                  fix = transfer$fix,
                                  index = index)

      main.parameters$alpha[[join$alpha[j]]] <- gcproc.update$main.parameters$alpha
      main.parameters$beta[[join$beta[j]]] <- gcproc.update$main.parameters$beta

      main.code <- gcproc.update$main.code

    }


    print_iterations(convergence.parameters,config)

    matrix.residuals <- prev_code$encode - main.code$encode

    total.mae <- mean(abs(matrix.residuals))

    # Check convergence
    convergence.parameters$score.vec <- c(convergence.parameters$score.vec, total.mae)
    convergence.parameters$MAE <- mean(tail(convergence.parameters$score.vec,accept_score))
    convergence.parameters$prev.MAE <- mean(tail(convergence.parameters$score.vec,score_lag)[1:accept_score])

    if (convergence.parameters$count > config$min_iter &  convergence.parameters$count > ( score_lag ) ){
      if ((convergence.parameters$count > config$max_iter ) | abs(convergence.parameters$prev.MAE - convergence.parameters$MAE) < config$tol){

        break

      }
    }

    convergence.parameters$count <- convergence.parameters$count + 1
    initial_setup = T

  }

  if (config$verbose){
    print("Learning has converged for gcproc, beginning prediction (if requested) and dimension reduction")
  }

  if (any(do.call('c',lapply(recover$design.list,function(X){!is.null(X)})))){

    recover_data <- recover_points(
      data_list,
      main.code = main.code,
      main.parameters = main.parameters,
      config = config,
      recover = recover,
      join = join
    )

    recover <- recover_data$recover
    data_list <- recover_data$data_list

  }


  dimension_reduction <- lapply(c(1:length(data_list)),function(Y){

    x <- as.matrix(data_list[[Y]])

    feature_x.dim_reduce.encode <- t(main.parameters$alpha[[join$alpha[Y]]]%*%x)
    sample_x.dim_reduce.encode <- x%*%main.parameters$beta[[join$beta[Y]]]

    feature_x.dim_reduce.code <- t(pinv(t(main.parameters$alpha[[join$alpha[Y]]]))%*%main.parameters$alpha[[join$alpha[Y]]]%*%x)
    sample_x.dim_reduce.code <- x%*%main.parameters$beta[[join$beta[Y]]]%*%pinv((main.parameters$beta[[join$beta[Y]]]))

    return(list(
      feature_x.dim_reduce.encode = feature_x.dim_reduce.encode,
      sample_x.dim_reduce.encode = sample_x.dim_reduce.encode,
      feature_x.dim_reduce.code = feature_x.dim_reduce.code,
      sample_x.dim_reduce.code = sample_x.dim_reduce.code
    ))
  })




  runtime.end <- Sys.time()


  if (config$verbose){
    print(paste("Done! Total runtime of   ", runtime.end - runtime.start ,sep=""))
  }

  return(list(

    main.parameters = main.parameters,

    main.code = main.code,

    recover =  recover,

    dimension_reduction = dimension_reduction,

    meta.parameters = list(
      config = config,
      join = join,
      main.index = main.index,
      runtime = list(
        runtime.start = runtime.start,
        runtime.end = runtime.end,
        runtime_total = runtime.end - runtime.start
      )
    ),

    convergence.parameters = convergence.parameters


  ))


}


update_set <- function(x,
                       main.parameters,
                       main.code,
                       index,
                       fix){

  internal.code <- Reduce('+',lapply(c(1:length(main.code$code)),function(X){
    index[X,2] * main.code$code[[X]]
  }))/sum(index[,2])

  main.parameters$alpha <- (t(x%*%t((internal.code)%*%t(main.parameters$beta))%*%pinv(t((internal.code)%*%t(main.parameters$beta)))))
  main.parameters$beta <- (t(pinv(((t(main.parameters$alpha)%*%(internal.code))))%*%t(t(main.parameters$alpha)%*%(internal.code))%*%x))

  main.code$encode <- (main.parameters$alpha%*%( x )%*%(main.parameters$beta))

  for (i in which(index[,2]==1)){

    if(!fix){
      main.code$code[[i]] <- pinv(t(main.parameters$alpha))%*%(main.code$encode)%*%pinv(main.parameters$beta)
    }

  }

  return(list(main.parameters = main.parameters,
              main.code = main.code
  ))

}


pinv <- function(X){
  MASS::ginv(t(X)%*%X)
}


chunk <- function(x,n){
  if (n==1){
    list(x)
  }
  else{
    split(x, cut(seq_along(x), n, labels = FALSE))
  }
}


print_iterations <- function(convergence.parameters,config){

  if (convergence.parameters$count <= 15){
    if (config$verbose){
      print(paste("Iteration: ",convergence.parameters$count,sep=""))
    }
  }
  if (convergence.parameters$count > 15 & convergence.parameters$count <= 100){
    if (convergence.parameters$count%%5 == 0 & config$verbose){
      print(paste("Iteration: ",convergence.parameters$count,sep=""))
    }
  }
  if (convergence.parameters$count > 100 & convergence.parameters$count <= 200){
    if (convergence.parameters$count%%10 == 0 & config$verbose){
      print(paste("Iteration: ",convergence.parameters$count,sep=""))
    }
  }
  if (convergence.parameters$count > 200 & convergence.parameters$count <= 300){
    if (convergence.parameters$count%%25 == 0 & config$verbose){
      print(paste("Iteration: ",convergence.parameters$count,sep=""))
    }
  }
  if (convergence.parameters$count > 300 & convergence.parameters$count <= 400){
    if (convergence.parameters$count%%50 == 0 & config$verbose){
      print(paste("Iteration: ",convergence.parameters$count,sep=""))
    }
  }
  if (convergence.parameters$count > 400){
    if (convergence.parameters$count%%100 == 0 & config$verbose){
      print(paste("Iteration: ",convergence.parameters$count,sep=""))
    }
  }
}
