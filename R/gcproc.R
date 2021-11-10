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

  convergence.parameters <- list(count=0,score.vec=c())

  recover$predict.list <- lapply(c(1:length(data_list)),function(X){NULL})

  initialise.model <- initialise.gcproc(data_list = data_list,
                                        config = config,
                                        transfer = transfer,
                                        join = join)

  main.parameters <- initialise.model$main.parameters
  main.code <- initialise.model$main.code

  if (config$verbose){
    print(paste("Beginning gcproc learning with:    Sample dimension reduction (config$i_dim): ",config$i_dim, "    Feature dimension reduction (config$j_dim): ", config$j_dim,"    Tolerance Threshold: ", config$tol, "   Maximum number of iterations: ", config$max_iter, "   Verbose: ", config$verbose, sep=""))
  }


  while (T){
    
    prev_code <- main.code
    
    for (i in 1:length(data_list)){
      
      internal.parameters <- list(alpha=main.parameters$alpha[[join$alpha[i]]],
                              beta=main.parameters$beta[[join$beta[i]]])
      
      return_update <- update_set(x = as.matrix(data_list[[i]]),
                                  main.parameters = internal.parameters,
                                  main.code = main.code
                                  )
      
      main.parameters$alpha[[join$alpha[i]]] <- return_update$main.parameters$alpha
      main.parameters$beta[[join$beta[i]]] <- return_update$main.parameters$beta
      
      main.code <- return_update$main.code
      
    }
    
    
    
    matrix.residuals <- main.code$encode - prev_code$encode
    
    total.mse <- mean(abs(matrix.residuals))
    
    # Check convergence
    convergence.parameters$score.vec <- c(convergence.parameters$score.vec, total.mse)
    MSE <- mean(tail(convergence.parameters$score.vec,2))
    prev.MSE <- mean(tail(convergence.parameters$score.vec,3)[1:2])
    
    if ( convergence.parameters$count > ( 3 ) ){
      if (config$verbose == T){
        print(paste("Iteration: ",convergence.parameters$count," with Tolerance of: ", abs(prev.MSE - MSE),sep=""))
      }
    } else {
      if (config$verbose){
        print(paste("Iteration: ",convergence.parameters$count," ... initialising ... ",sep=""))
      }
    }
    
    if (convergence.parameters$count > config$min_iter){
      if ((convergence.parameters$count > config$max_iter ) | abs(prev.MSE - MSE) < config$tol){
        break
      }
    }
    
    convergence.parameters$count = convergence.parameters$count + 1
    
  }
  
  
  if (convergence.parameters$count > 2 & any(do.call('c',lapply(recover$design.list,function(X){!is.null(X)})))){
    
    recover_data <- recover_points(
      data_list,
      main.code = main.code,
      main.parameters = main.parameters,
      config = config,
      recover = recover,
      join = join
    )
    
    recover <- recover_data
    data_list <- recover_data$predict.list
    
  }
  
  
  if (config$verbose){
    print("Learning has converged for gcproc, beginning (if requested) dimension reduction")
  }
  
  

  dimension_reduction <- lapply(c(1:length(data_list)),function(Y){

    x <- as.matrix(data_list[[Y]])

    feature_x.dim_reduce.encode <- t(main.parameters$alpha[[join$alpha[Y]]]%*%x)
    sample_x.dim_reduce.encode <- x%*%main.parameters$beta[[join$beta[Y]]]

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

    main.code = main.code,

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

    convergence.parameters = convergence.parameters


  ))

}


update_set <- function(x,
                       main.parameters,
                       main.code){

  main.code$encode <- (main.parameters$alpha%*%(x)%*%(main.parameters$beta))
  main.code$code <- pinv(t(main.parameters$alpha))%*%(main.code$encode)%*%pinv(main.parameters$beta)
  
  main.parameters$alpha <- (t((x)%*%t((main.code$code)%*%t(main.parameters$beta))%*%pinv(t((main.code$code)%*%t(main.parameters$beta)))))
  main.parameters$beta <- (t(pinv(((t(main.parameters$alpha)%*%(main.code$code))))%*%t(t(main.parameters$alpha)%*%(main.code$code))%*%(x)))
  
  
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