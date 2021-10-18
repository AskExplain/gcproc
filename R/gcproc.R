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

  set.seed(config$seed)

  convergence.parameters <- list()

  recover$predict.list <- lapply(c(1:length(data_list)),function(X){NULL})

  pivots <- list(
    alpha = chunk(sample(c(1:config$i_dim)),config$n_batch),
    beta = chunk(sample(c(1:config$j_dim)),config$n_batch)
  )

  batch_table <- cbind(rep(1:length(pivots$alpha),each=length(pivots$beta)),
                       rep(1:length(pivots$beta),times=length(pivots$alpha))
  )

  internal_pivots <- list(
    alpha=pivots$alpha[[1]],
    beta=pivots$beta[[1]]
  )


  initialise.model <- initialise.gcproc(data_list = data_list,
                                        config = config,
                                        covariate = covariate,
                                        transfer = transfer,
                                        join = join,
                                        pivots = internal_pivots)


  main.parameters <- initialise.model$main.parameters
  main.code <- initialise.model$main.code
  main.proportion <- initialise.model$main.proportion
  main.index <- initialise.model$main.index

  names(main.code$code) <- if(is.null(transfer$code)){unique(do.call('c',lapply(c(1:length(covariate$factor)),function(X){c(unique(colnames(covariate$factor[[X]])))})))}else{names(main.code$code)}
  names(main.parameters$alpha) <- unique(join$alpha)
  names(main.parameters$beta) <- unique(join$beta)

  if (config$verbose){
    print(paste("Beginning gcproc learning with:    Sample dimension reduction (config$i_dim): ",config$i_dim, "    Feature dimension reduction (config$j_dim): ", config$j_dim,"    Tolerance Threshold: ", config$tol, "   Maximum number of iterations: ", config$max_iter, "   Verbose: ", config$verbose, sep=""))
  }


  for (set.of.batch.id in c(0:(config$n_batch-1))){

    print(paste("Batch number :   ",set.of.batch.id,sep=""))

    mini.batch_table <- batch_table[sample(seq(c(1+config$n_batch*set.of.batch.id),(config$n_batch*(1+set.of.batch.id)),1)),]

    main_batches <-
      # parallel::mc
      lapply(X = c(1:dim(mini.batch_table)[1]),function(batch){
        pivots <- list(alpha = pivots$alpha[[mini.batch_table[batch,1]]],
                       beta = pivots$beta[[mini.batch_table[batch,2]]])

        for (iter.update.outer in c(1:config$n_epochs)){
          for (i in sample(1:length(data_list))){
            for (iter.update.inner in c(1:config$n_epochs)){

              internal.param <- list(
                alpha = main.parameters$alpha[[join$alpha[i]]],
                beta = main.parameters$beta[[join$beta[i]]]
              )

              return_update <- update_set(x = as.matrix(data_list[[i]]),
                                          main.parameters = internal.param,
                                          main.code = main.code,
                                          main.proportion = main.proportion[[i]],
                                          main.index = main.index[[i]],
                                          pivots = pivots,
                                          fix = transfer$fix)

              main.parameters$alpha[[join$alpha[i]]] <- return_update$main.parameters$alpha
              main.parameters$beta[[join$beta[i]]] <- return_update$main.parameters$beta

              main.code <- return_update$main.code
              main.proportion[[i]] <- return_update$main.proportion

              if (!covariate$fix){
                main.index[[i]] <- t(apply(return_update$main.proportion,1,function(X){X==min(X)}))
              }

            }
          }
        }

        return(list(pivots = pivots,
                    main.code = main.code,
                    main.parameters = main.parameters,
                    main.proportion = main.proportion,
                    main.index = main.index
        )
        )

      })
    # ,mc.cores = config$n_cores)


    for (X in (1:length(data_list))){
      for (Y in (1:length(main_batches))){

        main.proportion[[X]] <- main.proportion[[X]] + main_batches[[Y]]$main.proportion[[X]]
        if (!covariate$fix){
          main.index[[X]] <- t(apply(main.proportion[[X]],1,function(X){X==min(X)}))
        }
      }
    }

    for (batch.id in 1:length(main_batches)){

      for (join.id in c(1:length(data_list))){
        main.parameters$alpha[[join$alpha[join.id]]][main_batches[[batch.id]]$pivots$alpha,] <- main_batches[[batch.id]]$main.parameters$alpha[[join$alpha[join.id]]][main_batches[[batch.id]]$pivots$alpha,]
        main.parameters$beta[[join$beta[join.id]]][,main_batches[[batch.id]]$pivots$beta] <- main_batches[[batch.id]]$main.parameters$beta[[join$beta[join.id]]][,main_batches[[batch.id]]$pivots$beta]
      }

      main.code$encode[main_batches[[batch.id]]$pivots$alpha,main_batches[[batch.id]]$pivots$beta] <- main_batches[[batch.id]]$main.code$encode[main_batches[[batch.id]]$pivots$alpha,main_batches[[batch.id]]$pivots$beta]

      for (code.id in c(1:length(main.code$code))){
        main.code$code[[code.id]][main_batches[[batch.id]]$pivots$alpha,main_batches[[batch.id]]$pivots$beta] <- main_batches[[batch.id]]$main.code$code[[code.id]][main_batches[[batch.id]]$pivots$alpha,main_batches[[batch.id]]$pivots$beta]
      }


    }

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


  if (config$verbose){
    print("Learning has converged for gcproc, beginning prediction (if requested) and dimension reduction")
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

  main.labels <- list(probabilties = main.proportion,
                      labels = lapply(c(1:length(data_list)),function(X){apply(main.proportion[[X]],1,function(X){which(X==min(X))})})
  )

  runtime.end <- Sys.time()


  if (config$verbose){
    print(paste("Done! Total runtime of   ", runtime.end - runtime.start ,sep=""))
  }

  return(list(

    main.parameters = main.parameters,

    main.code = main.code,

    main.labels = main.labels,

    covariate = covariate,

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
                       main.code,
                       main.proportion,
                       main.index,
                       pivots,
                       fix){

  internal.code <- Reduce('+',lapply(main.code$code,function(X){X[pivots$alpha,pivots$beta]}))

  for (X in sample(c(1:dim(main.proportion)[2]))){
    main.parameters$alpha[pivots$alpha,which((main.index[,X]==1))] <- (t((x[which((main.index[,X]==1)),])%*%t((internal.code)%*%t(main.parameters$beta[,pivots$beta]))%*%pinv(t((internal.code)%*%t(main.parameters$beta[,pivots$beta])))))
    main.parameters$beta[,pivots$beta] <- (t(pinv(((t(main.parameters$alpha[pivots$alpha,which((main.index[,X]==1))])%*%(internal.code))))%*%t(t(main.parameters$alpha[pivots$alpha,which((main.index[,X]==1))])%*%(internal.code))%*%(x[which((main.index[,X]==1)),])))
  }



  for (X in 1:dim(main.proportion)[2]){

    main.code$encode[pivots$alpha,pivots$beta] <- (main.parameters$alpha[pivots$alpha,which((main.index[,X]==1))]%*%(x[which((main.index[,X]==1)),])%*%(main.parameters$beta[,pivots$beta]))

    if(!fix){
      main.code$code[[X]][pivots$alpha,pivots$beta] <- pinv(t(main.parameters$alpha[pivots$alpha,which(main.index[,X]==1)]))%*%(main.code$encode[pivots$alpha,pivots$beta])%*%pinv(main.parameters$beta[,pivots$beta])
    }

  }


  for (X in c(1:dim(main.proportion)[2])){
    x.alpha.code <- x%*%(main.parameters$beta[,pivots$beta])%*%MASS::ginv(t(main.parameters$beta[,pivots$beta])%*%(main.parameters$beta[,pivots$beta]))
    alpha.code <- t(main.parameters$alpha[pivots$alpha,])%*%main.code$code[[X]][pivots$alpha,pivots$beta]
    project.x <- MASS::ginv(t(alpha.code)%*%alpha.code)%*%t(alpha.code)%*%x.alpha.code

    main.proportion[,X] <- rowMeans(abs(x.alpha.code - alpha.code))
  }

  return(list(main.parameters = main.parameters,
              main.code = main.code,
              main.proportion = main.proportion
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
