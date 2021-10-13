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


  pivots <- list()
  pivots$alpha <- c(1:min(c(config$i_dim,config$extend_dim)))
  pivots$beta <- c(1:min(c(config$j_dim,config$extend_dim)))

  initialise.model <- initialise.gcproc(data_list = data_list,
                                        config = config,
                                        covariate = covariate,
                                        transfer = transfer,
                                        join = join,
                                        pivots = pivots)


  main.parameters <- initialise.model$main.parameters
  main.code <- initialise.model$main.code
  main.index <- initialise.model$main.index

  if (is.null(transfer$code)){
    # Update full dimensions iteratively in one pass
    if (config$verbose){
      print(paste("Construct parameter approximation"))
    }

    gcproc.update <- run_gcproc_single_pass(data_list = data_list,
                                            main.parameters = main.parameters,
                                            main.code = main.code,
                                            config = config,
                                            covariate = covariate,
                                            join = join,
                                            transfer = transfer,
                                            pivots = pivots,
                                            main.index = main.index)

    main.code <- gcproc.update$main.code
    main.parameters <- gcproc.update$main.parameters
  }


  if (config$verbose){
    print(paste("Beginning gcproc learning with:    Sample dimension reduction (config$i_dim): ",config$i_dim, "    Feature dimension reduction (config$j_dim): ", config$j_dim,"    Tolerance Threshold: ", config$tol, "   Maximum number of iterations: ", config$max_iter, "   Verbose: ", config$verbose, sep=""))
  }




  while (T){



    prev_code = main.code




    names(main.code$code) <- if(is.null(transfer$code)){do.call('c',lapply(c(1:dim(covariate$factor)[2]),function(X){c(unique(covariate$factor[,X]))}))}else{names(main.code$code)}
    names(main.parameters$alpha) <- unique(join$alpha)
    names(main.parameters$beta) <- unique(join$beta)



    set.seed(main_sample_seed)

    anchor <- list(alpha = sample(1:config$i_dim,min(config$extend_dim,config$i_dim)),
                   beta = sample(1:config$j_dim,min(config$extend_dim,config$j_dim)))

    pivots_i_dim.list <- chunk(sample(c(1:config$i_dim)[anchor$alpha]),length(anchor$alpha)/min(c(config$extend_dim,config$i_dim)))
    pivots_j_dim.list <- chunk(sample(c(1:config$j_dim)[anchor$beta]),length(anchor$beta)/min(c(config$extend_dim,config$j_dim)))

    pivots <- list(alpha = pivots_i_dim.list,
                   beta = pivots_j_dim.list)



    for (j in 1:length(data_list)){

      index <- main.index[[j]]

      internal.param <- list(
        alpha = main.parameters$alpha[[join$alpha[j]]],
        beta = main.parameters$beta[[join$beta[j]]]
      )

      gcproc.update <- run_gcproc_parallel(x = as.matrix(data_list[[j]]),
                                           main.parameters = internal.param,
                                           main.code = main.code,
                                           config = config,
                                           fix = transfer$fix,
                                           join = join,
                                           pivots = pivots,
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

        internal_anchor.update <- internal_anchor.update + 1

        if (internal_anchor.update > config$n_update){

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

            internal_bootstrap <- internal_bootstrap + 1

            if (internal_bootstrap > config$n_bootstrap){
              break
            }

          } else {
            break
          }

          internal_anchor.update <- 0

        }

        main_sample_seed <- main_sample_seed + 1

      }
    }

    convergence.parameters$count <- convergence.parameters$count + 1

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
                       pivots,
                       fix){

  internal.code <- Reduce('+',lapply(c(1:dim(index)[1]),function(X){
    index[X,2] * main.code$code[[index[X,1]]]
  }))/sum(index[,2])

  main.parameters$alpha[pivots$alpha,] <- (t(x%*%t((internal.code[pivots$alpha,pivots$beta])%*%t(main.parameters$beta[,pivots$beta]))%*%pinv(t((internal.code[pivots$alpha,pivots$beta])%*%t(main.parameters$beta[,pivots$beta])))))
  main.parameters$beta[,pivots$beta] <- (t(pinv(((t(main.parameters$alpha[pivots$alpha,])%*%(internal.code[pivots$alpha,pivots$beta]))))%*%t(t(main.parameters$alpha[pivots$alpha,])%*%(internal.code[pivots$alpha,pivots$beta]))%*%x))

  main.code$encode[pivots$alpha,pivots$beta] <- (main.parameters$alpha[pivots$alpha,]%*%( x )%*%(main.parameters$beta[,pivots$beta]))

  for (i in which(index[,2]==1)){

    if(!fix){
      main.code$code[[i]][pivots$alpha,pivots$beta] <- pinv(t(main.parameters$alpha[pivots$alpha,]))%*%(main.code$encode[pivots$alpha,pivots$beta])%*%pinv(main.parameters$beta[,pivots$beta])
    }

  }

  return(list(main.parameters = main.parameters,
              main.code = main.code
  ))

}



transform.data <- function(x,method="scale"){

  if (method == "scale"){
    center = T
    scale = T

    x <- as.matrix(x)
    nc <- ncol(x)
    if (is.logical(center)) {
      if (center) {
        center <- colMeans(x, na.rm=TRUE)
        x <- sweep(x, 2L, center, check.margin=FALSE)
      }
    }
    else if (is.numeric(center) && (length(center) == nc))
      x <- sweep(x, 2L, center, check.margin=FALSE)
    else
      stop("length of 'center' must equal the number of columns of 'x'")
    if (is.logical(scale)) {
      if (scale) {
        f <- function(v) {
          v <- v[!is.na(v)]
          sqrt(sum(v^2) / max(1, length(v) - 1L))
        }
        scale <- apply(x, 2L, f)
        scale <- sapply(scale,function(scale){if(scale==0|is.na(scale)){1}else{scale}})
        x <- sweep(x, 2L, scale, "/", check.margin=FALSE)
      }
    }
    else if (is.numeric(scale) && length(scale) == nc)
      x <- sweep(x, 2L, scale, "/", check.margin=FALSE)
    else
      stop("length of 'scale' must equal the number of columns of 'x'")
    if(is.numeric(center)) attr(x, "scaled:center") <- center
    if(is.numeric(scale)) attr(x, "scaled:scale") <- scale
    x
  }

  return(x)
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


run_gcproc_parallel <- function(x,
                                main.parameters,
                                main.code,
                                config,
                                join,
                                fix,
                                pivots,
                                index){

  batch_table <- cbind(rep(1:length(pivots$alpha),length(pivots$beta)),
                       rep(1:length(pivots$beta),length(pivots$alpha))
  )

  main_batches <-
    parallel::mclapply(X = c(1:dim(batch_table)[1]),function(batch){
      pivots <- list(alpha = pivots$alpha[[batch_table[batch,1]]],
                     beta = pivots$beta[[batch_table[batch,2]]])

      for (i in 1:3){

        return_update <- update_set(x = x,
                                    main.parameters = main.parameters,
                                    main.code = main.code,
                                    pivots = pivots,
                                    fix = fix,
                                    index = index)

        main.parameters <- return_update$main.parameters
        main.code <- return_update$main.code

      }



      return(list(pivots = pivots,
                  main.code = main.code,
                  main.parameters = main.parameters))

    },mc.cores = config$n_cores)

  for (i in 1:length(main_batches)){

    main.parameters$alpha[main_batches[[i]]$pivots$alpha,] <- main_batches[[i]]$main.parameters$alpha[main_batches[[i]]$pivots$alpha,]
    main.parameters$beta[,main_batches[[i]]$pivots$beta] <- main_batches[[i]]$main.parameters$beta[,main_batches[[i]]$pivots$beta]

    main.code$encode[main_batches[[i]]$pivots$alpha,main_batches[[i]]$pivots$beta] <- main_batches[[i]]$main.code$encode[main_batches[[i]]$pivots$alpha,main_batches[[i]]$pivots$beta]

    for (j in 1:length(main_batches[[i]]$main.code$code)){
      main.code$code[[j]][main_batches[[i]]$pivots$alpha,main_batches[[i]]$pivots$beta] <- main_batches[[i]]$main.code$code[[j]][main_batches[[i]]$pivots$alpha,main_batches[[i]]$pivots$beta]

    }

  }



  return(list(main.parameters = main.parameters,
              main.code = main.code))

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









run_gcproc_single_pass <- function(data_list,
                                   main.parameters,
                                   main.code,
                                   config,
                                   covariate,
                                   join,
                                   transfer,
                                   pivots,
                                   main.index){

  set.seed(config$seed)

  # Prepare convergence checking parameters
  count = 0
  update = 0
  score.vec <- c()

  score_lag <- 2 # How many previous scores kept track of
  accept_score <- 1 # How many scores used to calculate previous and current "mean score"



  sub.pivots <- pivots

  while (T){


    names(main.code$code) <- if(is.null(transfer$code)){do.call('c',lapply(c(1:dim(covariate$factor)[2]),function(X){c(unique(covariate$factor[,X]))}))}else{names(main.code$code)}
    names(main.parameters$alpha) <- unique(join$alpha)
    names(main.parameters$beta) <- unique(join$beta)


    prev_code <- main.code

    for (i in 1:length(data_list)){

      index <- main.index[[i]]

      internal.param <- list(
        alpha = main.parameters$alpha[[join$alpha[i]]],
        beta = main.parameters$beta[[join$beta[i]]]
      )

      return_update <- update_set(x = as.matrix(data_list[[i]]),
                                  main.parameters = internal.param,
                                  main.code = main.code,
                                  index = index,
                                  pivots = sub.pivots,
                                  fix = transfer$fix)

      main.parameters$alpha[[join$alpha[i]]] <- return_update$main.parameters$alpha
      main.parameters$beta[[join$beta[i]]] <- return_update$main.parameters$beta

      main.code <- return_update$main.code

    }

    matrix.residuals <- prev_code$encode - main.code$encode

    total.mae <- mean(abs(matrix.residuals))

    # Check convergence
    score.vec <- c(score.vec, total.mae)
    MAE <- mean(tail(score.vec,accept_score))
    prev.MAE <- mean(tail(score.vec,score_lag)[1:accept_score])

    if (count > config$min_iter &  count > ( score_lag ) ){
      if (count > config$max_iter | abs(prev.MAE - MAE) < config$tol){

        if ( length(pivots$alpha) == config$i_dim & length(pivots$beta) == config$j_dim){
          break
        }



        if (update ==  0 & length(pivots$alpha) != config$i_dim){
          for (i in 1:config$extend_dim){
            pivots$alpha <- c(pivots$alpha,tail(pivots$alpha,1)+1)
            main.parameters$alpha <- lapply(c(1:length(main.parameters$alpha)),function(X){
              internal_alpha <- main.parameters$alpha[[X]]
              internal_alpha[tail(pivots$alpha,1),] <- rnorm(dim(main.parameters$alpha[[X]])[2])

              return(internal_alpha)
            })
            main.code$encode[tail(pivots$alpha,1),1:tail(pivots$beta,1)] <- rnorm(tail(pivots$beta,1))


            for (j in 1:length(main.code$code)){
              main.code$code[[j]][tail(pivots$alpha,1),1:tail(pivots$beta,1)] <- rnorm(tail(pivots$beta,1))
            }

          }
        }
        if (update == 0 & length(pivots$beta) != config$j_dim){
          for (i in 1:config$extend_dim){
            pivots$beta <- c(pivots$beta,tail(pivots$beta,1)+1)
            main.parameters$beta <- lapply(c(1:length(main.parameters$beta)),function(X){
              internal_beta <- main.parameters$beta[[X]]
              internal_beta[,tail(pivots$beta,1)] <- rnorm(dim(main.parameters$beta[[X]])[1])

              return(internal_beta)
            })

            main.code$encode[1:tail(pivots$alpha,1),tail(pivots$beta,1)] <- rnorm(tail(pivots$alpha,1))


            for (j in 1:length(main.code$code)){
              main.code$code[[j]][1:tail(pivots$alpha,1),tail(pivots$beta,1)] <- rnorm(tail(pivots$alpha,1))
            }


          }
        }

        sub.pivots <- list(
          alpha = tail(pivots$alpha,config$extend_dim),
          beta = tail(pivots$beta,config$extend_dim)
        )

      }
    }

    print(paste("Approximation iteration of:   ",count,sep=""))

    count = count + 1

  }



  return(list(main.parameters = main.parameters,
              main.code = main.code))

}
