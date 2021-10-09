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

  set.seed(config$config$seed)

  initialise = TRUE

  # Prepare convergence checking parameters
  count = 0
  update = 0
  bootstrap_predict = 0
  score.vec <- c()

  convergence.parameters <- list(count = count, score.vec = score.vec)

  score_lag <- 2 # How many previous scores kept track of
  accept_score <- 1 # How many scores used to calculate previous and current "mean score"

  pivots <- list()
  pivots$alpha <- c(1:config$expand.dim)
  pivots$beta <- c(1:config$expand.dim)

  recover$predict.list <- lapply(c(1:length(data_list)),function(X){NULL})

  initialise.model <- initialise.gcproc(data_list = data_list,
                                        config = config,
                                        transfer = transfer,
                                        pivots = pivots)

  main.parameters <- initialise.model$main.parameters
  code <- initialise.model$code

  if (config$verbose){
    print(paste("Beginning gcproc learning with:    Sample dimension reduction (config$i_dim): ",config$i_dim, "    Feature dimension reduction (config$j_dim): ", config$j_dim,"    Tolerance Threshold: ", config$tol, "   Maximum number of iterations: ", config$max_iter, "   Verbose: ", config$verbose, sep=""))
  }

  while (T){



    if ("decode"%in%recover$method){

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




    prev_code <- code

    for (i in 1:length(data_list)){


      return_update <- update_set(x = as.matrix(data_list[[i]]),
                                  main.parameters = main.parameters[[i]],
                                  code = code,
                                  pivots = pivots)


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



    matrix.residuals <- prev_code$encode - code$encode

    total.mae <- mean(abs(matrix.residuals))

    # Check convergence
    convergence.parameters$score.vec <- c(convergence.parameters$score.vec, total.mae)
    convergence.parameters$MAE <- mean(tail(convergence.parameters$score.vec,accept_score))
    convergence.parameters$prev.MAE <- mean(tail(convergence.parameters$score.vec,score_lag)[1:accept_score])

    if ( convergence.parameters$count > ( score_lag ) ){
      if (config$verbose == T){
        print(paste("Iteration: ",convergence.parameters$count," with Tolerance of: ", abs(convergence.parameters$prev.MAE - convergence.parameters$MAE),sep=""))
      }
    } else {
      if (config$verbose){
        print(paste("Iteration: ",convergence.parameters$count," ... initialising ... ",sep=""))
      }
    }

    if (convergence.parameters$count > config$min_iter &  convergence.parameters$count > ( score_lag ) ){
      if (abs(convergence.parameters$prev.MAE - convergence.parameters$MAE) < config$tol){
        if ((convergence.parameters$count > config$max_iter ) | ( length(pivots$alpha) == config$i_dim & length(pivots$beta) == config$j_dim)){
          if (update == 0){
            main_pivot <- pivots
          }
          update = update + 1
          pivots$alpha <- sample(main_pivot$alpha,config$expand.dim)
          pivots$beta <- sample(main_pivot$beta,config$expand.dim)



          if (update > config$update){


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

              bootstrap_predict <- bootstrap_predict + 1
              if (bootstrap_predict > config$bootstrap){
                break
              }

            } else {
              break
            }

          }
        }
        if (update ==  0 & length(pivots$alpha) != config$i_dim){
          for (i in 1:config$expand.dim){
          pivots$alpha <- c(pivots$alpha,tail(pivots$alpha,1)+1)
          main.parameters <- lapply(c(1:length(data_list)),function(X){
            internal_alpha <- main.parameters[[X]]$alpha
            internal_alpha[tail(pivots$alpha,1),] <- rnorm(dim(data_list[[X]])[1])

            return(list(alpha = internal_alpha,
                        beta = main.parameters[[X]]$beta))
          })
          code$encode[tail(pivots$alpha,1),1:tail(pivots$beta,1)] <- code$code[tail(pivots$alpha,1),1:tail(pivots$beta,1)] <- rnorm(tail(pivots$beta,1))
          }
        }
        if (update == 0 & length(pivots$beta) != config$j_dim){
          for (i in 1:config$expand.dim){
          pivots$beta <- c(pivots$beta,tail(pivots$beta,1)+1)
          main.parameters <- lapply(c(1:length(data_list)),function(X){
            internal_beta <- main.parameters[[X]]$beta
            internal_beta[,tail(pivots$beta,1)] <- rnorm(dim(data_list[[X]])[2])

            return(list(alpha = main.parameters[[X]]$alpha,
                        beta = internal_beta))
          })
          code$encode[1:tail(pivots$alpha,1),tail(pivots$beta,1)] <- code$code[1:tail(pivots$alpha,1),tail(pivots$beta,1)] <- rnorm(tail(pivots$alpha,1))
          }
        }
      }
    }



    convergence.parameters$count = convergence.parameters$count + 1

  }



  if (config$verbose){
    print("Learning has converged for gcproc, beginning prediction (if requested) and dimension reduction")
  }









  dimension_reduction <- lapply(c(1:length(data_list)),function(Y){

    x <- as.matrix(data_list[[Y]])

    feature_x.dim_reduce.encode <- t(main.parameters[[Y]]$alpha%*%x)
    sample_x.dim_reduce.encode <- x%*%main.parameters[[Y]]$beta

    feature_x.dim_reduce.code <- t(pinv(t(main.parameters[[Y]]$alpha))%*%main.parameters[[Y]]$alpha%*%x)
    sample_x.dim_reduce.code <- x%*%main.parameters[[Y]]$beta%*%pinv((main.parameters[[Y]]$beta))

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

    code = code,

    recover =  recover,

    dimension_reduction = dimension_reduction,

    meta.parameters = list(
      config = config,
      join = join,
      pivots = pivots,
      update = update,
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
                       code,
                       pivots){

  main.parameters$alpha[tail(pivots$alpha,config$expand.dim),] <- (t(x%*%t((code$code[tail(pivots$alpha,config$expand.dim),tail(pivots$beta,config$expand.dim)])%*%t(main.parameters$beta[,tail(pivots$beta,config$expand.dim)]))%*%pinv(t((code$code[tail(pivots$alpha,config$expand.dim),tail(pivots$beta,config$expand.dim)])%*%t(main.parameters$beta[,tail(pivots$beta,config$expand.dim)])))))
  main.parameters$beta[,tail(pivots$beta,config$expand.dim)] <- (t(pinv(((t(main.parameters$alpha[tail(pivots$alpha,config$expand.dim),])%*%(code$code[tail(pivots$alpha,config$expand.dim),tail(pivots$beta,config$expand.dim)]))))%*%t(t(main.parameters$alpha[tail(pivots$alpha,config$expand.dim),])%*%(code$code[tail(pivots$alpha,config$expand.dim),tail(pivots$beta,config$expand.dim)]))%*%x))

  code$encode[tail(pivots$alpha,config$expand.dim),tail(pivots$beta,config$expand.dim)] <- (main.parameters$alpha[tail(pivots$alpha,config$expand.dim),]%*%( x )%*%(main.parameters$beta[,tail(pivots$beta,config$expand.dim)]))
  code$code[tail(pivots$alpha,config$expand.dim),tail(pivots$beta,config$expand.dim)] = MASS::ginv(main.parameters$alpha[tail(pivots$alpha,config$expand.dim),]%*%t(main.parameters$alpha[tail(pivots$alpha,config$expand.dim),]))%*%(code$encode[tail(pivots$alpha,config$expand.dim),tail(pivots$beta,config$expand.dim)])%*%MASS::ginv(t(main.parameters$beta[,tail(pivots$beta,config$expand.dim)])%*%main.parameters$beta[,tail(pivots$beta,config$expand.dim)])

  #
#   sample_beta.pivot <- sample(pivots$beta,config$expand.dim)
#   sample_alpha.pivot <- sample(pivots$alpha,config$expand.dim)
#
#   main.parameters$alpha[sample_alpha.pivot,] <- (t(x%*%t((code$code[sample_alpha.pivot,sample_beta.pivot])%*%t(main.parameters$beta[,sample_beta.pivot]))%*%pinv(t((code$code[sample_alpha.pivot,sample_beta.pivot])%*%t(main.parameters$beta[,sample_beta.pivot])))))
#   main.parameters$beta[,sample_beta.pivot] <- (t(pinv(((t(main.parameters$alpha[sample_alpha.pivot,])%*%(code$code[sample_alpha.pivot,sample_beta.pivot]))))%*%t(t(main.parameters$alpha[sample_alpha.pivot,])%*%(code$code[sample_alpha.pivot,sample_beta.pivot]))%*%x))
#
#   code$encode[sample_alpha.pivot,sample_beta.pivot] <- (main.parameters$alpha[sample_alpha.pivot,]%*%( x )%*%(main.parameters$beta[,sample_beta.pivot]))
#   code$code[sample_alpha.pivot,sample_beta.pivot] = MASS::ginv(main.parameters$alpha[sample_alpha.pivot,]%*%t(main.parameters$alpha[sample_alpha.pivot,]))%*%(code$encode[sample_alpha.pivot,sample_beta.pivot])%*%MASS::ginv(t(main.parameters$beta[,sample_beta.pivot])%*%main.parameters$beta[,sample_beta.pivot])


  return(list(main.parameters = main.parameters,
              code = code
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



# vi_update_set <- function(x,
#                           config,
#                           main.parameters,
#                           code,
#                           transfer,
#                           convergence.parameters
# ){
#
#   config$eta <- if(is.null(config$eta)){1e-2}else{config$eta}
#
#   b.a <- config$eta
#   a.b <- 1-b.a
#
#   internal.score <- c()
#
#   set.seed(config$seed+convergence.parameters$count)
#
#
#   config$dim_batches <- if(is.null(config$dim_batches)){5}else{config$dim_batches}
#
#   i.sample <- chunk(sample(c(1:config$i_dim)),config$dim_batches)
#   j.sample <- chunk(sample(c(1:config$j_dim)),config$dim_batches)
#
#   dim_batch_table <- c(1:config$dim_batches)
#   dim_batch_table <- sapply(c(1:2),function(X){sample(dim_batch_table)})
#
#   internal_list <- list()
#
#   if (T){
#
#     to_return <- parallel::mclapply(c(1:dim(dim_batch_table)[1]),function(i){
#       i.ids <- i.sample[[dim_batch_table[i,1]]]
#       j.ids <- j.sample[[dim_batch_table[i,2]]]
#
#       alpha <- main.parameters$alpha[i.ids,,drop=F]
#       beta <- main.parameters$beta[,j.ids,drop=F]
#
#       x.data <- (x)
#
#       x_code <- code$code[i.ids,j.ids]
#
#       alpha.param <- t(x.data%*%t((x_code)%*%t(beta))%*%pinv(t((x_code)%*%t(beta))))
#       beta.param <- t(pinv(((t(alpha.param)%*%(x_code))))%*%t(t(alpha.param)%*%(x_code))%*%x.data)
#
#       encode.param = MASS::ginv(alpha.param%*%t(alpha.param))%*%(x_encode)%*%MASS::ginv(t(beta.param)%*%beta.param)
#       code.param = MASS::ginv(alpha.param%*%t(alpha.param))%*%(x_encode)%*%MASS::ginv(t(beta.param)%*%beta.param)
#
#       internal_list$alpha = alpha.param
#       internal_list$beta = beta.param
#       internal_list$code = code.param
#
#       internal_list$i.ids = i.ids
#       internal_list$j.ids = j.ids
#
#       return(internal_list)
#     },mc.silent = config$verbose,mc.cores = config$n_cores)
#
#
#     for (i in 1:dim(dim_batch_table)[1]){
#
#       if (is.character(to_return[[i]])){
#         to_return[[i]] <- internal_list
#       }
#
#       alpha <- to_return[[i]]$alpha
#       beta <- to_return[[i]]$beta
#
#       i.ids <- to_return[[i]]$i.ids
#       j.ids <- to_return[[i]]$j.ids
#
#       main.parameters$alpha[i.ids,] <- a.b*main.parameters$alpha[i.ids,] + b.a*alpha
#       main.parameters$beta[,j.ids] <- a.b*main.parameters$beta[,j.ids] + b.a*beta
#
#       code$code[i.ids,j.ids] <- a.b*code$code[i.ids,j.ids] + b.a*to_return[[i]]$code
#
#     }
#
#
#   }
#
#
#   code$encode <- (main.parameters$alpha%*%( x )%*%(main.parameters$beta))
#
#
#
#
#
#
#
#
#
#   final_vi.update <- list(
#     main.parameters = main.parameters,
#     code = code
#   )
#
#   return(final_vi.update)
#
# }
