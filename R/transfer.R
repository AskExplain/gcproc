#' Transfer model from one to another
#'
#' Main function to transferring models 
#' 
#' @param model_list A list of gcproc models 
#' @param data_list A list of datasets (matrix or tensor etc.)
#' @param config Configuration parameters (required, default provided)
#' @param transfer Transferring pre-trained model parameters (not required)
#' @param recover Important information used for prediction or imputation (not required)
#' @param join Join parameters that share the same axes to be similar (not required)
#'
#' @return  All models to be learned or updated
#' @export
transfer <- function(model_list,
                     data_list,
                     config = gcproc::extract_config(verbose = F),
                     transfer = gcproc::extract_transfer_framework(verbose = F),
                     recover = gcproc::extract_recovery_framework(verbose = F),
                     join = gcproc::extract_join_framework(verbose=F)
){
  
  # run gcproc again, with new and old data, using old model
  
  # Options - 1 - make initialise optional, pass in parameters... (too cumbersome)
  #         - 2 - run multiple gcprocs, pass code between? (fast, scalable, reasonable - different from fixing gcproc)
  #       ***     Must ensure, it can take in features for beta (manually created features)
  #       ***     transfer beta
  # Following option 2 - how to update old model (run previous data with old code, pass in new code),
  
  convergence.parameters <- list(count=0,score.vec=c())
  
  # Initialise models
  gcproc.model_list <- list()
  
  # Set up transfer of model
  transfer_list <- list()
  recover_list <- list()
  join_list <- list()
  
  for (set.id in c(1:length(data_list))){
    
    transfer_list[[set.id]] <- gcproc::extract_transfer_framework(F)
    recover_list[[set.id]] <- gcproc::extract_recovery_framework(F)
    join_list[[set.id]] <- gcproc::extract_join_framework(F)
    
    
    if (!is.null(model_list[[set.id]])){
      
      transfer_list[[set.id]]$main.code <- model_list[[set.id]]$main.code
      transfer_list[[set.id]]$main.parameters <- model_list[[set.id]]$main.parameters

      recover_list[[set.id]] <- model_list[[set.id]]$recover
      join_list[[set.id]] <- model_list[[set.id]]$meta.parameters$join
      
    } else {
      
      transfer_list[[set.id]] <- transfer
      recover_list[[set.id]] <- recover
      join_list[[set.id]] <- join
      
    }

  }
  
  while (T){
    if (convergence.parameters$count > 0){
      prev.encode <- Reduce('+',lapply(c(1:length(model_list)),function(X){gcproc.model_list[[X]]$main.code$encode}))
    } else {
      prev.encode <- 0
    }

    for (set.id in c(1:length(data_list))){
      config$verbose <- F
      gcproc.model_list[[set.id]] <- gcproc(data_list = data_list[[set.id]],
                                            config = config,
                                            transfer = transfer_list[[set.id]],
                                            recover = recover_list[[set.id]],
                                            join = join_list[[set.id]])

      transfer_list[[set.id]]$main.code <- gcproc.model_list[[set.id]]$main.code
      transfer_list[[set.id]]$main.parameters <- gcproc.model_list[[set.id]]$main.parameters
      
    }

    config$verbose <- T
    
    mae <- mean(abs(prev.encode - Reduce('+',lapply(c(1:length(model_list)),function(X){gcproc.model_list[[X]]$main.code$encode}))))
    print(mae)
    
    # Check convergence
    convergence.parameters$score.vec <- c(convergence.parameters$score.vec, mae)
    MAE <- mean(tail(convergence.parameters$score.vec,2))
    prev.MAE <- mean(tail(convergence.parameters$score.vec,3)[1:2])

    
    if ( convergence.parameters$count > ( 3 ) ){
      if (config$verbose == T){
        print(paste("Iteration: ",convergence.parameters$count," with Tolerance of: ", abs(prev.MAE - MAE),sep=""))
      }
    } else {
      if (config$verbose){
        print(paste("Iteration: ",convergence.parameters$count," ... initialising ... ",sep=""))
      }
    }
    
    if (convergence.parameters$count > config$min_iter){
      if ((convergence.parameters$count > config$max_iter ) | abs(prev.MAE - MAE) < config$tol){
        break
      }
    }
    

    convergence.parameters$count = convergence.parameters$count + 1


  }



  return(gcproc.model_list)
  
}