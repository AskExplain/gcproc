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
  
  for (set.id in c(1:length(data_list))){
    
    internal.config <- config
    internal.config$verbose <- F
    
    gcproc.model_list[[set.id]] <- gcproc(data_list = data_list[[set.id]],
                                          config = internal.config,
                                          transfer = transfer_list[[set.id]],
                                          recover = recover_list[[set.id]],
                                          join = join_list[[set.id]])
    
    transfer_list[[set.id]]$main.code <- gcproc.model_list[[set.id]]$main.code
    transfer_list[[set.id]]$main.parameters <- gcproc.model_list[[set.id]]$main.parameters
    
  }


  return(gcproc.model_list)
  
}