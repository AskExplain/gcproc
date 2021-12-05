#' Multiple Imputation via Generalised Canonical Procrustes
#'
#'
#' @param data List of data matrices of varying dimensionality. Attempts to find similarities among all datasets with a core structure.
#' @param config Configuration parameters (required, default provided)
#' @param transfer Transferring pre-trained model parameters (not required)
#' @param recover Important information used for prediction or imputation (not required)
#' @param join Join parameters that share the same axes to be similar (not required)
#' 
#' @export
multiple_imputation <- function(data_list,
                                multi_impute = 5,
                                config = gcproc::extract_config(verbose = F),
                                transfer = gcproc::extract_transfer_framework(verbose = F),
                                recover = gcproc::extract_recovery_framework(verbose = F),
                                join = gcproc::extract_join_framework(verbose=F)
){
  
  
  main_predict.list <- list()
  for (seed_id in c(1:multi_impute)){
    config$seed <- seed_id
    config$verbose <- F
    gcproc.model <- gcproc::gcproc(data = data_list,
                                   config = config,
                                   transfer = transfer,
                                   recover = recover,
                                   join = join)
    main_predict.list <- c(main_predict.list,list(gcproc.model$recover$predict.list))
  }
  
  
  
  # final_scores <- lapply(c(1:length(recover$design.list)),function(Y){
  #   
  #   if (!is.null(recover$design.list[[Y]])){
  #     internal_predict.list <- lapply(c(1:length(main_predict.list)),function(X){
  #       data_list[[Y]]*(recover$design.list[[Y]]==0) + main_predict.list[[X]][[Y]]*(recover$design.list[[Y]]==1)
  #     })
  # 
  #     prediction <- array(Reduce("+",internal_predict.list) / multi_impute ,dim=dim(recover$design.list[[Y]]))
  #     
  #     standard.error <- array(sqrt(Reduce("+",lapply(c(1:length(main_predict.list)),function(X){
  #       ((prediction - (data_list[[Y]]*(recover$design.list[[Y]]==0) + main_predict.list[[X]][[Y]]*(recover$design.list[[Y]]==1)))^2)
  #     }))/multi_impute),dim=dim(recover$design.list[[Y]]))
  #     
  #     return(list(pred = prediction, pred.se = standard.error ))
  #   } else{
  #     return(list(pred = NULL, pred.se = NULL ))
  #   }
  #   
  # })
  # 
  return(main_predict.list)
  
}

