branch.gcproc <- function(gcproc.model,
                         y,
                         x,
                         branch_size = 1000){
  
  # Prepare
  gcproc.model <- gcproc::transfer.gcproc(
    gcproc.model = gcproc.model,
    x = x,
    y = y
  )
  
  branch_id <- chunk(c(1:dim(x)[1]),dim(x)[1]/branch_size)
  
  gcproc.model_list <- list()
  gcproc.rotation_parameter <- list(gcproc.model$main.parameters$u.beta)
  
  # Action
  for (i in c(1:length(branch_id))){
    
    if (i < length(branch_id)){
      internal_x1 <- x[branch_id[[i]],]
      internal_x2 <- x[branch_id[[i+1]],]
    }
    if (i == length(branch_id)){
      internal_x1 <- x[branch_id[[i]],]
      internal_x2 <- x[branch_id[[1]],]
    }
    
    
    anchors = list(
      anchor_y.sample = NULL,
      anchor_y.feature = gcproc.model$main.parameters$u.beta,
      anchor_x.sample = NULL,
      anchor_x.feature = NULL
    )
    
    gcproc.model_list[[i]] <- gcproc.model <- 
      transfer.gcproc(gcproc.model = gcproc.model,
                    y = internal_x1,
                    x = internal_x2,
                    anchors = anchors)
    
    gcproc.rotation_parameter[[i]] <- gcproc.model$main.parameters$u.beta
    
  }
  
  # Finalise
  final_rotation <- Reduce('*',gcproc.rotation_parameter)
  
  # Return
  return(list(gcproc.model_list = gcproc.model_list,
              final_rotation = final_rotation
              )
  )
}