# branch.gcproc <- function(
#                          y,
#                          x,
#                          config = NULL,
#                          branch_size = 3,
#                          anchors = NULL,
#                          pivots = NULL
# ){
#
#
#   if (is.null(anchors)){
#     method = "cv"
#   } else {
#     method = "transfer"
#   }
#
#   # Prepare
#   branch_id <- chunk(c(1:dim(x)[1]),branch_size)
#
#   gcproc.rotation_parameter <- list()
#
#   # Action
#   for (i in c(1:length(branch_id))){
#
#     if (method == "transfer"){
#       internal_gcproc.model <-  transfer.gcproc(
#         y = y,
#         x = x[branch_id[[i]],],
#         config = config,
#         anchors = anchors,
#         pivots = pivots
#       )
#     }
#     if (method == "cv"){
#       internal_gcproc.model <-  cv.gcproc(
#         y = y,
#         x =  x[branch_id[[i]],],
#         config = config,
#         pivots = pivots
#       )
#     }
#
#     gcproc.rotation_parameter[[i]] <- internal_gcproc.model$main.parameters$u.beta
#   }
#
#   # Finalise
#   final_rotation <- Reduce('+',gcproc.rotation_parameter)
#
#   # Return
#   return(
#     list(final_rotation = final_rotation)
#   )
# }
