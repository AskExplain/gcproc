#' Recover missing data points via imputation or prediction
#'
#' Main function to recover missing points. Used internally, but can be re-purposed by the user.
#'
#' @param data_list A list of datasets (matrix or tensor etc.)
#' @param code Code parameters from gcproc
#' @param main.parameters Main parameters from gcproc
#' @param config Configuration parameters from gcproc
#' @param recover Recover list from gcproc
#'
#' @return  Recovered data from imputation or prediction, with the design matrices and any user input parameters and functions
#' @export
recover_points <- function(data_list,
                           main.code,
                           main.parameters,
                           config,
                           join,
                           recover){
  
  for (task in recover$task){
    
    if ("regression" %in% task){
      
      for (method in recover$method){
        
        
        for (i in 1:length(data_list)){
          
          if (!is.null(recover$design.list[[i]])){
            
            if ("decode" %in% method){
              
              row_with_missing_points <- which((rowSums(recover$design.list[[i]])>0)==T,arr.ind = T)
              column_with_missing_points <- which((colSums(recover$design.list[[i]])>0)==T,arr.ind = T)
              
              x <- transform.data(as.matrix(data_list[[i]]), method = recover$link_function[1])
              
              alpha.code <-  cbind(1,t(main.parameters$alpha[[join$alpha[i]]]))
              beta.code <- main.code$code%*%t(main.parameters$beta[[join$beta[i]]])
              d <- MASS::ginv(t(alpha.code)%*%(alpha.code))%*%t(alpha.code)%*%x%*%t(beta.code)%*%MASS::ginv(beta.code%*%t(beta.code))
              
              alpha.code <-  cbind(1,t(main.parameters$alpha[[join$alpha[i]]]))%*%d%*%main.code$code
              beta.code <- t(main.parameters$beta[[join$beta[i]]])
              e <- MASS::ginv(t(alpha.code)%*%(alpha.code))%*%t(alpha.code)%*%x%*%t(beta.code)%*%MASS::ginv(beta.code%*%t(beta.code))
              
              alpha.code <- cbind(1,t(main.parameters$alpha[[join$alpha[i]]]))%*%d%*%main.code$code%*%e%*%t(main.parameters$beta[[join$beta[i]]])%*%(main.parameters$beta[[join$beta[i]]])
              f <- MASS::ginv(t(alpha.code)%*%(alpha.code))%*%t(alpha.code)%*%x

              pred <- cbind(1,t(main.parameters$alpha[[join$alpha[i]]]))%*%d%*%main.code$code%*%e%*%t(main.parameters$beta[[join$beta[i]]])%*%(main.parameters$beta[[join$beta[i]]])%*%f
              
              x[row_with_missing_points,column_with_missing_points]  <- pred[row_with_missing_points,column_with_missing_points]
              
              data_list[[i]] <- recover$predict.list[[i]] <- transform.data(x, method= recover$link_function[2])
              
            }
            if ("matrix.projection" %in% method){
              
              x <- transform.data(as.matrix(data_list[[i]]), method = recover$link_function[1])
              
              
              if (is.null(recover$encoded_covariate)){
                recover$encoded_covariate <- lapply(c(1:length(data_list))[-i],function(X){
                  transformed.data <- as.matrix(data_list[[X]])%*%(main.parameters$beta[[join$beta[X]]])
                  return(transformed.data)
                })
              }
              
              
              decoded_covariate <- cbind(1,
                                         transform.data(Reduce('+',lapply(c(1:length(recover$encoded_covariate)),function(X){
                                           recover$encoded_covariate[[X]]
                                         })))              )
              
              samples_with_missing_points <- which((rowSums(recover$design.list[[i]])>0)==T)
              covariate_predictors <-  decoded_covariate[-samples_with_missing_points,]
              test_predictors <- decoded_covariate[samples_with_missing_points,]
              
              elements_with_missing_points <- which((recover$design.list[[i]]>0)[samples_with_missing_points,]==T,arr.ind = T)
              x[samples_with_missing_points,][elements_with_missing_points]  <- (((test_predictors)%*%(MASS::ginv(t(covariate_predictors)%*%(covariate_predictors))%*%t(covariate_predictors)%*%(x[-samples_with_missing_points,]))))[elements_with_missing_points]
              
              data_list[[i]] <- recover$predict.list[[i]] <- transform.data(x, method = recover$link_function[2])
            }
            
            
            if ("knn" %in% method){
              
              
              x <- as.matrix(data_list[[i]])
              
              if (is.null(recover$encoded_covariate)){
                recover$encoded_covariate <- lapply(c(1:length(data_list))[-i],function(X){
                  transformed.data <- as.matrix(data_list[[X]])%*%(main.parameters$beta[[join$beta[X]]])
                  return(transformed.data)
                })
              }
              
              
              decoded_covariate <- cbind(1,
                                         transform.data(Reduce('+',lapply(c(1:length(recover$encoded_covariate)),function(X){
                                           recover$encoded_covariate[[X]]
                                         })))              )
              
              
              samples_with_missing_points <- which((rowSums(recover$design.list[[i]])>0)==T)
              covariate_predictors <-  decoded_covariate[-samples_with_missing_points,]
              test_predictors <- decoded_covariate[samples_with_missing_points,]
              
              knn_ix <- FNN::get.knnx(
                covariate_predictors,
                test_predictors,
                k = 20
              )$nn.index
              
              pred <- (x[-samples_with_missing_points,])[knn_ix[, 1], , drop = FALSE]
              if (20 > 1) {
                for (k in seq(2, 20)) {
                  pred <- pred + (x[-samples_with_missing_points,])[knn_ix[, k], , drop = FALSE]
                }
              }
              pred <- pred / 20
              
              elements_with_missing_points <- which((recover$design.list[[i]]>0)[samples_with_missing_points,]==T,arr.ind = T)
              x[samples_with_missing_points,][elements_with_missing_points]  <- (pred)[elements_with_missing_points]
              
              data_list[[i]] <- recover$predict.list[[i]] <- x
              
            }
            
            
            if (!is.null(recover$fn)){
              
              x <- transform.data(as.matrix(data_list[[i]]), method = recover$link_function[1])
              
              if (is.null(recover$encoded_covariate)){
                recover$encoded_covariate <- lapply(c(1:length(data_list))[-i],function(X){
                  transformed.data <- as.matrix(data_list[[X]])%*%(main.parameters$beta[[join$beta[X]]])
                  return(transformed.data)
                })
              }
              
              
              decoded_covariate <- cbind(1,
                                         transform.data(Reduce('+',lapply(c(1:length(recover$encoded_covariate)),function(X){
                                           recover$encoded_covariate[[X]]
                                         })))              )
              
              
              samples_with_missing_points <- which((rowSums(recover$design.list[[i]])>0)==T)
              covariate_predictors <-  decoded_covariate[-samples_with_missing_points,]
              test_predictors <- decoded_covariate[samples_with_missing_points,]
              
              pred <- recover$fn(train = covariate_predictors, test = test_predictors, y = x[-samples_with_missing_points,], parameters = recover$param)
              
              elements_with_missing_points <- which((recover$design.list[[i]]>0)[samples_with_missing_points,]==T,arr.ind = T)
              x[samples_with_missing_points,][elements_with_missing_points]  <- (pred)[elements_with_missing_points]
              
              data_list[[i]] <- recover$predict.list[[i]] <- transform.data(x, method = recover$link_function[2])
            }
            
            
          }
          
          
        }
      }
      
    }
    
    if ("classification" %in% task){
      
      
      label.projection <- c(recover$method=="label.projection")
      
      if (label.projection){
        
        for (j in which(recover$design.list==0)){
          
          label_code <- Reduce('+',lapply(c(covariate$factor[j,]),function(X){
            main.code$code[[X]]
          }))
          
          
          
          for (i in which(recover$design.list==1)){
            
            unlabel_code <- Reduce('+',lapply(c(covariate$factor[i,]),function(X){
              main.code$code[[X]]
            }))
            
            labels <- recover$labels
            
            recover$predict.list[[j]][[i]] <- apply((unlabel.decoded_covariate)%*%t(label.decoded_covariate),1,function(X){names(sort(table(labels[order(X,decreasing = T)[1]]))[1])})
            
          }
          
        }
      }
      
      
      
    }
    
  }
  
  return(list(recover=recover,
              data_list=data_list))
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
  if (method == "logistic"){
    x <- log((1e-9+x)/((1+1e-9+1e-10)-(1e-9+x)))
  }
  if (method == "logit"){
    x <- (exp(x)/(1+exp(x))) - 1e-9
  }
  if (method == "log"){
    x <- log(x+0.1)
  }
  if (method == "exp"){
    x <- exp(x)-0.1
  }
  if (method == "identity"){
    return(x)
  }
  return(x)
}
