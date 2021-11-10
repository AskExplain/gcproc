#' Recover missing data points via imputation or prediction
#'
#' Main function to recover missing points. Used internally, but can be re-purposed by the user.
#'
#' @param data_list List of Datasets
#' @param main.code Main code from gcproc
#' @param main.parameters Main parameters from gcproc
#' @param config Configuration parameters from gcproc
#' @param recover Recover list from gcproc
#' @param join Join parameters in gcproc
#' @param fixed Fixed parameters from gcproc
#'
#' @return  Recovered data from imputation or prediction, with the design matrices and any user input parameters and functions
#' @export
recover_points <- function(data_list,
                           main.code,
                           main.parameters,
                           config,
                           recover,
                           join,
                           fixed){
  
  
  
  for (task in recover$task){
    
    if ("regression" %in% task){
      
      for (method in recover$method){
        
        for (i in 1:length(data_list)){
          
          if (!is.null(recover$design.list[[i]])){
            
            
            y <- as.matrix(data_list[[i]])
            Y.y <- transform.data(y)
            
            if (is.null(recover$encoded_covariate)){
              recover$encoded_covariate <- lapply(c(1:length(data_list))[-i],function(X){
                transformed.data <- as.matrix(data_list[[X]])%*%(main.parameters$beta[[join$beta[X]]])
                return(transformed.data)
              })
            }
            
            
            X.x <- cbind(1,
                         transform.data(Reduce('+',lapply(c(1:length(recover$encoded_covariate)),function(X){
                           recover$encoded_covariate[[X]]
                         })))
            )
            
            
            y[,which((colSums(recover$design.list[[i]])>0)==T)]  <- do.call('cbind',parallel::mclapply(c(which((colSums(recover$design.list[[i]])>0)==T)),function(id_col){
              
              test_id <- as.logical(recover$design.list[[i]][,id_col])
              train_id <- as.logical(1 - recover$design.list[[i]][,id_col])
              sparse.y <- c(y[,id_col])
              
              if (any(test_id) & any(train_id)){
                
                x.covariate_predictors <- y.covariate_predictors <- cbind(1,Y.y[train_id,-id_col]%*%main.parameters$beta[-id_col,])
                x.test_predictors <- y.test_predictors <- cbind(1,Y.y[test_id,-id_col]%*%main.parameters$beta[-id_col,])
                
                b.a <- 0
                if (!identical(recover$x,recover$design.list[[i]]) & fixed$i_dim == T) {
                  x.covariate_predictors <- cbind(1,X.x[train_id,]%*%main.parameters$u.beta)
                  x.test_predictors <- cbind(1,X.x[test_id,]%*%main.parameters$u.beta)
                  b.a <- 0.5
                }
                a.b <- 1 - b.a
                
                
                if (method=="knn"){
                  
                  sparse.y[test_id] <-
                    a.b * FNN::knn.reg(
                      train = y.covariate_predictors,
                      test = y.test_predictors,
                      y = sparse.y[train_id],
                      k = 5
                    )$pred +
                    b.a * FNN::knn.reg(
                      train = x.covariate_predictors,
                      test = x.test_predictors,
                      y = sparse.y[train_id],
                      k = 5
                    )$pred
                  
                }
                if (method=="glmnet"){
                  
                  sparse.y[test_id] <- a.b * c(predict(glmnet::cv.glmnet(x=(y.covariate_predictors),y=sparse.y[train_id],type.measure = "mse"),(y.test_predictors), s = "lambda.min")) +
                    b.a * (predict(glmnet::cv.glmnet(x=(x.covariate_predictors),y=sparse.y[train_id],type.measure = "mse"),(x.test_predictors), s = "lambda.min"))
                  
                }
                
                if (method=="matrix.projection"){
                  
                  sparse.y[test_id] <- a.b*((y.test_predictors)%*%(MASS::ginv(t(y.covariate_predictors)%*%(y.covariate_predictors))%*%t(y.covariate_predictors)%*%((as.matrix(y[train_id,id_col]))))) +
                    b.a*((x.test_predictors)%*%(MASS::ginv(t(x.covariate_predictors)%*%(x.covariate_predictors))%*%t(x.covariate_predictors)%*%((as.matrix(y[train_id,id_col])))))
                  
                }
                if (!is.null(recover$fn)){
                  
                  sparse.y[test_id] <- a.b*recover$fn(train = y.covariate_predictors, test = y.test_predictors, y = sparse.y[train_id], parameters = recover$param) +
                    b.a*recover$fn(train = x.covariate_predictors, test = x.test_predictors, y = sparse.y[train_id], parameters = recover$param)
                  
                }
                
                
              }
              
              return(sparse.y)
            }))
            
            recover$predict.y[[i]] <- as.matrix(y)
            
            
          }
        }
      }
    }
  }
  
  return(recover)
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

