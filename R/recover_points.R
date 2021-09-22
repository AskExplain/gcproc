#' Recover missing data points via imputation or prediction
#'
#' Main function to recover missing points. Used internally, but can be re-purposed by the user.
#'
#' @param x Matrix of dataset x
#' @param y Matrix of dataset y
#' @param fixed Fixed parameters from gcproc
#' @param code Code parameters from gcproc
#' @param main.parameters Main parameters from gcproc
#' @param config Configuration parameters from gcproc
#' @param recover Recover list from gcproc
#'
#' @return  Recovered data from imputation or prediction, with the design matrices and any user input parameters and functions
#' @export
recover_points <- function(x,
                           y,
                           fixed,
                           code,
                           main.parameters,
                           config,
                           recover){


  for (method in recover$method){

    if (!is.null(recover$y)){

      X.x <- transform.data(x)%*%(main.parameters$u.beta)
      if (fixed$i_dim == T){
        Y.y <- transform.data(y)%*%(main.parameters$v.beta)
        Y.X <- Y.y + X.x
      } else {
        Y.X <- X.x
      }

      y[,which((colSums(recover$y)>0)==T)]  <- do.call('cbind',lapply(X = c(which((colSums(recover$y)>0)==T)), FUN = function(id_col){

        test_id.y <- as.logical(recover$y[,id_col])
        train_id.y <- as.logical(1 - recover$y[,id_col])

        if (min(y)==0){
          sparse.y <- log(y[,id_col]+1)
          to_exp <- T
        } else {
          sparse.y <- y[,id_col]
          to_exp <- F
        }

        if (any(test_id.y) & any(train_id.y)){

          covariate_predictors <- cbind(1,Y.X[train_id.y,])
          test_predictors <- cbind(1,Y.X[test_id.y,])

          if (method=="knn"){

            sparse.y[test_id.y] <-
              FNN::knn.reg(
                train = covariate_predictors,
                test = test_predictors,
                y = sparse.y[train_id.y],
                k = 5
              )$pred

          }

          if (method=="glmnet"){
            sparse.y[test_id.y] <-  c(predict(glmnet::cv.glmnet(x=(covariate_predictors),y=sparse.y[train_id.y],type.measure = "mse"),(test_predictors), s = "lambda.min"))
          }

          if (method=="matrix.projection"){
            sparse.y[test_id.y] <- ((test_predictors)%*%(MASS::ginv(t(covariate_predictors)%*%(covariate_predictors))%*%t(covariate_predictors)%*%((sparse.y[train_id.y]))))
          }

          if (!is.null(recover$fn)){
            sparse.y[test_id.y] <- recover$fn(train = covariate_predictors, test = test_predictors, y = sparse.y[train_id.y], parameters = recover$param)
          }


        }

        return(if(to_exp){exp(sparse.y)-1}else{sparse.y})
      }))

      y[y<0] <- 0
      y <- as.matrix(y)



      recover$predict.y <- y


    }

    if (!is.null(recover$x)){

      Y.y <- transform.data(y)%*%(main.parameters$v.beta)
      if (fixed$i_dim == T){
        X.x <- transform.data(x)%*%(main.parameters$u.beta)
        Y.X <- Y.y + X.x
      } else {
        Y.X <- Y.y
      }

      x[,which((colSums(recover$x)>0)==T)]  <- do.call('cbind',lapply(X = c(which((colSums(recover$x)>0)==T)),FUN = function(id_col){

        test_id.x <- as.logical(recover$x[,id_col])
        train_id.x <- as.logical(1 - recover$x[,id_col])

        if (min(x)==0){
          sparse.x <- log(x[,id_col]+1)
          to_exp <- T
        } else {
          sparse.x <- x[,id_col]
          to_exp <- F
        }

        if (any(test_id.x) & any(train_id.x)){

          covariate_predictors <- cbind(1,Y.X[train_id.x,])
          test_predictors <- cbind(1,Y.X[test_id.x,])

          if (method=="knn"){

            sparse.x[test_id.x] <-
              FNN::knn.reg(
                train = covariate_predictors,
                test = test_predictors,
                y = sparse.x[train_id.x],
                k = 5
              )$pred
          }
          if (method=="glmnet"){

            sparse.x[test_id.x] <- c(predict(glmnet::cv.glmnet(x=(covariate_predictors),y=sparse.x[train_id.x],type.measure = "mse"),(test_predictors), s = "lambda.min"))
          }

          if (method=="matrix.projection"){

            sparse.x[test_id.x] <- ((test_predictors)%*%(MASS::ginv(t(covariate_predictors)%*%(covariate_predictors))%*%t(covariate_predictors)%*%(sparse.x[train_id.x])))
          }

          if (!is.null(recover$fn)){

            sparse.x[test_id.x] <- recover$fn(train = covariate_predictors, test = test_predictors, y = sparse.x[train_id.x], parameters = recover$param)
          }


        }
        return(if(to_exp){exp(sparse.x)-1}else{sparse.x})
      }))

      x[x<0] <- 0
      x <- as.matrix(x)

      recover$predict.x <- x

    }


  }


  return(recover)
}


