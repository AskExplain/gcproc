#' Recover missing data points via imputation or prediction
#'
#' Main function to recover missing points. Used internally, but can be re-purposed by the user.
#'
#' @param x Matrix of dataset x
#' @param y Matrix of dataset y
#' @param main.parameters Main parameters from gcproc
#' @param config Configuration parameters from gcproc
#' @param recover Recover list from gcproc
#' @param fixed Fixed parameters from gcproc
#'
#' @return  Recovered data from imputation or prediction, with the design matrices and any user input parameters and functions
#' @export
recover_points <- function(x,
                           y,
                           main.parameters,
                           config,
                           recover,
                           fixed){

  for (method in recover$method){

    if (!is.null(recover$y)){


      Y.y <- scale(y)
      X.x <- scale(x)

      y[,which((colSums(recover$y)>0)==T)]  <- do.call('cbind',parallel::mclapply(c(which((colSums(recover$y)>0)==T)),function(id_col){

        test_id <- as.logical(recover$y[,id_col])
        train_id <- as.logical(1 - recover$y[,id_col])
        sparse.y <- c(y[,id_col])

        if (any(test_id) & any(train_id)){

          x.covariate_predictors <- y.covariate_predictors <- cbind(1,Y.y[train_id,-id_col]%*%main.parameters$v.beta[-id_col,])
          x.test_predictors <- y.test_predictors <- cbind(1,Y.y[test_id,-id_col]%*%main.parameters$v.beta[-id_col,])

          b.a <- 0
          if (!identical(recover$x,recover$y) & fixed$i_dim == T) {
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

      y <- as.matrix(y)


      recover$predict.y <- y


    }

    if (!is.null(recover$x)){

      Y.y <- scale(y)
      X.x <- scale(x)


      x[,which((colSums(recover$x)>0)==T)]  <- do.call('cbind',parallel::mclapply(c(which((colSums(recover$x)>0)==T)),function(id_col){

        test_id <- as.logical(recover$x[,id_col])
        train_id <- as.logical(1 - recover$x[,id_col])
        sparse.x <- c(x[,id_col])

        if (any(test_id) & any(train_id)){

          y.covariate_predictors <- x.covariate_predictors <- cbind(1,X.x[train_id,-id_col]%*%main.parameters$u.beta[-id_col,])
          y.test_predictors <- x.test_predictors <- cbind(1,X.x[test_id,-id_col]%*%main.parameters$u.beta[-id_col,])

          b.a <- 0
          if (!identical(recover$x,recover$y) & fixed$i_dim == T) {
            y.covariate_predictors <- cbind(1,Y.y[train_id,]%*%main.parameters$v.beta)
            y.test_predictors <- cbind(1,Y.y[test_id,]%*%main.parameters$v.beta)
          }
          a.b <- 1 - b.a


          if (method=="knn"){

            sparse.x[test_id] <-
              a.b * FNN::knn.reg(
                train = x.covariate_predictors,
                test = x.test_predictors,
                y = sparse.x[train_id],
                k = 5
              )$pred +
              b.a * FNN::knn.reg(
                train = y.covariate_predictors,
                test = y.test_predictors,
                y = sparse.x[train_id],
                k = 5
              )$pred

          }
          if (method=="glmnet"){

            sparse.x[test_id] <- a.b * c(predict(glmnet::cv.glmnet(x=(x.covariate_predictors),y=sparse.x[train_id],type.measure = "mse"),(x.test_predictors), s = "lambda.min")) +
              b.a * (predict(glmnet::cv.glmnet(x=(y.covariate_predictors),y=sparse.x[train_id],type.measure = "mse"),(y.test_predictors), s = "lambda.min"))

          }

          if (method=="matrix.projection"){

            sparse.x[test_id] <- a.b*((x.test_predictors)%*%(MASS::ginv(t(x.covariate_predictors)%*%(x.covariate_predictors))%*%t(x.covariate_predictors)%*%((as.matrix(x[train_id,id_col]))))) +
              b.a*((y.test_predictors)%*%(MASS::ginv(t(y.covariate_predictors)%*%(y.covariate_predictors))%*%t(y.covariate_predictors)%*%((as.matrix(x[train_id,id_col])))))

          }

          if (!is.null(recover$fn)){

            sparse.x[test_id] <- a.b*recover$fn(train = x.covariate_predictors, test = x.test_predictors, y = sparse.x[train_id], parameters = recover$param) +
              b.a*recover$fn(train = y.covariate_predictors, test = y.test_predictors, y = sparse.x[train_id], parameters = recover$param)

          }


        }
        return(sparse.x)
      }))

      x <- as.matrix(x)

      recover$predict.x <- x

    }


  }


  return(recover)
}


