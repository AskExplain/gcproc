#' Recover missing data points via imputation or prediction
#'
#' Main function to recover missing points. Used internally, but can be re-purposed by the user.
#'
#' @param x Matrix of dataset x
#' @param y Matrix of dataset y
#' @param fixed Fixed parameters from gcproc
#' @param main.parameters Main parameters from gcproc
#' @param config Configuration parameters from gcproc
#' @param recover Recover list from gcproc
#'
#' @return  Recovered data from imputation or prediction, with the design matrices and any user input parameters and functions
#' @export
recover_points <- function(x,
                           y,
                           fixed,
                           main.parameters,
                           code,
                           config,
                           recover){


  u.beta <- (main.parameters$u.beta)
  v.beta <- (main.parameters$v.beta)
  alpha.K <- (main.parameters$alpha.K)
  alpha.L <- (main.parameters$alpha.L)

  for (method in recover$method){

    if (!is.null(recover$y)){

      Y.y <- scale(y)
      X.x <- scale(t(alpha.K)%*%MASS::ginv((alpha.K)%*%t(alpha.K))%*%(alpha.L)%*%(x)%*%u.beta%*%MASS::ginv(t(v.beta)%*%(v.beta))%*%t(v.beta))

      y[,which((colSums(recover$y)>0)==T)]  <- do.call('cbind',parallel::mclapply(mc.cores = 4, mc.silent = config$verbose, X = c(which((colSums(recover$y)>0)==T)), FUN = function(id_col){

        test_id.y <- as.logical(recover$y[,id_col])
        train_id.y <- as.logical(1 - recover$y[,id_col])

        sparse.y <- c(y[,id_col])

        if (any(test_id.y) & any(train_id.y)){

          y.covariate_predictors <- cbind(1,Y.y[train_id.y,-id_col]%*%main.parameters$v.beta[-id_col,])
          y.test_predictors <- cbind(1,Y.y[test_id.y,-id_col]%*%main.parameters$v.beta[-id_col,])

          x.covariate_predictors <- X.x[train_id.y,-id_col]%*%main.parameters$v.beta[-id_col,]
          x.test_predictors <- X.x[test_id.y,-id_col]%*%main.parameters$v.beta[-id_col,]

          covariate_predictors <- cbind(x.covariate_predictors,y.covariate_predictors)
          test_predictors <- cbind(x.test_predictors,y.test_predictors)

          if (method=="knn"){

            sparse.y[test_id.y] <-
              FNN::knn.reg(
                train = covariate_predictors,
                test = test_predictors,
                y = sparse.y[train_id.y],
                k = 2
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

        return(sparse.y)
      }))

      y <- as.matrix(y)


      recover$predict.y <- y


    }

    if (!is.null(recover$x)){

      X.x <- scale(x)
      Y.y <- scale(t(alpha.L)%*%MASS::ginv((alpha.L)%*%t(alpha.L))%*%(alpha.K)%*%(y)%*%v.beta%*%MASS::ginv(t(u.beta)%*%(u.beta))%*%t(u.beta))

      x[,which((colSums(recover$x)>0)==T)]  <- do.call('cbind',parallel::mclapply(mc.cores = 4,mc.silent = config$verbose,X = c(which((colSums(recover$x)>0)==T)),FUN = function(id_col){

        test_id.x <- as.logical(recover$x[,id_col])
        train_id.x <- as.logical(1 - recover$x[,id_col])

        sparse.x <- c(x[,id_col])

        if (any(test_id.x) & any(train_id.x)){

          y.covariate_predictors <- cbind(1,Y.y[train_id.x,-id_col]%*%main.parameters$u.beta[-id_col,])
          y.test_predictors <- cbind(1,Y.y[test_id.x,-id_col]%*%main.parameters$u.beta[-id_col,])

          x.covariate_predictors <- X.x[train_id.x,-id_col]%*%main.parameters$u.beta[-id_col,]
          x.test_predictors <- X.x[test_id.x,-id_col]%*%main.parameters$u.beta[-id_col,]

          covariate_predictors <- cbind(x.covariate_predictors,y.covariate_predictors)
          test_predictors <- cbind(x.test_predictors,y.test_predictors)

          if (method=="knn"){

            sparse.x[test_id.x] <-
              FNN::knn.reg(
                train = covariate_predictors,
                test = test_predictors,
                y = sparse.x[train_id.x],
                k = 2
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
        return(sparse.x)
      }))

      x <- as.matrix(x)

      recover$predict.x <- x

    }


  }


  return(recover)
}


