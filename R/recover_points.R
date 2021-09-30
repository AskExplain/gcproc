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
recover_points <- function(data_list,
                           code,
                           main.parameters,
                           config,
                           recover){


  for (method in recover$method){

    for (i in 1:length(data_list)){

      x <- as.matrix(data_list[[i]])

      if (!is.null(recover$design.list[[i]])){
        print(c("check pre covariate",i))

        recover$covariate <- scale(Reduce("+",lapply(c(1:length(data_list)),function(X){
          transformed.data <- (t(main.parameters[[i]]$alpha)%*%MASS::ginv((main.parameters[[X]]$alpha)%*%t(main.parameters[[X]]$alpha))%*%(main.parameters[[X]]$alpha)%*%as.matrix(data_list[[X]])%*%(main.parameters[[X]]$beta)%*%MASS::ginv(t(main.parameters[[X]]$beta)%*%(main.parameters[[X]]$beta)))
        })))

        print(c("check post covariate",i))

        x[,which((colSums(recover$design.list[[i]])>0)==T)]  <- do.call('cbind',lapply(X = c(which((colSums(recover$design.list[[i]])>0)==T)),FUN = function(id_col){

          print(c("check pre design",i))

          test_id.x <- as.logical(recover$design.list[[i]][,id_col])
          train_id.x <- as.logical(1 - recover$design.list[[i]][,id_col])

          print(c("check post design",i))

          if (min(x)==0){
            sparse.x <- log(x[,id_col]+1)
            to_exp <- T
          } else {
            sparse.x <- x[,id_col]
            to_exp <- F
          }

          print(c("check pre predict inner",i))

          if (any(test_id.x) & any(train_id.x)){

            covariate_predictors <- cbind(1,recover$covariate[train_id.x,])
            test_predictors <- cbind(1,recover$covariate[test_id.x,])

            if (method=="knn"){

              sparse.x[test_id.x] <-
                FNN::knn.reg(
                  train = covariate_predictors,
                  test = test_predictors,
                  y = sparse.x[train_id.x],
                  k = 20
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

            print(c("check post predict inner",i))

          }
          return(if(to_exp){exp(sparse.x)-1}else{sparse.x})
        }))

        x <- as.matrix(x)

        recover$predict.list[[i]] <- x

        data_list[[i]] <- x

      }


    }


  }


  return(recover)
}


