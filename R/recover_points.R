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
                           code,
                           main.parameters,
                           config,
                           recover){

  for (task in recover$task){

    if ("regression" %in% task){

      for (method in recover$method){

        if ("matrix.projection" %in% method){

          covariate_predictors <-  cbind(1,scale(as.matrix(data_list$train.x)%*%(main.parameters[[which(names(data_list)=="train.x")]]$beta)))
          test_predictors <- cbind(1,scale(as.matrix(data_list$test.x)%*%(main.parameters[[which(names(data_list)=="test.x")]]$beta)))


          x  <- (((test_predictors)%*%(MASS::ginv(t(covariate_predictors)%*%(covariate_predictors))%*%t(covariate_predictors)%*%(data_list$train.y))))

          data_list$test.y <- recover$predict.list <- x
        }


        if ("knn" %in% method){


          covariate_predictors <-  cbind(1,scale(as.matrix(data_list$train.x)%*%(main.parameters[[which(names(data_list)=="train.x")]]$beta)))
          test_predictors <- cbind(1,scale(as.matrix(data_list$test.x)%*%(main.parameters[[which(names(data_list)=="test.x")]]$beta)))

          knn_ix <- FNN::get.knnx(
            covariate_predictors,
            test_predictors,
            k = 5
          )$nn.index

          pred <- (data_list$train.y)[knn_ix[, 1], , drop = FALSE]
          if (5 > 1) {
            for (k in seq(2, 5)) {
              pred <- pred + (data_list$train.y)[knn_ix[, k], , drop = FALSE]
            }
          }
          pred <- pred / 5

          data_list$test.y <- recover$predict.list <- as.matrix(pred)

        }

      }

    }

    if ("classification" %in% task){


      label.projection <- c(recover$method=="label.projection")

      if (label.projection){

        for (j in which(recover$design.list==0)){

          recover$encoded_covariate <- lapply(c(1:length(data_list)),function(X){
            transformed.data <- as.matrix(MASS::ginv((main.parameters[[X]]$alpha)%*%t(main.parameters[[X]]$alpha))%*%(main.parameters[[X]]$alpha)%*%as.matrix(data_list[[X]])%*%(main.parameters[[X]]$beta)%*%MASS::ginv(t((main.parameters[[X]]$beta))%*%(main.parameters[[X]]$beta)))
          })

          label.decoded_covariate <- scale(Reduce('+',lapply(c(1:length(recover$encoded_covariate)),function(X){
            t(main.parameters[[j]]$alpha)%*%recover$encoded_covariate[[X]]
          })))

          for (i in which(recover$design.list==1)){

            unlabel.decoded_covariate <- scale(Reduce('+',lapply(c(1:length(recover$encoded_covariate))[-j],function(X){
              t(main.parameters[[i]]$alpha)%*%recover$encoded_covariate[[X]]
            })))

            labels <- recover$labels

            recover$predict.list[[j]][[i]] <- apply((unlabel.decoded_covariate)%*%t(label.decoded_covariate),1,function(X){names(sort(table(labels[order(X,decreasing = T)[1]]))[1])})

          }

        }
      }



    }

  }

  return(recover)
}


