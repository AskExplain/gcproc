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


  if (recover$task == "regression"){


    matrix.projection <- c("matrix.projection"%in%recover$method)
    knn.reg <- c("knn.reg"%in%recover$method)

    for (i in 1:length(data_list)){

      if (!is.null(recover$design.list[[i]])){

        if (matrix.projection){


          if (min(data_list[[i]])==0){
            x <- log(as.matrix(data_list[[i]])+1)
            to_exp <- T
          } else {
            x <- as.matrix(data_list[[i]])
            to_exp <- F
          }

          recover$encoded_covariate <- lapply(c(1:length(data_list))[-i],function(X){
            transformed.data <- as.matrix(MASS::ginv((main.parameters[[X]]$alpha)%*%t(main.parameters[[X]]$alpha))%*%(main.parameters[[X]]$alpha)%*%as.matrix(data_list[[X]])%*%(main.parameters[[X]]$beta)%*%MASS::ginv(t((main.parameters[[X]]$beta))%*%(main.parameters[[X]]$beta)))
          })

          decoded_covariate <- cbind(1,scale(Reduce('+',lapply(c(1:length(recover$encoded_covariate)),function(X){
            t(main.parameters[[i]]$alpha)%*%recover$encoded_covariate[[X]]
          }))))

          samples_with_missing_points <- which((rowSums(recover$design.list[[i]])>0)==T)
          covariate_predictors <-  decoded_covariate[-samples_with_missing_points,]
          test_predictors <- decoded_covariate[samples_with_missing_points,]


          elements_with_missing_points <- which((recover$design.list[[i]]>0)==T,arr.ind = T)
          x[elements_with_missing_points]  <- (((test_predictors)%*%(MASS::ginv(t(covariate_predictors)%*%(covariate_predictors))%*%t(covariate_predictors)%*%(x[-samples_with_missing_points,]))))[elements_with_missing_points]

          data_list[[i]] <- recover$predict.list[[i]] <- x <- if(to_exp){exp(as.matrix(x))-1}else{as.matrix(x)}

        }


        if (knn.reg){

          if (min(data_list[[i]])==0){
            x <- log(as.matrix(data_list[[i]])+1)
            to_exp <- T
          } else {
            x <- as.matrix(data_list[[i]])
            to_exp <- F
          }

          recover$encoded_covariate <- lapply(c(1:length(data_list))[-i],function(X){
            transformed.data <- as.matrix(MASS::ginv((main.parameters[[X]]$alpha)%*%t(main.parameters[[X]]$alpha))%*%(main.parameters[[X]]$alpha)%*%as.matrix(data_list[[X]])%*%(main.parameters[[X]]$beta)%*%MASS::ginv(t((main.parameters[[X]]$beta))%*%(main.parameters[[X]]$beta)))
          })

          decoded_covariate <- cbind(1,scale(Reduce('+',lapply(c(1:length(recover$encoded_covariate)),function(X){
            t(main.parameters[[i]]$alpha)%*%recover$encoded_covariate[[X]]
          }))))

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

          elements_with_missing_points <- which((recover$design.list[[i]]>0)==T,arr.ind = T)
          pred.x <- x
          pred.x[samples_with_missing_points,] <- pred
          x[elements_with_missing_points] <- pred.x[elements_with_missing_points]

          data_list[[i]] <- recover$predict.list[[i]] <- x <- if(to_exp){exp(as.matrix(x))-1}else{as.matrix(x)}

        }


        if (!is.null(recover$fn)){

          if (min(x)==0){
            x <- log(as.matrix(data_list[[i]])+1)
            to_exp <- T
          } else {
            x <- as.matrix(data_list[[i]])
            to_exp <- F
          }

          recover$encoded_covariate <- lapply(c(1:length(data_list))[-i],function(X){
            transformed.data <- as.matrix(MASS::ginv((main.parameters[[X]]$alpha)%*%t(main.parameters[[X]]$alpha))%*%(main.parameters[[X]]$alpha)%*%as.matrix(data_list[[X]])%*%(main.parameters[[X]]$beta)%*%MASS::ginv(t((main.parameters[[X]]$beta))%*%(main.parameters[[X]]$beta)))
          })

          decoded_covariate <- cbind(1,scale(Reduce('+',lapply(c(1:length(recover$encoded_covariate)),function(X){
            t(main.parameters[[i]]$alpha)%*%recover$encoded_covariate[[X]]
          }))))


          samples_with_missing_points <- which((rowSums(recover$design.list[[i]])>0)==T)
          covariate_predictors <-  decoded_covariate[-samples_with_missing_points,]
          test_predictors <- decoded_covariate[samples_with_missing_points,]

          pred <- recover$fn(train = covariate_predictors, test = test_predictors, y = x[-samples_with_missing_points,], parameters = recover$param)

          elements_with_missing_points <- which((recover$design.list[[i]]>0)==T,arr.ind = T)
          pred.x <- x
          pred.x[samples_with_missing_points,] <- pred
          x[elements_with_missing_points] <- pred.x[elements_with_missing_points]

          data_list[[i]] <- recover$predict.list[[i]] <- x <- if(to_exp){exp(as.matrix(x))-1}else{as.matrix(x)}

        }


      }


    }


  }
  if (recover$task == "classification"){


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

          unlabel.decoded_covariate <- scale(Reduce('+',lapply(c(1:length(recover$encoded_covariate)),function(X){
            t(main.parameters[[i]]$alpha)%*%recover$encoded_covariate[[X]]
          })))

          labels <- recover$labels

          recover$predict.list[[j]][[i]] <- apply((unlabel.decoded_covariate)%*%t(label.decoded_covariate),1,function(X){names(sort(table(labels[order(X,decreasing = T)[1]]))[1])})

        }

      }
    }



  }

  return(recover)
}


