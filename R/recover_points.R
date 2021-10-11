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
                           recover){

  for (task in recover$task){

    if ("regression" %in% task){

      for (method in recover$method){


        for (i in 1:length(data_list)){

          if (!is.null(recover$design.list[[i]])){

            if ("decode" %in% method){

              x <- as.matrix(data_list[[i]])

              row_with_missing_points <- which((rowSums(recover$design.list[[i]])>0)==T,arr.ind = T)
              column_with_missing_points <- which((colSums(recover$design.list[[i]])>0)==T,arr.ind = T)

              # pred <- t(main.parameters[[i]]$alpha)%*%main.code$code%*%t(main.parameters[[i]]$beta)

              x[row_with_missing_points,column_with_missing_points]  <- pred[row_with_missing_points,column_with_missing_points]

              data_list[[i]] <- recover$predict.list[[i]] <- x


            }


            if ("matrix.projection" %in% method){

              x <- as.matrix(data_list[[i]])


              if (is.null(recover$encoded_covariate)){
                recover$encoded_covariate <- lapply(c(1:length(data_list)),function(X){
                  transformed.data <- as.matrix(data_list[[X]])%*%(main.parameters[[X]]$beta)
                  return(transformed.data)
                })
              }


              decoded_covariate <- cbind(1,
                                         scale(Reduce('+',lapply(c(1:length(recover$encoded_covariate)),function(X){
                                           recover$encoded_covariate[[X]]
                                         })))              )

              samples_with_missing_points <- which((rowSums(recover$design.list[[i]])>0)==T)
              covariate_predictors <-  decoded_covariate[-samples_with_missing_points,]
              test_predictors <- decoded_covariate[samples_with_missing_points,]

              elements_with_missing_points <- which((recover$design.list[[i]]>0)[samples_with_missing_points,]==T,arr.ind = T)
              x[samples_with_missing_points,][elements_with_missing_points]  <- (((test_predictors)%*%(MASS::ginv(t(covariate_predictors)%*%(covariate_predictors))%*%t(covariate_predictors)%*%(x[-samples_with_missing_points,]))))[elements_with_missing_points]

              data_list[[i]] <- recover$predict.list[[i]] <- x
            }


            if ("knn" %in% method){


              x <- as.matrix(data_list[[i]])

              if (is.null(recover$encoded_covariate)){
                recover$encoded_covariate <- lapply(c(1:length(data_list)),function(X){
                  transformed.data <- as.matrix(data_list[[X]])%*%(main.parameters[[X]]$beta)
                  return(transformed.data)
                })
              }


              decoded_covariate <- cbind(1,
                                         scale(Reduce('+',lapply(c(1:length(recover$encoded_covariate)),function(X){
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

              x <- as.matrix(data_list[[i]])

              if (is.null(recover$encoded_covariate)){
                recover$encoded_covariate <- lapply(c(1:length(data_list)),function(X){
                  transformed.data <- as.matrix(data_list[[X]])%*%(main.parameters[[X]]$beta)
                  return(transformed.data)
                })
              }


              decoded_covariate <- cbind(1,
                                         scale(Reduce('+',lapply(c(1:length(recover$encoded_covariate)),function(X){
                                           recover$encoded_covariate[[X]]
                                         })))              )


              samples_with_missing_points <- which((rowSums(recover$design.list[[i]])>0)==T)
              covariate_predictors <-  decoded_covariate[-samples_with_missing_points,]
              test_predictors <- decoded_covariate[samples_with_missing_points,]

              pred <- recover$fn(train = covariate_predictors, test = test_predictors, y = x[-samples_with_missing_points,], parameters = recover$param)

              elements_with_missing_points <- which((recover$design.list[[i]]>0)[samples_with_missing_points,]==T,arr.ind = T)
              x[samples_with_missing_points,][elements_with_missing_points]  <- (pred)[elements_with_missing_points]

              data_list[[i]] <- recover$predict.list[[i]] <- x
            }


          }


        }
      }

    }

    if ("classification" %in% task){


      label.projection <- c(recover$method=="label.projection")

      if (label.projection){

        for (j in which(recover$design.list==0)){

          label.decoded_covariate <- t(main.parameters[[j]]$alpha)

          for (i in which(recover$design.list==1)){

            unlabel.decoded_covariate <- t(main.parameters[[i]]$alpha)

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
