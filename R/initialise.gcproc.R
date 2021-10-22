#' @export
initialise.gcproc <- function(data_list,
                              config,
                              covariate,
                              transfer,
                              join,
                              pivots){

  index <- list()
  index$code_indicator <- unique(do.call('c',lapply(c(1:length(covariate$factor)),function(X){c(unique(colnames(covariate$factor[[X]])))})))



  if (config$verbose){
    print(paste("Initialising data with : ",config$init,sep=""))
  }

  main.code <- list(code=list(),encode=list())
  main.index <- list()
  main.proportion <- list()
  main.parameters <- list(alpha = list(), beta = list())

  for (i in 1:length(data_list)){
    alpha.list <- list()
    beta.list <- list()

    encode.list <- list()
    code.list <- list()

    for (j in c(1:length(index$code_indicator))){

      if (covariate$fix){
        main.proportion[[i]] <- covariate$factor[[i]]
      } else {
        main.proportion[[i]] <- array(runif(dim(data_list[[i]])[1]*length(index$code_indicator)),dim=c(dim(data_list[[i]])[1],length(index$code_indicator)))
      }
      main.proportion[[i]] <- main.proportion[[i]] / rowSums(main.proportion[[i]])


      initial.param <-initialise.parameters(x = as.matrix(data_list[[i]]),transfer = transfer, i_dim=config$i_dim,j_dim=config$j_dim,init=config$init,verbose=config$verbose)

      # Check anchoring parameters
      alpha <- initial.param$pivot_x.sample
      beta <- initial.param$pivot_x.feature


      if (transfer$fix){
        encode.d <- transfer$code$encode
        code.d <- transfer$code$code
      } else {

        encode.d <- (alpha%*%as.matrix(data_list[[i]])%*%(beta))
        code.d <- ((pinv(t(alpha))%*%(encode.d)%*%pinv((beta))))

      }

      # for (init_run in c(1:3)){
      #
      #   beta.encode_projection <- (beta%*%t(encode.d))
      #   beta.code_projection <- code.d%*%t(beta)
      #   beta.decode_projection <- beta.code_projection%*%beta.encode_projection
      #
      #   alpha <- (t(((as.matrix(data_list[[i]]))%*%beta.encode_projection%*%t(beta.decode_projection)%*%pinv(t(beta.decode_projection)))))
      #
      #   alpha.encode_projection <- (t(encode.d)%*%alpha)
      #   alpha.code_projection <- t(alpha)%*%code.d
      #   alpha.decode_projection <- alpha.encode_projection%*%alpha.code_projection
      #
      #   beta <- (t(pinv((alpha.decode_projection))%*%t(alpha.decode_projection)%*%(alpha.encode_projection)%*%(as.matrix(data_list[[i]]))))
      #
      #
      #   encode.d <- (alpha%*%as.matrix(data_list[[i]])%*%(beta))
      #   code.d <- ((pinv(t(alpha))%*%(encode.d)%*%pinv((beta))))
      #
      #
      #   if (!covariate$fix){
      #     pys <- main.proportion[[i]]
      #     for (X in c(1:dim(main.proportion[[i]])[2])){
      #       x.beta <- (as.matrix(data_list[[i]]))%*%(beta)
      #       x.decode <- t(alpha)%*%code.d%*%t(beta)%*%(beta)
      #
      #       pys[,X] <- log((pys[,X])) + mclust::dmvnorm(data = (x.beta - x.decode),sigma = diag(diag(t(x.decode)%*%(x.decode)/dim(x.decode)[1])) ,log = T)
      #     }
      #
      #     pys_max <- apply(pys, 1, max)
      #     pys <- sweep(pys, 1, pys_max, '-')
      #     pys <- exp(pys)
      #     main.proportion[[i]] <- sweep(pys, 1, rowSums(pys), '/')
      #     colnames(main.proportion[[i]]) <- index$code_indicator
      #
      #   }
      #
      # }


      alpha.list <- c(alpha.list,list(alpha))
      beta.list <- c(beta.list,list(beta))
      code.list <- c(code.list,list(code.d))
      encode.list <- c(encode.list,list(encode.d))

    }


    main.code = list(
      encode = encode.d,
      code = code.list
    )


    main.parameters$alpha[[i]] <- alpha.list
    main.parameters$beta[[i]] <- beta.list


  }


  return(
    list(
      main.parameters = main.parameters,
      main.code = main.code,
      main.index = main.index,
      main.proportion = main.proportion
    )
  )

}


#' @export
initialise.parameters <- function(x,transfer,i_dim,j_dim,init="svd",verbose=F){

  x <- Matrix::Matrix(x,sparse=T)

  set.seed(1)

  if (init=="random"){
    param.beta <- if(is.null(transfer$beta)){array(rnorm(config$j_dim),dim=c(dim(x)[2],config$j_dim))}else{transfer$beta}
    param.alpha = if(is.null(transfer$alpha)){array(rnorm(config$i_dim),dim=c(config$i_dim,dim(x)[1]))}else{transfer$alpha}
  } else {
    cov_x <- corpcor::cov.shrink(x,verbose = F)
    cov_tx <- corpcor::cov.shrink(Matrix::t(x),verbose = F)
  }

  if (init=="svd"){
    param.beta.svd <- irlba::irlba(
      cov_x,j_dim,verbose = F)
    rm(cov_x)

    param.beta <- if(is.null(transfer$beta)){param.beta.svd$v}else{transfer$beta}


    param.alpha.J.svd <- irlba::irlba(
      cov_tx,i_dim,verbose = F)
    rm(cov_tx)

    param.alpha = if(is.null(transfer$alpha)){t(param.alpha.J.svd$u)}else{transfer$alpha}

  }

  pivots <- list(
    pivot_x.sample = as.matrix(param.alpha),
    pivot_x.feature = as.matrix(param.beta)  )
  return(pivots)

}
