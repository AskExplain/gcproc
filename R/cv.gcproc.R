#' Run cross-validation on gcproc for regularisation
#'
#' Evaluate regularisation parameters, by running a search of the tuning parameter space
#'
#' @param x Matrix of dataset x
#' @param y Matrix of dataset y
#' @param fixed Fixed parameters from gcproc
#' @param config Configuration parameters from gcproc
#' @param recover Recover list from gcproc
#' @param anchors Transferring pre-trained model parameters (not required)
#' @param pivots Initialisation of model parameters (not required)
#'
#' @return  Recovered data from imputation or prediction, with the design matrices and any user input parameters and functions
#' @export
cv.gcproc <- function(x,
                      y,
                      reference = "y",
                      fixed = list(i_dim = F, j_dim = F),
                      config = gcproc::extract_config(verbose = F),
                      anchors = gcproc::extract_anchors_framework(verbose = F),
                      pivots = gcproc::extract_pivots_framework(verbose = F),
                      recover = gcproc::extract_recovery_framework(verbose = F)
)
{
  config.main <- config

  regularise.frame <- gcproc::extract_regularise_framework(F)

  lambda_sequence <- c(0.05,0.1,0.2)
  alpha_sequence <- c(0.5)

  score.list <- Inf
  main.score <- tail(score.list,1)

  config$min_iter <- 3
  config$max_iter <- 3
  config$verbose <- config.main$verbose

  for (j in 1:length(alpha_sequence)){

    regularise.frame$i_dim <- T
    regularise.frame$j_dim <- T

    for (i in 1:length(lambda_sequence)){

      regularise.frame$lambda <- lambda_sequence[i]
      regularise.frame$alpha <- alpha_sequence[j]

      gcproc.model <- gcproc::gcproc(
        x,
        y,
        reference = reference,
        fixed = fixed,
        config = config,
        anchors = anchors,
        pivots = pivots,
        recover = recover,
        regularise = regularise.frame
      )

      if (!is.infinite(gcproc.model$convergence.parameters$regularisation.penalty)){
        score.list <- c(score.list,
                        tail(gcproc.model$convergence.parameters$score.vec,1)
                        )
      } else {
        score.list <- c(score.list,Inf)
      }

      if (main.score > tail(score.list,1)){
        main.model <- gcproc.model
        main.score <- tail(score.list,1)
      }

    }
  }

  regularise.frame <- main.model$meta.parameters$regularise
  config <- config.main

  main.model <- gcproc::gcproc(
    x,
    y,
    reference = reference,
    fixed = fixed,
    config = config,
    anchors = anchors,
    pivots = pivots,
    recover = recover,
    regularise = regularise.frame
  )



  return(main.model)

}
