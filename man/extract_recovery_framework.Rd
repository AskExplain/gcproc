% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_config.R
\name{extract_recovery_framework}
\alias{extract_recovery_framework}
\title{Extract recovery framework to put into gcproc}
\usage{
extract_recovery_framework(verbose = T)
}
\arguments{
\item{task}{Allows user to specify either a regression or classification task}

\item{method}{The algorithm for the task (Options are regression: "knn.reg","matrix.projection", -- provide your own --   ;   classification: "label.projection")}

\item{design.list}{A list of design structures where each element is given a 1 to indicate the test set, 0 indicates the train set.}

\item{labels}{For classification, these are the pre-defined labels}

\item{predict.list}{This will be filled in by gpcroc with the predictions and return a prediction for indicated design matrices only. Leave as NULL to begin.}
}
\value{
Prediction framework for gcproc
}
\description{
Can recover data points by imputing or predicting missing values
}
