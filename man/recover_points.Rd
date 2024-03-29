% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recover_points.R
\name{recover_points}
\alias{recover_points}
\title{Recover missing data points via imputation or prediction}
\usage{
recover_points(data_list, code, main.parameters, config, recover)
}
\arguments{
\item{data_list}{A list of datasets (matrix or tensor etc.)}

\item{code}{Code parameters from gcproc}

\item{main.parameters}{Main parameters from gcproc}

\item{config}{Configuration parameters from gcproc}

\item{recover}{Recover list from gcproc}
}
\value{
Recovered data from imputation or prediction, with the design matrices and any user input parameters and functions
}
\description{
Main function to recover missing points. Used internally, but can be re-purposed by the user.
}
