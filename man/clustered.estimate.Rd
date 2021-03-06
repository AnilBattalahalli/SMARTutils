% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustered.R
\name{clustered.estimate}
\alias{clustered.estimate}
\title{Function to estimate beta, variance of beta_hat, and ICC}
\usage{
clustered.estimate(formula, data, est_var = T, weights = NULL)
}
\arguments{
\item{formula}{formula for the estimation}

\item{data}{a data.frame object containing cluster ID, participant ID, A1, R, A2, Y with
column names being i, j, A1, R, A2, Y}

\item{est_var}{T/F to estimate the variance or not}

\item{weights}{a data frame containing cluster level weights with columns `i` and `W` for cluster ID and weights respectively}
}
\value{
estimates object containing all the estimated parameters.
}
\description{
Function to estimate beta, variance of beta_hat, and ICC
}
