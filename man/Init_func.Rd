% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Init_func.R
\name{Init_func}
\alias{Init_func}
\title{Initialize model Judgement variables and so on.}
\usage{
Init_func(userdata, model, guess, method, plot)
}
\arguments{
\item{userdata}{a data frame containing all model variables. The "time" column must be included.}

\item{model}{a string of model.}

\item{guess}{a list or a string. Guess the coefficients or initial values.}

\item{method}{a list or a string. The available options are 'Nelder-Mead','L-BFGS-B','SANN' and 'BFGS'.}

\item{plot}{TRUE or FALSE.}
}
\value{
a list
}
\description{
Initialize model Judgement variables and so on.
}
