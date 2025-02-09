% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Solver_BiFirst_func.R
\name{Solver_BiFirst_func}
\alias{Solver_BiFirst_func}
\title{Solver of bivariate first-order differential equation}
\usage{
Solver_BiFirst_func(userdata, var_model, guess, method)
}
\arguments{
\item{userdata}{a data frame containing all model variables. The "time" column must be included.}

\item{var_model}{a dataframe containing equations.}

\item{guess}{a list or a string. Guess the coefficients or initial values.}

\item{method}{a list or a string. The available options are 'Nelder-Mead','L-BFGS-B','SANN' and 'BFGS'.}
}
\value{
a list
}
\description{
Solver of bivariate first-order differential equation
}
