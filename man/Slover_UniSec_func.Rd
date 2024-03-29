% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Solver_UniSec_func.R
\name{Slover_UniSec_func}
\alias{Slover_UniSec_func}
\title{Core code of univariable second-order differential equation}
\usage{
Slover_UniSec_func(data, model, guess, method)
}
\arguments{
\item{data}{User's data}

\item{model}{model's class is dataframe.}

\item{guess}{Guess values that contain coefficient and initial values.}

\item{method}{"Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN" and "Brent"}
}
\value{
The result of optimization, SE, RMSE, r-squared, users' data, predictor data and output table.
}
\description{
Core code of univariable second-order differential equation
}
