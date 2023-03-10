#' Fitting Differential Equations to Time Series Data
#'
#' @description Use numerical optimization to fit ordinary differential equations (ODEs) to time series data to examine the dynamic relationships between variables or the characteristics of a dynamical system. It can now be used to estimate the parameters of ODEs up to second order.
#'
#' @param data a data frame containing all model variables. The "time" column must be included.
#' @param model a string specifying the model to be used. The "=~" operator is used to define variables, with the name of the variable user defined on the left and the name of the variable in the data on the right. The '~' operator specifies a differential equation, with the dependent variable on the left and the independent variables on the right. See also ‘Details’.
#' @param guess an optional vector that allows the user to give starting values for the model parameters, including the model coefficients and variable initial states.
#' @param method an optional string indicating which optimizer to use. The default method is subject to the specific model. The available options are 'Nelder-Mead','L-BFGS-B','SANN' and 'BFGS'.
#' @param plot an optional TRUE or FALSE that TRUE will draw the plot of the raw data and the predicted values.
#' @details We suggest choosing the method by default. The guess values contain the coefficient of the model and initial values (the values of t0). Different models have different number of values.
#' @details Time(param) sequence for which output is wanted; the first value of times must be the initial time.
#' ```{r}
#' #eg1. An example of univariate second-order differential equation (damped oscillator model)
#' data('example1')
#' model1 <- '
#'    X =~ myX
#'    time =~ myTime
#'    X(2) ~ X(1) + X
#'   '
#' result1 <- defit(data = example1, model = model1)
#' # result1$table get the result
#' # names(result1) get all names of object
#' ```
#' @md
#' @returns object: directly type the defit object will print all results. The function summary is used to print the summary of all results, and the exact values of each result can be extracted by the "$" operator.
#' @returns userdata: the data that contains a sequence 'seq' starting from 1, the original time variable 'time', and all other variables user defined.
#' @returns parameter: the best set of parameters found, including parameter values, gradient, convergence, message and hessian matrix.
#' @returns predict: a dataframe of model predicted variable states at each time point.
#' @return r_squared: r_squared is the square of the correlation between the observed values and the predicted values, representing the proportion of variance explained by the model.
#' @return RMSE: RMSE (Root Mean Squared Error) is the standard deviation of the residuals.
#' @return SE: a symmetric matrix giving standard error of the model parameters.
#' @return equation: a string prints the estimated differential equations and initial states.
#' @return table: a summary table of parameter estimates and their corresponding SEs.
#' @return convergence: a message returns the result of the optimization convergence check.
#' @export defit
#'
#' @examples
#' #eg2. An example of bivariate first-order differential equation
#' data('example2')
#' model2 <- '
#'     # define variable
#'     X =~ myX
#'     Y =~ myY
#'
#'     # define time
#'     time =~ myTime
#'
#'     # define differential equation
#'     X(1) ~ X + Y
#'     Y(1) ~ Y + X
#'   '
#' result2 <- defit(data = example2, model = model2)
#' result2
#' # extract details and values
#' result2$summary()
#' result2$userdata
#' result2$parameter$par
#' result2$equation
#' result2$table
#'

defit <- function(data,model,guess=NULL,method=NULL,plot=FALSE){
  DeFit <- R6::R6Class("DeFit",
                       public=list(model = NULL,
                                   chooseModel = NULL,
                                   userdata = NULL,
                                   guess = NULL,
                                   method = NULL,
                                   parameter = NULL,
                                   predict = NULL,
                                   r_squared = NULL,
                                   RMSE = NULL,
                                   SE = NULL,
                                   equation = NULL,
                                   table = NULL,
                                   convergence = NULL,
                                   initialize = function(model = NULL) {
                                     self$model <- model
                                   },
                                   print = function(...){
                                     print(self$table)
                                   },
                                   summary = function(...){
                                     cat('-------Your Model-------\n')
                                     print(self$model)
                                     cat('-------Your Data (first 3 rows)-------\n')
                                     print(self$userdata[1:3,])
                                     cat('-------Your Results-------\n')
                                     cat('The fitted differential equation(s) \n')
                                     print(self$equation)
                                     cat('The optimization method is ',self$method,'\n',sep='')
                                     print(self$table)
                                     cat(self$convergence)
                                   }),
                       )
  model <- AdjustModel_func(model = model)
  outObject <- DeFit$new(model=model)
  # print('-------Your Model-------')
  # print(model)
  # print('-------Your Data (first 3 rows)-------')
  # print(data[1:3,])
  # print('-------Begin to estimate the parameters-------')
  chooseModel = JudgeModel_func(model=model)
  outObject$chooseModel <- chooseModel
  InitOption <- IntiOption_func(userdata=data,
                                model= model,
                                choosemodel=chooseModel,
                                guess=guess,
                                method=method)
  outObject$userdata <- InitOption$data
  outObject$guess <- InitOption$guess
  outObject$method <- InitOption$method
  calcDe <- CalcDe_func(data=InitOption$data,
                        model=model,
                        guess=InitOption$guess,
                        method=InitOption$method,
                        chooseModel=chooseModel)
  # print(calcDe$table)
  if(plot == TRUE){
    outPlot <- PlotDe_func(userdata = calcDe$userdata,calcdata = calcDe$Predictor)
    print(outPlot)
  }
  # print(calcDe$IsConvergence)
  outObject$parameter <- calcDe$Parameter
  outObject$predict <- calcDe$Predictor
  outObject$r_squared <- calcDe$Rsquared
  outObject$RMSE <- calcDe$Rmse
  outObject$SE <- calcDe$SE
  outObject$equation <- calcDe$DifferentialEquational
  outObject$table <- calcDe$table
  outObject$convergence <- calcDe$IsConvergence
  return(outObject)
}
utils::globalVariables(
  c("self")
)
