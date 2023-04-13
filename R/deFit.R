#' Fitting Differential Equations to Time Series Data
#'
#' @description Use numerical optimization to fit ordinary differential equations (ODEs) to time series data to examine the dynamic relationships between variables or the characteristics of a dynamical system. It can now be used to estimate the parameters of ODEs up to second order, and can also apply to multilevel systems.
#'
#' @param data a data frame containing all model variables. The "time" column must be included.
#' @param model a string specifying the model to be used. The "=~" operator is used to define variables, with the name of the variable user defined on the left and the name of the variable in the data on the right. The '~' operator specifies a differential equation, with the dependent variable on the left and the independent variables on the right. See also ‘Details’.
#' @param guess an optional vector that allows the user to give starting values for the model parameters, including the model coefficients and variable initial states.
#' @param method an optional string indicating which optimizer to use. The default method is subject to the specific model. The available options are 'Nelder-Mead','L-BFGS-B','SANN' and 'BFGS'.
#' @param plot an optional TRUE or FALSE that TRUE will draw the plot of the raw data and the predicted values.
#' @details We suggest choosing the method by default. The guess values contain the coefficient of the model and initial values (the values of t0). Different models have different number of values.
#' @details Time(param) sequence for which output is wanted; the first value of times must be the initial time.
#' ```{r}
#' # eg1. An example of the univariate second-order differential equation (damped oscillator model)
#' data('example1')
#' model1 <- '
#'    X =~ myX
#'    time =~ myTime
#'    X(2) ~ X(1) + X
#'   '
#' result1 <- defit(data = example1, model = model1)
#' # result1$table get the result
#' # names(result1) get all names of object
#'
#' #--------------
#' # eg3. An example of the multilevel univariate second-order differential equation
#' data('example3')
#' model3 <- '
#'    X =~ current
#'    time =~ myTime
#'    X(2) ~ X(1) + X + (1 + X(1) + X | year)
#'    '
#' example3_use <- example3[(example3["year"] >= 2015)&(example3["year"] <= 2018),] # Note: select a subset of the data as an example.
#' example3_c <- scale_within(example3_use, model3) # note: centering X variable by year
#' result3 <- defit(data=example3_c,model = model3,plot=FALSE)
#'
#' #--------------
#' # eg4. An example of the multilevel bivariate first-order differential equations
#' data('example3')
#' model4 <- '
#'    X =~ current
#'    Y =~ expected
#'    time =~ myTime
#'    X(1) ~ X + Y + (1 + X + Y | year)
#'    Y(1) ~ X + Y + (1 + X + Y | year)
#'    '
#' example4_use <- example3[(example3["year"] >= 2015)&(example3["year"] <= 2018),] # Note: select a subset of the data as an example.
#' example4_c <- scale_within(example4_use, model4) # centering X and Y variable by year
#' result4 <- defit(data=example4_c,model = model4,plot=FALSE)
#'
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
#' # result2$summary()
#' # result2$userdata
#' # result2$parameter$par
#' # result2$equation
#' # result2$table
#' # result2$plot()
#'

defit <- function(data,model,guess=NULL,method=NULL,plot=FALSE){
  DeFit <- R6::R6Class("DeFit",
                       public=list(model = NULL,
                                   chooseModel = NULL,
                                   userdata = NULL,
                                   guess = NULL,
                                   guess2 = NULL,
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
                                   plot = function(...){
                                     outPlot <- PlotDe_func(userdata = self$userdata,
                                                            predictor = self$predict,
                                                            model=self$model,
                                                            chooseModel=self$chooseModel)
                                   },
                                   summary = function(...){
                                     cat('-------Your Model-------\n')
                                     print(self$model)
                                     cat('-------Your Data (first 3 rows)-------\n')
                                     print(self$userdata[1:3,])
                                     cat('-------Your Results-------\n')
                                     cat('The fitted differential equation(s) \n')
                                     for (i_equ in self$equation){
                                       print(i_equ)
                                     }
                                     cat('-------Your method-------\n')
                                     cat('The optimization method is ',self$method,'\n',sep='')
                                     cat('-------Summary table-------\n')
                                     print(self$table)
                                     cat(self$convergence)
                                     cat('\nr_square:',self$r_squared)
                                     cat('\nRMSE:',self$RMSE)
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
                        chooseModel=chooseModel,
                        guess2=InitOption$guess2)
  # print(calcDe$table)
  # print(calcDe$IsConvergence)
  outObject$parameter <- calcDe$Parameter

  var_model = outObject$model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]

  outObject$predict <- calcDe$Predictor
  for (var_mode_i in var_model[,'field']) {
    names(outObject$predict)[names(outObject$predict) == paste('solver_',var_mode_i,sep = '')] = paste(var_mode_i,'_hat',sep='')
  }

  if(plot == TRUE){
    outPlot <- PlotDe_func(userdata = outObject$userdata,
                           predictor = outObject$predict,
                           model=outObject$model,
                           chooseModel=outObject$chooseModel)
    #print(outPlot)
  }

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
