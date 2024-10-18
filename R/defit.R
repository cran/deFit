#' Fitting Differential Equations to Time Series Data
#'
#' @description Use numerical optimization to fit ordinary differential equations (ODEs) to time series data to examine the dynamic relationships between variables or the characteristics of a dynamical system. It can now be used to estimate the parameters of ODEs up to second order, and can also apply to multilevel systems.
#'
#' @param data a data frame containing all model variables. The "time" column must be included.
#' @param model a string specifying the model to be used. The "=~" operator is used to define variables, with the name of the variable user defined on the left and the name of the variable in the data on the right. The '~' operator specifies a differential equation, with the dependent variable on the left and the independent variables on the right. See also ‘Details’.
#' @param guess an optional vector that allows the user to give starting values for the model parameters, including the model coefficients and variable initial states.
#' @param method an optional string indicating which optimizer to use. The default method is subject to the specific model. The available options are 'Nelder-Mead','L-BFGS-B','SANN' and 'BFGS'.
#' @param plot an optional TRUE or FALSE that TRUE will draw the plot of the raw data and the predicted values.
#' @param fixed_first an optional TRUE or FALSE that TRUE will estimate the multilevel model parameters using a two-step approach.
#' @details We suggest choosing the method by default. The guess values contain the coefficient of the model and initial values (the values of t0). Different models have different number of values.
#' @details Time(param) sequence for which output is wanted; the first value of times must be the initial time.
#' \preformatted{
#' # eg1. An example of the univariate second-order differential equation (damped oscillator model)
#' data('example1')
#' model1 <- '
#'    X =~ myX
#'    time =~ myTime
#'    X(2) ~ X + X(1)
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
#'    X(2) ~ X + X(1) + (1 + X + X(1) | year)
#'    '
#' example3_use <- example3[(example3["year"] >= 2015)&(example3["year"] <= 2018),] # Note: select a subset of the data as an example.
#' example3_c <- Scale_within(example3_use, model3) # note: centering X variable by year
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
#' example4_c <- Scale_within(example4_use, model4) # centering X and Y variable by year
#' result4 <- defit(data=example4_c,model = model4,plot=FALSE)
#' }
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
#'     Y(1) ~ X + Y
#'   '
#' result2 <- defit(data = example2, model = model2,method='Nelder-Mead')
#' # Note: the method argument will override the default "L-BFGS-B" method
#' result2
#' # #extract details and values
#' # result2$summary()
#' # result2$userdata
#' # result2$parameter$par
#' # result2$equation
#' # result2$table
#' # result2$plot()
#'
defit <- function(data,model,guess=NULL,method=NULL,plot=FALSE,fixed_first=TRUE){
  DeFit <- R6::R6Class("DeFit",
                       public=list(model = NULL,
                                   finalModel = NULL,
                                   userdata = NULL,
                                   guess = NULL,
                                   method = NULL,
                                   parameter = NULL,
                                   predict = NULL,
                                   predict_fixed = NULL,
                                   r_squared = NULL,
                                   RMSE = NULL,
                                   SE = NULL,
                                   Hessian = NULL,
                                   equation = NULL,
                                   table = NULL,
                                   var_model = NULL,
                                   field_model = NULL,
                                   order_model = NULL,
                                   multi_model = NULL,
                                   random_effects = NULL,
                                   obj = NULL,
                                   initialize = function(model = NULL) {
                                     # if the model is undefined,
                                     model = '
                                        X =~ myX
                                        Y =~ myY
                                        time =~ myTime
                                        X(1) ~ X + Y
                                        Y(1) ~ X + Y
                                      '
                                     self$model <- model
                                   },
                                   print = function(...){
                                     print(self$table)
                                   },
                                   plot = function(...){
                                     outPlot <- Plot_func(userdata = self$userdata,
                                                          predict = self$predict,
                                                          modelDF=self$model,
                                                          var_model=self$var_model,
                                                          field_model=self$field_model,
                                                          order_model=self$order_model,
                                                          multi_model=self$multi_model)
                                   },
                                   summary = function(...){
                                     cat('-------Your Model-------\n')
                                     print(self$model)
                                     cat('-------Your Data (first 3 rows)-------\n')
                                     print(self$userdata[1:3,])
                                     cat('-------Your Results-------\n')
                                     # cat('The fitted differential equation(s) \n')
                                     # for (i_equ in self$equation){
                                     #   print(i_equ)
                                     # }
                                     # print(self$equation)
                                     cat('-------Your method-------\n')
                                     cat('The optimization method is ',self$method,'\n',sep='')
                                     cat('-------Summary table-------\n')
                                     print(self$table)
                                     cat('\nr_square:',self$r_squared)
                                     cat('\nRMSE:',self$RMSE)
                                     try({
                                       cat('\nNeg2LL:',self$obj$Neg2LL[[1]][1])
                                     },
                                         silent = FALSE)
                                   }),
  )
  init_list = Init_func(userdata=data,
                        model,
                        guess,
                        method,
                        plot)
  if(init_list$multi_model){
    if(length(init_list$field_model) == 2){
      if(init_list$order_model ==1){
        # multilevel bivariate first-order differential equation
        init_list$fixed_first = fixed_first
        out_solve = Solver_MultiBiFirst_func(init_list)
      }else{
        stop("********************\n",
             "Not supported now!.\n",
             "Two variables model is only supported by multilevel bivariate first-order differential model.\n",
             "********************")
      }
    }else if(length(init_list$field_model) == 1){
      if(init_list$order_model ==2){
        # multilevel univariate second-order differential equation
        init_list$fixed_first = fixed_first
        out_solve = Solver_MultiUniSecond_func(init_list)
      }else{
        stop("********************\n",
             "Not supported now!.\n",
             "One variable is only supported by multilevel univariate second-order differential model.\n",
             "********************")
      }
    }else{
      stop("********************\n",
           "Not supported now!.\n",
           "********************")
    }
  }else{
    # print('not a multilevel model')
    if(length(init_list$field_model) == 2){
      if(init_list$order_model ==1){
        # bivariate first-order differential equational
        out_solve = Solver_BiFirst_func(init_list$userdata,
                                        init_list$var_model,
                                        init_list$guess,
                                        init_list$method)
      }else{
        stop("********************\n",
             "Not supported now!.\n",
             "Two variables model is only supported by bivariate first-order differential model.\n",
             "********************")
      }
    }else if(length(init_list$field_model) == 1){
      if(init_list$order_model ==2){
        out_solve = Solver_UniSecond_func(init_list$userdata,
                                        init_list$var_model,
                                        init_list$guess,
                                        init_list$method)
      }else{
        stop("********************\n",
             "Not supported now!.\n",
             "One variables model is only supported by univariate second-order differential model.\n",
             "********************")
      }
    }else{
      stop("********************\n",
           "Not supported now!.\n",
           "********************")
    }
  }
  out_info = Info_func(solve_data = out_solve['solve_data'],
                       userdata = out_solve['userdata'][[1]],
                       predict = out_solve['predict_data'][[1]],
                       var_model = init_list$var_model,
                       table = out_solve['table'][[1]])
  Obj_defit = DeFit$new(model = model)
  Obj_defit$userdata = out_solve$userdata
  Obj_defit$model = init_list$modelDF
  Obj_defit$guess = init_list$guess[[1]]
  Obj_defit$method = init_list$method[[1]]
  Obj_defit$parameter = out_solve[1][[1]]
  Obj_defit$equation = unlist(out_solve[5][[1]])
  Obj_defit$r_squared = out_info[1][[1]]
  Obj_defit$RMSE = out_info[2][[1]]
  Obj_defit$Hessian = out_info[3][[1]]
  Obj_defit$table = out_info[4][[1]]
  Obj_defit$var_model = init_list$var_model
  Obj_defit$predict = out_solve['predict_data'][[1]]
  Obj_defit$predict_fixed = out_solve['predict_fixed_data'][[1]]


  Obj_defit$field_model = init_list$field_model
  Obj_defit$order_model = init_list$order_model
  Obj_defit$multi_model = init_list$multi_model
  Obj_defit$random_effects = out_solve['random_effects'][[1]]
  Obj_defit$obj = out_solve

  if(plot){
    # plot is true will draw graph
    outPlot <- Plot_func(userdata = out_solve$userdata,
                         predict = out_solve['predict_data'][[1]],
                         modelDF=init_list$modelDF,
                         var_model=init_list$var_model,
                         field_model=init_list$field_model,
                         order_model=init_list$order_model,
                         multi_model=init_list$multi_model)
  }
  return(Obj_defit)
}
utils::globalVariables(
  c("self")
)
