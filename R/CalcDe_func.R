#' Choose the core code by the model.
#'
#' @param data Users' data
#' @param model The original string users defined.
#' @param guess The guess values users input.
#' @param method The method users selected.
#' @param chooseModel c(2,1): Bivariate first-order differential equation; c(1,2): Univariable second-order differential equation.
#'
#' @return The result of optimization,SE,RMSE,r-squared,users' data,predictor data and output table.
CalcDe_func <- function(data,model,guess=c(0,0,0,0,0,0),method,chooseModel=NULL){

  if(all(chooseModel == c(2,1))){
    # print('Bivariate first-order differential equation')
    solverdata = Slover_BinFirst_func(data=data,model,guess,method)
  }else if(all(chooseModel==c(1,2))){
    # print('Univariate second-order differential equation')
    solverdata = Slover_UniSec_func(data=data,model,guess,method)
  }else{
    warning('Your model is not supported!')
  }
  return(solverdata)
}
