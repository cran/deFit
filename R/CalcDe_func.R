#' Choose the core code by the model.
#'
#' @param data Users' data
#' @param model The original string users defined.
#' @param guess The guess values users input.
#' @param method The method users selected.
#' @param chooseModel c(2,1): Bivariate first-order differential equation; c(1,2): Univariable second-order differential equation.
#' @param guess2 Use for Multilevel
#'
#' @return The result of optimization,SE,RMSE,r-squared,users' data,predictor data and output table.
CalcDe_func <- function(data,model,guess=c(0,0,0,0,0,0),method,chooseModel=NULL,guess2){

  if(all(chooseModel == c(2,1))){
    # print('Bivariate first-order differential equation')
    if(all(is.na(model[,'fixRand']) & all(is.na(model[,'subject'])))){
      # cat('Random effects and fixed effects are not defined\n')
      solverdata = Solver_BinFirst_func(data=data,model,guess,method)
    }else{
      cat('Random effects and fixed effects have been defined\n')
      solverdata = Solver_MultiBiFirst_func(data,model,guess,method,guess2)
    }
  }else if(all(chooseModel==c(1,2))){
    if(all(is.na(model[,'fixRand']) & all(is.na(model[,'subject'])))){
      # cat('Random effects and fixed effects are not defined\n')
      solverdata = Slover_UniSec_func(data=data,model,guess,method)
    }else{
      cat('Random effects and fixed effects have been defined\n')
      solverdata = Solver_MultiUniSec_func(data=data,model,guess,method,guess2)
    }
  }else{
    warning('Your model is not supported!')
  }
  return(solverdata)
}
