#' Solver of univariate second-order differential equation
#'
#' @param userdata a data frame containing all model variables. The "time" column must be included.
#' @param var_model a dataframe containing equations
#' @param guess a list or a string. Guess the coefficients or initial values.
#' @param method a list or a string. The available options are 'Nelder-Mead','L-BFGS-B','SANN' and 'BFGS'.
#'
#' @return a list
Solver_UniSecond_func <- function(userdata,var_model,guess,method){
  message('Program will fit the data with a univariate second-order differential equation.')
  message('The differential equation is:')
  message("x(2) = beta1 * x + beta2 * x(1)")
  message('Optimizing...')

  ## identify variables and information
  names(userdata)[names(userdata) == var_model[var_model$field == 'time','variable']] <- 'time'
  var_notime_model = var_model[var_model$field != 'time','field']
  var_notime_hat = paste(var_notime_model,'_hat',sep = "")
  var_notime_data = var_model[var_model$field != 'time','variable']
  var_model_new = var_model[,'field']
  var_userdata = var_model[,'variable']

  ## add the columns of model in userdata
  userdata[,var_model_new[1]] = userdata[,var_userdata[1]]

  calc_data = calc_UniSec_func(userdata,guess,method,var_notime_model)
  predict_data = calc_data[2][[1]]
  colnames(predict_data)[2] = var_notime_hat
  ## equation
  equation1 = paste(var_notime_model[1],"(2) = ",calc_data[1][[1]]$par[1]," * ",var_notime_model[1]," + ",calc_data[1][[1]]$par[2]," * ",var_notime_model[1],"(1)",sep = "")
  equation2 = paste("Init t0_",var_notime_model[1],": ",calc_data[1][[1]]$par[3],sep = "")
  equation = list(equation1,equation2)
  ## table
  table <- data.frame()
  table[seq(1,4,by=1),'parameter'] = c(paste(var_notime_model[1],"(0) to ",var_notime_model[1],"(2)",sep = ""),
                                       paste(var_notime_model[1],"(1) to ",var_notime_model[1],"(2)",sep = ""),
                                       paste("init t0_",var_notime_model[1],sep = ""),
                                       paste("init t0_d",var_notime_model[1],sep = ""))
  table[,'value'] = c(calc_data[1][[1]]$par[1],
                      calc_data[1][[1]]$par[2],
                      calc_data[1][[1]]$par[3],
                      calc_data[1][[1]]$par[4])
  return(list(solve_data = calc_data[1][[1]],
              userdata=userdata,
              predict_data=predict_data,
              predict_fixed_data=NULL,
              table=table,
              equation=equation))
}

calc_UniSec_func <- function(userdata,guess,method,var_notime_model){
  mid_UniSec_data = stats::optim(guess[[1]],
                                  min_UniSec_func,
                                  NULL,
                                  method = method,
                                  userdata=userdata,
                                  var_notime_model=var_notime_model,
                                  hessian = TRUE)
  message('Finishing optimization...')

  ## calculate predict data
  times <- seq(min(userdata$time),max(userdata$time),length.out=length(userdata$time))
  initial_conditions <- c(y1 = mid_UniSec_data$par[3], y2 = mid_UniSec_data$par[4])  # initial state,initial speed
  parameters <- c(beta1=mid_UniSec_data$par[1],beta2=mid_UniSec_data$par[2])  # beta1,beta2
  predict_data <- deSolve::ode(y = initial_conditions,
                               times = times,
                               func = damped_oscillator,
                               parms = parameters)
  return(list(mid_UniSec_data,predict_data))

}
min_UniSec_func <- function(x0,userdata,var_notime_model){
  times <- seq(min(userdata$time),max(userdata$time),length.out=length(userdata$time))
  initial_conditions <- c(y1 = x0[3], y2 = x0[4])  # initial state,initial speed
  parameters <- c(beta1=x0[1],beta2=x0[2])  # beta1,beta2
  mid_solve_data <- deSolve::ode(y = initial_conditions,
                                 times = times,
                                 func = damped_oscillator,
                                 parms = parameters)
  res = sum((mid_solve_data[,'y1'] - userdata[,var_notime_model[1]]) **2 )

  ## hhhh so happy
  if(!(is.finite(res))){
    mid_merge_df = merge(userdata,mid_solve_data, by = "time", all.x = TRUE)
    cat('waiting...\n')
    mid_merge_df[is.na(mid_merge_df)] = 0
    res = sum(mid_merge_df[,'y1'] - mid_merge_df[,var_notime_model[1]] **2)
  }
  return(res)
}
damped_oscillator <- function(t, y, parms) {
  with(as.list(parms), {
    dy1 <- y[2]
    dy2 <- beta1 * y[1] + beta2 * y[2]
    return(list(c(dy1, dy2)))
  })
}
