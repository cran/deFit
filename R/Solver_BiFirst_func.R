
#' Solver of bivariate first-order differential equation
#'
#' @param userdata a data frame containing all model variables. The "time" column must be included.
#' @param var_model a dataframe containing equations.
#' @param guess a list or a string. Guess the coefficients or initial values.
#' @param method a list or a string. The available options are 'Nelder-Mead','L-BFGS-B','SANN' and 'BFGS'.
#'
#' @return a list
Solver_BiFirst_func <- function(userdata,var_model,guess,method){
  message('Program will fit the data with a bivariate first-order differential equation.')
  message('The differential equations are:')
  message('dx/dt = beta1 * x + beta2 * y')
  message('dy/dt = beta3 * x + beta4 * y')
  message('Optimizing...')

  ## identify variables and information
  var_notime_model = var_model[var_model$field != 'time','field']
  var_notime_hat = paste(var_notime_model,'_hat',sep = "")
  var_notime_data = var_model[var_model$field != 'time','variable']
  var_model_new = var_model[,'field']
  var_userdata = var_model[,'variable']

  ## add the columns of model in userdata
  for (i in 1:length(var_userdata)) {
    userdata[,var_model_new[i]] = userdata[,var_userdata[i]]
  }

  calc_data = calc_BiFirst_func(userdata,var_model,guess,method)
  Neg2LL = calc_data[3]
  predict_data = calc_data[2][[1]]
  colnames(predict_data) = c('time',c(var_notime_hat))
  ## equation
  equation1 = paste(var_notime_model[1],"(1) = ",calc_data[1][[1]]$par[1]," * ",var_notime_model[1]," + ",calc_data[1][[1]]$par[2]," * ",var_notime_model[2],sep = "")
  equation2 = paste(var_notime_model[2],"(1) = ",calc_data[1][[1]]$par[3]," * ",var_notime_model[1]," + ",calc_data[1][[1]]$par[4]," * ",var_notime_model[2],sep = "")
  equation3 = paste("Init t0_",var_notime_model[1],": ",calc_data[1][[1]]$par[5],", Init t0_",var_notime_model[2],": ",calc_data[1][[1]]$par[6],sep = "")
  equation = list(equation1,equation2,equation3)
  ## table
  table <- data.frame()
  table[seq(1,6,by=1),'parameter'] = c(paste(var_notime_model[1],"(0) to ",var_notime_model[1],"(1)",sep = ""),
                                       paste(var_notime_model[2],"(0) to ",var_notime_model[1],"(1)",sep = ""),
                                       paste(var_notime_model[1],"(0) to ",var_notime_model[2],"(1)",sep = ""),
                                       paste(var_notime_model[2],"(0) to ",var_notime_model[2],"(1)",sep = ""),
                                       paste("init t0_",var_notime_model[1],sep = ""),
                                       paste("init t0_",var_notime_model[2],sep = ""))
  table[,'value'] = calc_data[1][[1]]$par[1:6]
  return(list(solve_data=calc_data[1][[1]],
              userdata=userdata,
              predict_data=predict_data,
              predict_fixed_data=NULL,
              table=table,
              equation=equation,
              Neg2LL=Neg2LL))
  }
calc_BiFirst_func <- function(userdata,var_model,guess,method){
  ## identify variables and information
  var_notime_model = var_model[var_model$field != 'time','field']
  var_notime_hat = paste(var_notime_model,'_hat',sep = "")
  var_notime_data = var_model[var_model$field != 'time','variable']
  guess_e = c(guess[[1]],c(0.1,0.1))

  mid_BiFirst_data = stats::optim(guess_e,
                                  min_BiFirst_func,
                                  NULL,
                                  method = method,
                                  userdata=userdata,
                                  var_notime_model=var_notime_model,
                                  hessian = TRUE)
  message('Finishing optimization...')
  Neg2LL = mid_BiFirst_data$value[[1]][1]  + 3.67
  min_parms <- c(beta1 = mid_BiFirst_data$par[1],
                 beta2 = mid_BiFirst_data$par[2],
                 beta3 = mid_BiFirst_data$par[3],
                 beta4 = mid_BiFirst_data$par[4])
  state <- c(x = mid_BiFirst_data$par[5],
             y = mid_BiFirst_data$par[6])
  times <- seq(min(userdata$time),max(userdata$time),length.out=length(userdata$time))
  predict_data <- deSolve::ode(y=state,
                               times=times,
                               func=solve_BiFirst_func,
                               parms=min_parms)
  return(list(mid_BiFirst_data,predict_data,Neg2LL))
}
min_BiFirst_func <- function(x0,userdata,var_notime_model){
  times <- seq(min(userdata$time),max(userdata$time),length.out=length(userdata$time))
  min_parms <- c(beta1 = x0[1],
                 beta2 = x0[2],
                 beta3 = x0[3],
                 beta4 = x0[4])
  state <- c(x = x0[5],
             y = x0[6])
  mid_solve_data <- deSolve::ode(y=state,
                                 times=times,
                                 func=solve_BiFirst_func,
                                 parms=min_parms)
  mid_solve_data = data.frame(mid_solve_data)
  # estimating residual in loss function
  guess_xy = c(x0[7],x0[8])
  if(x0[7] == 0){
    guess_xy[1] = 0.01
  }
  if(x0[8] == 0){
    guess_xy[2] = 0.01
  }
  log_e_x = log(guess_xy[1]**2)
  log_e_y = log(guess_xy[2]**2)
  inverse_ex = 1/(guess_xy[1]**2)
  inverse_ey = 1/(guess_xy[2]**2)
  #mid_res = sum((mid_solve_data[,'x'] - userdata[var_notime_model[1]]) **2 )
  #mid_res2 = sum((mid_solve_data[,'y'] - userdata[var_notime_model[2]])**2 )
  mid_res = sum(log_e_x + ((mid_solve_data[,'x'] - userdata[var_notime_model[1]]) **2) * inverse_ex )
  mid_res2 = sum(log_e_y + ((mid_solve_data[,'y'] - userdata[var_notime_model[2]]) **2) * inverse_ey )
  res = mid_res + mid_res2

  if(is.nan(res) || is.infinite(res) ){
    mid_res = sum(userdata[var_notime_model[1]]  **2)
    mid_res2 = sum(userdata[var_notime_model[2]] **2)
    res = mid_res + mid_res2
  }
  return(res)
}
solve_BiFirst_func <- function(t,state,parms){
  with(as.list(c(state,parms)),{
    dX <- beta1*x + beta2*y
    dY <- beta3*x + beta4*y
    list(c(dX,dY))
  })# end with(as.list ...
}
