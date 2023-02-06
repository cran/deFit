#' Core code of Binary first-order differential equational
#'
#' @param data Users' data
#' @param model model's class is dataframe.
#' @param guess Guess values that contain coefficient and initial values.
#' @param method "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN" and "Brent"
#' @importFrom stats cor
#'
#' @return The result of optimization,SE,RMSE,r-squared,users's data,predictor data and output table.

Slover_BinFirst_func <- function(data,model,guess,method){
  # print('------Begin to estimate the parameters -------')
  message('Program will fit the data with a bivariate first-order differential equation.')
  message('The differential equations are:')
  message('dx/dt = beta1 * x + beta2 * y')
  message('dy/dt = beta3 * x + beta4 * y')
  message('Optimizing...')
  userdata = data #Suggest the time step by 0.1
  if (all(is.na(guess))){
    # If you did not guess the number, the program will guess parms from 0 . And the yini will guess from the data you have upload
    userdata_field <- unlist(c(names(userdata)))
    #remove 'seq' and 'time'
    delSeq <- which(userdata_field == 'seq')
    userdata_field = userdata_field[-delSeq]
    deltime <- which(userdata_field == 'time')
    userdata_field = userdata_field[-deltime]
    # We can use the average of documents' values in head 5
    guess = c(0,0,0,0,userdata[1,userdata_field[1]],userdata[1,userdata_field[2]])
    # message(paste("The optimization method is",method,"."))
    # print("The initial guess values are")
    # print(paste('beta1:',guess[1],', beta2:',guess[2],', beta3:',guess[3],', beta4:',guess[4],', x0:',guess[5],', y0:',guess[6]))
  }else{
    # message('The initial guess values you input are')
    # message(guess)
  }
  user_calc_data <- calc_BinFirst_func(userdata,guess,method=method)
  Predictor_data <- Predictor_BinFirst_func(userdata,user_calc_data)
  squareR_data <- SquareR_BinFirst_func(userdata,Predictor_data)
  user_rmse_data <- RMSE_BinFirst_func(userdata,Predictor_data)
  SE_data <- Hessian_BinFirst_func(user_calc_data$hessian)
  outputDE1 <- paste(userdata_field[1],'(1)=',user_calc_data$par[1],'*',userdata_field[1],'+',user_calc_data$par[2],' * ',userdata_field[2],sep = '')
  outputDE2 <- paste(userdata_field[2],'(1)=',user_calc_data$par[3],'*',userdata_field[1],'+',user_calc_data$par[4],' * ',userdata_field[2],sep = '')
  outputDE3 <- paste('Init t0_x:',user_calc_data$par[5],', Init t0_y:',user_calc_data$par[6])
  outtable <- data.frame()
  outtable[seq(1,6,by=1),'parameter'] = c(paste(userdata_field[1],'(0) to ',userdata_field[1],'(1)',sep = ""),
                                         paste(userdata_field[2],'(0) to ',userdata_field[1],'(1)',sep = ""),
                                         paste(userdata_field[1],'(0) to ',userdata_field[2],'(1)',sep = ""),
                                         paste(userdata_field[2],'(0) to ',userdata_field[2],'(1)',sep = ""),
                                         'init01',
                                         'init02')
  outtable[,'value'] = user_calc_data$par
  outtable[,'SE'] = c(SE_data[1,1],
                      SE_data[2,2],
                      SE_data[3,3],
                      SE_data[4,4],
                      SE_data[5,5],
                      SE_data[6,6])
  IsConvergence = user_calc_data$convergence
  if(IsConvergence == 0){
    IsConvergence = 'successful:successful completion:(which is always the case for SANN and Brent)'
  }else if(IsConvergence == 1){
    IsConvergence = 'successful:indicates that the iteration limit maxit had been reached.'
  }else if(IsConvergence == 10){
    IsConvergence = 'indicates degeneracy of the Nelder-Mead simplex.")'
  }else if(IsConvergence == 51){
    IsConvergence = 'indicates a warning from the "L-BFGS-B" method; see component message for further details.'
  }else if(IsConvergence == 52){
    IsConvergence = 'indicates an error from the "L-BFGS-B" method; see component message for further details.'
  }else{
    IsConvergence = 'Something error")'
  }
  return(list(userdata = userdata,
              Parameter = user_calc_data,
              DifferentialEquational = c(outputDE1,outputDE2,outputDE3),
              Predictor = Predictor_data,
              Rsquared = squareR_data,
              Rmse = user_rmse_data,
              SE = SE_data,
              table = outtable,
              IsConvergence = IsConvergence))
}


solve_BinFirst_func <- function(t,state,parms){
  with(as.list(c(state,parms)),{
    dX <- a*x + b*y
    dY <- c*x + d*y
    list(c(dX,dY))
  })# end with(as.list ...
}

min_BinFirst_func <- function(x0,userdata){
  times <- userdata[,'time']
  #remove 'seq' and 'time'
  user_data_field <- names(userdata)
  delSeq <- which(user_data_field == 'seq')
  user_data_field = user_data_field[-delSeq]
  deltime <- which(user_data_field == 'time')
  user_data_field = user_data_field[-deltime]
  min_parms <- c(a = x0[1],
                 b = x0[2],
                 c = x0[3],
                 d = x0[4])
  state <- c(x = x0[5],
             y = x0[6])
  mini_BinFirst_data <- deSolve::ode(y=state,
                                     times=times,
                                     func=solve_BinFirst_func,
                                     parms=min_parms)
  mid_res = sum((mini_BinFirst_data[,'x'] - userdata[user_data_field[1]]) **2 )
  mid_res2 = sum((mini_BinFirst_data[,'y'] - userdata[user_data_field[2]])**2 )
  res = mid_res + mid_res2
  return(res)
}

calc_BinFirst_func <- function(userdata,guessdata,method){
  message('Finishing optimization...')
  mid_calc_BinFirst_data = stats::optim(c(guessdata[1],
                                          guessdata[2],
                                          guessdata[3],
                                          guessdata[4],
                                          guessdata[5],
                                          guessdata[6]),
                                        min_BinFirst_func,
                                        NULL,
                                        method = method,
                                        userdata = userdata,
                                        hessian=TRUE)
  return(mid_calc_BinFirst_data)

}
Predictor_BinFirst_func <- function(userdata,calcdata){
  mid_times = userdata[,'time']
  mid_a = calcdata['par'][[1]][1]
  mid_b = calcdata['par'][[1]][2]
  mid_c = calcdata['par'][[1]][3]
  mid_d = calcdata['par'][[1]][4]
  mid_initX = calcdata['par'][[1]][5]
  mid_initY = calcdata['par'][[1]][6]
  mid_parms <- c(a = mid_a,
                 b = mid_b,
                 c = mid_c,
                 d = mid_d)
  mid_state <- c(x = mid_initX,
                 y = mid_initY)
  mid_usersol_data <- deSolve::ode(y=mid_state,
                                   times=mid_times,
                                   func=solve_BinFirst_func,
                                   parms=mid_parms)
  #remove 'seq' and 'time
  user_data_field <- names(userdata)
  delSeq <- which(user_data_field == 'seq')
  user_data_field = user_data_field[-delSeq]
  deltime <- which(user_data_field == 'time')
  user_data_field = user_data_field[-deltime]
  # Change the field to "solver_***"
  PredictorColumns <- unlist(c(names(userdata)))
  solverNum1 <- which(PredictorColumns == user_data_field[1])
  mid_usersol_data = data.frame(mid_usersol_data)
  names(mid_usersol_data)[2] <- paste('solver_',user_data_field[1],sep = '')
  solverNum2 <- which(PredictorColumns == user_data_field[2])
  names(mid_usersol_data)[3] <- paste('solver_',user_data_field[2],sep = '')
  return(mid_usersol_data)

}
SquareR_BinFirst_func <- function(userdata,Predictor_data){
  #
  userDF <- merge(userdata,Predictor_data,by='time')
  # SSxe <- sum((userDF[,3] - userDF[,5])**2)
  # SSx <- sum((userDF[,3] - mean(userDF[,3]))**2)
  # #Square_x <- SSxe / SSx
  # SSye <- sum((userDF[,4] - userDF[,6])**2)
  # SSy <- sum((userDF[,4] - mean(userDF[,4]))**2)
  # #Square_y <- SSye / SSy
  # SquareR_ <- (SSxe + SSye) / (SSx + SSy)
  # return(1 - SquareR_)

  r_squared01 = stats::cor(userDF[,3],userDF[,5]) **2
  r_squared02 = stats::cor(userDF[,4],userDF[,6]) **2
  name1 = paste('r_squared',names(userDF)[3],': ',r_squared01,sep = "")
  name2 = paste('r_squared',names(userDF)[4],': ',r_squared02,sep = "")
  message('Estimating R_squared')
  return(c(name1,name2))
}

RMSE_BinFirst_func <- function(userdata,Predictor_data){
  message('Estimating RMSE')
  userDF <- merge(userdata,Predictor_data,by='time')
  SSxe <- sum((userDF[,3] - userDF[,5])**2)
  RMSEx <- sqrt(SSxe / NROW(userDF[,'time']))
  SSye <- sum((userDF[,4] - userDF[,6])**2)
  RMSEy <- sqrt(SSye / NROW(userDF[,'time']))
  return(c(RMSEx,RMSEy))
}

Hessian_BinFirst_func <- function(hessian){
  message('Estimating Hessian')
  hessian[hessian<0] <- NaN
  mid_SE <- sqrt(1/hessian)
  return(mid_SE)
}
