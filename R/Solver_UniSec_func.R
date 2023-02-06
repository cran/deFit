#' Core code of univariable second-order differential equation
#'
#' @param data Users' data
#' @param model model's class is dataframe.
#' @param guess Guess values that contain coefficient and initial values.
#' @param method "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN" and "Brent"
#'
#' @return The result of optimization, SE, RMSE, r-squared, users' data, predictor data and output table.

Slover_UniSec_func <- function(data,model,guess,method){
  # print('------Begin to estimate the parameters -------')
  message('Program will fit the data with a univariate second-order differential equation.')
  message('The differential equation is:')
  message("x(2) = beta1 * x + beta2 * x(1)")
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
    guess = c(0,0,userdata[1,userdata_field[1]],0,0,0)
    # message(paste("The optimization method is",method,"."))
    # message("The initial guess values are")
    # message(paste('beta1:',guess[1],', beta2:',guess[2],', init01(x_0):',guess[3],', init02(dx_0):',guess[4],', init03(dx_0):',guess[5],', init04(d2x_0):',guess[6]))
  }else{
    # message('The initial guess values you input are')
    # message(guess)
  }
  user_CalcUniSec_data <- calc_UniSec_func(userdata,guess,method=method)
  Predictor_UniSec_data <- Predictor_UniSec_func(userdata,user_CalcUniSec_data)
  squareR_UniSec_data <- SquareR_UniSec_func(userdata,Predictor_UniSec_data)
  rmse_UniSec_data <- RMSE_UniSec_func(userdata,Predictor_UniSec_data)
  SE_UniSec_data <- Hessian_UniSec_func(user_CalcUniSec_data$hessian)
  outputDE1 <- paste(userdata_field[1],'(2)=',user_CalcUniSec_data$par[1],'*',userdata_field[1],'+',user_CalcUniSec_data$par[2],' * ',userdata_field[1],'(1)')
  outputDE2 <- paste('Init t0:',user_CalcUniSec_data$par[3])
  outtable <- data.frame()
  outtable[seq(1,6,by=1),'parameter'] = c(paste(userdata_field[1],'(0) to ',userdata_field[1],'(2)',sep = ""),
                                         paste(userdata_field[1],'(1) to ',userdata_field[1],'(2)',sep = ""),
                                         'init01','init02','init03','init04')
  outtable['value'] = user_CalcUniSec_data$par
  outtable['SE'] = c(SE_UniSec_data[1,1],SE_UniSec_data[2,2],SE_UniSec_data[3,3],SE_UniSec_data[4,4],SE_UniSec_data[5,5],SE_UniSec_data[6,6])
  IsConvergence = user_CalcUniSec_data$convergence
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
              Parameter = user_CalcUniSec_data,
              DifferentialEquational = c(outputDE1,outputDE2),
              Predictor = Predictor_UniSec_data,
              Rsquared = squareR_UniSec_data,
              Rmse = rmse_UniSec_data,
              SE = SE_UniSec_data,
              table = outtable,
              IsConvergence = IsConvergence))
}

solve_UniSec_func <- function(t, y, dy, parms){
  res1 <- dy[1] - y[2]
  res2 <- dy[2] - parms[1] * y[2] - parms[2] * y[1]
  list(c(res1, res2))
}

min_UniSec_func <- function(userdata,x0){
  times <- userdata[,'time']
  #remove 'seq' and 'time'
  user_data_field <- names(userdata)
  delSeq <- which(user_data_field == 'seq')
  user_data_field = user_data_field[-delSeq]
  deltime <- which(user_data_field == 'time')
  user_data_field = user_data_field[-deltime]
  min_parms <- c(a = x0[1],
                 b = x0[2])
  y0 <- c(y01 = x0[3],
             y02 = x0[4])
  dy0 <- c(dy01 = x0[5],
           dy02 = x0[6])
  mini_UniSec_data <- deSolve::daspk(y=y0,
                                    dy = dy0,
                                     times=times,
                                     res=solve_UniSec_func,
                                     parms=min_parms)
  # print(mini_UniSec_data)
  mid_res = sum((mini_UniSec_data[,2] - userdata[user_data_field[1]]) **2 )
  return(mid_res)
}
calc_UniSec_func <- function(userdata,guessdata,method){
  message('Finishing optimization...')
  mid_calc_BinFirst_data = stats::optim(c(guessdata[1],
                                          guessdata[2],
                                          guessdata[3],
                                          guessdata[4],
                                          guessdata[5],
                                          guessdata[6]),
                                        min_UniSec_func,
                                        NULL,
                                        method = method,
                                        userdata = userdata,
                                        hessian=TRUE)
  return(mid_calc_BinFirst_data)
}
Predictor_UniSec_func <- function(userdata,calcdata){
  mid_times = userdata[,'time']
  mid_a = calcdata['par'][[1]][1]
  mid_b = calcdata['par'][[1]][2]
  mid_initY01 = calcdata['par'][[1]][3]
  mid_initY02 = calcdata['par'][[1]][4]
  mid_initDY01 = calcdata['par'][[1]][5]
  mid_initDY02 = calcdata['par'][[1]][6]
  mid_usersol_data <- deSolve::daspk(y=c(y1=mid_initY01,y2=mid_initY02),
                                     dy = c(dy1=mid_initDY01,dy2=mid_initDY02),
                                     times=mid_times,
                                     res=solve_UniSec_func,
                                     parms=c(mid_a,mid_b))
  #remove 'seq' and 'time
  user_data_field <- names(userdata)
  delSeq <- which(user_data_field == 'seq')
  user_data_field = user_data_field[-delSeq]
  deltime <- which(user_data_field == 'time')
  user_data_field = user_data_field[-deltime]
  # Change the field(columns) to "solver_***"
  # because of the equational is Univariable Second order differential the values belong of second column, so the value is 2.
  mid_usersol_data = data.frame(mid_usersol_data)
  names(mid_usersol_data)[2] <- paste('solver_',user_data_field[1],sep = '')
  return(mid_usersol_data[,1:2])
}
SquareR_UniSec_func <- function(userdata,Predictor_data){
  message('Estimating R_squared')
  userDF <- merge(userdata,Predictor_data,by='time')
  # SSxe <- sum((userDF[,3] - userDF[,4])**2)
  # SSx <- sum((userDF[,3] - mean(userDF[,3]))**2)
  # #Square_x <- SSxe / SSx
  # SquareR_ <- SSxe / SSx
  # return(1 - SquareR_)
  r_squared01 = cor(userDF[,3],userDF[,4]) **2
  name1 = paste('r_squared',names(userDF)[3],': ',r_squared01,sep = "")
  return(name1)
}
RMSE_UniSec_func <- function(userdata,Predictor_data){
  message('Estimating RMSE')
  userDF <- merge(userdata,Predictor_data,by='time')
  SSxe <- sum((userDF[,3] - userDF[,4])**2)
  RMSEx <- sqrt(SSxe / NROW(userDF[,'time']))
  return(RMSEx)
}
Hessian_UniSec_func <- function(hessian){
  message('Estimating Hessian')
  hessian[hessian<0] <- NaN
  mid_SE <- sqrt(1/hessian)
  return(mid_SE)
}
