#' Core code of univariable second-order differential equation
#'
#' @param data User's data
#' @param model model's class is dataframe.
#' @param guess Guess values that contain coefficient and initial values.
#' @param method "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN" and "Brent"
#'
#' @return The result of optimization, SE, RMSE, r-squared, users' data, predictor data and output table.

Slover_UniSec_func <- function(data,model,guess,method){
  var_model = model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]

  # print('------Begin to estimate the parameters -------')
  message('Program will fit the data with a univariate second-order differential equation.')
  message('The differential equation is:')
  message("x(2) = beta1 * x + beta2 * x(1)")
  message('Optimizing...')
  userdata = data #Suggest the time step by 0.1
  # If you did not guess the number, the program will guess parms from 0 . And the yini will guess from the data you have upload
  userdata_field <- unlist(c(names(userdata)))
  #remove 'seq' and 'time'
  delSeq <- which(userdata_field == 'seq')
  userdata_field = userdata_field[-delSeq]
  deltime <- which(userdata_field == 'time')
  userdata_field = userdata_field[-deltime]
  if (all(is.na(guess))){
    # We can use the average of documents' values in head 5
    guess = c(0,0,userdata[1,userdata_field[1]],0,0,0)
    # message(paste("The optimization method is",method,"."))
    # message("The initial guess values are")
    # message(paste('beta1:',guess[1],', beta2:',guess[2],', init01(x_0):',guess[3],', init02(dx_0):',guess[4],', init03(dx_0):',guess[5],', init04(d2x_0):',guess[6]))
  }else{
    # message('The initial guess values you input are')
    if(length(guess) == 6){
      guess = guess
    }else{
      stop('Your model is a univariate second-order differential equation that the guess values need 6 numbers. like c(0,0,0,0,0,0)')
    }

    # message(guess)
  }
  user_CalcUniSec_data <- calc_UniSec_func(userdata,guess,method=method,model=model)
  Predictor_UniSec_data <- Predictor_UniSec_func(user_CalcUniSec_data$detailtable,model)
  squareR_UniSec_data <- SquareR_UniSec_func(Predictor_UniSec_data,model)
  rmse_UniSec_data <- RMSE_UniSec_func(Predictor_UniSec_data,model)
  rmse_UniSec_data <- paste('RMSE_',var_model[1,'field'],' = ',rmse_UniSec_data,sep='')
  SE_UniSec_data <- Hessian_UniSec_func(user_CalcUniSec_data$calc_res$hessian)
  outputDE1 <- paste(userdata_field[1],'(2)=',user_CalcUniSec_data$calc_res$par[1],'*',userdata_field[1],'+',user_CalcUniSec_data$calc_res$par[2],'*',userdata_field[1],'(1)')
  outputDE2 <- paste('Init t0:',user_CalcUniSec_data$calc_res$par[3])
  outtable <- data.frame()
  outtable[seq(1,3,by=1),'parameter'] = c(paste(userdata_field[1],'(0) to ',userdata_field[1],'(2)',sep = ""),
                                         paste(userdata_field[1],'(1) to ',userdata_field[1],'(2)',sep = ""),
                                         paste('init_',var_model[1,'field'],sep='')) #,'init02','init03','init04')
  outtable['value'] = user_CalcUniSec_data$calc_res$par[1:3]
  outtable['SE'] = c(SE_UniSec_data[1,1],SE_UniSec_data[2,2],SE_UniSec_data[3,3])#,SE_UniSec_data[4,4],SE_UniSec_data[5,5],SE_UniSec_data[6,6])
  IsConvergence = user_CalcUniSec_data$calc_res$convergence
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
              Parameter = user_CalcUniSec_data$calc_res,
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
calc_UniSec_func <- function(userdata,guessdata,method,model){
  message('Finishing optimization...')

  # ###the variable of model
  var_model = model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]

  mid_detail_table = data.frame(seq=1:nrow(userdata))
  mid_detail_table['time'] = userdata[,'time']
  mid_detail_table[var_model[1,'field']] = userdata[,var_model[1,'field']]
  mytry <- tryCatch({
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
  },
  warning = function(war){
    message('Waring @ ')
  },
  error = function(err){
    message('Error @  ',err,'You can try method ="Nelder-Mead"')
  }
  # ,finally = {
  #   message('next...subject...')
  # }
  )
  mid_detail_table['beta1'] = mid_calc_BinFirst_data['par'][[1]][1]
  mid_detail_table['beta2'] = mid_calc_BinFirst_data['par'][[1]][2]
  mid_detail_table['init01'] = mid_calc_BinFirst_data['par'][[1]][3]
  mid_detail_table['init02'] = mid_calc_BinFirst_data['par'][[1]][4]
  mid_detail_table['init03'] = mid_calc_BinFirst_data['par'][[1]][5]
  mid_detail_table['init04'] = mid_calc_BinFirst_data['par'][[1]][6]
  # print(mid_detail_table)
  mid_detail_table[is.na(mid_detail_table)] <- 0

  return(list(detailtable=mid_detail_table,calc_res=mid_calc_BinFirst_data))
}
Predictor_UniSec_func <- function(detailtable,model){
  mid_detailtable = detailtable

  var_model = model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]

  detailtable[paste('solver_',var_model[1,'field'],sep='')]= NA

  times = sort(unique(mid_detailtable[,'time']),decreasing = FALSE)
  mid_a = mean(mid_detailtable[,'beta1'])
  mid_b = mean(mid_detailtable[,'beta2'])
  mid_initY01 = mean(mid_detailtable[,'init01'])
  mid_initY02 = mean(mid_detailtable[,'init02'])
  mid_initDY01 = mean(mid_detailtable[,'init03'])
  mid_initDY02 = mean(mid_detailtable[,'init04'])
  mid_usersol_data <- deSolve::daspk(y=c(y1=mid_initY01,y2=mid_initY02),
                                     dy = c(dy1=mid_initDY01,dy2=mid_initDY02),
                                     times=times,
                                     res=solve_UniSec_func,
                                     parms=c(mid_a,mid_b))
  #remove 'seq' and 'time
  mid_usersol_data <- data.frame(mid_usersol_data)
  detailtable[,paste('solver_',var_model[1,'field'],sep='')] = mid_usersol_data[mid_usersol_data['time'] == detailtable[,'time'],'y1']
  return(detailtable)
}
SquareR_UniSec_func <- function(Predictor_data,model){
  message('Estimating R_squared')

  var_model = model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]

  # SSxe <- sum((userDF[,3] - userDF[,4])**2)
  # SSx <- sum((userDF[,3] - mean(userDF[,3]))**2)
  # #Square_x <- SSxe / SSx
  # SquareR_ <- SSxe / SSx
  # return(1 - SquareR_)
  r_squared01 = stats::cor(Predictor_data[,var_model[1,'field']],Predictor_data[,paste('solver_',var_model[1,'field'],sep='')]) **2
  name1 = paste('r_squared_',var_model[1,'field'],' = ',r_squared01,sep = "")
  return(name1)
}
RMSE_UniSec_func <- function(Predictor_data,model){
  message('Estimating RMSE')

  var_model = model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]
  SSxe <- sum((Predictor_data[,var_model[1,'field']] - Predictor_data[,paste('solver_',var_model[1,'field'],sep = '')])**2)
  RMSEx <- sqrt(SSxe / NROW(Predictor_data[,'time']))
  return(RMSEx)
}
Hessian_UniSec_func <- function(hessian){
  message('Estimating Hessian')
  hessian[hessian<0] <- NaN
  mid_SE <- sqrt(1/hessian)
  return(mid_SE)
}
