
#' Title
#'
#' @param data User's data
#' @param model Model's class is dataframe.
#' @param guess Guess values that contain coefficient and initial values.
#' @param method "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN" and "Brent"
#' @param guess2 Guess values of multilevel that contain coefficient and initial values.
#'
#' @return The result of optimization,SE,RMSE,r-squared,users's data,predictor data and output table.

Solver_MultiBiFirst_func <- function(data,model,guess,method,guess2){
  # cat('### Program will fit the data with multilevel bivariate first-order differential equations')
  # print('### Firstly, we should estimate the population parameter of differential equations')
  # ###the variable of model
  var_model = model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]
  sub_model = model
  sub_model = sub_model[which(sub_model['operator'] == '~'),]
  if(all(is.na(guess))){
    guess = c(0,0,0,0,data[1,var_model[1,'field']],data[1,var_model[2,'field']])
  }
  x0 = guess
  cat('Estimating parameters \n')
  cat('Waiting...\n')
  if(method == 'bayesian'){
    min_population_data <- bayesian_population_func(userdata = data,
                                                    model = model,
                                                    x0 = x0)
  }else{
    min_population_data <- stats::optim(x0,
                                        min_MultiPopulation_func,
                                        NULL,
                                        method = method,
                                        userdata = data,
                                        model = model,
                                        hessian=TRUE
                                        )
  }
  # print(min_population_data)
  fixvalues = min_population_data$par
  detailtable = min_multiBiFirst_IndiOnebyOne_func(guess2,
                                        userdata=data,
                                        model=model,
                                        fixvalues=fixvalues,
                                        method = method
                                        )
  outputDE0 <- 'Fixed effects:'
  outputDE1 <- paste(var_model[1,'field'],'(1)=',min_population_data$par[1],'*',var_model[1,'field'],'+',min_population_data$par[2],' * ',var_model[2,'field'],sep = '')
  outputDE2 <- paste(var_model[2,'field'],'(1)=',min_population_data$par[3],'*',var_model[1,'field'],'+',min_population_data$par[4],' * ',var_model[2,'field'],sep = '')
  outputDE3 <- paste('Init t0_x:',min_population_data$par[5],', Init t0_y:',min_population_data$par[6])

  Predictor_data = Predictor_MultiBiFirst_func(detailtable,model)
  squareR_data = SquareR_MultiBiFirst_func(Predictor_data,model)
  user_rmse_data = RMSE_MultiBiFirst_func(Predictor_data,model)
  user_rmse_data <- paste('RMSE_',var_model[1,'field'],'=',user_rmse_data[1],
                          ' & RMSE_',var_model[2,'field'],'=',user_rmse_data[2],sep = '')
  SE_data = Hessian_MultiBiFirst_func(min_population_data$hessian)

  outputDE4 <- 'Random effects:'
  outputDE5 <- data.frame(seq=1:6,
                          Groups = NA,
                          Name = NA,
                          Variance = NA,
                          Std.Dev. = NA)
  outputDE5[1,'Groups'] = 'Subject'
  outputDE5[1,'Name'] = paste(var_model[1,'field'],'(0) to ',var_model[1,'field'],'(1)',sep='')
  outputDE5[1,'Variance'] = stats::var(Predictor_data[Predictor_data['time'] == Predictor_data[1,'time'],'etaI1'])
  outputDE5[1,'Std.Dev.'] = sqrt(outputDE5[1,'Variance'])
  outputDE5[2,'Groups'] = 'Subject'
  outputDE5[2,'Name'] = paste(var_model[2,'field'],'(0) to ',var_model[1,'field'],'(1)',sep='')
  outputDE5[2,'Variance'] = stats::var(Predictor_data[Predictor_data['time'] == Predictor_data[1,'time'],'etaI2'])
  outputDE5[2,'Std.Dev.'] = sqrt(outputDE5[2,'Variance'])
  outputDE5[3,'Groups'] = 'Subject'
  outputDE5[3,'Name'] = paste(var_model[1,'field'],'(0) to ',var_model[2,'field'],'(1)',sep='')
  outputDE5[3,'Variance'] = stats::var(Predictor_data[Predictor_data['time'] == Predictor_data[1,'time'],'etaI3'])
  outputDE5[3,'Std.Dev.'] = sqrt(outputDE5[3,'Variance'])
  outputDE5[4,'Groups'] = 'Subject'
  outputDE5[4,'Name'] = paste(var_model[2,'field'],'(0) to ',var_model[2,'field'],'(1)',sep='')
  outputDE5[4,'Variance'] = stats::var(Predictor_data[Predictor_data['time'] == Predictor_data[1,'time'],'etaI4'])
  outputDE5[4,'Std.Dev.'] = sqrt(outputDE5[4,'Variance'])
  outputDE5[5,'Groups'] = 'Subject'
  outputDE5[5,'Name'] = paste('init_',var_model[1,'field'],sep='')
  outputDE5[5,'Variance'] = stats::var(Predictor_data[Predictor_data['time'] == Predictor_data[1,'time'],paste('etaInit',var_model[1,'field'],'i1',sep='')])
  outputDE5[5,'Std.Dev.'] = sqrt(outputDE5[5,'Variance'])
  outputDE5[6,'Groups'] = 'Subject'
  outputDE5[6,'Name'] = paste('init_',var_model[2,'field'],sep='')
  outputDE5[6,'Variance'] = stats::var(Predictor_data[Predictor_data['time'] == Predictor_data[1,'time'],paste('etaInit',var_model[2,'field'],'i1',sep='')])
  outputDE5[6,'Std.Dev.'] = sqrt(outputDE5[6,'Variance'])
  outputDE6 <- paste('Number of obs:',nrow(Predictor_data),', groups:',sub_model[1,'subject'],',',length(unique(Predictor_data[,'subject'])))

  outtable <- data.frame()
  outtable[seq(1,6,by=1),'parameter'] = c(paste(var_model[1,'field'],'(0) to ',var_model[1,'field'],'(1)',sep = ""),
                                          paste(var_model[2,'field'],'(0) to ',var_model[1,'field'],'(1)',sep = ""),
                                          paste(var_model[1,'field'],'(0) to ',var_model[2,'field'],'(1)',sep = ""),
                                          paste(var_model[2,'field'],'(0) to ',var_model[2,'field'],'(1)',sep = ""),
                                          paste('init_',var_model[1,'field'],sep=''),
                                          paste('init_',var_model[2,'field'],sep=''))
  outtable[,'value'] = min_population_data$par
  outtable[,'SE'] = c(SE_data[1,1],
                      SE_data[2,2],
                      SE_data[3,3],
                      SE_data[4,4],
                      SE_data[5,5],
                      SE_data[6,6])
  IsConvergence = min_population_data$convergence
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
  return(list(userdata = data,
              Parameter = fixvalues,
              DifferentialEquational = list(outputDE0,outputDE1,outputDE2,outputDE3,outputDE4,outputDE5,outputDE6),
              Predictor = Predictor_data,
              Rsquared = squareR_data,
              Rmse = user_rmse_data,
              SE = SE_data,
              table = outtable,
              IsConvergence = IsConvergence))
  }
min_MultiPopulation_func <- function(x0,userdata,model){
  BetaXToX = x0[1] # beta1
  BetaYToX = x0[2] # beta2
  BetaXToY = x0[3] # beta3
  BetaYToY = x0[4] # beta4
  BetaXt0 = x0[5] # init X
  BetaYt0 = x0[6] # init Y
  # Get the time and sort it in descending order of time.
  # Now we have not considered the time of NA. to-do
  times = sort(unique(userdata[,'time']),decreasing = FALSE)
  mid_deEstimate_data <- deSolve::ode(y=c(x = BetaXt0,
                                          y = BetaYt0),
                                      times=times,
                                      func=solve_BinFirst_func,
                                      parms=c(a = BetaXToX,
                                              b = BetaYToX,
                                              c = BetaXToY,
                                              d = BetaYToY))
  mid_df_desolve <- data.frame(mid_deEstimate_data)
  # Get the variables of the model
  var_model = model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]
  # print(userdata)
  names(mid_df_desolve)[names(mid_df_desolve) == "x"] = 'solver_x'
  names(mid_df_desolve)[names(mid_df_desolve) == "y"] = 'solver_y'
  mid_df_res <- data.frame(seq=1:nrow(userdata),
                           time = userdata[,'time'],
                           x = userdata[,var_model[1,'field']],
                           y = userdata[,var_model[2,'field']],
                           subject = userdata[,'subject'],
                           solver_x = NA,
                           solver_y = NA
                           )
  # When the time is the same, userdata and solverdata will be merged.
  for (i in 1:nrow(mid_df_res)){
    mid_df_res[i,'solver_x'] = mid_df_desolve[mid_df_desolve['time'] == mid_df_res[i,'time'],'solver_x']
    mid_df_res[i,'solver_y'] = mid_df_desolve[mid_df_desolve['time'] == mid_df_res[i,'time'],'solver_y']
  }
  mid_df_res['e_y'] = (mid_df_res['y'] - mid_df_res['solver_y']) ** 2
  mid_df_res['e_x'] = (mid_df_res['x'] - mid_df_res['solver_x']) ** 2
  mid_sum_res = sum(mid_df_res['e_y'] + mid_df_res['e_x'])
  return(mid_sum_res)
}

min_multiBiFirst_IndiOnebyOne_func <- function(x0,userdata,model,fixvalues,method){
  # ###the variable of model
  var_model = model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]
  # ###deal the data
  mid_onebyone_df = data.frame(seq=1:nrow(userdata))
  mid_onebyone_df['time'] = userdata[,'time']
  mid_onebyone_df[var_model[1,'field']] = userdata[,var_model[1,'field']]
  mid_onebyone_df[var_model[2,'field']] = userdata[,var_model[2,'field']]
  mid_onebyone_df['subject'] = userdata[,'subject']
  mid_onebyone_df['beta1'] = fixvalues[1]
  mid_onebyone_df['beta2'] = fixvalues[2]
  mid_onebyone_df['beta3'] = fixvalues[3]
  mid_onebyone_df['beta4'] = fixvalues[4]
  mid_onebyone_df['initXt0'] = fixvalues[5]
  mid_onebyone_df['initYt0'] = fixvalues[6]
  # ###calculate/estimate the random effect
  # Identify the fix random variable
  mid_Identify_fixrand = model
  mid_Identify_fixrand = mid_Identify_fixrand[which(mid_Identify_fixrand['operator'] == '~'),]
  # ###Judge the eta guess values
  etaX0 =rep(0.01,times=6) #Not in use.
  for(i in unique(userdata[,'subject'])){
    mytry <- tryCatch({
      mid_individual_data <- stats::optim(x0,
                                          min_onebyone_BiFirst_Eta_func,
                                          NULL,
                                          method = method,
                                          userdata= userdata[userdata['subject']==i,],
                                          model = model,
                                          fixvalues = fixvalues)
      mid_coef = mid_individual_data$par
    },
    warning = function(war){
      message('Waring @ ')
      mid_coef = c(0,0,0,0,0,0)
      return(mid_coef)
    },
    error = function(err){
      # message('Error @  ',err)
      mid_coef = c(0,0,0,0,0,0)
      return(mid_coef)
    },
    finally = {
      message('next...subject...',i)
    })
    if(!(exists('mid_coef'))){
      mid_coef = c(0,0,0,0,0,0)
    }
    fixRand_model = model
    fixRand_model = fixRand_model[which(fixRand_model['operator'] == '~'),]
    if(!(grepl('1',fixRand_model[1,'fixRand']))){
      mid_coef[5] = 0
    }
    if(!(grepl('1',fixRand_model[2,'fixRand']))){
      mid_coef[6] = 0
    }
    if(!(grepl(var_model[1,'field'],fixRand_model[1,'fixRand']))){
      mid_coef[1] = 0
    }
    if(!(grepl(var_model[2,'field'],fixRand_model[1,'fixRand']))){
      mid_coef[2] = 0
    }
    if(!(grepl(var_model[1,'field'],fixRand_model[2,'fixRand']))){
      mid_coef[3] = 0
    }
    if(!(grepl(var_model[2,'field'],fixRand_model[2,'fixRand']))){
      mid_coef[4] = 0
    }

    mid_onebyone_df[userdata['subject']==i,'etaI1'] = mid_coef[1]
    mid_onebyone_df[userdata['subject']==i,'etaI2'] = mid_coef[2]
    mid_onebyone_df[userdata['subject']==i,'etaI3'] = mid_coef[3]
    mid_onebyone_df[userdata['subject']==i,'etaI4'] = mid_coef[4]
    mid_onebyone_df[userdata['subject']==i,'etaInitXi1'] = mid_coef[5]
    mid_onebyone_df[userdata['subject']==i,'etaInitYi1'] = mid_coef[6]
  }
  return(mid_onebyone_df)
}

min_onebyone_BiFirst_Eta_func <- function(x0,userdata,model,fixvalues){
  EtaXToX = 1 # eta1
  EtaYToX = 2 # eta2
  EtaXToY = 3 # eta3
  EtaYToY = 4 # eta4
  EtaXt0 = 5 # etaXit0
  EtaYt0 = 6 # etayit0

  var_model = model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]

  fixRand_model = model
  fixRand_model = fixRand_model[which(fixRand_model['operator'] == '~'),]
  if(!(grepl('1',fixRand_model[1,'fixRand']))){
    x0[5] = 0
  }
  if(!(grepl('1',fixRand_model[2,'fixRand']))){
    x0[6] = 0
  }
  if(!(grepl(var_model[1,'field'],fixRand_model[1,'fixRand']))){
    x0[1] = 0
  }
  if(!(grepl(var_model[2,'field'],fixRand_model[1,'fixRand']))){
    x0[2] = 0
  }
  if(!(grepl(var_model[1,'field'],fixRand_model[2,'fixRand']))){
    x0[3] = 0
  }
  if(!(grepl(var_model[2,'field'],fixRand_model[2,'fixRand']))){
    x0[4] = 0
  }

  mid_df = data.frame(seq=1:nrow(userdata),
                      time = userdata[,'time'],
                      x = userdata[,var_model[1,'field']],
                      y = userdata[,var_model[2,'field']],
                      subject = userdata[,'subject'],
                      solver_x = NA,
                      solver_y = NA
                      )
  times = sort(unique(userdata[,'time']),decreasing = FALSE)
  mid_deEstimate_data <- deSolve::ode(y=c(x = fixvalues[5] + x0[EtaXt0],
                                          y = fixvalues[6] + x0[EtaYt0]),
                                      times=times,
                                      func=solve_MultiBiFirst_func,
                                      parms=c(a = fixvalues[1] + x0[EtaXToX],
                                              b = fixvalues[2] + x0[EtaYToX],
                                              c = fixvalues[3] + x0[EtaXToY],
                                              d = fixvalues[4] + x0[EtaYToY])
  )
  mid_deEstimate_data = data.frame(mid_deEstimate_data)
  for(i_mid in rownames(mid_df)){
    mid_df[i_mid,'solver_x'] = mid_deEstimate_data[mid_deEstimate_data['time'] == mid_df[i_mid,'time'],'x']
    mid_df[i_mid,'solver_y'] = mid_deEstimate_data[mid_deEstimate_data['time'] == mid_df[i_mid,'time'],'y']
  }
  mid_df['e_y'] = (mid_df['y'] - mid_df['solver_y']) **2
  mid_df['e_x'] = (mid_df['x'] - mid_df['solver_x']) **2
  res_sum = sum(mid_df['e_y'] + mid_df['e_x'])
  return(res_sum)
}

solve_MultiBiFirst_func <- function(t,state,parms){
  with(as.list(c(state,parms)),{
    # error1 <- rnorm(mean = 0, sd = 0.0005,n=1)
    # error2 <- rnorm(mean = 0, sd = 0.0005,n=1)
    dX <- a*x + b*y
    dY <- c*x + d*y
    list(c(dX,dY))
  })# end with(as.list ...
}

Predictor_MultiBiFirst_func <- function(detailtable,model){
  mid_detailtable = detailtable

  var_model = model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]

  detailtable[paste('solver_',var_model[1,'field'],sep='')]= NA
  detailtable[paste('solver_',var_model[2,'field'],sep='')]= NA

  times = sort(unique(mid_detailtable[,'time']),decreasing = FALSE)

  for(i in unique(detailtable[,'subject'])){
    mid_indi_table = mid_detailtable[which(mid_detailtable['subject'] == i),]
    # There are same,so use mean
    mid_solver_df = deSolve::ode(y=c(x = mean(mid_indi_table[,'initXt0']) + mean(mid_indi_table[,'etaInitXi1']),
                                    y = mean(mid_indi_table[,'initYt0']) + mean(mid_indi_table[,'etaInitYi1'])),
                                times=times,
                                func=solve_MultiBiFirst_func,
                                parms=c(a = mean(mid_indi_table[,'beta1']) + mean(mid_indi_table[,'etaI1']),
                                        b = mean(mid_indi_table[,'beta2']) + mean(mid_indi_table[,'etaI2']),
                                        c = mean(mid_indi_table[,'beta3']) + mean(mid_indi_table[,'etaI3']),
                                        d = mean(mid_indi_table[,'beta4']) + mean(mid_indi_table[,'etaI4']))
                                )
    mid_solver_df = data.frame(mid_solver_df)
    for(i_mid in rownames(mid_indi_table)){
      detailtable[i_mid,paste('solver_',var_model[1,'field'],sep='')] = mid_solver_df[mid_solver_df['time'] == mid_indi_table[i_mid,'time'],'x']
      detailtable[i_mid,paste('solver_',var_model[2,'field'],sep='')] = mid_solver_df[mid_solver_df['time'] == mid_indi_table[i_mid,'time'],'y']
    }
  }
  return(detailtable)
}

SquareR_MultiBiFirst_func <- function(Predictor_data,model){
  message('Estimating R_squared')

  var_model = model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]

  r_squared01 = stats::cor(Predictor_data[,var_model[1,'field']],Predictor_data[,paste('solver_',var_model[1,'field'],sep = '')]) **2
  r_squared02 = stats::cor(Predictor_data[,var_model[2,'field']],Predictor_data[,paste('solver_',var_model[2,'field'],sep = '')]) **2
  name1 = paste('r_squared_',var_model[1,'field'],' = ',r_squared01,sep = "")
  name2 = paste('r_squared_',var_model[2,'field'],' = ',r_squared02,sep = "")
  return(c(name1,name2))
}
RMSE_MultiBiFirst_func <- function(Predictor_data,model){
  message('Estimating RMSE')

  var_model = model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]

  SSxe <- sum((Predictor_data[,var_model[1,'field']] - Predictor_data[,paste('solver_',var_model[1,'field'],sep = '')])**2)
  RMSEx <- sqrt(SSxe / NROW(Predictor_data[,'time']))
  SSye <- sum((Predictor_data[,var_model[2,'field']] - Predictor_data[,paste('solver_',var_model[2,'field'],sep = '')])**2)
  RMSEy <- sqrt(SSye / NROW(Predictor_data[,'time']))
  return(c(RMSEx,RMSEy))
}
Hessian_MultiBiFirst_func <- function(hessian){
  message('Estimating Hessian')
  hessian[hessian<0] <- NaN
  mid_SE <- sqrt(1/hessian)
  return(mid_SE)
}

bayesian_population_func <- function(userdata,model,x0){
  BetaXToX = x0[1] # beta1
  BetaYToX = x0[2] # beta2
  BetaXToY = x0[3] # beta3
  BetaYToY = x0[4] # beta4
  BetaXt0 = x0[5] # init X
  BetaYt0 = x0[6] # init Y

  times = sort(unique(userdata[,'time']),decreasing = FALSE)
  mid_deEstimate_data <- deSolve::ode(y=c(x = BetaXt0,
                                          y = BetaYt0),
                                      times=times,
                                      func=solve_BinFirst_func,
                                      parms=c(a = BetaXToX,
                                              b = BetaYToX,
                                              c = BetaXToY,
                                              d = BetaYToY))
  mid_df_desolve <- data.frame(mid_deEstimate_data)

  var_model = model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]
  # print(userdata)
  names(mid_df_desolve)[names(mid_df_desolve) == "x"] = 'solver_x'
  names(mid_df_desolve)[names(mid_df_desolve) == "y"] = 'solver_y'
  mid_df_res <- data.frame(seq=1:nrow(userdata),
                           time = userdata[,'time'],
                           x = userdata[,var_model[1,'field']],
                           y = userdata[,var_model[2,'field']],
                           subject = userdata[,'subject'],
                           solver_x = NA,
                           solver_y = NA
  )

  for (i in 1:nrow(mid_df_res)){
    mid_df_res[i,'solver_x'] = mid_df_desolve[mid_df_desolve['time'] == mid_df_res[i,'time'],'solver_x']
    mid_df_res[i,'solver_y'] = mid_df_desolve[mid_df_desolve['time'] == mid_df_res[i,'time'],'solver_y']
  }
  mid_df_res['e_y'] = (mid_df_res['y'] - mid_df_res['solver_y']) ** 2
  mid_df_res['e_x'] = (mid_df_res['x'] - mid_df_res['solver_x']) ** 2
  mid_sum_res = sum(mid_df_res['e_y'] + mid_df_res['e_x'])

  mid_gpfit_X = GPfit::GP_fit(X=mid_df_res[,'x'],
                              Y=mid_df_res[,'solver_x'],
                              corr = list(type = "exponential", power = 1.95))
  mid_gpfit_Y = GPfit::GP_fit(X=mid_df_res[,'y'],
                              Y=mid_df_res[,'solver_y'],
                              corr = list(type = "exponential", power = 1.95))
  x_new = seq(0, 1, length.out = nrow(mid_df_res))
  mid_gppred_X = GPfit::predict.GP(mid_gpfit_X,xnew = data.frame(x = x_new))
  mid_gppred_Y = GPfit::predict.GP(mid_gpfit_Y,xnew = data.frame(x = x_new))
  mid_mu = mid_gppred_X$Y_hat + mid_gppred_Y$Y_hat
  mid_sigma <- sqrt(mid_gppred_X$MSE) + sqrt(mid_gppred_Y$MSE)

  mid_y_best <- min(mid_df_res[,'solver_x']) + min(mid_df_res[,'solver_y'])
  eps <- 0.01
  expected_improvement <- numeric()

  ei_calc <- function(m, s,mid_y_best) {
    if (s == 0) {
      return(0)
    }
    Z <- (m - mid_y_best - eps)/s
    expected_imp <- (m - mid_y_best - eps) * pnorm(Z) + s * dnorm(Z)
    return(expected_imp)
  }

  for (i in 1:length(mid_mu)) {
    expected_improvement[i] <- ei_calc(m = mid_mu[i],s =  mid_sigma[i],mid_y_best=mid_y_best)
  }

  return(mid_sum_res)
}

