
#' Title
#'
#' @param data User's data
#' @param model Model's class is dataframe.
#' @param guess Guess values that contain coefficient and initial values.
#' @param method "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN" and "Brent"
#' @param guess2 Guess values of multilevel that contain coefficient and initial values.
#' @param method2 "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN" and "Brent"
#'
#' @return The result of optimization,SE,RMSE,r-squared,users's data,predictor data and output table.

Solver_MultiUniSec_func <- function(data,model,guess,method,guess2,method2){
  # message(method,method2)
  # cat('### Program will fit the data with multilevel bivariate first-order differential equations')
  # print('### Firstly, we should estimate the population parameter of differential equations')
  # ###the variable of model
  var_model = model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]
  sub_model = model
  sub_model = sub_model[which(sub_model['operator'] == '~'),]
  if(all(is.na(guess))){
    guess = c(0,0,data[1,var_model[1,'field']],0,0,0)
  }
  x0 = guess
  cat('Estimating population parameter \n')
  cat('Waiting...\n')
  if(method == 'bayesian'){
    min_population_data <- bayesian_population_func(userdata = data,
                                                    model = model,
                                                    x0 = x0)
  }else{
    min_population_data <- stats::optim(x0,
                                        min_MultiUniSec_Population_func,
                                        NULL,
                                        method = method,
                                        userdata = data,
                                        model = model,
                                        hessian=TRUE
                                        )
  }
  fixvalues = min_population_data$par
  detailtable = min_multiUniSec_IndiOnebyOne_func(guess2,
                                           userdata=data,
                                           model=model,
                                           fixvalues=fixvalues,
                                           method = method2
                                           )

  outputDE0 <- 'Fixed effects:'
  outputDE1 <- paste(var_model[1,'field'],'(2)=',min_population_data$par[1],'*',var_model[1,'field'],'+',min_population_data$par[2],' * ',var_model[1,'field'],'(1)',sep = '')
  outputDE2 <- paste('init01:',min_population_data$par[3])

  Predictor_data = Predictor_MultiUniSec_func(detailtable,model)
  squareR_data = SquareR_MultiUniSec_func(Predictor_data,model)
  user_rmse_data = RMSE_MultiUniSec_func(Predictor_data,model)
  user_rmse_data <- paste('RMSE_',var_model[1,'field'],' = ',user_rmse_data,sep='')
  SE_data = Hessian_MultiUniSec_func(min_population_data$hessian)

  outputDE3 <- 'Random effects:'
  outputDE4 <- data.frame(seq=1:4,
                          Groups = NA,
                          Name = NA,
                          Variance = NA,
                          Std.Dev. = NA)
  outputDE4[1,'Groups'] = 'Subject'
  outputDE4[1,'Name'] = paste(var_model[1,'field'],'(0) to ',var_model[1,'field'],'(1)',sep = "")
  outputDE4[1,'Variance'] = stats::var(Predictor_data[Predictor_data['time'] == Predictor_data[1,'time'],'etaI1'])
  outputDE4[1,'Std.Dev.'] = sqrt(outputDE4[1,'Variance'])
  outputDE4[2,'Groups'] = 'Subject'
  outputDE4[2,'Name'] = paste(var_model[1,'field'],'(1) to ',var_model[1,'field'],'(1)',sep = "")
  outputDE4[2,'Variance'] = stats::var(Predictor_data[Predictor_data['time'] == Predictor_data[1,'time'],'etaI2'])
  outputDE4[2,'Std.Dev.'] = sqrt(outputDE4[2,'Variance'])
  outputDE4[3,'Groups'] = 'Subject'
  outputDE4[3,'Name'] = paste('init_',var_model[1,'field'],sep='')
  outputDE4[3,'Variance'] = stats::var(Predictor_data[Predictor_data['time'] == Predictor_data[1,'time'],'etaInitXi1'])
  outputDE4[3,'Std.Dev.'] = sqrt(outputDE4[3,'Variance'])
  outputDE4[4,'Groups'] = 'Residual'
  outputDE4[4,'Name'] = ''
  outputDE4[4,'Variance'] = stats::var(Predictor_data[,var_model[1,'field']] - Predictor_data[,paste('solver_',var_model[1,'field'],sep='')])
  outputDE4[4,'Std.Dev.'] = sqrt(outputDE4[4,'Variance'])
  outputDE5 <- paste('Number of obs:',nrow(Predictor_data),', groups:',sub_model[1,'subject'],',',length(unique(Predictor_data[,'subject'])))

  outtable <- data.frame()
  outtable[seq(1,3,by=1),'parameter'] = c(paste(var_model[1,'field'],'(0) to ',var_model[1,'field'],'(2)',sep = ""),
                                          paste(var_model[1,'field'],'(1) to ',var_model[1,'field'],'(2)',sep = ""),
                                          paste('init_',var_model[1,'field'],sep='')
                                          #'init02',
                                          #'init03',
                                          #'init04'
                                          )
  outtable[,'value'] = min_population_data$par[1:3]
  outtable[,'SE'] = c(SE_data[1,1],
                      SE_data[2,2],
                      SE_data[3,3]
                      # SE_data[4,4],
                      # SE_data[5,5],
                      # SE_data[6,6])
                      )
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
              DifferentialEquational = list(outputDE0,outputDE1,outputDE2,outputDE3,outputDE4,outputDE5),
              Predictor = Predictor_data,
              Rsquared = squareR_data,
              Rmse = user_rmse_data,
              SE = SE_data,
              table = outtable,
              IsConvergence = IsConvergence))
}

solve_UniSec_func <- function(t, y, dy, parms){
  res1 <- dy[1] - y[2]
  res2 <- dy[2] - parms[1] * y[2] - parms[2] * y[1]
  list(c(res1, res2))
}

min_MultiUniSec_Population_func <- function(x0,userdata,model){
  Beta1 = x0[1] # beta1
  Beta2 = x0[2] # beta2
  init01 = x0[3] # y
  init02 = x0[4] # dy
  init03 = x0[5] # dy
  init04 = x0[6] # d2y
  # Get the time and sort it in descending order of time.
  # Now we have not considered the time of NA. to-do
  times = sort(unique(userdata[,'time']),decreasing = FALSE)
  mid_deEstimate_data <- deSolve::daspk(y=c(y01 = init01,
                                            y02 = init02),
                                        dy = c(dy01 = init03,
                                               dy02 = init04),
                                        times=times,
                                        res=solve_UniSec_func,
                                        parms=c(a = Beta1,
                                                b = Beta2))
  mid_df_desolve <- data.frame(mid_deEstimate_data)
  # Get the variables of the model
  var_model = model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]
  # print(userdata)
  names(mid_df_desolve)[names(mid_df_desolve) == "y01"] = 'solver_y'
  mid_df_res <- data.frame(seq=1:nrow(userdata),
                           time = userdata[,'time'],
                           y = userdata[,var_model[1,'field']],
                           subject = userdata[,'subject'],
                           solver_y = NA
  )
  # When the time is the same, userdata and solverdata will be merged.
  for (i in 1:nrow(mid_df_res)){
    mid_df_res[i,'solver_y'] = mid_df_desolve[mid_df_desolve['time'] == mid_df_res[i,'time'],'solver_y']
  }
  mid_df_res['e_y'] = (mid_df_res['y'] - mid_df_res['solver_y']) ** 2
  mid_sum_res = sum(mid_df_res['e_y'])
  return(mid_sum_res)
}
min_multiUniSec_IndiOnebyOne_func <- function(x0,userdata,model,fixvalues,method){
  # ###the variable of model
  var_model = model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]
  # ###deal the data
  mid_onebyone_df = data.frame(seq=1:nrow(userdata))
  mid_onebyone_df['time'] = userdata[,'time']
  mid_onebyone_df[var_model[1,'field']] = userdata[,var_model[1,'field']]
  mid_onebyone_df['subject'] = userdata[,'subject']
  mid_onebyone_df['beta1'] = fixvalues[1]
  mid_onebyone_df['beta2'] = fixvalues[2]
  mid_onebyone_df['init01'] = fixvalues[3]
  mid_onebyone_df['init02'] = fixvalues[4]
  mid_onebyone_df['init03'] = fixvalues[5]
  mid_onebyone_df['init04'] = fixvalues[6]
  # ###calculate/estimate the random effect
  # Identify the fix random variable
  mid_Identify_fixrand = model
  mid_Identify_fixrand = mid_Identify_fixrand[which(mid_Identify_fixrand['operator'] == '~'),]
  # ###Judge the eta guess values
  etaX0 =rep(0.01,times=6) #Not in use.
  for(i in unique(userdata[,'subject'])){
    mytry <- tryCatch({
      mid_individual_data <- stats::optim(x0,
                                          min_onebyone_UniSec_Eta_func,
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
      message('Error @  ',err)
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
    fixRand_model = fixRand_model[which(fixRand_model['operator'] == '~'),'fixRand']
    fixRand_list = unlist(strsplit(fixRand_model,'[+]'))
    fixRand_list = trimws(fixRand_list,which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
    if(!('1' %in% fixRand_list)){
      mid_coef[3] = 0

      mid_coef[4] = 0
      mid_coef[5] = 0
      mid_coef[6] = 0
    }
    if(!(var_model[1,'field'] %in% fixRand_list)){
      mid_coef[1] = 0

      mid_coef[4] = 0
      mid_coef[5] = 0
      mid_coef[6] = 0
    }
    if(!(paste(var_model[1,'field'],'(1)',sep = '') %in% fixRand_list)){
      mid_coef[2] = 0

      mid_coef[4] = 0
      mid_coef[5] = 0
      mid_coef[6] = 0
    }

    mid_onebyone_df[userdata['subject']==i,'etaI1'] = mid_coef[1]
    mid_onebyone_df[userdata['subject']==i,'etaI2'] = mid_coef[2]
    mid_onebyone_df[userdata['subject']==i,'etaInitXi1'] = mid_coef[3]
    mid_onebyone_df[userdata['subject']==i,'etaInitXi2'] = mid_coef[4]
    mid_onebyone_df[userdata['subject']==i,'etaInitXi3'] = mid_coef[5]
    mid_onebyone_df[userdata['subject']==i,'etaInitXi4'] = mid_coef[6]
  }
  return(mid_onebyone_df)
}
min_onebyone_UniSec_Eta_func <- function(x0,userdata,model,fixvalues){
  etaI1 = 1 # eta1
  etaI2 = 2 # eta2
  etaInitXi1 = 3 # etaXit0
  etaInitXi2 = 4 # etaXit1
  etaInitXi3 = 5 # etaXit2
  etaInitXi4 = 6 # etaXit3

  var_model = model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]

  fixRand_model = model
  fixRand_model = fixRand_model[which(fixRand_model['operator'] == '~'),'fixRand']
  fixRand_list = unlist(strsplit(fixRand_model,'[+]'))
  fixRand_list = trimws(fixRand_list,which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
  if(!('1' %in% fixRand_list)){
    x0[3] = 0
  }
  if(!(var_model[1,'field'] %in% fixRand_list)){
    x0[1] = 0
  }
  if(!(paste(var_model[1,'field'],'(1)',sep = '') %in% fixRand_list)){
    x0[2] = 0
  }

  mid_df = data.frame(seq=1:nrow(userdata),
                      time = userdata[,'time'],
                      y = userdata[,var_model[1,'field']],
                      subject = userdata[,'subject'],
                      solver_y = NA
  )
  times = sort(unique(userdata[,'time']),decreasing = FALSE)
  mid_deEstimate_data <- deSolve::daspk(y=c(y01 =fixvalues[3] +  x0[etaInitXi1],
                                            y02 =fixvalues[4] + x0[etaInitXi2]),
                                        dy = c(dy01 =fixvalues[5] + x0[etaInitXi3],
                                               dy02 =fixvalues[6] + x0[etaInitXi4]),
                                        times=times,
                                        res=solve_UniSec_func,
                                        parms=c(a =fixvalues[1] + x0[etaI1],
                                                b =fixvalues[2] + x0[etaI2]))
  mid_deEstimate_data = data.frame(mid_deEstimate_data)
  for(i_mid in rownames(mid_df)){
    mid_df[i_mid,'solver_y'] = mid_deEstimate_data[mid_deEstimate_data['time'] == mid_df[i_mid,'time'],'y01']
  }
  mid_df['e_y'] = (mid_df['y'] - mid_df['solver_y']) **2
  res_sum = sum(mid_df['e_y'])
  return(res_sum)
}
Predictor_MultiUniSec_func <- function(detailtable,model){
  mid_detailtable = detailtable

  var_model = model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]

  detailtable[paste('solver_',var_model[1,'field'],sep='')]= NA

  times = sort(unique(mid_detailtable[,'time']),decreasing = FALSE)

  for(i in unique(detailtable[,'subject'])){
    mid_indi_table = mid_detailtable[which(mid_detailtable['subject'] == i),]
    # There are same,so use mean
    mid_solver_df = deSolve::daspk(y=c(y01 =mean(mid_indi_table[,'init01']) + mean(mid_indi_table[,'etaInitXi1']),
                                       y02 =mean(mid_indi_table[,'init02']) + mean(mid_indi_table[,'etaInitXi2'])),
                                   dy = c(dy01 =mean(mid_indi_table[,'init03']) + mean(mid_indi_table[,'etaInitXi3']),
                                          dy02 =mean(mid_indi_table[,'init04']) + mean(mid_indi_table[,'etaInitXi4'])),
                                   times=times,
                                   res=solve_UniSec_func,
                                   parms=c(a =mean(mid_indi_table[,'beta1']) + mean(mid_indi_table[,'etaI1']),
                                           b =mean(mid_indi_table[,'beta2']) + mean(mid_indi_table[,'etaI2']))
                                   )
    mid_solver_df = data.frame(mid_solver_df)
    for(i_mid in rownames(mid_indi_table)){
      detailtable[i_mid,paste('solver_',var_model[1,'field'],sep='')] = mid_solver_df[mid_solver_df['time'] == mid_indi_table[i_mid,'time'],'y01']
    }
  }
  return(detailtable)
}
SquareR_MultiUniSec_func <- function(Predictor_data,model){
  message('Estimating R_squared')

  var_model = model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]

  r_squared01 = stats::cor(Predictor_data[,var_model[1,'field']],Predictor_data[,paste('solver_',var_model[1,'field'],sep = '')]) **2
  name1 = paste('r_squared',var_model[1,'field'],' = ',r_squared01,sep = "")
  return(name1)
}
RMSE_MultiUniSec_func <- function(Predictor_data,model){
  message('Estimating RMSE')

  var_model = model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]

  SSxe <- sum((Predictor_data[,var_model[1,'field']] - Predictor_data[,paste('solver_',var_model[1,'field'],sep = '')])**2)
  RMSEx <- sqrt(SSxe / NROW(Predictor_data[,'time']))
  return(RMSEx)
}
Hessian_MultiUniSec_func <- function(hessian){
  message('Estimating Hessian')
  hessian[hessian<0] <- NaN
  mid_SE <- sqrt(1/hessian)
  return(mid_SE)
}
bayesian_population_func <- function(userdata,model,x0){
  message('Not supported now!')
}
