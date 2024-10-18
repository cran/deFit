
#' Solver of Multilevel univariate second-order differential equation
#'
#' @param init_list a list of init_func
#'
#' @return a list
#' @export
Solver_MultiUniSecond_func <- function(init_list){
  userdata = init_list$userdata
  var_model = init_list$var_model
  guess = init_list$guess
  method = init_list$method
  subject_model = init_list$subject_model
  modelDF = init_list$modelDF
  fixed_first = init_list$fixed_first

  message('Program will fit the data with multilevel univariate second-order differential equation.')
  message('The multilevel differential equations are:')
  message('X(2) = (beta1 + etaI1) * X + (beta2 + etaI2) * X(1)')
  message('Optimizing...')

  ## identify variables and information
  var_notime_model = var_model[var_model$field != 'time','field']
  var_notime_hat = paste(var_notime_model,'_hat',sep = "")
  var_notime_data = var_model[var_model$field != 'time','variable']
  fix_model = modelDF[modelDF['operator'] == '~',]

  ############################
  # optimization
  #--------------------------
  if(all(method == 'bayesian')){
    warning('Bayesian Not supported')
  }else if(fixed_first == FALSE){
    calc_data = calc_MultiUniSec_United_func(userdata,var_model,guess,method,subject_model,fix_model)
  }else{
    calc_data = calc_MultiUniSec_func(userdata,var_model,guess,method,subject_model,fix_model)
  }
  predict_data = calc_data['predict_data'][[1]]
  predict_fixed_data = calc_data['predict_fixed_data'][[1]]
  ## equation
  equation1 = paste(var_notime_model[1],"(1) = ",calc_data['fixed_effects'][[1]]$par[1]," * ",var_notime_model[1]," + ",calc_data['fixed_effects'][[1]]$par[2]," * ",var_notime_model[1],sep = "")
  equation2 = paste("Init t0_",var_notime_model[1],": ",calc_data['fixed_effects'][[1]]$par[3],sep = "")
  equation = list(equation1,equation2)
  ## table
  table <- data.frame()
  table[seq(1,4,by=1),'parameter'] = c(paste(var_notime_model[1],"(0) to ",var_notime_model[1],"(2)",sep = ""),
                                       paste(var_notime_model[1],"(1) to ",var_notime_model[1],"(2)",sep = ""),
                                       paste("init t0_",var_notime_model[1],sep = ""),
                                       paste("init dt0_",var_notime_model[1],sep = ""))
  table[1,'value'] = calc_data['fixed_effects'][[1]]$par[1]
  table[2,'value'] = calc_data['fixed_effects'][[1]]$par[2]
  table[3,'value'] = calc_data['fixed_effects'][[1]]$par[3]
  table[4,'value'] = calc_data['fixed_effects'][[1]]$par[4]
  return(list(solve_data=calc_data['fixed_effects'][[1]],
              userdata=userdata,
              predict_data=predict_data,
              predict_fixed_data=predict_fixed_data,
              table=table,
              equation=equation,
              random_effects = calc_data['random_effects'][[1]]))
}

############################
# calculate coefficents of differential equation
#--------------------------
calc_MultiUniSec_func <- function(userdata,var_model,guess,method,subject_model,fix_model){
  ## identify variables and information
  var_notime_model = var_model[var_model$field != 'time','field']
  var_notime_hat = paste(var_notime_model,'_hat',sep = "")
  var_notime_fixed_hat = paste(var_notime_model,'_fixed_hat',sep = "")
  var_notime_data = var_model[var_model$field != 'time','variable']

  ## fix effect
  calcFix_data = stats::optim(guess[[1]],
                              calcFix_MultiUniSec_func,
                              NULL,
                              method = method[[1]],
                              userdata = userdata,
                              var_notime_model = var_notime_model,
                              subject_model =subject_model,
                              hessian=TRUE
  )
  #--------------------------
  ## random effects
  randEffect_data = data.frame(Subject=NULL,
                               EtaI_1=NULL,
                               EtaI_2=NULL,
                               EtaI_y1=NULL,
                               EtaI_y2=NULL)
  for (i_subject in 1:length(unique(userdata[,subject_model]))) {
    mid_subject_name = unique(userdata[,subject_model])
    cat('Estimating random effects ',mid_subject_name[i_subject],'\n')
    # print(randEffect_data)
    mid_userdata = userdata[userdata[subject_model] == mid_subject_name[i_subject],]
    randEffect_list = list(Subject=mid_subject_name[i_subject],
                           EtaI_1=0,
                           EtaI_2=0,
                           EtaI_y1=0,
                           EtaI_y2=0)
    tryCatch(
      {
        mid_calcRandom_data = stats::optim(guess[[2]],
                                           calcRandom_MultiUniSec_func,
                                           NULL,
                                           method = method[[2]],
                                           userdata = mid_userdata,
                                           var_notime_model = var_notime_model,
                                           subject_model = subject_model,
                                           calcFix_data = calcFix_data,
                                           fix_model = fix_model,
                                           hessian=TRUE)
        randEffect_list = list(Subject=mid_subject_name[i_subject],
                               EtaI_1=mid_calcRandom_data$par[1],
                               EtaI_2=mid_calcRandom_data$par[2],
                               EtaI_y1=mid_calcRandom_data$par[3],
                               EtaI_y2=mid_calcRandom_data$par[4])
      },
      error = function(e){
        cat("Error occurred: ", conditionMessage(e), "\n")
      }
    )

    # binding random effects into dataframe
    randEffect_data = rbind(randEffect_data,randEffect_list)
  }
  #--------------------------
  ## predict fixed effects data

  times = seq(min(userdata$time),max(userdata$time),length.out=length(unique(userdata[,'time'])))
  parms_fix = c(beta1 = calcFix_data$par[1],
            beta2 = calcFix_data$par[2])
  state_fix = c(y1 = calcFix_data$par[3],
            y2 = calcFix_data$par[4])
  predict_fixed_data <- deSolve::ode(y=state_fix,
                                   times=times,
                                   func=solve_MultiUniSec_func,
                                   parms=parms_fix)
  predict_fixed_data = data.frame(predict_fixed_data)
  names(predict_fixed_data) = c('time',var_notime_fixed_hat)

  ## predict data
  predict_data = data.frame(Subject=NULL,
                            time=NULL,
                            x=NULL,
                            y=NULL)
  for (i_subject in 1:NROW(randEffect_data)) {
    i_subject_parameters = unlist(randEffect_data[i_subject,c(2,3,4,5)])
    parms = c(beta1 = calcFix_data$par[1] + i_subject_parameters[1][[1]],
              beta2 = calcFix_data$par[2] + i_subject_parameters[2][[1]])
    state = c(y1 = calcFix_data$par[3] + i_subject_parameters[3][[1]],
              y2 = calcFix_data$par[4] + i_subject_parameters[4][[1]])
    mid_predict_data <- deSolve::ode(y=state,
                                     times=times,
                                     func=solve_MultiUniSec_func,
                                     parms=parms)
    mid_predictDF_data = data.frame(Subject=rep(randEffect_data[i_subject,1],length(mid_predict_data[,'time'])),
                                    time=mid_predict_data[,'time'],
                                    y1 =mid_predict_data[,'y1']
    )
    predict_data = rbind(predict_data,mid_predictDF_data)
  }
  names(predict_data) = c('Subject','time',var_notime_hat)
  return(list(fixed_effects = calcFix_data,
              random_effects = randEffect_data,
              predict_fixed_data = predict_fixed_data,
              predict_data = predict_data))
}

############################
# Fix effects
#--------------------------
calcFix_MultiUniSec_func <- function(x0,userdata,var_notime_model,subject_model){
  times <- seq(min(userdata$time),max(userdata$time),length.out=length(unique(userdata[,'time'])))
  parms = c(beta1 = x0[1],
            beta2 = x0[2])
  state = c(y1 = x0[3],
            y2 = x0[4])
  mid_solve_data <- deSolve::ode(y=state,
                                 times=times,
                                 func=solve_MultiUniSec_func,
                                 parms=parms)

  mid_solve_data = data.frame(mid_solve_data)
  mid_df_res <- data.frame(seq=1:nrow(userdata),
                           time = times,
                           y1 = userdata[,var_notime_model[1]],
                           subject = userdata[,subject_model],
                           solver_y1 = NA
                           )

  # When the time is the same, userdata and solverdata will be merged.
  for (i in 1:nrow(mid_df_res)){
    mid_df_res[i,'solver_y1'] = mid_solve_data[mid_solve_data['time'] == mid_df_res[i,'time'],'y1']
  }
  mid_df_res['e_y'] = (mid_df_res['y1'] - mid_df_res['solver_y1']) ** 2
  mid_sum_res = sum(mid_df_res['e_y'])

  ## The sum of squared may be infinitely
  if(!(is.finite(mid_sum_res))){
    mid_df_res[,'e_y'] = (mid_df_res[,'y1']) ** 2
    print('waiting')
    mid_sum_res = sum(mid_df_res['e_y'])
  }

  return(mid_sum_res)
}

############################
# Random effects
#--------------------------
calcRandom_MultiUniSec_func <- function(x0,userdata,var_notime_model,subject_model,calcFix_data,fix_model){

  fixRand_string = fix_model[1,'fixRand']
  fixRand_list = unlist(strsplit(fixRand_string,'\\+'))
  fixRand_list = gsub(" ", "", fixRand_list)
  if(!(var_notime_model[1] %in% fixRand_list)){
    x0[1] = 0
  }
  if(!(paste(var_notime_model[1],"(1)",sep = "") %in% fixRand_list)){
    x0[2] = 0
  }
  if(!("1" %in% fixRand_list)){
    x0[3] = 0
  }
  times <- seq(min(userdata$time),max(userdata$time),length.out=length(unique(userdata[,'time'])))
  parms = c(beta1 = calcFix_data$par[1] + x0[1], # beta1 is etaI1
            beta2 = calcFix_data$par[2] + x0[2] # beta2 is etaI2
            )
  state = c(y1 = calcFix_data$par[3] + x0[3],
            y2 = calcFix_data$par[4] + x0[4])
  mid_solve_data <- deSolve::ode(y=state,
                                 times=times,
                                 func=solve_MultiUniSec_func,
                                 parms=parms)

  mid_solve_data = data.frame(mid_solve_data)
  mid_solve_data = stats::setNames(mid_solve_data,c('time','solver_y1','solver_y2'))
  mid_df_res <- data.frame(seq=1:nrow(userdata),
                           time = userdata[,'time'],
                           y1 = userdata[,var_notime_model[1]],
                           subject = userdata[,subject_model]
                           )

  # When the time is the same, userdata and solverdata will be merged.
  mid_df_res = merge(mid_df_res,mid_solve_data,by='time',all = TRUE)
  for (col in names(mid_df_res)) {
    mid_df_res[is.na(mid_df_res[[col]]), col] <- 0
  }
  mid_df_res['e_y'] = (mid_df_res['y1'] - mid_df_res['solver_y1']) ** 2
  mid_sum_res = sum(mid_df_res['e_y'])

  ## hhhh so happy
  if(!(is.finite(mid_sum_res))){
    mid_df_res[,'e_y'] = (mid_df_res[,'y']) ** 2
    print('waiting')
    mid_sum_res = sum(mid_df_res['e_y'])
  }
  return(mid_sum_res)
}

############################
# Differential equations
#--------------------------
solve_MultiUniSec_func <- function(t,y,parms){
  with(as.list(parms), {
    dy1 <- y[2]
    dy2 <- beta1 * y[1] + beta2 * y[2]
    return(list(c(dy1, dy2)))
  })
}

############################
# calc_MultiUniSec_United_func
#--------------------------
calc_MultiUniSec_United_func <- function(userdata,var_model,guess,method,subject_model,fix_model){
  ## identify variables and information
  var_notime_model = var_model[var_model$field != 'time','field']
  var_notime_hat = paste(var_notime_model,'_hat',sep = "")
  var_notime_fixed_hat = paste(var_notime_model,'_fixed_hat',sep = "")
  var_notime_data = var_model[var_model$field != 'time','variable']

  fixed_random_guess <- guess[[1]]
  len_subject = unique(userdata[subject_model])
  len_subject = as.list(len_subject)
  len_subject = length(len_subject[[1]])
  fix_rand_guess <- c(fixed_random_guess,rep(0, len_subject * 4))

  ## fixed and random effect
  calcFix_data = stats::optim(fix_rand_guess,
                              calcFix_MultiUniSec_United_func,
                              NULL,
                              method = method[[1]],
                              # control = list(reltol = 1e-6),
                              userdata = userdata,
                              var_notime_model = var_notime_model,
                              subject_model =subject_model,
                              hessian=TRUE
  )
  ## random effects
  randEffect_data = data.frame(Subject=NULL,
                               EtaI_1=NULL,
                               EtaI_2=NULL,
                               EtaI_x=NULL,
                               EtaI_y=NULL)
  rand_li = calcFix_data$par[5:length(calcFix_data$par)]
  mid_subject_name = unique(userdata[subject_model])
  mid_subject_li = as.list(mid_subject_name)
  for(i_subject in seq(0,len_subject-1)){
    randEffect_list = list(Subject=mid_subject_li[[1]][i_subject + 1],
                           EtaI_1=rand_li[1 + i_subject*4],
                           EtaI_2=rand_li[2 + i_subject*4],
                           EtaI_x=rand_li[3 + i_subject*4],
                           EtaI_y=rand_li[4 + i_subject*4])
    randEffect_data = rbind(randEffect_data,do.call("data.frame", randEffect_list))
  }
  ## predict
  times = sort(unique(userdata[,'time']),decreasing = FALSE)

  fixed_parms = c(beta1 = calcFix_data$par[1],
                  beta2 = calcFix_data$par[2])
  fixed_state = c(x = calcFix_data$par[3],
                  y = calcFix_data$par[4])
  predict_fixed_data = deSolve::ode(y=fixed_state,
                                    times=times,
                                    func=solve_MultiUniSec_func,
                                    parms=fixed_parms)
  predict_fixed_data = data.frame(predict_fixed_data)
  names(predict_fixed_data) = c('time',var_notime_fixed_hat)

  ## predict random effects data
  predict_data = data.frame(Subject=NULL,
                            time=NULL,
                            x=NULL,
                            y=NULL)
  for (i_subject in 1:NROW(randEffect_data)) {
    i_subject_parameters = unlist(randEffect_data[i_subject,c(2,3,4,5)])
    parms = c(beta1 = calcFix_data$par[1] + i_subject_parameters[1][[1]],
              beta2 = calcFix_data$par[2] + i_subject_parameters[2][[1]])
    state = c(x = calcFix_data$par[3] + i_subject_parameters[3][[1]],
              y = calcFix_data$par[4] + i_subject_parameters[4][[1]])
    mid_predict_data <- deSolve::ode(y=state,
                                     times=times,
                                     func=solve_MultiUniSec_func,
                                     parms=parms)
    mid_predictDF_data = data.frame(Subject=rep(randEffect_data[i_subject,1],length(mid_predict_data[,'time'])),
                                    time=mid_predict_data[,'time'],
                                    x =mid_predict_data[,'x'],
                                    y =mid_predict_data[,'y']
    )
    predict_data = rbind(predict_data,mid_predictDF_data)
  }
  names(predict_data) = c('Subject','time',var_notime_hat)


  return(list(fixed_effects = calcFix_data,
              random_effects = randEffect_data,
              predict_data = predict_data,
              predict_fixed_data = predict_fixed_data))
}
#--------------------------
calcFix_MultiUniSec_United_func <- function(x0,userdata,var_notime_model,subject_model){
  # times <- sort(unique(userdata[,'time']),decreasing = FALSE)
  times <- seq(min(userdata$time),max(userdata$time),length.out=length(unique(userdata[,'time'])))
  parms = c(beta1 = x0[1],
            beta2 = x0[2])
  state = c(x = x0[3],
            y = x0[4])
  mid_fixed_rand_df = data.frame(rand_beta1 = numeric(),
                                 rand_beta2 = numeric(),
                                 rand_init1 = numeric(),
                                 rand_init2 = numeric(),
                                 stringsAsFactors = FALSE)
  subjects <- unique(userdata[[subject_model]])
  for (idx in seq_along(subjects)) {
    mid_idx <- idx * 4
    mid_i_key <- data.frame(
      rand_beta1 = x0[mid_idx + 1],
      rand_beta2 = x0[mid_idx + 2],
      rand_init1 = x0[mid_idx + 3],
      rand_init2 = x0[mid_idx + 4],
      stringsAsFactors = FALSE
    )
    mid_fixed_rand_df <- rbind(mid_fixed_rand_df, mid_i_key)
  }
  mid_fixed_rand_df[['fixed_beta1']] = x0[1]
  mid_fixed_rand_df[['fixed_beta2']] = x0[2]
  mid_fixed_rand_df[['fixed_init1']] = x0[3]
  mid_fixed_rand_df[['fixed_init2']] = x0[4]
  mid_idx = 1
  res_sum = 0
  for(i in subjects){
    mid_parms = c(beta1 = x0[1] + mid_fixed_rand_df[mid_idx,'rand_beta1'],
                  beta2 = x0[2] + mid_fixed_rand_df[mid_idx,'rand_beta2'])
    mid_state = c(x = x0[3] + mid_fixed_rand_df[mid_idx,'rand_init1'],
                  y = x0[4] + mid_fixed_rand_df[mid_idx,'rand_init2'])
    mid_solve_data <- deSolve::ode(y=mid_state,
                                   times=times,
                                   func=solve_MultiUniSec_func,
                                   parms=mid_parms)
    mid_df_res = data.frame(mid_solve_data)
    mid_df_res['mid_subject'] = mid_idx
    colnames(mid_df_res)[colnames(mid_df_res) == "x"] <- "solver_x"
    mid_userdata = userdata[userdata[subject_model] == i,]
    mid_res_df = merge(mid_userdata, mid_df_res, by = "time", all.x = TRUE)
    solver_li = c('solver_x','solver_y')
    solver_idex = 1
    mid_sum_res = 0
    for(i in var_notime_model){
      mid_col_name = paste('e_',i,sep = '')
      mid_res_df[mid_col_name] = (mid_res_df[i] - mid_res_df[solver_li[solver_idex]]) ** 2
      mid_sum_res = mid_sum_res + sum(mid_res_df[mid_col_name])
      solver_idex = solver_idex + 1
    }
    ## hhhh
    if(!(is.finite(mid_sum_res))){
      mid_sum_res = 0
      for (i in var_notime_model){
        mid_col_name = paste('e_',i,sep = '')
        mid_res_df[mid_col_name] = (mid_res_df[i]) ** 2
        mid_sum_res = mid_sum_res + sum(mid_res_df[mid_col_name])
        solver_idex = solver_idex + 1
      }
    }
    res_sum = res_sum + mid_sum_res
    mid_idx = mid_idx + 1
  }

  return(res_sum)
}
