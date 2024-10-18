
#' Solver of Multilevel bivariate first-order differential equation
#'
#' @param init_list a list of init_func
#'
#' @return a list
#' @export
Solver_MultiBiFirst_func <- function(init_list){
  ###################################
  # init_list come from init_func
  userdata = init_list$userdata
  var_model = init_list$var_model
  guess = init_list$guess
  method = init_list$method
  subject_model = init_list$subject_model
  modelDF = init_list$modelDF
  fixed_first = init_list$fixed_first

  message('Program will fit the data with multilevel bivariate first-order differential equations.')
  message('The multilevel differential equations are:')
  message('dx/dt = (beta1 + etaI1) * x + (beta2 + etaI2) * y')
  message('dy/dt = (beta3 + etaI3) * x + (beta4 + etaI4) * y')
  message('Optimizing...')

  ## identify variables and information
  var_notime_model = var_model[var_model$field != 'time','field']
  var_notime_hat = paste(var_notime_model,'_hat',sep = "")
  var_notime_data = var_model[var_model$field != 'time','variable']
  fix_model = modelDF[modelDF['operator'] == '~',]

  ############################
  # optimization
  #--------------------------
  # bayesian is not supported. fixed_fist=TRUE, Random effects will be estimate by two steps
  if(all(method == 'bayesian')){
    warning('Bayesian Not supported')
  }else if(fixed_first == FALSE){
    calc_data = calc_MultiBiFirst_United_func(userdata,var_model,guess,method,subject_model,fix_model)
  }else{
    calc_data = calc_MultiBiFirst_func(userdata,var_model,guess,method,subject_model,fix_model)
  }

  predict_data = calc_data['predict_data'][[1]]
  predict_fixed_data = calc_data['predict_fixed_data'][[1]]

  ## equation
  equation1 = paste(var_notime_model[1],"(1) = ",calc_data['fixed_effects'][[1]]$par[1]," * ",var_notime_model[1]," + ",calc_data['fixed_effects'][[1]]$par[2]," * ",var_notime_model[2],sep = "")
  equation2 = paste(var_notime_model[2],"(1) = ",calc_data['fixed_effects'][[1]]$par[3]," * ",var_notime_model[1]," + ",calc_data['fixed_effects'][[1]]$par[4]," * ",var_notime_model[2],sep = "")
  equation3 = paste("Init t0_",var_notime_model[1],": ",calc_data['fixed_effects'][[1]]$par[5],", Init t0_",var_notime_model[2],": ",calc_data['fixed_effects'][[1]]$par[6],sep = "")
  equation = list(equation1,equation2,equation3)
  ## table
  table <- data.frame()
  table[seq(1,6,by=1),'parameter'] = c(paste(var_notime_model[1],"(0) to ",var_notime_model[1],"(1)",sep = ""),
                                       paste(var_notime_model[2],"(0) to ",var_notime_model[1],"(1)",sep = ""),
                                       paste(var_notime_model[1],"(0) to ",var_notime_model[2],"(1)",sep = ""),
                                       paste(var_notime_model[2],"(0) to ",var_notime_model[2],"(1)",sep = ""),
                                       paste("init t0_",var_notime_model[1],sep = ""),
                                       paste("init t0_",var_notime_model[2],sep = ""))
  table[1,'value'] = calc_data['fixed_effects'][[1]]$par[1]
  table[2,'value'] = calc_data['fixed_effects'][[1]]$par[2]
  table[3,'value'] = calc_data['fixed_effects'][[1]]$par[3]
  table[4,'value'] = calc_data['fixed_effects'][[1]]$par[4]
  table[5,'value'] = calc_data['fixed_effects'][[1]]$par[5]
  table[6,'value'] = calc_data['fixed_effects'][[1]]$par[6]
  return(list(solve_data=calc_data['fixed_effects'][[1]],
              userdata=userdata,
              predict_data=predict_data,
              predict_fixed_data=predict_fixed_data,
              table=table,
              equation=equation,
              random_effects = calc_data['random_effects'][[1]],
              Neg2LL = calc_data['Neg2LL']))
}

############################
# calculate coefficents of differential equation
#--------------------------
calc_MultiBiFirst_func <- function(userdata,var_model,guess,method,subject_model,fix_model){
  ## identify variables and information
  var_notime_model = var_model[var_model$field != 'time','field']
  var_notime_hat = paste(var_notime_model,'_hat',sep = "")
  var_notime_fixed_hat = paste(var_notime_model,'_fixed_hat',sep = "")
  var_notime_data = var_model[var_model$field != 'time','variable']

  ## fix effect
  guess_e = c(guess[[1]],c(1,1))
  calcFix_data = stats::optim(guess_e,
                              calcFix_MultiBiFirst_func,
                              NULL,
                              method = method[[1]],
                              # control = list(reltol = 1e-6),
                              userdata = userdata,
                              var_notime_model = var_notime_model,
                              subject_model =subject_model,
                              hessian=TRUE
                              )
  Neg2LL = calcFix_data$value[[1]][1] + 3.67
  #--------------------------
  ## random effects
  randEffect_data = data.frame(Subject=NULL,
                               EtaI_1=NULL,
                               EtaI_2=NULL,
                               EtaI_3=NULL,
                               EtaI_4=NULL,
                               EtaI_x=NULL,
                               EtaI_y=NULL)
  for (i_subject in 1:length(unique(userdata[,subject_model]))) {
    mid_subject_name = unique(userdata[,subject_model])
    cat('Estimating random effects ',mid_subject_name[i_subject],'\n')
    mid_userdata = userdata[userdata[subject_model] == mid_subject_name[i_subject],]
    randEffect_list = list(Subject=mid_subject_name[i_subject],
                           EtaI_1=0,
                           EtaI_2=0,
                           EtaI_3=0,
                           EtaI_4=0,
                           EtaI_x=0,
                           EtaI_y=0)
    tryCatch({
      guess_e_rand = c(guess[[2]],c(1,1))
      mid_calcRandom_data = stats::optim(guess_e_rand,
                                         calcRandom_MultiBiFirst_func,
                                         NULL,
                                         method = method[[2]],
                                         # control = list(reltol = 1e-6),
                                         # lower = -Inf,
                                         # upper = Inf,
                                         userdata = mid_userdata,
                                         var_notime_model = var_notime_model,
                                         subject_model = subject_model,
                                         calcFix_data = calcFix_data,
                                         fix_model = fix_model,
                                         hessian=TRUE)
      # binding random effects into dataframe
      randEffect_list = list(Subject=mid_subject_name[i_subject],
                             EtaI_1=mid_calcRandom_data$par[1],
                             EtaI_2=mid_calcRandom_data$par[2],
                             EtaI_3=mid_calcRandom_data$par[3],
                             EtaI_4=mid_calcRandom_data$par[4],
                             EtaI_x=mid_calcRandom_data$par[5],
                             EtaI_y=mid_calcRandom_data$par[6])
    },
    error = function(e){
      # binding random effects into dataframe
      randEffect_list = list(Subject=mid_subject_name[i_subject],
                             EtaI_1=0,
                             EtaI_2=0,
                             EtaI_3=0,
                             EtaI_4=0,
                             EtaI_x=0,
                             EtaI_y=0)
      cat("Error occurred: ", conditionMessage(e), "\n")
    })
    randEffect_data = rbind(randEffect_data,randEffect_list)
  }
  #--------------------------
  ## predict fixed effects data
  times = sort(unique(userdata[,'time']),decreasing = FALSE)

  fixed_parms = c(beta1 = calcFix_data$par[1],
                 beta2 = calcFix_data$par[2],
                 beta3 = calcFix_data$par[3],
                 beta4 = calcFix_data$par[4])
  fixed_state = c(x = calcFix_data$par[5],
            y = calcFix_data$par[6])
  predict_fixed_data = deSolve::ode(y=fixed_state,
                                    times=times,
                                    func=solve_MultiBiFirst_func,
                                    parms=fixed_parms)
  predict_fixed_data = data.frame(predict_fixed_data)
  names(predict_fixed_data) = c('time',var_notime_fixed_hat)

  ## predict random effects data
  predict_data = data.frame(Subject=NULL,
                            time=NULL,
                            x=NULL,
                            y=NULL)
  for (i_subject in 1:NROW(randEffect_data)) {
    i_subject_parameters = unlist(randEffect_data[i_subject,c(2,3,4,5,6,7)])
    i_time = userdata[userdata[subject_model] == randEffect_data[i_subject,1],'time']
    i_time = sort(i_time)
    parms = c(beta1 = calcFix_data$par[1] + i_subject_parameters[1][[1]],
              beta2 = calcFix_data$par[2] + i_subject_parameters[2][[1]],
              beta3 = calcFix_data$par[3] + i_subject_parameters[3][[1]],
              beta4 = calcFix_data$par[4] + i_subject_parameters[4][[1]])
    state = c(x = calcFix_data$par[5] + i_subject_parameters[5][[1]],
              y = calcFix_data$par[6] + i_subject_parameters[6][[1]])
    mid_predict_data <- deSolve::ode(y=state,
                                   times=i_time,
                                   func=solve_MultiBiFirst_func,
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
              predict_fixed_data = predict_fixed_data,
              Neg2LL = Neg2LL))
  }

############################
# Fix effects
#--------------------------
calcFix_MultiBiFirst_func <- function(x0,userdata,var_notime_model,subject_model){
  # times <- sort(unique(userdata[,'time']),decreasing = FALSE)
  times <- seq(min(userdata$time),max(userdata$time),length.out=length(unique(userdata[,'time'])))
  parms = c(beta1 = x0[1],
            beta2 = x0[2],
            beta3 = x0[3],
            beta4 = x0[4])
  state = c(x = x0[5],
            y = x0[6])
  mid_solve_data <- deSolve::ode(y=state,
                                 times=times,
                                 func=solve_MultiBiFirst_func,
                                 parms=parms)

  mid_solve_data = data.frame(mid_solve_data)
  mid_df_res <- data.frame(seq=1:nrow(userdata),
                           time = times,
                           x = userdata[,var_notime_model[1]],
                           y = userdata[,var_notime_model[2]],
                           subject = userdata[,subject_model]
                           )
  # estimating residual in loss function
  guess_xy = c(x0[7],x0[8])
  if(x0[7] == 0){
    guess_xy[1] = 0.00001
  }
  if(x0[8] == 0){
    guess_xy[2] = 0.00001
  }
  log_e_x = log(guess_xy[1]**2)
  log_e_y = log(guess_xy[2]**2)
  inverse_ex = 1/(guess_xy[1]**2)
  inverse_ey = 1/(guess_xy[2]**2)

  # When the time is the same, userdata and solverdata will be merged.
  tryCatch({
    mid_df_res = merge(mid_solve_data,mid_df_res,by='time',all.x=TRUE)
    mid_df_res['e_x'] = log_e_x + ((mid_df_res['x.x'] - mid_df_res['x.y']) ** 2) * inverse_ex
    mid_df_res['e_y'] = log_e_y + ((mid_df_res['y.x'] - mid_df_res['y.y']) ** 2) * inverse_ey
    mid_sum_res = sum(mid_df_res['e_y'] + mid_df_res['e_x'])
  },error={
    ## The sum of squared may be infinitely
    mid_df_res[,'e_y'] = (mid_df_res[,'y']) ** 2
    mid_df_res[,'e_x'] = (mid_df_res[,'x']) ** 2
    mid_sum_res = sum(mid_df_res['e_y'] + mid_df_res['e_x'])
  },finally={
    # nothing
  })


  return(mid_sum_res)
}

############################
# Random effects
#--------------------------
calcRandom_MultiBiFirst_func <- function(x0,userdata,var_notime_model,subject_model,calcFix_data,fix_model){
  # If the user defines a random effect that does not need to be estimated, it is judged by the following code.
  if(!(grepl('1',fix_model[1,'fixRand']))){
    x0[5] = 0
  }
  if(!(grepl('1',fix_model[2,'fixRand']))){
    x0[6] = 0
  }
  if(!(grepl(var_notime_model[1],fix_model[1,'fixRand']))){
    x0[1] = 0
  }
  if(!(grepl(var_notime_model[2],fix_model[1,'fixRand']))){
    x0[2] = 0
  }
  if(!(grepl(var_notime_model[1],fix_model[2,'fixRand']))){
    x0[3] = 0
  }
  if(!(grepl(var_notime_model[2],fix_model[2,'fixRand']))){
    x0[4] = 0
  }
  times <- seq(min(userdata$time),max(userdata$time),length.out=length(unique(userdata[,'time'])))
  parms = c(beta1 = calcFix_data$par[1] + x0[1], # beta1 is etaI1
            beta2 = calcFix_data$par[2] + x0[2], # beta2 is etaI2
            beta3 = calcFix_data$par[3] + x0[3], # beta3 is etaI3
            beta4 = calcFix_data$par[4] + x0[4]) # beta4 is etaI4
  state = c(x = calcFix_data$par[5] + x0[5],
            y = calcFix_data$par[6] + x0[6])
  mid_solve_data <- deSolve::ode(y=state,
                                 times=times,
                                 func=solve_MultiBiFirst_func,
                                 parms=parms)

  mid_solve_data = data.frame(mid_solve_data)
  mid_solve_data = stats::setNames(mid_solve_data,c('time','solver_x','solver_y'))
  times <- seq(min(userdata$time),max(userdata$time),length.out=length(unique(userdata[,'time'])))
  mid_df_res <- data.frame(seq=1:nrow(userdata),
                           time = times,
                           x = userdata[,var_notime_model[1]],
                           y = userdata[,var_notime_model[2]],
                           subject = userdata[,subject_model]
                           # solver_x = NA,
                           # solver_y = NA
                           )
  # estimating residual in loss function
  guess_xy = c(x0[7],x0[8])
  if(x0[7] == 0){
    guess_xy[1] = 0.00001
  }
  if(x0[8] == 0){
    guess_xy[2] = 0.00001
  }
  log_e_x = log(guess_xy[1]**2)
  log_e_y = log(guess_xy[2]**2)
  inverse_ex = 1/(guess_xy[1]**2)
  inverse_ey = 1/(guess_xy[2]**2)

  # When the time is the same, userdata and solverdata will be merged.s
  mid_df_res = merge(mid_df_res,mid_solve_data,by='time',all = TRUE)
  for (col in names(mid_df_res)) {
    mid_df_res[is.na(mid_df_res[[col]]), col] <- 0
  }
  mid_df_res[,'e_x'] = log_e_x + ((mid_df_res[,'x'] - mid_df_res[,'solver_x']) ** 2) * inverse_ex
  mid_df_res[,'e_y'] = log_e_y + ((mid_df_res[,'y'] - mid_df_res[,'solver_y']) ** 2) * inverse_ey
  mid_sum_res = sum(mid_df_res[,'e_y'] + mid_df_res[,'e_x'])

  ## The sum of squared may be infinitely
  if(!(is.finite(mid_sum_res))){
    mid_df_res[,'e_y'] = (mid_df_res[,'y']) ** 2
    mid_df_res[,'e_x'] = (mid_df_res[,'x']) ** 2
    print('waiting')
    mid_sum_res = sum(mid_df_res['e_y'] + mid_df_res['e_x'])
  }
  return(mid_sum_res)
}

############################
# Differential equations
#--------------------------
solve_MultiBiFirst_func <- function(t,state,parms){
  with(as.list(c(state,parms)),{
    dX <- beta1*x + beta2*y
    dY <- beta3*x + beta4*y
    list(c(dX,dY))
  })# end with(as.list ...
}

############################
# calc_MultiBiFirst_United_func
#--------------------------
calc_MultiBiFirst_United_func <- function(userdata,var_model,guess,method,subject_model,fix_model){
  ## identify variables and information
  var_notime_model = var_model[var_model$field != 'time','field']
  var_notime_hat = paste(var_notime_model,'_hat',sep = "")
  var_notime_fixed_hat = paste(var_notime_model,'_fixed_hat',sep = "")
  var_notime_data = var_model[var_model$field != 'time','variable']

  fixed_random_guess <- guess[[1]]
  len_subject = unique(userdata[subject_model])
  len_subject = as.list(len_subject)
  len_subject = length(len_subject[[1]])
  fix_rand_guess <- c(fixed_random_guess,rep(0, len_subject * 6))
  fix_rand_guess <- c(fix_rand_guess,rep(1,len_subject * 2))

  ## fixed and random effect
  calcFix_data = stats::optim(fix_rand_guess,
                              calcFix_MultiBiFirst_United_func,
                              NULL,
                              method = method[[1]],
                              # control = list(reltol = 1e-6),
                              userdata = userdata,
                              var_notime_model = var_notime_model,
                              subject_model =subject_model,
                              hessian=TRUE
                              )
  Neg2LL = calcFix_data$value[[1]][1] + 3.67
  ## random effects
  randEffect_data = data.frame(Subject=NULL,
                               EtaI_1=NULL,
                               EtaI_2=NULL,
                               EtaI_3=NULL,
                               EtaI_4=NULL,
                               EtaI_x=NULL,
                               EtaI_y=NULL)
  rand_li = calcFix_data$par[7:length(calcFix_data$par)]
  mid_subject_name = unique(userdata[subject_model])
  mid_subject_li = as.list(mid_subject_name)
  for(i_subject in seq(0,len_subject-1)){
    randEffect_list = list(Subject=mid_subject_li[[1]][i_subject + 1],
                           EtaI_1=rand_li[1 + i_subject*6],
                           EtaI_2=rand_li[2 + i_subject*6],
                           EtaI_3=rand_li[3 + i_subject*6],
                           EtaI_4=rand_li[4 + i_subject*6],
                           EtaI_x=rand_li[5 + i_subject*6],
                           EtaI_y=rand_li[6 + i_subject*6])
    randEffect_data = rbind(randEffect_data,do.call("data.frame", randEffect_list))
  }
  ## predict
  times = sort(unique(userdata[,'time']),decreasing = FALSE)

  fixed_parms = c(beta1 = calcFix_data$par[1],
                  beta2 = calcFix_data$par[2],
                  beta3 = calcFix_data$par[3],
                  beta4 = calcFix_data$par[4])
  fixed_state = c(x = calcFix_data$par[5],
                  y = calcFix_data$par[6])
  predict_fixed_data = deSolve::ode(y=fixed_state,
                                    times=times,
                                    func=solve_MultiBiFirst_func,
                                    parms=fixed_parms)
  predict_fixed_data = data.frame(predict_fixed_data)
  names(predict_fixed_data) = c('time',var_notime_fixed_hat)

  ## predict random effects data
  predict_data = data.frame(Subject=NULL,
                            time=NULL,
                            x=NULL,
                            y=NULL)
  for (i_subject in 1:NROW(randEffect_data)) {
    i_subject_parameters = unlist(randEffect_data[i_subject,c(2,3,4,5,6,7)])
    parms = c(beta1 = calcFix_data$par[1] + i_subject_parameters[1][[1]],
              beta2 = calcFix_data$par[2] + i_subject_parameters[2][[1]],
              beta3 = calcFix_data$par[3] + i_subject_parameters[3][[1]],
              beta4 = calcFix_data$par[4] + i_subject_parameters[4][[1]])
    state = c(x = calcFix_data$par[5] + i_subject_parameters[5][[1]],
              y = calcFix_data$par[6] + i_subject_parameters[6][[1]])
    mid_predict_data <- deSolve::ode(y=state,
                                     times=times,
                                     func=solve_MultiBiFirst_func,
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
              predict_fixed_data = predict_fixed_data,
              Neg2LL = Neg2LL))
  }
#--------------------------
calcFix_MultiBiFirst_United_func <- function(x0,userdata,var_notime_model,subject_model){
  # times <- sort(unique(userdata[,'time']),decreasing = FALSE)
  times <- seq(min(userdata$time),max(userdata$time),length.out=length(unique(userdata[,'time'])))
  parms = c(beta1 = x0[1],
            beta2 = x0[2],
            beta3 = x0[3],
            beta4 = x0[4])
  state = c(x = x0[5],
            y = x0[6])
  mid_fixed_rand_df = data.frame(rand_beta1 = numeric(),
                                 rand_beta2 = numeric(),
                                 rand_beta3 = numeric(),
                                 rand_beta4 = numeric(),
                                 rand_init1 = numeric(),
                                 rand_init2 = numeric(),
                                 stringsAsFactors = FALSE)
  subjects <- unique(userdata[[subject_model]])
  for (idx in seq_along(subjects)) {
    mid_idx <- idx * 6
    mid_i_key <- data.frame(
      rand_beta1 = x0[mid_idx + 1],
      rand_beta2 = x0[mid_idx + 2],
      rand_beta3 = x0[mid_idx + 3],
      rand_beta4 = x0[mid_idx + 4],
      rand_init1 = x0[mid_idx + 5],
      rand_init2 = x0[mid_idx + 6],
      stringsAsFactors = FALSE
    )
    mid_fixed_rand_df <- rbind(mid_fixed_rand_df, mid_i_key)
  }
  mid_fixed_rand_df[['fixed_beta1']] = x0[1]
  mid_fixed_rand_df[['fixed_beta2']] = x0[2]
  mid_fixed_rand_df[['fixed_beta3']] = x0[3]
  mid_fixed_rand_df[['fixed_beta4']] = x0[4]
  mid_fixed_rand_df[['fixed_init1']] = x0[5]
  mid_fixed_rand_df[['fixed_init2']] = x0[6]
  mid_idx = 1
  res_sum = 0
  e_guess = x0[6 + length(subjects)*6+1:length(x0)]
  e_guess_two = 1
  for(i in subjects){
    mid_parms = c(beta1 = x0[1] + mid_fixed_rand_df[mid_idx,'rand_beta1'],
                  beta2 = x0[2] + mid_fixed_rand_df[mid_idx,'rand_beta2'],
                  beta3 = x0[3] + mid_fixed_rand_df[mid_idx,'rand_beta3'],
                  beta4 = x0[4] + mid_fixed_rand_df[mid_idx,'rand_beta4'])
    mid_state = c(x = x0[5] + mid_fixed_rand_df[mid_idx,'rand_init1'],
                  y = x0[6] + mid_fixed_rand_df[mid_idx,'rand_init2'])
    mid_solve_data <- deSolve::ode(y=mid_state,
                                   times=times,
                                   func=solve_MultiBiFirst_func,
                                   parms=mid_parms)
    mid_df_res = data.frame(mid_solve_data)
    mid_df_res['mid_subject'] = mid_idx
    colnames(mid_df_res)[colnames(mid_df_res) == "x"] <- "solver_x"
    colnames(mid_df_res)[colnames(mid_df_res) == "y"] <- "solver_y"
    mid_userdata = userdata[userdata[subject_model] == i,]
    mid_res_df = merge(mid_userdata, mid_df_res, by = "time", all.x = TRUE)
    solver_li = c('solver_x','solver_y')
    solver_idex = 1
    # estimating residual in loss function
    guess_xy = c(e_guess[e_guess_two],e_guess[e_guess_two + 1])
    e_guess_two = e_guess_two + 2
    if(guess_xy[1] == 0){
      guess_xy[1] = 0.00001
    }
    if(guess_xy[2] == 0){
      guess_xy[2] = 0.00001
    }
    mid_sum_res = 0
    mid_i_num = 1
    for(i in var_notime_model){
      log_e = log(guess_xy[mid_i_num] ** 2)
      inverse_e = 1 / (guess_xy[mid_i_num] ** 2)
      mid_col_name = paste('e_',i,sep = '')
      mid_res_df[mid_col_name] = log_e + ((mid_res_df[i] - mid_res_df[solver_li[solver_idex]]) ** 2) * inverse_e
      mid_sum_res = mid_sum_res + sum(mid_res_df[mid_col_name])
      solver_idex = solver_idex + 1
      mid_i_num = mid_i_num + 1
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
