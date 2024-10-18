
#' Calculate R-squared RMSE SE and so on.
#'
#' @param solve_data a list data, that is a solve of differential equations.
#' @param userdata a data frame containing all model variables. The "time" column must be included.
#' @param predict predict data.
#' @param var_model a dataframe containing equations.
#' @param table outtable.
#'
#' @return a list
Info_func <- function(solve_data,userdata,predict,var_model,table){

  var_notime_data = var_model[var_model$field != 'time','field']
  var_notime_hat = paste(var_notime_data,"_hat",sep = "")
  rsquared_data <- Rsquared_func(userdata,predict,var_notime_data,var_notime_hat)
  RMSE_data <- RMSE_func(userdata,predict,var_notime_data,var_notime_hat)
  Hessian_BiFirst_data <- Hessian_BiFirst_func(solve_data[[1]]$hessian)
  SE_vector = c()
  for (i in 1:NROW(Hessian_BiFirst_data)) {
    SE_vector[i] = Hessian_BiFirst_data[i,i]
  }
  ## Identify univariate second-order. Changeing the SE
  if(length(var_notime_data) == 1){
    SE_vector = c(SE_vector[1],SE_vector[2],SE_vector[3],SE_vector[4])
  }else{
    SE_vector = c(SE_vector[1],SE_vector[2],SE_vector[3],SE_vector[4],SE_vector[5],SE_vector[6])
  }
  table[,'SE'] = SE_vector
  return(list(rsquared_data = rsquared_data,
              RMSE = RMSE_data,
              SE = SE_vector,
              table = table))
}

############################
# R-squared
#--------------------------
Rsquared_func <- function(userdata,predict,var_notime_data,var_notime_hat){
  message('Estimating R_squared')
  # merge_df = merge(userdata,predict,by.x=)
  cor_list = c()
  for (i in 1:length(var_notime_data)) {
    mid_cor = stats::cor(userdata[,var_notime_data[i]],predict[,var_notime_hat[i]])
    cor_list = c(cor_list,mid_cor)
  }
  cor_list = cor_list **2
}

############################
# RMSE
#--------------------------
RMSE_func <- function(userdata,predict,var_notime_data,var_notime_hat){
  message('Estimating RMSE')
  RMSE_res = c()
  for (i in 1:length(var_notime_hat)) {
    SSxe <- sum((userdata[,var_notime_data[i]] - predict[,var_notime_hat[i]])**2)
    RMSEx <- sqrt(SSxe / NROW(userdata[,'time']))
    RMSE_res[i] = RMSEx
  }
  return(RMSE_res)
}

############################
# Hessian
#--------------------------
Hessian_BiFirst_func <- function(hessian){
  message('Estimating Hessian')
  hessian[hessian<0] <- NaN
  mid_SE <- sqrt(1/hessian)
  return(mid_SE)
}
