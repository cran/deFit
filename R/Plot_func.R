#' Draw the diagram of differential equation
#'
#' @param userdata a data frame containing all model variables. The "time" column must be included.
#' @param predict predict data.
#' @param modelDF a dataframe of full model.
#' @param var_model a dataframe containing equations.
#' @param field_model the user's data columns.
#' @param order_model N-order differential equation.
#' @param multi_model TRUE or FALSE
#'
#' @return plot
Plot_func <- function(userdata,predict,modelDF,var_model,field_model,order_model,multi_model){
  ## identify variables and information
  names(userdata)[names(userdata) == var_model[var_model$field == 'time','variable']] <- 'time'
  var_notime_model = var_model[var_model$field != 'time','field']
  var_notime_hat = paste(var_notime_model,'_hat',sep = "")
  var_notime_data = var_model[var_model$field != 'time','variable']

  combine_data = cbind(userdata,predict)
  combine_data <- combine_data[, !duplicated(colnames(combine_data))]

  if(length(field_model) == 2){
    if(order_model == 1){
      if(multi_model){
        # multilevel bivariate first-order differential equation
        outplot = ggplot2::ggplot() +
          ggplot2::geom_line(data=combine_data,ggplot2::aes(x=time,y=combine_data[,var_notime_model[1]]),color='red', linetype = "dashed") +
          ggplot2::geom_line(data=combine_data,ggplot2::aes(x=time,y=combine_data[,var_notime_model[2]]),color='green', linetype = "dashed") +
          ggplot2::geom_line(data=combine_data,ggplot2::aes(x=time,y=combine_data[,var_notime_hat[1]]),color='red') +
          ggplot2::geom_line(data=combine_data,ggplot2::aes(x=time,y=combine_data[,var_notime_hat[2]]),color='green') +
          ggplot2::facet_wrap(.~Subject)+
          ggplot2::labs(y = "values",
                        # title = 'Multilevel Bivariate first-order differential equation',
                        caption =  paste("Raw data (Dashed Lines)  & Predict values (Solid Lines)\n",
                                         var_model[1,'field'],'(Red) & ',var_model[2,'field'],'(Green)',
                                         sep=''))
      }else{
        outplot = ggplot2::ggplot() +
          ggplot2::geom_line(data=combine_data,ggplot2::aes(x=time,y=combine_data[,var_notime_model[1]]),color='red', linetype = "dashed") +
          ggplot2::geom_line(data=combine_data,ggplot2::aes(x=time,y=combine_data[,var_notime_model[2]]),color='green', linetype = "dashed") +
          ggplot2::geom_line(data=combine_data,ggplot2::aes(x=time,y=combine_data[,var_notime_hat[1]]),color='red') +
          ggplot2::geom_line(data=combine_data,ggplot2::aes(x=time,y=combine_data[,var_notime_hat[2]]),color='green') +
          ggplot2::labs(y = "values",
                        # title = 'Bivariate first-order differential equation',
                        caption =  paste("Raw data (Dashed Lines)  & Predict values (Solid Lines)\n",
                                         var_model[1,'field'],'(Red) & ',var_model[2,'field'],'(Green)',
                                         sep=''))
      }
    }
  }else if(length(field_model) == 1){
    if(order_model == 2){
      if(multi_model){
        cat('Plotting')
        outplot = ggplot2::ggplot() +
          ggplot2::geom_line(data=combine_data,ggplot2::aes(x=time,y=combine_data[,var_notime_model[1]]),color=2, linetype = "dashed") +
          ggplot2::geom_line(data=combine_data,ggplot2::aes(x=time,y=combine_data[,var_notime_hat[1]]),color=2)  +
          ggplot2::facet_wrap(.~Subject)+
          ggplot2::labs(y = "values",
                        # title = 'Univariate Second order differential equation',
                        caption =  paste("Raw Data (Dashed Lines)  & Predict Values (Solid Lines) \n",
                                         var_notime_model[1],'(Red)',
                                         sep=''))
      }else{
        outplot = ggplot2::ggplot() +
          ggplot2::geom_line(data=combine_data,ggplot2::aes(x=time,y=combine_data[,var_notime_model[1]]),color='red', linetype = "dashed") +
          ggplot2::geom_line(data=combine_data,ggplot2::aes(x=time,y=combine_data[,var_notime_hat[1]]),color='red') +
          ggplot2::labs(y = "values",
                        # title = 'Bivariate first-order differential equation',
                        caption =  paste("Raw data (Dashed Lines)  & Predict values (Solid Lines)\n",
                                         var_model[1,'field'],'(Red) & ',
                                         sep=''))
      }
    }
  }else{
    stop("********************\n",
         "Someting error\n",
         "Contact me",
         "********************")
  }
  print(outplot)
}
utils::globalVariables(
  c("time")
)
