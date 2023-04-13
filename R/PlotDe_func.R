#' Plot all data
#'
#' @param predictor The data of differential equation model
#' @param model The model of standardization
#' @param chooseModel c(2,1) or c(2,1)
#' @param userdata User's data.
#'
#' @return The plot of all data.
PlotDe_func <- function(userdata,predictor,model,chooseModel){

  mid_predict_data = predictor

  var_model = model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]
  limy_solvervar <- c(paste(var_model[,'field'],'_hat',sep=''))
  limy_var <- c(var_model[,'field'])
  limy_min <- min(mid_predict_data[,c(limy_var,limy_solvervar)])
  limy_max <- max(mid_predict_data[,c(limy_var,limy_solvervar)])
  times = sort(unique(mid_predict_data[,'time']),decreasing = FALSE)
  if (all(chooseModel == c(1,2))){
    # Univariable Second order differential equation
    if(all(is.na(model[,'fixRand']) & all(is.na(model[,'subject'])))){
      # cat('Random effects and fixed effects are not defined\n')
      cat('Plotting')
      outplot = ggplot2::ggplot() +
        ggplot2::geom_line(data=mid_predict_data,ggplot2::aes(x=time,y=mid_predict_data[,var_model[1,'field']]),color=2, linetype = "dashed") +
        ggplot2::geom_line(data=mid_predict_data,ggplot2::aes(x=time,y=mid_predict_data[,paste(var_model[1,'field'],'_hat',sep='')]),color=2) +
        ggplot2::labs(y = "values",
                      # title = 'Univariable Second order differential equation',
                      caption =  paste("Raw Data (Dashed Lines)  & Predict Values (Solid Lines) \n",
                                       var_model[1,'field'],'(Red)',
                                       sep='')) +
        ggplot2::coord_cartesian(ylim = c(limy_min,
                                          limy_max))
    }else{
      # cat('Random effects and fixed effects have been defined\n')
      cat('Plotting')
      outplot = ggplot2::ggplot() +
        ggplot2::geom_line(data=mid_predict_data,ggplot2::aes(x=time,y=mid_predict_data[,var_model[1,'field']]),color=2, linetype = "dashed") +
        ggplot2::geom_line(data=mid_predict_data,ggplot2::aes(x=time,y=mid_predict_data[,paste(var_model[1,'field'],'_hat',sep='')]),color=2) +
        ggplot2::facet_wrap(.~subject) +
        ggplot2::labs(y = "values",
                      # title = 'Multilevel Univariable Second order differential equation',
                      caption =  paste("Raw data (Dashed Lines) & Predict values (Solid Lines)\n",
                                       var_model[1,'field'],'(Red)',
                                       sep='')) +
        ggplot2::coord_cartesian(ylim = c(limy_min,
                                          limy_max))
    }

  }else if (all(chooseModel == c(2,1))){
    # Binary First order differential equation
    # print('Bivariate first-order differential equation')
    if(all(is.na(model[,'fixRand']) & all(is.na(model[,'subject'])))){
      # cat('Random effects and fixed effects are not defined\n')
      outplot = ggplot2::ggplot() +
        ggplot2::geom_line(data=mid_predict_data,ggplot2::aes(x=time,y=mid_predict_data[,var_model[1,'field']]),color='red', linetype = "dashed") +
        ggplot2::geom_line(data=mid_predict_data,ggplot2::aes(x=time,y=mid_predict_data[,var_model[2,'field']]),color='green', linetype = "dashed") +
        ggplot2::geom_line(data=mid_predict_data,ggplot2::aes(x=time,y=mid_predict_data[,paste(var_model[1,'field'],'_hat',sep='')]),color='red') +
        ggplot2::geom_line(data=mid_predict_data,ggplot2::aes(x=time,y=mid_predict_data[,paste(var_model[2,'field'],'_hat',sep='')]),color='green') +
        ggplot2::labs(y = "values",
                      # title = 'Bivariate first-order differential equation',
                      caption =  paste("Raw data (Dashed Lines)  & Predict values (Solid Lines)\n",
                                       var_model[1,'field'],'(Red) & ',var_model[2,'field'],'(Green)',
                                       sep='')) +
        ggplot2::coord_cartesian(ylim = c(limy_min,
                                          limy_max))
    }else{
      # cat('Random effects and fixed effects have been defined\n')
      cat('Plotting')
      outplot = ggplot2::ggplot() +
        ggplot2::geom_line(data=mid_predict_data,ggplot2::aes(x=time,y=mid_predict_data[,var_model[1,'field']]),color=2, linetype = "dashed") +
        ggplot2::geom_line(data=mid_predict_data,ggplot2::aes(x=time,y=mid_predict_data[,var_model[2,'field']]),color=3, linetype = "dashed") +
        ggplot2::geom_line(data=mid_predict_data,ggplot2::aes(x=time,y=mid_predict_data[,paste(var_model[1,'field'],'_hat',sep='')]),color=2) +
        ggplot2::geom_line(data=mid_predict_data,ggplot2::aes(x=time,y=mid_predict_data[,paste(var_model[2,'field'],'_hat',sep='')]),color=3) +
        ggplot2::facet_wrap(.~subject) +
        ggplot2::labs(y = "values",
                      # title = 'Multilevel Bivariate first-order differential equation',
                      caption =  paste("Raw Data (Dashed Lines)  & Predict Values (Solid Lines)\n",
                                       var_model[1,'field'],'(Red) & ',var_model[2,'field'],'(Green)',
                                       sep='')) +
        ggplot2::coord_cartesian(ylim = c(limy_min,
                                          limy_max))
    }

  }
  else{
    cat('Plotting')
    outplot = ggplot2::ggplot() +
      ggplot2::geom_line(data=mid_predict_data,ggplot2::aes(x=time,y=mid_predict_data[,var_model[1,'field']]),color=2, linetype = "dashed") +
      ggplot2::geom_line(data=mid_predict_data,ggplot2::aes(x=time,y=mid_predict_data[,var_model[2,'field']]),color=3, linetype = "dashed") +
      ggplot2::geom_line(data=mid_predict_data,ggplot2::aes(x=time,y=mid_predict_data[,paste(var_model[1,'field'],'_hat',sep='')]),color=4) +
      ggplot2::geom_line(data=mid_predict_data,ggplot2::aes(x=time,y=mid_predict_data[,paste(var_model[2,'field'],'_hat',sep='')]),color=6) +
      ggplot2::facet_wrap(.~subject) +
      ggplot2::labs(y = "values",
                    # title = 'Multilevel Bivariate first-order differential equation',
                    caption =  paste("Raw Data (Dashed Lines)  & Predict Values (Solid Lines)\n",
                                     var_model[1,'field'],'(Red) &',var_model[2,'field'],'(Green)',
                                     sep='')) +
      ggplot2::coord_cartesian(ylim = c(limy_min,
                                        limy_max))
  }
  print(outplot)

}
utils::globalVariables(
  c("time")
)

