#' Plot all data
#'
#' @param userdata Users's data.
#' @param calcdata Predictor's data.
#'
#' @return The plot of all data.
PlotDe_func <- function(userdata,calcdata){
  PlotDataFrame = merge(userdata,calcdata,by='time')
  #remove 'seq'
  PlotDF = PlotDataFrame[,which(names(PlotDataFrame) != 'seq')]
  ###################################################################
  # differential equational
  # ----------------------------------
  # Calculate how many variable of data
  columnsVariable = (length(names(PlotDF)) - 1) / 2
  limy_data = PlotDF[,which(names(PlotDF) != 'time')]
  limy_max = max(limy_data)
  limy_min = min(limy_data)
  if (columnsVariable == 1){
    # Univariable Second order differential equation
    outplot = ggplot2::ggplot(PlotDF,ggplot2::aes(x=time)) +
      ggplot2::geom_point(y =  PlotDF[,2],color= 2) +
      ggplot2::geom_line(y =  PlotDF[,3],color = 3) +
      ggplot2::labs(y = "y",color="group",title = 'Univariate second-order differential equation') +
      ggplot2::coord_cartesian(ylim = c(limy_min,limy_max))
    outplot = outplot +
      ggplot2::annotate("text",x = max(PlotDF[,'time'])/2, y = limy_max*0.9, color=2, label = paste(names(PlotDF)[2],':point'))+
      ggplot2::annotate("text",x = max(PlotDF[,'time'])/2, y = limy_max*0.96, color=3, label = paste(names(PlotDF)[3],':line'))
  }
  if (columnsVariable == 2){
    # Binary First order differential equation
    outplot = ggplot2::ggplot(PlotDF,ggplot2::aes(x=time)) +
      ggplot2::geom_point(y = PlotDF[,2],color= 2) +
      ggplot2::geom_point(y = PlotDF[,3],color = 3) +
      ggplot2::geom_line(y = PlotDF[,4],color = 4) +
      ggplot2::geom_line(y = PlotDF[,5],color = 5) +
      ggplot2::labs(y = "y",group="color",title = 'Bivariate first-order differential equation') +
      # ggplot2::scale_y_continuous(name = "y",
      #                             sec.axis = ggplot2::sec_axis( trans=~.*1, name="y2"))+
      ggplot2::coord_cartesian(ylim = c(limy_min,limy_max))
    outplot = outplot +
      ggplot2::annotate("text", x = max(PlotDF[,'time'])/2, y = limy_max*0.83, color=2, label = paste(names(PlotDF)[2],':point'))+
      ggplot2::annotate("text", x = max(PlotDF[,'time'])/2, y = limy_max*0.88, color=3, label = paste(names(PlotDF)[3],":point"))+
      ggplot2::annotate("text", x = max(PlotDF[,'time'])/2, y = limy_max*0.93, color=4, label = paste(names(PlotDF)[4],":line"))+
      ggplot2::annotate("text", x = max(PlotDF[,'time'])/2, y = limy_max*0.98, color=5, label = paste(names(PlotDF)[5],":line"))

  }
  print(outplot)

}
utils::globalVariables(
  c("time")
)
