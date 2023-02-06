#' Judge the model belongs of N-variable and N-order differential equation.
#'
#' @param model The original string users defined.
#'
#' @return c(2,1): Bivariate first-order DE; c(1,2): Univariable second-order DE.
JudgeModel_func <- function(model){
  varnum <- NROW(model[which(model$operator == '=~'),]) #the variables of differential equation
  varnum <- varnum - 1 # The 'time' occupy one field,  so need to reduce 1
  orderDF = model[which(model$operator == '~'),]
  ordernum <- orderDF['field']
  ordernum <- as.numeric(gsub("\\D", "", ordernum))
  ordernum <- strsplit(as.character(ordernum),"")
  ordernum <- as.numeric(unlist(ordernum)) # the order of derivative
  ordernum <- max(ordernum)
  #print(c(varnum,ordernum))
  return(c(varnum,ordernum))
}
