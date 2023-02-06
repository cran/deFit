
IntiOption_func <- function(userdata,model,choosemodel,guess=NULL,method=NULL){

  # standard the userdata by model
  # judge the userdata whether include 'time' and variable
  data <- data.frame()
  usersTime <- model[which(model$field == 'time'),'variable']
  usersVariable <- model[which(model$operator == '=~'),]
  usersVariable <- usersVariable[-which(usersVariable$field== 'time'),]
  data[seq(1:NROW(userdata)),'seq'] <- seq(1:NROW(userdata))
  if(!(usersTime %in% names(userdata))){
    # print(paste('"',usersTime,'"is not found in columns of data.'))
    stop(paste('"',usersTime,'"is not found in columns of data.'))
  }
  data[,'time'] <- userdata[,usersTime]

  if(!all(c(usersVariable$variable) %in% c(names(userdata)))){
    # print(paste('"',usersVariable$variable,'"are not found in columns of data.'))
    stop(paste('"',usersVariable$variable,'"are not found in columns of data.'))
  }
  data[,usersVariable$field] <- userdata[,usersVariable$variable]

  # choose method by choosemodel
  if(is.null(method)){
    if(all(choosemodel == c(1,2))){
      method = 'L-BFGS-B'
    }else if(all(choosemodel == c(2,1))){
      method = 'Nelder-Mead'
    }else{
      method = 'BFGS'
    }
  }

  # keep to modify guess values
  if (is.null(guess)){
   # if the values of guess are NA, to modify in Solver_***_function
    guess <- c(NA,NA,NA,NA,NA,NA)
  }else{
    # message('The values of your guess are')
    # message(guess)
  }

  reslist <- list(data= data,choosemodel=choosemodel,guess=guess,method=method)
  return(reslist)
}
