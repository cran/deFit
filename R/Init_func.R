#' Initialize model Judgement variables and so on.
#'
#' @param userdata a data frame containing all model variables. The "time" column must be included.
#' @param model a string of model.
#' @param guess a list or a string. Guess the coefficients or initial values.
#' @param method a list or a string. The available options are 'Nelder-Mead','L-BFGS-B','SANN' and 'BFGS'.
#' @param plot TRUE or FALSE.
#'
#' @return a list
Init_func <- function(userdata,model,guess,method,plot){
  ############################
  # Model
  #--------------------------
  ## Judgement Model
  if (is.null(model)){
    message('Model can not be NULL.')
    model = '
            x =~ x
            y =~ y
            time =~ time
            x(1) ~ x + y
            y(1) ~ y + x
    '
  }
  # Model standardization
  # The model syntax: delete the useless \t\r\n and split field the by \n
  model_str <- model
  model <- trimws(model,which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
  modellist <- strsplit(model,'\n')
  # Split the variable and data's field. Save them as dataframe
  modelDF <- data.frame(field=NULL,
                        order = NULL,
                        operator=NULL,
                        variable=NULL,
                        fixRand=NULL,
                        subject=NULL)
  modelunlist = unlist(modellist)
  for (i in seq(along = modelunlist)){
    if (grepl('=~',modelunlist[i])){ #if the syntax contain "=~"
      mid_operator <- unlist(strsplit(modelunlist[i],'=~'))
      modelDF[i,'field'] <- trimws(mid_operator[1],which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
      modelDF[i,'order'] <- NA
      modelDF[i,'operator'] <- '=~'
      modelDF[i,'variable'] <- trimws(mid_operator[2],which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
      modelDF[i,'fixRand'] <- NA
      modelDF[i,'subject'] <- NA
    }else if (grepl('~~',modelunlist[i])){ #if the syntax contain "~~"
      mid_operator <- unlist(strsplit(modelunlist[i],'~~'))
      modelDF[i,'field'] <- trimws(mid_operator[1],which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
      modelDF[i,'order'] <- NA
      modelDF[i,'operator'] <- '~~'
      modelDF[i,'variable'] <- trimws(mid_operator[2],which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
      modelDF[i,'fixRand'] <- NA
      modelDF[i,'subject'] <- NA
    }else if (grepl('~',modelunlist[i])){ #if the syntax contain "~"
      mid_operator <- unlist(strsplit(modelunlist[i],'~'))
      mid_order <- unlist(strsplit(mid_operator[1],")"))
      mid_order <- unlist(strsplit(mid_order[1],"\\("))
      mid_splitoperator <- strsplit(mid_operator[2],"\\+\\s*\\(")
      mid_fixRand <- unlist(mid_splitoperator)
      mid_subject <- unlist(strsplit(mid_fixRand[2],'\\|'))
      modelDF[i,'field'] <- trimws(mid_operator[1],which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
      modelDF[i,'order'] <- mid_order[2]
      modelDF[i,'operator'] <- '~'
      modelDF[i,'variable'] <- trimws(mid_fixRand[1],which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
      modelDF[i,'fixRand'] <- trimws(mid_subject[1],which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
      modelDF[i,'subject'] <- trimws(mid_subject[2],which = c("both", "left", "right"), whitespace = "[ \t\r\n)]")
    }
  }
  # The operator are not contained in dataframe that will not display.
  modelDF <- modelDF[which(modelDF$field != '<NA>'),]

  # Model information
  var_model = modelDF[modelDF$operator == '=~',]
  var_notime_model = var_model[var_model$field != 'time','variable']
  var_data = var_model[,'variable']
  multi_model = !all(is.na(modelDF[,'fixRand'])) #multi_model is true. Multilevel DE model is required.
  field_model = modelDF[modelDF$operator == '~','field']
  order_model = max(modelDF[modelDF$operator == '~','order'])
  subject_model = modelDF[modelDF$operator == '~','subject']

  ############################
  # userdata
  #--------------------------
  ## Judgement userdata
  colname_data = names(userdata)
  if (!all(var_data %in% colname_data)){
    stop("********************\n",
         "Model variable must be included in your data.\n",
         "********************")
  }
  ## Judgement time
  # must be contain time variable
  if(!any('time' %in% var_model[,'field'])){
    stop("********************\n",
         "Model must contain 'time'.\n",
         "********************")
  }
  ## Judgement multilevel
  # subject will be covert to string
  if (isTRUE(multi_model)){
    userdata[[subject_model[1]]] = as.character(userdata[[subject_model[1]]])
  }

  ############################
  # guess
  #--------------------------
  ## guess values can be a list or a string. When it is a list, the model must be multilevel model.
  ## the guess values now is support string, list(multilevel), dict(bayesian).
  if(is.null(guess)){
    # if guess values is null, Then, we need to give guess values base on model.
    if(isTRUE(multi_model)){
      # multi_model is true
      if(length(field_model) == 2){
        # multilevel bivariate
        if(order_model == 1){
          #first-order
          guess = list(c(0.0001,0.0001,0.0001,0.0001,0.0001,0.0001),
                       c(0,0,0,0,0,0))
        }
      }else if(length(field_model) == 1){
        if(order_model == 2){
          #seond-order
          # multilevel univariate second-order
          guess = list(c(0.0001,0.0001,0.0001,0.0001),
                       c(0,0,0,0))
        }
      }
    }else{
      # multi_model is false.
      if(length(field_model) == 2){
        # bivariate
        if(order_model == 1){
          #first-order
          guess = list(c(0.0001,0.0001,0.0001,0.0001,userdata[[var_notime_model[1]]][1],userdata[[var_notime_model[1]]][2]))
          # guess2 = NULL
        }
      }else{
        # univariate second-order
        guess = list(c(0.0001,0.0001,userdata[[var_notime_model[1]]][1],0.0001))
      }
    }
  }else{
    if(!is.list(guess)){
      if(isTRUE(multi_model)){
        ## if is multilevel the guess must contain two step
        guess = list(guess,guess)
      }
      guess = list(guess)
    }
  }
  ############################
  # method
  #--------------------------
  if(is.null(method)){
    # if guess values is null, Then, we need to give guess values base on model.
    if(isTRUE(multi_model)){
      # multilevel
      if(length(field_model) == 2){
        # bivariate
        if(order_model == 1){
          # multilevel bivariate first-order
          method = list('L-BFGS-B','L-BFGS-B')
        }
      }else{
        # multilevel univariate second-order
        method = list('L-BFGS-B','L-BFGS-B')
      }
    }else if(!isTRUE(multi_model)){
      if(length(field_model) == 2){
        # bivariate
        if(order_model == 1){
          #first-order
          method = c('L-BFGS-B')
        }
      }else if(length(field_model) == 1){
        # univariate
        if(order_model == 2){
          # univariate second order
          method = c('L-BFGS-B')
        }
      }
    }else{
      stop("********************\n",
           "Error 001\n",
           "********************")
    }
  }else{
    if(!is.list(method)){
      if(isTRUE(multi_model)){
        ## if is multilevel the guess must contain two step
        method = list(method,method)
      }
    }
  }

  ############################
  # multilevel
  #--------------------------
  if(isTRUE(multi_model)){
    # subject_model
    #--------------------------
    subject_model = unique(subject_model)
    # if(all(subject_model == subject_model[1])){
    #   subject_model = unique(subject_model)
    # }else{
    #   stop("********************\n",
    #        "ALL subject must be the same.\n",
    #        "********************")
    # }

    name_subject = modelDF[modelDF$operator == '~','subject']
    if(!all(unique(name_subject) %in% names(userdata))){
      stop("********************\n",
           "Subject must be contained in columns of your data.\n",
           "********************")
    }
  }

  ############################
  # return
  #--------------------------
  return(list(userdata=userdata,
              modelDF=modelDF,
              field_model=field_model,
              multi_model=multi_model,
              order_model=order_model,
              var_model = var_model,
              guess=guess,
              method=method,
              subject_model = subject_model))
}
