#' Adjust the model function
#'
#' @param model The original string of users' defined.
#'
#' @return dataframe data.frame(model)
AdjustModel_func <- function(model){
  if (is.null(model)){
    # print('Model is not defined. So the program will fit the default bivariate first-order differential equation. And the field of your data must contain "seq","time","x","y".')
    model = '
            x =~ x
            y =~ y
            time =~ time
            x(1) ~ x + y
            y(1) ~ y + x
    '
  }
  modelDF <- data.frame()
  # The model syntax: delate the useless \t\r\n and split field the by \n
  model <- trimws(model,which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
  modellist <- strsplit(model,'\n')
  # Split the variable and data's field. Save them as dataframe
  modelDF <- data.frame()
  modelunlist = unlist(modellist)
  for (i in seq(along = modelunlist)){
    if (grepl('=~',modelunlist[i])){ #if the syntax contain "=~"
      mid_operator <- unlist(strsplit(modelunlist[i],'=~'))
      modelDF[i,'field'] <- trimws(mid_operator[1],which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
      modelDF[i,'operator'] <- '=~'
      modelDF[i,'variable'] <- trimws(mid_operator[2],which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
    }else if (grepl('~~',modelunlist[i])){ #if the syntax contain "~~"
      mid_operator <- unlist(strsplit(modelunlist[i],'~~'))
      modelDF[i,'field'] <- mid_operator[1]
      modelDF[i,'operator'] <- '~~'
      modelDF[i,'variable'] <- mid_operator[2]
    }else if (grepl('~',modelunlist[i])){ #if the syntax contain "~"
      mid_operator <- unlist(strsplit(modelunlist[i],'~'))
      modelDF[i,'field'] <- mid_operator[1]
      modelDF[i,'operator'] <- '~'
      modelDF[i,'variable'] <- mid_operator[2]
    }
  }

  # The operator are not contained in dataframe that will not display.
  modelDF <- modelDF[which(modelDF$field != '<NA>'),]
  return(modelDF)
}
