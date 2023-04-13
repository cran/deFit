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
  modelDF <- data.frame(field=NULL,
                        operator=NULL,
                        variable=NULL,
                        fixRand=NULL,
                        subject=NULL)
  modelunlist = unlist(modellist)
  for (i in seq(along = modelunlist)){
    if (grepl('=~',modelunlist[i])){ #if the syntax contain "=~"
      mid_operator <- unlist(strsplit(modelunlist[i],'=~'))
      modelDF[i,'field'] <- trimws(mid_operator[1],which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
      modelDF[i,'operator'] <- '=~'
      modelDF[i,'variable'] <- trimws(mid_operator[2],which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
      modelDF[i,'fixRand'] <- NA
      modelDF[i,'subject'] <- NA
    }else if (grepl('~~',modelunlist[i])){ #if the syntax contain "~~"
      mid_operator <- unlist(strsplit(modelunlist[i],'~~'))
      modelDF[i,'field'] <- trimws(mid_operator[1],which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
      modelDF[i,'operator'] <- '~~'
      modelDF[i,'variable'] <- trimws(mid_operator[2],which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
      modelDF[i,'fixRand'] <- NA
      modelDF[i,'subject'] <- NA
    }else if (grepl('~',modelunlist[i])){ #if the syntax contain "~"
      mid_operator <- unlist(strsplit(modelunlist[i],'~'))
      mid_splitoperator <- strsplit(mid_operator[2],"\\+\\s*\\(")
      mid_fixRand <- unlist(mid_splitoperator)
      mid_subject <- unlist(strsplit(mid_fixRand[2],'\\|'))
      modelDF[i,'field'] <- trimws(mid_operator[1],which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
      modelDF[i,'operator'] <- '~'
      modelDF[i,'variable'] <- trimws(mid_fixRand[1],which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
      modelDF[i,'fixRand'] <- trimws(mid_subject[1],which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
      modelDF[i,'subject'] <- trimws(mid_subject[2],which = c("both", "left", "right"), whitespace = "[ \t\r\n)]")
    }
  }

  # The operator are not contained in dataframe that will not display.
  modelDF <- modelDF[which(modelDF$field != '<NA>'),]
  return(modelDF)
}
