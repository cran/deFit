
#' Title
#'
#' @param userdata users' data
#' @param model a string specifying the model to be used. The "=~" operator is used to define variables, with the name of the variable user defined on the left and the name of the variable in the data on the right. The '~' operator specifies a differential equation, with the dependent variable on the left and the independent variables on the right. See also ‘Details’.
#' @param center TRUE or FALSE
#' @param scale TRUE or FALSE
#'
#' @return dataframe
#' @export scale_within
#'
#' @examples
#' #eg1.
#' data('example3')
#' multi_model <- '
#'   X =~ current
#'   time =~ myTime
#'   X(2) ~ X(1) + X + (1 + X(1) + X | year)
#'   '
#' scale_mydata <- scale_within(example3[(example3["year"] >= 2015)&(example3["year"] <= 2018),]
#' ,multi_model
#' ,center=TRUE)
scale_within <- function(userdata,model=NA,center=FALSE,scale=FALSE){
  Adjust_model = AdjustModel_func(model)
  mid_data = userdata
  mid_ColumnsNames = names(userdata)
  # Get the variables of the model
  var_model = Adjust_model
  var_model = var_model[which(var_model['operator']=='=~'),]
  var_model = var_model[(which(var_model['field'] != 'time')),]
  sub_model = Adjust_model
  sub_model = sub_model[which(sub_model['operator'] == '~'),]

  # if there is no model. This function will scale all columns.
  if(is.na(model)){
    cat('The model is not defined, all columns will be scaled.\n')
    cat('The NA will not be filled in with 0.')
    for(i in mid_ColumnsNames){
      mid_data[paste('scale_',i,sep='')] = userdata[,i] - mean(userdata[,i])
      #Centering is the data minus the average of the data.
    }
  } else if (all(is.na(Adjust_model[,'subject']))){
    message('There is no "subject" in your model')
    print(mid_ColumnsNames)
    for(i in var_model[,'variable']){
      mid_everyItem <- userdata[,i]
      mid_data[,i] = mid_data[,i] - mean(mid_everyItem)
    }
  } else if (nrow(sub_model)==1){
    sub_model = sub_model[,'subject']
    for(i in var_model[,'variable']){
      # print(unique(userdata[,sub_model[1]]))
      for(i_subject in unique(userdata[,sub_model[1]])){
        mid_everyItem <- userdata[which(userdata[,sub_model[1]]==i_subject),i]
        mid_data[which(userdata[,sub_model[1]]==i_subject),i] = mid_data[which(userdata[,sub_model[1]]==i_subject),i] - mean(mid_everyItem)
      }
    }
  }else if (nrow(sub_model)==2){
    sub_model = sub_model[,'subject']
    # sub_model = gsub("\\s+", "", sub_model[,'subject'])

    # Determine if the subject is the same.
    if( sub_model[1] != sub_model[2]){
      warning('The subject must be kept the same.')
    }
    for(i in var_model[,'variable']){
      # print(unique(userdata[,sub_model[1]]))
      for(i_subject in unique(userdata[,sub_model[1]])){
        mid_everyItem <- userdata[which(userdata[,sub_model[1]]==i_subject),i]
        mid_data[which(userdata[,sub_model[1]]==i_subject),i] = mid_data[which(userdata[,sub_model[1]]==i_subject),i] - mean(mid_everyItem)
      }
    }
  }else{
    message('There are some errors')
  }
  return(mid_data)
}
