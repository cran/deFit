
#' Center the data according to model
#'
#' @param userdata users' data
#' @param model a string specifying the model to be used. The "=~" operator is used to define variables, with the name of the variable user defined on the left and the name of the variable in the data on the right. The '~' operator specifies a differential equation, with the dependent variable on the left and the independent variables on the right. See also ‘Details’.
#' @param center TRUE or FALSE
#' @param scale TRUE or FALSE
#'
#' @return dataframe
#' @export Scale_within
#'
#' @examples
#' #eg1.
#' data('example3')
#' multi_model <- '
#'   X =~ current
#'   time =~ myTime
#'   X(2) ~ X(1) + X + (1 + X(1) + X | year)
#'   '
#' scale_mydata <- Scale_within(example3[(example3["year"] >= 2015)&(example3["year"] <= 2018),]
#' ,multi_model
#' ,center=TRUE)
Scale_within <- function(userdata,model=NA,center=FALSE,scale=FALSE){
  # Get the variables of the model
  init_list = Init_func(userdata,model,guess=NA,method=NA)
  subject_model = unique(init_list$subject_model)
  var_userdata = init_list$var_model[,'variable']
  var_model_new = init_list$var_model[,'field']
  var_model_notime = var_model_new[var_model_new!= 'time']
  var_notime_userdata = init_list$var_model[init_list$var_model[,'field'] != 'time','variable']

  ## add the columns of model in userdata
  for (i in 1:length(var_userdata)) {
    userdata[,var_model_new[i]] = userdata[,var_userdata[i]]
  }

  ## the subject data subtract the means of same subject
  uni_subject = unique(userdata[,subject_model])
  for (i_subject in uni_subject) {
    for (i_var in var_model_notime) {
      mid_mean = mean(userdata[userdata[,subject_model] == i_subject,i_var])
      mid_center_data = userdata[userdata[,subject_model] == i_subject,i_var] - mid_mean
      userdata[userdata[,subject_model] == i_subject,i_var] = mid_center_data
    }
  }
  return(userdata)
}
