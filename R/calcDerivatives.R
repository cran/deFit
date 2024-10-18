
#' Calculating the Derivatives
#'
#' @param data a data frame.
#' @param column names of variables in the long format that correspond to multiple variables in the wide format.
#' @param groupby Character vector. Only used if the data is in a data.frame.
#' @param time A variable name in the data frame containing sampling time information.
#' @param order integer scalar.
#' @param window integer scalar.Must be an odd number
#'
#' @return a data frame contains derivatives.
#' @export calcDerivatives
#'
#' @details examples
#' \preformatted{
#' #eg1.
#' derivatives1 <- calcDerivatives(data=example3,column='expected',groupby='year',order=2,window=5,interval=1)
#' #eg2.
#' derivatives2 <- calcDerivatives(data=example3,column=c('expected','current'),groupby='year',time='myTime',order=2,window=5,interval=1)
#' #eg3.
#' derivativese3 <- calcDerivatives(data=example3,column=c('expected','current'),groupby='year',order=2,window=5,interval=1)
#'
#' }
#' @md
calcDerivatives <- function(data,column,groupby,time=NA,order=2,window=5){
  deltaT=1
  interval=1

  res_df = data.frame()
  for (i_group in unique(data[[groupby]])){
    mid_df = data.frame()
    mid_group <- data[data[[groupby]] == i_group,]
    for (i_column in column){
      mid_x <- mid_group[,i_column]
      tEmbedded <- gllaEmbed(mid_x,embed=window,tau=interval)
      wMatrix <- gllaWMatrix(embed=window, tau=interval, deltaT=deltaT, order=order)
      if (nrow(mid_df) == 0){
        cbind_group <- mid_group[,groupby]
        cbind_tE_wM <- tEmbedded[,2:ncol(tEmbedded)] %*% wMatrix
        cbind_group = cbind_group[1:nrow(cbind_tE_wM)]
        glladata <- data.frame(cbind(cbind_group,cbind_tE_wM))
        ggla_name <- c('ID',paste0(i_column,"_", 0:order))
        names(glladata) <- ggla_name

        if(is.na(time)){
          mid_df <- data.frame(row.names = 1:nrow(glladata))
          mid_start_time = (window+1)/2
          mid_end_time = (length(mid_x)-(window-1)/2)
          mid_df$time <- seq(mid_start_time,mid_end_time)
        }else{
          mid_df <- data.frame(row.names = 1:nrow(glladata))
          mid_start_time = (window+1)/2
          mid_end_time = (length(mid_x)-(window-1)/2)
          mid_df$time <- mid_group[[time]][mid_start_time:mid_end_time]
        }

        mid_df <- cbind(mid_df,glladata)
      }else{
        glladata <- data.frame(tEmbedded[,2:ncol(tEmbedded)]%*%wMatrix)
        ggla_name <- c(paste0(i_column,"_", 0:order))
        names(glladata) <- ggla_name
        mid_df <- cbind(mid_df,glladata)
      }

    }
    res_df <- rbind(res_df,mid_df)
  }
  #######################
  # merge data
  if(is.na(time)){
    mid_data <- data.frame()
    for(i_group in unique(data[[groupby]])){
      mid_i_group_df = data[data[[groupby]] == i_group,]
      mid_df <- data.frame(ID=i_group,
                           time=seq(1,nrow(mid_i_group_df)))
      mid_df[column] <- mid_i_group_df[column]
      mid_data <- rbind(mid_data,mid_df)
    }
    res_df <- merge(mid_data,res_df,by.x = c("ID", "time"), by.y = c("ID", "time"), all.x = TRUE, all.y = FALSE)
  }else{
    mid_data <- data[c(groupby,time,column)]
    res_df <- merge(mid_data,res_df,by.x = c(groupby, time), by.y = c("ID", "time"), all.x = TRUE, all.y = FALSE)
  }
  return(res_df)
}

gllaWMatrix <- function(embed=5, tau=1, deltaT=1, order=2) {
  L <- rep(1,embed)
  for(i in 1:order) {
    L <- cbind(L,(((c(1:embed)-mean(1:embed))*tau*deltaT)^i)/factorial(i))
  }
  return(L%*%solve(t(L)%*%L))
}


gllaEmbed <- function(x, embed=2, tau=1, groupby=NA, label="_", idColumn=TRUE) {

  minLen <- (tau + 1 + ((embed - 2) * tau))
  if (!is.vector(groupby) | length(groupby[!is.na(groupby[])])<1) {
    groupby <- rep(1,length(x))
  }
  x <- x[!is.na(groupby[])]
  groupby <- groupby[!is.na(groupby[])]
  if (embed < 2 | is.na(embed) | tau < 1 | is.na(tau) |
      !is.vector(x) | length(x) < minLen)
    return(NA)
  if (length(groupby) != length(x))
    return(NA)
  embeddedMatrix <- matrix(NA, length(x) + (embed*tau), embed+1)
  colNames <- c("ID", paste(label, "0", sep=""))
  for (j in 2:embed) {
    colNames <- c(colNames, paste(label, (j-1)*tau, sep=""))
  }
  dimnames(embeddedMatrix) <- list(NULL, colNames)
  tRow <- 1
  for (i in unique(groupby)) {
    tx <- x[groupby==i]
    if (length(tx) < minLen)
      next
    tLen <- length(tx) - minLen
    embeddedMatrix[tRow:(tRow+tLen), 1] <- i
    for (j in 1:embed) {
      k <- 1 + ((j-1)*tau)
      embeddedMatrix[tRow:(tRow+tLen), j+1] <- tx[k:(k+tLen)]
    }
    tRow <- tRow + tLen + 1
  }
  if (idColumn==TRUE) {
    return(embeddedMatrix[1:(tRow-1),])
  }
  return(embeddedMatrix[1:(tRow-1), 2:(embed+1)])
}
