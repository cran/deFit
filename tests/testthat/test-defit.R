test_that("Univariable second-order differential test", {
  data('example1')
  model1 <- '
      X =~ myX
      time =~ myTime
      X(2) ~ X(1) + X
      '
  result1 <- defit(data = example1, model = model1,method = 'Nelder-Mead')
  if(is.null(result1)){
    test_out = 'error'
  }else{
    test_out = 'pass'
  }
  expect_equal(test_out,'pass')
})
test_that("bivariate first-order differential test", {
  data('example2')
  model2 <- '
    # define variable
    X =~ myX
    Y =~ myY

    # define time
    time =~ myTime

    # define differential equation
    X(1) ~ X + Y
    Y(1) ~ Y + X
  '
  result2 <- defit(data = example2, model = model2)
  if(is.null(result2)){
    test_out2 = 'error'
  }else{
    test_out2 = 'pass'
  }
  expect_equal(test_out2,'pass')
})
test_that("Guess values test", {
  data('example2')
  model2 <- '
    # define variable
    X =~ myX
    Y =~ myY

    # define time
    time =~ myTime

    # define differential equation
    X(1) ~ X + Y
    Y(1) ~ Y + X
  '
  result3 <- defit(data = example2, model = model2,guess=c(0,0,0,0,0,0))
  if(is.null(result3)){
    test_out3 = 'error'
  }else{
    test_out3 = 'pass'
  }
  expect_equal(test_out3,'pass')
})
test_that("Multilevel univariate differential euqation", {
  data('example3')
  model3 <- '
   X =~ current
   time =~ myTime
   X(2) ~ X(1) + X + (1 + X(1) + X | year)
   '
  example3_use <- example3[(example3["year"] >= 2015)&(example3["year"] <= 2018),] # Note: select a subset of the data as an example.
  example3_c <- Scale_within(example3_use, model3) # note: centering X variable by year
  result4 <- defit(data=example3_c,model = model3,plot=FALSE)
  if(is.null(result4)){
    test_out4 = 'error'
  }else{
    test_out4 = 'pass'
  }
  expect_equal(test_out4,'pass')
})
test_that("Plot: Multilevel univariate differential euqation", {
  data('example3')
  model3 <- '
   X =~ current
   time =~ myTime
   X(2) ~ X(1) + X + (1 + X(1) + X | year)
   '
  example3_use <- example3[(example3["year"] >= 2015)&(example3["year"] <= 2018),] # Note: select a subset of the data as an example.
  example3_c <- Scale_within(example3_use, model3) # note: centering X variable by year
  result5 <- defit(data=example3_c,model = model3,plot=FALSE)
  if(is.null(result5)){
    test_out5 = 'error'
  }else{
    test_out5 = 'pass'
  }
  expect_equal(test_out5,'pass')
}
)
test_that("guess and plot: Multilevel univariate differential euqation", {
  data('example3')
  model3 <- '
   X =~ current
   time =~ myTime
   X(2) ~ X(1) + X + (1 + X(1) + X | year)
   '
  example3_use <- example3[(example3["year"] >= 2015)&(example3["year"] <= 2018),] # Note: select a subset of the data as an example.
  example3_c <- Scale_within(example3_use, model3) # note: centering X variable by year
  result6 <- defit(data=example3_c,model = model3,guess=list(c(0,0,0,0),c(0,0,0,0)),plot=TRUE)
  if(is.null(result6)){
    test_out6 = 'error'
  }else{
    test_out6 = 'pass'
  }
  expect_equal(test_out6,'pass')
}
)
test_that("guess, plot and method: Multilevel univariate differential euqation", {
  data('example3')
  model3 <- '
   X =~ current
   time =~ myTime
   X(2) ~ X(1) + X + (1 + X(1) + X | year)
   '
  example3_use <- example3[(example3["year"] >= 2015)&(example3["year"] <= 2018),] # Note: select a subset of the data as an example.
  example3_c <- Scale_within(example3_use, model3) # note: centering X variable by year
  result7 <- defit(data=example3_c,model = model3,guess=list(c(0,0,0,0),c(0,0,0,0),mthod=list('L-BFGS-B','BFGS')),plot=TRUE)
  if(is.null(result7)){
    test_out7 = 'error'
  }else{
    test_out7 = 'pass'
  }
  expect_equal(test_out7,'pass')
}
)
# test_that("example1 simulation", {
  # damped_oscillator <- function(t, y, parms) {
  #   with(as.list(parms), {
  #     dy1 <- y[2]
  #     dy2 <- beta1 * y[1] + beta2 * y[2]
  #     return(list(c(dy1, dy2)))
  #   })
  # }
  #
  # mid_solve_data <- deSolve::ode(y = c(1,0.1),
  #                                times = seq(1,30),
  #                                func = damped_oscillator,
  #                                parms = c(beta1=-0.3,beta2=-0.1))
  #
  # mid_solve_data
  # res = mid_solve_data
  # res[,2] = mid_solve_data[,2] + rnorm(30,0,sd(mid_solve_data[,2])/2)
  # res[,3] = mid_solve_data[,3] + rnorm(30,0,sd(mid_solve_data[,3])/2)
  # plot(res)
  # res[,'myX'] = res[,2]
  #
  # out = data.frame()
  # out$myX = res[,2]
  # out
  #
  # library(readr)
  # mydata <- read_csv("data/example1.csv")
  # mydata$myTime <- mydata$myTime
  # str(mydata)
# }
# )
