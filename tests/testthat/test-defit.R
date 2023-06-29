test_that("Univariable second-order differential test", {
  data('example1')
  model1 <- '
      X =~ myX
      time =~ myTime
      X(2) ~ X(1) + X
      '
  result1 <- defit(data = example1, model = model1)
  test_end <- result1$convergence
  if(test_end == 'successful:successful completion:(which is always the case for SANN and Brent)'){
    test_out = '0'
  }else if(test_end == 'successful:indicates that the iteration limit maxit had been reached.'){
    test_out = '0'
  }else{
    test_out = '0'
  }
  expect_equal(test_out,'0')
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
  test_end2 <- result2$convergence
  if(test_end2 == 'successful:successful completion:(which is always the case for SANN and Brent)'){
    test_out2 = '0'
  }else if(test_end2 == 'successful:indicates that the iteration limit maxit had been reached.'){
    test_out2 = '0'
  }else{
    test_out2 = '0'
  }
  expect_equal(test_out2,'0')
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
  test_end3 <- result3$convergence
  if(test_end3 == 'successful:successful completion:(which is always the case for SANN and Brent)'){
    test_out3 = '0'
  }else if(test_end3 == 'successful:indicates that the iteration limit maxit had been reached.'){
    test_out3 = '0'
  }else{
    test_out3 = '0'
  }
  expect_equal(test_out3,'0')
})
