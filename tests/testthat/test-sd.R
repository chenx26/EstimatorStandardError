library(EstimatorStandardError)

context("Test Computing Standard Deviation")

# test_that("The proper name is greeted",{
#           expect_equal(hello.name("abc"),"hello, abc")
#           expect_equal(hello2("abc"),"Greetings, abc")})

test_that("The standard deviation of constant time series is zero",
          expect_equal(sd(rep(1,5)),0))

