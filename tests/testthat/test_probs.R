#Test probabilities
context("Test probabilities")

test_that("Logged and unlogged probabilities are consistent.", {
  #Poisson
  expect_equal(log(pois_inclusive(5,10)),pois_inclusive(5,10,log=T))
  expect_equal(log(pois_exclusive(5,8,10)),pois_exclusive(5,8,10,log=T))

  #Negative-Binomial
  expect_equal(log(negbin_inclusive(5,10,0.5)),negbin_inclusive(5,10,0.5,log=T))
})

test_that("Inclusive and exclusive probabilities are equal if n=k.", {
  #Poisson
  expect_equal(pois_inclusive(5,10),pois_exclusive(5,5,10))
})
