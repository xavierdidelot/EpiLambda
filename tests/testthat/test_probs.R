#Test probabilities
context("Test probabilities")

test_that("p2t and pnt are consistent.", {
  #Poisson
  expect_equal(pois_p2t(10),pois_pnt(2,10))

  #Negative-Binomial
  expect_equal(negbin_p2t(10,0.5),negbin_pnt(2,10,0.5))
})

test_that("Logged and unlogged probabilities are consistent.", {
  #Poisson
  expect_equal(log(pois_p2t(10)),pois_p2t(10,log=T))
  expect_equal(log(pois_pnt(5,10)),pois_pnt(5,10,log=T))

  #Negative-Binomial
  expect_equal(log(negbin_p2t(10,0.5)),negbin_p2t(10,0.5,log=T))
  expect_equal(log(negbin_pnt(5,10,0.5)),negbin_pnt(5,10,0.5,log=T))
})
