#Test probabilities
context("Test probabilities")

test_that("p2t and pnt are consistent.", {
  #Poisson
  expect_equal(pois_p2t(10),pois_pnt(2,10))

  #Negative-Binomial
  expect_equal(negbin_p2t(10,0.5),negbin_pnt(2,10,0.5))
})

