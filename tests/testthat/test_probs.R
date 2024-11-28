#Test probabilities
context("Test probabilities")

test_that("Logged and unlogged probabilities are consistent.", {
  #Poisson
  expect_equal(log(pois_inclusive(k=5,nt=10)),pois_inclusive(k=5,nt=10,log=T))
  expect_equal(log(pois_exclusive(k=5,n=8,nt=10)),pois_exclusive(k=5,n=8,nt=10,log=T))

  #Negative-Binomial
  expect_equal(log(negbin_inclusive(k=5,nt=10,r=0.5)),negbin_inclusive(k=5,10,0.5,log=T))
  expect_equal(log(negbin_exclusive(k=5,n=8,nt=10,r=0.5)),negbin_exclusive(k=5,n=8,nt=10,r=0.5,log=T))
})

test_that("Inclusive probability is one for k=1.", {
  expect_equal(pois_inclusive(k=1,nt=10),1)
  expect_equal(negbin_inclusive(k=1,nt=10,r=5),1)
})

test_that("Inclusive and exclusive probabilities are equal if n=k.", {
  expect_equal(pois_inclusive(k=5,nt=10),pois_exclusive(k=5,n=5,nt=10))
  expect_equal(negbin_inclusive(k=5,nt=10,r=5),negbin_exclusive(k=5,n=5,nt=10,r=5))
})
