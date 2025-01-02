#Test simulation
context("Test MLE")

test_that("MLE code works.", {
  set.seed(0)
  expect_silent(t<-beta_simtree(n=50,alpha=1))
  expect_silent(beta_mle(t))
  expect_silent(t<-omega_simtree(n=50,nt=100,r=0.1))
  expect_silent(omega_mle(t,nt=100))
  expect_silent(omega_mle(t,r=0.1))
  expect_silent(omega_mle(t))
})
