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

test_that("Exclusive probabilities add up to one.", {
  n=10;nt=20;k=1:n;expect_equal(sum(pois_exclusive(n=n,k=k,nt=nt)*choose(n-1,k-1)),1)
  n=10;nt=20;k=1:n;expect_equal(sum(negbin_exclusive(n=n,k=k,nt=nt,r=0.1)*choose(n-1,k-1)),1)
})

test_that("Exclusive probabilities are consistent.", {
  n=10;k=3;nt=15
  expect_equal(pois_exclusive(n,k,nt),pois_exclusive(n+1,k,nt)+pois_exclusive(n+1,k+1,nt))
  expect_equal(negbin_exclusive(n,k,nt,r=0.5),negbin_exclusive(n+1,k,nt,r=0.5)+negbin_exclusive(n+1,k+1,nt,r=0.5))
})
