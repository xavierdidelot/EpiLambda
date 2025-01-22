#Test beta coalescent
context("Test beta coalescent")

test_that("Beta coalescent can be simulated.", {
  expect_silent(tree<-beta_simtree(100,1.5))
  expect_is(tree,'phylo')
  expect_silent(tree<-beta_simtree(100,0))
  expect_equal(tree$Nnode,1)
  expect_silent(tree<-beta_simtree(100,2))
  expect_equal(tree$Nnode,99)
})

test_that("Beta coalescent has same likelihood in simulation and computation.", {
  tree<-beta_simtree(100,0.5)
  expect_equal(tree$ll,beta_loglik(tree,0.5),tolerance=1e-10)
  tree<-beta_simtree(100,1)
  expect_equal(tree$ll,beta_loglik(tree,1),tolerance=1e-10)
  tree<-beta_simtree(100,1.5)
  expect_equal(tree$ll,beta_loglik(tree,1.5),tolerance=1e-10)
})

test_that("Beta coalescent rates are consistent.", {
  n=10;k=3;alpha=0.5
  expect_equal(beta_lambda(n,k,alpha),beta_lambda(n+1,k,alpha)+beta_lambda(n+1,k+1,alpha))
  n=5;k=2;alpha=1
  expect_equal(beta_lambda(n,k,alpha),beta_lambda(n+1,k,alpha)+beta_lambda(n+1,k+1,alpha))
  n=15;k=5;alpha=1.5
  expect_equal(beta_lambda(n,k,alpha),beta_lambda(n+1,k,alpha)+beta_lambda(n+1,k+1,alpha))
})
