#Test omega-coalescent
context("Test new omega-coalescent")

test_that("Omega-coalescent can be simulated.", {
  expect_silent(tree<-omega_simtree(20,40,0.2))
  expect_is(tree,'phylo')
})

test_that("Omega-coalescent has same likelihood in simulation and computation.", {
  tree<-omega_simtree(20,40,0.2)
  expect_equal(tree$ll,omega_loglik(tree,40,0.2),tolerance=1e-10)
  tree<-omega_simtree(5,10,0.5)
  expect_equal(tree$ll,omega_loglik(tree,10,0.5),tolerance=1e-10)
  tree<-omega_simtree(20,30,1.1)
  expect_equal(tree$ll,omega_loglik(tree,30,1.1),tolerance=1e-10)
})
