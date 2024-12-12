#Test new lambda-coalescent
context("Test new lambda-coalescent")

test_that("New lambda-coalescent can be simulated.", {
  expect_silent(tree<-new_simtree(20,40,0.2))
  expect_is(tree,'phylo')
})

test_that("New lamba-coalescent has same likelihood in simulation and computation.", {
  tree<-new_simtree(20,40,0.2)
  expect_equal(tree$ll,new_loglik(tree,40,0.2),tolerance=1e-10)
  tree<-new_simtree(5,10,0.5)
  expect_equal(tree$ll,new_loglik(tree,10,0.5),tolerance=1e-10)
  tree<-new_simtree(20,30,1.1)
  expect_equal(tree$ll,new_loglik(tree,30,1.1),tolerance=1e-10)
})
