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
