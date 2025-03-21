#Test simulation
context("Test simulation")

test_that("Simulation code works.", {
  expect_silent(simul_pcoal(rt=1.5,nt=10,nrep=10))
  expect_silent(simul_pcoal(rt=1.5,vt=10,nt=10,nrep=10))
  expect_silent(simul_pcoal(rt=1.5,nt=10,nrep=10,method='multinomial'))
  expect_silent(simul_pcoal(rt=1.5,vt=10,nt=10,nrep=10,method='multinomial'))
  expect_silent(simul_pcoal(rt=1.5,nt=10,n=8,type='exclusive',nrep=10))
  expect_silent(simul_pcoal(rt=1.5,vt=10,n=8,type='exclusive',nt=10,nrep=10))
  expect_silent(simul_pcoal(rt=1.5,nt=10,n=8,type='exclusive',nrep=10,method='multinomial'))
  expect_silent(simul_pcoal(rt=1.5,vt=10,n=8,type='exclusive',nt=10,nrep=10,method='multinomial'))
})
