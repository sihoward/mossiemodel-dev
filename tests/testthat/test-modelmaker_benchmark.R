context("match ModelMaker outputs to mospopn")

# Description: ------------------------------------------------------------
#
# Re-create mosquito model to match saved inputs & outputs from ModelMaker 
#  software
#  - inputs and outputs were saved to csv from original ModelMaker runs
#  - sources mosqpopn() function
#
# Built under R version 4.0.3 (2020-10-10)
# Simon Howard; howards@landcareresearch.co.nz | si.w.howard@gmail.com
#-------------------------------------------------------------------------#

test_that("mospopn matches ModelMaker output", {
  
  # load benchmark temperature time-series (Penguin 89-90) ------------------
  d <- read.csv(file = system.file("benchmark/benchmark_temperatures.csv", package = "mosqmod"))
  temp_seq <- rep(d$temp, 100)     # repeat 100 times for burn-in
  
  # run model using ModelMaker (MM) parameters ------------------------------
  
  out <- mosqmod::mosqpopn(temp_seq,
                           b = 100,           # number of female eggs per clutch
                           alpha = 0.073,     # adult mortality rate (1/days)
                           beta = 0.0315,     # larval mortality rate (1/days)
                           K_L = 2710200,     # larval carrying capacity (numbers/km^2)
                           M_max = 1000060,   # max adult density
                           MTD = 7.783,       # minimum temperature for mosquito development)
                           L_1 = 0, L_2 = 0,  L_3 = 0, L_4 = 0, L_5 = 0,
                           M = 100)
  
  # compare to MM results ---------------------------------------------------
  
  # read saved MM results
  popn <- read.csv(system.file("benchmark/benchmark_population.csv", package = "mosqmod"))
  
  # tests
  
  # NOTE: exported ModelMaker results get rounded for large estimates so
  # tolerance value set lower when predicted values are large (i.e. for adult
  # mosquito - M estimates))
  expect_equal(popn$M, out[,"M"], tolerance = 1e-3)
  expect_equal(popn$L5, out[,"L_5"], tolerance = 1e-6)
  
})












