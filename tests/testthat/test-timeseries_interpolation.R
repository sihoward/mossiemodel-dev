context("Check interpolation of temperature time series")


test_that("time series formatting interpolates missing points", {

  # data frame to pass to formatTempSeq()
  d.missing <- d <- data.frame(`Date.local.` = sprintf("200001%02.f:0900", 1:25), `Tmin.C.` = 1:25, `Tmax.C.` = 25:1)

  d.missing$Tmin.C.[seq(2,24,2)] <- NA    # create gaps in time series
  d.missing$Tmax.C.[seq(2,24,2)] <- NA

  res <- mosqmod::formatTempSeq(d.missing)

  # tests
  expect_equal(res$Tmin, d$Tmin.C., tolerance = 1e-6)
  expect_equal(res$Tmax, d$Tmax.C., tolerance = 1e-6)

})












