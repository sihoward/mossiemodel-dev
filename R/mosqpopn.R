
# Description: ------------------------------------------------------------
#
# function to simulate mosquito populations as a function of temperature
#  - see '' script for parameters matching Niebuhr(2016) PhD thesis
#  - note2
#  - note3
# Built under R version 4.0.3 (2020-10-10)
# Simon Howard; howards@landcareresearch.co.nz | si.w.howard@gmail.com
#-------------------------------------------------------------------------#

mosqpopn <- function(temp_ts, # temperature time series
                     b,       # number of female eggs per clutch
                     alpha,   # adult mortality rate (1/days)
                     beta,    # larval mortality rate (1/days)
                     K_L,     # larval carrying capacity (numbers/km^2)
                     M_max,   # max adult density
                     MTD,     # minimum temperature for mosquito development)
                     L_1 = 0,     # instar densities (1-5)
                     L_2 = 0,
                     L_3 = 0,
                     L_4 = 0,
                     L_5 = 0,
                     M = 0,
                     Mfloor = 100){      # adult density

  # mosquito population model as ordinary differential eqns -----------------
  #
  # (see Niebuhr, 2016. Avian malaria transmission dynamics in New Zealand:
  # investigating host and vector relationships along an elevational gradient. PhD
  # Thesis, U of Otago)

  f <- function(t, y, parms, printProgress = F){
    # debugonce(delta_T_fun)
    # t <- 1; parms <- c(b= b)

    with(as.list(c(y, parms)), {

      # attach(as.list(c(yini, pars))); t <- 0

      temp = temp_ts[floor(t+1)]

      # length gonotrophic cycle
      g_T = if(temp < 1) 241 else 241 * temp^-1.11
      # proportion of adults ovipositing (1/days)
      delta_T = (alpha*exp(-alpha*g_T)) / (1 - exp(-alpha*g_T))
      # larval maturation rate (1/days)
      d_T = if(temp < MTD) 0 else (temp - MTD)/179

      # M <- c(M = 0)
      # L_tot = sum(L)
      dL <- rep(0,5)
      dL[1] = b * delta_T * M * (1-L/K_L) - (1 * d_T + beta) * L_1
      # for(i in 2:5){
      #   # i <- 2
      #   dL[i] = i * d_T * L_i[i-1] - (i*d_T + beta) * L_i[i]
      # }; rm(i)

      dL[2] = 1 * d_T * L_1 - (2 * d_T + beta) * L_2
      dL[3] = 2 * d_T * L_2 - (3 * d_T + beta) * L_3
      dL[4] = 3 * d_T * L_3 - (4 * d_T + beta) * L_4
      dL[5] = 4 * d_T * L_4 - (5 * d_T + beta) * L_5

      L <- dL[1] + dL[2] + dL[3] + dL[4] + dL[5]


      if(printProgress){
        print(sprintf("t = %0.2f; temp = %0.6f; d_T = %0.3f; g_T = %0.3f, L5 = %0.3f; M = %0.3f", t, temp, d_T, g_T, L_5, M))
      }

      dM = if((M + (5 * d_T * L_5 - alpha * M)) < Mfloor) Mfloor-M else 5 * d_T * L_5 - alpha * M

      return(list(c(dL, dM, L)))
    })
  }

  # combine parameters
  pars <- c(b=b, alpha=alpha, beta=beta,
            K_L=K_L, M_Max = M_max, MTD = MTD)
  # combine starting values
  yini <- c(L_1 = L_1, L_2 = L_2, L_3 = L_3,
            L_4 = L_4, L_5 = L_5, M = M)
  yini <- c(yini, L = sum(yini[c("L_1","L_2","L_3","L_4","L_5")]))

  # solve ordinary differential equations over temperature sequence
  out <- deSolve::ode(y = yini, func = f,
             t = c(0, seq_along(temp_ts)),
             parms = pars, method = "euler")

  return(out)
}

#' Run mosquito population model.
#'
#' @param burnin.dates Trim temperature time series to these dates for burn-in.
#' @param burnin.reps Repeat burn-in temperature time series n times.
#' @param run.dates Run temperature time series at these dates after burn-in.
#' @param temp_seq temperature time series
#' @param b number of female eggs per clutch
#' @param alpha adult mortality rate (1/days)
#' @param beta larval mortality rate (1/days)
#' @param K_L larval carrying capacity (numbers/km^2)
#' @param M_max max adult density
#' @param MTD minimum temperature for mosquito development)
#' @param L_1 Initial instar densities (1-5)
#' @param L_2 see above
#' @param L_3 see above
#' @param L_4 see above
#' @param L_5 see above
#' @param M Initial adult density
#' @param Mfloor Minimum adult density
#'
#' @export
#'
#' @examples
#' library(mosqmod)
#'
#' # load stored temperatures
#' temperatures <- data.frame(mosqmod::saved_station_temps)
#' temperatures <- subset(temperatures, Station == "Dunedin, Musselburgh Ews")
#'
#' # append temperatures (requires 'cliflo_usrid' & 'cliflo_pwd' environment
#' # variables to be set if not using request = FALSE)
#' temp_seq <- mosqmod::append_TempSeq(temp_stored = temperatures, request = FALSE)
#'
#' # format sequence (fill gaps, format dates etc. )
#' temp_seq <- formatTempSeq(d = temp_seq)
#'
#' # get projected temperatures from calendar day mean temperatures
#'
#' calendar_day_mean_temps <- getCalendarDayMeans(temp_seq)
#'
#' temp_projected <- project_TempSeq(temp_seq = temp_seq, extend_days = 90,
#'                                   lookback_days = 90,
#'                                   calendar_day_mean_temps)
#'
#' # add projected temperatures
#' temp_seq <- dplyr::bind_rows(temp_seq, temp_projected)
#'
#' # run model
#' res <- runModel(temp_seq = temp_seq,
#'                 burnin.dates = seq(as.Date("2016-07-01"), as.Date("2017-06-30"), 1),
#'                 run.dates = seq(as.Date("2020-07-01"), max(temp_seq$Date), 1))
#' \dontrun{
#' mosqmod::plot_popn(resdf = res, include_temp = TRUE, selectPopn = "M")
#' }
runModel <- function(# burn-in range
  burnin.dates = seq(as.Date("2019-07-01"), as.Date("2020-06-30"), 1),
  burnin.reps = 100,
  # run model between dates (after burn-in)
  run.dates = seq(as.Date("2020-07-01"), max(temp_seq$Date), 1),
  # temperature sequence (Date, Tmean)
  temp_seq = temp_seq[c("Date","Tmean")],
  b = 100,           # number of female eggs per clutch
  alpha = 0.073,     # adult mortality rate (1/days)
  beta = 0.0315,     # larval mortality rate (1/days)
  K_L = 2710200,     # larval carrying capacity (numbers/km^2)
  M_max = 1000060,   # max adult density
  MTD = 7.783,
  L_1 = 0, L_2 = 0,  L_3 = 0, L_4 = 0, L_5 = 0,
  M = 100, Mfloor = 100){

  # repeat burn-in temperatures and append with run dates
  temp_ts <- c(rep(temp_seq$Tmean[temp_seq$Date %in% burnin.dates], burnin.reps),
               temp_seq$Tmean[temp_seq$Date %in% run.dates])

  # run popn model
  modOut <- mosqpopn(temp_ts = temp_ts,
                              b = b,
                              alpha = alpha,
                              beta = beta,
                              K_L = K_L,
                              M_max = M_max,
                              MTD = MTD,
                              L_1 = L_1, L_2 = L_2,  L_3 = L_3, L_4 = L_4, L_5 = L_5,
                              M = M,
                              Mfloor = Mfloor)
  modOut <- data.frame(modOut)
  # add temperature time series
  modOut$Tmean <- c(NA, temp_ts)
  # keep only run dates
  modOut <- subset(modOut, modOut$time > (length(burnin.dates) * burnin.reps))
  # add run dates
  modOut$Date <- run.dates
  # add extra columns from temp_seq
  modOut <- dplyr::left_join(modOut, temp_seq)

  # add attributes to modOut
  attr(modOut, "burnin.range") <- range(burnin.dates)
  attr(modOut, "burnin.reps") <- burnin.reps
  attr(modOut, "run.daterange") <- range(run.dates)
  attr(modOut, "temp_seq") <- temp_seq
  attr(modOut, "params") <- c(b = b, alpha = alpha, beta = beta, K_L = K_L, M_max = M_max, MTD = MTD,
                              L_1 = L_1, L_2 = L_2, L_3 = L_3, L_4 = L_4, L_5 = L_5, M = M, Mfloor = Mfloor)

  return(modOut)
}



#' Temperature requirements for plasmodium development
#'
#' Applies thermal requirements for plasmodium development to temperature time
#' series.
#'
#' Whether the thermal requirements for plasmodium developent are met is
#' controlled by three values. Minimum threshold temperature (\code{MTT}) for
#' development is the temperature above which sporogenic development occurs.
#' Time to development is measured in degree days and is the cumulative sum of
#' the difference between daily temperatures and the MTT. If temperatures remain
#' below the MTT longer than the value of \code{timeout} then the cumulative
#' degree days reset back to zero. Any time the temperature drops below zero
#' (degress Celsius) the cumulative degree days are also reset to zero.
#'
#' Development degree days (\code{devel_degdays}) is the degree days required to
#' complete development, after which the thermal requirements will be met.
#'
#' Defaults for minimum threshold temperature and development degree days are
#' from LaPointe et al. (2010) and the logoc for resetting degree days follows
#' Fortini et al. (2020).
#'
#' Lapointe DA, Goff ML, Atkinson CT 2010. Thermal constraints to the sporogonic
#' development and altitudinal distribution of avian malaria plasmodium relictum
#' in Hawai’i. Journal of Parasitology. 96(2):318–324.
#'
#' Berio Fortini L, Kaiser LR, Lapointe DA, Berio L, Kaiser LR, Lapointe DA
#' 2020. Fostering real-time climate adaptation: Analyzing past, current, and
#' forecast temperature to understand the dynamic risk to Hawaiian honeycreepers
#' from avian malaria. Global Ecology and Conservation. 23:e01069.
#'
#' @param envtemp temperature time series
#' @param MTT minimum threshold temperature (MTT) for sporogenic development in plasmodium
#' @param devel_degdays degree days above MTT to complete development
#' @param timeout reset cumulative degree days if temperature below MTT for consecutive values
#'
#' @return \code{plasmod_devel()} returns a \code{\link[base]{data.frame}} with
#'   the initial temperature time series ('envtemp'), the physiologically
#'   effective temperature (temperature minus the MTT: 'phystemp'), days since
#'   last zero physiologically effective temperature ('lastzero'), degree days
#'   ('degdays') and whether the thermal requirements for development are met
#'   ('thermal_req').
#' @export
#'
#' @examples
#' d <- plasmod_devel(envtemp = (saved_station_temps$`Tmin(C)`[1:(3*365)] +
#'                      saved_station_temps$`Tmax(C)`[1:(3*365)])/2)
#' \dontrun{
#' plot(y = d$envtemp, x = seq_along(d$envtemp), type = "l", ylim = c(0, max(d$envtemp)))
#' points(y = d$degdays / max(d$degdays) * max(d$envtemp), x = seq_along(d$envtemp),
#'        type = "l", col = "red")
#' abline(h = 12.97, col = "blue")
#' }
plasmod_devel <- function(envtemp, MTT = 12.97, devel_degdays = 86.21, timeout = 30L){

  if(any(is.na(envtemp))) stop("NAs present in envtemp argument passed to plasmod_devel function")

  # get degree days above minimum threshold temperature
  out <- getDegreeDays(envtemp, developtemp = MTT, timeout)

  # degrees days meet devel_degdays
  thermal_req <- out$degdays >= devel_degdays

  d <- data.frame(out, thermal_req)

  return(d)
}


#' Run shiny app locally
#'
#' Runs shiny app locally using the 'app/app.R' file included in 'mosqmod'
#' package (i.e. 'Program Files/R/R.X.X.X/library/mosqmod').
#'
#' @param cliflo_requests Make requests to CliFlo database to update
#'   temperatures (default = TRUE). Set to FALSE to disable when working on
#'   restricted networks for example.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' mosqmod::runMosqModApp(cliflo_requests = FALSE)
#' }
runMosqModApp <- function(cliflo_requests = TRUE){
  cliflo_requests <<- cliflo_requests   # add to global environment
  shiny::runApp(appDir = system.file("app", package = "mosqmod"),
                launch.browser = TRUE)
}






