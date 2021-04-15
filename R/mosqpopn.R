
# Description: ------------------------------------------------------------
#
# function to simulate mosquito populations as a function of temperature
#  - see '' script for parameters matching Niebuhr(2016) PhD thesis
#  - note2
#  - note3
# Built under R version 4.0.3 (2020-10-10)
# Simon Howard; howards@landcareresearch.co.nz | si.w.howard@gmail.com
#-------------------------------------------------------------------------#

#' Title
#'
#' @param temp_ts temperature time series
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
#'
#' @return
#' @export
#'
#' @examples

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
                     M = 0){      # adult density

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

      dM = if((M + (5 * d_T * L_5 - alpha * M)) < 0.01) -M + 0.01 else 5 * d_T * L_5 - alpha * M

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
#'
#' @return
#' @export
#'
#' @examples
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
  M = 100){

  # repeat burn-in temperatures and append with run dates
  temp_ts <- c(rep(temp_seq$Tmean[temp_seq$Date %in% burnin.dates], burnin.reps),
               temp_seq$Tmean[temp_seq$Date %in% run.dates])

  # run popn model
  modOut <- mosqmod::mosqpopn(temp_ts = temp_ts,
                              b = b,
                              alpha = alpha,
                              beta = beta,
                              K_L = K_L,
                              M_max = M_max,
                              MTD = MTD,
                              L_1 = L_1, L_2 = L_2,  L_3 = L_3, L_4 = L_4, L_5 = L_5,
                              M = M)
  modOut <- data.frame(modOut)
  # add temperature time series
  modOut$Tmean <- c(NA, temp_ts)
  # keep only run dates
  modOut <- subset(modOut, modOut$time > (length(burnin.dates) * burnin.reps))
  # add run dates
  modOut$Date <- run.dates
  # add extra columns from temp_seq
  modOut <- dplyr::left_join(modOut, temp_seq)

  return(modOut)
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
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' mosqmod::runMosqModApp()
#' }

runMosqModApp <- function(cliflo_requests = TRUE){
  shiny::runApp(appDir = system.file("app", package = "mosqmod"),
                launch.browser = TRUE)
}






