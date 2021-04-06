



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
#' @param L_1 instar densities (1-5)
#' @param L_2 see above
#' @param L_3 see above
#' @param L_4 see above
#' @param L_5 see above
#' @param M adult density
#'
#' @return 1
#' @export 1
#'
#' @examples 1

mosqpopn <- function(temp_ts, # temperature time series
                     b,       # number of female eggs per clutch
                     alpha,   # adult mortality rate (1/days)
                     beta,    # larval mortality rate (1/days)
                     K_L,     # larval carrying capacity (numbers/km^2)
                     M_max,   # max adult density
                     MTD,     # minimum temperature for mosquito development)
                     L_1,     # instar densities (1-5)
                     L_2,
                     L_3,
                     L_4,
                     L_5,
                     M){      # adult density

  require(deSolve)

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
  out <- ode(y = yini, func = f,
             t = c(0, seq_along(temp_ts)),
             parms = pars, method = "euler")

  return(out)
}
