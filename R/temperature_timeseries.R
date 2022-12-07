#' Format cliflo temperature time series
#'
#' Formats raw cliflo temperature time series for use in later functions. A
#' major operation is filling gaps in the timeseries which are present in the
#' cliflo database. Gaps are filled using a straight line interpolation between
#' the points adjoining the gap in the record.
#'
#' Interpolated values are also stored in a separate column 'interpolated' so
#' gaps can be identified (dates without gaps have NA in this column).
#'
#' Mean temperatures are calculated by dividing the daily min and max
#' temperatures by 2.
#'
#' @param d data.frame from \code{\link[mosqmod]{append_TempSeq}}.
#'
#' @export
formatTempSeq <- function(d){

  d <- dplyr::rename(d, Date = .data$Date.local., Tmax = .data$Tmax.C., Tmin = .data$Tmin.C.)

  # convert date strings
  d$Date_orig <- d$Date
  d$datetime <- as.POSIXlt(strptime(d$Date_orig, format = "%Y%m%d:%H%M", tz = "NZ"))
  d$Date <- as.Date(d$datetime)

  date_seq <- seq(min(d$Date), max(d$Date), by = "day")

  # check leap years
  # date_seq[grepl("02-29", format(date_seq, "%m-%d"))]

  # check same time each day
  # table(format(date_seq, "%H:%M"))

  # match each measurement to contiguous sequence
  d <- dplyr::left_join(data.frame(Date = date_seq), d)

  # check for missing measurements
  if(any(is.na(d$Tmin) | is.na(d$Tmax))){
    message(sprintf("%s of %s rows missing temperatures across %s separate gaps at rows %s",
                    sum(is.na(d$Tmin) | is.na(d$Tmax)),
                    nrow(d),
                    sum(diff(which(is.na(d$Tmin)|is.na(d$Tmin))) > 1) + 1,
                    paste(which(is.na(d$Tmin)|is.na(d$Tmax)), collapse = ",")))
  }

  # interploate points between
  # - creates sequence of all dates to save adding logic to identify gaps
  # - only replaces missing values
  # - keeps a record of interpolated points with 'interpolated' column
  d$interpolated <- is.na(d$Tmin)|is.na(d$Tmax)
  Tmin.est <- stats::approx(y = d$Tmin[!is.na(d$Tmin)], x = which(!is.na(d$Tmin)),
                     xout = seq_len(nrow(d)))$y
  d$Tmin[is.na(d$Tmin)] <- Tmin.est[is.na(d$Tmin)]

  Tmax.est <- stats::approx(y = d$Tmax[!is.na(d$Tmax)], x = which(!is.na(d$Tmax)),
                     xout = seq_len(nrow(d)))$y
  d$Tmax[is.na(d$Tmax)] <- Tmax.est[is.na(d$Tmax)]

  # check for remaining NAs
  if(any(is.na(d$Tmin)|is.na(d$Tmax))) stop("NAs remain in temperature series")

  #-------------------------------------------------------------------------#
  # illustrative plot showing interpolation between points
  # plot(Tmin ~ Date, data = subset(d[5800:6200,], !interpolated), type = "l")
  # lines(y = Tmin.est[d$interpolated], x = d$Date[d$interpolated], col = "green", type = "b", pch = 19)
  #-------------------------------------------------------------------------#

  # get mean temperature
  d$Tmean <- (d$Tmin + d$Tmax)/2

  d$calday <- format(d$Date, "%m-%d")

  return(d)
}

#' Append temperature times series with latest temperature from NIWA's CliFlo database.
#'
#' @param temp_stored Stored temperatures to append.
#' @param date_append Get updated temperatures after this date. Defaults to current system clock.
#' @param request Make CliFlo request (default = TRUE).
#' @param username CliFlo login username (see 'https://cliflo.niwa.co.nz')
#' @param password CliFlo login password.
#'
#' @export
#'
#' @examples
#' # get stored data
#' dat <- data.frame(mosqmod::saved_station_temps)
#' # subset selected temperature series
#' subdat <- subset(dat, Station == "Dunedin, Musselburgh Ews")
#'
#' temp_hist <- mosqmod::append_TempSeq(temp_stored = subdat,
#'                                      request = FALSE,
#'                                      date_append = as.Date(Sys.Date()))
append_TempSeq <- function(temp_stored,
                           date_append = as.Date(Sys.time()),
                           request = TRUE,
                           username = NULL, password = NULL){

  # get latest date from
  date_last_stored <- max(as.Date(temp_stored$Date.local., format = "%Y%m%d"))

  ## CliFro request latest temperatures


  if (request) {

    # check CliFlo credentials
    if(is.null(username) | is.null(username)) stop("Include CliFlo login and password when calling append_TempSeq(); see 'https://cliflo.niwa.co.nz/'")

    # add user & stations
    me <- clifro::cf_user(username = username,
                          password = password)
    my.dts <- clifro::cf_datatype(select_1 = 4, select_2 = 2, check_box = 1)
    # find matching station ID from stored name
    match.StationID <- clifro::cf_find_station(temp_stored$Station[1], search = "name")[["agent"]]
    # my.stations = clifro::cf_station(15752)       # 15752 - Musselburgh station ID
    my.stations <- clifro::cf_station(match.StationID)

    if(!all(temp_stored$Station %in% my.stations$name)){
      stop("Stored temperatures contain station IDs missing from CliFro request")
    }

    # fetch data
    cf.temp_latest = clifro::cf_query(user = me,
                                      datatype = my.dts,
                                      station = my.stations,
                                      start_date = strftime(date_last_stored + 1, format = "%Y-%m-%d 00"),
                                      end_date = strftime(date_append, format = "%Y-%m-%d 00"))
    # append updated data to stored
    temp_append <- rbind(data.frame(temp_stored, source = "stored"),
                         data.frame(cf.temp_latest, source = "retreived"))
  } else {
    temp_append <- data.frame(temp_stored, source = "stored")
  }

  return(temp_append)
}


#' Projected temperatures from current temperature difference between calendar
#' day means
#'
#' @param temp_seq Formatted temperature time series from
#'   \code{\link[mosqmod]{formatTempSeq}}.
#' @param extend_days Number of days to extend time series.
#' @param lookback_days Number of preceeding days to calculate difference
#'   between calendar day means and recent temperatures.
#' @param calday_means Calendar day means used to calculate differences
#'   between historical and recent temperatures. Generated using from
#'   \code{\link[mosqmod]{getCalendarDayMeans}}.
#'
#' @export
project_TempSeq <- function(temp_seq = temp_seq,
                            extend_days = 30,
                            lookback_days = 90,
                            calday_means = calday_means){

  # latest date in recorded temperatures
  date_latest <- max(temp_seq$Date)

  # projected temps
  temp_proj <- data.frame(Date = seq(date_latest + 1, date_latest + extend_days, 1),
                          Tmean = NA, calday = NA, source = "projected")
  temp_proj$calday <- format(temp_proj$Date, "%m-%d")

  # add calendar day means
  temp_seq <- dplyr::left_join(temp_seq, calday_means)
  temp_proj <- dplyr::left_join(temp_proj, calday_means)

  # adjust historical mean +/- difference between historical and current temps
  # - average difference from calendar day mean for last recorded n days
  # - difference to calendar day mean
  temp_lookback <- utils::tail(temp_seq, lookback_days)
  temp_delta <- mean(temp_lookback$Tmean - temp_lookback$Tmean_calday)

  message(sprintf("%s-day difference in mean and mean calendar day temperatures = %0.3f degreesC", lookback_days, temp_delta))

  temp_proj$Tmean <- temp_proj$Tmean_calday + temp_delta

  return(temp_proj)
}


#' Projected degree days temperatures from current difference between degree
#' days and calendar degree days
#'
#' @param temp_seq Formatted temperature time series from
#'   \code{\link[mosqmod]{formatTempSeq}}.
#' @param extend_days Number of days to extend time series.
#' @param developtemp TODO
#' @param timeout TODO
#' @param cal_degday_means Calendar degree day means used to calculate
#'   differences between historical and recent degree days. Generated using
#'   from \code{\link[mosqmod]{getCalendarDegreeDays}}.
#'
#' @export
#'
#' @examples
#' data("saved_station_temps")
#' temp_seq <- formatTempSeq(subset(data.frame(saved_station_temps),
#'                                  Station == "Nugget Point Aws"))
#' # calculate mean degree days for each celendar day on last ten years
#' cal_degday_means <- getCalendarDegreeDays(temp_hist = temp_seq,
#'                                           lookback = "-10 year",
#'                                           developtemp = 12.97, timeout = 30)
#' # pass last year of temperature time series and project degree days
#' project_DegreeDays(temp_seq = tail(temp_seq, 365),
#'                    extend_days = 30,
#'                    developtemp = 12.97, timeout = 30,
#'                    cal_degday_means = cal_degday_means)
project_DegreeDays <- function(temp_seq = temp_seq,
                               extend_days = 30,
                               developtemp,
                               timeout,
                               cal_degday_means){

  # check current temperature sequence covers > 1 year
  if(nrow(temp_seq) < 365){
    stop("project_DegreeDays(): Temperature time series < 1 year.
         1+ years of temperatures (covering a winter) needed to project degree days")
  }

  # latest date in recorded temperatures
  date_latest <- max(temp_seq$Date)

  # projected degree days
  degday_proj <- data.frame(Date = seq(date_latest + 1, date_latest + extend_days, 1),
                            degdays = NA, calday = NA, source = "projected")
  degday_proj$calday <- format(degday_proj$Date, "%m-%d")

  # add calendar day means
  temp_seq <- dplyr::left_join(temp_seq, cal_degday_means)
  degday_proj <- dplyr::left_join(degday_proj, cal_degday_means)

  # adjust historical mean +/- difference between historical and current degree days
  # - average difference from calendar day mean for last recorded n days
  # - difference to calendar day mean

  # get degree days from entire temperature time series up to present
  degdays_past <- getDegreeDays(temp_seq$Tmean, developtemp, timeout)
  temp_seq$degdays <- degdays_past[["degdays"]]

  # check that degree days are reset via timeout sometime in the last year
  if(all(utils::tail(degdays_past$lastzero, 365) == 0L)){
    stop("No timeout triggered in degree day timeseries; temperatures never < development temperature in previous 365 days")
  }

  # get difference between current and annual mean degree days on the same calendar day
  degdays_delta <- temp_seq$degdays[temp_seq$Date == date_latest] - temp_seq$meanDegDays[temp_seq$Date == date_latest]

  # if(abs(degdays_delta) > 10)
  if(as.numeric(format(date_latest, "%m")) %in% 5:7){
    warning("degree days typically reset to zero from May-Aug. ",
            "Adjusting mean calendar day degree-days using current temperatures is unreliable around these months...\n",
            "using unadjusted mean calendar day degree-days")
    degdays_delta <- 0
  }
  if(abs(degdays_delta) > 50){
    warning("difference in mean calendar day degree-days and current degree-days is > +/- 50 degree-days ...\nusing unadjusted mean calendar day degree-days")
    degdays_delta <- 0
  }

  message(sprintf("difference in degree days between current and past degree days for latest date = %0.3f degreesC", degdays_delta))

  # add difference to calendar day degree days
  degday_proj$degdays <- degday_proj$meanDegDays + degdays_delta

  #------------------------------------------------------------------------------#
  # plots for troubleshooting
  #
  # plotddays <- c(tail(temp_seq$degdays, 100), degday_proj$degdays)
  # plot(plotddays, type = "l", ylim = range(c(degday_proj$meanDegDays, plotddays)))
  # lines(c(tail(temp_seq$meanDegDays, 100), degday_proj$meanDegDays), type = "l", col = "red")
  # lines(c(rep(0, 100), degday_proj$degdays), type = "l", col = "blue")

  return(degday_proj)
}



#' Calculate calendar day mean temperatures from historical series
#'
#' @param temp_hist Temperature time series formatted using
#'   \code{\link[mosqmod]{formatTempSeq}}.
#' @param lookback Subset the historical time series starting at a past date.
#'   Must match the format for the 'by' argument in \code{\link{seq.Date}}.
#' @param series_ending Subset the historical time series ending date. Must be a
#'   valid \code{\link{Dates}} object. When set to NULL (default) the most
#'   recent date in the temperature time series is used.
#'
#' @export
#'
#' @examples
#' ## get past temperature series from append_TempSeq example
#' example("append_TempSeq")
#' temp_hist <- mosqmod::formatTempSeq(temp_hist)
#' # using default - 10 years from most recent temperature
#' getCalendarDayMeans(temp_hist = temp_hist)
# \dontrun{
# # customised - 5 years from end-2015
# getCalendarDayMeans(temp_hist = temp_hist, lookback = "-5 year",
#                     series_ending = as.Date("2015-12-31"))
# }
getCalendarDayMeans <-
  function(temp_hist = temp_hist,
           lookback = "-10 year",
           series_ending = NULL){

  Date <- NULL

  # use last date in temp series if no series_ending date given
  if(is.null(series_ending)){
    series_ending <- max(temp_hist$Date)
  }

  # get dates starting the specified range
  start_end_dates <- seq(series_ending, by = lookback, length.out = 2)

  # check supplied historical date range matches period to average over
  if(min(temp_hist$Date) > start_end_dates[2]){
    stop("Historical temperatures start after date range for calendar dat averages")
  }

  # only keep past temperatures in specified date range
  subdat <- subset(temp_hist, Date >= start_end_dates[2])

  # calculate calendar day means
  calday_means <- stats::aggregate(Tmean ~ calday, FUN = mean, data = subdat)

  # rename Tmean column
  calday_means$Tmean_calday <- calday_means$Tmean
  calday_means$Tmean <- NULL

  # replace leap day with mean for 28th Feb
  calday_means$Tmean_calday[calday_means$calday == "02-29"] <-
    calday_means$Tmean_calday[calday_means$calday == "02-28"]

  attr(calday_means, "start_end_dates") <- start_end_dates

  return(calday_means)
}

#' Degree days above developmental temperature from temperature time series
#'
#' Environmental temperatures are input and the developmental temperature is
#' subtracted to get the physiologically effective temperature, and degrees days
#' are the cumulative sum of the physiologically effective temperature.
#'
#' A timeout of consecutive days below physiologically effective temperature can
#' be used and after which the degree days are reset back to zero. This
#' represents that an animal cannot stay below physiological temperature
#' indefinitely.
#' see Trudgill et al. (2005) for details:
#'
#' Trudgill et al. 2005. Thermal time - Concepts and utility.
#' Annals of Applied Biology. 146(1):1–14.
#'
#' @param envtemp TODO
#' @param developtemp TODO
#' @param timeout TODO
#'
#' @export
getDegreeDays <- function(envtemp, developtemp, timeout){

  # physiologically effective temperature (see Trudgill et al. 2005 p2)
  phystemp <- (envtemp - developtemp) * ((envtemp - developtemp) > 0)

  # days since last zero degree day
  lastzero <- cumsum(phystemp == 0) - cummax(cumsum(phystemp == 0) * (phystemp > 0))

  # cumulative degree days (using physiologically effective temperature)
  degdays <- cumsum(phystemp)

  # reset degree days when consecutive days below developmental temperature is
  # greater than the timeout
  # if > timeout days subtract cumulative phystemp from last day < developmental temperature
  # if envtemp < 0 subtract cumulative phystemp from last day above zero
  degdays <- degdays - cummax((lastzero > timeout | envtemp <= 0) * cumsum(phystemp))

  # combine temp series into dataframe
  d <- data.frame(envtemp, phystemp, lastzero, degdays)
  # add developmental temperature as attribute
  attributes(d) <- c(attributes(d), developtemp = developtemp)

  return(d)
}
#-------------------------------------------------------------------------#
# getDegreeDays() example
# get fomatted temp series (no gaps, interpolated mean temps)
# d <- formatTempSeq(subset(data.frame(saved_station_temps), Station == "Nugget Point Aws"))
# envtemp <- d$Tmean
# # quick checks
# all(cumsum(envtemp) == getDegreeDays(envtemp, developtemp = 0, timeout = 0L)[["degdays"]])
# all(getDegreeDays(envtemp, developtemp = 100, timeout = 0L)[["degdays"]] == 0L)
# # check timeout - should reset from 5 degree days after being at or below
# # developtemp for five consecutive days
# getDegreeDays(c(15, rep(10, 9)), developtemp = 10, timeout = 5)
# getDegreeDays(c(15, rep(10, 4), 15, rep(10, 4)), developtemp = 10, timeout = 5)
#-------------------------------------------------------------------------#

#' TODO
#'
#' @param temp_hist TODO
#' @param lookback TODO
#' @param series_ending TODO
#' @param developtemp TODO
#' @param timeout TODO
#'
#' @export
#'
#' @examples
#' temp_hist <- formatTempSeq(subset(data.frame(saved_station_temps),
#'                                   Station == "Nugget Point Aws"))
#' meanDegDays <- getCalendarDegreeDays(temp_hist, developtemp = 12.97, timeout = 30L)
#' head(meanDegDays, 15)
#'
#' \dontrun{
#'   # plot historical temperatures, degree days and degree day means
#'   # - blue lines are mean temperatures with a trace for each calendar year
#'   # - black lines are degree days based on historical mean temperatures
#'   # - red lines are the calendar day means for degree days
#'   # - NOTE: the scales for degree days are adjusted to plot alongside mean
#'   #   temperatures and do not correspond to the y axis
#'   library(ggplot2)
#'   temp_hist$Year <- as.numeric(format(temp_hist$Date, "%Y"))
#'   temp_hist <- subset(temp_hist, Year %in% 2016:2020)
#'   temp_hist$degdays <-
#'     getDegreeDays(temp_hist$Tmean, developtemp = 12.97,
#'                   timeout = 30L)[["degdays"]]
#'
#'   ggplot(data = temp_hist, aes(x = as.numeric(factor(calday)), group = Year)) +
#'     geom_line(aes(y = Tmean), color = "blue") +
#'     geom_line(aes(y = degdays/365*25)) +
#'     geom_line(data = meanDegDays, aes(x = as.numeric(factor(calday)),
#'                                       y = meanDegDays/365*25),
#'               color = "red", inherit.aes = FALSE) +
#'     geom_hline(aes(yintercept = 12.97)) + labs(x = "Calendar day") # +
#'   # facet_wrap(~ Year)           ## uncomment to show separate years
#' }
getCalendarDegreeDays <-
  function(temp_hist,
           lookback = "-10 year",
           series_ending = NULL,
           developtemp,
           timeout){

    # temp_hist <- formatTempSeq(subset(data.frame(saved_station_temps), Station == "Dunedin, Musselburgh Ews"))
    # lookback <- "-10 year"
    # series_ending <- NULL
    # developtemp <- 12.97
    # timeout <- 30L

    Date <- NULL

    # use last date in temp series if no series_ending date given
    if(is.null(series_ending)){
      series_ending <- max(temp_hist$Date)
    }

    # get dates starting the specified range
    start_end_dates <- seq(series_ending, by = lookback, length.out = 2)

    # check supplied historical date range matches period to average over
    if(min(temp_hist$Date) > start_end_dates[2]){
      stop("Historical temperatures start after date range for calendar day averages")
    }

    # only keep past temperatures in specified date range
    subdat <- subset(temp_hist, Date >= start_end_dates[2])

    # get degree days
    subdat$degdays <- getDegreeDays(subdat$Tmean, developtemp, timeout)[["degdays"]]

    # calculate calendar day mean degree days
    calday_means <- stats::aggregate(degdays ~ calday, FUN = mean, data = subdat)

    # rename Tmean column
    calday_means$meanDegDays <- calday_means$degdays
    calday_means$degdays <- NULL

    # replace leap day with mean for 28th Feb
    calday_means$meanDegDays[calday_means$calday == "02-29"] <-
      calday_means$meanDegDays[calday_means$calday == "02-28"]

    attr(calday_means, "start_end_dates") <- start_end_dates

    return(calday_means)
  }

