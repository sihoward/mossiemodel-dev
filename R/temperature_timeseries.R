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

  # subset the 10 years of past temperatures
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





