#' @importFrom rlang .data
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
  Tmin.est <- stats::approx(y = d$Tmin, x = seq_len(nrow(d)),
                     xout = seq_len(nrow(d)), na.rm = TRUE)$y
  d$Tmin[is.na(d$Tmin)] <- Tmin.est[is.na(d$Tmin)]

  Tmax.est <- stats::approx(y = d$Tmax, x = seq_len(nrow(d)),
                     xout = seq_len(nrow(d)), na.rm = TRUE)$y
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
#' @return
#' @export
#'
#' @examples
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
    my.stations = clifro::cf_station(15752)       # 15752 - Musselburgh station ID

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
  temp_delta <- mean(temp_lookback$Tmean - temp_lookback$Tmean_10yr)

  message(sprintf("%s-day difference in mean and mean calendar day temperatures = %0.3f degreesC", lookback_days, temp_delta))

  temp_proj$Tmean <- temp_proj$Tmean_10yr + temp_delta

  return(temp_proj)
}



