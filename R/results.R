plot_popn <- function(resdf,
                      selectPopn = c("L", "L_1", "L_2", "L_3", "L_4", "L_5", "M"),
                      include_temp = TRUE){

  if(include_temp){
    selectPopn <- c(selectPopn, "Tmean")
  }

  # append last recorded temperature (removes gap when drawing geom_line)
  if(any(resdf$source != "projected")){
    last.recorded <- resdf[max(which(resdf$source != "projected")),]
    last.recorded$source <- "projected"
    resdf <- rbind(resdf, last.recorded)
  }

  # stack selected columns
  d <- tidyr::pivot_longer(resdf, cols = dplyr::all_of(selectPopn))
  # respect plotting order of selectPopn
  d$name <- factor(d$name, levels = unique(d$name))
  # create factor for recorded vs projected temperature estimates
  d$Temperature <- factor(ifelse(d$source %in% "projected", "projected", "recorded"),
                          levels = c("recorded", "projected"))


  gg <-
    ggplot2::ggplot(d, ggplot2::aes(y = .data$value, x = .data$Date,
                                    color = .data$name, size = .data$Temperature)) +
    ggplot2::geom_line() +
    ggplot2::scale_size_manual(breaks = c("recorded", "projected"), values = c(0.8,0.3))

  if(!include_temp){
    gg <- gg + ggplot2::labs(y = "Number per x")
  } else {
    gg <- gg + ggplot2::facet_grid(name == "Tmean" ~ . , scale = "free_y",
                                   labeller = ggplot2::as_labeller(c('TRUE' = "Temperature", 'FALSE' = "Population size")),
                                   switch = "y") +
      ggplot2::labs(y = NULL)
  }

  return(gg)
}


plot_popn_years <- function(resdf,
                       selectPopn = c("L", "L_1", "L_2", "L_3", "L_4", "L_5", "M")){

  # append last recorded temperature (removes gap when drawing geom_line)
  if(any(resdf$source != "projected")){
    last.recorded <- resdf[max(which(resdf$source != "projected")),]
    last.recorded$source <- "projected"
    resdf <- rbind(resdf, last.recorded)
  }

  # stack selected columns
  d <- tidyr::pivot_longer(resdf, cols = dplyr::all_of(selectPopn))
  # respect plotting order of selectPopn
  d$name <- factor(d$name, levels = unique(d$name))
  # calendar year
  d$Year <- as.numeric(format(d$Date, "%Y"))
  # Start years in July
  d$Year <- d$Year - (as.numeric(format(d$Date, "%m")) < 7)
  # drop any leap days
  d <- d[!grepl("02-29", format(d$Date, "%m-%d")),]
  # align start of each year with first date in sequence
  d <- dplyr::mutate(dplyr::group_by(d, .data$Year),
                     Date = min(d$Date) + (.data$Date - min(.data$Date)))
  # make year discrete
  d$Year <- factor(d$Year)
  # create factor for recorded vs projected temperature estimates
  d$Temperature <- factor(ifelse(d$source %in% "projected", "projected", "recorded"),
                          levels = c("recorded", "projected"))

  gg <- ggplot2::ggplot(d, ggplot2::aes(y = .data$value, x = .data$Date,
                                        color = .data$Year, group = .data$Year,
                                        size = .data$Temperature)) +
    ggplot2::geom_line() +
    ggplot2::labs(y = "Number per x") +
    ggplot2::scale_x_date(date_labels = "%b") +
    ggplot2::scale_color_brewer(palette = "Dark2") +
    ggplot2::scale_size_manual(breaks = c("recorded", "projected"), values = c(0.8,0.3))


  if(length(selectPopn) == 1){
    return(gg)
  } else {
    return(gg + ggplot2::facet_wrap(~name, scales = "free_y"))
  }
}
