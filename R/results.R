plotModOut <- function(resdf,
                       selectPopn = c("L", "L_1", "L_2", "L_3", "L_4", "L_5", "M")){

  # convert to long format
  d <- tidyr::pivot_longer(resdf[c("time", "Date", "Tmean", selectPopn)], cols = all_of(selectPopn))
  # respect plotting order of selectPopn
  d$name <- factor(d$name, levels = unique(d$name))

  gg <- ggplot2::ggplot(d, ggplot2::aes(y = value, x = Date, color = name)) +
    ggplot2::geom_line() +
    ggplot2::labs(y = "Number per x")

  return(gg)
}
