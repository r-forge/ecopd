# Quick plot of some of the different PD values
caterpillar <- function(data, center=c("mean", "median"), sort=TRUE,
  ...) {

  method <- attr(data, "method")
  center <- match.arg(center)
  if(class(data)=="matrix") {
    data <- as.data.frame(data)
  }
  if(sort) {
    data <- data[order(data$pd.obs),]
  }

  # set up plotting area
  ind <- seq_len(nrow(data))
  xlab <- paste("PD",
    if(!is.null(method)) paste(" (", method, ")", sep="") else NULL,
    sep="")
  plot(data[, "0.025"], ind, type = "n", yaxt = "n", xlab = xlab,
    ylab = NA, xlim = range(c(data[,"0.025"], data[, "0.975"])),
    ...)
  axis(2, at = ind, labels = rownames(data), las = 2)

  # plot 95% bounds
  arrows(data[, "0.025"], ind, data[, "0.975"], ind, code = 3,
    angle = 90, length = 0.01)

  # plot mean or median 
  if (center == "mean") {
    points(data[, "mean"], ind, pch = 3)
  } else if (center == "median") {
    points(data[, "median"], ind, pch = 3)
  }

  # plot actual pd values
  points(data[, "pd.obs"], ind, pch = 1, col = "red")

}
