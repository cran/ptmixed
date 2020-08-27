#' Probability mass function plot for a discrete variable
#'
#' This function produces a simple plot of the probability mass function
#' of a discrete variable
#' 
#' @param x the (discrete) variable of interest
#' @param absolute logical. If \code{TRUE} (default) absolute frequencies 
#' are plotted, if \code{FALSE} relative frequencies
#' @param xlim limits for the x axis
#' @param lwd line width
#' @param col color used for the vertical frequency bars
#' @param title plot title
#' @param xlab label for the x axis
#' @param bty box type (default is \code{bty="l"})
#' @param cex.title title font size
#' @param cex.axis font size for the axes
#' @import graphics
#' @export
#' @author Mirko Signorelli
#' @references Signorelli, M., Spitali, P., Tsonaka, R. (2020). Poisson-Tweedie 
#' mixed-effects model: a flexible approach for the analysis of longitudinal RNA-seq
#' data. Statistical Modelling. URL: https://doi.org/10.1177/1471082X20936017
#' @examples
#' pmf(cars$speed)
#' pmf(cars$speed, absolute = FALSE)
#' pmf(cars$speed, lwd = 2, col = 'blue')

pmf = function(x, absolute = T, 
                xlim = NULL, lwd = 1, col = 'black',
                title = NULL, xlab = NULL, bty = 'l',
                cex.title = NULL, cex.axis = NULL) {
  table = table(x)
  x = as.numeric(row.names(table))
  freq = as.numeric(table)
  ylab = 'Absolute frequency'
  if (!absolute) {
    freq = freq / sum(freq)
    ylab = 'Relative frequency'
  }
  if (is.null(xlim)) xlim = c(min(x), max(x))
  plot(NULL, xlim = xlim, ylim = c(0, max(freq)), bty = bty,
       cex.axis = cex.axis, cex.lab = cex.axis, cex.main = cex.title,
       xlab = xlab, ylab = ylab, main = title)
  segments(x0 = x, x1 = x, y0 = 0, y1 = freq, lwd = lwd, col = col)
}
