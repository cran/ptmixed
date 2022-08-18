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
#' @examples
#' pmf(cars$speed)
#' pmf(cars$speed, absolute = FALSE)
#' pmf(cars$speed, lwd = 2, col = 'blue')

pmf = function(x, absolute = T, 
                xlim = NULL, lwd = 1, col = 'black',
                title = NULL, xlab = 'x', bty = 'l',
                cex.title = NULL, cex.axis = NULL) {
  table = table(x)
  xvals = as.numeric(row.names(table))
  freq = as.numeric(table)
  ylab = 'Absolute frequency'
  if (!absolute) {
    freq = freq / sum(freq)
    ylab = 'Relative frequency'
  }
  if (is.null(xlim)) xlim = c(min(xvals), max(xvals))
  plot(NULL, xlim = xlim, ylim = c(0, max(freq)), bty = bty,
       cex.axis = cex.axis, cex.lab = cex.axis, cex.main = cex.title,
       xlab = xlab, ylab = ylab, main = title)
  segments(x0 = xvals, x1 = xvals, y0 = 0, y1 = freq, lwd = lwd, col = col)
}
