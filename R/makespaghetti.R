#' Generate a spaghetti plot to visualize longitudinal data
#'
#' A spaghetti plot, or trajectory plot, is a plot that allows
#' to compare across individuals or groups the trajectories 
#' of a longitudinal outcome
#' 
#' @param x the time variable (numeric vector)
#' @param y the longitudinal outcome (numeric vector)
#' @param id the subject indicator
#' @param group the group that each subject belongs to (optional, 
#' do not specify if not relevant)
#' @param data a data frame containing x, y, id and optionally group
#' @param col a vector of colors (optional)
#' @param lty line type
#' @param lwd line width
#' @param title plot title
#' @param xlab label for the x axis
#' @param ylab label for the y axis
#' @param pch dot type
#' @param legend.title legend title
#' @param ylim limits for the y axis
#' @param cex.axis font size for the axes
#' @param cex.title title font size
#' @param cex.lab font size for axis labels
#' @param cex.leg font size for the legend
#' @param margins use this argument if you want to overwrite the default 
#' function margins
#' @param legend.inset moves legend more to the left / right (default is -0.3)
#' @param legend.space interspace between lines in the legend (default is 1)
#' @import graphics
#' @export
#' @author Mirko Signorelli
#' @references Signorelli, M., Spitali, P., Tsonaka, R. (2020). Poisson-Tweedie 
#' mixed-effects model: a flexible approach for the analysis of longitudinal RNA-seq
#' data. Statistical Modelling. URL: https://doi.org/10.1177/1471082X20936017
#' @examples
#' \donttest{
#' # generate example data
#' set.seed(123)
#' n = 12; t = 6
#' id = rep(1:n, each = t)
#' rand.int = rep(rnorm(n, sd = 0.5), each = t)
#' group = rep(c(0,1), each = n*t/2)
#' time = rep(0:(t-1), n)
#' offset = rnorm(n*t, sd = 0.3)
#' beta = c(3, 0, 0.1, 0.3)
#' X = model.matrix(~group + time + group*time)
#' mu = 2^(X %*% beta + rand.int + offset)
#' y = rpois(n*t, lambda = mu)
#' group = ifelse(group == 0, 'control', 'treatment')
#' data.long = data.frame(y, group, time, id, offset)
#' rm(list = setdiff(ls(), 'data.long'))
#' 
#' # create plot
#' make.spaghetti(x = time, y, id, group, 
#' data = data.long, title = 'spaghetti plot')
#' }

make.spaghetti = function(x, y, id, group = NULL, data,
                      col = NULL, pch = 16,
                      lty = 1, lwd = 1,
                      title = '', xlab = NA, ylab = NA,
                      legend.title = '', ylim = NULL,
                      cex.axis = 1, cex.title = 1, 
                      cex.lab = 1, cex.leg = 1,
                      margins = NULL,
                      legend.inset = -0.3, legend.space = 1) {
  if (is.na(xlab)) xlab = deparse(substitute(x))
  if (is.na(ylab)) ylab = deparse(substitute(y))
  sort.x = order(data[ , deparse(substitute(x))])
  data = data[sort.x, ]
  x = data[ , deparse(substitute(x))]
  y = data[ , deparse(substitute(y))]
  id = data[ , deparse(substitute(id))]
  if (is.null(ylim)) ylim = range(y)
  if (length(ylim) !=2) warning('ylim should be a vector of length 2')
  check = missing(group)
  if (check) group = NULL
  if (!check) group = data[ , deparse(substitute(group))]
  par(bty = 'l')
  if (is.null(group)) {
    par(mar = c(4, 4, 2, 2))
    if (is.null(col)) col = 'deepskyblue3'
  }
  if (!is.null(group)) {
    par(mar = c(4, 4, 2, 7))
    nlevs = length(unique(group))
    if (is.null(col)) {
      palette = c("#1F78B4", "#33A02C", "#FB9A99", "#E31A1C","#A6CEE3", "#B2DF8A", 
      "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928")
      # this palette is a reordered version of brewer.pal(12, 'Paired') 
      # from RColorBrewer; the yellow was removed because it is hardly visible
      if (nlevs > 11) palette = rep(palette, 50)
      palette = palette[1:nlevs]
    }
    if (!is.null(margins)) par(mar = margins)
    if (!is.null(col)) palette = col
    group = as.factor(group)
    group.names = levels(group)
    col = group
    levels(col) = palette
    col = as.character(col)
  }
  plot(y ~ x, col = col, pch = pch, main = title, 
       xlab = xlab, ylab = ylab, ylim = ylim,
       cex.axis = cex.axis, cex.lab = cex.lab)
  data.long = data.frame(x, y, id, col, stringsAsFactors = F)
  seg.list = names(which(table(data.long$id) >= 2))
  nsegm = length(seg.list)
  for (i in 1:nsegm) {
    subset = data.long[data.long$id == seg.list[i],]
    for (j in 1:(nrow(subset)-1)) {
      segments(x0 = subset$x[j], x1 = subset$x[j+1], 
               y0 = subset$y[j], y1 = subset$y[j+1], 
               col = subset$col[j], lwd = lwd)
    }
  }
  # add legend if multiple groups are present
  if (!is.null(group)) {
    par(xpd = T)
    legend(x = "right", inset=c(legend.inset, 0), levels(group), 
           title = legend.title, pch = pch,
           y.intersp = legend.space,
           col = palette, bty = 'n', cex = cex.leg)
  }
}

