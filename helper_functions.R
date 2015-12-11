require(MASS)

foilPlot <- function(foil, ...) {
  eqscplot(foil$xu[, 1], foil$zu[, 1], type = 'l', bty = 'n', xlab = "x", ylab = "z", ...)
  lines(foil$xl, foil$zl, ...)
  lines(foil$x, foil$zc, col = '#FF1300', ...)
  lines(foil$x, foil$zt, col = '#34AADC', ...)
}