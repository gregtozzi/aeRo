require(MASS)
require(pracma)

foilPlot <- function(foil, ...) {
  eqscplot(foil$xu[, 1], foil$zu[, 1], type = 'l', bty = 'n', xlab = "x", ylab = "z", ...)
  lines(foil$xl, foil$zl, ...)
  lines(foil$x, foil$zc, col = '#FF1300', ...)
  lines(foil$x, foil$zt, col = '#34AADC', ...)
}

flowPlot <- function(x, z, u, w, New = TRUE) {
  if(New == TRUE) {
    plot.new()
    plot.window(xlim = c(min(x), max(x)), ylim = c(min(z), max(z)))
  }
  if(sum(is.nan(u)) > 0) u[which(is.nan(u))] <- 0
  if(sum(is.nan(w)) > 0) w[which(is.nan(w))] <- 0
  quiver(x, z, u, w)
}

pressPlot <- function(x, z, u, w, ...) {
  if(sum(is.nan(u)) > 0) u[which(is.nan(u))] <- 0
  if(sum(is.nan(w)) > 0) w[which(is.nan(w))] <- 0
  P <- sqrt((u ^ 2) + (w ^ 2))
  contour(x = x, y = z, z = P)
}

thinDisc <- function(X, Z, N) {
  # Discretize the thin foil into N panels
  discFoil <- approx(X, Z, n = N + 1)
  dx <- diff(discFoil$x)
  dz <- diff(discFoil$y)
  panelLength <- sqrt(dx ^ 2 + dz ^ 2)
  
  # Determine the angle and length of each panel
  theta <- asin(dz / panelLength)
  pct <- panelLength * cos(theta)
  pst <- panelLength * sin(theta)
  
  # Set the quarter chord and three-quarter chord points
  x <- 0.25 * pct + discFoil$x[1:N]
  z <- 0.25 * pst + discFoil$y[1:N]
  xc <- 0.75 * pct + discFoil$x[1:N]
  zc <- 0.75 * pst + discFoil$y[1:N]
  
  # Compute the normal vectors
  nx <- - sin(theta)
  nz <- cos(theta)
  
  return(list(x = x, z = z, xc = xc, zc = zc, nx = nx, nz = nz))
}

thinPlot <- function(X, Z, foil) {
  eqscplot(X, Z, type = 'l', bty = "n",
           col = "#DBDDDE", lwd = 2)
  points(foil$x, foil$z, pch = 19, col = "#FF3B30")
  points(foil$xc, foil$zc, pch = 19, col = "#34AADC")
  arrows(foil$xc, foil$zc, foil$xc + foil$nx, foil$zc + foil$nz,
         col = "#34AADC")
}