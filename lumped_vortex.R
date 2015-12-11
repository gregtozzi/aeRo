thinDisc <- function(X, Z, N) {
  # Discretize the thin foil into N panels
  discFoil <- approx(X, Z, n = N + 1)
  dx <- diff(discFoil$x)
  dz <- diff(discFoil$y)
  panelLength <- sqrt(dx ^ 2 + dz ^ 2)
  theta <- asin(dz / panelLength)
  pct <- panelLength * cos(theta)
  pst <- panelLength * sin(theta)
  x <- 0.25 * pct + discFoil$x[1:N]
  z <- 0.25 * pst + discFoil$y[1:N]
  xc <- 0.75 * pct + discFoil$x[1:N]
  zc <- 0.75 * pst + discFoil$y[1:N]
}

X <- seq(0, pi, length.out = 100)
Z <- 0.8 * sin(X)

eqscplot(X, Z, type = "l")
points(x, z, pch = 19, col = "red")
points(xc, zc, pch = 19, col = "blue")
points(discFoil$x, discFoil$y, pch = 19, cex = 0.5)
