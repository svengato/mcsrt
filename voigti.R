# --------------------------------------------------------------
# Monte Carlo Simulation of Radiation Trapping
# Copyright (c) 1988-2018 svengato. All rights reserved.
# --------------------------------------------------------------
# Voigt profile and related functions

# constants
k.sqrtpi <- sqrt(pi)

# --------------------------------------------------------------

# Dawson integral
# 						          x
# dawson(x) = exp(-x^2) ∫ exp(t^2) dt
#                       0
# by the Taylor expansion of Lohmander & Rittsten,
# Kungl. Fysiogr. Sällsk. Lund Förh. [28], 45 (1958).
# -> n = 64 is more stable than n = 50
dawson <- function(x, n = 64) {
  a0 <- 0
  a1 <- h <- min(x, 2)
  dawson_x <- a0 + a1
  for (s in 2:n) {
    a2 <- -2*h/s * h*a0
    dawson_x <- dawson_x + a2
    a0 <- a1
    a1 <- a2
  }
  
  if (x > 2) {
    a0 <- dawson_x
    h <- x - 2
    a1 <- h*(1 - 4*dawson_x)
    dawson_x <- a0 + a1
    for (s in 2:n) {
      a2 <- -2*h/s * (2*a1 + h*a0)
      dawson_x <- dawson_x + a2
      a0 <- a1
      a1 <- a2
    }
  }
  
  dawson_x
}

# --------------------------------------------------------------

# Voigt pdf
# y = (frequency detuning from line center)/(Doppler line width)
# b = (natural line width)/(Doppler line width)
#             b     ∞  exp(-x^2) dx
# V(y,b) = -------  ∫ -------------
#          π^(3/2) -∞ b^2 + (y-x)^2

# There are two different regions of integration:
# (1) Taylor expansion in y (for y < 4)
#     F. Hjerting, Astrophys. J. [88], 508 (1938).
# (2) 8-term Gauss-Hermite quadrature (for y ≥ 4)
#     Abramowitz & Stegun, pp. 890 & 924.
# These approximations are accurate for b < 0.1 (TODO: quantify how accurate they are.)

voigt <- function(yy, b) {
  # x[i] are zeroes of the Hermite polynomial H8(x)
  x <- c(0.381187, 1.157194, 1.981657, 2.930637)
  x <- c(x, -x)
  # and w[i] are its Gaussian weighting factors
  w <- c(0.661147, 0.207802, 0.017078, 0.000199604)
  w <- rep(w, 2)
  
  # vectorize it in order to enable integration
  sapply(yy, FUN = function(y) {
    y <- abs(y) # sign-symmetric
    if (y < 4) {
      # Taylor series expansion
      voigt_yb <- exp(-y*y) - 2*b/k.sqrtpi * (1 - 2*y*dawson(y))
    } else {
      # Gauss-Hermite quadrature
      voigt_yb <- sum(w/(b^2 + (y - x)^2)) * b/pi
    }
    voigt_yb
  })/k.sqrtpi
}

# integrate(voigt, -Inf, -0.5, b = 0.006)
# #0.2413713 with absolute error < 3.8e-07
# integrate(voigt, -Inf, 0, b = 0.006)
# #0.5 with absolute error < 1.8e-06
# integrate(voigt, -Inf, Inf, b = 0.006)
# #0.9999999 with absolute error < 3.7e-06

# --------------------------------------------------------------

# Voigt integrand, for the parallel velocity distribution:
# the idea is to find vz given the frequency detuning y and a
# random number Q in [0,1), where
#       1       b     vz  exp(-x^2) dx
# Q = ------ -------  ∫  -------------
#     V(y,b) π^(3/2) -∞  b^2 + (y-x)^2	
# We must numerically invert this equation.

# Probability density/distribution function (PDF) of *non*-normalized Voigt integrand
dvoigti <- function(x, y, b) {
  exp(-x^2)/(b^2 + (y - x)^2) # *b/(pi*sqrt(pi)*voigt(y, b))
}

# Cumulative distribution function (CDF) of *normalized* Voigt integrand
pvoigti <- function(x, y, b) {
  # fnorm = normalization factor (outside of dvoigti for speed)
  fnorm <- b/(pi*k.sqrtpi*voigt(y, b))
  rfnorm <- 1/fnorm

  pv <- sapply(x, FUN = function(xi) {
    # Since the Lorentzian contribution near x = y is so spiky,
    # approach it from the left xor right side
    if (xi < y) {
      integrate(dvoigti, -Inf, xi, y = y, b = b)$value
    } else {
      rfnorm - integrate(dvoigti, xi, Inf, y = y, b = b)$value
    }
    # changing "subdivisions =" in integrate() has no noticeable effect
  })
  pv*fnorm
}

# Quantile function = inverse of the CDF
# TODO: recover if f() of interval end points not of opposite sign, etc
qvoigti <- function(q, y, b) {
  qv <- uniroot(function(x) pvoigti(x, y, b) - q, interval = c(-2.5, 8.5))
  qv$root
}

# Random deviates from this distribution
rvoigti <- function(n, y, b) {
  # TODO: vectorize
  # rv <- qvoigti(runif(n), y, b)
  qq <- runif(n)
  rv <- numeric(n)
  for (i in 1:n) rv[i] <- qvoigti(qq[i], y, b)
  rv
}

# vz <- rvoigti(512, 3, 0.006)
# plot(density(vz), xlim = c(-3, 5), main = "")
# rug(vz)

# --------------------------------------------------------------
