# --------------------------------------------------------------
# Monte Carlo Simulation of Radiation Trapping
# Copyright (c) 1988-2018 svengato. All rights reserved.
# --------------------------------------------------------------
# Physical constants and methods

k.2pi <- 2*pi
k.amu <- 1.67e-24	# atomic mass unit, in g
k.c <- 3.00e+10	# speed of light, in cm/sec
k.osc <- 1/3 # D1 oscillator strength, roughly
k.kB <- 1.38e-16 # Boltzmann's constant, in erg/K
k.re <- 2.81e-13 # classical electron radius, in cm
k.v.thermal <- sqrt(0.5)

# import the Voigt integrand probability functions
source("voigti.R")
# --------------------------------------------------------------

atom.K <- list(
  symbol = "K",
  isotope = 39,
  lifetime = 26.0e-9, # excited state lifetime, in seconds
  wavelength = 770.1e-7, # or 769.9e-7 # D1 absorption wavelength, in cm
  concentration = function(temp) 10^(11.57 - 4964/temp - log10(k.kB*temp)) # in atoms/cm^3
)
atom.Rb <- list(
  symbol = "Rb",
  isotope = 85,
  lifetime = 28.0e-9,
  wavelength = 794.8e-7,
  concentration = function(temp) 10^(10.55 - 4132/temp - log10(k.kB*temp))
)
atom.Cs <- list(
  symbol = "Cs",
  isotope = 133,
  lifetime = 34.0e-9,
  wavelength = 894.4e-7,
  concentration = function(temp) 10^(14.178 - 4041/temp - 1.35*log10(temp) - log10(k.kB*temp))
)

alkali.atoms <- list("Potassium (K)" = atom.K, "Rubidium (Rb)" = atom.Rb, "Cesium (Cs)" = atom.Cs)

alkaliVapor <- function(atom.name, temperature, polarization) {
  atom <- alkali.atoms[[atom.name]]

  # Recompute all vapor quantities at the specified temperature
  mass <- atom$isotope*k.amu
  v.doppler <- sqrt(2*k.kB*temperature/mass)
  w.doppler <- v.doppler/atom$wavelength # Doppler width

  vapor <- list(
    symbol = atom$symbol,
    # atomic concentration, in atoms/cm^3
    concentration = atom$concentration(temperature),
    # absorption cross section, in cm^2
    sigma0 = k.re*k.c*k.osc*k.sqrtpi/w.doppler,
    # ratio of natural and Doppler linewidths (dimensionless)
    w.natural = 1/(4*pi*atom$lifetime*w.doppler),
    # vapor polarization (2*<Jz>)
    polarization = polarization,
    # display color for the vapor
    color = color.vapor[[atom$symbol]]
  )
  vapor
}

# --------------------------------------------------------------

newPhoton <- function(freq, f.up, x, y, z, theta, phi) {
  s <- ifelse(runif(1) < f.up, 1, -1)
  photon <- list(frequency = freq, helicity = s,
    # "position" where created or scattered
    x = x, y = y, z = z,
    # direction
    theta = theta, phi = phi
  )
  photon
}

# Translate the current photon to a new location
# (where it either gets absorbed or escapes from the cell)
propagate <- function(photon, vapor) {
  cos.theta <- cos(photon$theta)
  sin.theta <- sin(photon$theta)
  cos.phi <- cos(photon$phi)
  sin.phi <- sin(photon$phi)
  # sigma = absorption cross section
  sigma <- vapor$sigma0*voigt(photon$frequency, vapor$w.natural)*
    (1 - photon$helicity*vapor$polarization*cos.theta)

  if (length(photon$x) == 1) {
    # always absorb the initial photon
    q.max <- 1 - exp(-vapor$concentration*sigma*cell.length)
    distance <- -log(1 - runif(1, min = 0, max = q.max))/(vapor$concentration*sigma)
    z.new <- distance*cos.theta
  } else {
    # subsequent scattered photons may escape
    distance <- rexp(1, rate = vapor$concentration*sigma)
    z.new <- tail(photon$z, 1) + distance*cos.theta
  }
  x.new <- tail(photon$x, 1) + distance*sin.theta*cos.phi
  y.new <- tail(photon$y, 1) + distance*sin.theta*sin.phi
  # append() is faster than c()?
  # To do: replace append() or c() with a vector initialized to all NA?
  photon$x <- append(photon$x, x.new)
  photon$y <- append(photon$y, y.new)
  photon$z <- append(photon$z, z.new)

  photon
}

# Return the velocity components of the atom that scatters this photon
getAtomVelocity <- function(photon, vapor) {
  # v.zeta is the velocity component in the direction of the photon's momentum,
  # v.xi and v.eta are transverse to it
  v.transverse <- rnorm(2, mean = 0, sd = k.v.thermal)
  v.xi <- v.transverse[1]
  v.eta <- v.transverse[2]
  v.zeta <- rvoigti(1, photon$frequency, vapor$w.natural)

  # convert to the cell's coordinate system
  cos.theta <- cos(photon$theta)
  sin.theta <- sin(photon$theta)
  cos.phi <- cos(photon$phi)
  sin.phi <- sin(photon$phi)
  vx <- v.xi*cos.theta*cos.phi - v.eta*sin.phi + v.zeta*sin.theta*cos.phi
  vy <- v.xi*cos.theta*sin.phi + v.eta*cos.phi + v.zeta*sin.theta*sin.phi
  vz <- -v.xi*sin.theta + v.zeta*cos.theta

  velocity <- list(vx = vx, vy = vy, vz = vz, v.zeta = v.zeta, v = sqrt(v.xi^2 + v.eta^2 + v.zeta^2))
  velocity
}

scatter <- function(photon, vapor) {
  # determine the atom's velocity
  atom <- getAtomVelocity(photon, vapor)

  # transfer angular momentum
  Pz <- vapor$polarization
  s <- photon$helicity
  cos.theta <- cos(photon$theta)
  sin.theta <- sin(photon$theta)
  cos.phi <- cos(photon$phi)
  sin.phi <- sin(photon$phi)
  f.abs <- 1 - s*Pz*cos.theta
  # delta.Jz = angular momentum transferred from the photon to the vapor
  delta.Jz <- (2*s*cos.theta - Pz*(3 - cos.theta^2))/(6*f.abs)
  # Pz.exchanged = ...
  Pz.exchanged <- (s - Pz*cos.theta)/f.abs
  
  # compute the scattered photon's trajectory
  theta.new <- acos(runif(1, min = -1, max = 1))
  phi.new <- runif(1, min = 0, max = k.2pi)
  cos.theta.new <- cos(theta.new)
  sin.theta.new <- sin(theta.new)
  cos.phi.new <- cos(phi.new)
  sin.phi.new <- sin(phi.new)
  # cosine of the angle between the atom and the scattered photon
  cos.chi.a.new <- (atom$vx*sin.theta.new*cos.phi.new + atom$vy*sin.theta.new*sin.phi.new + atom$vz*cos.theta.new)/atom$v
  # cosine of the angle between the incident and scattered photon (TODO: confirm this)
  cos.chi.i.new <- sin.theta*sin.theta.new*(cos.phi*cos.phi.new + sin.phi*sin.phi.new) + cos.theta*cos.theta.new

  # update all values
  photon$frequency <- photon$frequency + (atom$v*cos.chi.a.new - atom$v.zeta)
  f.up.new <- 0.5*(1 + Pz.exchanged*cos.chi.i.new)
  photon$helicity <- ifelse(runif(1) < f.up.new, 1, -1)
  photon$theta <- theta.new
  photon$phi <- phi.new

  # besides the photon, return the angular momentum transfer from this scatter
  list(photon = photon, delta.Jz = delta.Jz)
}

# --------------------------------------------------------------
