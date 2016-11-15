# --------------------------------------------------------------
# Monte Carlo Simulation of Radiation Trapping
# Copyright (c) 1988-2016 svengato. All rights reserved.
# --------------------------------------------------------------
library(hexbin) # for hexbinplot()
library(plotrix) # for draw.circle()
# --------------------------------------------------------------

initializeDisplay <- function() {
  # Cell dimensions
  cell.radius <<- 2.5 # cm
  cell.length <<- 7.5 # cm
  # z offset of front view from side view
  # (this leaves one cell radius of space between the two views)
  z.offset <<- cell.length + 2*cell.radius

  # Distances for displaying the photon outside the cell
  z.initial <<- -1.0 # cm, apparent "distance" of initial photon from cell (negative = left)
  d.pretty <<- abs(z.initial) # apparent "distance" of photon after escaping the cell

  # Plot colors
  color.vapor <<- list("K" = "lavender", "Rb" = "mistyrose", "Cs" = "lightblue")
  color.laser <<- "red"
  color.delta.Jz <<- "#FF8000" # orange
  color.num.scatters <<- "#0080FF" # aqua
  color.bluescale <<- colorRampPalette(c("#F0F0FF", "#0000FF"))
}

# --------------------------------------------------------------

# Define the server logic
shinyServer(function(input, output, session) {
  rv <- reactiveValues()
  initializeDisplay()

  # Reset the simulation, clear the statistics
  resetSimulation <- function() {
    # Simulation speed requested by user
    rv$run.state <- "reset"

    # Statistics:
    # n = number of incident photons so far
    # delta.Jz = vector of spin transfers by each incident photon
    # num.scatters = vector of number of scatters for each incident photon
    rv$statistics <- list(n = 0, delta.Jz = c(), num.scatters = c())

    # Alkali vapor at the specified temperature
    atom.name <- isolate(input$alkali.atom)
    temp <- isolate(input$temperature)
    polz <- 0.01*isolate(input$polarization)
    rv$vapor <- alkaliVapor(atom.name, temp, polz)

    output$optical.depth.0 <- renderText(sprintf("Optical depth at line center: %3.2e",
      rv$vapor$concentration*rv$vapor$sigma0*voigt(0, rv$vapor$w.natural)*cell.radius))
  }

  resetSimulation()

  # Determine where to display the photon after it escapes the cell
  # (= d.pretty from the cell wall, in the direction it was traveling)
  getPrettyPosition <- function() {
    x <- tail(rv$photon$x, 2)
    y <- tail(rv$photon$y, 2)
    z <- tail(rv$photon$z, 2)

    if (x[2] > cell.radius + d.pretty) {
      f <- (cell.radius + d.pretty - x[1])/(x[2] - x[1])
      x[2] <- cell.radius + d.pretty
      y[2] <- y[1] + f*(y[2] - y[1])
      z[2] <- z[1] + f*(z[2] - z[1])
    } else if (-x[2] > cell.radius + d.pretty) {
      f <- -(cell.radius + d.pretty + x[1])/(x[2] - x[1])
      x[2] <- -(cell.radius + d.pretty)
      y[2] <- y[1] + f*(y[2] - y[1])
      z[2] <- z[1] + f*(z[2] - z[1])
    }

    if (y[2] > cell.radius + d.pretty) {
      f <- (cell.radius + d.pretty - y[1])/(y[2] - y[1])
      x[2] <- x[1] + f*(x[2] - x[1])
      y[2] <- cell.radius + d.pretty
      z[2] <- z[1] + f*(z[2] - z[1])
    } else if (-y[2] > cell.radius + d.pretty) {
      f <- -(cell.radius + d.pretty + y[1])/(y[2] - y[1])
      x[2] <- x[1] + f*(x[2] - x[1])
      y[2] <- -(cell.radius + d.pretty)
      z[2] <- z[1] + f*(z[2] - z[1])
    }
    
    if (z[2] > cell.length + d.pretty) {
      f <- (cell.length + d.pretty - z[1])/(z[2] - z[1])
      x[2] <- x[1] + f*(x[2] - x[1])
      y[2] <- y[1] + f*(y[2] - y[1])
      z[2] <- cell.length + d.pretty
    } else if (-z[2] > d.pretty) {
      f <- -(d.pretty + z[1])/(z[2] - z[1])
      x[2] <- x[1] + f*(x[2] - x[1])
      y[2] <- y[1] + f*(y[2] - y[1])
      z[2] <- -d.pretty
    }

    c(x[2], y[2], z[2])
  }

  # Test whether the photon is inside the cell
  isInsideCell <- function() {
    if (rv$statistics$n == 0) return(FALSE)

    x <- tail(rv$photon$x, 1)
    y <- tail(rv$photon$y, 1)
    z <- tail(rv$photon$z, 1)
    (z > 0 && z < cell.length && (x^2 + y^2 < cell.radius^2))
  }

  # Begin tracking the next incident photon
  nextPhoton <- function() {
    r0 <- cell.radius*sqrt(runif(1))
    phi0 <- runif(1, min = 0, max = k.2pi)
    rv$photon <- newPhoton(freq = input$detuning, f.up = 1,
      x = r0*cos(phi0), y = r0*sin(phi0), z = z.initial, theta = 0, phi = 0)
      # note phi is undefined because theta = 0

    rv$statistics$delta.Jz <- append(rv$statistics$delta.Jz, 0)
    rv$statistics$num.scatters <- append(rv$statistics$num.scatters, 0)
    rv$statistics$n <- rv$statistics$n + 1
  }

  # Propagate the current photon (decide where it gets absorbed)
  propagatePhoton <- function() {
    rv$photon <- propagate(rv$photon, rv$vapor)

    # test whether the photon escaped from the cell
    if (!isInsideCell()) {
      # clean up "position" of escaped photon
      xyz.pretty <- getPrettyPosition()
      ns <- 2 + tail(rv$statistics$num.scatters, 1)
      rv$photon$x[ns] <- xyz.pretty[1]
      rv$photon$y[ns] <- xyz.pretty[2]
      rv$photon$z[ns] <- xyz.pretty[3]
    }
  }

  # Scatter the current photon
  scatterPhoton <- function() {
    results <- scatter(rv$photon, rv$vapor)
    rv$photon <- results$photon

    # update the statistics
    rv$statistics$delta.Jz[rv$statistics$n] <- rv$statistics$delta.Jz[rv$statistics$n] + results$delta.Jz
    rv$statistics$num.scatters[rv$statistics$n] <- rv$statistics$num.scatters[rv$statistics$n] + 1
  }

  # Track one scatter
  scatterOnce <- function() {
    if (isInsideCell()) {
      scatterPhoton()
    } else {
      nextPhoton()
    }
    propagatePhoton()
  }

  # Track all scatters for an incident photon until it escapes
  scatterUntilEscaped <- function() {
    # alternatively, scatterOnce()
    if (!isInsideCell()) {
      nextPhoton()
      propagatePhoton()
    }
    while (isInsideCell()) {
      scatterPhoton()
      propagatePhoton()
    }
  }

  # Reactive event handlers
  observeEvent(input$alkali.atom, resetSimulation())
  observeEvent(input$polarization, resetSimulation())
  observeEvent(input$temperature, resetSimulation())

  observeEvent(input$detuning, {
#    resetSimulation()
    output$optical.depth.y <- renderText(sprintf("Optical depth at laser detuning: %3.2e",
      rv$vapor$concentration*rv$vapor$sigma0*voigt(input$detuning, rv$vapor$w.natural)*cell.radius))
    rv$run.state <- "reset" # TODO: comment this out? (just tuning the laser)
  })

  observeEvent(input$run3, {
    rv$run.state <- "pause"
    rv$run.state <- "run3"
  })

  observeEvent(input$run2, {
    rv$run.state <- "pause"
    rv$run.state <- "run2"
  })

  observeEvent(input$run1, {
    rv$run.state <- "pause"
    rv$run.state <- "run1"
  })

  observeEvent(input$pause, {
    rv$run.state <- "pause"
  })

  observeEvent(input$reset, resetSimulation())

  # Visualization of results
  drawCell <- function() {
    plot(x = c(0, z.offset + cell.radius), y = c(-2*cell.radius, 2*cell.radius), type = "n", asp = 1,
      xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    text(0.5*(z.offset + cell.radius), 1.5*cell.radius,
      paste("Cylindrical cell containing", rv$vapor$symbol, "vapor"))
    # side view
    rect(0, -cell.radius, cell.length, cell.radius, col = rv$vapor$color)
    text(0.5*cell.length, -1.5*cell.radius, "Side View")
    # front view
    draw.circle(z.offset, 0, radius = cell.radius, col = rv$vapor$color)
    text(z.offset, -1.5*cell.radius, "Front View")
  }

  drawPhotonTrajectories <- function() {
    # side view
    # TODO: replace rv$photon$z with is.na(rv$photon$z), etc?
    lines(rv$photon$z, rv$photon$y, col = color.laser)
    # front view
    lines(rv$photon$x + z.offset, rv$photon$y, col = color.laser)
  }

  output$plot.cell <- renderPlot({
    drawCell()

    if (rv$run.state == "reset") {
      return()
    } else if (rv$run.state == "pause") {
      drawPhotonTrajectories()
    } else if (rv$run.state == "run1") {
      isolate(scatterOnce())
      drawPhotonTrajectories()
    } else if (rv$run.state == "run2") {
      isolate(scatterUntilEscaped())
      drawPhotonTrajectories()
    } else if (rv$run.state == "run3") {
      isolate({
        drawCell()
        scatterUntilEscaped()
        drawPhotonTrajectories()
      })
      invalidateLater(100, session)
    }
  })

  output$label.statistics <- renderUI(
    # TODO: density plot of ∆Jz for this photon
    HTML(
      ifelse(rv$statistics$n == 0, "This photon:",
        sprintf("This photon: &Delta;Jz = %3.2f, # Scatters = %d",
          tail(rv$statistics$delta.Jz, 1), tail(rv$statistics$num.scatters, 1)
        )
      )
    )
  )

  output$label.mean.statistics <- renderUI(HTML(
    ifelse(rv$statistics$n == 0 || tail(rv$statistics$num.scatters, 1) == 0, "",
      sprintf("All %d photons: &lt;&Delta;Jz&gt; = %3.2f, &lt;# Scatters&gt; = %2.1f",
        rv$statistics$n, mean(rv$statistics$delta.Jz), mean(rv$statistics$num.scatters)
      # sprintf("All %d photons: Median &Delta;Jz = %3.2f, Median # Scatters = %2.1f",
      #   rv$statistics$n, median(rv$statistics$delta.Jz), median(rv$statistics$num.scatters)
      )
    )
  ))

  output$plot.absorption.profile <- renderPlot({
    x <- seq(-0.5, 5.5, by = 0.01)
    y <- input$detuning
    b <- rv$vapor$w.natural
    pdf_x <- b/(pi*k.sqrtpi*voigt(y, b)) * exp(-x^2)/((x - y)^2 + b^2)
    plot(x, pdf_x, type = "l", log = "y", xlab = "Frequency (in Doppler widths)", ylab = "Absorption PDF")
  })

  # Density plots of ∆Jz and number of scatters for all photons, on same chart with different axes
  # TODO: allow user to set the bandwidths?
  # (not currently used)
  output$plot.delta.Jz.and.num.scatters <- renderPlot({
    if (rv$statistics$n > 1) {
      plot(density(rv$statistics$delta.Jz), xlab = expression(paste(Delta, "Jz", sep = "")),
        ylab = "", main = "", col = color.delta.Jz, axes = FALSE, ann = FALSE)
      at <- axTicks(1)
      axis(1, at = at, labels = at, col = color.delta.Jz)
      
      par(new = TRUE)
      plot(density(rv$statistics$num.scatters), xlab = "Scatters", ylab = "Density", main = "",
        col = color.num.scatters, axes = FALSE, ann = FALSE)
      at <- axTicks(1)
      axis(3, at = at, labels = at, col = color.num.scatters)
    }
  })

  # Density plot of ∆Jz for all photons
  # TODO: allow user to set the bandwidth?
  output$plot.delta.Jz <- renderPlot({
    if (rv$statistics$n > 1) {
      # no need for ylab = "Relative Frequency" as the ∆Jz plot y axis already has that label
      plot(density(rv$statistics$delta.Jz), main = "", xlab = expression(paste(Delta, "Jz", sep = "")),
        ylab = "", yaxt = "n", col = color.delta.Jz)
      # rug(jitter(rv$statistics$delta.Jz, amount = 0.01), col = color.delta.Jz)
      abline(v = 0.3333, col = color.delta.Jz, lty = "dashed")
    }
  })

  # Density plot of number of scatters for all photons
  # TODO: allow user to set the bandwidth?
  output$plot.num.scatters <- renderPlot({
    if (rv$statistics$n > 1) {
      plot(density(rv$statistics$num.scatters), main = "", xlab = "Scatters", ylab = "Frequency",
        yaxt = "n", col = color.num.scatters)
      # plot(density(rv$statistics$num.scatters, bw = 1.0), main = "", xlab = "Scatters", ylab = "Frequency", col = color.num.scatters)
      # rug(jitter(rv$statistics$num.scatters, amount = 0.2), col = color.num.scatters)
      abline(v = 1, col = color.num.scatters, lty = "dashed")
    }
  })

  # Scatter plot of ∆Jz v. number of scatters for all photons
  # (not currently used)
  output$plot.delta.Jz.v.num.scatters <- renderPlot({
    if (rv$statistics$n > 1) {
      hexbinplot(rv$statistics$delta.Jz ~ rv$statistics$num.scatters, main = "", xlab = "# Scatters",
        ylab = expression(paste(Delta, "Jz", sep = "")),
        aspect = 1, colramp = color.bluescale, colorkey = FALSE)
    }
  })
})

# --------------------------------------------------------------
