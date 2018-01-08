# --------------------------------------------------------------
# Monte Carlo Simulation of Radiation Trapping
# Copyright (c) 1988-2018 svengato. All rights reserved.
# --------------------------------------------------------------
library(shiny)
library(markdown)

source("physics.R")
# --------------------------------------------------------------

# Define UI for application
shinyUI(fluidPage(
  theme = "mcsrt.css",
  tags$head(
    tags$style("#run3 { background-color: lime }"),
    tags$style("#run2 { background-color: palegreen }"),
    tags$style("#run1 { background-color: gold }"),
    tags$style("#pause { background-color: red; margin-left: 32px }"),
    tags$style("#reset { background-color: cornflowerblue }")
  ),

  # Application title
  titlePanel(HTML("Monte Carlo Simulation of Radiation Trapping<br>in a Dense Optically Pumped Alkali Vapor")),

  tabsetPanel(
    # Simulation tab
    tabPanel("Simulation",
      # Sidebar: various inputs
      sidebarPanel(
        selectInput(inputId = "alkali.atom", label = "Alkali Vapor", choices = names(alkali.atoms),
          selected = "Rubidium (Rb)"),
        sliderInput(inputId = "polarization", label = "Vapor Polarization",
          # TODO: make sure tick values are { 0, 10, 20, ..., 80, 90, 99 }
          min = 0, max = 99, value = 0, step = 1, post = "%"),
        sliderInput(inputId = "temperature", label = "Vapor Temperature",
          # TODO: what should the maximum temperature be?
          min = 293, max = 333, value = 293, step = 5, post = "&deg;K"),
        textOutput("optical.depth.0"),
        sliderInput(inputId = "detuning", label = "Laser Detuning (in Doppler widths)",
          # TODO: what should the initial and maximum detuning be?
          min = 0, max = 5, value = 0, step = 0.1),
        textOutput("optical.depth.y"),
        plotOutput("plot.absorption.profile", height = "256px")
      ),

      # Main panel: display the simulation results
      mainPanel(
        fluidRow(
          actionButton(inputId = "run3", label = "Run"),
          actionButton(inputId = "run2", label = "One Photon"),
          actionButton(inputId = "run1", label = "One Scatter"),
          actionButton(inputId = "pause", label = "Pause"),
          actionButton(inputId = "reset", label = "Reset")
        ),

        plotOutput("plot.cell"),
        htmlOutput("label.statistics"),
        htmlOutput("label.mean.statistics"),
        flowLayout(
          plotOutput("plot.delta.Jz", width = "224px", height = "256px"),
          plotOutput("plot.num.scatters", width = "224px", height = "256px")
        )
        # plotOutput("plot.delta.Jz.and.num.scatters")
        # plotOutput("plot.delta.Jz.v.num.scatters")
      )
    ),

    # About tab: just display README.md
    tabPanel("About",
      includeMarkdown("README.md")
    )
  )
))

# --------------------------------------------------------------
