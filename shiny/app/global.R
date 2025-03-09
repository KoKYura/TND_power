library(dplyr)
library(tidyverse)
library(shiny)
library(shinydashboard)
library(DT)
library(JuliaCall)

julia_setup()
julia_command(paste0("include(\"", getwd(), "/shiny/app/TND_power.jl\")"))







