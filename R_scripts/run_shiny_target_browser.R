################################################################################
# run_shiny_target_browser.R
#
# Convenience launcher for local shell use.
#
# Usage:
#   Rscript R_scripts/run_shiny_target_browser.R
################################################################################

source("R_scripts/shiny_target_browser.R")

shiny::runApp(
  app,
  host = Sys.getenv("SHINY_HOST", unset = "127.0.0.1"),
  port = as.integer(Sys.getenv("SHINY_PORT", unset = "3838")),
  launch.browser = FALSE
)
