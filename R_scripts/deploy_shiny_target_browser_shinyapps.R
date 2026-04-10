################################################################################
# deploy_shiny_target_browser_shinyapps.R
#
# Deploy the 3'-tRF target explorer to shinyapps.io using a minimal bundle that
# includes only a deployment wrapper, the Shiny app script, and a clean DuckDB
# database file.
#
# One-time setup:
#   install.packages("rsconnect")
#   rsconnect::setAccountInfo(name = "<ACCOUNT>", token = "<TOKEN>", secret = "<SECRET>")
#
# Interactive use:
#   source("R_scripts/deploy_shiny_target_browser_shinyapps.R")
#   deploy_trf_target_browser_shinyapps()
#
# Non-interactive use:
#   Rscript R_scripts/deploy_shiny_target_browser_shinyapps.R
################################################################################

deploy_trf_target_browser_shinyapps <- function(
  app_name = "trf-target-browser-poc",
  app_title = "3'-tRF target explorer",
  app_visibility = "public",
  launch_browser = interactive()
) {
  if (!requireNamespace("rsconnect", quietly = TRUE)) {
    stop("The 'rsconnect' package is required. Install it with install.packages('rsconnect').")
  }

  app_visibility <- tolower(app_visibility)
  if (!app_visibility %in% c("public", "private")) {
    stop("app_visibility must be either 'public' or 'private'.")
  }

  project_dir <- normalizePath(".", mustWork = TRUE)
  app_source <- file.path(project_dir, "R_scripts", "shiny_target_browser.R")
  db_source <- file.path(project_dir, "import", "miranda", "target_browser.duckdb")
  wal_source <- file.path(project_dir, "import", "miranda", "target_browser.duckdb.wal")

  required_files <- c(app_source, db_source)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0L) {
    stop(
      "Missing required deployment files:\n",
      paste(sprintf("- %s", missing_files), collapse = "\n")
    )
  }

  if (file.exists(wal_source) && file.info(wal_source)$size > 0L) {
    stop(
      "The local DuckDB database still has an active WAL file:\n",
      wal_source, "\n\n",
      "Please stop the local Shiny app (or any other process using the database),\n",
      "wait for the WAL to clear, and then redeploy. This avoids shipping an\n",
      "incomplete live database snapshot to shinyapps.io."
    )
  }

  deploy_dir <- file.path(tempdir(), "trf_target_browser_shinyapps_bundle")
  if (dir.exists(deploy_dir)) {
    unlink(deploy_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(deploy_dir, recursive = TRUE, showWarnings = FALSE)

  app_wrapper <- c(
    "Sys.setenv(",
    "  TRF_TARGET_BROWSER_DB = 'target_browser.duckdb',",
    "  TRF_TARGET_BROWSER_CSV = '',",
    "  TRF_TARGET_BROWSER_READ_ONLY = 'true'",
    ")",
    "source('shiny_target_browser.R', local = TRUE)",
    "app"
  )
  writeLines(app_wrapper, file.path(deploy_dir, "app.R"))

  file.copy(app_source, file.path(deploy_dir, "shiny_target_browser.R"), overwrite = TRUE)
  file.copy(db_source, file.path(deploy_dir, "target_browser.duckdb"), overwrite = TRUE)

  tryCatch(
    rsconnect::deployApp(
      appDir = deploy_dir,
      appFiles = c("app.R", "shiny_target_browser.R", "target_browser.duckdb"),
      appPrimaryDoc = "app.R",
      appMode = "shiny",
      appName = app_name,
      appTitle = app_title,
      appVisibility = app_visibility,
      launch.browser = launch_browser
    ),
    error = function(err) {
      if (
        identical(app_visibility, "private") &&
        grepl("application\\.visibility|Feature is not enabled", conditionMessage(err))
      ) {
        stop(
          paste(
            "This shinyapps.io account does not support private app visibility.",
            "Deploy with app_visibility = 'public' for a free-tier proof of concept,",
            "or upgrade the shinyapps.io plan before requesting a private deployment."
          ),
          call. = FALSE
        )
      }

      stop(err)
    }
  )
}

show_trf_target_browser_logs <- function(
  app_name = "trf-target-browser-poc",
  account = NULL,
  entries = 200
) {
  if (!requireNamespace("rsconnect", quietly = TRUE)) {
    stop("The 'rsconnect' package is required. Install it with install.packages('rsconnect').")
  }

  rsconnect::showLogs(
    appName = app_name,
    account = account,
    entries = entries,
    streaming = FALSE
  )
}

if (!interactive()) {
  deploy_trf_target_browser_shinyapps(
    app_name = Sys.getenv("TRF_TARGET_BROWSER_APP_NAME", unset = "trf-target-browser-poc"),
    app_title = Sys.getenv("TRF_TARGET_BROWSER_APP_TITLE", unset = "3'-tRF target explorer"),
    app_visibility = Sys.getenv("TRF_TARGET_BROWSER_APP_VISIBILITY", unset = "public"),
    launch_browser = FALSE
  )
}
