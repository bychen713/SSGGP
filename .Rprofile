local({
  # Activate the project's renv library from any working directory. The sandbox is
  # off so that packages provided by an HPC R module (e.g. sf, ncdf4) stay visible.
  Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = Sys.getenv("RENV_CONFIG_SANDBOX_ENABLED", "FALSE"))
  proj <- Sys.getenv("RENV_PROJECT", unset = getwd())
  activate <- file.path(proj, "renv", "activate.R")
  if (file.exists(activate)) source(activate)
})
