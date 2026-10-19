# Paper scripts run from paper/; use the repository-root .Rprofile and renv library.
local({
  root <- normalizePath("..", mustWork = TRUE)
  Sys.setenv(RENV_PROJECT = Sys.getenv("RENV_PROJECT", unset = root))
  profile <- file.path(root, ".Rprofile")
  if (file.exists(profile)) source(profile)
})
