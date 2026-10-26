#!/bin/bash
# Set up the project's R environment (renv) and install SSGGP from pkg/.
#   bash setup_env.sh               # packages for the demo, simulations, and fits
#   bash setup_env.sh --with-merra  # also the MERRA-2 preprocessing packages (ncdf4, sf,
#                                    # patchwork, maps; need GDAL/GEOS/PROJ/udunits/netCDF)
#   bash setup_env.sh --with-iox    # also the bundled spiox (paper/iox/spiox) for IOX
#   bash setup_env.sh --no-system   # do not install system libraries or load modules
set -uo pipefail
cd "$(dirname "$0")"

# Activate renv through RENV_PROJECT, including worker sessions.
if ! grep -q "RENV_PROJECT" .Rprofile 2>/dev/null; then
  cat > .Rprofile <<'RPROFILE'
# Activate the project library using RENV_PROJECT.
local({
  # Allow packages supplied by the R module; avoid sandbox startup overhead.
  Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = Sys.getenv("RENV_CONFIG_SANDBOX_ENABLED", "FALSE"))
  proj <- Sys.getenv("RENV_PROJECT", unset = getwd())
  activate <- file.path(proj, "renv", "activate.R")
  if (file.exists(activate)) source(activate)
})
RPROFILE
  echo "Wrote .Rprofile (renv activation)."
fi

SYS=1; WANT_MERRA=0; WANT_IOX=0
for a in "$@"; do case "$a" in
  --with-merra)     WANT_MERRA=1;;
  --with-iox)       WANT_IOX=1;;
  --no-system)      SYS=0;;
  *) echo "Unknown option: $a"; exit 2;; esac; done

# ---- R packages (dependencies are pulled in automatically from renv.lock) ---
ROOTS_CORE="fields glasso Matrix this.path Rcpp RcppEigen RhpcBLASctl"
ROOTS_PLOT="ggplot2 ggh4x gridExtra reshape2 forcats dplyr tidyr purrr stringr magrittr scales"
ROOTS="$ROOTS_CORE $ROOTS_PLOT"
VERIFY="Matrix Rcpp fields spam glasso ggplot2"
# MERRA-2 preprocessing packages are not pinned in renv.lock: sf's dependencies
# (s2, units) are hard to build on some clusters, so a copy provided by the R
# module is used when present, with CRAN as the fallback.
ROOTS_MERRA=""
[ $WANT_MERRA -eq 1 ] && ROOTS_MERRA="ncdf4 sf patchwork maps"

# ---- system prerequisites, grouped by purpose -----------------------------
if [ $SYS -eq 1 ]; then
  if command -v module >/dev/null 2>&1; then                          # ---- HPC (modules) ----
    # default: R (compilers are available by default on the node; R also stores
    #          absolute compiler paths, so sourceCpp works without a gcc module)
    # iox: build spiox with the gcc whose C++ runtime R loads (the IOX launchers
    #      load the same module), avoiding "GLIBCXX ... not found" at load time
    if [ $WANT_IOX -eq 1 ]; then
      module load gcc/8.5.0 2>/dev/null || echo "note: no 'gcc/8.5.0' module here"
    fi
    module load r/4.4.1 2>/dev/null || module load r 2>/dev/null \
      || echo "note: could not module-load R; ensure R 4.4.x is on PATH"
    # SSGGP requires GSL at build and load time.
    module load gsl/2.8 2>/dev/null || module load gsl 2>/dev/null \
      || echo "note: could not module-load GSL; SSGGP will not build without it"
    # plots: no extra system modules needed
    # merra: spatial + netCDF libraries (module names vary by cluster; non-fatal)
    if [ $WANT_MERRA -eq 1 ]; then
      for m in sqlite geos proj gdal udunits netcdf-c; do
        module load $m 2>/dev/null || echo "note: no '$m' module here -- fine if R already provides sf/ncdf4"
      done
    fi
  elif [ "$(uname -s)" = "Darwin" ]; then                             # ---- macOS (Homebrew) ----
    # compiler + OpenMP runtime (libomp fixes the spam load error) + GSL for SSGGP
    xcode-select -p >/dev/null 2>&1 || xcode-select --install 2>/dev/null
    brew install gcc libomp gsl 2>/dev/null
    [ $WANT_MERRA -eq 1 ] && brew install pkg-config gdal geos proj sqlite udunits netcdf 2>/dev/null
  elif command -v apt-get >/dev/null 2>&1; then                       # ---- Debian/Ubuntu ----
    P="build-essential gfortran libgsl-dev"
    [ $WANT_MERRA -eq 1 ] && P="$P pkg-config libgdal-dev libgeos-dev libproj-dev libsqlite3-dev libudunits2-dev libnetcdf-dev"
    sudo apt-get update -qq && sudo apt-get install -y $P || echo "Install manually: sudo apt-get install -y $P"
  elif command -v dnf >/dev/null 2>&1; then                           # ---- Fedora/RHEL ----
    P="gcc gcc-c++ gfortran gsl-devel"
    [ $WANT_MERRA -eq 1 ] && P="$P pkgconf-pkg-config gdal-devel geos-devel proj-devel sqlite-devel udunits2-devel netcdf-devel"
    sudo dnf install -y $P || echo "Install manually: sudo dnf install -y $P"
  elif uname -s | grep -qE "MINGW|MSYS|CYGWIN"; then                  # ---- Windows (Git Bash) ----
    echo "Windows: install R 4.4.x + Rtools44 (https://cran.r-project.org/bin/windows/Rtools/) and GSL before running -- SSGGP compiles C++ at install time."
  fi
fi

# ---- restore pinned packages + verify + repair spam + check compiler -------
# Install serially: renv's parallel (socket-based) installer crashes on some
# Linux setups (socketSelect error in renv_graph_install).
export RENV_CONFIG_INSTALL_JOBS=1
export RENV_CONFIG_SANDBOX_ENABLED=FALSE   # see .Rprofile
export ROOTS VERIFY ROOTS_MERRA WANT_IOX
Rscript -e '
# Activate renv for this session without renv::activate(), which would rewrite
# .Rprofile with a relative path that breaks parallel workers.
source("renv/activate.R")
options(renv.config.install.transactional = FALSE)   # keep successes if some fail
roots <- strsplit(Sys.getenv("ROOTS")," ")[[1]]; roots <- roots[nzchar(roots)]
if (!isTRUE(tryCatch({renv::restore(packages=roots, prompt=FALSE); TRUE}, error=function(e) FALSE)))
  for (p in roots) tryCatch(renv::restore(packages=p, prompt=FALSE),
                            error=function(e) message("  failed: ", p))

bad <- Filter(function(p) !isTRUE(tryCatch({suppressMessages(loadNamespace(p)); TRUE}, error=function(e) FALSE)),
              strsplit(Sys.getenv("VERIFY")," ")[[1]])
if ("spam" %in% bad) {                       # spam: rebuild against local toolchain (macOS libomp / missing gfortran)
  message("Rebuilding spam/dotCall64 from source...")
  tryCatch(renv::install(c("dotCall64","spam"), type="source", rebuild=TRUE, prompt=FALSE), error=function(e) NULL)
  if (isTRUE(tryCatch({suppressMessages(loadNamespace("spam")); TRUE}, error=function(e) FALSE))) bad <- setdiff(bad,"spam")
}
# Generic retry: reinstall anything still failing to load. Covers dependents
# stranded when their dependency was fixed/installed only later (e.g. fields
# after spam), or a transient build failure -- on any platform.
if (length(bad)) {
  tryCatch(renv::restore(packages=bad, prompt=FALSE), error=function(e) NULL)
  bad <- Filter(function(p) !isTRUE(tryCatch({suppressMessages(loadNamespace(p)); TRUE}, error=function(e) FALSE)), bad)
}

# ---- MERRA-2 preprocessing packages (--with-merra) --------------------------
merra <- strsplit(Sys.getenv("ROOTS_MERRA"), " ")[[1]]; merra <- merra[nzchar(merra)]
if (length(merra)) {
  have <- function(p) isTRUE(tryCatch({suppressMessages(loadNamespace(p)); TRUE}, error=function(e) FALSE))
  need <- Filter(Negate(have), merra)
  if (length(need)) {
    message("Installing MERRA packages from CRAN: ", paste(need, collapse=", "))
    tryCatch(renv::install(need, prompt=FALSE), error=function(e) message("  ", conditionMessage(e)))
  }
  bad <- c(bad, Filter(Negate(have), merra))
}

# ---- the SSGGP package itself (pkg/), into the renv library -----------------
# R CMD INSTALL with an explicit --library, rather than install.packages(), so
# renv cannot mistake the local path for a CRAN package name.
if (!length(bad)) {
  message("Installing SSGGP from pkg/ ...")
  # Zip downloads can drop the executable bits that R CMD INSTALL requires.
  Sys.chmod(c("pkg/configure", "pkg/cleanup"), "755")
  Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
  # A failed install keeps any previous copy, so check the exit status too.
  status <- system2(file.path(R.home("bin"), "R"),
                    c("CMD", "INSTALL", "--preclean", paste0("--library=", shQuote(.libPaths()[1])), "pkg"))
  if (status != 0L ||
      !isTRUE(tryCatch({suppressMessages(loadNamespace("SSGGP")); TRUE}, error=function(e) FALSE))) {
    message("SSGGP (pkg/) failed to install -- see the compiler output above. ",
            "The usual cause is GSL not found: check that gsl-config is on PATH.")
    bad <- c(bad, "SSGGP")
  }
}

# ---- the bundled, patched spiox (paper/iox/spiox), for the IOX baseline -----
if (!length(bad) && Sys.getenv("WANT_IOX") == "1") {
  message("Installing the bundled spiox from paper/iox/spiox ...")
  tryCatch(renv::restore(packages=c("FNN", "GpGp", "RcppArmadillo"), prompt=FALSE),
           error=function(e) message("  ", conditionMessage(e)))
  status <- system2(file.path(R.home("bin"), "R"),
                    c("CMD", "INSTALL", "--preclean", "--clean",
                      paste0("--library=", shQuote(.libPaths()[1])), "paper/iox/spiox"))
  if (status != 0L ||
      !isTRUE(tryCatch(utils::packageVersion("spiox") == "0.8.1.9001", error=function(e) FALSE))) {
    message("spiox (paper/iox/spiox) failed to install -- see the compiler output above.")
    bad <- c(bad, "spiox")
  }
}

if (length(bad)) { message("Not loading: ", paste(bad, collapse=", "), " (see README.MD)"); quit(status=1) }
message("OK. Next: bash run_demo.sh")
'
