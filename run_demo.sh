#!/bin/bash
# Run demo.R in the project environment. Run bash setup_env.sh first.
set -e
cd "$(dirname "$0")"

if [ ! -f renv.lock ] || [ ! -d renv/library ]; then
  echo "Project environment not found. Run: bash setup_env.sh first."
  exit 1
fi

# R and GSL from environment modules, where available
if command -v module >/dev/null 2>&1; then
  module load r/4.4.1 2>/dev/null || module load r 2>/dev/null || true
  module load gsl/2.8 2>/dev/null || module load gsl 2>/dev/null
fi

if ! Rscript -e 'quit(status = !requireNamespace("SSGGP", quietly = TRUE))' >/dev/null 2>&1; then
  echo "The SSGGP package is not installed in the project library."
  echo "Run: bash setup_env.sh   (it installs pkg/ as its last step)"
  exit 1
fi

Rscript demo.R
echo "Demo completed successfully."
