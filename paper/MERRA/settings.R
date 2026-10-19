# MERRA-2 configuration, sourced by every script in MERRA/. Run from paper/.

# raw data and preprocessing (preprocess_and_eda.R)
# Native-coordinate P = 40 configuration. Aerosol values remain on native model
# layers; no vertical interpolation.
merra_source_file <- Sys.getenv("MERRA_SOURCE_FILE", file.path(
  "MERRA", "raw", "MERRA2_400.inst3_3d_aer_Nv.20220715.nc4"
))

merra_lon_limit <- c(-95, -65)
merra_lat_limit <- c(25, 50)
merra_time_index <- 1L                 # 00 UTC in R's one-based indexing
merra_mean_bandwidth_km <- 1000

merra_variables <- c(
  "DU001", "DU003", "DU005",
  "SS001", "SS005",
  "BCPHOBIC", "BCPHILIC",
  "OCPHOBIC", "OCPHILIC",
  "SO4"
)

# Selected by closest regional-mean native mid-layer pressure at 00 UTC.
# These labels describe approximate regional means only; pressure varies by cell.
merra_layers <- c(50L, 52L, 55L, 62L)
merra_target_pressure_hpa <- c(500, 600, 700, 850)
merra_regional_mean_pressure_hpa <- c(506.461, 581.530, 687.881, 843.032)

stopifnot(length(merra_variables) * length(merra_layers) == 40L)

# Select the first P processes in variable-major, layer-minor order.
merra_p <- as.integer(Sys.getenv("MERRA_P", "40"))
stopifnot(length(merra_p) == 1L, !is.na(merra_p), merra_p >= 2L, merra_p <= 40L)
merra_data_root <- Sys.getenv("MERRA_DATA_ROOT", file.path("MERRA", "data"))
merra_output_base <- Sys.getenv("MERRA_OUTPUT_ROOT", file.path("MERRA", "output"))
merra_preprocessed <- file.path(merra_data_root, sprintf("preprocessed_P%d.rds", merra_p))
merra_source <- file.path(merra_data_root, sprintf("preprocessed_P%d_uncentered.rds", merra_p))

# cross-validation design (prepare_folds.R)
merra_seed <- 20260929L
merra_n_folds <- 3L
merra_n_train <- 500L
merra_n_test <- 200L
merra_misaligned_prop <- 0.40          # share of training values masked per process
merra_grid_scale <- 10

# which design a fitting script runs on
# misaligned: 40% of training values masked per process; colocated: none.
merra_design <- Sys.getenv("MERRA_DESIGN", "misaligned")
stopifnot(merra_design %in% c("misaligned", "colocated"))
merra_missing_prop <- if (merra_design == "misaligned") merra_misaligned_prop else 0
merra_data_dir <- Sys.getenv("MERRA_DATA_DIR", file.path(merra_data_root, merra_design))
merra_output_root <- file.path(merra_output_base, merra_design)
merra_output_dir <- Sys.getenv("MERRA_OUTPUT_DIR", file.path(merra_output_root, "ssggp"))

# SSGGP (main/)
merra_L <- as.integer(Sys.getenv("MERRA_L", "20000"))
merra_MC <- as.integer(Sys.getenv("MERRA_MC", "40"))
merra_burnin <- merra_MC %/% 2L
merra_lag <- as.integer(Sys.getenv("MERRA_LAG", "30"))
merra_iter <- as.integer(Sys.getenv("MERRA_ITER", "700"))
merra_tol <- 1e-3
merra_rho <- c(0, 0.0001, 0.00025, 0.0004, 0.00055, 0.0007,
               0.00085, 0.001, 0.00125, 0.00175, 0.0025, 0.003,
               0.005, 0.01, 0.02)
merra_L_ll <- 1000L
merra_L_pred <- as.integer(Sys.getenv("MERRA_L_PRED", "30000"))
merra_MC_pred <- as.integer(Sys.getenv("MERRA_MC_PRED", "20000"))
merra_pred_suffix <- Sys.getenv("MERRA_PRED_SUFFIX", "")
merra_pred_single_method <- identical(Sys.getenv("MERRA_PRED_SINGLE_METHOD", "0"), "1")
merra_threads <- as.integer(Sys.getenv("MERRA_THREADS", "16"))
stopifnot(merra_L >= 1L, merra_MC >= 2L, merra_lag >= 1L, merra_iter >= 1L,
          merra_rho[1L] == 0, all(diff(merra_rho) > 0),
          merra_L_pred >= 1L, merra_MC_pred >= 2L,
          grepl("^(_[A-Za-z0-9]+)?$", merra_pred_suffix))
merra_fold_file <- function(fold) file.path(merra_data_dir, sprintf("fold%d.rds", fold))
merra_ind_file <- function(fold) file.path(merra_output_dir, sprintf("ind_%d.rds", fold))
merra_est_file <- function(fold, ri) file.path(
  merra_output_dir, sprintf("est_fold%d_rho%d.rds", fold, ri))
