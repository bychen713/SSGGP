#!/usr/bin/env Rscript
# Build the cross-validation folds used by every method, from the output of
# preprocess_and_eda.R. Run from paper/:  Rscript MERRA/prepare_folds.R
#
#   1. uncentered source   data/preprocessed_P40_uncentered.rds
#   2. misaligned splits   data/baseline/misaligned/fold[1-3].rds
#   3. colocated splits    data/baseline/colocated/fold[1-3].rds
#   4. final folds         data/{misaligned,colocated}/fold[1-3].rds
#
# Steps 2-3 fix the train/test locations and the misaligned mask; step 4
# re-standardizes each process by the SD of its observed training residuals in
# that fold. The data/baseline/ folds are intermediate.
source("MERRA/settings.R")
baseline_root <- file.path(merra_data_root, "baseline")
final_root <- merra_data_root

# 1. uncentered source
# Preserve local-plane residuals without subtracting a global process mean.
local({
  old_file <- merra_preprocessed
  out_file <- merra_source
  stopifnot(file.exists(old_file), !file.exists(out_file))
  x <- readRDS(old_file)
  stopifnot(identical(dim(x$residual_log), c(2499L, merra_p)),
            all(is.finite(x$residual_log)),
            all(is.finite(x$scale)), all(x$scale > 0))

  old_cov <- cov(x$normalized)
  x$residual_mean <- colMeans(x$residual_log)
  x$center <- setNames(rep(0, ncol(x$residual_log)), colnames(x$residual_log))
  x$center_applied <- FALSE
  x$normalized <- sweep(x$residual_log, 2L, x$scale, "/")
  x$preprocessing <- "offset/log10 -> 1000-km local-plane removal -> divide by residual SD; no global centering"
  stopifnot(isTRUE(all.equal(unname(cov(x$normalized)), unname(old_cov),
                             tolerance = 1e-10)))
  x$sample_correlation <- cor(x$normalized)
  x$sample_precision <- solve(cov(x$normalized))

  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, out_file)
  cat("Saved", out_file, "\n")
  cat("Standardized residual mean range:", range(colMeans(x$normalized)), "\n")
})

# 2. misaligned splits
local({
  x <- readRDS(merra_source)
  stopifnot(identical(dim(x$normalized), c(2499L, merra_p)),
            identical(dim(x$xy_km), c(2499L, 2L)))
  data_dir <- file.path(baseline_root, "misaligned")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  targets <- file.path(data_dir, sprintf("fold%d.rds", seq_len(merra_n_folds)))
  if (any(file.exists(targets))) stop("refusing to overwrite prepared MERRA folds")
  xy0 <- sweep(x$xy_km, 2L, colMeans(x$xy_km))
  grid <- xy0 / max(diff(range(xy0[, 1L])), diff(range(xy0[, 2L]))) * merra_grid_scale
  colnames(grid) <- c("x", "y")
  # Allocate the test locations once so that evaluation sets are mutually
  # disjoint.  All training samples come from the remaining locations; training
  # sets may overlap each other, but no test location appears in any training set.
  set.seed(merra_seed)
  test_pool <- sample.int(nrow(x$normalized), merra_n_folds * merra_n_test)
  test_sets <- split(test_pool, rep(seq_len(merra_n_folds), each = merra_n_test))
  train_pool <- setdiff(seq_len(nrow(x$normalized)), test_pool)

  for (fold in seq_len(merra_n_folds)) {
    set.seed(merra_seed + 1000L * fold)
    tr <- sample(train_pool, merra_n_train, replace = FALSE)
    te <- test_sets[[fold]]
    set.seed(merra_seed + fold)
    masked <- matrix(FALSE, merra_n_train, merra_p,
                     dimnames = list(NULL, x$process_names))
    for (j in seq_len(merra_p))
      masked[sample.int(merra_n_train, merra_misaligned_prop * merra_n_train), j] <- TRUE
    Y <- x$normalized[tr, , drop = FALSE]
    Y[masked] <- NA_real_
    ans <- list(Y_mat = Y, Y_mat_test = x$normalized[te, , drop = FALSE],
      grid = grid[tr, , drop = FALSE], grid_test = grid[te, , drop = FALSE],
      process_names = x$process_names, masked_train = masked,
      missing_prop = merra_misaligned_prop, train_index = tr, test_index = te,
      lonlat_train = x$coords[tr, ], lonlat_test = x$coords[te, ],
      grid_scale = merra_grid_scale, source_md5 = unname(tools::md5sum(merra_source)),
      seed = merra_seed, mask_seed = merra_seed + fold,
      preprocessing = sprintf("offset/log10 -> full-grid %g-km local-plane removal -> divide by residual SD without global centering",
                              x$bandwidth_km))
    saveRDS(ans, targets[fold])
    span <- max(apply(grid[c(tr, te), , drop = FALSE], 2L, function(z) diff(range(z))))
    cat(sprintf("fold %d: %d x %d train, %d test; span %.1f\n",
                fold, merra_n_train, merra_p, merra_n_test, span))
  }
})

# 3. colocated splits
# Same train/test locations as the misaligned splits, with nothing masked.
local({
  full <- readRDS(merra_source)
  source_dir <- file.path(baseline_root, "misaligned")
  target_dir <- file.path(baseline_root, "colocated")
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  targets <- file.path(target_dir, sprintf("fold%d.rds", 1:3))
  if (any(file.exists(targets))) stop("refusing to overwrite colocated folds")
  for (fold in 1:3) {
    old <- readRDS(file.path(source_dir, sprintf("fold%d.rds", fold)))
    old$Y_mat <- full$normalized[old$train_index, , drop = FALSE]
    old$Y_mat_test <- full$normalized[old$test_index, , drop = FALSE]
    old$masked_train <- matrix(FALSE, nrow(old$Y_mat), ncol(old$Y_mat),
                               dimnames = dimnames(old$Y_mat))
    old$missing_prop <- 0
    old$design <- "colocated"
    old$source_misaligned_file <- file.path(source_dir, sprintf("fold%d.rds", fold))
    old$source_misaligned_md5 <- unname(tools::md5sum(old$source_misaligned_file))
    stopifnot(!anyNA(old$Y_mat), identical(dim(old$Y_mat), c(merra_n_train, merra_p)),
              identical(dim(old$Y_mat_test), c(merra_n_test, merra_p)))
    saveRDS(old, targets[fold])
  }
  z <- lapply(targets, readRDS)
  tests <- lapply(z, `[[`, "test_index")
  stopifnot(length(unique(unlist(tests))) == 600L,
            !any(unlist(tests) %in% unlist(lapply(z, `[[`, "train_index"))))
  cat("Prepared three colocated folds with 500 complete training locations each\n")
})

# 4. final folds: standardize by the observed training SD of each fold
# Reuse the uncentered residuals and the baseline train/test splits. Only the
# process-wise SD changes: estimate it from values actually observed in each
# training fold, then apply that same SD to its training and test data.
local({
  source_file <- merra_source
  full <- readRDS(source_file)
  stopifnot(identical(dim(full$residual_log), c(2499L, merra_p)),
            isFALSE(full$center_applied), all(is.finite(full$residual_log)),
            all(is.finite(full$scale)), all(full$scale > 0))

  targets <- unlist(lapply(c("misaligned", "colocated"), function(design)
    file.path(final_root, design, sprintf("fold%d.rds", 1:3))))
  if (any(file.exists(targets)))
    stop("refusing to overwrite existing fold-SD outputs")

  for (design in c("misaligned", "colocated")) {
    out_dir <- file.path(final_root, design)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    for (fold in 1:3) {
      old_file <- file.path(baseline_root, design, sprintf("fold%d.rds", fold))
      old <- readRDS(old_file)
      stopifnot(identical(dim(old$Y_mat), c(merra_n_train, merra_p)),
                identical(dim(old$Y_mat_test), c(merra_n_test, merra_p)),
                length(old$train_index) == 500L,
                length(old$test_index) == 200L,
                identical(dim(old$masked_train), c(merra_n_train, merra_p)))
      r_train <- full$residual_log[old$train_index, , drop = FALSE]
      r_test <- full$residual_log[old$test_index, , drop = FALSE]
      old_expected <- sweep(r_train, 2L, full$scale, "/")
      old_expected[old$masked_train] <- NA_real_
      stopifnot(isTRUE(all.equal(unname(old$Y_mat), unname(old_expected),
                                 tolerance = 1e-12)),
                isTRUE(all.equal(unname(old$Y_mat_test),
                                 unname(sweep(r_test, 2L, full$scale, "/")),
                                 tolerance = 1e-12)))
      if (design == "colocated") stopifnot(!any(old$masked_train))
      if (design == "misaligned")
        stopifnot(all(colSums(!old$masked_train) == 300L))

      observed <- r_train
      observed[old$masked_train] <- NA_real_
      fold_sd <- apply(observed, 2L, sd, na.rm = TRUE)
      stopifnot(length(fold_sd) == merra_p,
                all(is.finite(fold_sd)), all(fold_sd > 0))
      y_train <- sweep(r_train, 2L, fold_sd, "/")
      y_train[old$masked_train] <- NA_real_
      y_test <- sweep(r_test, 2L, fold_sd, "/")
      stopifnot(all(abs(apply(y_train, 2L, sd, na.rm = TRUE) - 1) < 1e-10),
                all(is.finite(y_test)))

      new <- old
      new$Y_mat <- y_train
      new$Y_mat_test <- y_test
      new$training_observed_sd <- fold_sd
      new$full_grid_sd_not_used <- full$scale
      new$source_fold_file <- old_file
      new$source_fold_md5 <- unname(tools::md5sum(old_file))
      new$source_md5 <- unname(tools::md5sum(source_file))
      new$preprocessing <- paste(
        "process-specific offset/log10 -> full-grid 1000-km local-plane removal",
        "-> divide by process SD of observed training residuals in this fold;",
        "no global centering")
      if (design == "colocated") {
        new$source_misaligned_file <- file.path(final_root, "misaligned",
                                                 sprintf("fold%d.rds", fold))
        new$source_misaligned_md5 <- unname(tools::md5sum(new$source_misaligned_file))
      }
      saveRDS(new, file.path(out_dir, sprintf("fold%d.rds", fold)))
      cat(design, "fold", fold, "observed/process",
          unique(colSums(!old$masked_train)), "SD range", range(fold_sd), "\n")
    }
  }
})
