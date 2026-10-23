#!/usr/bin/env Rscript
# Build the cross-validation folds used by every method, from the output of
# preprocess_and_eda.R. Run from paper/:  Rscript MERRA/prepare_folds.R
#
#   1. uncentered source   data/preprocessed_P40_uncentered.rds
#   2. misaligned splits   data/baseline/misaligned/fold[1-3].rds
#   3. colocated splits    data/baseline/colocated/fold[1-3].rds
#   4. final folds         data/{misaligned,colocated}/fold[1-3].rds
#
# Steps 2-3 fix the train/test locations and the misaligned mask. Step 4 builds
# each fold from the raw values: offset/log10, 1000-km local-plane trend, and
# division by the residual SD of each process, all estimated from the observed
# training values of that fold. The data/baseline/ folds are intermediate.
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

# 4. final folds from the raw values
# Within each fold: offset/log10, a local-plane trend within merra_mean_bandwidth_km
# of each location, and the residual SD of each process, all fitted to the
# observed training values and then applied to the training and test locations.
fit_local_plane <- function(target_xy, train_xy, train_log, observed, radius_km) {
  fitted <- matrix(NA_real_, nrow(target_xy), ncol(train_log))
  counts <- matrix(0L, nrow(target_xy), ncol(train_log))
  for (i in seq_len(nrow(target_xy))) {
    delta <- sweep(train_xy, 2L, target_xy[i, ], "-")
    nearby <- rowSums(delta^2) < radius_km^2
    for (j in seq_len(ncol(train_log))) {
      use <- nearby & observed[, j]
      counts[i, j] <- sum(use)
      if (counts[i, j] < 5L)
        stop("too few observed training neighbors at target ", i, ", process ", j)
      fit <- lm.fit(cbind(1, delta[use, 1L], delta[use, 2L]), train_log[use, j])
      if (fit$rank != 3L || !is.finite(fit$coefficients[1L]))
        stop("invalid local plane at target ", i, ", process ", j)
      fitted[i, j] <- fit$coefficients[1L]
    }
  }
  list(fitted = fitted, counts = counts)
}

local({
  source_file <- merra_source
  full <- readRDS(source_file)
  stopifnot(identical(dim(full$raw), c(2499L, merra_p)),
            identical(dim(full$xy_km), c(2499L, 2L)))

  targets <- unlist(lapply(c("misaligned", "colocated"), function(design)
    file.path(final_root, design, sprintf("fold%d.rds", 1:3))))
  if (any(file.exists(targets)))
    stop("refusing to overwrite existing fold outputs")

  for (design in c("misaligned", "colocated")) {
    out_dir <- file.path(final_root, design)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    for (fold in 1:3) {
      old_file <- file.path(baseline_root, design, sprintf("fold%d.rds", fold))
      old <- readRDS(old_file)
      observed <- !old$masked_train
      train_raw <- full$raw[old$train_index, , drop = FALSE]
      test_raw <- full$raw[old$test_index, , drop = FALSE]
      stopifnot(identical(dim(train_raw), c(merra_n_train, merra_p)),
                identical(dim(test_raw), c(merra_n_test, merra_p)),
                identical(dim(observed), c(merra_n_train, merra_p)),
                !any(old$train_index %in% old$test_index))
      if (design == "colocated") stopifnot(all(observed))

      offset <- vapply(seq_len(merra_p), function(j) {
        positive <- train_raw[observed[, j], j]
        positive <- positive[positive > 0]
        if (!length(positive)) stop("no positive training value for process ", j)
        min(positive) / 2
      }, numeric(1))
      train_log <- log10(sweep(train_raw, 2L, offset, "+"))
      test_log <- log10(sweep(test_raw, 2L, offset, "+"))
      train_xy <- full$xy_km[old$train_index, , drop = FALSE]
      test_xy <- full$xy_km[old$test_index, , drop = FALSE]
      train_trend <- fit_local_plane(train_xy, train_xy, train_log, observed,
                                     merra_mean_bandwidth_km)
      test_trend <- fit_local_plane(test_xy, train_xy, train_log, observed,
                                    merra_mean_bandwidth_km)
      train_resid <- train_log - train_trend$fitted
      test_resid <- test_log - test_trend$fitted
      observed_resid <- train_resid
      observed_resid[!observed] <- NA_real_
      fold_sd <- apply(observed_resid, 2L, sd, na.rm = TRUE)
      stopifnot(all(is.finite(fold_sd)), all(fold_sd > 0))
      y_train <- sweep(train_resid, 2L, fold_sd, "/")
      y_train[!observed] <- NA_real_
      y_test <- sweep(test_resid, 2L, fold_sd, "/")
      stopifnot(all(is.finite(y_train[observed])), all(is.finite(y_test)),
                max(abs(apply(y_train, 2L, sd, na.rm = TRUE) - 1)) < 1e-10)

      new <- old
      new$Y_mat <- y_train
      new$Y_mat_test <- y_test
      new$training_offset <- offset
      new$training_observed_sd <- fold_sd
      new$training_local_trend <- train_trend$fitted
      new$test_local_trend <- test_trend$fitted
      new$train_neighbor_count <- train_trend$counts
      new$test_neighbor_count <- test_trend$counts
      new$source_fold_file <- old_file
      new$source_fold_md5 <- unname(tools::md5sum(old_file))
      new$source_md5 <- unname(tools::md5sum(source_file))
      new$preprocessing <- paste(
        "offset/log10, 1000-km local-plane trend and process SD, each estimated",
        "from the observed training values of this fold; no global centering")
      saveRDS(new, file.path(out_dir, sprintf("fold%d.rds", fold)))
      cat(design, "fold", fold, "min neighbors",
          min(c(train_trend$counts, test_trend$counts)), "SD range",
          range(fold_sd), "\n")
    }
  }
})
