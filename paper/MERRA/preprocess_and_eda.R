#!/usr/bin/env Rscript

source("MERRA/settings.R")
suppressPackageStartupMessages({
  library(ncdf4)
  library(sf)
})

data_dir <- merra_data_root
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# Read the 10 selected aerosol tracers at four native MERRA-2 layers.
nc <- nc_open(merra_source_file)
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
lev <- ncvar_get(nc, "lev")
time <- ncvar_get(nc, "time")
ii <- which(lon >= merra_lon_limit[1L] & lon <= merra_lon_limit[2L])
jj <- which(lat >= merra_lat_limit[1L] & lat <= merra_lat_limit[2L])
kk <- match(merra_layers, lev)
stopifnot(!anyNA(kk), merra_time_index >= 1L, merra_time_index <= length(time))

coords <- expand.grid(lon = lon[ii], lat = lat[jj], KEEP.OUT.ATTRS = FALSE)
read_field <- function(variable) {
  dim_names <- vapply(nc$var[[variable]]$dim, function(z) z$name, character(1))
  start <- count <- rep(1L, length(dim_names))
  selection <- list(
    lon = ii,
    lat = jj,
    lev = seq.int(min(kk), max(kk)),
    time = merra_time_index
  )
  for (name in names(selection)) {
    d <- match(name, dim_names)
    start[d] <- min(selection[[name]])
    count[d] <- max(selection[[name]]) - min(selection[[name]]) + 1L
  }
  value <- ncvar_get(nc, variable, start = start, count = count,
                     collapse_degen = FALSE)
  value <- aperm(value, match(c("lon", "lat", "lev", "time"), dim_names))
  value[, , kk - min(kk) + 1L, 1L, drop = FALSE]
}

p <- length(merra_variables) * length(merra_layers)
raw <- matrix(NA_real_, nrow(coords), p)
process_names <- character(p)
q <- 0L
for (variable in merra_variables) {
  value <- read_field(variable)
  for (layer_index in seq_along(merra_layers)) {
    q <- q + 1L
    raw[, q] <- as.vector(value[, , layer_index, 1L])
    process_names[q] <- sprintf(
      "%s_L%02d(~%.0fhPa)", variable, merra_layers[layer_index],
      merra_regional_mean_pressure_hpa[layer_index]
    )
  }
}
nc_close(nc)
colnames(raw) <- process_names
valid <- apply(raw, 1L, function(z) all(is.finite(z)))
raw <- raw[valid, , drop = FALSE]
coords <- coords[valid, , drop = FALSE]
process_indices <- seq_len(merra_p)
raw <- raw[, process_indices, drop = FALSE]
process_names <- process_names[process_indices]
p <- merra_p
if (any(raw < 0)) stop("negative aerosol mixing ratios detected")

# Match Argo-style preprocessing:
# raw -> offset/log10 -> full-grid local-plane removal -> residual SD.
offset <- vapply(seq_len(p), function(j) {
  positive <- raw[, j][raw[, j] > 0]
  if (!length(positive)) stop("process has no positive observations: ", process_names[j])
  0.5 * min(positive)
}, numeric(1))
log_value <- log10(sweep(raw, 2L, offset, "+"))

points <- st_as_sf(coords, coords = c("lon", "lat"), crs = 4326,
                   remove = FALSE)
xy_km <- st_coordinates(st_transform(points, 5070)) / 1000
local_mean_log <- matrix(NA_real_, nrow(log_value), p)
neighbor_count <- integer(nrow(log_value))
for (i in seq_len(nrow(log_value))) {
  delta <- sweep(xy_km, 2L, xy_km[i, ], "-")
  use <- which(rowSums(delta^2) < merra_mean_bandwidth_km^2)
  if (length(use) < 5L) stop("fewer than five local neighbors at row ", i)
  design <- cbind(1, delta[use, 1L], delta[use, 2L])
  fit <- lm.fit(design, log_value[use, , drop = FALSE])
  if (fit$rank < 3L) stop("rank-deficient local plane at row ", i)
  local_mean_log[i, ] <- fit$coefficients[1L, ]
  neighbor_count[i] <- length(use)
}
colnames(local_mean_log) <- process_names
residual_log <- log_value - local_mean_log
residual_mean <- colMeans(residual_log)
scale <- apply(residual_log, 2L, sd)
if (any(!is.finite(scale) | scale <= 0)) stop("invalid residual scale")
normalized <- sweep(residual_log, 2L, scale, "/")
colnames(normalized) <- process_names

sample_correlation <- cor(normalized)
sample_covariance <- cov(normalized)
eigen_cov <- eigen(sample_covariance, symmetric = TRUE)
tolerance <- max(eigen_cov$values) * 1e-8
retained <- eigen_cov$values > tolerance
sample_precision <- eigen_cov$vectors[, retained, drop = FALSE] %*%
  (t(eigen_cov$vectors[, retained, drop = FALSE]) /
     eigen_cov$values[retained])
rownames(sample_correlation) <- colnames(sample_correlation) <- process_names
rownames(sample_precision) <- colnames(sample_precision) <- process_names

saveRDS(list(
  coords = coords, xy_km = xy_km, raw = raw, offset = offset,
  log = log_value, local_mean_log = local_mean_log,
  residual_log = residual_log, residual_mean = residual_mean,
  center_applied = FALSE, scale = scale,
  normalized = normalized, process_names = process_names,
  process_indices = process_indices,
  variables = merra_variables, layers = merra_layers,
  target_pressure_hpa = merra_target_pressure_hpa,
  regional_mean_pressure_hpa = merra_regional_mean_pressure_hpa,
  bandwidth_km = merra_mean_bandwidth_km, neighbor_count = neighbor_count,
  sample_correlation = sample_correlation,
  sample_precision = sample_precision,
  precision_rank = sum(retained), covariance_eigenvalues = eigen_cov$values,
  source_md5 = unname(tools::md5sum(merra_source_file))
), file.path(data_dir, sprintf("preprocessed_P%d.rds", p)))

cat(sprintf("Created MERRA-2 P=%d preprocessing outputs\n", p))
cat("locations:", nrow(raw), " processes:", p, " precision rank:",
    sum(retained), "\n")
cat("neighbor count min/median/max:", min(neighbor_count),
    median(neighbor_count), max(neighbor_count), "\n")
