# Same report as sim_result.R; kept as the entry point for the all-Matern design.
candidates <- c("sim/result/sim_result.R", "result/sim_result.R", "sim_result.R")
report <- candidates[file.exists(candidates)][1L]
if (is.na(report)) stop("run from paper/, paper/sim/, or paper/sim/result/")
source(report)


