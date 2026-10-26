# Simulation results for the active setting: prediction and graph-recovery
# tables, estimation error, computation time, precision/correlation figures,
# and BIC-selected SAECM diagnostic traces.
#   Rscript sim/result/sim_result.R
# Tables go to sim/output/tables/<setting>/ and figures to sim/output/figure/<setting>/.
# LMC columns are included when its results exist. For mixed families the
# family-specific estimation tables of sim_result_est.R are produced as well.

rm(list=ls())
library(Matrix)
if (!file.exists("settings.R")) {  # work from sim/, where settings.R lives
  if (file.exists(file.path("sim", "settings.R"))) setwd("sim") else
  if (file.exists(file.path("..", "settings.R"))) setwd("..") else
    stop("run this script from paper/")
}
source("settings.R")

family_names <- vapply(f_fun_true, function(x) {
  if (is.character(x) && length(x) == 1L) x else NA_character_
}, character(1))
run_family_result <- anyNA(family_names) || !all(family_names == "matern")
est_result_script <- here("sim", "result", "sim_result_est.R")

summarize_results <- function() {
  helper_names <- c("trimat", "get_result", "build_matrix", "compute_sparsity")
  sys.source(here("sim", "result", "helpers.R"), envir = environment())
  missing_helpers <- helper_names[!vapply(helper_names, exists, logical(1),
                                          envir=environment(), inherits=FALSE)]
  if (length(missing_helpers))
    stop("could not load result helpers: ", paste(missing_helpers, collapse=", "))

  all_results <- list()
  all_settings <- expand.grid(prec_structure = c('tridiagonal', 'random'),
                              offdiag_sign = c('neg', 'pos'),
                              missing_TF = c('coloc','miss'),
                              stringsAsFactors = FALSE
  )
  plot_set <- which(all_settings$prec_structure == prec_structure &
                    all_settings$offdiag_sign == offdiag_sign &
                    all_settings$missing_TF == missing_TF)
  if (length(plot_set) != 1L)
    stop("could not match the active setting in all_settings: ", setting_name)
  all_skip <- list()
  for(setting_ind in plot_set){
    
    prec_structure <- as.character( all_settings[setting_ind, 'prec_structure'])
    offdiag_sign <- as.character( all_settings[setting_ind, 'offdiag_sign'])
    missing_TF <- as.character( all_settings[setting_ind, 'missing_TF'])
    
    # The penalty grid and the true precision come from settings.R.
    param_grid <- unique(subset(grid_params, select = -rho))
    param_grid$prec_struc <- rep(prec_structure, n_sim)
    Sigma_inv <- true_precision()
    
    path <- out_dir
    setwd(path)
    
    n_rhos <- length(rho_vals)
    n_settings <- nrow(grid_params)/n_rhos/n_sim
    n_per_setting <- nrow(grid_params)/n_settings
    
    files_all <- list.files(path)
    files_all <- files_all[!file.info(list.files(path,full.names = T))$isdir]
    
    # independent process files
    files_ind_exist <- files_all[substr(files_all,1,3)=='ind']
    files_number_ind <- if (length(files_ind_exist))
      as.integer(sub("^ind_([0-9]+)\\.RData$", "\\1", files_ind_exist)) else
      integer(0)
    if (anyNA(files_number_ind))
      stop("cannot read a replicate index from: ",
           paste(files_ind_exist[is.na(files_number_ind)], collapse=", "))
    files_ind_exist <- files_ind_exist[order(files_number_ind)]
    files_number_ind <- files_number_ind[order(files_number_ind)]
    files_ind <- rep(NA,n_settings*n_sim)
    files_ind[files_number_ind] <- files_ind_exist
    if (!length(files_ind_exist))
      stop("no ind_*.RData files found in ", path,
           "; run sim/main/sim_ind.R for setting ", setting_name,
           " or use sim_result_est.R to summarise the available SAEM fits")
    
    # LMC files (lmcscov_*) are optional; without them the LMC columns are omitted.
    files_lmc_exist <- files_all[substr(files_all,1,4)=='lmcs']
    have_lmc <- length(files_lmc_exist) > 0
    if (have_lmc) {
      files_number_lmc <- as.integer(
        sapply(strsplit(sapply(strsplit(files_lmc_exist, '\\.'), function(x) x[[1]]), '_'), function(x) x[[2]])
      )
      files_lmc_exist <- files_lmc_exist[order(files_number_lmc)]
      files_number_lmc <- files_number_lmc[order(files_number_lmc)]
      files_lmc <- rep(NA,n_settings*n_sim)
      files_lmc[files_number_lmc] <- files_lmc_exist
    } else {
      message("no lmcs_*.RData in ", path, " -- LMC columns will be omitted")
      files_lmc <- rep(NA, n_settings*n_sim)
    }
    
    # SAEM files
    files_dense <- files_BIC <- rep(NA,n_settings*n_sim)
    files_dense_exist <- files_all[substr(files_all,1,3)=='den']
    files_BIC_exist <- files_all[substr(files_all,1,3)=='BIC']
    
    # Extract the replicate index from the filename's final field.
    index_from <- function(f) {
      j <- as.integer(sub("^.*_([0-9]+)\\.RData$", "\\1", f))
      if (anyNA(j)) stop("cannot read a replicate index from: ",
                         paste(f[is.na(j)], collapse = ", "))
      list(i = rep(1L, length(f)), j = j)
    }
    # dense
    ix <- index_from(files_dense_exist)
    ord <- order(ix$i, ix$j)
    files_dense_exist <- files_dense_exist[ord]
    files_number_i <- ix$i[ord]; files_number_j <- ix$j[ord]
    files_dense[n_sim*(files_number_i-1)+files_number_j] <- files_dense_exist
    # BIC
    ix <- index_from(files_BIC_exist)
    ord <- order(ix$i, ix$j)
    files_BIC_exist <- files_BIC_exist[ord]
    files_number_i <- ix$i[ord]; files_number_j <- ix$j[ord]
    files_BIC[n_sim*(files_number_i-1)+files_number_j] <- files_BIC_exist
    
    results_list <- BIC_list <- llmcem_list <- n_par_list <- list()
    for(i in 1:n_settings){
      print(i)
      BIC_all <- llMCEM_all <- n_par <- c()
      for(j in 1:n_sim){
        print(j)
        array_id <- (n_sim*(i-1)+j)
        # independent process
        load(file.path(path, files_ind[array_id]))
        
        n <- nrow(result_ind$grid)
        n_test <- nrow(result_ind$grid_test)
        n_all <- n + n_test
        result_ind <- get_result(result_ind,"ind")

        # LMC
        if (have_lmc) {
          load(file.path(path, files_lmc[array_id]))
          result_lmc <- get_result(result_lmc,"lmc")
        }

        # dense precision, rho=0
        load(file.path(path, files_dense[array_id]))
        result_dense <- get_result(result_final,"dense")
        rm(result_final)
        time_dense <- time_approx
        
        load(file.path(path, files_BIC[array_id]))
        BIC_all <- rbind(BIC_all, result_final$BICs)
        result_BIC <- get_result(result_final,"BIC")
        time_BIC <- time_approx
        rm(result_final)

        # Fitting times are stored as seconds or difftime objects; report minutes.
        as_mins <- function(x, numeric_units = "mins") {
          if (inherits(x, "difftime")) { units(x) <- "mins"; as.numeric(x) }
          else as.numeric(x) / if (numeric_units == "secs") 60 else 1
        }
        time_true  <- as_mins(time_true)
        time_BIC   <- as_mins(time_BIC, "secs")
        time_dense <- as_mins(time_dense, "secs")
        if (have_lmc) time_lmc <- as_mins(time_lmc)

        lmc_block <- if (have_lmc) c(unlist(result_lmc$data_plot), time_lmc) else NULL

        results_list <- c(results_list,
                          list(c(i,
                                 unlist(result_ind$data_plot),
                                 time_true,

                                 lmc_block,

                                 unlist(result_dense$data_plot),
                                 time_dense,

                                 unlist(result_BIC$data_plot),
                                 time_BIC
                          )
                          )
        )
      }
      BIC_list <- c(BIC_list, list(BIC_all))
    }
    
    library(tidyr)
    library(dplyr)
    library(ggplot2)
    library(forcats)
    library(stringr)
    library(purrr)
    
    param_grid <- unique(subset(grid_params, select = -rho))
    param_grid$row <- 1:nrow(param_grid)
    results_df <- cbind(param_grid %>% filter((row %in% (n_sim*(files_number_i-1)+files_number_j) )),
                        do.call(rbind, results_list))
    ind <- which(upper.tri(matrix(nrow=P, ncol=P),diag=T),arr.ind=T)
    ind <- ind[order(ind[,2]-ind[,1]),,drop=F]
    ind <- paste(ind[,1],ind[,2],sep="")
    sigma_inv_1_ind <- paste0(1:P,1:P)
    
    {
      grid_names <- c('sim', 'L','lag', 'MC', 'row', 'number'
      )
      
      ind_names <- c(paste0('prec_true', 1:P),paste0('a_true',1:P), 
                     paste0('nu_true',1:P),paste0('nugget_true', 1:P),
                     'norm_true', paste0('rmse_true',1:P),
                     paste0('coverage_true',1:P),
                     paste0('var_true',1:P), 
                     'rmsesum_true',"coveragesum_true", "varsum_true",
                     'time_true'
      )
      lmc_names <- if (!have_lmc) NULL else
               c(paste0('sigma_lmc', 1:(P*q)),paste0('a_lmc',1:q),
                     paste0('nu_lmc',1:q),paste0('nugget_lmc', 1:P),
                     'norm_lmc', paste0('rmse_lmc',1:P),
                     paste0('coverage_lmc',1:P),
                     paste0('var_lmc',1:P),
                     'rmsesum_lmc',"coveragesum_lmc", "varsum_lmc",
                     'time_lmc'
      )
      approx_names <- c('prec', 'nugget', 'a', 'nu',
                        'rho',
                        'll',
                        'norm',
                        'rmse',
                        'coverage',
                        'var',
                        'rmsesum',
                        'coveragesum', 
                        'varsum'
      )
      get_names <- function(approx_names,method){
        out <- c()
        for(s in approx_names){
          if(s=='prec'){
            out <- c(out, paste0(paste0(s,'_',method),ind))
          }else if(substr(s,1,4)=='norm' || s %in% c('nugget','ll','rho') || 
                   any(grepl('sum', s))){
            out <- c(out, paste0(s,'_',method))
          }else{
            out <- c(out, paste0(paste0(s,'_',method),1:P))
          }
        }
        out 
      }
      
      colnames(results_df) <- c(grid_names, ind_names, lmc_names,
                                get_names(approx_names, 'dense'), 'time_dense',
                                get_names(approx_names, 'BIC'),'time_BIC'
                                
      )
      # Convert squared prediction errors to RMSE once, before reporting.
      is_rmse <- grepl('^rmse', colnames(results_df))
      results_df[is_rmse] <- lapply(results_df[is_rmse], sqrt)
    }
    
    # Keep the columns used in the report.
    results_df <- results_df %>% dplyr::select(contains(c("lmc","BIC","dense","true","P","L","prec","a_","nu_","nugget_","time"),
                                                        ignore.case=F))
    
    results_df <- results_df %>%
      mutate(across(starts_with("var"), ~ sqrt(.x) * 1.96 * 2,  .names = "tmp_{.col}")) %>%
      rename_with(~ sub("^tmp_var", "len", .x), starts_with("tmp_var"))
    
    # Graph-recovery rates against the true precision.
    off_diag_mask <- row(Sigma_inv) != col(Sigma_inv)
    true_pos_mask <- (Sigma_inv != 0) & off_diag_mask
    true_neg_mask <- (Sigma_inv == 0) & off_diag_mask
    
    n_true_pos    <- sum(true_pos_mask)
    n_true_neg    <- sum(true_neg_mask)
    
    results_df <- results_df %>%
      rowwise() %>%
      mutate(
        md   = list(list(compute_sparsity(pick(everything()), P, "prec_dense", true_pos_mask, true_neg_mask, n_true_pos, n_true_neg))),
        ma   = list(list(compute_sparsity(pick(everything()), P, "prec_BIC",   true_pos_mask, true_neg_mask, n_true_pos, n_true_neg))),
        FDR_dense = md[[1]][1], FPR_dense = md[[1]][2], TPR_dense = md[[1]][3], F2_dense = md[[1]][4],
        FDR_BIC   = ma[[1]][1], FPR_BIC   = ma[[1]][2], TPR_BIC   = ma[[1]][3], F2_BIC = ma[[1]][4]
      ) %>%
      ungroup() %>%
      select(-md, -ma)
    
    
    results_df_ori <- results_df
    all_results[[setting_ind]] <- results_df
    
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    setwd(fig_dir)
    all_ind <- rep(1:n_settings, each = n_sim)
    for(j in 1:n_settings){
      print(j)
      i <- n_sim*(j-1)+1
      results_df <- results_df_ori[all_ind==j,]
      Sigma <- solve(Sigma_inv)
      
      library(reshape2)
      library(gridExtra)
      library(ggh4x)
      
      # Median estimated precision matrices
      quan <- 0.5
      precision_types <- c("prec_BIC", "prec_dense")
      quant_row <- results_df %>%
        dplyr::select(!contains('prec_struc')) %>%
        summarise(across(everything(), \(x) quantile(x, probs = quan, na.rm = TRUE)))
      
      lapply(precision_types, function(prefix) {
        mat <- build_matrix(quant_row, P, prefix)
        mat_melt <- melt(mat)
        colnames(mat_melt) <- c("Row", "Col", "Value")
        
        n_show <- 5
        axis_labels <- levels(factor(mat_melt$Col))
        thin_labels <- ifelse(seq_along(axis_labels) %% n_show == 1, axis_labels, "")
        
        p <- ggplot(mat_melt, aes(x = factor(Col), y = factor(Row), fill = Value)) +
          geom_raster() +
          scale_fill_gradient2(
            low      = "#2166ac",
            mid      = "white",
            high     = "#b2182b",
            midpoint = 0
          ) +
          scale_x_discrete(labels = thin_labels) +
          scale_y_discrete(labels = thin_labels) +
          labs(x = 'P', y = 'P', fill = 'Precision') +
          theme_bw() +
          theme(
            axis.title   = element_text(size = 14, face = "bold"),
            axis.text.x  = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
            axis.text.y  = element_text(size = 12),
            legend.title = element_text(size = 10),
            legend.text  = element_text(size = 9),
            panel.grid   = element_blank()
          ) +
          coord_fixed() +
          force_panelsizes(rows = unit(10, "cm"), cols = unit(10, "cm"))
        
        file_name <- paste0("sigmainv_", prefix, ".png")
        ggsave(file_name, p, height = 5, width = 5.5, device = "png")
        return(p)
      })
      
      # True correlation matrix
      corr <- diag(1/sqrt(diag(Sigma))) %*% Sigma %*% diag(1/sqrt(diag(Sigma)))
      mat_melt <- melt(corr)
      colnames(mat_melt) <- c("Row", "Col", "Value")
      threshold <- 0
      mat_melt$TextColor <- ifelse(mat_melt$Value > threshold, "black", "white")
      
      n_show <- 5
      axis_labels <- levels(factor(mat_melt$Col))
      thin_labels <- ifelse(seq_along(axis_labels) %% n_show == 1, axis_labels, "")
      
      ggplot(mat_melt, aes(x = factor(Col), y = factor(Row), fill = Value)) +
        geom_raster() +
        scale_fill_gradient2(
          low      = "#2166ac",
          mid      = "white",
          high     = "#b2182b",
          midpoint = 0
        ) +
        scale_color_identity() +
        scale_x_discrete(labels = thin_labels) +
        scale_y_discrete(labels = thin_labels) +
        labs(x = 'P', y = 'P', fill = 'Correlation') +
        theme_bw() +
        theme(
          axis.title   = element_text(size = 14, face = "bold"),
          axis.text.x  = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
          axis.text.y  = element_text(size = 12),
          legend.title = element_text(size = 10),
          legend.text  = element_text(size = 9),
          panel.grid   = element_blank()
        ) +
        coord_fixed() +
        force_panelsizes(rows = unit(10, "cm"), cols = unit(10, "cm"))
      
      ggsave('correlation.png', height = 5, width = 5.5, device = "png", bg = "white")
      
      # True precision matrix
      corr <- Sigma_inv
      mat_melt <- melt(corr)
      colnames(mat_melt) <- c("Row", "Col", "Value")
      mat_melt$Value <- factor(mat_melt$Value)
      
      unique_vals <- unique(mat_melt$Value)
      other_val <- unique_vals[!unique_vals %in% c(0, 2)]
      color_mapping <- c("2" = "firebrick3", "0" = "white")
      if (length(other_val) > 0) {
        color_mapping[as.character(other_val)] <- ifelse(offdiag_sign == "neg", "steelblue", "pink")
      }
      
      n_show <- 5
      axis_labels <- levels(factor(mat_melt$Col))
      thin_labels <- ifelse(seq_along(axis_labels) %% n_show == 1, axis_labels, "")
      
      p <- ggplot(mat_melt, aes(x = factor(Col), y = factor(Row), fill = factor(Value))) +
        geom_tile(color = "white", linewidth = 0.1) +
        scale_fill_manual(
          values = color_mapping,
          breaks = names(color_mapping)
        ) +
        scale_x_discrete(labels = thin_labels) +
        scale_y_discrete(labels = thin_labels) +
        labs(x = 'P', y = 'P', fill = "Value") +
        theme_bw() +
        theme(
          axis.title   = element_text(size = 14, face = "bold"),
          axis.text.x  = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
          axis.text.y  = element_text(size = 12),
          legend.title = element_blank(),
          legend.text  = element_text(size = 9),
          panel.grid   = element_blank()
        ) +
        coord_fixed() +
        force_panelsizes(rows = unit(10, "cm"), cols = unit(10, "cm"))
      
      ggsave('precision.png', p, height = 5, width = 5.5, dpi = 300)
      
    }
  }

  all_results_df <- cbind(do.call(rbind,all_results), 
                          all_settings[rep(which(lengths(all_results) > 0L),
                                           times = vapply(all_results[lengths(all_results) > 0L], nrow, integer(1))), ]
  )

  variables <- c("a","nu","nugget","norm", "coverage", "rmsesum",
                 "coveragesum","len","lensum","time","FPR","TPR","FDR","F2")
  all_results_list <- list()
  for (variable in variables) {
    for (setting_ind in plot_set) {
      prec_structure <- as.character(all_settings[setting_ind, 'prec_structure'])
      offdiag_sign   <- as.character(all_settings[setting_ind, 'offdiag_sign'])
      missing_TF     <- as.character(all_settings[setting_ind, 'missing_TF'])
      
      ind <- all_results_df$prec_structure == prec_structure &
        all_results_df$offdiag_sign == offdiag_sign &
        all_results_df$missing_TF == missing_TF
      results_df <- all_results_df[ind, ]
      
      out <- results_df %>%
        pivot_longer(cols = starts_with(paste0(variable, "_"))) %>%
        mutate(
          method = case_when(
            str_starts(name, paste0(variable, "_BIC"))   ~ "BIC",
            str_starts(name, paste0(variable, "_dense")) ~ "Dense",
            str_starts(name, paste0(variable, "_true"))        ~ "Ind",
            str_starts(name, paste0(variable, "_lmc"))         ~ "LMC",
            TRUE                                               ~ name
          )
        ) %>%
        group_by(method) %>%
        summarise(
          avg_len = mean(value, na.rm = TRUE),
          sd_len  = sd(value, na.rm = TRUE),
          n       = n(),
          .groups = "drop"
        ) %>%
        mutate(
          variable       = variable,
          prec_structure = prec_structure,
          offdiag_sign   = offdiag_sign,
          missing_TF     = missing_TF
        )
      all_results_list[[length(all_results_list) + 1]] <- out
    }
  }

  combined_df <- bind_rows(all_results_list) %>%
    arrange(desc(variable), desc(missing_TF), desc(prec_structure), desc(offdiag_sign), desc(method))

  combined_df <- combined_df %>%
    mutate(summary_val = avg_len)

  methods <- c("BIC", "Dense","LMC","Ind")
  ordered_cols <- as.vector(c(
    paste0(methods, "_summary_val"),
    paste0(methods, "_sd_val")
  ))

  wide_df <- combined_df %>%
    mutate(
      summary_val = avg_len,
      sd_val      = sd_len
    ) %>%
    unite("setting", missing_TF, prec_structure, offdiag_sign, sep = "_") %>%
    select(variable, setting, method, summary_val, sd_val) %>%
    pivot_wider(
      id_cols     = c(variable, setting),
      names_from  = method,
      values_from = c(summary_val, sd_val),
      names_glue  = "{method}_{.value}"
    ) %>%
    separate(setting, into = c("missing_TF", "prec_structure", "offdiag_sign"), sep = "_") %>%
    arrange(variable, missing_TF, desc(prec_structure), desc(offdiag_sign)) %>%
    select(variable, missing_TF, prec_structure, offdiag_sign, any_of(ordered_cols))

  # Estimation error against each process's true parameters.
  truth_by_process <- function(par) vapply(seq_len(P), function(j) {
    v <- f_par_true[[j]]
    if (!is.null(v) && par %in% names(v)) as.numeric(v[[par]]) else NA_real_
  }, numeric(1))

  par_truth <- list(a      = truth_by_process("a"),
                    nu     = truth_by_process("nu"),
                    nugget = rep_len(nugget_true, P))

  err_list <- list()
  for (variable in names(par_truth)) {
    for (setting_ind in plot_set) {
      prec_structure <- as.character(all_settings[setting_ind, 'prec_structure'])
      offdiag_sign   <- as.character(all_settings[setting_ind, 'offdiag_sign'])
      missing_TF     <- as.character(all_settings[setting_ind, 'missing_TF'])

      ind <- all_results_df$prec_structure == prec_structure &
        all_results_df$offdiag_sign == offdiag_sign &
        all_results_df$missing_TF == missing_TF
      results_df <- all_results_df[ind, ]

      # Read the truth before mutate() creates a column named `variable`.
      truth_vec <- par_truth[[variable]]

      # A single estimate (e.g. a shared nugget) is compared only when all process
      # truths agree; LMC indices are latent factors, not processes.
      shared_truth <- if (length(unique(truth_vec[!is.na(truth_vec)])) == 1L)
        truth_vec[!is.na(truth_vec)][1] else NA_real_

      truth_lab <- if (!any(is.finite(truth_vec))) NA_character_ else
        if (!is.na(shared_truth)) sprintf("%.3g", shared_truth) else
          sprintf("%.2f-%.2f", min(truth_vec, na.rm = TRUE),
                  max(truth_vec, na.rm = TRUE))

      err_list[[length(err_list) + 1]] <- results_df %>%
        pivot_longer(cols = starts_with(paste0(variable, "_"))) %>%
        mutate(
          method = case_when(
            str_starts(name, paste0(variable, "_BIC"))   ~ "BIC",
            str_starts(name, paste0(variable, "_dense")) ~ "Dense",
            str_starts(name, paste0(variable, "_true"))        ~ "Ind",
            str_starts(name, paste0(variable, "_lmc"))         ~ "LMC",
            TRUE                                               ~ NA_character_
          ),
          # trailing integer = process; none = a single shared estimate
          process = suppressWarnings(as.integer(sub("^[^0-9]*", "", name))),
          err = value - ifelse(is.na(process) | method == "LMC",
                               shared_truth, truth_vec[process])
        ) %>%
        filter(!is.na(method), is.finite(err)) %>%
        group_by(method) %>%
        summarise(
          bias = mean(err),
          # pooled over processes, whose true values differ
          sd   = stats::sd(err),
          .groups = "drop"
        ) %>%
        mutate(
          variable = variable,
          truth    = truth_lab,
          setting  = paste(missing_TF, prec_structure, offdiag_sign, sep = "_")
        )
    }
  }

  ordered_err <- unlist(lapply(c("BIC", "Dense", "LMC", "Ind"),
                               function(m) paste0(c("bias_", "sd_"), m)))

  err_wide <- bind_rows(err_list) %>%
    pivot_wider(
      id_cols     = c(variable, setting, truth),
      names_from  = method,
      values_from = c(bias, sd),
      names_glue  = "{.value}_{method}"
    ) %>%
    separate(setting, into = c("missing_TF", "prec_structure", "offdiag_sign"), sep = "_") %>%
    arrange(variable, missing_TF, desc(prec_structure), desc(offdiag_sign)) %>%
    select(variable, missing_TF, prec_structure, offdiag_sign, truth,
           any_of(ordered_err))

  cat("\n================ error against the truth: bias and SD, pooled over ",
      "processes and replicates ================\n\n", sep = "")
  old_width <- options(width = 200)
  print(as.data.frame(err_wide %>% mutate(across(where(is.numeric), \(z) round(z, 4)))),
        row.names = FALSE, na.print = "-")
  options(old_width)

  # Prediction table
  pred_rows <- c(
    norm        = "Overall: residual norm ||y - yhat||_2",
    coverage    = "Individual processes: 95% coverage",
    len         = "Individual processes: 95% interval length",
    rmsesum     = "Sum of processes: RMSE",
    coveragesum = "Sum of processes: 95% coverage",
    lensum      = "Sum of processes: 95% interval length")

  cat("\n================ prediction and interval summaries: ",
      "mean (sd) over replicates ================\n\n", sep = "")
  pred_tab <- do.call(rbind, lapply(names(pred_rows), function(v) {
    r <- wide_df[wide_df$variable == v, , drop = FALSE]
    if (!nrow(r)) {
      out <- data.frame(metric = unname(pred_rows[[v]]), stringsAsFactors = FALSE)
      for (m in methods) out[[m]] <- NA_character_
      return(out)
    }
    out <- data.frame(metric = unname(pred_rows[[v]]), stringsAsFactors = FALSE)
    for (m in methods) {
      sv <- paste0(m, "_summary_val"); sdv <- paste0(m, "_sd_val")
      out[[m]] <- if (sv %in% names(r) && !is.na(r[[sv]][1]))
        sprintf("%.3f (%.3f)", r[[sv]][1], r[[sdv]][1]) else NA_character_
    }
    out
  }))
  old_width <- options(width = 200)
  print(pred_tab, row.names = FALSE, na.print = "NA")
  options(old_width)
  if (all(is.na(pred_tab[["BIC"]]))) {
    message("  prediction rows are empty -- run sim_ll.R, sim_BIC.R and ",
            "sim_dense.R for this setting first")
  }

  # Save the tables
  table_dir <- here("sim", "output", "tables", setting_name)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(pred_tab, file.path(table_dir, "prediction.csv"), row.names = FALSE)
  write.csv(err_wide, file.path(table_dir, "estimation_error.csv"), row.names = FALSE)
  time_rows <- wide_df[wide_df$variable == "time", , drop = FALSE]
  time_table <- do.call(rbind, lapply(methods, function(method) {
    mn <- paste0(method, "_summary_val"); sd <- paste0(method, "_sd_val")
    if (!nrow(time_rows) || !mn %in% names(time_rows)) return(NULL)
    data.frame(stage = "estimation", method = method,
               mean_minutes = time_rows[[mn]][1],
               sd_minutes = time_rows[[sd]][1])
  }))
  for (method in c("BIC", "dense")) {
    minutes <- vapply(seq_len(n_sim), function(sim) {
      file <- file.path(out_dir, paste0(method, "_", setting_name, "_", sim, ".RData"))
      if (!file.exists(file)) return(NA_real_)
      e <- new.env(parent = emptyenv())
      load(file, envir = e)
      if (is.null(e$time_pred)) return(NA_real_)
      as_mins(e$time_pred)
    }, numeric(1))
    time_table <- rbind(time_table, data.frame(
      stage = "prediction", method = if (method == "dense") "Dense" else method,
      mean_minutes = mean(minutes, na.rm = TRUE),
      sd_minutes = sd(minutes, na.rm = TRUE)))
  }
  cat("\n================ computation time (minutes): mean and SD ================\n")
  print(time_table, row.names = FALSE)
  write.csv(time_table, file.path(table_dir, "computation_time.csv"), row.names = FALSE)

  wide_df

}
wide_df <- summarize_results()

# Trace diagnostics for the BIC-selected SAECM fits. Microergodic and precision
# plots are made directly from the saved chains. Observed likelihood evaluation
# is substantially more expensive, so missing cache points are computed only
# when SSGGP_COMPUTE_DIAGNOSTIC_LOGLIK=true; completed points are reusable.
diagnostic_script <- here("sim", "result", "sim_diagnostic_plots.R")
sys.source(diagnostic_script, envir = environment())
compute_diagnostic_loglik <- tolower(trimws(
  Sys.getenv("SSGGP_COMPUTE_DIAGNOSTIC_LOGLIK", "false"))) %in%
  c("1", "true", "yes")
diagnostic_files <- plot_simulation_diagnostics(
  out_dir = out_dir,
  fig_dir = file.path(fig_dir, "diagnostics"),
  setting_name = setting_name,
  grid_params = grid_params,
  n_sim = n_sim,
  P = P,
  family_names = family_names,
  n_iter = iter,
  L_ll = L_ll,
  seed_ll = seed_ll,
  threads = threads,
  compute_observed_ll = compute_diagnostic_loglik)

if (run_family_result) {
  message("non-Matern or mixed families detected (",
          paste(family_names, collapse = ", "),
          "); calculating family-specific parameter summaries")
  source(est_result_script, chdir = TRUE)
}

as.data.frame(wide_df[wide_df$variable == "a", ])
