make_figures_replicates <- function(d, exposure, outcome, dat, reps, location, name, covariates) {
  
  sum_stats_dat <- d %>%
    group_by(strata, strata_method) %>%
    summarise_all(., mean) %>%
    ungroup()
  
  get_pvals = function(df){
    p_Q <- 1 - pchisq(metafor::rma(beta_mr, vi=(se_mr)^2, data = df, control=list(maxiter=1000))$QE, df=(9))
    tibble(Q_p = p_Q)
  }
  
  p_res <- d %>%
    filter(strata_method != "full") %>%
    group_by(strata_method, replicate) %>%
    nest() %>%
    mutate(res = map(data, get_pvals)) %>%
    unnest(res) %>%
    ungroup()
  
  ranked_ps <- p_res %>%
    filter(strata_method == "ranked") %>%
    mutate(Q_p = if_else(Q_p ==0 , 0, Q_p))
  resid_ps <- p_res %>%
    filter(strata_method == "residual") %>%
    mutate(Q_p = if_else(Q_p == 0, 0, Q_p))
  
  r2_gene_exposure <- round(sum_stats_dat$r2_gene_exposure[1], digits = 4)
  r2_gene_outcome <- round(sum_stats_dat$r2_gene_outcome[1], digits = 4)
  r2_exposure_outcome <- round(sum_stats_dat$r2_exposure_outcome[1], digits = 4)
  
  p_gene_exposure <- signif(sum_stats_dat$p_gx[1], digits = 4)
  p_gene_outcome <- signif(sum_stats_dat$p_gy[1], digits = 4)
  p_exposure_outcome <- signif(sum_stats_dat$p_xy[1], digits = 4)
  
  # Prepare data with overall effect as first stratum
  d_with_overall <- d %>%
    mutate(strata = as.numeric(as.character(strata))) %>%
    bind_rows(
      d %>% 
        filter(strata_method == "full") %>% 
        mutate(strata = -1, strata_method = "ranked")
    ) %>%
    bind_rows(
      d %>% 
        filter(strata_method == "full") %>% 
        mutate(strata = -1, strata_method = "residual")
    ) %>%
    mutate(strata = factor(strata, levels = c(-1, 0:10), 
                           labels = c("Overall", as.character(0:10))))
  
  # Update plots to include overall effect as first stratum
  p1 <- d_with_overall %>%
    filter(strata_method != "full") %>% 
    ggplot(aes(x = strata, y = beta_gx)) +
    geom_boxplot() +
    theme_bw() +
    geom_hline(yintercept = sum_stats_dat$beta_gx[1], linetype = "dashed", color = "blue") +
    facet_wrap(~strata_method) +
    labs(title = "Gene-exposure effect estimates across strata",
         x = "Strata",
         y = "Betas from replicates (interquartile ranges)")
  
  p2 <- d_with_overall %>%
    filter(strata_method != "full") %>% 
    ggplot(aes(x = strata, y = beta_gy)) +
    geom_boxplot() +
    theme_bw() +
    geom_hline(yintercept = sum_stats_dat$beta_gy[1], linetype = "dashed", color = "blue") +
    facet_wrap(~strata_method) +
    labs(title = "Gene-outcome effect estimates across strata",
         x = "Strata",
         y = "Betas from replicates (interquartile ranges)")
  
  # Calculate MSE for MR estimates
  mse_results <- d_with_overall %>%
    filter(strata_method != "full") %>%
    group_by(strata_method, strata) %>%
    summarise(
      mse_mr = mean((beta_mr - 0)^2),  # 0 is the true value for MR
      .groups = 'drop'
    ) %>%
    mutate(label_mse_mr = ifelse(mse_mr < 0.001, "<0.001", sprintf("%.3f", mse_mr)))
  
  
  # Function to add MSE labels to MR plot
  add_mse_labels <- function(p, mse_data) {
    p + scale_x_discrete(labels = setNames(
      paste0(mse_data$strata, "\n", mse_data$label_mse_mr),
      mse_data$strata
    ))
  }
  
  p3 <- d_with_overall %>%
    filter(strata_method != "full") %>% 
    ggplot(aes(x = strata, y = beta_mr)) +
    geom_boxplot() +
    theme_bw() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
    facet_wrap(~strata_method) +
    labs(title = "MR effect estimates across strata",
         x = "Strata (MSE below)",
         y = "Betas from replicates (interquartile ranges)",
         caption = paste0("Cochrane FDR at p = 0.05; Ranked: ",  
                          signif(mean(ranked_ps$Q_p <0.05),3)," Residual: ",
                          signif(mean(resid_ps$Q_p < 0.05),3)))
  
  # Add MSE labels to MR plot
  p3 <- add_mse_labels(p3, mse_results)
  
  # Create the base filename and caption
  base_filename <- paste0(location, name, "_", reps, "_reps_with_overall_mse")
  base_caption <- name
  
  # Modify filename and caption if covariates are provided
  if (!is.null(covariates)) {
    base_filename <- paste0(base_filename, "_adjusted_for_", paste(covariates, collapse = "_"))
    base_caption <- paste0(base_caption, " (adjusted for ", paste(covariates, collapse = ", "), ")")
  }
  
  # Update master plot with potentially modified caption
  master_plot <- (p1 / p2 / p3) + 
    plot_annotation(tag_levels = "A", caption = base_caption)
  
  # Save the updated master plot with potentially modified filename
  ggsave(master_plot, filename = paste0(base_filename, ".pdf"), width = 10, height = 14)
  
  # Return the MSE results
  return(mse_results)
}
  # Return the MSE results
