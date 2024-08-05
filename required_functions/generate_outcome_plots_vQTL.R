generate_outcome_plots_vQTL<- function(exposure, outcome, seed = "123", reps = 10,
                                       g_sd = 0.1, u_sd = 0,

                                                                bux = 0, bx = 0, buy = 0, covariates = NULL) {
  

  replicate_internally <- function(n_rep) {
  seed = seed + n_rep
  dat <-  sim_mydata_outcomes_vQTL(n = 100000, seed = seed,
                                   g_sd = g_sd,

                                         bux= bux, bx = bx, buy = buy)
  dat_with_both <- generate_all_sumstats(data = dat, exposure = exposure, outcome = outcome, k = 10, covariates = covariates)
  
  ## redefine data set
  
  
  ## define sumstats tibble
  sum_stats_dat = dat_with_both$ss
  
  sum_stats_dat %>%
    mutate(replicate = n_rep)
  }
  
  d <- map_dfr(1:reps, replicate_internally)
 
  dat <-  sim_mydata_outcomes_vQTL(n = 100000, seed = seed)
  
  ## make figures


  make_figures_replicates(d,exposure,outcome,dat,reps, "figures/linear/", paste0("variance_alone_","bx=",bx, "_bux=",bux, "_buy=", buy, "_b_var=",g_sd), covariates = covariates) 
}





