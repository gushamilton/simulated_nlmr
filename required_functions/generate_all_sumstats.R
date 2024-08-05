generate_all_sumstats = function(
    datain = NULL,
    exposure = "x",
    outcome = "y",
    covariates = NULL,
    k = 10){
  
  #####################
  ## make working data frame
  #####################
  wdata = datain %>% mutate(iv_free_x = iv_free_exposure( wdata = .,
                                                          exposure = exposure,
                                                          instrument = "g",
                                                          covariates = covariates,
                                                          exposure_mean_normalize = TRUE) ) %>% 
    mutate( obs_strata = generate_residualized_strata(iv_free_exposure = .data[[exposure]], k = k ) ) %>%
    mutate( residual = generate_residualized_strata(iv_free_exposure = iv_free_x, k = k ) ) %>% 
    mutate( ranked = generate_ranked_strata(instrument = g, exposure = .data[[exposure]], k = k) )
  
  # Function to create formula with covariate
  create_formula <- function(dependent, independent) {
    if (is.null(covariates)) {
      as.formula(paste(dependent, "~", independent))
    } else {
      as.formula(paste(dependent, "~", independent, "+", covariates))
    }
  }
  
  #####################
  ## Summary Statistics
  ## for the whole data set
  #####################
  f = create_formula(exposure, "g")
  f2 = create_formula(outcome, "g")
  f3 = create_formula(outcome, exposure)
  fit = lm(f, data = wdata)
  fit2 = lm(f2, data = wdata)
  fit3 = lm(f3, data = wdata)
  ###
  full_ss = data.frame( 
    strata = 0,
    v = var( wdata[, exposure] ),
    sd = sd( wdata[, exposure] ),
    skew = psych::skew( wdata[, exposure] ),
    kurtosis = psych::kurtosi( wdata[, exposure] ),
    min = min( wdata[, exposure] ),
    max = max( wdata[, exposure] ),
    mean = mean( wdata[, exposure]),
    range = max( wdata[, exposure] ) - min( wdata[, exposure] ),
    r2_gene_exposure = summary(fit)$r.sq,
    r2_gene_outcome = summary(fit2)$r.sq,
    r2_exposure_outcome = summary(fit3)$r.sq,
    beta_gx = summary(fit)$coef[2,1] , 
    se_gx = summary(fit)$coef[2,2] , 
    p_gx = summary(fit)$coef[2,4] ,
    beta_gy = summary(fit2)$coef[2,1] , 
    se_gy = summary(fit2)$coef[2,2] , 
    p_gy = summary(fit2)$coef[2,4] ,
    p_xy = summary(fit3)$coef[2,4] ,
    strata_method = "full"
  )
  
  ####################################
  ## Ranked Strata Summary Statistics
  ####################################
  
  compute_stats <- function(df, exposure, outcome) {
    f1 <- create_formula(exposure, "g")
    f2 <- create_formula(outcome, "g")
    f3 <- create_formula(outcome, exposure)
    
    fit1 <- lm(f1, data = df)
    fit2 <- lm(f2, data = df)
    fit3 <- lm(f3, data = df)
    
    list(
      v = var(df[[exposure]]),
      sd = sd(df[[exposure]]),
      skew = psych::skew(df[[exposure]]),
      kurtosis = psych::kurtosi(df[[exposure]]),
      min = min(df[[exposure]]),
      max = max(df[[exposure]]),
      mean = mean(df[[exposure]]),
      range = max(df[[exposure]]) - min(df[[exposure]]),
      r2_gene_exposure = summary(fit1)$r.sq,
      r2_gene_outcome = summary(fit2)$r.sq,
      r2_exposure_outcome = summary(fit3)$r.sq,
      beta_gx = coef(fit1)["g"],
      beta_gy = coef(fit2)["g"],
      se_gx = summary(fit1)$coef["g", "Std. Error"],
      se_gy = summary(fit2)$coef["g", "Std. Error"],
      p_gx = summary(fit1)$coef["g", "Pr(>|t|)"],
      p_gy = summary(fit2)$coef["g", "Pr(>|t|)"],
      p_xy = summary(fit3)$coef[exposure, "Pr(>|t|)"]
    )
  }
  
  ranked_ss <- wdata %>%
    group_by(ranked) %>%
    group_modify(~ as_tibble(compute_stats(., exposure, outcome))) %>%
    ungroup()
  
  ranked_ss$strata <- as.integer(ranked_ss$ranked)
  ranked_ss$strata_method <- "ranked"
  
  residual_ss <- wdata %>%
    group_by(residual) %>%
    group_modify(~ as_tibble(compute_stats(., exposure, outcome))) %>%
    ungroup()
  
  residual_ss$strata <- as.integer(residual_ss$residual)
  residual_ss$strata_method <- "residual"
  
  #  GENERATE Q and P
  
  ss_out <- bind_rows(full_ss, ranked_ss, residual_ss)
  
  ss_out <- ss_out %>%
    mutate(beta_mr = beta_gy/beta_gx, se_mr = se_gy/beta_gx, p_mr = pnorm(abs(beta_mr/se_mr), lower.tail=F))
  
  ss_out$strata <- factor(ss_out$strata, levels = 0:max(ss_out$strata))
  
  # Reorder columns to match original output
  ss_out <- ss_out %>%
    select(strata, v, sd, skew, kurtosis, min, max, mean, range, r2_gene_exposure, r2_gene_outcome, 
           r2_exposure_outcome, beta_gx, se_gx, p_gx, beta_gy, se_gy, p_gy, p_xy, strata_method, 
           beta_mr, se_mr, p_mr)
  
  out <- list(wdata = wdata, 
              ss = ss_out)
  return(out)
}
