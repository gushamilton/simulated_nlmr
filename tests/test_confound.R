library(DRMR)
library(ggforestplot)
library(tidyverse)

# Example demonstration of confounding in stratified analysis
n <- 1e5

# Generate example data with effect heterogeneity
dat <- tibble(Z = rnorm(n, 2, 1),
              U = rnorm(n, 2, 1),
              X = rnorm(n) + 2 * Z * U + 2 * Z,
              Y = 0.2 * X + rnorm(n) + 2 * U)

# Analyze stratified effects
rdat <- Stratify(dat)
res <- getSummaryInf(rdat)
plot <- res$DRres |> 
  as_tibble(rownames = "strat") |>
  mutate(strat = factor(strat, levels = paste("Stratum", 1:10))) |>
  forestplot(name = strat, estimate = est, se = se) +
  xlab("Causal effect")

plot

# Demonstrate confounding by simulating null Y with confounder
dat_test <- dat |> 
  mutate(V = rnorm(n),
         X = X + V, 
         Y = rnorm(n) + V)

rdat_conf <- Stratify(dat_test)
res_conf <- getSummaryInf(rdat_conf)
conf_plot <- res_conf$DRres |> 
  as_tibble(rownames = "strat") |>
  mutate(strat = factor(strat, levels = paste("Stratum", 1:10))) |>
  forestplot(name = strat, estimate = est, se = se) +
  xlab("Causal effect (all estimates should be null)")

conf_plot

# Function to test confounding with arbitrary inputs
test_confound <- function(Z, X, K = 100) {
  n <- length(Z)
  if(length(X) != n) stop("Z and X must be same length")
  
  # Store results across iterations
  all_results <- map(1:K, function(i) {
    # Generate confounder and null outcome
    dat <- tibble(Z = Z,
                 X = X,
                 V = rnorm(n),
                 Y = rnorm(n))
    
    # Add confounding
    dat <- dat |> 
      mutate(X = X + V,
             Y = Y + V)
    
    # Get stratified results
    rdat <- Stratify(dat)
    res <- getSummaryInf(rdat)
    
    # Return results for this iteration
    res$DRres |> 
      as_tibble(rownames = "strat")
  })
  
  # Combine and summarize results
  results_df <- bind_rows(all_results, .id = "iteration")
  
  # Calculate mean effects across iterations
  summary_plot <- results_df |> 
    group_by(strat) |> 
    summarise(
      estimate = mean(est),
      se = sd(est)/sqrt(K)
    ) |>
    mutate(strat = factor(strat, levels = paste("Stratum", 1:10))) |>
    forestplot(name = strat, estimate = estimate, se = se) +
    xlab("Mean causal effect across iterations (should be null)")
  
  return(summary_plot)
}

# Example usage:
# test_confound(Z = rnorm(1000), X = rnorm(1000), K = 100)

test_confound(Z = dat$Z, X = dat$X, K = 10)
