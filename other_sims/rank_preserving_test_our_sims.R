# Clear workspace and load required libraries
rm(list=ls())
library(tidyverse)
library(broom)
library(skedastic)
library(furrr)
library(future)

# Set up parallel processing
plan(multisession, workers = parallel::detectCores() - 1)

n <- 5e4  # Sample size
n.sims <- 100  # Number of simulated datasets
res <- NULL

# Define all models with correct naming for order
# Define all models with correct naming for order and descriptions
sim_models <- list(
  model1 = list(func = function(G.in, U.in, eX.in) 0.3 * G.in + U.in + eX.in, bux = 1, bgux = 0, description = "Linear model with G and U"),
  model2 = list(func = function(G.in, U.in, eX.in) 0.3 * G.in + 0.1 * U.in + eX.in, bux = 0.1, bgux = 0, description = "Linear model with reduced U effect"),
  model3 = list(func = function(G.in, U.in, eX.in) 0.3 * G.in + 1 * U.in + 0.1 * G.in * U.in + eX.in, bux = 1, bgux = 0.1, description = "Interaction model with small GxU effect"),
  model4 = list(func = function(G.in, U.in, eX.in) 0.3 * G.in + 1 * U.in + 0.2 * G.in * U.in + eX.in, bux = 1, bgux = 1, description = "Interaction model with large GxU effect"),
  model5 = list(func = function(G.in, U.in, eX.in) 0.3 * G.in + 0.1 * G.in * U.in + eX.in, bux = 0, bgux = 0.1, description = "Interaction model without main U effect"),
  model6 = list(func = function(G.in, U.in, eX.in) 0.3 * G.in + -0.2 * G.in * U.in + U.in + eX.in, bux = 1, bgux = 1, description = "Interaction model with negative GxU effect"),
  model7 = list(func = function(G.in, U.in, eX.in) 0.3 * G.in + 0.1 * U.in + 0.5 * G.in * U.in + eX.in, bux = 0.1, bgux = 0.5, description = "Interaction model with large GU effect but small U effect"),
  model8 = list(func = function(G.in, U.in, eX.in) 0.3 * G.in + U.in + 0.05 * G.in * U.in + eX.in, bux = 1, bgux = 0.05, description = "Interaction model with positive small GxU effect"),
  model9 = list(func = function(G.in, U.in, eX.in) {
    bg <- 1  # Assuming bg = 1, adjust if needed
    bux <- 1 # Assuming bux = 1, adjust if needed
    X <- bg * G.in + bux * U.in
    eX <- rnorm(length(G.in), 0, 0.2 * sqrt(abs(G.in) + 1))  # Variance effect
    return(X + eX)
  }, bux = 1, bgux = NA, description = "Variance Effect Model"),
  model10 = list(func = function(G.in, U.in, eX.in) {
    0.3 * G.in + 0.1 * G.in^2 + U.in + eX.in
  }, bux = 1, bgux = NA, description = "G^2 Model"),
  model11 = list(func = function(G.in, U.in, eX.in) {
    -10 + (1.5 + 0.4 * U.in) * (G.in + 5) + U.in + eX.in
  }, bux = NA, bgux = NA, description = "PLOS Genetics"),
  model12 = list(func = function(G.in, U.in, eX.in) {
    alpha <- rnorm(length(G.in), mean = 0.3 + 0.1 * eX.in, sd = 0.1)
    alpha * G.in + U.in + eX.in
  }, bux = NA, bgux = NA, description = "HH: X only"),
  model13 = list(func = function(G.in, U.in, eX.in) {
    alpha <- rnorm(length(G.in), mean = 0.3 + 0.1/sqrt(2) * eX.in + 0.1/sqrt(2) * U.in, sd = 0.1)
    alpha * G.in + U.in + eX.in
  }, bux = NA, bgux = NA, description = "HH: X and U equally"),
  model14 = list(func = function(G.in, U.in, eX.in) {
    alpha <- rnorm(length(G.in), mean = 0.3 + 0.1 * U.in, sd = 0.1)
    alpha * G.in + U.in + eX.in
  }, bux = NA, bgux = NA, description = "HH: U only")
)

# The rest of the script remains the same

# Function to assess GxU interaction
assess_GxU_interaction <- function(model, n = 10000) {
  if (model$description == "PLOS Genetics") {
    G <- rnorm(n, mean = 0, sd = 0.25)  # SD of 0.5^2 for PLOS Genetics
  } else {
    G <- rnorm(n, mean = 0, sd = 1)  # SD of 1 for all other models
  }
  
  U <- rnorm(n, mean = 0, sd = 1)
  eX <- rnorm(n, mean = 0, sd = 1)
  
  X <- model$func(G, U, eX)
  
  fit_main <- lm(X ~ G + U)
  fit_interaction <- lm(X ~ G * U)
  
  r2_main <- summary(fit_main)$r.squared
  r2_interaction <- summary(fit_interaction)$r.squared
  var_explained_by_interaction <- r2_interaction - r2_main
  
  # Use broom::tidy to extract coefficients
  coefficients <- broom::tidy(fit_interaction)
  interaction_coefficient <- coefficients$estimate[coefficients$term == "G:U"]
  
  # Check if interaction_coefficient is empty and handle it
  if (length(interaction_coefficient) == 0) {
    interaction_coefficient <- 0  # or any other appropriate value
  }
  
  return(c(var_explained_by_interaction = var_explained_by_interaction,
           interaction_coefficient = interaction_coefficient))
}
# Function to generate data and calculate results for one simulation
run_simulation <- function(model, n) {
  # Generate G with SD of 1 for all models except PLOS Genetics
  if (model$description == "PLOS Genetics") {
    G <- rnorm(n, mean = 0, sd = 0.25)  # SD of 0.5^2 for PLOS Genetics
  } else {
    G <- rnorm(n, mean = 0, sd = 1)  # SD of 1 for all other models
  }
  
  U <- rnorm(n, mean = 0, sd = 1)
  eX <- rnorm(n, mean = 0, sd = 1)
  
  X <- model$func(G, U, eX)
  
  # Calculate X at -1SD, mean, and +1SD of G
  if (model$description == "PLOS Genetics") {
    X_neg <- model$func(rep(-0.25, n), U, eX)  # -1SD for PLOS Genetics
    X0 <- model$func(rep(0, n), U, eX)         # Mean
    X_pos <- model$func(rep(0.25, n), U, eX)   # +1SD for PLOS Genetics
  } else {
    X_neg <- model$func(rep(-1, n), U, eX)     # -1SD for other models
    X0 <- model$func(rep(0, n), U, eX)         # Mean
    X_pos <- model$func(rep(1, n), U, eX)      # +1SD for other models
  }
  
  X_neg.rank <- rank(X_neg)
  X0.rank <- rank(X0)
  X_pos.rank <- rank(X_pos)
  
  rank_diff_neg_0 <- mean(X_neg.rank != X0.rank)
  rank_diff_neg_pos <- mean(X_neg.rank != X_pos.rank)
  rank_diff_0_pos <- mean(X0.rank != X_pos.rank)
  
  mean_rank_diff_neg_0 <- mean(abs(X_neg.rank - X0.rank)) / n
  mean_rank_diff_neg_pos <- mean(abs(X_neg.rank - X_pos.rank)) / n
  mean_rank_diff_0_pos <- mean(abs(X0.rank - X_pos.rank)) / n
  
  gxu_results <- assess_GxU_interaction(model, 1e5)
  
  tibble(
    p_het = glejser(lm(X ~ G))[2] %>% as.numeric(),
    rank_diff_neg_0 = rank_diff_neg_0,
    rank_diff_neg_pos = rank_diff_neg_pos,
    rank_diff_0_pos = rank_diff_0_pos,
    mean_rank_diff_neg_0 = mean_rank_diff_neg_0,
    mean_rank_diff_neg_pos = mean_rank_diff_neg_pos,
    mean_rank_diff_0_pos = mean_rank_diff_0_pos,
    var_explained_by_interaction = gxu_results["var_explained_by_interaction"],
    interaction_coefficient = gxu_results["interaction_coefficient"]
  )
}

# Run simulations using future_map
results <- crossing(
  sim_model = names(sim_models),
  iteration = 1:n.sims
) %>%
  mutate(
    data = future_map2(sim_model, n, ~run_simulation(sim_models[[.x]], .y), .options = furrr_options(seed = TRUE))
  ) %>%
  unnest(data)


results$interaction_coefficient
# Prepare data for plotting
results_long <- results %>%
  pivot_longer(cols = -c(sim_model, iteration), names_to = "metric", values_to = "value")

# Define more descriptive labels for the plots
metric_labels <- c(
  "p_het" = "Heteroscedasticity p-value",
  "rank_diff_neg_0" = "Proportion Rank Change: G=-1 vs G=0",
  "rank_diff_neg_pos" = "Proportion Rank Change: G=-1 vs G=1",
  "rank_diff_0_pos" = "Proportion Rank Change: G=0 vs G=1",
  "mean_rank_diff_neg_0" = "Mean Normalized Rank Difference: G=-1 vs G=0",
  "mean_rank_diff_neg_pos" = "Mean Normalized Rank Difference: G=-1 vs G=1",
  "mean_rank_diff_0_pos" = "Mean Normalized Rank Difference: G=1 vs G=1",
  "var_explained_by_interaction" = "Variance Explained by GxU Interaction",
  "interaction_coefficient" = "GxU Interaction Coefficient"
)

# Create a factor with levels in the desired order
model_order <- paste0("M", 1:14)
results_long$sim_model <- factor(results_long$sim_model, levels = names(sim_models), labels = model_order)

# Create the plot
ggplot(results_long, aes(x = sim_model, y = value, fill = sim_model)) +
  geom_boxplot() +
  facet_wrap(~ metric, scales = "free_y", labeller = as_labeller(metric_labels), ncol = 3) +
  theme_bw() +
  labs(title = "Distribution of Metrics Across All Models",
       x = "Model", y = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave("metric_distributions_across_all_models_updated.pdf", width = 20, height = 10, bg = "white")

# Generate CSV for Word
model_descriptions <- tibble(
  Model = paste0("M", 1:14),
  Description = c(
    "Linear model with G and U",
    "Linear model with reduced U effect",
    "Interaction model with small GxU effect",
    "Interaction model with large GxU effect",
    "Interaction model without main U effect",
    "Interaction model with negative GxU effect",
    "Interaction model with large GU effect but small U effect",
    "Interaction model with positive small GxU effect",
    "Variance effect model",
    "Quadratic G model",
    "PLOS Genetics model",
    "Human Heredity model (X only)",
    "Human Heredity model (X and U equally)",
    "Human Heredity model (U only)"
  ),
  Equation = c(
    "X = 0.3G + U + ε",
    "X = 0.3G + 0.1U + ε",
    "X = 0.3G + U + 0.1GU + ε",
    "X = 0.3G + U + 0.2GU + ε",
    "X = 0.3G + 0.1GU + ε",
    "X = 0.3G - 0.2GU + ε",
    "X = 0.3G + 0.1U + 0.5GU + ε",
    "X = 0.3G + U + 0.05GU + ε",
    "X = G + U + ε, where ε ~ N(0, (0.2√(|G|+1))^2)",
    "X = 0.3G + 0.1G^2 + U + ε",
    "X = -10 + (1.5 + 0.4U)(G + 5) + U + ε",
    "X = αG + U + ε, where α ~ N(0.3 + 0.1ε, 0.1^2)",
    "X = αG + U + ε, where α ~ N(0.3 + 0.1/√2(ε + U), 0.1^2)",
    "X = αG + U + ε, where α ~ N(0.3 + 0.1U, 0.1^2)"
  )
)

model_descriptions %>%
  gt::gt()

# Write CSV file
write_csv(model_descriptions, "model_descriptions.csv")