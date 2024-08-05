library(furrr)
library(tidyverse)
library(SUMnlmr)

# Set up parallel processing
plan(multisession, workers = availableCores() - 1)


plos_g_function =  function(G.in, U.in, eX.in) {
  -10 + (1.5 + 0.4 * U.in) * (G.in + 5) + U.in + eX.in
}
replicate_plos_gen <- function(n_rep, n) {
  set.seed(n_rep)
  t <- tibble(
    g = rnorm(n, sd = 0.25),
    u = rnorm(n),
    e = rnorm(n),
    x = plos_g_function(g,u,e),
    y = u + rnorm(n)
  ) 
  summary <- SUMnlmr::create_nlmr_summary(t$y, t$x, t$g, strata_method = "ranked", q = 10)
  summary$summary %>%
    mutate(mr = by/bx, mr_se = byse/bx) %>%
    as_tibble() %>%
    mutate(strata = as.factor(row_number())) %>%
    mutate(replicate = n_rep, sample_size = n)
}



# Use future_map_dfr instead of map_dfr - this will take some time!
replicates_1e4 <- future_map_dfr(1:100, ~replicate_plos_gen(.x, n = 1e4), .options = furrr_options(seed = TRUE))
replicates_1e5 <- future_map_dfr(1:100, ~replicate_plos_gen(.x, n = 1e5), .options = furrr_options(seed = TRUE))
replicates_5e5 <- future_map_dfr(1:100, ~replicate_plos_gen(.x, n = 3.5e5), .options = furrr_options(seed = TRUE))

# For the 1000 replications

# The plotting code remains the same
replicates_1e4 %>%
    bind_rows(replicates_1e5) %>%
    bind_rows(replicates_5e5) %>%
  mutate(sample_size = case_when(
    sample_size == "10000" ~ "10k",
    sample_size == "1e+05" ~ "100k",
    sample_size == "350000" ~ "350k",
    TRUE ~ as.character(sample_size)
  )) %>%
  mutate(sample_size = factor(sample_size, levels = c("10k", "100k", "350k"))) %>%
  ggplot(aes(x = strata, y = mr)) +
  geom_hline(aes(yintercept = 0), col = "red", lty = "dashed") +
  geom_boxplot(outliers = F) +
  facet_wrap( ~ sample_size, ncol = 1, scales = "free") +
  ylab("Estimates from MR (interquartile range)") +
  xlab("Doubly-ranked strata") +
  theme_bw()
