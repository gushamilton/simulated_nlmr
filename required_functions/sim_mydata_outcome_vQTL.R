sim_mydata_outcomes_vQTL <- function(n = 100000,
                                     mu = c(0.5, 0.5, 0.5, 0.5, 0.5), 
                                     sd = c(1, 1, 1, 1, 1), 
                                     r = 0, 
                                     varnames = c("g", "ex", "ey", "u", "v"),
                                     bux = 0,
                                     bx = 0,
                                     buy = 0,
                                     g_sd = 0,
                                     u_sd = 0,
                                     seed = seed,
                                     covariates = NULL
){ 
  
  ########################
  ## Set seed and define helper function
  ########################
  set.seed(seed)
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  
  ########################
  ## Simulate the data 
  ########################
  mydata <- faux::rnorm_multi(n = n, 
                              mu = mu,
                              sd = sd,
                              r = r, 
                              varnames = varnames,
                              empirical = FALSE) %>%
    mutate(g_rank = range01(g),
           u_rank = range01(u),
           x = bx * g + rnorm(n, mean = 0, sd = u_rank*u_sd + g_rank*g_sd + 1) + bux * v,
           y = 0 * x + buy * v + ey)
  
  ########################
  ## return the sims
  ########################
  return(mydata)
}
