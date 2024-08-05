sim_mydata_outcomes_confounder <- function(n = 100000,
                                           mu = c(2, 2, 2, 2, 2), 
                                           sd = c(1, 1, 1, 1, 1), 
                                           r = 0, 
                                           varnames = c("g", "u", "v", "ex", "ey"),
                                           bux = 0.3,
                                           buy = 0.3,
                                           bgux = -0.1,
                                           bv = 0,
                                           seed = NULL,
                                           u_cor=0.8
){ 
  if (!is.null(seed)) set.seed(seed)
  
  mydata <- faux::rnorm_multi(
    n = n, 
    mu = mu,
    sd = sd,
    r = r, 
    varnames = varnames,
    empirical = FALSE
  )
  
  mydata$x <- with(mydata, 0.3 * g + bux * u + bgux * g * u + bv * v + ex)
  mydata$y <- with(mydata, buy * u + bv * v + ey)
  mydata$u_incomplete = faux::rnorm_pre(mydata$u,2,1, r= u_cor)
  return(mydata)
}

# Example usage
sim_mydata_outcomes_confounder(n = 100, u_cor = 0.1)

