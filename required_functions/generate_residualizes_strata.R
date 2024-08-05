generate_residualized_strata <- function(iv_free_exposure, k){
  
  if(length(k) == 1){
    q = quantile( iv_free_exposure, probs = seq(0, 1, 1/k), na.rm = TRUE )
  } else {
    if(class(k) == "numeric" & length(k) >=3 ){
      q = k
    } else {
      stop("please check strata parameter. Acceptable values are a single numeric value indicating the number of quantiles, or a numeric vector of at least length 3 to define strata boundries.")
    }
  }
  
  ## identify samples of each strata
  ## and add a strata label to the model data frame
  strata = as.factor( cut(iv_free_exposure, q, include.lowest=TRUE, labels=FALSE) )
  
  ## Return to user
  return( strata )
}





