############################
## D.Hughes function
## to make an IV free exposure
## vector
############################
iv_free_exposure = function( wdata,
                             exposure,
                             instrument,
                             covariates = NULL,
                             exposure_mean_normalize = TRUE){
  
  ######################
  ### I. Linear model
  ######################
  if( !is.null(covariates) ){
    form = formula(paste0(exposure, " ~ ", paste0( covariates, collapse = " + ") , " + ", instrument ))
  } else {
    form = formula(paste0(exposure, " ~ ", instrument ))
  }
  
  #########################
  ## II. RUN the LINEAR MODEL
  #########################
  index = rownames(wdata)
  lm_mod = lm(form, data = wdata)
  
  #########################
  ## III. Extract the residuals
  #########################
  ## IV free exposures a.k.a the residuals
  res = residuals(lm_mod)
  
  ## accounting for NAs in model, to insure length res == nrow(wdata)
  m = match(index, names(res))
  res = res[m]; names(res) = index
  
  #########################
  ## IV. mean normalize the
  ##     residuals
  #########################
  if(exposure_mean_normalize == TRUE){
    exposure_mean = mean(wdata[, exposure], na.rm = TRUE)
    res = res + exposure_mean
  }
  
  #########################
  ## IV. add IV free exposure
  ##     to working data frame
  #########################
  wdata$iv_free_exposure = res
  
  #########################
  ## VI. Return to user
  #########################
  # return(wdata)
  return(res)
  
}
