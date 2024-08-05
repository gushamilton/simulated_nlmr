

generate_ivs <- function(all_snps, exposure, covariates = NULL, p_thresh) {
  
  # run and get the betas
  get_betas<- function(name,exposure = exposure, covariates = NULL){
    
    
    just_bmi <- d %>%
      select(exposure, IID, covariates) %>%
      distinct(.keep_all = T)
    
    d <- all_snps %>%
      select(IID, contains(name)) %>%
      select(!contains("HET")) %>%
      rename("SNP" = 2) %>%
      right_join(just_bmi) %>%
      left_join(snp)
    
    if (is.null(covariates) || length(covariates) == 0) {
      f <- as.formula(paste(exposure, "~ SNP"))
    } else {
      f <- as.formula(paste(exposure, paste(covariates, collapse = " + "), sep = " ~ SNP + "))
    }
    
    m <- lm(f, data = d)
    tidy(m) %>%
      bind_cols(skedastic::glejser(m) %>%
                  transmute(gles = p.value)) %>%
      filter(term == "SNP") %>%
      mutate(SNP = name)
    
  }


  get_betas("rs12", exposure = "body_mass_index_bmi", covariates = c("age_at_recruitment"))
  
  all_snp_names <- colnames(all_snps) %>%
    str_subset("rs") %>%
    str_subset("HET", negate = T)
  
  betas <- map_dfr(all_snp_names, get_betas, exposure, covariates)  
  
  # generate PRS
  
  generate_prs_internal <- function(p_thresh) {
    
    bs <- betas %>%
    filter(gles > p_thresh) 
  
  # Ensure SNP columns in all_snps have the same order as in betas
  snps_data <- all_snps %>% select(IID, one_of(bs$SNP))
  
  # Convert SNP data and betas to matrices
  snp_matrix <- as.matrix(snps_data %>% select(-IID))
  beta_vector <- as.numeric(bs$estimate)
  
  # Calculate the PRS for each participant using matrix multiplication
  prs_values <- snp_matrix %*% beta_vector
  working <- tibble(prs = prs_values[,1],
                    IID = all_snps$IID) %>%
    
                       rename_at(vars(prs), ~paste0("prs_", p_thresh))

  

  return(working)
  }

  map(p_thresh, generate_prs_internal) %>%
    reduce(., left_join, by = "IID") %>%
    janitor::clean_names() %>%
    rename(IID = iid) %>%
    relocate(IID)
  
}


