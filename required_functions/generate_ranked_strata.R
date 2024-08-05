# generate_ranked_strata <- function(g, x, q, seed = 123){
#   # haodong ranked strata method
#   z = rank(g, ties.method = "random")
#   strata1 = floor((z-1)/q)+1
#   
#   set.seed(seed)
#   id = seq(x)
#   temp<- data.frame(x=x,strata1=strata1,id=id, g=g)
#   temp<- arrange(.data=temp, x)
#   temp<- group_by(.data=temp, strata1)
#   temp<- mutate(.data=temp, x0q= rank(x, ties.method = "random"))
#   temp<- arrange(.data=temp, id)
#   x0q <- temp$x0q
#   return(x0q)
# }



generate_ranked_strata <- function(instrument, exposure, k, seed = 310723){
  # haodong ranked strata method
  set.seed(seed)
  ###
  wdata = tibble(instrument = instrument, exposure = exposure) 
  ###
  wdata = wdata %>% mutate(z = rank(instrument, ties.method = "random")) %>%
    mutate(strata1 = floor((z-1)/k)+1) %>% 
    mutate(id = seq(exposure)) %>%
    arrange(exposure) %>%
    group_by(strata1) %>%
    mutate(x0q = rank(exposure, ties.method = "random") ) %>%
    arrange(id)
  ###
  bins = wdata %>% pull(x0q)
  
  return( as.factor(bins) )
}





