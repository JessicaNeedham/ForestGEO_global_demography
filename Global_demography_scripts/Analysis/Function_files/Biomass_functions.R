################## 
# BIOMASS FUNCTIONS 
###################
# Get Basal Area in m2 
get.ba <- function(dbh){
  dbh <- dbh/1000  # get it in m
  ba <- pi*((dbh/2)^2)
  return(ba)
}

# temperate equation for bioground
J.agb.model <- function(dbh, param1, param2){
  agb.est <- exp(param1 + param2*log(dbh))
  return(agb.est) 
}

# tropical equation for bioground
agb.est <- function(D, E, wsg){
  agb.est <- exp(-1.803 - 0.976*E + 0.976*log(wsg) 
                 + 2.673*log(D) - 0.0299*(log(D)^2))
  return(agb.est)
}

make.transp <- function(tmpcol, transp = 50) {
  col.out <- paste(tmpcol, transp, sep = "")
  return(col.out)
}