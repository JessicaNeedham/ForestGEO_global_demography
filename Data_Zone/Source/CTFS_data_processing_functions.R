#######################################################################################
### Data processing functions ###

#### correct dbh measurements for all censuses at once
correct.dbh <- function(dbhs,
                        times = 5, 
                        sda = 0.927, sdb = 0.0038, 
                        give.up = 20000){
  
 
  # need to correct everything negative AND positive since
  # just correcting negative introduces a bias
  # so to start everything gets jittered according to the error distribution
  # things that are still negative get repeatedly jittered until positive
  neg <- seq(nrow(dbhs))
  dim.col <- ncol(dbhs)
  
  dbh.stars <- matrix(NA, nrow = nrow(dbhs), ncol = dim.col)
  diffs <- rep(NA, nrow(dbhs))
  
  ct <- 0
  # while anything is still negative keep going
  # until give.up to try and correct them
  while(length(neg)>1 & ct < give.up) {
    
    # propose new sizes for negative stems/increments
    suppressWarnings(dbh.stars[neg, ] <- 
                       matrix(rnorm(length(dbhs[neg ,]), dbhs[neg ,], 
                                    (sda + sdb * dbhs[neg ,])), 
                              ncol = dim.col, byrow = FALSE))
    
    # new increments for the new sizes
    diffs[neg] <- t(apply(dbh.stars[neg, ], 1, diff))
    # which ones are still negative
    neg.incrs <- which(diffs <= 0)
    neg.stems <- which(apply(dbh.stars, 1, function(x) any(x <= 0)) == TRUE)
    neg <- unique(c(neg.incrs, neg.stems))
    ct <- ct + 1
    
    if(ct %% 1000 == 0){
      print(sprintf('%d: Negative stems/increments: %d', ct, length(neg)))
    }
  }
  # at give up... if there are still negative stems/increments
  if(length(neg) > 0) {  
    
    # first correct negative sizes
    dbh.stars[which(dbh.stars < 0)] <- rnorm(length(which(dbh.stars < 0)), 10, 0.5)
    
    # now get new diffs
    diffs <- apply(dbh.stars, 1, diff)
    
    # which rows still have negative growths
    neg <- which(diffs <= 0)
    
    if(length(neg) > 0){
      diffs[neg] <- sample(diffs[diffs > 0 & !is.na(diffs)], length(neg), replace = TRUE)
      # adjust the sizes to take account of the new increments
      dbh.stars[ ,2] <- dbh.stars[ ,1] + diffs
      
    }
  } 
  n.sampled <- length(neg)  # how many were not corrected before give.up
  return(list(dbh.stars = dbh.stars, n.sampled = n.sampled))
}


