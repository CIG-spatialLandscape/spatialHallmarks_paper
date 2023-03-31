##################################################
## Project: Cancer Hallmarks
## Script purpose: Function to compute the real distance between sub-spots (utils)
## Author: Sergi Cervilla* & Mustafa Sibai*
##################################################

#Function to compute real distances at sub-spot level
subspot_coord <- function(coord)  {
  for (subspot in rownames(coord)) {
    realrow <-  coord[paste0(substr(subspot, 0, nchar(subspot)-2), ".5"), "row"]
    realrow <- (realrow+1)*100
    realcol <- trunc(coord[paste0(substr(subspot, 0, nchar(subspot)-2), ".3"), "col"])
    realcol <- (realcol+1)*100
    n <- substr(subspot, nchar(subspot), nchar(subspot))
    
    factor <- 55/2/3
    if (n == 1) {
      realrow <- realrow + factor
      realcol <- realcol + factor
    }  else if (n == 2) {
      realrow <- realrow + factor
      realcol <- realcol - factor
    }  else if (n == 3) {
      realrow <- realrow - factor
      realcol <- realcol + factor
    }  else if (n == 4) {
      realrow <- realrow - factor
      realcol <- realcol - factor
    }  else if (n == 5) {
      realrow <- realrow
      realcol <- realcol + 2 * factor
    }  else if (n == 6) {
      realrow <- realrow
      realcol <- realcol - 2 * factor
    }
    coord[subspot, c("realrow", "realcol")] <- c(realrow, realcol)
  }
  return(coord)
}

