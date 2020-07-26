create_time <- function(year, years, status){
  start <- list()
  stop <- list()
  event <- list()
  index <- 1
  
  for (i in 1:(length(year) - 1)) {
    if ((year[i] == 0 && year[i+1] != 0) || (year[i] != 0 && year[i+1] != 0)) {
      start[index] <- year[i]
      stop[index] <- year[i+1]
      
      if (status[i] == 1) {
        event[index] <- 0
      }
      else {
        event[index] <- status[i]
      }
    }
    
    else if ((year[i] == 0 && year[i+1] == 0) || (year[i] != 0 && year[i+1] == 0)) {
      start[index] <- year[i]
      stop[index] <- years[i]
      event[index] <- status[i]
    }
    index <- index + 1
  }
  
  start[index] <- year[i+1]
  stop[index] <- years[i+1]
  event[index] <- status[i+1]
  return(list(start, stop, event))
}
