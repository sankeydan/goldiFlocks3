runningmean = function ( x , lag){
  
  #x = solo.speed[5:nrow(solo.speed),1]
  #lag = 8
  len = length(x)
  sq = (lag+1):(len-lag)
  
  vec = rep(NA, len)
  for ( i in sq){
    vec[i] = mean(x[(i-lag):(i+lag)] , na.rm = T )  
  }
  for ( i in lag:1){
    
    # i = 2
    vec[i] = mean(x[((i-i)+1):(i+(i-1))],na.rm = T)
    vec[length(vec)-i] = mean( x[((len-i-i)+1):((len-i)+(i-1))],na.rm = T)
  }

return(   vec )
}
