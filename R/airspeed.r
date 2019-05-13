#' Airspeed
#'
#' Methods rewritten from: Safi et al. 2013 
#' https://movementecologyjournal.biomedcentral.com/articles/10.1186/2051-3933-1-4
#' 
#' 
#' @export
#' 
airspeed = function ( ground.speed, heading, wind.dir , wind.speed ){
  
  x = abs(atan2(sin(wind.dir - heading), cos(wind.dir - heading)))
  
  if ( abs(x) > pi/2){
    
    x = abs(x) - (pi/2)
    cross = cos(x) * wind.speed
    support  = -abs(sin(x)) * wind.speed
    
  } else{
    cross = sin(x) * wind.speed
    support  = abs(cos(x)) * wind.speed
    
  }
  return( sqrt ( (ground.speed - support)^2 + cross^2))
}
