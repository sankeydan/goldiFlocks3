#' @export
#' 
split_times = function ( data, pigeon, fission.matrices.num , gr.num, fl.num){
  
centroid = data.frame( x = apply(data[,"lon",], 1, mean), y = apply(data[,"lat",], 1, mean)) # centroid of whole group 

dist2cent = matrix(NA, nrow(data[,,1]), length(pigeon)) # empty matrix for distance to centroid

for ( j in 1:nrow(dist2cent)){ # for each timestep
  for ( k in 1:length(pigeon)){ # for each pigeon
    dist2cent[j,k] = get_dist( data[j,"lon",k], data[j, "lat", k], centroid$x[j], centroid$y[j] , method = "distance") # get the distance between the individual and the centroid
  }
}

fis.data = data
omit.list = list() # a list will store omitted individuals. 

for ( j in 1:nrow(dist2cent)){ # for each timestamp 
  
  vec = vector() # set up a vector to record ommitted individuals
  while( length( which(dist2cent[j,]  > fis.dist)) >0){ # while there is at least one individual over the fission distance threshold
    
    omit = which(dist2cent[j,]== max(dist2cent[j,], na.rm = T) ) # which individual needs to be removed?
    fis.data[j,c("lat","lon"), omit] = NA # make NA the latitude and the longitude of the furthest individual
    centroid$x[j] = mean(fis.data[j,"lon",], na.rm = T) # calcutate new centroid
    centroid$y[j] = mean(fis.data[j,"lat",], na.rm = T)
    
    for( k in 1:length(pigeon)){
      dist2cent[j,k] = get_dist( fis.data[j,"lon",k], fis.data[j, "lat", k], centroid$x[j], centroid$y[j] , method = "distance") # distance to centroid - same as above
    }
    
    vec = c(vec, omit) # build the vector of omitted individuals
  }
  omit.list[[j]] = vec # store in the list
  
}


length.omit = rep(NA, nrow(dist2cent)) # total omitted individuals per time stamp- using this object later for unbiased centroid speed

for ( j in 1:nrow(dist2cent)){  # for each timestamp 
  
  if ( length( omit.list[[j]] > 0 )) { # if there were any ommitted individuals 
    
    fission.matrices.num[[as.numeric(gr.num)]][ pigeon[omit.list[[j]]],fl.num] =    
      fission.matrices.num[[as.numeric(gr.num)]][ pigeon[omit.list[[j]]],fl.num] +1 # add one for each timestep
  }
  
  length.omit[j] = length( omit.list[[j]])
  
}



return(list(fission.matrices.num, nrow(dist2cent), length.omit, centroid))

}
