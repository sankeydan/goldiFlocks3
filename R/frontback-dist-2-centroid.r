##### Distance front / back

############


fb_lr_d2c = function ( data, 
                       fis.dist = 50,
                       hz = 5,
                       plot = T){
  
  # OBJECTS 
  
  # fis.dist = 50
  # data = data
  # hz = 5
  # plot = T

  # centroid of whole group regardless of fission
  centroid = data.frame( x = apply(data[,"lon",], 1, mean), y = apply(data[,"lat",], 1, mean)) 
  
  # distance to centroid
  dist2cent = matrix(NA, dim(data)[1], dim(data)[3] ) 
  for ( j in 1:nrow(dist2cent)){ # for each timestep
    for ( k in 1:dim(data)[3]){ # for each pigeon
      dist2cent[j,k] = get_dist( data[j,"lon",k], data[j, "lat", k], centroid$x[j], centroid$y[j] , method = "distance") # get the distance between the individual and the centroid
    }
  }
  
  # Fission 
  fis.data = data
  omit.list = list() # a list will store omitted individuals. 
  for ( j in 1:nrow(dist2cent)){ # for each timestamp 
    vec = vector() # set up a vector to record ommitted individuals
    while( length( which(dist2cent[j,]  > fis.dist)) >0){ # while there is at least one individual over the fission distance threshold
      omit = which(dist2cent[j,]== max(dist2cent[j,], na.rm = T) ) # which individual needs to be removed?
      fis.data[j,c("lat","lon"), omit] = NA # make NA the latitude and the longitude of the furthest individual
      centroid$x[j] = mean(fis.data[j,"lon",], na.rm = T) # calcutate new centroid
      centroid$y[j] = mean(fis.data[j,"lat",], na.rm = T)
      for( k in 1:dim(data)[3]){
        dist2cent[j,k] = get_dist( fis.data[j,"lon",k], fis.data[j, "lat", k], centroid$x[j], centroid$y[j] , method = "distance") # distance to centroid - same as above
      }
      vec = c(vec, omit) # build the vector of omitted individuals
    }
    omit.list[[j]] = vec # store in the list
  }
  
  # remove fissioned individuals from the group
  for ( i in 1:dim(data)[1]){
    data[i,1:2,][,omit.list[[i]]] = NA
  }
  
  # centroid of whole group with  fission
  centroid = data.frame( x = apply(data[,"lon",], 1,function( x){ mean(x,na.rm=T)}),
                         y = apply(data[,"lat",], 1,function( x){ mean(x,na.rm=T)})) 
 
  #  centroid speed and heading
  centroid = cbind(centroid,
                   head = get_heading( lon1 =  centroid[,"x"], lat1 = centroid[,"y"], indivs = 1),
                   speed =   get_dist( lon1 =  centroid[,"x"], lat1 = centroid[,"y"], method = "speed" , hz = hz))

  # frontback / leftright dist 2 centroid
  frontback.res.cent = array( NA, c(dim(data)[1] , dim(data)[3]), dimnames = list(NULL,names(data[1,1,]) ))
  leftright.res.cent = array( NA, c(dim(data)[1] , dim(data)[3]), dimnames = list(NULL,names(data[1,1,]) ))
  
  for ( j in 1:dim(data)[3]){
    ta <- centroid[,"x"]
    tb <- data[,1, j]
    dl <- data[,2, j] - centroid[, "y" ]
    X <- cos(tb)*sin(dl)
    Y <- (cos(ta)*sin(tb))-(sin(ta)*cos(tb)*cos(dl))
    focal.pos <- atan2(Y,X)
    for ( k in 1:dim(data)[1]){
      if ( !is.na(focal.pos[k])){
    focal.rel.pos.centroid <- ifelse (focal.pos[k]-centroid[k,"head"]<=-pi,
                                      focal.pos[k]+(2*pi)-centroid[k,"head"],
                                      ifelse(focal.pos[k]-centroid[k,"head"]>=pi,
                                             focal.pos[k]-(2*pi)-centroid[k,"head"],
                                             focal.pos[k]-centroid[k,"head"]))
    frontback.res.cent[k,j] <- cos(focal.rel.pos.centroid)*dist2cent[k,j]
    leftright.res.cent[k,j] <- sin(focal.rel.pos.centroid)*dist2cent[k,j]
      }
    }
  }
  
  #plot
  if ( plot){
  boxplot( frontback.res.cent)
  boxplot( leftright.res.cent)
  }
  
  # prepare for return
  frontback.res.cent = as.data.frame(frontback.res.cent)
  leftright.res.cent = as.data.frame(leftright.res.cent)

  
  return( list(frontback.res.cent, leftright.res.cent))
}

