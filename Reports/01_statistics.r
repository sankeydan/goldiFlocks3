# HOUSEKEEPING

rm(list=ls()) # clear workspace


# LIBRARIES
{
  library(ggplot2)
  library(lme4)
  library(nlme)
  library(Goldilocks)
  library(car) 
  library(broom)
}

# PARAMS
{
  models = c("anova_slope", 
             "anova_intercept",
             "lm_mass_speed",
             "lm_massresid_speedcomp",
             "x2_poly_workrate_vs_mass",
             "lmm_mass_speed",
             "anchormod",
             "linearmod")
  statistics = c("R2","T","pval")
  params = c("standard.methods", "use.cent", "rm.cooks","omit.first.6","wind.support","med4mean")
  base = c(F,T,T,F,F,F)
  stopifnot(length(base)==length(params))
  arr = array( NA, c(length(models), length(statistics), length(base)),
               dimnames = list( models,statistics,params))
  arr2 = array( NA,c( 3,4,length(params)), dimnames = list( c("intercept", "x.minus","x.plus"),
                                                            c("Estimate", "Std. Error", "t value", "Pr(>|t|)"), params   ))
}

#### LOOOP ####
for ( i in 1:length(base)){
  #i=2
  
  
  # params
  base2   = base
  if ( i != 0){
    base2[i]= !base[i] 
  }
  for ( j in 1:length(params)){
    assign(params[j], base2[j])
  }
  
  # objects
  group.vec = c(7,9) # This is how many individuals were group 1 and 2 respecively
  group.num = rep(1:2, group.vec)
  g1 = which(group.num == 1)
  g2 = which(group.num == 2) 
  
  ### DATA
  # group speed
  if ( med4mean){
    load( file.path(PROJHOME, "Output" , "cent-speed-array-median.rda" ))
  } else {
    load( file.path(PROJHOME, "Output" , "cent-speed-array.rda" ))
  }
  # solo speeds and other metrics
  if ( med4mean){
    load(file.path(PROJHOME , "Output" , "metrics-midpoint-median.rda"))
    metrics = metrics_midpoint
  } else {
    load(file.path(PROJHOME , "Output" , "metrics-midpoint.rda"))
    metrics = metrics_midpoint
  }
  metrics.all = metrics
  rm(metrics_midpoint)
  # morphometrics
  load( file.path(PROJHOME , "Output" , "morpho-finishers.rda"))
  mass = morpho$mean.mass
  # Release site coordinates
  load(file.path ( PROJHOME , "Output" , "release-site-coordinates.rda"))
  # wind support
  if ( wind.support){
    load( file.path( PROJHOME , "Output", "wind-support.rda"))
  }
  # groupspeed
  load( file.path (PROJHOME , "Output", "observed-speeds.rda"))
  group.speed= observed_speeds[[2]]
  
  ## 2. TRANSFORM SPEEDS FOR ANALYSIS
  solo.speed = matrix(NA, ncol = length(unique(metrics$pigeon)),
                      nrow = 30 )
  
  for (j in 1:length(unique(metrics$pigeon))){ 
    #j = 1
    solo.speed[,j] = metrics$air.speed[metrics$pigeon == unique(metrics$pigeon)[j] & metrics$g_s == "Sol" ]
  }
  #group speed
  if ( use.cent){
    stat.group = rep(c(median(cent.speed.array[,1,"air.speed"]),
                       median(cent.speed.array[,2,"air.speed"])),group.vec)
  } else{ 
    stat.group =rep(c(median( as.vector(group.speed[,g1]),na.rm = T),
                      median( as.vector(group.speed[,g2]),na.rm = T)),group.vec)
  }
  # individual speed 
  if ( omit.first.6 ){
    stat.ind = apply(solo.speed[-c(1:6, 13:17),], 2, median, na.rm = T)
  } else {
    stat.ind = apply(solo.speed, 2, median, na.rm = T)
  }
  
  # wind support? 
  if ( wind.support){
    stat.ind   = stat.ind   - wind_support$sol
    stat.group = stat.group - wind_support$gro
  }
  
  # delta - individual minus group speed
  delta = ( stat.ind - stat.group)
  delta_m.s = delta
  delta = ((delta)/stat.group )*100
  
  # mass residual 
  h.l = mass- c(rep(mean(mass[g1]),group.vec[1]), rep(mean(mass[g2]), group.vec[2]))
  sub = 1:length(h.l)
  mass.resid = abs(h.l)
  
  ## 3. STATISTICAL ANALYSIS
  
  ## Cooks distance
  if( rm.cooks ){
    
    # Which pigeon speeds are above threshold cooks distance? 
    cooks.dist =  c( cooks.distance(lm( delta[g1] ~ mass[g1])) ,
                     cooks.distance(lm( delta[g2] ~ mass[g2]))) 
    # {plot(cooks.dist ~ c(1:nrow(morpho)), xlab = "Individual")
    #   abline( h = 3 * mean(cooks.dist) , col  = "red", lty = 2)
    #   abline( h = 1 , lty = 2 )
    #   legend ( "topleft" , legend = c( "Threshold = 1" , "Threshold = 3 * mean"), lty = 2, col = c(1,2))
    # }
    cook = which(cooks.dist >( 3 * mean(cooks.dist) ))
    
    if( length(cook) >= 1){
      # remove pigeons above cooks dist
      metrics = metrics[ metrics$pigeon != morpho$Pigeon[cook],] 
      morpho = morpho[ -cook,]
      mass = morpho$mean.mass
      
      # modify objects 
      group.vec = c( length(which(morpho$group.num == 1)) , length(which(morpho$group.num == 2))) # This is how many individuals were group 1 and 2 respecively
      group.num = rep(1:2 , group.vec)
      g1 = which(group.num == 1)
      g2 = which(group.num == 2)
      delta = delta[-cook]
      mass = mass[-cook]
      mass.resid = mass.resid[-cook]
      stat.ind = stat.ind[-cook]
      h.l = h.l[-cook]
    }
  } else {
    allind = list ( delta, h.l)
    save( allind , file = file.path (PROJHOME, "Output", "allind-rmcooks.rda"))
  }
  plot(abs(delta) ~ h.l)
  if ( standard.methods){
    y = abs(delta)
    x = h.l
    xy = list(y,x)
    save ( xy , file =file.path (PROJHOME, "Output", "kneemod.rda"))
  }
  # Statistics
  
  ### ANCOVA to test for significant differences between slopes
  a1 = summary(aov( lm ( stat.ind ~ mass * group.num)))
  a3 = unlist(a1,recursive = F)
  a2 = aov( lm ( stat.ind ~ mass * group.num))
  tidy_aov <- tidy(a2)
  sum_squares_regression <- tidy_aov$sumsq[1]
  sum_squares_residuals <- tidy_aov$sumsq[2]
  R_squared <- sum_squares_regression /
    (sum_squares_regression + sum_squares_residuals)
  
  arr[ "anova_slope","pval",params[i]] = round( a3[[5]][3],3)
  arr[ "anova_slope","T"   ,params[i]] =   round(a3$`F value`[3],3)
  arr[ "anova_slope","R2",params[i]] =   round( R_squared , 3)
  
  ### ANOVA test for differences between intercepts
  a1 = Anova( lm( stat.ind ~ mass + as.character(group.num)))
  a3 = unlist(a1,recursive = F)
  a2 = aov( lm ( stat.ind ~ mass  + as.character(group.num)))
  tidy_aov <- tidy(a2)
  sum_squares_regression <- tidy_aov$sumsq[1]
  sum_squares_residuals <- tidy_aov$sumsq[2]
  R_squared <- sum_squares_regression /
    (sum_squares_regression + sum_squares_residuals)
  arr[ "anova_intercept","pval",params[i]] = round( a3['Pr(>F)2'] , 3)
  arr[ "anova_intercept","T"   ,params[i]] =   round( a1$`F value`[2], 3)
  arr[ "anova_intercept","R2",params[i]] = round(R_squared,3)
  
  # LINEAR MODEL - statistics
  lm1 = summary(lm(delta~ mass))
  arr["lm_mass_speed","R2"  ,params[i]] = round( lm1$r.squared,3)
  arr["lm_mass_speed","T"   ,params[i]] = round( lm1$coeff[2,3],3)
  arr["lm_mass_speed","pval",params[i]] = round( lm1$coeff[2,4],3)
  
  # LM2 delta
  lm2 = summary(lm(abs(delta)~mass.resid))
  arr["lm_massresid_speedcomp","R2"  ,params[i]] = round( lm2$r.squared,3)
  arr["lm_massresid_speedcomp","T"   ,params[i]] = round( lm2$coeff[2,3],3)
  arr["lm_massresid_speedcomp","pval",params[i]] = round( lm2$coeff[2,4],3)
  
  
  ## 6. THEORETICAL DELTA POWER DIFFERENCES 
  load( file.path(PROJHOME , "Output", "data-theo-power.rda"))
  dtp = data.theo.power
  if ( rm.cooks){
    if ( length(cook >= 1)){
      dtp = dtp[-16,]
    }
  }
  y = dtp$power
  x=dtp$mass-mean(dtp$mass)
  mod <- lm(y ~ 1 + x + I(x^2))
  m3 = summary(mod)
  arr["x2_poly_workrate_vs_mass" ,"R2",params[i]] = round( m3$r.squared,3)
  arr["x2_poly_workrate_vs_mass" ,"T",params[i]] = round( m3$coeff[3,3],3)
  arr["x2_poly_workrate_vs_mass" ,"pval",params[i]] = round( m3$coeff[3,4],3)
  
  
  ##  LARGER MODEL TO INCLUDE ALL FLIGHTS 
  
  # just solo metrics 
  metrics.sol = metrics[ metrics$g_s =="Sol",]#need support wind
  x = abs(atan2(sin(metrics.sol$wind.dir.circ - metrics.sol$av.head), cos(metrics.sol$wind.dir.circ - metrics.sol$av.head)))
  metrics.sol$cross = NA
  metrics.sol$support=NA
  for ( j in 1:length(x)){
    #j=71
    if( !is.na(x[j])){
      if ( abs(x[j]) > pi/2){
        x[j] = abs(x[j]) - (pi/2)
        metrics.sol$cross[j] = cos(x[j] ) * metrics.sol$wind.speed[j] 
        metrics.sol$support[j]   = -abs(sin(x[j] )) * metrics.sol$wind.speed[j]
      } else{
        metrics.sol$cross[j]  = sin(x[j] ) * metrics.sol$wind.speed[j] 
        metrics.sol$support[j]   = abs(cos(x[j] )) * metrics.sol$wind.speed[j] 
      }
    }
  }
  # Speed 
  stat.group = c(median(cent.speed.array[,1,"air.speed"]),
                 median(cent.speed.array[,2,"air.speed"]))
  del = ifelse ( metrics.sol$group.num == 1,((metrics.sol$air.speed - stat.group[1])/stat.group[1])*100,((metrics.sol$air.speed-stat.group[2])/stat.group[2])*100 )
  metrics.sol$abs.rel.speed = abs(del)
  metrics.sol$rel.speed = del
  # Mass
  h.l = mass-
    c(rep(mean(mass[g1]),group.vec[1]), rep(mean(mass[g2]), group.vec[2]))
  mass.resid = abs(h.l)
  metrics.sol$mass.resid = NA
  metrics.sol$mass = NA
  for ( j in 1:nrow(morpho) ){
    #j=1
    metrics.sol$mass.resid [which( metrics.sol$pigeon == morpho$Pigeon[j])] = mass.resid[j]
    metrics.sol$mass       [which( metrics.sol$pigeon == morpho$Pigeon[j])] = mass[j]
  }
  # last flights count as flights 13,14,15 not 1,2,3
  len.pid = length(unique(metrics$pigeon))
  last = (24*len.pid+1)
  metrics.sol$flight[last:nrow(metrics.sol)] = rep(rep(13:15,each = len.pid),2)
  # pigeon as random variable, so 
  metrics.sol$pigeon = as.character(metrics.sol$pigeon)
  
  ##### MODELS 
  m.mixed1 = lmer(ground.speed~ mass + cross +support+ flight + (1 | group.num ) + (1|pigeon), data = metrics.sol)
  coefs <- data.frame(coef(summary(m.mixed1)))
  cf = confint(m.mixed1)
  cf2 = round( cf["mass",],3)
  arr["lmm_mass_speed", "pval",params[i]] = paste0( c("[", cf2[1], " ", cf2[2] , "]"),collapse = "")
  arr["lmm_mass_speed", "T",params[i]] = round( coefs["mass","t.value"],3)
  arr["lmm_mass_speed", "R2",params[i]] = round( MuMIn::r.squaredGLMM(m.mixed1)[2],3)
  
  ## elbow/anchor mod
  y = abs(delta)
  x = h.l
  x.0 = 0
  x.plus <- ifelse( x> 0, x,0)
  X <- data.frame(x=x, x.plus=x.plus)
  mod = lm(formula = y ~ x + x.plus, data = X)
  sum = summary(mod)
  arr2[,,i] = sum$coefficients[,]
  library(broom)
  gl = glance(mod)
  
  arr["anchormod","R2",params[i]] =   round(sum$r.squared    ,3) 
  arr["anchormod","T",params[i]] =    round(sum$fstatistic[1],3)
  arr["anchormod","pval",params[i]] = round(gl$p.value       ,3)
  
  
  #### SEX
  sex = c("M", "M", "M" ,"F" ,"F" ,"M" ,"F" ,"M" ,"M" ,"F" ,"F" ,"F" ,"F" ,"M" ,"F" ,"M")
  if ( rm.cooks){
    sex = sex[-cook]
  } 
  x.save = x
  y.save = y
  x = x[sex=="M"]
  y = y[sex=="M"]
  x.plus <- ifelse( x> 0, x,0)
  X <- data.frame(x=x, x.plus=x.plus)
  mod = lm(formula = y ~ x + x.plus, data = X)
  sum = summary(mod)
  x = x.save
  y = y.save
  x = x[sex=="F"]
  y = y[sex=="F"]
  x.plus <- ifelse( x> 0, x,0)
  X <- data.frame(x=x, x.plus=x.plus)
  mod = lm(formula = y ~ x + x.plus, data = X)
  sum = summary(mod)
  
  
  #### linear mod
  
  y = abs(delta)
  x = h.l
  sum = summary(lm(y~x))
  
  arr["linearmod","R2",params[i]]   = round(sum$r.squared        ,3)
  arr["linearmod","T",params[i]]    = round(sum$coefficients[2,3],3)
  arr["linearmod","pval",params[i]] = round(sum$coefficients[2,4],3)
  
  # END OF SCRIPT
  # TAKE STOCK
  print(i)
  
}

# save
write.csv(arr ,file = file.path(PROJHOME ,"Output", "paramsensitivitytables.csv"))
write.csv(arr2,file = file.path(PROJHOME ,"Output", "anchormod-output.csv"))

