#' @export
#' 

plot_and_stats = function (delta, mass, group.num , tarsus,  LMM.plot = T, qq.plot = F){

  g1 = which(group.num == 1)
  g2 = which(group.num == 2)
  
  # PLOT
  
  if ( qq.plot ){
    par(mfrow = c(2,1))
  }

  # LMM
  delta_speed = data.frame( group.num, delta , mass)
  m.mixed = lmer(delta ~ mass + (1 | group.num ), data = delta_speed)
  coefs <- data.frame(coef(summary(m.mixed)))
  # use normal distribution to approximate p-value
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))


  if ( LMM.plot){
    p = ggplot( delta_speed, aes( x = mass, y = delta)  )
    p = p + geom_point(aes(color=as.factor(group.num))) +
      geom_smooth(aes(color=as.factor(group.num)),
                  method=lm, se=FALSE, fullrange=TRUE) +
       geom_abline(intercept = 0, slope = 0)
    print(p)
  } else {
    plot( delta ~ mass)
    co = summary(lm(delta~mass))$coeff
    abline( co[1,1], co[2,1])
    
    plot( lm(delta~ mass) , which = 2)
  }


  # LM
  mass.mod = summary(lm(delta~ mass))
  tars.mod = summary(lm(delta~ tarsus))

  group2.mod = summary(lm(delta[g2]~ mass[g2]))
  group1.mod = summary(lm(delta[g1]~ mass[g1]))

  return ( list( list(paste("Linear model for only group1"), group1.mod),
                 list(paste("linear model for only group2"), group2.mod),
                 list( paste("Linear model, delta over mass"), mass.mod) ,
                 list(paste("Linear model, delta over tarsus length"), tars.mod) ,
                 list(paste("mixed model with group as random"), paste("P value = "), coefs$p.z[2], paste("summary "), summary(m.mixed))))


}
