#' metrics_row
#'
#'   Which row in the metricws dataframe corresponds to our data?
#' @export
#' 
metrics_row = function( file , input = c("GPS" , "ACC" , "ufi&id"), 
                        output = c("row_num" , "pigeon", "group.num", "f_u", "site", "g_s", "dir", "miles",
                                          "flight", "accuracy", "mean.speed", "var.speed", "total.dist",
                                          "time.taken", "time.stationary")) {


  load( file.path( PROJHOME, "Data" , "metrics.rda"))

  vec = rep(NA, length( file))
  for ( metrics_row_i in 1:length(file)){

    v1 = str_split_fixed( file[metrics_row_i] , "\\.", 16)

    if ( input[1] == "GPS"){
      pigeon = as.numeric(v1[,11])
    }

    if ( input[1] == "ACC"){
      pigeon= as.numeric(v1[,12])
    }

    if ( input[1] == "ufi&id"){
      pigeon = as.numeric(v1[,10])
    }

    w.1 = which(metrics$pigeon == pigeon)
    w.2 = which(metrics$g_s == v1[,4])
    w.3 = which(metrics$site == v1[,6])
    w.4 = which(metrics$flight == v1[,9])

    if ( output[1] == "row_num"){
      vec[metrics_row_i] = Reduce(intersect, list(w.1,w.2,w.3,w.4)) # this will be a vector of length 1. The row which you will output new metrics.
    }

    if ( output[1] != "row_num"){
      vec[metrics_row_i] = eval(parse(text = paste0( "metrics$" , output , "[" ,
                                        Reduce(intersect, list(w.1,w.2,w.3,w.4,w.5)) , "]")))
    }
  }

  return(vec)

}

