#' plot_map
#'
#' A handy and interactive plotter for solo or group GPS tracks. Uses ggmaps to produce images or videos of animal trajectories on a googlemaps background.
#' @param data data frame or list of data frames with 2 dimensions: longitude & latitude (in that order; columns) For group plots, ensure data frames in the list have the same number of rows
#' @param type Solo or Group plots.
#' @param image.or.vid "image" or "video" (sequence of images) to be produced
#' @param tail.size animals represented by current position and the previous n positions (video only)
#' @param frame.size plot every row of data, or every 1/frame size for shorter running speeds and output files (video only)
#' @param zoom zoom = 19 is good for video. Other plot types are worth playing around with, try 13 for approx 5 km square
#' @param lines.or.points On the plotted map, animal trajectory represented by either "lines" or "points"
#' @param maptype "satellite" is default, though check ?get_map in the ggmap package for more options
#' @param wd Specify the working directory that you would like to save the plots into.
#' @param col Colour of plotted trajectory (solo only)
#' @param lwd Thickness of line
#' @param n.indiv Number of individuals. Only use in group context.
#' @param centre Choose your own centre points for the mapform: c(lon, lat)
#' @param file.to.folder TRUE if you want to make videos when happy with other params, but for testing videos set to FALSE
#' @export



plot_map = function( data,  type = c("solo", "group") , image.or.vid = "image",
                     tail.size = 20, frame.size = 1, wd = getwd(), plot_name = "my_plot" ,
                     lines.or.points = "lines", col = "red" , axis.nums.size = 12, axis.labs.size = 14, n.indiv = NULL,
                     lwd = 2, file.to.folder = T, centre = NA , maptype = "satellite", zoom = 19 ){


  if ( type[1] == "solo"){

    if ( class(data) == "data.frame"){
      dat = as.data.frame(data)
      names(dat) = c("lon", "lat")
      num.of.files = 1
      lis = F
    }

    if ( class(data) == "list"){
      num.of.files = length(data)
      dat = as.data.frame(data[[1]])
      names(dat) = c("lon", "lat")
      lis = T
    }

    if ( image.or.vid == "image"){

      if ( is.na(centre[1]) ){ # if no centre is specified.
        centre = c ( (max( dat[ ,1] ) + min( dat [ ,1]) ) /2 , (max( dat[ ,2] ) + min( dat [ ,2]) ) /2 ) # find centre for map
      }


      map = ggmap::get_map( location = centre,
                            zoom = zoom,
                            maptype = maptype) # load map
      p = ggmap::ggmap(map)

      for ( i in 1:num.of.files){

        if ( lis == T){
          dat = as.data.frame(data[[i]])
          names(dat) = c("lon", "lat")
        }
        p = p + ggplot2::geom_path (data= dat, ggplot2::aes(x=lon, y=lat),
                                    color = col[i], size = lwd)
      }

      scale_x_longitude <- function(xmin=-180, xmax=180, step=1, ...) {
        xbreaks <- seq(xmin,xmax,step)
        xlabels <- unlist(lapply(xbreaks, function(x) ifelse(x < 0, parse(text=paste0(x,"o", "*W")), ifelse(x > 0, parse(text=paste0(x,"^o", "*E")),x))))
        return(scale_x_continuous("Longitude", breaks = xbreaks, labels = xlabels, expand = c(0, 0), ...))
      }

      scale_y_latitude <- function(ymin=-90, ymax=90, step=0.5, ...) {
        ybreaks <- seq(ymin,ymax,step)
        ylabels <- unlist(lapply(ybreaks, function(x) ifelse(x < 0, parse(text=paste0(x,"o", "*S")), ifelse(x > 0, parse(text=paste0(x,"^o", "*N")),x))))
        return(scale_y_continuous("Latitude", breaks = ybreaks, labels = ylabels, expand = c(0, 0), ...))
      }

      p = p  + ggplot2::xlab("Longitude")+
        ggplot2::ylab("Latitude") +
        ggplot2::theme(axis.text= ggplot2::element_text(size= axis.nums.size),
                       axis.title= ggplot2::element_text(size= axis.labs.size ,face="bold"))
      print(p)
    }

    if ( image.or.vid == "video"){

      for ( j in seq(tail.size,nrow(dat),frame.size)){ # from the start (restricted by tail size), to the end of the data set, by frame size

        dat2 = dat[(j-tail.size+1):j,]
        centroid = data.frame ( lon = mean(dat[j,1]),
                                lat = mean(dat[j,2]))# A more advanced centriod algorithm is currently being updated for potential publication, email me and I can send it over though if you're interested.
        names(dat2) = c("lon" , "lat")

        if( file.to.folder){
          png(file = file.path( wd , paste0( plot_name,  j , ".png")))
        }

        map = ggmap::get_map( location = c(centroid$lon,
                                           centroid$lat),
                              zoom = zoom,
                              maptype = maptype)

        p = ggmap::ggmap(map, extent = "device")  +
          ggplot2::theme(axis.line = element_blank(),
                axis.text  = element_blank(),
                axis.ticks = element_blank(),
                plot.margin = unit(c(0, 0, -1, -1), 'lines')) +
          xlab('') +
          ylab('')
        p = p + ggplot2::geom_path (data= dat2, aes(x=lon, y=lat), color= col, size = lwd)
        print(p)

        if(file.to.folder){
          dev.off()
        }
      }
    }
  }

  if ( type[1] == "group"){

    if ( image.or.vid == "video"){

      for ( j in seq(tail.size,nrow(data[[1]]),frame.size)){ # from the start (restricted by tail size), to the end of the data set, by frame size

        data2 = data.frame(ldply(data, data.frame),
                           rep(1:n.indiv, each =  length(dataG[[1]][,1])))
        names(data2) = c( "lon", "lat", "id") # give the columns names
        data3 = simplify2array(by(data2, data2$id, as.matrix))
        data4 = plyr::adply(data3[(j-tail.size+1):j,1:2,],3)

        centroid = data.frame ( lon = mean(data3[j,1,]),
                                lat = mean(data3[j,2,]))# A more advanced centriod algorithm is currently being updated for potential publication, email me and I can send it over though if you're interested.

        if( file.to.folder){
          png(file = file.path( wd , paste0( plot_name,  j , ".png")))
        }


        map = ggmap::get_map( location = c(centroid$lon,
                                           centroid$lat),
                              zoom = zoom,
                              maptype = maptype)

        p = ggmap::ggmap(map, extent = "device")  +
          ggplot2::theme(axis.line = element_blank(),
                axis.text  = element_blank(),
                axis.ticks = element_blank(),
                plot.margin = unit(c(0, 0, -1, -1), 'lines')) +
          xlab('') +
          ylab('')

        rain = rainbow(dim(data3)[3])

        for ( k in 1:dim(data3)[3]){
          sub = (((k-1)* tail.size)+1):(k * tail.size)
          p = p + ggplot2::geom_path (data= data4[sub,], aes(x=lon, y=lat), color= rain[k], size = 2)

        }

        print( paste ( j , "/" , nrow(data3) ))
        print(suppressWarnings(p))

        if(file.to.folder){
          dev.off()
        }
      }
    }
  }
}






