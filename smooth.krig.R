
#' Regularized kriging, optionally with externally supplied variogram
#'
#' Smoothes data using Gausssian kernel, then performs universal Krige interpolation
#'
#' @param latlon dataframe of latitutde and longitude for datapoints. Must have columns named "lon" and "lat".
#' @param dset dataframe of environmental variables
#' @param V name of the variable to interpolate
#' @param polygon dataframe of points (lon, lat) defining a polygon. The first row must be the same as the last row.
#' @param variogram_model fixed variogram model, made with variogram() and fit.variogram(). If omitted, will be fit on the fly.
#' @param grid.res grid resolution, in meters. 
#' @param smooth.range in meters. Set to some small value to turn off smoothing.
#' @param utm.zone UTM zone number
#' @param verbose whether to plot original data, smoothed data, and kriging result
#' @return dataframe of lat-lon (grid within polygon), and projected variable
#' @import sf
#' @import ggplot2
#' @import raster
#' @export
smooth.krig=function(latlon,dset,V,polygon,variogram_model=NULL,grid.res=1000,smooth.range=5000,utm.zone=17,verbose=TRUE){
#  latlon=ll;dset=dataset;V="NOX.B_range";variogram_model=variogram_model;polygon=keys.polygon;grid.res=1000;smooth.range=1;utm.zone=17;verbose=TRUE
  require(spatstat)
  require(automap)
  require(sp)
  require(sf)
  require(raster)
  require(gstat)
  
  # reformatting polygon and making grid  
  polygon <- st_polygon(list(polygon)) %>%
    st_sfc(crs = 4326) %>%
    st_sf()
  
  # Transform polygon to UTM Zone 17N (EPSG:32617)
  polygon_utm <- st_transform(polygon, crs = 32617)
  
  # Create a raster grid within the polygon's extent
  r <- raster(
    extent(st_bbox(polygon_utm)),
    res = grid.res,  # Resolution in meters
    crs = CRS(paste0("+proj=utm +zone=",utm.zone," +datum=WGS84 +units=m +no_defs"))
  )

  # Convert raster to SpatialPixelsDataFrame for interpolation
  grid <- as(r, "SpatialPixelsDataFrame")

  data <- na.omit(data.frame(cbind(x=latlon$lon,y=latlon$lat,z=dset[,V])))
  names(data)=c("x","y","z")
  # Convert to SpatialPointsDataFrame
  coordinates(data) <- ~x + y
  proj4string(data) <- CRS("+proj=longlat +datum=WGS84")  # WGS 84 CRS
  
  # Transform point data to UTM Zone 17N
  data_utm <- spTransform(data, CRS(paste0("+proj=utm +zone=",utm.zone," +datum=WGS84 +units=m +no_defs")))
#  if(!is.null(smooth.range)) {
          # Create a point pattern (ppp) object
    pp <- ppp(
      x = coordinates(data_utm)[, 1],
      y = coordinates(data_utm)[, 2],
      window = owin(
        xrange = range(coordinates(data_utm)[, 1]),
        yrange = range(coordinates(data_utm)[, 2])
      ),
      marks = data_utm$z
    )
    data_utm_sf <- st_as_sf(data_utm)
    # if (verbose==TRUE) { 
    #   plot(data_utm_sf["z"], pch=16,main="original data")
    # }
    
    # Apply kernel smoothing
    smoothed <- Smooth(pp, sigma = smooth.range,at="points",leaveoneout=FALSE)
#  }   else { smoothed=data_utm }
  # # Convert smoothed raster to points
  #   plot(smoothed)
  # smoothed_raster <- raster(smoothed)
  # smoothed_values <- raster::extract(smoothed_raster, coordinates(data_utm))
  # data_utm$z <- smoothed_values
  data_utm$z <- smoothed
  data_utm$z[is.na(data_utm$z)]=mean(data_utm$z,na.rm=T)
  data_utm_sf <- st_as_sf(data_utm)
#  if (verbose==TRUE) { plot(data_utm_sf["z"], pch=16,main="smoothed data")}
  if(is.null(variogram_model)) { 
    kriging_result <- autoKrige(z ~ x+y, input_data = data_utm, new_data = grid)
    kriging_raster <- raster(kriging_result$krige_output)
  } else {
    kriging_result <- krige(z ~ x+y, locations = data_utm, newdata = grid, model = variogram_model)
    kriging_raster <- raster(kriging_result)
  }
  
  # Extract interpolated raster
  
  # Mask raster to polygon
  kriging_clipped <- mask(kriging_raster, as(polygon_utm, "Spatial"))
  
  # Plot the result
  if (verbose==TRUE) { 
    plot(kriging_clipped, main = paste(V,smooth.range)) 
#    par(mfrow=c(1,1))
  }
  
  # Convert raster to SpatialPointsDataFrame
   rpoints <- rasterToPoints(kriging_clipped, spatial = TRUE)

  # Transform SpatialPointsDataFrame to WGS84 (Lat/Lon)
  rpoints_wgs84 <- spTransform(rpoints, CRS("+proj=longlat +datum=WGS84"))
  
  # Convert to data frame
  rdf <- as.data.frame(rpoints_wgs84)

  # Rename columns
  colnames(rdf) <- c(V,"lon","lat")
  return(rdf)
}

