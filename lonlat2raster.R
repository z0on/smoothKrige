
#' Create a multilayer SpatRaster out of lonlat and values dataframes
#'
#' ..to extract data with coordinates from it: as.data.frame(RASTERname,xy=TRUE)
#' 
#' 
#' @param lonlat dataframe of latitutde and longitude corresponding to rows in `values`. Must have columns named "lon" and "lat".
#' @param values dataframe of environmental variables
#' @param res desired resolution, in degrees latitude
#' @param polygon dataframe of points (lon, lat) defining a polygon. The first row must be the same as the last row.
#' @param variogram_model fixed variogram model, made with variogram() and fit.variogram(). If omitted, will be fit on the fly.
#' @param grid.res grid resolution, in meters. 
#' @param smooth.range in meters
#' @param utm.zone UTM zone number
#' @param verbose whether to plot original data, smoothed data, and kriging result
#' @return dataframe of lat-lon (grid within polygon), and projected variable
#' @import terra
#' @export
lonlat2raster=function(lonlat,values,res=0.01){
  dff=data.frame(cbind(lonlat,values))
  lon_range <- range(lonlat$lon, na.rm = TRUE)
  lat_range <- range(lonlat$lat, na.rm = TRUE)
  
  # Create an empty raster template
  template_raster <- rast(
    xmin = lon_range[1], xmax = lon_range[2],
    ymin = lat_range[1], ymax = lat_range[2],
    resolution = res,
    crs = "EPSG:4326"  # WGS84
  )
  
  # Convert dataframe to SpatVector
  data_vect <- vect(dff, geom = c("lon", "lat"), crs = "EPSG:4326")
  
  # Rasterize each layer
  layers <- lapply(names(values), function(layer_name) {
    rasterize(data_vect, template_raster, field = layer_name)
  })
  
  # Combine all layers into a multi-layer raster
  mrast <- rast(layers)
  
  # Assign names to layers
  names(mrast) <- names(values)
  
  return(mrast)
}

