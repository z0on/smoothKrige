rm(list=ls())
library()

ll=load("smooth.krig_testData.RData")
ll
# [1] "keys.polygon" "latlon.pre.s" "s.pre"  

# keys.polygon: simple data frame of polygon corners. 
# generate it using locator() by clicking on a map (see appendix below)
keys.polygon
# lon      lat
# [1,] -82.95499 24.82131
# [2,] -83.15715 24.64552
# [3,] -83.01652 24.44337
# [4,] -82.71768 24.30274
# [5,] -81.80359 24.35547
# [6,] -80.81919 24.64552
# [7,] -80.26546 24.99709
# [8,] -80.07210 25.32230
# [9,] -80.01936 25.62992
# [10,] -80.13362 25.68266
# [11,] -80.39730 25.14651
# [12,] -80.76645 24.89162
# [13,] -81.02134 24.99709
# [14,] -81.88270 24.72462
# [15,] -82.95499 24.82131

# s.pre : SERC surface variables, 1995-2008
names(s.pre)
# [1] "DEPTH_median"   "NOX.S_median"   "NO3_S_median"   "NO2.S_median"  
# [5] "NH4.S_median"   "TN.S_median"    "DIN.S_median"   "TON.S_median"  
# [9] "TP.S_median"    "SRP.S_median"   "CHLA.S_median"  "TOC.S_median"  
# [13] "SiO2.S_median"  "TURB.S_median"  "SAL.S_median"   "TEMP.S_median" 
# [17] "DO.S_median"    "TN.TP_median"   "N.P_median"     "DIN.TP_median" 
# [21] "Si.DIN_median"  "X.SAT.S_median" "DSIGT_median"   "DEPTH_range"   
# [25] "NOX.S_range"    "NO3_S_range"    "NO2.S_range"    "NH4.S_range"   
# [29] "TN.S_range"     "DIN.S_range"    "TON.S_range"    "TP.S_range"    
# [33] "SRP.S_range"    "CHLA.S_range"   "TOC.S_range"    "SiO2.S_range"  
# [37] "TURB.S_range"   "SAL.S_range"    "TEMP.S_range"   "DO.S_range"    
# [41] "TN.TP_range"    "N.P_range"      "DIN.TP_range"   "Si.DIN_range"  
# [45] "X.SAT.S_range"  "DSIGT_range" 

# latlon.pre.s : lat, lon of measurements in s.pre

# simple map of our area
plot(keys.polygon,type="l",asp=1)
maps::map(add=T)

# variable to interpolate (change 8 to something else, too)
V=names(s.pre)[7]

# here we go
source("smooth.krig.R")
sk=smooth.krig(latlon=latlon.pre.s,dset=s.pre,V=V,keys.polygon,
            grid.res=1000,smooth.range=5000)
grid.latlon=sk[,1:2]
save(grid.latlon,file="grid.latlon.nov2024.RData")
# also try with smooth.range=10 to see how it looks without smoothing

#----------- "documentaion"-------

#' Smoothed kriging
#'
#' smoothes data, then performs krige interpolation using automap::autoKrige
#'
#' @param latlon dataframe of latitutde and longitude for datapoints. Must have columns named "lon" and "lat".
#' @param dset dataframe of environmental variables
#' @param V name of the variable to interpolate
#' @param polygon dataframe of points (lon, lat) defining a polygon. The first row must be the same as the last row.
#' @param variogram_model fixed variogram model, made with variogram() and fit.variogram(). If omitted, will be fit on the fly.
#' @param grid.res grid resolution, in meters. 
#' @param smooth.range in meters
#' @param utm.zone UTM zone number
#' @param verbose whether to plot original data, smoothed data, and kriging result

#-----appendix 1: drawing polygon using locator()

# margin around our study area (in degrees lat)
marg=0.2
ll=latlon.pre.s
plot(lat~lon,ll,pch=16,cex=0.5,col="red",
     xlim=c(min(ll$lon)-marg,max(ll$lon)+marg),
     ylim=c(min(ll$lat)-marg,max(ll$lat)+marg),
     asp=1
)
map(add=T, interior=F)

# drawing polygon using locator function
# click around the area, going around in sequence; 
# the last click must be near the first click
keys1=locator(n=18,type="l")

# closing the polygon
keys.polygon=cbind("lon"=keys1[["x"]],"lat"=keys1[["y"]])
keys.polygon[nrow(keys.polygon),]=keys.polygon[1,]

# adding the polygon to the figue
polygon(keys.polygon,border="cyan")

#----- appendix 2: defining fixed variogram model based on satellite data

# getting PC1 and PC2 scores for satellite data
load("satellite.RData")
rr=rda(sats~1,scale=T)
# plot(rr,scaling=1,display="sites")
# plot(rr$CA$eig/sum(rr$CA$eig))
scc=data.frame(scores(rr,scaling=1,display='sites'))

# making a colorful map plot based on pcs (skip if just need variogram)
b= scc$PC2
g = scc$PC1-scc$PC2
r = scc$PC1
r <- (r - min(r))/(max(r) - min(r)) * 255
g <- (g - min(g))/(max(g) - min(g)) * 255
b <- (b - min(b))/(max(b) - min(b)) * 255
colors=rgb(r,g,b,max=255)
plot(XY,pch=16,cex=0.3,col=colors,asp=1)
# plotting PC1
ggplot(cbind(XY,PC1=scc$PC1),aes(lon,lat,color=PC1))+
  geom_point(cex=0.03)+scale_color_viridis()+coord_equal()+theme_minimal()


# converting data to UTM coordinates
variogram_data <- data.frame(cbind(XY,z=scc$PC1))
coordinates(variogram_data) <- ~lon + lat
proj4string(variogram_data) <- CRS("+proj=longlat +datum=WGS84")
vdata_utm <- spTransform(variogram_data, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
dd=extent(vdata_utm)
# max distance from corner to corner
maxd=dist(rbind(c(dd@xmin,dd@ymin),c(dd@xmax,dd@ymax)))
maxd
# 350224

# empirical variogram
variogram_empirical = variogram(z ~ 1, data = vdata_utm,cutoff = 30000, width = 500)
plot(variogram_empirical)

# fitting a Matern model. It is quite robust to vgm() parameter settings.
variogram_model <- fit.variogram(variogram_empirical, model = vgm(psill = 1, model = "Mat", range = 10000, nugget = 10,kappa=1))
plot(variogram_empirical, variogram_model, main = "Fitted Variogram")
save(variogram_empirical, variogram_model,file="variogram_model_satPC1.RData")

# the object variogram_model can be passed to smooth.krig() as argument variogram_model, 
# to force it not to fit its own model but use the supplied one.
