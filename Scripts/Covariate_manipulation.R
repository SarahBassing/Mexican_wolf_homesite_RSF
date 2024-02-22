  #'  -----------------------------
  #'  Generate spatial covariates 
  #'  January 2024
  #'  -----------------------------
  #'  Load DEM, seasonal NDVI, and annual forest cover. Generate new terrain variables
  #'  and format NDVI/forest data for easier extraction.
  #'  -----------------------------
  
  #'  Load libraries
  library(sf)
  library(terra)
  library(spatialEco)
  library(tidyverse)
  
  #'  --------------------
  ####  Load raster data  ####
  #'  --------------------
  #'  Load spatial data
  dem <- terra::rast("./Shapefiles/Terrain_variables/Mosaic_DEM.tif"); res(dem); st_crs(dem)
  
  #####  Generate terrain data from DEM  #####
  #'  -----------------------------------
  slope <- terrain(dem, v = "slope", neighbors = 8, unit = "degrees")
  aspect <- terrain(dem, v = "aspect", neighbors = 8, unit = "degrees")
  
  res(slope); crs(slope)
  res(aspect); crs(aspect)
  
  # writeRaster(slope, "./Shapefiles/Terrain_variables/slope.tif")
  # writeRaster(aspect, "./Shapefiles/Terrain_variables/aspect.tif")
  
  #####  Vector Ruggedness Measure (VRM)  #####
  #'  -------------------------------------
  #'  Using VRM to represent surface roughness
  #'  1) Decompose DEM into x, y, & z components using trigonometric operators 
  z <- 1 * cos(slope)
  xy <- 1 * sin(slope)
  x <- xy * sin(aspect)
  y <- xy * cos(aspect)
  
  #'  2) Define neighborhood for moving-window analysis
  neighborhood <- 3
  
  #'  3) Sum x, y, and z values within each focal pixel's neighborhood
  x.sum <- focal(x, neighborhood, sum) 
  y.sum <- focal(y, neighborhood, sum) 
  z.sum <- focal(z, neighborhood, sum) 
  
  writeRaster(x.sum, "./Shapefiles/Terrain_variables/VRM parts/x.sum.tif")
  writeRaster(y.sum, "./Shapefiles/Terrain_variables/VRM parts/y.sum.tif")
  writeRaster(z.sum, "./Shapefiles/Terrain_variables/VRM parts/z.sum.tif")
  
  #'  4) Calculate the magnitude of the resultant vector for each pixel based on the
  #'  x, y, and z values within each pixel's neighborhood
  #'  |r| = sqrt((sum(x)^2) + (sum(y)^2) + (sum(z)^2))
  r <- sqrt((x.sum^2) + (y.sum^2) + (z.sum^2))
  writeRaster(r, "./Shapefiles/Terrain_variables/VRM parts/r.tif")
  
  #'  5) Calculate standardized ruggedness value per pixel
  #'  ruggedness = 1 - (|r|/n) where n = number of pixels used to estimate |r|
  #'  Values can range 0 (flat) - 1 (most rugged)
  ncells <- neighborhood * neighborhood
  vrm <- 1 - (r/ncells)
  writeRaster(vrm, "./Shapefiles/Terrain_variables/VRM.tif", overwrite = TRUE)
  
  #'  Review vrm attributes (min, max, mean, and sd of pixel values)
  minmax(vrm)
  app(vrm, mean)
  app(vrm, sd)
  plot(vrm)
  
  #'  -------------------------------------------------------
  ####  Martinez-Meyer et al. 2021 Habitat suitability data  ####
  #'  -------------------------------------------------------
  #'  Load data downloaded from https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/MartinezMeyer_et_al_DivDist_2021
  #'  niche.txt = Climatic suitability raster
  #'  human_population.txt = Human population density suitability raster
  #'  road_density.txt = Road density suitability raster
  niche <- terra::rast("./Shapefiles/Martinez_Meyer_2021_layers/niche.txt")
  #'  Review raster details and visualize
  niche
  plot(niche)
  
  #'  Convert raster to polygon outlining spatial extent of suitable habitat
  suitable_habitat <- as.polygons(niche)
  suitable_habitat_sf <- st_as_sf(suitable_habitat)
  plot(suitable_habitat_sf)
  st_write(suitable_habitat_sf, "./Shapefiles/Martinez_Meyer_2021_layers/suitable_habitat.shp")
  
  
  
  
  
  
  #' #####  Surface curvatures  #####
  #' #'  ------------------------
  #' #'  Gaussian curvature and Profile curvature
  #' #'  Code adapted from spatialEco curvature function (https://rdrr.io/cran/spatialEco/src/R/curvature.R)
  #' #'  and curvature equations defined in Minar et al. 2020
  #' my_curves <- function(x, type = c("planform", "profile", "total", "gaussian", "mcnab", "bolstad"), ...) { 
  #'   if(!inherits(x, "SpatRaster"))
  #'     stop(deparse(substitute(x)), " must be a terra SpatRaster object")
  #'   m <- matrix(1, nrow=3, ncol=3)
  #'   type = type[1] 
  #'   if(!any(c("planform", "profile", "total", "gaussian", "mcnab", "bolstad") %in% type)  )
  #'     stop("Not a valid curvature option")	
  #'   zt.crv <- function(m, method = type, res = terra::res(x)[1]) {
  #'     p=(m[6]-m[4])/(2*res)
  #'     q=(m[2]-m[8])/(2*res)
  #'     r=(m[4]+m[6]-2*m[5])/(2*(res^2))
  #'     s=(m[3]+m[7]-m[1]-m[9])/(4*(res^2))
  #'     tx=(m[2]+m[8]-2*m[5])/(2*(res^2))
  #'     if(type == "planform") {
  #'       crv <- round( -(q^2*r-2*p*q*s+p^2*tx)/((p^2+q^2)*sqrt(1+p^2+q^2)),6) 
  #'     } else if(type == "profile") {
  #'       crv <- round( -(p^2*r+2*p*q*s+q^2*tx)/((p^2+q^2)*sqrt(1+p^2+q^2)^3),6 )
  #'     } else if(type == "total") {
  #'       crv <- round( -(q^2*r-2*p*q*s+p^2*tx)/((p^2+q^2)*sqrt(1+p^2+q^2)),6) + 
  #'         round( -(p^2*r+2*p*q*s+q^2*tx)/((p^2+q^2)*sqrt(1+p^2+q^2)^3),6 )  
  #'     } else if(type == "gaussian"){
  #'       crv <- round( -(r*tx-s^2)/((1+p^2+q^2)^2),6 )
  #'     }
  #'     return(crv)
  #'   }  	
  #'   if(type == "bolstad") {
  #'     return( 10000 * ((x - terra::focal(x, w=m, fun=mean, ...)) / 1000 /
  #'                        (terra::res[1] + terra::res[1]/2) ) )
  #'   } else if(type == "mcnab") {
  #'       mcnab <- function(x) {
  #'       m <- ceiling(length(x)/2)
  #'       sum(x[m] - x[-m], na.rm=TRUE) / length(x)-1
  #'     }
  #'     tidx <- mask(terra::focal(x, w=m, fun=mcnab, ...), x)
  #'     return( tidx / as.numeric(terra::global(tidx, "max", na.rm=TRUE)) )
  #'   }
  #'   else {
  #'     return( terra::focal(x, w=m, fun = zt.crv, fillvalue = 0, ...) )
  #'   }
  #' }	
  #' gaus_curv <- my_curves(dem, type = "gaussian")
  #' plot(gaus_curv)
  #' writeRaster(gaus_curv, "./Shapefiles/Terrain_variables/Gaussian_curvature.tif", overwrite = TRUE)
  #' 
  #' # Converting Minar parameters to spatialEco parameters 
  #' # p = x
  #' # q = y
  #' # r = xx
  #' # s = xy
  #' # tx = yy
  #' 
  #' #'  Profile curvature: second derivative of elevation surface (slope of the slope)
  #' #'  i.e., direction of the maximum slope where (-) values indicate surface is 
  #' #'  upwardly convex and (+) values indicate surface is upwardly concave. Values
  #' #'  of 0 indicate flat surface.
  #' profcurv <- curvature(dem, type = "profile")
  #' writeRaster(profcurv, "./Shapefiles/Terrain_variables/Profile_curvature.tif", overwrite = TRUE)
  
  
  
  
  
  
  
