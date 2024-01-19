  #'  ----------------
  #'  Generate spatial covariates 
  #'  January 2024
  #'  ----------------
  
  #'  Load libraries
  # install.packages("remotes")
  # remotes::install_github("ailich/MultiscaleDEM")
  # library(MultiscaleDEM)
  library(sf)
  library(terra)
  library(tidyverse)
  
  #'  Load spatial data
  dem <- terra::rast("./Shapefiles/GEE/Mosaic_DEM.tif"); res(dem); st_crs(dem)
  
  #'  Generate terrain data from DEM
  slope <- terrain(dem, v = "slope", neighbors = 8, unit = "degrees")
  aspect <- terrain(dem, v = "aspect", neighbors = 8, unit = "degrees")
  
  # writeRaster(slope, "./Shapefiles/Terrain_variables/slope.tif")
  # writeRaster(aspect, "./Shapefiles/Terrain_variables/aspect.tif")
  
  
  ####  Vector Ruggedness Measure (VRM)  ####
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
  
  writeRaster(x.sum, "./Shapefiles/Terrain_variables/x.sum.tif")
  writeRaster(y.sum, "./Shapefiles/Terrain_variables/y.sum.tif")
  writeRaster(z.sum, "./Shapefiles/Terrain_variables/z.sum.tif")
  
  #'  4) Calculate the magnitude of the resultant vector for each pixel based on the
  #'  x, y, and z values within each pixel's neighborhood
  #'  |r| = sqrt((sum(x)^2) + (sum(y)^2) + (sum(z)^2))
  r <- sqrt((x.sum^2) + (y.sum^2) + (z.sum^2))
  writeRaster(r, "./Shapefiles/Terrain_variables/r.tif")
  
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
  
  ####  Surface curvatures  ####
  #'  Profile curvature: second derivative of elevation surface (slope of the slope)
  #'  i.e., direction of the maximum slope where (-) values indicate surface is 
  #'  upwardly convex and (+) values indicate surface is upwardly concave. Values
  #'  of 0 indicate flat surface.
  profcurv <- curvature(dem, type = "profile")
  
  #'  Total curvature: sigma of the profile and planform curvatures (planform is 
  #'  perpendicular to direction of maximum slope).
  meancurv <- terra::curvature(dem, type = "total")
  
  
