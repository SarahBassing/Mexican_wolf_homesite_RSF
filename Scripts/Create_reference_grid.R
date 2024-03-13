  #'  -----------------------------------------
  #'  Create reference grid for suitable WMEPA
  #'  Mexican wolf project
  #'  Sarah B. Bassing
  #'  December 2023
  #'  -----------------------------------------
  #'  Script to create a study area-wide grid. Read in empty raster that spans
  #'  the entire Mexican wolf experimental population area, mask out large water
  #'  bodies and habitat identified as "unsuitable" for Mexican gray wolves. Then
  #'  extract the centroid of each remaining pixel and save as a series of shapefiles
  #'  to be used to extract covariate data from numerous rasters and Google Earth Engine.
  #'  -----------------------------------------
  
  #'  Load libraries
  library(sf)
  library(terra)
  library(tidyverse)
  
  #'  ------------------------------------------
  ####  Load and format rasters and shapefiles  ####
  #'  ------------------------------------------
  #'  Load a raster with desired extent, resolution, and coordinate system
  elev <- terra::rast("./Shapefiles/Terrain_variables/Mosaic_DEM.tif"); res(elev); crs(elev)
  
  #'  Define WGS84 and NAD83 coordinate systems
  wgs84 <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
  nad83 <- st_crs(elev)
  
  #'  Load experimental population area and suitable habitat shapefiles and transform
  suitable_habitat <- st_read("./Shapefiles/Martinez_Meyer_2021_layers/suitable_habitat.shp") %>% st_transform(wgs84)
  wmepa <- st_read("./Shapefiles/MWEPA Layers Zone 1-3 & Boundary/MWEPA Final.shp"); crs(wmepa)
  wmepa_wgs84 <- st_transform(wmepa, wgs84)
  
  #'  Crop suitable habitat to experimental population area and reproject
  wmepa_suitable <- st_intersection(wmepa_wgs84, suitable_habitat)
  plot(wmepa_suitable[1])
  wmepa_suitable_nad83 <- st_transform(wmepa_suitable, crs = nad83)
  
  #'  Load NHD waterbodies shapefile
  waterbodies <- st_read("./Shapefiles/National Hydrography Dataset (NHD)/AZ_NM_waterbodies_NHD.shp"); crs(waterbodies)
  #'  Identify large bodies of water (anything larger than 1 sq-km in size)
  bigwater <- waterbodies[waterbodies$areasqkm > 1,]
  bigwater_nad83 <- st_transform(bigwater, crs = nad83)

  #'  -------------------------------------------------
  ####  Generate reference raster and centroid points  ####
  #'  -------------------------------------------------
  #####  Generate reference raster  #####
  #'  ------------------------------
  #'  Create empty raster based on extent, resolution, & coord system of elevation raster
  empty_rast <- rast(elev, nlyrs = nlyr(elev), names = "CellID", vals = 1:ncell(elev), keeptime = FALSE,
                     keepunits = FALSE, props = FALSE)
  #'  Make sure it checks out
  res(empty_rast); crs(empty_rast); st_bbox(empty_rast); st_bbox(elev)
  
  #'  Save that sucker
  writeRaster(empty_rast, filename  = "./Shapefiles/Reference_grid_30m.tif", overwrite = TRUE)
  
  #####  Mask out unsuitable areas to reduce pixels  #####
  #"  -----------------------------------------------
  #' #'  Load emtpy raster
  #' empty_rast <- terra::rast("./Shapefiles/Reference_grid_30m.tif")
  
  #'  Mask out waterbodies
  masked_water <- mask(empty_rast, bigwater_nad83, inverse = TRUE)
  plot(masked_water)
  
  #'  Mask out unsuitable areas (so areas NOT in the wmep_suitable polygon)
  wmepa_grid <- mask(masked_water, wmepa_suitable_nad83)
  plot(wmepa_grid); plot(homesites_nad83[1], add = TRUE)
  
  #'  Save!
  writeRaster(wmepa_grid, "./Shapefiles/WMEPA_masked_grid.tif", overwrite = TRUE)
  
  #####  Extract pixel centroids  #####
  #'  ----------------------------
  #' #'  Load masked reference grid
  #' wmepa_grid <- terra::rast("./Shapefiles/WMEPA_masked_grid.tif")
  
  #'  Convert grid to a polygon and sf object #### THIS TAKES FOREVER!!!
  wmepa_poly <- as.polygons(wmepa_grid)    
  wmepa_poly_sf <- st_as_sf(wmepa_poly)
  
  #'  Save so I never have to do that again!
  st_write(wmepa_poly_sf, "./Shapefiles/WMEPA_masked_polygon.shp")
  
  #'  Grab centroid of each pixel
  wmepa_grid_centroid <- st_centroid(wmepa_poly_sf)
  
  #'  Grab coordinates of pixel centroids
  wmepa_grid_coords <- st_coordinates(wmepa_grid_centroid)
  
  #'  Grab cellID data for each pixel 
  wmepa_grid_cellID <- wmepa_grid_centroid$CellID
  wmepa_grid_cellID <- as.data.frame(wmepa_grid_cellID)
  names(wmepa_grid_cellID) <- "cellID"
  
  #'  Merge cellID and centroid coordinates to create one dataframe 
  wmepa_grid_pts <- bind_cols(wmepa_grid_cellID, wmepa_grid_coords)
  
  #'  Save!
  write_csv(wmepa_grid_pts, file = "./Data/WMEPA_suitable_grid_points.csv")
  
  #'  Convert dataframe to shapefiles
  wmepa_grid_pts <- read_csv("./Data/WMEPA_suitable_grid_points.csv")
  
  #'  Convert to an sf object
  wmepa_grid_pts <- st_as_sf(wmepa_grid_pts, coords = c("X", "Y"), crs = nad83) %>%
    #'  Transform to WGS84 (needed for Google Earth Engine)
    st_transform(wgs84) %>%
    #'  Add a unique identifier that's ordered sequentially
    mutate(newID = 1:nrow(.))
  
  #'  Save the new suitable MWEPA grid as a shapefile
  st_write(wmepa_grid_pts, "./Shapefiles/MWEPA_suitable_reference_grid.shp")
  
  #'  Slit massive sf object into more manageable chunks (needed for Google Earth Engine)
  #'  Some of these are split into even smaller pieces b/c GEE had issues with the
  #'  larger datasets
  nrow(wmepa_grid_pts)
  wmepa_grid_pts1aa <- wmepa_grid_pts[1:250000,]
  wmepa_grid_pts1ab <- wmepa_grid_pts[250001:500000,]
  wmepa_grid_pts1b <- wmepa_grid_pts[500001:1000000,]
  wmepa_grid_pts2a <- wmepa_grid_pts[1000001:1500000,]
  wmepa_grid_pts2ba <- wmepa_grid_pts[1500001:1750000,]
  wmepa_grid_pts2bb <- wmepa_grid_pts[1750001:2000000,]
  wmepa_grid_pts3aa <- wmepa_grid_pts[2000001:2250000,]
  wmepa_grid_pts3ab <- wmepa_grid_pts[2250001:2500000,]
  wmepa_grid_pts3b <- wmepa_grid_pts[2500001:3000000,]
  wmepa_grid_pts4a <- wmepa_grid_pts[3000001:3500000,]
  wmepa_grid_pts4b <- wmepa_grid_pts[3500001:4000000,]
  wmepa_grid_pts5a <- wmepa_grid_pts[4000001:4500000,]
  wmepa_grid_pts5b <- wmepa_grid_pts[4500001:5000000,]
  wmepa_grid_pts6a <- wmepa_grid_pts[5000001:5500000,]
  wmepa_grid_pts6b <- wmepa_grid_pts[5500001:6000000,]
  wmepa_grid_pts7a <- wmepa_grid_pts[6000001:6500000,]
  wmepa_grid_pts7b <- wmepa_grid_pts[6500001:7000000,]
  wmepa_grid_pts8a <- wmepa_grid_pts[7000001:7500000,]
  wmepa_grid_pts8b <- wmepa_grid_pts[7500001:8000000,]
  wmepa_grid_pts9a <- wmepa_grid_pts[8000001:nrow(wmepa_grid_pts),]

  #'  Save individual shapefiles
  st_write(wmepa_grid_pts, "./Shapefiles/WMEPA_grid_points.shp")
  st_write(wmepa_grid_pts1aa, "./Shapefiles/WMEPA_grid_points1aa.shp")
  st_write(wmepa_grid_pts1ab, "./Shapefiles/WMEPA_grid_points1ab.shp")
  st_write(wmepa_grid_pts1b, "./Shapefiles/WMEPA_grid_points1b.shp")
  st_write(wmepa_grid_pts2a, "./Shapefiles/WMEPA_grid_points2a.shp")
  st_write(wmepa_grid_pts2ba, "./Shapefiles/WMEPA_grid_points2ba.shp")
  st_write(wmepa_grid_pts2bb, "./Shapefiles/WMEPA_grid_points2bb.shp")
  st_write(wmepa_grid_pts3aa, "./Shapefiles/WMEPA_grid_points3aa.shp")
  st_write(wmepa_grid_pts3ab, "./Shapefiles/WMEPA_grid_points3ab.shp")
  st_write(wmepa_grid_pts3b, "./Shapefiles/WMEPA_grid_points3b.shp")
  st_write(wmepa_grid_pts4a, "./Shapefiles/WMEPA_grid_points4a.shp")
  st_write(wmepa_grid_pts4b, "./Shapefiles/WMEPA_grid_points4b.shp")
  st_write(wmepa_grid_pts5a, "./Shapefiles/WMEPA_grid_points5a.shp")
  st_write(wmepa_grid_pts5b, "./Shapefiles/WMEPA_grid_points5b.shp")
  st_write(wmepa_grid_pts6a, "./Shapefiles/WMEPA_grid_points6a.shp")
  st_write(wmepa_grid_pts6b, "./Shapefiles/WMEPA_grid_points6b.shp")
  st_write(wmepa_grid_pts7a, "./Shapefiles/WMEPA_grid_points7a.shp")
  st_write(wmepa_grid_pts7b, "./Shapefiles/WMEPA_grid_points7b.shp")
  st_write(wmepa_grid_pts8a, "./Shapefiles/WMEPA_grid_points8a.shp")
  st_write(wmepa_grid_pts8b, "./Shapefiles/WMEPA_grid_points8b.shp")
  st_write(wmepa_grid_pts9a, "./Shapefiles/WMEPA_grid_points9a.shp")