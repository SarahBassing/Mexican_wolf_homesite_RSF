  #'  -----------------------------
  #'  Generate/reformat covariates
  #'  Mexican wolf project
  #'  Sarah B. Bassing
  #'  January 2024
  #'  -----------------------------
  #'  Load DEM, habiat suitability, and GEE extracted seasonal NDVI and annual canopy 
  #'  cover data. Generate vector ruggedness measurement (VRM) from elevation data,
  #'  create shapefile of area considered "suitable habitat" by Martinez-Meyer et 
  #'  al. 2021, reformat NDVI data, and calculate annual canopy cover.
  #'  -----------------------------
  
  #'  Load libraries
  library(sf)
  library(terra)
  library(stars)
  library(spatialEco)
  library(tidyverse)
  
  #'  -------------------------
  ####  Load elevation raster  ####
  #'  -------------------------
  #'  Load spatial data
  dem <- terra::rast("./Shapefiles/Terrain_variables/Mosaic_DEM.tif"); res(dem); st_crs(dem)
  nad83 <- crs(dem)
  
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
  #'  Raster extent represents the approximated historical distribution of Mexican 
  #'  wolves based on climatically suitable areas 
  #'  niche.txt = Climatic suitability raster
  #'  human_population.txt = Human population density suitability raster (-1 to 1, 
  #'  where -1 is most unsuitable, reflecting highest human densities in the region)
  niche <- terra::rast("./Shapefiles/Martinez_Meyer_2021_layers/niche.txt")
  humans <- terra::rast("./Shapefiles/Martinez_Meyer_2021_layers/human_population.txt")
  #'  Review raster details and visualize
  niche
  plot(niche)
  plot(humans)
  
  #'  Convert rasters to polygons outlining spatial extent of suitable habitat features
  suitable_niche <- as.polygons(niche)
  suitable_niche_sf <- st_as_sf(suitable_niche)
  human_pop <- as.polygons(humans)
  human_pop_sf <- st_as_sf(human_pop)
  
  #'  Visualize
  plot(suitable_niche_sf)
  plot(human_pop_sf)
  #'  Focus on only area with low human densities
  plot(human_pop_sf[human_pop_sf$human_population >0.5,])
  low_human_pop <- human_pop_sf[human_pop_sf$human_population >0.5,]
  plot(low_human_pop)
  
  #'  Remove areas of high human density from suitable_niche polygon to create
  #'  a final "suitable habitat" polygon 
  suitable_low_ppl <- st_intersection(suitable_niche_sf, low_human_pop)
  plot(suitable_low_ppl[1])
  
  #'  Save
  st_write(suitable_niche_sf, "./Shapefiles/Martinez_Meyer_2021_layers/suitable_climatic.shp")
  st_write(low_human_pop, "./Shapefiles/Martinez_Meyer_2021_layers/low_human_density.shp")
  st_write(suitable_low_ppl, "./Shapefiles/Martinez_Meyer_2021_layers/suitable_habitat.shp")
  
  #'  ------------------------
  ####  Resample all rasters  ####
  #'  ------------------------
  elev <- terra::rast("./Shapefiles/Terrain_variables/Mosaic_DEM.tif"); res(elev); crs(elev); st_bbox(elev)
  slope <- terra::rast("./Shapefiles/Terrain_variables/slope.tif"); res(slope); crs(slope); st_bbox(slope)
  rough <- terra::rast("./Shapefiles/Terrain_variables/VRM.tif"); res(rough); crs(rough); st_bbox(rough)
  curve <- terra::rast("./Shapefiles/Terrain_variables/Gaussian_curvature.tif"); res(curve); crs(curve); st_bbox(curve)
  water <- terra::rast("./Shapefiles/National Hydrography Dataset (NHD)/Mosaic_Dist2Water.tif"); res(water); crs(water); st_bbox(water)
  human_mod <- terra::rast("./Shapefiles/Human_variables/mosaic_global_Human_Modification.tif"); res(human_mod); crs(human_mod); st_bbox(human_mod)
  roads <- terra::rast("./Shapefiles/Human_variables/mosaic_dist2road.tif"); res(roads); crs(roads); st_bbox(roads)
  # ext(water) <- ext(elev)
  # ext(human_mod) <- ext(elev)
  # ext(roads) <- ext(elev)
  # st_bbox(elev); st_bbox(slope); st_bbox(water); st_bbox(roads)
  # dim(elev); dim(slope); dim(water); dim(human_mod); dim(roads)
  # water_resampled <- resample(water, elev, method = "bilinear", filename = "./Shapefiles/National Hydrography Dataset (NHD)/Mosaic_Dist2Water_resampled.tif")
  # human_mod_resampled <- resample(human_mod, elev)
  # roads_resampled <- resample(roads, elev)
  # st_bbox(elev); st_bbox(slope); st_bbox(water_resampled); st_bbox(roads_resampled); st_bbox(human_mod_resampled)
  # dim(elev); dim(slope); dim(water_resampled); dim(human_mod_resampled); dim(roads_resampled)
  # 
  # writeRaster(water_resampled, "./Shapefiles/National Hydrography Dataset (NHD)/Mosaic_Dist2Water_resampled.tif")
  # writeRaster(human_mod_resampled, "./Shapefiles/Human_variables/mosaic_global_Human_Modification_resampled.tif")
  # writeRaster(roads_resampled, "./Shapefiles/Human_variables/mosaic_dist2road_resampled.tif")
  
  #'  -------------------------------------
  ####  Exported Google Earth Engine data  ####
  #'  -------------------------------------
  #'  Global forest change data set to estimate annual canopy cover (Hansen Global Forest Change v1.10 (2000-2022))
  #'  NDVI data set to estimate mean seasonal greenness (MOD13Q1.061 Terra Vegetation Indices 16-Day Global 250m (2000 - 2023))
  
  #####  Annual canopy cover  ##### 
  #'  ------------------------
  #'  Derived from Global Forest Change (2000 - 2022) data
  #'  Average 2000 canopy cover within 250 of each location (used & available)
  mean_canopy_den <- read_csv("./Data/GEE extracted data/GEE_meanCanopyCover_2000_den.csv") %>%
    dplyr::select(c(ID, Pack_year, Site_Type, mean, used)) %>%
    arrange(ID) %>%
    rename("Mean_canopy_cover_2000" = "mean")
  mean_canopy_rnd <- read_csv("./Data/GEE extracted data/GEE_meanCanopyCover_2000_rnd.csv") %>%
    dplyr::select(c(ID, Pack_year, Site_Type, mean, used)) %>%
    arrange(ID) %>%
    rename("Mean_canopy_cover_2000" = "mean")
  
  #'  Sum of area (sq. m) of canopy loss over time within 250m of each location
  canopy_loss_den <- read_csv("./Data/GEE extracted data/GEE_accumulated_canopy_loss_area_den.csv") %>%
    dplyr::select(c(ID, Pack_year, used, BufferArea_sq_m, CanopyLossArea_sq_m, Year_add1)) %>%
    group_by(Pack_year, Year_add1) %>%
    arrange(ID, .by_group = TRUE) %>%
    mutate(Year = Year_add1 + 2001,
           LossYear = paste0("lossYr_", Year)) %>%
    ungroup() %>%
    dplyr::select(-Year_add1) 
  canopy_loss_rnd <- read_csv("./Data/GEE extracted data/GEE_accumulated_canopy_loss_area_rnd.csv") %>%
    dplyr::select(c(ID, Pack_year, used, BufferArea_sq_m, CanopyLossArea_sq_m, Year_add1))  %>%
    group_by(Pack_year, Year_add1) %>%
    arrange(ID, .by_group = TRUE) %>%
    mutate(Year = Year_add1 + 2001,
           LossYear = paste0("lossYr_", Year)) %>%
    ungroup() %>%
    dplyr::select(-Year_add1) 
  
  #'  Brief look
  head(mean_canopy_den); head(canopy_loss_den)
  bufferedArea <- max(canopy_loss_den$BufferArea_sq_m)
  
  #'  Reformat and update canopy cover based on years when canopy loss occurred
  canopy_loss_den_reformat <- canopy_loss_den %>%
    #'  Each column represents area (sq. m) of canopy loss as it accumulates across years
    pivot_wider(!Year, names_from = LossYear, values_from = CanopyLossArea_sq_m) %>%
    #'  Join with mean canopy cover from 2000
    full_join(mean_canopy_den, by = c("ID", "Pack_year", "used")) %>%
    relocate(Site_Type, .after = Pack_year) %>%
    relocate(Mean_canopy_cover_2000, .after = used) %>%
    #'  Put mean % canopy cover in a real percentage
    mutate(Mean_2000_CC_percent = Mean_canopy_cover_2000 /100,
           #'  Calculate proportion of area within buffer that was lost over time
           across(lossYr_2001:lossYr_2022, function(x)(x/BufferArea_sq_m)),
           #'  Calculate percent of area in buffer lost as it accumulates over time
           across(lossYr_2001:lossYr_2022, function(x)(x * Mean_2000_CC_percent)),
           #'  Adjust mean % canopy cover from 2000 by percentage of area lost across years
           across(lossYr_2001:lossYr_2022, function(x)(Mean_2000_CC_percent - x)),
           #'  Remove all characters before year in pack_year
           site_year = as.numeric(gsub(".*_", " ", Pack_year)),
           #'  Change site_year if from before 2000 or after 2022 (need to apply
           #'  2000 or 2022 to sites outside this time window)
           site_year = ifelse(site_year < 2000, 2000, site_year),
           site_year = ifelse(site_year > 2022, 2022, site_year)) %>%
    rename("canopy_cov_2000" = "Mean_2000_CC_percent") %>%
    relocate(canopy_cov_2000, .after = Mean_canopy_cover_2000) %>%
    dplyr::select(-BufferArea_sq_m)

  canopy_loss_rnd_reformat <- canopy_loss_rnd %>%
    #'  Each column represents area (sq. m) of canopy loss as it accumulates across years
    pivot_wider(!Year, names_from = LossYear, values_from = CanopyLossArea_sq_m) %>%
    #'  Join with mean canopy cover from 2000
    full_join(mean_canopy_rnd, by = c("ID", "Pack_year", "used")) %>%
    relocate(Site_Type, .after = Pack_year) %>%
    relocate(Mean_canopy_cover_2000, .after = used) %>%
    #'  Put mean % canopy cover in a real percentage
    mutate(Mean_2000_CC_percent = Mean_canopy_cover_2000 /100,
           #'  Calculate proportion of area within buffer that was lost over time
           across(lossYr_2001:lossYr_2022, function(x)(x/BufferArea_sq_m)),
           #'  Calculate percent of area in buffer lost as it accumulates over time
           across(lossYr_2001:lossYr_2022, function(x)(x * Mean_2000_CC_percent)),
           #'  Adjust mean % canopy cover from 2000 by percentage of area lost across years
           across(lossYr_2001:lossYr_2022, function(x)(Mean_2000_CC_percent - x)),
           #'  Remove all characters before year in pack_year
           site_year = as.numeric(gsub(".*_", " ", Pack_year)),
           #'  Change site_year if from before 2000 or after 2022 (need to apply
           #'  2000 or 2022 to sites outside this time window)
           site_year = ifelse(site_year < 2000, 2000, site_year),
           site_year = ifelse(site_year > 2022, 2022, site_year)) %>%
    rename("canopy_cov_2000" = "Mean_2000_CC_percent") %>%
    relocate(canopy_cov_2000, .after = Mean_canopy_cover_2000) %>%
    dplyr::select(-BufferArea_sq_m)
  
  #'  Re-pivot wide datasets into long format and filter to only years where lossYr and site_year match
  canopy_cover_den <- canopy_loss_den_reformat %>%
    dplyr::select(-c(Site_Type, Mean_canopy_cover_2000)) %>%
    pivot_longer(!c(ID, Pack_year, site_year, used), names_to = "Canopy_year", values_to = "Mean_percent_canopy")
  
  canopy_cover_rnd <- canopy_loss_rnd_reformat %>%
    dplyr::select(-c(Site_Type, Mean_canopy_cover_2000)) %>%
    pivot_longer(!c(ID, Pack_year, site_year, used), names_to = "Canopy_year", values_to = "Mean_percent_canopy")
  
  #'  Function to relabel the canopy_year for each data set
  relabel_canopy_year <- function(dat) {
    relabeled_dat <- dat %>%
      mutate(Canopy_year = ifelse(Canopy_year == "canopy_cov_2000", "2000", Canopy_year),
             Canopy_year = ifelse(Canopy_year == "lossYr_2001", "2001", Canopy_year), 
             Canopy_year = ifelse(Canopy_year == "lossYr_2002", "2002", Canopy_year),
             Canopy_year = ifelse(Canopy_year == "lossYr_2003", "2003", Canopy_year),
             Canopy_year = ifelse(Canopy_year == "lossYr_2004", "2004", Canopy_year),
             Canopy_year = ifelse(Canopy_year == "lossYr_2005", "2005", Canopy_year),
             Canopy_year = ifelse(Canopy_year == "lossYr_2006", "2006", Canopy_year),
             Canopy_year = ifelse(Canopy_year == "lossYr_2007", "2007", Canopy_year),
             Canopy_year = ifelse(Canopy_year == "lossYr_2008", "2008", Canopy_year),
             Canopy_year = ifelse(Canopy_year == "lossYr_2009", "2009", Canopy_year),
             Canopy_year = ifelse(Canopy_year == "lossYr_2010", "2010", Canopy_year),
             Canopy_year = ifelse(Canopy_year == "lossYr_2011", "2011", Canopy_year),
             Canopy_year = ifelse(Canopy_year == "lossYr_2012", "2012", Canopy_year),
             Canopy_year = ifelse(Canopy_year == "lossYr_2013", "2013", Canopy_year),
             Canopy_year = ifelse(Canopy_year == "lossYr_2014", "2014", Canopy_year),
             Canopy_year = ifelse(Canopy_year == "lossYr_2015", "2015", Canopy_year),
             Canopy_year = ifelse(Canopy_year == "lossYr_2016", "2016", Canopy_year),
             Canopy_year = ifelse(Canopy_year == "lossYr_2017", "2017", Canopy_year),
             Canopy_year = ifelse(Canopy_year == "lossYr_2018", "2018", Canopy_year),
             Canopy_year = ifelse(Canopy_year == "lossYr_2019", "2019", Canopy_year),
             Canopy_year = ifelse(Canopy_year == "lossYr_2020", "2020", Canopy_year), 
             Canopy_year = ifelse(Canopy_year == "lossYr_2021", "2021", Canopy_year),
             Canopy_year = ifelse(Canopy_year == "lossYr_2022", "2022", Canopy_year))
    return(relabeled_dat)
  }
  canopy_cover_den <- relabel_canopy_year(canopy_cover_den) 
  canopy_cover_rnd <- relabel_canopy_year(canopy_cover_rnd) 
  
  #'  Annual canopy cover, averaged across entire study area
  annual_cover_den <- canopy_cover_den %>%
    filter(used == 0) %>%
    group_by(Canopy_year) %>%
    summarise(avg_MCP_canopycover = mean(Mean_percent_canopy)) %>%
    ungroup()
  annual_cover_rnd <- canopy_cover_rnd %>%
    filter(used == 0) %>%
    group_by(Canopy_year) %>%
    summarise(avg_MCP_canopycover = mean(Mean_percent_canopy)) %>%
    ungroup()
  
  #'  Bind site-level and study area average canopy cover data per site and year
  percent_canopy_den <- canopy_cover_den %>%
    filter(site_year == Canopy_year) %>%
    left_join(annual_cover_den, by = "Canopy_year")
  percent_canopy_rnd <- canopy_cover_rnd %>%
    filter(site_year == Canopy_year) %>%
    left_join(annual_cover_rnd, by = "Canopy_year")
    
  #'  Save
  write_csv(percent_canopy_den, "./Data/GEE extracted data/percent_canopy_denSeason.csv")
  write_csv(percent_canopy_rnd, "./Data/GEE extracted data/percent_canopy_rndSeason.csv")
  
  
  ######  2022 canopy cover MWEPA grid  ######
  #'  ----------------------------------
  #' #'  Merge pieces of 2000 canopy cover data for entire MWEPA suitable grid
  #' canopy_grid <-list.files(path = "./Data/GEE extracted data/MWEPA suitable grid 2000 canopy/", pattern = "\\.csv$", full.names = T) %>%
  #'   map_df(~read_csv(., col_types = cols(.default = "c"))) %>%
  #'   dplyr::select(-c(`system:index`,`.geo`))
  #' write_csv(canopy_grid, "./Data/GEE extracted data/GEE_2000_canopy_cover_grid.csv")
  #' 
  #' #'  Merge pieces of 2022 canopy loss data for entire MWEPA suitable grid
  #' loss_grid <-list.files(path = "./Data/GEE extracted data/MWEPA suitable grid 2022 canopy loss/", pattern = "\\.csv$", full.names = T) %>%
  #'   map_df(~read_csv(., col_types = cols(.default = "c"))) %>%
  #'   dplyr::select(-c(`system:index`,`.geo`))
  #' write_csv(loss_grid, "./Data/GEE extracted data/GEE_accu_canopy_loss_2022_grid.csv")
  
  #'  Load MWEPA-wide 2000 canopy cover and accumulated area lost through 2022
  canopy_cover_mwepa <- terra::rast("./Shapefiles/Vegetation_variables/Resampled_CanopyCover_2000_movingwindow.tif"); crs(canopy_cover_mwepa); st_bbox(canopy_cover_mwepa)
  acculoss_2022_mwepa <- terra::rast("./Shapefiles/Vegetation_variables/Resampled_lossAccu2022_reduceNeighborhood.tif"); crs(acculoss_2022_mwepa); st_bbox(acculoss_2022_mwepa)
  dim(acculoss_2022_mwepa); dim(canopy_cover_mwepa)
  
  #'  Load reference grid and extract canopy data at each point
  grid_pts <- st_read("./Shapefiles/WMEPA_grid_clip_pts.shp"); crs(grid_pts)
  xy <- st_coordinates(grid_pts)
  # grid_pts <- st_read("./Shapefiles/MWEPA_suitable_reference_grid.shp") %>%
  #   st_transform(crs = nad83); crs(grid_pts)
  canopy_grid <- terra::extract(canopy_cover_mwepa, grid_pts)
  loss_grid <- terra::extract(acculoss_2022_mwepa, grid_pts)
  
  head(canopy_grid); head(loss_grid)
  
  #'  Calculate 2022 canopy cover using 2000 baseline and accumulated loss area
  percent_canopy_2022_grid <- full_join(canopy_grid, loss_grid, by = "ID") %>% #by = c("cellID", "newID")
    mutate(BufferArea_sq_m = bufferedArea) %>%
    rename("lossYr_2022" = "Resampled_lossAccu2022_reduceNeighborhood") %>%
    rename("Mean_canopy_cover_2000" = "Resampled_CanopyCover_2000_movingwindow") %>%
    # rename("lossYr_2022" = "sum") %>% 
    # rename("Mean_canopy_cover_2000" = "mean") %>%
    # relocate(newID, .after = cellID) %>%
    #'  Make sure canopy cover and area loss are numeric
    mutate(Mean_canopy_cover_2000 = as.numeric(Mean_canopy_cover_2000),
           lossYr_2022 = as.numeric(lossYr_2022),
           lossYr_2022 = ifelse(is.na(lossYr_2022), 0, lossYr_2022),
           #'  Put mean % canopy cover in a real percentage
           Mean_2000_CC_percent = Mean_canopy_cover_2000/100,
           #'  Calculate proportion of area within buffer that was lost over time
           PropLost = lossYr_2022/BufferArea_sq_m,
           #'  Calculate accumulated percentage of area in buffer lost over time
           PercentLost = PropLost * Mean_2000_CC_percent,
           #'  Adjust mean % canopy cover from 2000 by percentage of area lost up to 2022
           Mean_percent_canopy = Mean_2000_CC_percent - PercentLost) %>%
    rename("canopy_cov_2000" = "Mean_2000_CC_percent") %>%
    relocate(canopy_cov_2000, .after = Mean_canopy_cover_2000) %>%
    dplyr::select(ID, Mean_percent_canopy)
  
  #'  2022 canopy cover, averaged across entire study area
  mean_cover_2022_grid <- mean(percent_canopy_2022_grid$Mean_percent_canopy) 
  
  #'  Add 2022 study area average canopy cover to 2022 canopy cover data set
  percent_canopy_grid <- percent_canopy_2022_grid %>%
    mutate(avg_MWEPA_canopycover = mean_cover_2022_grid) %>%
    cbind(xy)
  
  write_csv(percent_canopy_grid, "./Data/GEE extracted data/Resampled_percent_canopy_2022_grid.csv")
  
  #####  Mean Seasonal Greenness  #####
  #'  ---------------------------
  #'  Derived from MODIS Terra 16-day (2000 - 2023) data
  #'  Load and format seasonal mean NDVI data, averaged within 250 m radius of each location
  mean_denNDVI <- read_csv("./Data/GEE extracted data/GEE_annual_mean_NDVI_den.csv") %>%
    dplyr::select(-`.geo`) %>%
    #'  Remove all characters after year identifier in system:index column
    mutate(Year_ID = gsub("_.*", " ", `system:index`),
           #'  Turn Year_ID into an actual year (2000 is first year of NDVI data)
           NDVI_year = as.numeric(Year_ID) + 2000,
           #'  Remove all characters before year in pack_year
           site_year = as.numeric(gsub(".*_", " ", Pack_year)),
           #'  Change site_year if from before 2000 (need to apply 2000 data to 
           #'  sites from before this time window)
           site_year = ifelse(site_year < 2000, 2000, site_year)) %>%
    #'  Rearrange data by pack year so used location is listed first, then all available sites
    group_by(Pack_year, Year_ID) %>%
    arrange(ID, .by_group = TRUE) %>%
    ungroup() %>%
    rename("meanNDVI" = "mean") %>%
    arrange(Pack_year, NDVI_year)
  
  mean_rndNDVI <- read_csv("./Data/GEE extracted data/GEE_annual_mean_NDVI_rnd.csv") %>%
    dplyr::select(-`.geo`) %>%
    #'  Remove all characters after year identifier in system:index column
    mutate(Year_ID = gsub("_.*", " ", `system:index`),
           #'  Turn Year_ID into an actual year (2000 is first year of NDVI data)
           NDVI_year = as.numeric(Year_ID) + 2000,
           #'  Remove all characters before year in pack_year
           site_year = as.numeric(gsub(".*_", " ", Pack_year)),
           #'  Change site_year if from before 2000 (need to apply 2000 data to 
           #'  sites from before this time window)
           site_year = ifelse(site_year < 2000, 2000, site_year)) %>%
    #'  Rearrange data by pack year so used location is listed first, then all available sites
    group_by(Pack_year, Year_ID) %>%
    arrange(ID, .by_group = TRUE) %>%
    ungroup() %>%
    rename("meanNDVI" = "mean") %>%
    arrange(Pack_year, NDVI_year)
  
  #'  Reduce to only observations where the site year matches the meanNDVI year
  #'  This reflects what was used and available each year
  NDVI_denSeason <- mean_denNDVI %>% 
    filter(site_year == NDVI_year) %>%
    dplyr::select(c(ID, Pack_year, used, NDVI_year, meanNDVI))
      
  NDVI_rndSeason <- mean_rndNDVI %>%
    filter(site_year == NDVI_year) %>%
    dplyr::select(c(ID, Pack_year, used, NDVI_year, meanNDVI))
  
  #'  Format average NDVI based on what's available across entire study area
  avg_NDVI_denSeason <- read_csv("./Data/GEE extracted data/GEE_annual_mean_NDVI_MCPavg_den.csv") %>%
    dplyr::select(-`.geo`) %>%
    #'  Remove all characters after year identifier in system:index column
    mutate(Year_ID = gsub("_.*", " ", `system:index`),
           #'  Turn Year_ID into an actual year (2000 is first year of NDVI data)
           NDVI_year = as.numeric(Year_ID) + 2000) %>%
    dplyr::select(c(NDVI_year, mean)) %>%
    rename("avg_MCP_meanNDVI" = "mean")
  
  avg_NDVI_rndSeason <- read_csv("./Data/GEE extracted data/GEE_annual_mean_NDVI_MCPavg_rnd.csv") %>%
    dplyr::select(-`.geo`) %>%
    #'  Remove all characters after year identifier in system:index column
    mutate(Year_ID = gsub("_.*", " ", `system:index`),
           #'  Turn Year_ID into an actual year (2000 is first year of NDVI data)
           NDVI_year = as.numeric(Year_ID) + 2000) %>%
    dplyr::select(c(NDVI_year, mean)) %>%
    rename("avg_MCP_meanNDVI" = "mean")
  
  #'  Merge site-level NDVI with overall mean available NDVI
  ndvi_den <- left_join(NDVI_denSeason, avg_NDVI_denSeason, by = "NDVI_year") 
  ndvi_rnd <- left_join(NDVI_rndSeason, avg_NDVI_rndSeason, by = "NDVI_year")
  
  #'  Save
  write_csv(ndvi_den, "./Data/GEE extracted data/mean_NDVI_denSeason.csv")
  write_csv(ndvi_rnd, "./Data/GEE extracted data/mean_NDVI_rndSeason.csv")
  
  ######  2023 Rendezvous NDVI MWEPA grid  ######
  #'  -------------------------------------
  #'  Merge pieces of 2023 rendezvous NDVI data for entire MWEPA suitable grid
  ndvi_grid <-list.files(path = "./Data/GEE extracted data/MWEPA suitable grid NDVI/", pattern = "\\.csv$", full.names = T) %>%
    map_df(~read_csv(., col_types = cols(.default = "c"))) %>%
    dplyr::select(-c(`system:index`,`.geo`))
  write_csv(ndvi_grid, "./Data/GEE extracted data/GEE_mean_NDVI_grid_rnd_2023.csv")
  names(ndvi_grid) <- c("CellID", "meanNDVI", "newID")
  
  #'  Load reference grid & centroids
  ref_grid <- rast("./Shapefiles/WMEPA_masked_grid.tif"); dim(ref_grid)
  # ref_grid <- setValues(ref_grid, 1:ncell(ref_grid))
  grid_poly <- st_read("./Shapefiles/WMEPA_masked_polygon.shp"); crs(grid_poly)
  
  #'  Append NDVI to polygon sf object
  ndvi_poly <- full_join(grid_poly, ndvi_grid, by = "CellID"); crs(ndvi_poly)
  #'  Rename for st_rasterize
  names(ndvi_poly) <- c("cellID", "value", "newID", "geometry")
  
  #'  Rasterize NDVI data
  #'  Use MWEPA masked grid as the template for rasterizing so the resolution, extent, and coordinate system are correct
  meanNDVI_2023rnd_raster <- st_rasterize(ndvi_poly %>% dplyr::select(value, geometry), template = read_stars("./Shapefiles/WMEPA_masked_grid.tif"), align = TRUE)
  #'  Save
  write_stars(meanNDVI_2023rnd_raster, "./Shapefiles/meanNDVI_2023rnd_raster2.tif")
  
  ndvi_raster <- rast("./Shapefiles/meanNDVI_2023rnd_raster2.tif"); res(ndvi_raster); crs(ndvi_raster)
  plot(ndvi_raster)
  
  