  #'  ------------------------
  #'  Data prepartion
  #'  Mexican wolf project
  #'  Sarah B. Bassing
  #'  December 2023
  #'  ------------------------
  #'  Script to visualize and explore wolf homesite location data, randomly sample
  #'  "available" points and extract covariate data for each "used" and "available"
  #'  location. Final data set to be used as input for RSF analyses.
  #'  ------------------------
  
  #'  Clear memory
  rm(list = ls())
  
  #'  Load libraries
  library(readr)
  library(sf)
  library(terra)
  library(adehabitatHR)
  library(tidyverse)
  
  #'  Load used locations
  homesites <- read_csv("./Data/MexWolf_dens_rend_sites_1998_2023.csv") %>%
    mutate(Pack_year = paste0(Pack, "_", Year)) %>%
    relocate(Pack_year, .before = Year)
  
  #'  Define projection
  wgs84 <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
  #' ESRI:102003 USA_Contiguous_Albers_Equal_Area_Conic
  aea <- st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
  
  
  #'  Load spatial data
  usa <- st_read("./Shapefiles/tl_2012_us_state/tl_2012_us_state.shp")
  hwys <- st_read("./Shapefiles/GEE/PrimaryRoads_AZ_NM.shp") %>% st_transform(wgs84)
  dem <- terra::rast("./Shapefiles/GEE/DEM_Arizona_NewMexico.tif")
  slope1 <- terra::rast("./Shapefiles/GEE/Slope_Arizona_NewMexico-1.tif")
  slope2 <- terra::rast("./Shapefiles/GEE/Slope_Arizona_NewMexico-2.tif")
  
  #'  Filter and create new sf objects
  #'  AZ & NM polygon
  az_nm <- filter(usa, NAME == "Arizona" | NAME == "New Mexico") %>% st_transform(wgs84) 
  st_crs(az_nm); st_bbox(az_nm)
  #'  Merge into single "southwest" polygon
  southwest <- st_union(az_nm[az_nm$NAME == "Arizona",], az_nm[az_nm$NAME == "New Mexico",]) %>%
    st_cast("POLYGON")
  #'  Create a single LINESTRING for I40
  I40 <- filter(hwys, fullname == "I- 40") %>%
    st_union()
  I10 <- filter(hwys, fullname == "I- 10") %>%
    st_union()
  #'  Join highways that boarder experimental population area
  I40_I10 <- st_union(I40, I10)
  #'  Split southwest polygon by I40
  split_southwest <- lwgeom::st_split(southwest, I40_I10) %>%
    st_collection_extract("POLYGON")
  split_southwest$include <- factor(seq_len(nrow(split_southwest)))
  ggplot(split_southwest) + geom_sf(aes(fill = include))
  #'  Filter to experimental population area and anything south of I40
  exp_pop <- filter(split_southwest, include == 2); st_crs(exp_pop); st_bbox(exp_pop)
  southI40 <- filter(split_southwest, include != 1) %>% st_union(); st_crs(southI40); st_bbox(southI40)
  ggplot(exp_pop) + geom_sf();  ggplot(az_nm) + geom_sf() + geom_sf(data = southI40)
  
  #'  Save select features
  # st_write(az_nm, "./Shapefiles/tl_2012_us_state/Arizona_NewMexico.shp")
  # st_write(az_nm, "./Shapefiles/tl_2012_us_state/Arizona_NewMexico.kml", driver = "kml", delete_dsn = TRUE)
  # st_write(exp_pop, "./Shapefiles/experimental_pop_poly.shp")
  # st_write(southI40, "./Shapefiles/I40_south_poly.shp")
  
  #'  Create a sf object for locations
  spatial_locs <- function(locs, proj) {
    locs <- arrange(locs, Pack_year)
    sf_locs <- st_as_sf(locs, coords = c("Longitude", "Latitude"), crs = wgs84) %>%
      mutate(obs = seq(1:nrow(.))) %>%
      relocate(obs, .before = Pack_year) #%>%
      # st_transform(proj)
    return(sf_locs)
  }
  homesites_wgs84 <- spatial_locs(homesites, wgs84)

  plot(homesites_wgs84[homesites_wgs84$Site_Type == "Den",])  
  plot(homesites_wgs84[homesites_wgs84$Site_Type == "Rendezvous",])  
  
  #'  Save
  # st_write(homesites_wgs84[homesites_wgs84$Site_Type == "Den",], "./Shapefiles/Homesites/homesites_d.kml", driver = "kml", delete_dsn = TRUE)
  # st_write(homesites_wgs84[homesites_wgs84$Site_Type == "Den",], "./Shapefiles/Homesites/homesites_den.shp")
  # st_write(homesites_wgs84[homesites_wgs84$Site_Type == "Rendezvous",], "./Shapefiles/Homesites/homesites_r.kml", driver = "kml", delete_dsn = TRUE)
  # st_write(homesites_wgs84[homesites_wgs84$Site_Type == "Rendezvous",], "./Shapefiles/Homesites/homesites_rendezvous.shp")
  
  #'  Explore homesites by year, pack, and type (den vs rendezvous)
  
  
  
  
  
  
  
  
  