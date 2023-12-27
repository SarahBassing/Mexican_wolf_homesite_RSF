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
  az_nm <- filter(usa, NAME == "Arizona" | NAME == "New Mexico") %>% st_transform(wgs84) 
  st_crs(az_nm); st_bbox(az_nm)
  # st_write(az_nm, "./Shapefiles/tl_2012_us_state/Arizona_NewMexico.shp")
  # st_write(az_nm, "./Shapefiles/tl_2012_us_state/Arizona_NewMexico.kml", driver = "kml", delete_dsn = TRUE)
  dem <- terra::rast("./Shapefiles/GEE/DEM_Arizona_NewMexico.tif")
  dem2 <- terra::rast("./Shapefiles/GEE/DEM_Arizona_NewMexico-2.tif")
  
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
  
  #'  Explore homesites by year, pack, and type (den vs rendezvous)
  
  
  
  
  
  
  
  
  