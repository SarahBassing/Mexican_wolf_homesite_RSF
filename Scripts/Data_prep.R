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
  # slope1 <- terra::rast("./Shapefiles/GEE/Slope_Arizona_NewMexico-1.tif")
  # slope2 <- terra::rast("./Shapefiles/GEE/Slope_Arizona_NewMexico-2.tif")
  
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
  #'  Join highways that border experimental population area
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

  ggplot(exp_pop) + geom_sf() + geom_sf(data = hwys) + 
    geom_sf(data = homesites_wgs84[homesites_wgs84$Site_Type == "Den",], aes(color = Year), shape = 16, size = 3) 
  ggplot(exp_pop) + geom_sf() + geom_sf(data = hwys) + 
    geom_sf(data = homesites_wgs84[homesites_wgs84$Site_Type == "Rendezvous",], aes(color = Year), shape = 16, size = 3) #+
    #' #'  Constrain plot to bbox of the experimental population area
    #' coord_sf(xlim = c(-114.53588, -103.04233), ylim = c(31.96210, 35.53438), expand = FALSE)
  
  #'  Save
  # st_write(homesites_wgs84[homesites_wgs84$Site_Type == "Den",], "./Shapefiles/Homesites/homesites_d.kml", driver = "kml", delete_dsn = TRUE)
  # st_write(homesites_wgs84[homesites_wgs84$Site_Type == "Den",], "./Shapefiles/Homesites/homesites_den.shp")
  # st_write(homesites_wgs84[homesites_wgs84$Site_Type == "Rendezvous",], "./Shapefiles/Homesites/homesites_r.kml", driver = "kml", delete_dsn = TRUE)
  # st_write(homesites_wgs84[homesites_wgs84$Site_Type == "Rendezvous",], "./Shapefiles/Homesites/homesites_rendezvous.shp")
  
  #'  Explore data
  #'  Which den site falls way outside experimental population area?
  home_exp_intersection <- st_intersection(homesites_wgs84, exp_pop)
  ggplot(exp_pop) + geom_sf() + geom_sf(data = hwys) + geom_sf(data = home_exp_intersection, aes(color = Year), shape = 16, size = 3)
  far_away_home <- subset(homesites_wgs84, !(obs %in% home_exp_intersection$obs))
  #'  2023 den of the Manada del Arroyo pack
  #'  This pair was released in the state of Chihuahua, Mexico so excluding from analyses
  homesites_wgs84_usa <- filter(homesites_wgs84, Pack != "Manada del Arroyo")
  
  #'  Split out dens and rendezvous sites
  dens <- filter(homesites_wgs84_usa, Site_Type == "Den") %>%
    #'  Flag packs with multiple den sites in the same year and if any were noted 
    #'  as the natal den
    mutate(NewDen_SameYear = ifelse(Pack_year == lead(Pack_year), "TRUE", "FALSE"),
           NewDen_SameYear = ifelse(Pack_year == lag(Pack_year), "TRUE", NewDen_SameYear),
           NewDen_SameYear = ifelse(is.na(NewDen_SameYear), "FALSE", NewDen_SameYear)#,
           #Comments = ifelse(grepl("natal", Comments), "Natal den", Comments),
           ) %>%
    relocate(NewDen_SameYear, .after = Comments)
  rnds <- filter(homesites_wgs84_usa, Site_Type == "Rendezvous")
  
  #'  Save pack-years with >1 den site
  multiple_dens <- filter(dens, NewDen_SameYear == "TRUE")
  den_moves <- subset(homesites, (Pack_year %in% multiple_dens$Pack_year)) %>%
    filter(Site_Type == "Den")
  write_csv(den_moves, "./Data/Multiple_den_sites.csv")
  
  #'  Create buffer around each homesite (width in meters)
  st_bbox(homesites_wgs84_usa)
  homesite_buffers <- homesites_wgs84_usa %>%
    vect() %>%
    terra::buffer(width = 1000) %>%
    st_as_sf()
  ggplot(exp_pop) + geom_sf() + 
    geom_sf(data = dens, aes(color = Year), shape = 16) +
    geom_sf(data = homesite_buffers[homesite_buffers$Site_Type == "Den",], colour = "black", fill = NA) +
    coord_sf(xlim = c(-109.8417, -107.3039), ylim = c(33.0215, 34.2875), expand = FALSE)
  ggplot(exp_pop) + geom_sf() + 
    geom_sf(data = rnds, aes(color = Year), shape = 16) +
    geom_sf(data = homesite_buffers[homesite_buffers$Site_Type == "Rendezvous",], colour = "black", fill = NA) +
    coord_sf(xlim = c(-109.8417, -107.3039), ylim = c(33.0215, 34.2875), expand = FALSE)
  
  #'  Which homesites are within 250m of each other?  
  close_sites <- function(sites) {
    #'  Create a sequential observation for each unique site
    pack_yr <- dplyr::select(sites, "Pack_year") %>% 
      mutate(obs = seq(1:nrow(.))) %>% st_drop_geometry()
    print(nrow(pack_yr))
    #'  Calculate pairwise distances between all sites
    dist_btwn_sites <- sites$geometry %>% 
      st_distance()
    #'  Identify which pairings were within 100m of each other
    indices <- which(dist_btwn_sites <= units::set_units(250, "m"), arr.ind = TRUE)
    #'  Simplify distance matrix into a dataframe where 1st column = matrix rows,
    #'  2nd column = matrix columns, and 3rd column = distance IF < 100m between
    #'  sites indexed by matrix row and column
    indices_df <- as.data.frame(indices)
    indices_df$dist <- as.numeric(dist_btwn_sites[indices])
    #'  Append pack_year as two new columns corresponding to the sites indexed 
    #'  by the matrix's "row" and "column" 
    pairwise_distances <- full_join(indices_df, pack_yr, by = c("row" = "obs")) %>%
      full_join(pack_yr, by = c("col" = "obs")) %>%
      relocate(Pack_year.x, .before = dist) %>%
      relocate(Pack_year.y, .before = dist) %>%
      #'  Pull out just the pack name per site
      mutate(pack_A = sub("\\_.*", "", Pack_year.x),
             pack_B = sub("\\_.*", "", Pack_year.y),
             pack_pair = paste0(Pack_year.x, "_", Pack_year.y)) %>%
      filter(row != col)
    
    return(pairwise_distances)
  }
  close_dens <- close_sites(dens)
  close_rnd <- close_sites(rnds)
  
  #'  Reduce to packs that used approximately the same site over multiple years 
  #'  (i.e., recurring use of a site across) or relocated to a site <250m away 
  #'  in the same year
  close_dens_skinny <- filter(close_dens, pack_A == pack_B) %>%
    dplyr::select(c("Pack_year.x", "Pack_year.y", "pack_pair", "dist")) %>%
    mutate(Repeat_DenYear = duplicated(Pack_year.x))
  # Need to be able to group dens by pack_years that are spatially overlapping but
  # in some cases a pack will use the same general area repeatedly for several
  # years (e.g, 2012-2014 cluster) then move to a new area and repeatedly use
  # same general area for several more years (e.g., 2015-2016 new cluster). Need
  # to differentiate these two clusters and then randomly draw one location per
  # cluster for that pack.
  
  #'  Filter spatial data to the above sites, then
  #'  https://stackoverflow.com/questions/77179007/group-spatial-points-by-distance-in-r-how-to-group-cluster-spatial-points-so-gr
  adj <- st_distance(points)
  adj <- matrix(as.numeric(as.numeric(adj)) < 11000, nrow = nrow(adj))
  library(igraph)
  g <- graph_from_adjacency_matrix(adj)
  plot(g)
  df$group <- factor(components(g)$membership)
  ggplot(rnaturalearth::ne_countries(scale = 10, returnclass = 'sf')) +
    geom_sf() +
    geom_point(data = df, aes(x = longitude, y = latitude, color = group)) +
    coord_sf(xlim = c(71.5, 72.5), ylim = c(-6, -5))
  