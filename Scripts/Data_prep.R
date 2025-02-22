  #'  ------------------------
  #'  Used and Available locations
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
  library(tidyterra)
  library(spdep)
  library(adehabitatHR)
  library(paletteer)
  library(patchwork)
  library(tidyverse)
  #'  NOTE: dplyr arrange() orders pack names differently depending on package version. 
  #'  Use dplyr 1.1.4 for consistency
  
  #'  ------------------------
  ####  Load and format data  ####
  #'  ------------------------
  #'  Used locations
  homesites <- read_csv("./Data/MexWolf_dens_rend_sites_1998_2023_updated_01.19.24.csv") %>%
    mutate(Pack_year = paste0(Pack, "_", Year)) %>%
    relocate(Pack_year, .before = Year)
    
  #'  Spatial data (rasters)
  elev <- terra::rast("./Shapefiles/Terrain_variables/Mosaic_DEM.tif"); res(elev); crs(elev); st_bbox(elev)
  slope <- terra::rast("./Shapefiles/Terrain_variables/slope.tif"); res(slope); crs(slope); st_bbox(slope)
  rough <- terra::rast("./Shapefiles/Terrain_variables/VRM.tif"); res(rough); crs(rough); st_bbox(rough)
  curve <- terra::rast("./Shapefiles/Terrain_variables/Gaussian_curvature.tif"); res(curve); crs(curve); st_bbox(curve)
  water <- terra::rast("./Shapefiles/National Hydrography Dataset (NHD)/Mosaic_Dist2Water.tif"); res(water); crs(water); st_bbox(water)
  human_mod <- terra::rast("./Shapefiles/Human_variables/mosaic_global_Human_Modification.tif"); res(human_mod); crs(human_mod); st_bbox(human_mod)
  # roads <- terra::rast("./Shapefiles/Human_variables/mosaic_dist2road.tif"); res(roads); crs(roads); st_bbox(roads)
  roads <- terra::rast("./Shapefiles/Human_variables/dist_to_road_updated_NAD83.tif"); res(roads); crs(roads); st_bbox(roads)
  
  #'  Define WGS84 coordinate systems
  wgs84 <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
  
  #'  Recovery zones and suitable habitat (identified by Martinez-Meyer et al. 2021)
  mwepa <- st_read("./Shapefiles/MWEPA Layers Zone 1-3 & Boundary/MWEPA Final.shp"); crs(mwepa)
  mwz1 <- st_read("./Shapefiles/MWEPA Layers Zone 1-3 & Boundary/Final_MWEPA_Zone_1.shp") %>% st_transform(wgs84)
  mwz2 <- st_read("./Shapefiles/MWEPA Layers Zone 1-3 & Boundary/Final_MWEPA_Zone_2.shp") %>% st_transform(wgs84)
  mwz3 <- st_read("./Shapefiles/MWEPA Layers Zone 1-3 & Boundary/Final_MWEPA_Zone_3.shp") %>% st_transform(wgs84)
  suitable_habitat <- st_read("./Shapefiles/Martinez_Meyer_2021_layers/suitable_habitat.shp") %>% st_transform(wgs84)
  mwepa_wgs84 <- st_transform(mwepa, wgs84)
  crs(mwz1); crs(mwepa_wgs84); crs(suitable_habitat)
  
  #'  Define NAD27 and NAD83 projected coordinate system
  nad27_12N <- st_crs(mwepa)
  nad83 <- st_crs(elev)
  
  #'  Crop suitable habitat to experimental population area and reproject
  mwepa_suitable <- st_intersection(mwepa_wgs84, suitable_habitat)
  plot(mwepa_suitable[1])
  mwepa_suitable_nad27 <- st_transform(mwepa_suitable, crs = nad27_12N)
  mwepa_suitable_nad83 <- st_transform(mwepa_suitable, crs = nad83)
  
  #'  State, highway, and waterbody shapefiles
  usa <- st_read("./Shapefiles/tl_2012_us_state/tl_2012_us_state.shp")
  # hwys <- st_read("./Shapefiles/GEE/PrimaryRoads_AZ_NM.shp") %>% st_transform(wgs84)
  az_nm <- filter(usa, NAME == "Arizona" | NAME == "New Mexico") %>% st_transform(wgs84) 
  crs(az_nm); st_bbox(az_nm)
  #'  Merge into single "southwest" polygon
  southwest <- st_union(az_nm[az_nm$NAME == "Arizona",], az_nm[az_nm$NAME == "New Mexico",]) %>%
    st_cast("POLYGON")
  waterbodies <- st_read("./Shapefiles/National Hydrography Dataset (NHD)/AZ_NM_waterbodies_NHD.shp"); crs(waterbodies)
  #'  Identify large bodies of water (anything larger than 1 sq-km in size)
  bigwater <- waterbodies[waterbodies$areasqkm > 1,]
  bigwater_nad27 <- st_transform(bigwater, crs = nad27_12N)
  bigwater_nad83 <- st_transform(bigwater, crs = nad83)
  
  #'  Review each zones & suitable habitat within context of larger experimental 
  #'  population area boundary
  ggplot(mwepa_wgs84) + geom_sf() + geom_sf(data = mwz1) #'  Recovery zone 1
  ggplot(mwepa_wgs84) + geom_sf() + geom_sf(data = mwz2) #'  Recovery zone 2
  ggplot(mwepa_wgs84) + geom_sf() + geom_sf(data = mwz3) #'  Recovery zone 3
  ggplot(mwepa_wgs84) + geom_sf() + geom_sf(data = mwepa_suitable, fill = "blue") + geom_sf(data = mwz1, color = "orange", fill = NA)
  
  #'  Recovery zone 1 bounding box
  st_bbox(mwz1)
  mwz1_bbox <- st_as_sfc(st_bbox(mwz1))
  mwz1_extent_wgs84 <- st_transform(mwz1_bbox, wgs84)
  
  #'  Create a sf object for locations
  #'  ####  NOTE: dplyr arrange() orders pack names differently depending on package version. Use dplyr 1.1.4 to be consistent
  spatial_locs <- function(locs, proj) {
    locs <- arrange(locs, Pack_year)      
    sf_locs <- st_as_sf(locs, coords = c("Longitude", "Latitude"), crs = wgs84) %>%
      mutate(obs = seq(1:nrow(.))) %>%
      relocate(obs, .before = Pack_year) %>%
      st_transform(proj)
    return(sf_locs)
  }
  homesites_wgs84 <- spatial_locs(homesites, wgs84)
  homesites_nad27 <- spatial_locs(homesites, nad27_12N)
  homesites_nad83 <- spatial_locs(homesites, nad83)

  #'  Save
  # st_write(homesites_wgs84[homesites_wgs84$Site_Type == "Den",], "./Shapefiles/Homesites/homesites_d.kml", driver = "kml", delete_dsn = TRUE)
  # st_write(homesites_wgs84[homesites_wgs84$Site_Type == "Den",], "./Shapefiles/Homesites/homesites_den.shp")
  # st_write(homesites_wgs84[homesites_wgs84$Site_Type == "Rendezvous",], "./Shapefiles/Homesites/homesites_r.kml", driver = "kml", delete_dsn = TRUE)
  # st_write(homesites_wgs84[homesites_wgs84$Site_Type == "Rendezvous",], "./Shapefiles/Homesites/homesites_rendezvous.shp")
  
  
  #'  ------------------------------
  ####  Explore & filter homesites  ####
  #'  ------------------------------
  #'  Visualize homesites within Recovery Zone 1
  ggplot(mwz1_bbox) + geom_sf() + geom_sf(data = mwz1) +
    geom_sf(data = homesites_wgs84[homesites_wgs84$Site_Type == "Den",], aes(color = Year), shape = 16, size = 3) 
  ggplot(mwz1_bbox) + geom_sf() + geom_sf(data = mwz1) + 
    geom_sf(data = homesites_wgs84[homesites_wgs84$Site_Type == "Rendezvous",], aes(color = Year), shape = 16, size = 3) #+
    #' #'  Constrain plot to bbox of the experimental population area
    #' coord_sf(xlim = c(-114.53588, -103.04233), ylim = c(31.96210, 35.53438), expand = FALSE)
    
  #'  Which den site falls way outside experimental population area?
  home_exp_intersection <- st_intersection(homesites_wgs84, mwepa_wgs84)
  # ggplot(mwepa_wgs84) + geom_sf() + geom_sf(data = hwys) + geom_sf(data = home_exp_intersection, aes(color = Year), shape = 16, size = 3)
  # far_away_home <- subset(homesites_wgs84, !(obs %in% home_exp_intersection$obs))
  #'  2023 den of the Manada del Arroyo pack
  #'  This pair was released in the state of Chihuahua, Mexico so excluding from analyses
  homesites_wgs84_usa <- filter(homesites_wgs84, Pack != "Manada del Arroyo")
  homesites_nad27_usa <- st_transform(homesites_wgs84_usa, nad27_12N)
  homesites_nad83_usa <- st_transform(homesites_wgs84_usa, nad83)
  
  #'  What's the min and max elevation (meters) of homesites? Using 30m res DEM
  homesite_elev <- terra::extract(elev, vect(homesites_wgs84_usa)) 
  names(homesite_elev) <- c("ID", "elevation")
  range(homesite_elev$elevation) # min = 1718, max = 3218
  
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
  
  #'  Unique packs
  length(unique(dens$Pack))
  length(unique(rnds$Pack))
  
  #'  Unique pack-years
  length(unique(dens$Pack_year))
  length(unique(rnds$Pack_year))
  
  #'  Unique homesites
  nrow(dens)
  nrow(rnds)
  
  #'  ---------------------------------
  ##### Group and rarify repeat sites  ####
  #'  ---------------------------------
  #'  Pack-years with >1 den site
  multiple_dens <- filter(dens, NewDen_SameYear == "TRUE")
  den_moves <- subset(homesites, (Pack_year %in% multiple_dens$Pack_year)) %>%
    filter(Site_Type == "Den")
  # write_csv(den_moves, "./Data/Multiple_den_sites.csv")
  
  #'  Create buffer around each homesite (width in meters)
  st_bbox(homesites_wgs84_usa)
  st_bbox(homesites_nad27_usa)
  homesite_buffers <- homesites_wgs84_usa %>% #homesites_nad27_usa
    vect() %>%
    terra::buffer(width = 1000) %>%
    st_as_sf()
  ggplot(mwepa_wgs84) + geom_sf() + 
    geom_sf(data = dens, aes(color = Year), shape = 16) +
    geom_sf(data = homesite_buffers[homesite_buffers$Site_Type == "Den",], colour = "black", fill = NA) +
    coord_sf(xlim = c(-109.8417, -107.3039), ylim = c(33.0215, 34.2875), expand = FALSE) # mwepa_sbb_defined
    #coord_sf(xlim = c(606823.5, 840562.3), ylim = c(3656509.8, 3798894.4), expand = FALSE) # mwepa
  ggplot(mwepa_wgs84) + geom_sf() + 
    geom_sf(data = rnds, aes(color = Year), shape = 16) +
    geom_sf(data = homesite_buffers[homesite_buffers$Site_Type == "Rendezvous",], colour = "black", fill = NA) +
    coord_sf(xlim = c(-109.8417, -107.3039), ylim = c(33.0215, 34.2875), expand = FALSE) # mwepa_sbb_defined
    #coord_sf(xlim = c(606823.5, 840562.3), ylim = c(3656509.8, 3798894.4), expand = FALSE) # mwepa
    
  #'  Identify pairwise combinations of sites that are within 250m of each other
  close_sites <- function(sites) {
    #'  Create a sequential observation for each unique site
    pack_yr <- dplyr::select(sites, "Pack_year") %>% 
      mutate(obs = seq(1:nrow(.))) %>% st_drop_geometry()
    print(nrow(pack_yr))
    #'  Calculate pairwise distances between all sites
    dist_btwn_sites <- sites$geometry %>% 
      st_distance()
    #'  Identify which pairings were within 250m of each other
    indices <- which(dist_btwn_sites <= units::set_units(250, "m"), arr.ind = TRUE)
    #'  Simplify distance matrix into a data frame where 1st column = matrix rows,
    #'  2nd column = matrix columns, and 3rd column = distance IF < 250m between
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
  #'  These dfs an be used as comparison to the cluster results below
  
  #'  Group sites based on their proximity and assign a unique cluster ID per group
  #'  Proximity based on 250m (0.25km)
  homesite_clusters <- function(sites) {
    site_groups <- sites %>%
      #'  Identify neighbors of region points by Euclidean distance
      #'  Distances measured in km if coordinates are geographic lat/long
      mutate(cluster_id = dnearneigh(geometry, 0, 0.25) %>% 
               #'  Find the number of distjoint connected subgraphs (point clusters)
               n.comp.nb() %>% 
               #' Extract the comp.id generated with n.comp.nb and append to df
               getElement("comp.id")) %>%
      #'  Count the number of sites in each group
      add_count(cluster_id) %>%
      relocate(cluster_id, .after = "Pack") %>%
      relocate(n, .after = "cluster_id")
    return(site_groups)
  }
  den_clusters <- homesite_clusters(dens)
  rnd_clusters <- homesite_clusters(rnds)
  
  #'  Assign a confidence weight to sites based on comments in data or pers. comm.
  #'  from J. Oakleaf & A. Greenleaf
  #'  Prior to 2012: located dens via flights (higher error)
  #'  2012 - 2015: identified dens via den visits (lower error)
  #'  2015-present: identified dens via GPS clusters & den visits if cross-fostering occurred
  #'  Most rendezvous sites identified via GPS clusters and not visited
  den_weights <- den_clusters %>%
    #'  Highest weight goes to sites that were confirmed or have detailed information about the site (visited)
    mutate(wgts = ifelse(grepl("Confirmed", Comments), 4, NA),
           wgts = ifelse(grepl("natal", Comments), 4, wgts),
           wgts = ifelse(grepl("Natal", Comments), 4, wgts),
           wgts = ifelse(grepl("Den Description:", Comments), 4, wgts),
           wgts = ifelse(grepl("Failed", Comments), 4, wgts),
           wgts = ifelse(grepl("failed", Comments), 4, wgts), 
           #'  Downgrade sites with den descriptions that were not the natal den 
           #'  (moved from original den so selection process was slightly different)
           wgts = ifelse(Pack_year == "Iron Creek_2015" & !grepl("Natal den", Comments), 3, wgts),
           wgts = ifelse(Pack_year == "Maverick_2014" & !grepl("natal", Comments), 3, wgts),
           # wgts = ifelse(Pack_year == "Rim_2014" & !grepl("Natal den", Comments), 3, wgts),
           #'  Middle weight goes to sites with strong evidence, often identified via GPS clustering
           wgts = ifelse(grepl("Strong", Comments), 3, wgts),
           wgts = ifelse(grepl("based on GPS", Comments), 3, wgts),
           wgts = ifelse(grepl("cluster", Comments), 3, wgts),
           wgts = ifelse(Pack_year == "Willow Creek_2022", 3, wgts), # assuming ID'd by GPS locations
           #'  Lowest weight goes to sites with weak evidence or based on approximate triangulation
           wgts = ifelse(grepl("Weak", Comments), 1, wgts),
           wgts = ifelse(Pack_year == "Luna_2019" | Pack_year == "Luna_2020", 1, wgts),
           #'  Medium-low weight to assumed natal den site but had weak evidence for location
           wgts = ifelse(Pack_year == "Hawks Nest_2001" & grepl("natal den", Comments), 2, wgts),
           #'  Assign low to medium weights for sites with no comments based on 
           #'  year and type of monitoring at that time (<2012 flights, 2012-2015 
           #'  den visits, 2015+ GPS clusters)
           wgts = ifelse(is.na(wgts) & Year < 2012, 1, wgts),
           wgts = ifelse(is.na(wgts) & Year >= 2015, 2, wgts),
           wgts = ifelse(is.na(wgts), 3, wgts))  
  
  rnd_weights <- rnd_clusters %>%
    #'  Highest weight goes to sites that were confirmed or have detailed information about pups/the site
    mutate(wgts = ifelse(grepl("Confirmed", Comments), 4, NA),
           wgts = ifelse(grepl("visual", Comments), 4, wgts),
           wgts = ifelse(grepl("pups", Comments), 4, wgts),
           wgts = ifelse(grepl("Rendezvous", Comments), 4, wgts),
           #'  Middle weight goes to sites with strong evidence, often identified via GPS clustering
           wgts = ifelse(grepl("Strong", Comments), 3, wgts),
           wgts = ifelse(grepl("Moved", Comments), 3, wgts), 
           #'  Lowest weight goes to sites with weak evidence
           wgts = ifelse(grepl("Weak", Comments), 1, wgts),
           #'  Assign low to medium weights for sites with no comments based on 
           #'  year and type of monitoring at that time (<2012 flights, 2012-2015 
           #'  den visits, 2015+ GPS clusters)
           wgts = ifelse(is.na(wgts) & Year < 2012, 1, wgts),
           wgts = ifelse(is.na(wgts) & Year >= 2015, 2, wgts),
           wgts = ifelse(is.na(wgts), 3, wgts))
  
  #'  Select one site per pack per cluster to keep for RSF analyses
  #'  Need to exclude repeat use of sites by the same pack to keep each "used" 
  #'  observation in the RSF as independent as possible (i.e., some packs re-use
  #'  same den or rendezvous every year - multiple used points at this location 
  #'  is not representative of what the population is selecting for)
  reduce_sites <- function(sites) {
    #' #'  Set seed so sampling is reproducible
    #' set.seed(2024)
    site_random_sample <- sites %>%
      #'  Grouping by pack and cluster allows same site to be used by different
      #'  packs but prevents same site from being used by one pack multiple times
      group_by(Pack, cluster_id) %>%
      #'  Arrange sites by confidence level in descending order
      arrange(-wgts) %>%
      #'  Retain the site with the highest level of confidence in each cluster
      slice(1L) %>%
      #' #'  Weight observations by strength of confidence and draw 1 sample
      #' slice_sample(n = 1, weight_by = wgts) %>%
      ungroup() %>%
      arrange(Pack, cluster_id)
    return(site_random_sample)
  }
  den_sample <- reduce_sites(den_weights) %>%
    dplyr::select(-NewDen_SameYear)
  rnd_sample <- reduce_sites(rnd_weights)

  #'  Merge sampled den and rendezvous sites into single df
  used_homesites <- den_sample %>% bind_rows(rnd_sample)
  used_homesites_nad27 <- st_transform(used_homesites, nad27_12N)
  used_homesites_nad83 <- st_transform(used_homesites, nad83)
  
  #'  Final homesite count
  nrow(used_homesites[used_homesites$Site_Type == "Den",])
  nrow(used_homesites[used_homesites$Site_Type == "Rendezvous",])

  #'  --------------------------------
  ####  Generate available locations  ####
  #'  --------------------------------
  #####  Create buffered sampling extent  #####
  #'  ------------------------------------
  #'  Define extent of "available" habitat for wolves to den/rendezvous in
  #'  Create single MCP that includes den and rendezvous sites
  #'  NOTE the coordinate system - buffer is smoother in projected coord system
  used_homesites_sp <- as(used_homesites_nad27, "Spatial"); proj4string(used_homesites_sp)
  homesite_mcp <- mcp(used_homesites_sp, percent = 100) # ignore warnings 
  homesite_mcp_sf <- st_as_sf(homesite_mcp)
  
  #'  Hacky way to estimate "radius" of polygon (needed as a buffering distance)
  (mcp_radius_sqr <- as.numeric(sqrt(st_area(homesite_mcp_sf))/2)) #if polygon was a perfect square
  (mcp_radius_crl <- as.numeric(sqrt(st_area(homesite_mcp_sf)/pi)/2)) #if polygon was a perfect circle
  
  #'  Buffer MCP by its ~radius so the available habitat extends beyond known used sites
  homesite_mcp_buff <- st_buffer(homesite_mcp_sf, mcp_radius_crl)
  #'  How much bigger is the buffered MCP compared to the original?
  st_area(homesite_mcp_buff)/st_area(homesite_mcp_sf)
  
  #'  Mask out large water bodies so not available when drawing random locations
  homesite_mcp_buff_watermask <- st_difference(homesite_mcp_buff, st_union(bigwater_nad27))
  plot(homesite_mcp_buff_watermask[1])
  
  #'  Mask out unsuitable habitat so not available
  #'  First generate polygon of unsuitable habitat
  homesite_mcp_buff_UNsuitablemask <- st_difference(homesite_mcp_buff, st_union(mwepa_suitable_nad27))
  #'  Mask out unsuitable habitat from buffer with large water bodies already masked out
  homesite_mcp_buff_suitablemask <- st_difference(homesite_mcp_buff_watermask, st_union(homesite_mcp_buff_UNsuitablemask))
  plot(homesite_mcp_buff_suitablemask[1])
  
  buff_bbox <- st_bbox(homesite_mcp_buff)
  
  #'  Visualize (note the coordinate system!)
  #'  100% MCP and buffered MCP
  ggplot(homesite_mcp_buff) + geom_sf() + geom_sf(data = homesite_mcp_sf)
  #'  Den/rendezvous sites within maksed & buffered MCP and Zone 1 for context within Exp. Pop. Area
  ggplot(st_transform(mwepa, nad27_12N)) + geom_sf(data = mwepa_suitable_nad27, fill = "gray65") + 
    geom_sf(data = homesite_mcp_buff_suitablemask, color = "red", fill = NA) + 
    # geom_sf(data = homesite_mcp_sf, color = "blue", fill = NA) + 
    geom_sf(data = st_transform(mwz1, nad27_12N), color = "gray15", fill = NA, size = 0.7) + #fill = "gray15", alpha = 0.50) +
    geom_sf(data = homesites_nad27_usa[homesites_nad27_usa$Site_Type == "Den",], aes(color = Year), shape = 16, size = 1.5) +
    coord_sf(xlim = c(buff_bbox[1], buff_bbox[3]), ylim = c(buff_bbox[2], buff_bbox[4]), expand = TRUE) +
    ggtitle("Den sites, MWZ1, and buffered MCP (excluding unsuitable habitat)")
  ggplot(st_transform(mwepa, nad27_12N)) + geom_sf(data = mwepa_suitable_nad27, fill = "gray65") + 
    geom_sf(data = homesite_mcp_buff_suitablemask, color = "red", fill = NA) + 
    # geom_sf(data = homesite_mcp_sf, color = "blue", fill = NA) + 
    geom_sf(data = st_transform(mwz1, nad27_12N), color = "gray15", fill = NA, size = 0.7) + #, fill = "gray15", alpha = 0.50) +
    geom_sf(data = homesites_nad27_usa[homesites_nad27_usa$Site_Type == "Rendezvous",], aes(color = Year), shape = 16, size = 1.5) +
    coord_sf(xlim = c(buff_bbox[1], buff_bbox[3]), ylim = c(buff_bbox[2], buff_bbox[4]), expand = TRUE) +
    ggtitle("Rendezvous sites, MWZ1, and buffered MCP (excluding unsuitable habitat)")
  
  #'  Save as shapefiles
  homesite_mcp_buff_wgs84 <- st_transform(homesite_mcp_buff, wgs84); st_bbox(homesite_mcp_buff_wgs84)
  homesite_mcp_buff_suitablemask_wgs84 <- st_transform(homesite_mcp_buff_suitablemask, wgs84); st_bbox(homesite_mcp_buff_suitablemask_wgs84)
  # st_write(homesite_mcp_buff_wgs84, "./Shapefiles/Homesites/Homesite_buffered_MCP.shp")
  # st_write(homesite_mcp_buff_suitablemask, "./Shapefiles/Homesites/Homesite_buffered_MCP_suitableHabitat.shp")
  
  #'  Area of buffer (without unsuitable habitat masked out)
  buff_area <- st_area(homesite_mcp_buff_wgs84)
  require(units)
  set_units(buff_area, km^2)
  
  #'  Area of final buffer with unsuitable habitat excluded
  suitable_buff_area <- st_area(homesite_mcp_buff_suitablemask)
  set_units(suitable_buff_area, km^2)
  
  #'  ------------------------------------------------
  #####  Generate random points within buffered area  #####
  #'  ------------------------------------------------
  #'  Number of available locations to generate per used location 
  #'  Draw 1000 to use for sensitivity analysis
  # avail_pts <- 1000
  #'  Draw 100 to use for final analyses based on results of sensitivity analysis
  avail_pts <- 100
  
  #'  Function to generate random available locations based on number of used locations
  sample_avail_locs <- function(locs, newcrs, navail, seed, sitetype) {
    #'  Reproject locations
    locs <- st_transform(locs, newcrs)
    
    #'  Number of used locations
    nused <- nrow(locs)
    print(nused)
    
    #'  Total number of available locations to generate
    navailable <- nused * navail
    
    #'  Sequence of pack_year repeated avail_pts times
    #'  (need to know which available points correspond to each site)
    navailable_packyear <- as.vector(rep(locs$Pack_year, each = avail_pts))
     
    #'  Set seed for reproducibility
    set.seed(seed)
    
    #'  Draw random sample of locations within the buffered MCP (excluding large waterbodies)
    rndpts <- st_sample(homesite_mcp_buff_suitablemask, size = navailable, type = "random", exact = TRUE) %>%
      #'  Reformat to a normal sf object
      st_as_sf() %>%
      mutate(obs = seq(1:nrow(.)),
             Site_Type = sitetype,
             used = 0,
             #'  Add weights to used/available locations (used = 1, available = 5000 
             #'  per Fieberg et al. 2021)
             wgts = 5000) %>%
      bind_cols(navailable_packyear) %>%
      relocate(x, .after = last_col()) %>%
      relocate(used, .after = x) %>%
      relocate(Site_Type, .before = used) %>%
      relocate(wgts, .after = used) %>%
      rename(geometry = x)
    names(rndpts) <- c("obs", "Pack_year", "Site_Type", "used", "wgts", "geometry")
    
    #'  Plot available points within buffered MCP
    avail_locs_plot <- ggplot(mwepa) + geom_sf() + 
      geom_sf(data = homesite_mcp_buff, color = "red", size = 1.2) + 
      geom_sf(data = mwz1, fill = "gray25", alpha = 0.30) +
      geom_sf(data = rndpts, shape = 16, size = 0.5, color = "blue") 
    plot(avail_locs_plot)
  
    return(rndpts)
  }
  avail_locs_den <- sample_avail_locs(den_sample, newcrs = nad27_12N, navail = avail_pts, seed = 108, sitetype = "Den")
  avail_locs_rnd <- sample_avail_locs(rnd_sample, newcrs = nad27_12N, navail = avail_pts, seed = 211, sitetype = "Rendezvous")
  
  #'  Reproject available locations
  avail_locs_den_nad83 <- st_transform(avail_locs_den, nad83)
  avail_locs_den_wgs84 <- st_transform(avail_locs_den, wgs84)
  avail_locs_rnd_nad83 <- st_transform(avail_locs_rnd, nad83)
  avail_locs_rnd_wgs84 <- st_transform(avail_locs_rnd, wgs84)
  
  #' #'  Save available locations 
  #' st_write(avail_locs_den_wgs84, "./Shapefiles/Homesites/Available_locations_den.shp")
  #' st_write(avail_locs_rnd_wgs84, "./Shapefiles/Homesites/Available_locations_rnd.shp")
  
  #'  Combine used and available locations per homesite type 
  all_locs <- function(used, avail) {
    #'  Add "used" classification to used locations
    used_and_avail <- used %>%
      mutate(used = 1,
             #'  Add weights to used/available locations (used = 1, available = 5000 
             #'  per Fieberg et al. 2021)
             wgts = 1) %>%
      dplyr::select(-c(cluster_id, n, Year, Pack, Comments)) %>%
      #'  Bind available locations to the used sites
      bind_rows(avail) %>%
      #'  Reogranize so all used & available locations per pack_year are sequential
      arrange(Pack_year) %>%
      #'  Renumber all used and available locations sequentially
      mutate("ID" = seq(1:nrow(.))) %>%
      relocate(used, .after = "Site_Type") %>%
      relocate(ID, .before = obs) %>%
      dplyr::select(-obs)
    return(used_and_avail)
  }
  den_locs_wgs84 <- all_locs(den_sample, avail = avail_locs_den_wgs84) 
  rnd_locs_wgs84 <- all_locs(rnd_sample, avail = avail_locs_rnd_wgs84)
  
  #'  Transform projection for covariate extraction
  den_locs_nad83 <- st_transform(den_locs_wgs84, nad83)
  rnd_locs_nad83 <- st_transform(rnd_locs_wgs84, nad83)
  
  #' #'  Save all locations
  #' st_write(den_locs_wgs84, "./Shapefiles/Homesites/Used_Available_locations_den_updated_12.16.24.shp")
  #' st_write(rnd_locs_wgs84, "./Shapefiles/Homesites/Used_Available_locations_rnd_updated_12.16.24.shp")
  
  #'  ---------------------
  ####  Gather covariates  ####
  #'  ---------------------
  #'  Read in covariates extracted from Google Earth Engine - this requires uploading
  #'  used/avail locations to GEE and extracting data
  #'  Apply 2000 data to sites older than 2000 and 2022 data to sites from 2023
  #'  due to temporal extent of MODIS and Hansen data
  canopy_den <- read_csv("./Data/GEE extracted data/percent_canopy_denSeason_updated121624.csv") %>%
    dplyr::select(-site_year)
  canopy_rnd <- read_csv("./Data/GEE extracted data/percent_canopy_rndSeason_updated121624.csv") %>%
    dplyr::select(-site_year)
  ndvi_den <- read_csv("./Data/GEE extracted data/mean_NDVI_denSeason_updated121624.csv")
  ndvi_rnd <- read_csv("./Data/GEE extracted data/mean_NDVI_rndSeason_updated121624.csv")
  
  #'  Create raster stack of all terrain covariates (must be same grid & res)
  terrain_stack <- c(elev, slope, rough, curve) 
  
  #'  Extract covariate values at each used and available location
  get_covs <- function(locs_nad83, locs_wgs84, canopy, ndvi) {
    #'  Extract covariate values at each location (use correct projections!)
    terrain <- terra::extract(terrain_stack, locs_nad83) 
    names(terrain) <- c("ID", "Elevation_m", "Slope_degrees", "Roughness_VRM", "Gaussian_curvature")
    rds <- terra::extract(roads, locs_nad83) %>% rename("Nearest_road_m" = "distance") #"mosaic_dist2road
    gHM <- terra::extract(human_mod, locs_nad83) %>% rename("Human_mod_index" = "mosaic_global_Human_Modification")
    h20 <- terra::extract(water, locs_wgs84) %>% rename("Nearest_water_m" = "Mosaic_Dist2Water")
    
    #'  Grab homesite info
    homesite <- locs_wgs84 %>%
      dplyr::select(c("Pack_year", "Site_Type", "used", "wgts")) %>%
      rename("Homesite" = "Site_Type") %>%
      mutate(ID = seq(1:nrow(.))) %>%
      relocate(ID, .before = "Pack_year")
    
    #'  Combine homesite info and covariates into single data frame
    df <- full_join(homesite, terrain, by = "ID") %>%
      full_join(h20, by = "ID") %>%
      full_join(gHM, by = "ID") %>%
      full_join(rds, by = "ID") %>%
      full_join(canopy, by = c("ID", "Pack_year", "used")) %>%
      full_join(ndvi, by = c("ID", "Pack_year", "used")) %>%
      dplyr::select(-c(Canopy_year, NDVI_year))
    
    #'  Review data
    print(summary(df))
    
    return(df)
  }
  all_data_den <- get_covs(den_locs_nad83, locs_wgs84 = den_locs_wgs84, canopy = canopy_den, ndvi = ndvi_den)
  all_data_rnd <- get_covs(rnd_locs_nad83, locs_wgs84 = rnd_locs_wgs84, canopy = canopy_rnd, ndvi = ndvi_rnd)
  
  #'  Save covariate data for all used and available locations
  # write_csv(all_data_den, "./Data/all_data_den_updated121624.csv")
  # write_csv(all_data_rnd, "./Data/all_data_rnd_updated121624.csv")
  # write_csv(all_data_den, "./Data/all_data_den_1to1000ratio.csv")
  # write_csv(all_data_rnd, "./Data/all_data_rnd_1to1000ratio.csv")
  
  #'  -----------------------------------------
  ####  Distance btwn den and rendezvous sites ####
  #'  -----------------------------------------
  #'  Identify pairwise combinations of den and rendezvous sites 
  close_sites <- function(sites) {
    #'  Create a sequential observation for each unique site
    pack_yr <- sites %>%
      dplyr::select(c("Pack_year", "Pack", "Year", "Site_Type")) %>% 
      mutate(Pack_yr_site = paste0(Pack_year, "_", Site_Type),
             obs = seq(1:nrow(.))) %>% st_drop_geometry() 
    print(nrow(pack_yr))
    #'  Calculate pairwise distances between all sites
    dist_btwn_sites <- sites$geometry %>% 
      st_distance()
    #'  Identify which pairings were within 1000000m of each other
    indices <- which(dist_btwn_sites <= units::set_units(1000000, "m"), arr.ind = TRUE)
    #'  Simplify distance matrix into a data frame where 1st column = matrix rows,
    #'  2nd column = matrix columns, and 3rd column = distance IF < 1000000m between
    #'  sites indexed by matrix row and column
    indices_df <- as.data.frame(indices)
    indices_df$dist <- as.numeric(dist_btwn_sites[indices])
    #'  Append pack_year as two new columns corresponding to the sites indexed 
    #'  by the matrix's "row" and "column" 
    pairwise_distances <- full_join(indices_df, pack_yr, by = c("row" = "obs")) %>%
      full_join(pack_yr, by = c("col" = "obs")) %>%
      relocate(Pack_yr_site.x, .before = dist) %>%
      relocate(Pack_yr_site.y, .before = dist) %>%
      #'  Pull out just the pack name per site
      mutate(pack_A = sub("\\_.*", "", Pack_yr_site.x),
             pack_B = sub("\\_.*", "", Pack_yr_site.y),
             pack_pair = paste0(Pack_yr_site.x, "_", Pack_yr_site.y)) %>%
      filter(row != col) %>%
      #'  Reduce to just pairs of den and rendezvous sites from the same pack in the same year
      filter(Pack.x == Pack.y) %>%
      filter(Site_Type.x != Site_Type.y) %>%
      filter(Year.x == Year.y) %>%
      #'  Focus on data from 2015 or later when GPS collars were used and rendezvous sites were 
      #'  sampled more frequently
      filter(Year.x >= 2015) %>%
      #'  Reduce to just one observation per pair
      distinct(dist, .keep_all = TRUE)
    
    #'  Report summary stats and visualize distribution of data
    print(summary(pairwise_distances$dist))
    boxplot(pairwise_distances$dist)
    hist(pairwise_distances$dist)
    
    return(pairwise_distances)
  }
  close_sites_same_yr <- close_sites(homesites_wgs84_usa)
  quantile(close_sites_same_yr$dist, 0.8)
  (sites_within_2.5km <- length(close_sites_same_yr$dist[close_sites_same_yr$dist < 2500]))
  (sites_within_5km <- length(close_sites_same_yr$dist[close_sites_same_yr$dist < 5000]))
  (sites_within_100km <- length(close_sites_same_yr$dist))
  sites_within_2.5km/sites_within_100km
  sites_within_5km/sites_within_100km
  
  #'  -----------------------
  ####  Explore covaraites  ####
  #'  -----------------------
  #'  Maximum used elevations
  max(all_data_den$Elevation_m[all_data_den$used == 1])
  max(all_data_rnd$Elevation_m[all_data_rnd$used == 1])
  
  #'  Maximum available elevations
  max(all_data_den$Elevation_m[all_data_den$used == 0])
  max(all_data_rnd$Elevation_m[all_data_rnd$used == 0])
  
  #'  Summarize available covariate data
  summary(all_data_den)
  summary(all_data_rnd) 
  
  #'  Compare spread of covaraite values between use and available locations
  all_data_den <- mutate(all_data_den, used = ifelse(used == 0, "Available", "Used")) %>%
    mutate(used = factor(used))
  all_data_rnd <- mutate(all_data_rnd, used = ifelse(used == 0, "Available", "Used")) %>%
    mutate(used = factor(used))
  
  #'  ------------------------
  #####  Visualize covariates #####
  #'  ------------------------
  ######  Calculate mean of each grouped covaraite  ######
  cov_means <- function(dat) {
    means <- dat %>% 
      group_by(used) %>%
      summarise(across(where(is.numeric), ~round(mean(.x), 2))) %>%
      ungroup() %>%
      dplyr::select(-c(ID, wgts, geometry))
    return(means)
  }
  
  ######  Density plots: usd vs available  ######
  density_plots <- function(dat, cov, col_id, x1, x2, y1, y2, legend_title, xtitle) {
    #'  Call cov_means function to cacluate covariate mean by group
    means <- cov_means(dat)
    print(means)
    
    #'  Create a data frame to hold relevant values and xy coordinates
    annot <- data.frame(loc_type = c("Available", "Used"), 
                        means = c(means[,col_id]),
                        x = c(x1, x2),
                        y = c(y1, y2)) %>%
      dplyr::select(-means.geometry) 
    names(annot) <- c("loc_type", "used", "x", "y")
    annot <- mutate(annot, used = as.character(used))
    
    d_plot <- ggplot(dat, aes(x = cov, group = used, color = used, fill = used)) +
      geom_density(alpha = 0.4, linewidth = 0.25) + theme_classic() +
      scale_color_manual(values = c("Available" = "#6699CC", "Used" = "#993300"), 
                         name = paste0("Sampled locations \n", legend_title, " RSF")) +
      scale_fill_manual(values = c("Available" = "#6699CC", "Used" = "#993300"), 
                         name = paste0("Sampled locations \n", legend_title, " RSF")) +
      geom_text(data = annot, aes(x = x, y = y, label = paste0(loc_type, " (mean = ", used, ")"), color = loc_type), 
                hjust = 0, size = 2, show.legend  = FALSE) + #
      xlab(xtitle) + ylab("Density") + 
      theme(legend.position = "none",
            axis.line = element_line(color = 'black', linewidth = 0.2),
            axis.ticks = element_line(colour = "black", linewidth = 0.2)) 
    print(d_plot)
    return(d_plot)
  }
  elev_den <- density_plots(all_data_den, cov = all_data_den$Elevation_m, col_id = 2, x1 = 1000, x2 = 1000, y1 = 0.00155, y2 = 0.00185, 
                            legend_title = "Den", xtitle = "Elevation (m)")
  elev_rnd <- density_plots(all_data_rnd, cov = all_data_rnd$Elevation_m, col_id = 2, x1 = 1000, x2 = 1500, y1 = 0.0016, y2 = 0.0025, 
                            legend_title = "Rendezvous site", xtitle = "Elevation (m)")
  slope_den <- density_plots(all_data_den, cov = all_data_den$Slope_degrees, col_id = 3, x1 = 5, x2 = 27, y1 = 0.08, y2 = 0.037, 
                             legend_title = "Den", xtitle = "Slope (degrees)")
  slope_rnd <- density_plots(all_data_rnd, cov = all_data_rnd$Slope_degrees, col_id = 3, x1 = 5, x2 = 12, y1 = 0.08, y2 = 0.055, 
                             legend_title = "Rendezvous site", xtitle = "Slope (degrees)")
  rough_den <- density_plots(all_data_den, cov = all_data_den$Roughness_VRM, col_id = 4, x1 = 0, x2 = 0, y1 = 2.2, y2 = 3.05, 
                             legend_title = "Den", xtitle = "Terrain roughness (VRM)")
  rough_rnd <- density_plots(all_data_rnd, cov = all_data_rnd$Roughness_VRM, col_id = 4, x1 = 0, x2 = 0, y1 = 2.22, y2 = 2, 
                             legend_title = "Rendezvous site", xtitle = "Terrain roughness (VRM)")
  curve_den <- density_plots(all_data_den, cov = all_data_den$Gaussian_curvature, col_id = 5, x1 = -0.0002, x2 = 0.00005, y1 = 150000, y2 = 20000, 
                             legend_title = "Den", xtitle = "Surface curvature")
  curve_rnd <- density_plots(all_data_rnd, cov = all_data_rnd$Gaussian_curvature, col_id = 5, x1 = -0.00045, x2 = -0.00045, y1 = 150000, y2 = 100000, 
                             legend_title = "Rendezvous site", xtitle = "Surface curvature")
  water_den <- density_plots(all_data_den, cov = all_data_den$Nearest_water_m, col_id = 6, x1 = 1500, x2 = 1500, y1 = 0.00045, y2 = 0.00052, 
                             legend_title = "Den", xtitle = "Distance to nearest water (m)")
  water_rnd <- density_plots(all_data_rnd, cov = all_data_rnd$Nearest_water_m, col_id = 6, x1 = 1500, x2 = 1500, y1 = 0.00047, y2 = 0.00054, 
                             legend_title = "Rendezvous site", xtitle = "Distance to nearest water (m)")
  human_den <- density_plots(all_data_den, cov = all_data_den$Human_mod_index, col_id = 7, x1 = 0.1, x2 = 0.1, y1 = 33, y2 = 40, 
                             legend_title = "Den", xtitle = "Human modification index")
  human_rnd <- density_plots(all_data_rnd, cov = all_data_rnd$Human_mod_index, col_id = 7, x1 = 0.1, x2 = 0.1, y1 = 33, y2 = 50, 
                             legend_title = "Rendezvous site", xtitle = "Human modification index")
  road_den <- density_plots(all_data_den, cov = all_data_den$Nearest_road_m, col_id = 8, x1 = 1300, x2 = 1700, y1 = 0.00065, y2 = 0.00045, 
                            legend_title = "Den", xtitle = "Distance to nearest road (m)")
  road_rnd <- density_plots(all_data_rnd, cov = all_data_rnd$Nearest_road_m, col_id = 8, x1 = 1300, x2 = 1700, y1 = 0.00065, y2 = 0.0005, 
                            legend_title = "Rendezvous site", xtitle = "Distance to nearest road (m)")
  canopy_den <- density_plots(all_data_den, cov = all_data_den$Mean_percent_canopy, col_id = 9, x1 = 0.05, x2 = 0.3, y1 = 6.75, y2 = 1.75, 
                              legend_title = "Den", xtitle = "Average percent canopy cover")
  ndvi_rnd <- density_plots(all_data_rnd, cov = all_data_rnd$meanNDVI, col_id = 11, x1 = 0.01, x2 = 0.20, y1 = 2.9, y2 = 4.9, 
                            legend_title = "Rendezvous site", xtitle = "Average NDVI")
  

  ######  Histograms: used vs available  ######
  ggplot(all_data_den, aes(x = Elevation_m, color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Elevation m", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Den \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Den \nLocations") +
    geom_vline(xintercept = mean(all_data_den$Elevation_m[all_data_den$used == "used"]), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(all_data_den$Elevation_m[all_data_den$used == "available"]), linetype = "dashed", color = "#7570b3")
  ggplot(all_data_rnd, aes(x = Elevation_m, color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Elevation m", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Rnd \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Rnd \nLocations") +
    geom_vline(xintercept = mean(all_data_rnd$Elevation_m[all_data_rnd$used == "used"]), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(all_data_rnd$Elevation_m[all_data_rnd$used == "available"]), linetype = "dashed", color = "#7570b3")
  
  ######  Mapping covariates  ######
  #'  Load boundaries
  buff <- st_read("./Shapefiles/Homesites/Homesite_buffered_MCP_suitableHabitat.shp") %>% st_transform(nad83)
  az_nm <- filter(usa, NAME == "Arizona" | NAME == "New Mexico") %>% st_transform(nad83) 
  mwepa <- st_read("./Shapefiles/MWEPA Layers Zone 1-3 & Boundary/MWEPA Final.shp") %>% st_transform(nad83)
  mwz1 <- st_read("./Shapefiles/MWEPA Layers Zone 1-3 & Boundary/Final_MWEPA_Zone_1.shp") %>% st_transform(nad83)
  mwz2 <- st_read("./Shapefiles/MWEPA Layers Zone 1-3 & Boundary/Final_MWEPA_Zone_2.shp") %>% st_transform(nad83)
  mwz3 <- st_read("./Shapefiles/MWEPA Layers Zone 1-3 & Boundary/Final_MWEPA_Zone_3.shp") %>% st_transform(nad83)
  
  #'  Load additional covaraite rasters
  canopy <- rast("./Shapefiles/Vegetation_variables/percent_canopy_2022_raster.tif")
  #'  Change NAs to 0 (represents canopy loss that year)
  canopy_cov <- canopy[[2]]
  canopy_cov[is.na(canopy_cov[])] <- 0
  ndvi <- rast("./Shapefiles/Vegetation_variables/Resampled_meanNDVI_rnd_movingwindow.tif")
  water_5km <- rast("./Shapefiles/National Hydrography Dataset (NHD)/Dist2Water_5km_res_NAD83.tif")
  curve_5km <- rast("./Shapefiles/Terrain_variables/Gaussian_curvature_5km_res_NAD83.tif")
  roads_5km <- rast("./Shapefiles/Human_variables/dist_to_road_5km_res_NAD83.tif")
  human_mod_5km <- rast("./Shapefiles/Human_variables/global_Human_Modification_5km_res_NAD83.tif")
  
  #'  Visualize
  plot(elev)
  plot(az_nm[1], color = NA, add = T); plot(mwz3, add = T); plot(mwz2, add = T); plot(mwz1, add = T); plot(buff[1], add = T)
  plot(curve_5km)
  plot(az_nm[1], color = NA, add = T); plot(mwz3, add = T); plot(mwz2, add = T); plot(mwz1, add = T); plot(buff[1], add = T)
  
  
  #'  Crop covariates and plot
  crop_covs <- function(cov, shp, mycolors, covname) {
    #'  Reproject shapefile to match raster projection
    shp_reproj <- st_transform(shp, crs(cov))
    #'  Crop raster to shapefile
    crop_cov <- crop(cov, shp_reproj)
    #'  Mask raster to shapefile
    mask_cov <- mask(crop_cov, shp_reproj)
    #'  Plot raster
    mapped_cov <- ggplot() +
      geom_spatraster(data = mask_cov) +
      geom_sf(data = shp_reproj, fill = NA) +
      coord_sf(crs = nad83) +
      scale_fill_gradientn(colours = mycolors, na.value = NA) +
      labs(fill = covname) + xlab("Longitude") + ylab("Latitude") +
      theme_bw() + 
      theme(text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.line = element_line(color = 'black', linewidth = 0.2),
            axis.ticks = element_line(colour = "black", linewidth = 0.2),
            legend.justification = c(1, 0),
            plot.margin = unit(c(0.1,0,0,0.1), "cm"))
    print(mapped_cov)
    return(mapped_cov)
  }
  crop_elev <- crop_covs(elev, shp = buff, covname = "Elevation \n(m)", mycolors = terrain.colors(6))
  crop_slope <- crop_covs(slope, shp = buff, covname = "Slope \n(degrees)", mycolors = terrain.colors(6))
  crop_rough <- crop_covs(rough, shp = buff, covname = "Terrain \nroughness \n(VRM)", mycolors = terrain.colors(6))
  crop_curve <- crop_covs(curve, shp = buff, covname = "Surface \ncurvature", mycolors = c("#26456E", "#2B5C8A", "#508BB8", "#96B7D1", "#F4F9FC", "#FFC382", "#F18932", "#CA5621", "#9E3D22", "#7B3014")) #c("darkblue", "blue", "lightblue", "red", "firebrick", "darkred")
  crop_water <- crop_covs(water, shp = buff, covname = "Nearest \nwater (m)", mycolors = paletteer_c("ggthemes::Classic Blue", 30))
  crop_human <- crop_covs(human_mod, shp = buff, covname = "Human \nmodification \nindex", mycolors = paletteer_c("ggthemes::Orange-Blue Diverging", 30))
  crop_roads <- crop_covs(roads, shp = buff, covname = "Nearest \nroad (m)", mycolors = paletteer_c("ggthemes::Classic Area-Brown", 30))
  crop_canopy <- crop_covs(canopy_cov, shp = buff, covname = "2022 Percent \ncanopy cover", mycolors = paletteer_c("grDevices::Greens", 30, direction = -1))
  crop_ndvi <- crop_covs(ndvi, shp = buff, covname = "2022 Mean \nNDVI", mycolors = paletteer_c("ggthemes::Green-Gold", 30))

  crop_curve_5km <- crop_covs(curve_5km, shp = buff, covname = "Surface \ncurvature", mycolors = c("#26456E", "#2B5C8A", "#508BB8", "#96B7D1", "#F4F9FC", "#FFC382", "#F18932", "#CA5621", "#9E3D22", "#7B3014")) #c("darkblue", "blue", "lightblue", "red", "firebrick", "darkred")
  crop_water_5km <- crop_covs(water_5km, shp = buff, covname = "Nearest \nwater (m)", mycolors = paletteer_c("ggthemes::Classic Blue", 30))
  crop_human_5km <- crop_covs(human_mod_5km, shp = buff, covname = "Human \nmodification \nindex", mycolors = paletteer_c("ggthemes::Orange-Blue Diverging", 30))
  crop_roads_5km <- crop_covs(roads_5km, shp = buff, covname = "Nearest \nroad (m)", mycolors = paletteer_c("ggthemes::Classic Area-Brown", 30))
  
  #'  Customize individual plots and patchwork together
  elev_den <- elev_den + ggtitle("Sampled locations")
  elev_rnd <- elev_rnd + ggtitle("Sampled locations")
  
  crop_elev <- crop_elev + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) #axis.ticks.x = element_blank(), 
  crop_slope <- crop_slope + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) #axis.ticks.x = element_blank(), 
  crop_rough <- crop_rough + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) #axis.ticks.x = element_blank(), 
  crop_human_5km <- crop_human_5km + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) #axis.ticks.x = element_blank(), 
  crop_water <- crop_water + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) #axis.ticks.x = element_blank(), 
  # crop_canopy <- crop_canopy + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) #axis.ticks.x = element_blank(), 
  crop_ndvi <- crop_ndvi + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) #axis.ticks.x = element_blank(), 
  
  ######  Patchwork covariate maps & use/available data  ######
  #'  Significant covaraites from top den RSF
  top_den_covs <- (crop_elev + elev_den) / #+ plot_layout(width = c(1.5, 1)) 
    (crop_slope + slope_den) / 
    (crop_water + water_den) /
    (crop_rough + rough_den) / 
    (crop_roads + road_den) & 
    theme(text = element_text(size = 6),
          #'  Make legend smaller
          legend.key.size = unit(0.5, "lines"),
          legend.justification = c(1, 0))
  #'  Annotate figure 
  top_den_covs[[1]] <- top_den_covs[[1]] + plot_layout(tag_level = 'new')
  top_den_covs[[2]] <- top_den_covs[[2]] + plot_layout(tag_level = 'new')
  top_den_covs[[3]] <- top_den_covs[[3]] + plot_layout(tag_level = 'new')
  top_den_covs[[4]] <- top_den_covs[[4]] + plot_layout(tag_level = 'new')
  top_den_covs[[5]] <- top_den_covs[[5]] + plot_layout(tag_level = 'new')
  top_den_covs <- top_den_covs + plot_annotation(tag_levels = c('a', '1'))
  
  ggsave("./Outputs/Figures/Den_top_mod_signif_covs_all.tiff", top_den_covs, units = "cm", 
         height = 18, width = 11, dpi = 600, device = 'tiff', compression = 'lzw')
  
  #'  Non-significant covariates from top den RSF
  top_den_covs_nonsig <- (crop_human_5km + human_den) / 
    (crop_canopy + canopy_den) & 
    theme(text = element_text(size = 6),
          #'  Make legend smaller
          legend.key.size = unit(0.5, "lines"),
          legend.justification = c(1, 0))
  #'  Annotate figure 
  top_den_covs_nonsig[[1]] <- top_den_covs_nonsig[[1]] + plot_layout(tag_level = 'new')
  top_den_covs_nonsig[[2]] <- top_den_covs_nonsig[[2]] + plot_layout(tag_level = 'new')
  top_den_covs_nonsig <- top_den_covs_nonsig + plot_annotation(tag_levels = c('a', '1'))
  
  ggsave("./Outputs/Figures/Den_top_mod_nonsignif_covs_human_5km.tiff", top_den_covs_nonsig, units = "cm", 
         height = 7, width = 10, dpi = 600, device = 'tiff', compression = 'lzw')
  
  #'  All covariates from top rendezvous site RSF
  top_rnd_covs <- (crop_elev + elev_rnd) / 
    (crop_ndvi + ndvi_rnd) /
    (crop_water + water_rnd) /
    (crop_rough + rough_rnd) / 
    (crop_curve_5km + curve_rnd) &
    theme(text = element_text(size = 6),
          #'  Make legend smaller
          legend.key.size = unit(0.5, "lines"),
          legend.justification = c(1, 0))
  #'  Annotate figure 
  top_rnd_covs[[1]] <- top_rnd_covs[[1]] + plot_layout(tag_level = 'new')
  top_rnd_covs[[2]] <- top_rnd_covs[[2]] + plot_layout(tag_level = 'new')
  top_rnd_covs[[3]] <- top_rnd_covs[[3]] + plot_layout(tag_level = 'new')
  top_rnd_covs[[4]] <- top_rnd_covs[[4]] + plot_layout(tag_level = 'new')
  top_rnd_covs[[5]] <- top_rnd_covs[[5]] + plot_layout(tag_level = 'new')
  top_rnd_covs <- top_rnd_covs + plot_annotation(tag_levels = c('a', '1'))
  
  ggsave("./Outputs/Figures/Rendezvous_top_mod_covs_final.tiff", top_rnd_covs, units = "cm", 
         height = 18, width = 11, dpi = 600, device = 'tiff', compression = 'lzw')
  
  #' #'  Selection of significant covariates from top den RSF
  #' top_den_covs <- (crop_elev + elev_den) / #+ plot_layout(width = c(1.5, 1))
  #'   (crop_slope + slope_den) /
  #'   (crop_water + water_den) &
  #'   theme(text = element_text(size = 6),
  #'         #'  Make legend smaller
  #'         legend.key.size = unit(0.5, "lines"),
  #'         legend.justification = c(1, 0))
  #' #'  Annotate figure
  #' top_den_covs[[1]] <- top_den_covs[[1]] + plot_layout(tag_level = 'new')
  #' top_den_covs[[2]] <- top_den_covs[[2]] + plot_layout(tag_level = 'new')
  #' top_den_covs[[3]] <- top_den_covs[[3]] + plot_layout(tag_level = 'new')
  #' top_den_covs <- top_den_covs + plot_annotation(tag_levels = c('a', '1'))
  #' 
  #' ggsave("./Outputs/Figures/Den_top_mod_signif_covs.tiff", top_den_covs, units = "cm",
  #'        height = 12, width = 12, dpi = 600, device = 'tiff', compression = 'lzw')
  #' 
  #' #'  Additional significant covariates from top den RSF
  #' top_den_covs_b <- (crop_rough + rough_den) /
  #'   (crop_roads + road_den) &
  #'   theme(text = element_text(size = 6),
  #'         #'  Make legend smaller
  #'         legend.key.size = unit(0.5, "lines"),
  #'         legend.justification = c(1, 0))
  #' #'  Annotate figure
  #' top_den_covs_b[[1]] <- top_den_covs_b[[1]] + plot_layout(tag_level = 'new')
  #' top_den_covs_b[[2]] <- top_den_covs_b[[2]] + plot_layout(tag_level = 'new')
  #' top_den_covs_b <- top_den_covs_b + plot_annotation(tag_levels = c('a', '1'))
  #' 
  #' ggsave("./Outputs/Figures/Den_top_mod_signif_covs_secondset.tiff", top_den_covs_b, units = "cm",
  #'        height = 7, width = 11, dpi = 600, device = 'tiff', compression = 'lzw')
  #' 
  #' #'  Selection of significant covariates from top rendezvous site RSF
  #' top_rnd_covs <- (crop_elev + elev_rnd) / 
  #'   (crop_ndvi + ndvi_rnd) /
  #'   (crop_water + water_rnd) &
  #'   theme(text = element_text(size = 6),
  #'         #'  Make legend smaller
  #'         legend.key.size = unit(0.5, "lines"),
  #'         legend.justification = c(1, 0))
  #' #'  Annotate figure 
  #' top_rnd_covs[[1]] <- top_rnd_covs[[1]] + plot_layout(tag_level = 'new')
  #' top_rnd_covs[[2]] <- top_rnd_covs[[2]] + plot_layout(tag_level = 'new')
  #' top_rnd_covs[[3]] <- top_rnd_covs[[3]] + plot_layout(tag_level = 'new')
  #' top_rnd_covs <- top_rnd_covs + plot_annotation(tag_levels = c('a', '1'))
  #' 
  #' ggsave("./Outputs/Figures/Rendezvous_top_mod_signif_covs.tiff", top_rnd_covs, units = "cm", 
  #'        height = 12, width = 12, dpi = 600, device = 'tiff', compression = 'lzw')
  #' 
  #' #'  Selection of nonsignificant covariates from top rendezvous site RSF
  #' top_rnd_covs_nonsig <- (crop_rough + rough_rnd) / 
  #'   (crop_curve_5km + curve_rnd) &
  #'   theme(text = element_text(size = 6),
  #'         #'  Make legend smaller
  #'         legend.key.size = unit(0.5, "lines"),
  #'         legend.justification = c(1, 0))
  #' #'  Annotate figure 
  #' top_rnd_covs_nonsig[[1]] <- top_rnd_covs_nonsig[[1]] + plot_layout(tag_level = 'new')
  #' top_rnd_covs_nonsig[[2]] <- top_rnd_covs_nonsig[[2]] + plot_layout(tag_level = 'new')
  #' top_rnd_covs_nonsig <- top_rnd_covs_nonsig + plot_annotation(tag_levels = c('a', '1'))
  #' 
  #' ggsave("./Outputs/Figures/Rendezvous_top_mod_nonsignif_covs_curve_5km.tiff", top_rnd_covs_nonsig, units = "cm", 
  #'        height = 8, width = 12, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  #'  -----------------------------------------
  ####  Extract covariates for reference grid  ####
  #'  -----------------------------------------
  #'  Load reference grid centroids
  grid_pts <- st_read("./Shapefiles/WMEPA_grid_clip_pts.shp"); crs(grid_pts)
  xy <- st_coordinates(grid_pts)
  
  #' #'  Reproject grid points
  #' grid_pts_wgs84 <- st_transform(grid_pts, crs = wgs84)
  
  #'  Load resampled covariate data
  elev <- terra::rast("./Shapefiles/Terrain_variables/Mosaic_DEM.tif"); res(elev); crs(elev); st_bbox(elev)
  slope <- terra::rast("./Shapefiles/Terrain_variables/slope.tif"); res(slope); crs(slope); st_bbox(slope)
  rough <- terra::rast("./Shapefiles/Terrain_variables/VRM.tif"); res(rough); crs(rough); st_bbox(rough)
  curve <- terra::rast("./Shapefiles/Terrain_variables/Gaussian_curvature.tif"); res(curve); crs(curve); st_bbox(curve)
  water <- terra::rast("./Shapefiles/National Hydrography Dataset (NHD)/Resampled_Dist2Water.tif"); res(water); crs(water); st_bbox(water)
  human_mod <- terra::rast("./Shapefiles/Human_variables/Resampled_global_Human_Modication.tif"); res(human_mod); crs(human_mod); st_bbox(human_mod)
  roads <- terra::rast("./Shapefiles/Human_variables/Resampled_Dist2Road.tif"); res(roads); crs(roads); st_bbox(roads)
  ndvi <- rast("./Shapefiles/Vegetation_variables/Resampled_meanNDVI_rnd_movingwindow.tif"); res(ndvi); crs(ndvi); st_bbox(ndvi)
  canopy <- rast("./Shapefiles/Vegetation_variables/percent_canopy_2022_raster.tif"); res(canopy); crs(canopy); st_bbox(canopy)
  
  #'  Extract covariates
  terrain_stack <- c(elev, slope, rough, curve)
  grid_terrain <- terra::extract(terrain_stack, grid_pts) 
  names(grid_terrain) <- c("ID", "Elevation_m", "Slope_degrees", "Roughness_VRM", "Gaussian_curvature")
  grid_rds <- terra::extract(roads, grid_pts); grid_rds <- grid_rds %>% rename("Nearest_road_m" = "Resampled_Dist2Road")
  grid_gHM <- terra::extract(human_mod, grid_pts); grid_gHM <- grid_gHM %>% rename("Human_mod_index" = "Resampled_global_Human_Modication")
  grid_h20 <- terra::extract(water, grid_pts); grid_h20 <- grid_h20 %>% rename("Nearest_water_m" = "Resampled_Dist2Water")
  grid_ndvi <- terra::extract(ndvi, grid_pts) 
  avg_MWEPA_meanNDVI <- mean(grid_ndvi$Resampled_meanNDVI_rnd_movingwindow)
  grid_ndvi <- grid_ndvi %>% rename("meanNDVI" = "Resampled_meanNDVI_rnd_movingwindow") %>%
    cbind(AvgSeasonalNDVI)
  grid_canopy <- terra::extract(canopy, grid_pts)
  
  
  #'  Combine covariates into single data frame
  grid_covs <- full_join(grid_terrain, grid_h20, by = "ID") %>%
    full_join(grid_gHM, by = "ID") %>%
    full_join(grid_rds, by = "ID") %>% 
    full_join(grid_ndvi, by = "ID") %>% 
    full_join(grid_canopy, by = "ID") 
  summary(grid_covs)
  
  #'  Save
  write_csv(grid_covs, "./Data/MWEPA_suitable_grid_covs.csv")
  
  
  #'  Explore elevation data - are there elevations that should be considered "unavailable"?
  max(grid_covs$Elevation_m)
  summary(grid_covs$Elevation_m)
  quantile(grid_covs$Elevation_m, probs = 0.99995)
  
  