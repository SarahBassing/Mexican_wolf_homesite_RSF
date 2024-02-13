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
  library(spdep)
  library(adehabitatHR)
  library(tidyverse)
  
  #'  -------------
  ####  Load data  ####
  #'  -------------
  #'  Used locations
  homesites <- read_csv("./Data/MexWolf_dens_rend_sites_1998_2023_updated_01.19.24.csv") %>%
    mutate(Pack_year = paste0(Pack, "_", Year)) %>%
    relocate(Pack_year, .before = Year)
    
  #'  Spatial data (rasters)
  elev <- terra::rast("./Shapefiles/Terrain_variables/Mosaic_DEM.tif"); res(elev); crs(elev)
  slope <- terra::rast("./Shapefiles/Terrain_variables/slope.tif"); res(slope); crs(slope)
  rough <- terra::rast("./Shapefiles/Terrain_variables/VRM.tif"); res(rough); crs(rough)
  curve <- terra::rast("./Shapefiles/Terrain_variables/Gaussian_curvature.tif"); res(curve); crs(curve)
  water <- terra::rast("./Shapefiles/National Hydrography Dataset (NHD)/Mosaic_Dist2Water.tif"); res(water); crs(water)
  ndvi_den <- terra::rast("./Shapefiles/Vegetation_variables/Mosaic_Avg_NDVI_Mar_June.tif"); res(ndvi_den); crs(ndvi_den)
  ndvi_rnd <- terra::rast("./Shapefiles/Vegetation_variables/Mosaic_Avg_NDVI_June_Aug.tif"); res(ndvi_rnd); crs(ndvi_rnd)
  ndvi_den0 <- terra::rast("./Shapefiles/Vegetation_variables/meanNDVI_den_2000.tif"); res(ndvi_den0); crs(ndvi_den0)
  
  #' #'  Empty raster to extract covariates for each pixel for predicting phase of study
  #' s <- elev
  #' values(s) <- NA
  #' s
  #' s_gridID <- as.data.frame(s)
  #' s_coord <- coordinates(s)
  #' s_grid <- cbind(s_gridID, s_coord)
  #' #'  Mask large bodies of water from raster
  #' s_mask <- mask(s, bigwater, inverse = TRUE)
  
  #'  Define WGS84 coordinate systems
  wgs84 <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
  
  #'  Recovery zones
  exp_pop <- st_read("./Shapefiles/MWEPA Layers Zone 1-3 & Boundary/MWEPA Final.shp"); crs(exp_pop)
  zone1 <- st_read("./Shapefiles/MWEPA Layers Zone 1-3 & Boundary/Final_MWEPA_Zone_1.shp") %>% st_transform(wgs84)
  zone2 <- st_read("./Shapefiles/MWEPA Layers Zone 1-3 & Boundary/Final_MWEPA_Zone_2.shp") %>% st_transform(wgs84)
  zone3 <- st_read("./Shapefiles/MWEPA Layers Zone 1-3 & Boundary/Final_MWEPA_Zone_3.shp") %>% st_transform(wgs84)
  exp_pop_wgs84 <- st_transform(exp_pop, wgs84)
  crs(zone1); crs(exp_pop_wgs84)
  
  #'  Define NAD27 and NAD83 projected coordinate system
  nad27_12N <- st_crs(exp_pop)
  nad83 <- st_crs(elev)
  
  #'  State, highway, and waterbody shapefiles
  usa <- st_read("./Shapefiles/tl_2012_us_state/tl_2012_us_state.shp")
  hwys <- st_read("./Shapefiles/GEE/PrimaryRoads_AZ_NM.shp") %>% st_transform(wgs84)
  az_nm <- filter(usa, NAME == "Arizona" | NAME == "New Mexico") %>% st_transform(wgs84) 
  crs(az_nm); st_bbox(az_nm)
  #'  Merge into single "southwest" polygon
  southwest <- st_union(az_nm[az_nm$NAME == "Arizona",], az_nm[az_nm$NAME == "New Mexico",]) %>%
    st_cast("POLYGON")
  water <- st_read("./Shapefiles/National Hydrography Dataset (NHD)/AZ_NM_waterbodies_NHD.shp"); crs(water)
  #'  Identify large bodies of water (anything larger than 1 sq-km in size)
  bigwater <- water[water$areasqkm > 1,]
  bigwater_nad27 <- st_transform(bigwater, crs = nad27_12N)
  #' #'  Create a single LINESTRING for I40 and I10
  #' I40 <- filter(hwys, fullname == "I- 40") %>%
  #'   st_union()
  #' I10 <- filter(hwys, fullname == "I- 10") %>%
  #'   st_union()
  #' #'  Join highways that border experimental population area
  #' I40_I10 <- st_union(I40, I10)
  #' #'  Split southwest polygon by I40
  #' split_southwest <- lwgeom::st_split(southwest, I40_I10) %>%
  #'   st_collection_extract("POLYGON")
  #' split_southwest$include <- factor(seq_len(nrow(split_southwest)))
  #' ggplot(split_southwest) + geom_sf(aes(fill = include))
  #' #'  Filter to SBB defined experimental population area and anything south of I40
  #' exp_pop_sbb_defined <- filter(split_southwest, include == 2); st_crs(exp_pop_sbb_defined); st_bbox(exp_pop_sbb_defined)
  #' southI40 <- filter(split_southwest, include != 1) %>% st_union(); st_crs(southI40); st_bbox(southI40)
  #' ggplot(exp_pop_sbb_defined) + geom_sf();  ggplot(az_nm) + geom_sf() + geom_sf(data = southI40)
  
  #'  Review each zones within context of larger experimental population area boundary
  ggplot(exp_pop) + geom_sf() + geom_sf(data = zone1) #'  Recovery zone 1
  ggplot(exp_pop) + geom_sf() + geom_sf(data = zone2) #'  Recovery zone 2
  ggplot(exp_pop) + geom_sf() + geom_sf(data = zone3) #'  Recovery zone 3
  
  #'  Recovery zone 1 bounding box
  st_bbox(zone1)
  zone1_bbox <- st_as_sfc(st_bbox(zone1))
  zone1_extent_wgs84 <- st_transform(zone1_bbox, wgs84)
  
  #'  Save select features
  # st_write(az_nm, "./Shapefiles/tl_2012_us_state/Arizona_NewMexico.shp")
  # st_write(az_nm, "./Shapefiles/tl_2012_us_state/Arizona_NewMexico.kml", driver = "kml", delete_dsn = TRUE)
  # st_write(exp_pop_sbb_defined, "./Shapefiles/experimental_pop_poly.shp")
  # st_write(exp_pop_wgs84, "./Shapefiles/experimental_pop_poly_wgs84.shp")
  # st_write(southI40, "./Shapefiles/I40_south_poly.shp")
  # st_write(zone1_extent_wgs84, "./Shapefiles/RecoveryZone1_bbox_wgs84.shp")
  
  #'  Create a sf object for locations
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
  ggplot(zone1_bbox) + geom_sf() + geom_sf(data = zone1) + 
    geom_sf(data = homesites_wgs84[homesites_wgs84$Site_Type == "Den",], aes(color = Year), shape = 16, size = 3) 
  ggplot(zone1_bbox) + geom_sf() + geom_sf(data = zone1) + 
    geom_sf(data = homesites_wgs84[homesites_wgs84$Site_Type == "Rendezvous",], aes(color = Year), shape = 16, size = 3) #+
    #' #'  Constrain plot to bbox of the experimental population area
    #' coord_sf(xlim = c(-114.53588, -103.04233), ylim = c(31.96210, 35.53438), expand = FALSE)
    
  #'  Which den site falls way outside experimental population area?
  home_exp_intersection <- st_intersection(homesites_wgs84, exp_pop_wgs84)
  ggplot(exp_pop_wgs84) + geom_sf() + geom_sf(data = hwys) + geom_sf(data = home_exp_intersection, aes(color = Year), shape = 16, size = 3)
  far_away_home <- subset(homesites_wgs84, !(obs %in% home_exp_intersection$obs))
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
  
  #'  Save pack-years with >1 den site
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
  ggplot(exp_pop_wgs84) + geom_sf() + 
    geom_sf(data = dens, aes(color = Year), shape = 16) +
    geom_sf(data = homesite_buffers[homesite_buffers$Site_Type == "Den",], colour = "black", fill = NA) +
    coord_sf(xlim = c(-109.8417, -107.3039), ylim = c(33.0215, 34.2875), expand = FALSE) # exp_pop_sbb_defined
    #coord_sf(xlim = c(606823.5, 840562.3), ylim = c(3656509.8, 3798894.4), expand = FALSE) # exp_pop
  ggplot(exp_pop_wgs84) + geom_sf() + 
    geom_sf(data = rnds, aes(color = Year), shape = 16) +
    geom_sf(data = homesite_buffers[homesite_buffers$Site_Type == "Rendezvous",], colour = "black", fill = NA) +
    coord_sf(xlim = c(-109.8417, -107.3039), ylim = c(33.0215, 34.2875), expand = FALSE) # exp_pop_sbb_defined
    #coord_sf(xlim = c(606823.5, 840562.3), ylim = c(3656509.8, 3798894.4), expand = FALSE) # exp_pop
    
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
           wgts = ifelse(grepl("Den Description", Comments), 4, wgts),
           wgts = ifelse(grepl("Failed", Comments), 4, wgts),
           wgts = ifelse(grepl("failed", Comments), 4, wgts), 
           #'  Downgrade sites with den descriptions that were not the natal den 
           #'  (moved from original den so selection process was slightly different)
           wgts = ifelse(Pack_year == "Iron Creek_2015" & !grepl("Natal den", Comments), 3, wgts),
           wgts = ifelse(Pack_year == "Maverick_2014" & !grepl("natal", Comments), 3, wgts),
           wgts = ifelse(Pack_year == "Rim_2014" & !grepl("Natal den", Comments), 3, wgts),
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
  
  #'  Randomly select one site per pack per cluster to keep for RSF analyses
  #'  Need to exclude repeat use of sites by the same pack to keep each "used" 
  #'  observation in the RSF as independent as possible (i.e., some packs re-use
  #'  same den or rendezvous every year - multiple used points at this location 
  #'  is not representative of what the population is selecting for)
  reduce_sites <- function(sites) {
    #'  Set seed so sampling is reproducible
    set.seed(2024)
    site_random_sample <- sites %>%
      #'  Grouping by pack and cluster allows same site to be used by different
      #'  packs but prevents same site from being used by one pack multiple times
      group_by(Pack, cluster_id) %>%
      #'  Weight observations by strength of confidence and draw 1 sample
      slice_sample(n = 1, weight_by = wgts) %>%
      ungroup()
    return(site_random_sample)
  }
  den_sample <- reduce_sites(den_weights)
  rnd_sample <- reduce_sites(rnd_weights)

  #'  Merge sampled den and rendezvous sites into single df
  used_homesites <- den_sample %>%
    dplyr::select(-NewDen_SameYear) %>%
    bind_rows(rnd_sample)
  used_homesites_nad27 <- st_transform(used_homesites, nad27_12N)

  #'  --------------------------------
  ####  Generate available locations  ####
  #'  --------------------------------
  #'  Define extent of "available" habitat for wolves to den/rendezvous in
  #'  Create single MCP that includes den and rendezvous sites
  #'  NOTE the coordinate sytem - buffer is smoother in projected coord system
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
  
  #'  Mask out large waterbodies so not available when drawing random locations
  homesite_mcp_buff_watermask <- st_difference(homesite_mcp_buff, st_union(bigwater_nad27))
  
  
  #### MASK OUT LOW ELEVATION OR UNSUITABLE HABITAT & URBAN CENTERS (CREATE A SHAPEFILE FIRST)  ####
  
  
  #'  Visualize (note the coordinate system!)
  #'  100% MCP and buffered MCP
  ggplot(homesite_mcp_buff) + geom_sf() + geom_sf(data = homesite_mcp_sf)
  #'  Den/rendezvous sites within buffered MCP and Zone 1 for context within Exp. Pop. Area
  ggplot(st_transform(exp_pop, nad27_12N)) + geom_sf() + geom_sf(data = homesite_mcp_buff, color = "red") + 
    geom_sf(data = homesite_mcp_sf, color = "blue") + geom_sf(data = st_transform(zone1, nad27_12N), fill = "gray25", alpha = 0.30) +
    geom_sf(data = homesites_nad27_usa[homesites_nad27_usa$Site_Type == "Den",], aes(color = Year), shape = 16, size = 1.5) +
    ggtitle("Mexican wolf den sites, 100% MCP, and buffered MCP")
  ggplot(st_transform(exp_pop, nad27_12N)) + geom_sf() + geom_sf(data = homesite_mcp_buff, color = "red") + 
    geom_sf(data = homesite_mcp_sf, color = "blue") + geom_sf(data = st_transform(zone1, nad27_12N), fill = "gray25", alpha = 0.30) +
    geom_sf(data = homesites_nad27_usa[homesites_nad27_usa$Site_Type == "Rendezvous",], aes(color = Year), shape = 16, size = 1.5) +
    ggtitle("Mexican wolf rendezvous sites, 100% MCP, and buffered MCP")
  
  #' #'  Save as shapefile and kml
  #' homesite_mcp_buff_wgs84 <- st_transform(homesite_mcp_buff, wgs84); st_bbox(homesite_mcp_buff_wgs84)
  #' st_write(homesite_mcp_buff_wgs84, "./Shapefiles/Homesites/Homesite_buffered_MCP.shp")
  #' st_write(homesite_mcp_buff_wgs84, "./Shapefiles/Homesites/Homesite_buffered_MCP.kml", driver = "kml", delete_dsn = TRUE)
  
  #'  Number of available locations to generate per used location 
  avail_pts <- 1000
  
  #'  Function to generate random available locations based on number of used locations
  sample_avail_locs <- function(locs, newcrs, navail) {
    #'  Reproject locations
    locs <- st_transform(locs, newcrs)
    
    #'  Number of used locations
    nused <- nrow(locs)
    print(nused)
    
    #'  Total number of available locations to generate
    navailable <- nused * navail
    
    #'  Set seed for reproducibility
    set.seed(108)
    
    #'  Draw random sample of locations within the buffered MCP (excluding large waterbodies)
    rndpts <- st_sample(homesite_mcp_buff_watermask, size = navailable, type = "random", exact = TRUE) %>%
      #'  Reformat to a normal sf object
      st_as_sf() %>%
      mutate(obs = seq(1:nrow(.))) %>%
      relocate(obs, .before = x) %>%
      rename(geometry = x)
    
    #'  Plot available points within buffered MCP
    avail_locs_plot <- ggplot(exp_pop) + geom_sf() + 
      geom_sf(data = homesite_mcp_buff, color = "red") + 
      geom_sf(data = zone1, fill = "gray25", alpha = 0.30) +
      geom_sf(data = rndpts, shape = 16, size = 0.5, color = "blue") 
    plot(avail_locs_plot)
  
    return(rndpts)
  }
  avail_locs_den <- sample_avail_locs(den_sample, newcrs = nad27_12N, navail = avail_pts)
  avail_locs_rnd <- sample_avail_locs(rnd_sample, newcrs = nad27_12N, navail = avail_pts)
  
  #'  Reproject available locations
  avail_locs_den_wgs84 <- st_transform(avail_locs_den, wgs84)
  avail_locs_rnd_wgs84 <- st_transform(avail_locs_rnd, wgs84)
  
  #'  Save available locations
  # st_write(avail_locs_den_wgs84, "./Shapefiles/Homesites/Available_locations_den.shp")
  # st_write(avail_locs_den_wgs84, "./Shapefiles/Homesites/Available_locations_den.kml", driver = "kml", delete_dsn = TRUE)
  # st_write(avail_locs_rnd_wgs84, "./Shapefiles/Homesites/Available_locations_rnd.shp")
  # st_write(avail_locs_rnd_wgs84, "./Shapefiles/Homesites/Available_locations_rnd.kml", driver = "kml", delete_dsn = TRUE)
  
  
  #'  ---------------------
  ####  Gather covariates  ####
  #'  ---------------------
  #'  Create raster stacks of all relevant covariates (must be same grid & res)
  terrain_stack <- c(elev, slope, rough, curve) 
  ndvi_stack <- c(ndvi_den, ndvi_rnd) #water, 
  
  #'  DON'T FORGET TO LOAD WATERBODIES
  #'  Identify large bodies of water (anything larger than 1 sq-km in size)
  bigwater <- waterbody[waterbody$AreSqKm > 1,]
  #'  Mask large bodies of water from raster
  grid_mask <- mask(s, bigwater, inverse = TRUE)
  
  #'  Extract covariate values at each used and available location
  get_covs <- function(locs) {
    tst <- terra::extract(locs, terrain_stack)
  }
  
  
  
  