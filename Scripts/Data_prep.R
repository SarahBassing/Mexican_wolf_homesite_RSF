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
  
  #'  ------------------------
  ####  Load and format data  ####
  #'  ------------------------
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
  human_mod <- terra::rast("./Shapefiles/Human_variables/mosaic_global_Human_Modification.tif"); res(human_mod); crs(human_mod)
  roads <- terra::rast("./Shapefiles/Human_variables/mosaic_dist2road.tif"); res(roads); crs(roads)
  
  #' #'  Grab dimensions of elev raster
  #' rast_dim <- dim(elev); rast_rows <- rast_dim[1]; rast_cols <- rast_dim[2]
  #' new_cell_id <- seq(1:(rast_rows * rast_cols))
  #' #'  Create empty raster based on extent, resolution, & coord system of elevation raster
  #' empty_rast <- rast(elev, nlyrs = nlyr(elev), names = "CellID", vals = new_cell_id, keeptime = FALSE, 
  #'                    keepunits = FALSE, props = FALSE) 
  #' #'  Make sure it checks out
  #' res(empty_rast); crs(empty_rast); st_bbox(empty_rast); st_bbox(elev)
  #' #'  Save that sucker
  #' writeRaster(empty_rast, filename  = "./Shapefiles/Reference_grid_30m.tif", overwrite = TRUE)
  empty_rast <- terra::rast("./Shapefiles/Reference_grid_30m.tif")
  
  #'  Define WGS84 coordinate systems
  wgs84 <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
  
  #'  Recovery zones and suitable habitat (identified by Martinez-Meyer et al. 2021)
  wmepa <- st_read("./Shapefiles/MWEPA Layers Zone 1-3 & Boundary/MWEPA Final.shp"); crs(wmepa)
  wmz1 <- st_read("./Shapefiles/MWEPA Layers Zone 1-3 & Boundary/Final_MWEPA_Zone_1.shp") %>% st_transform(wgs84)
  wmz2 <- st_read("./Shapefiles/MWEPA Layers Zone 1-3 & Boundary/Final_MWEPA_Zone_2.shp") %>% st_transform(wgs84)
  wmz3 <- st_read("./Shapefiles/MWEPA Layers Zone 1-3 & Boundary/Final_MWEPA_Zone_3.shp") %>% st_transform(wgs84)
  suitable_habitat <- st_read("./Shapefiles/Martinez_Meyer_2021_layers/suitable_habitat.shp") %>% st_transform(wgs84)
  wmepa_wgs84 <- st_transform(wmepa, wgs84)
  crs(wmz1); crs(wmepa_wgs84); crs(suitable_habitat)
  
  #'  Define NAD27 and NAD83 projected coordinate system
  nad27_12N <- st_crs(wmepa)
  nad83 <- st_crs(elev)
  
  #'  Crop suitable habitat to experimental population area and reproject
  wmepa_suitable <- st_intersection(wmepa_wgs84, suitable_habitat)
  plot(wmepa_suitable[1])
  wmepa_suitable_nad27 <- st_transform(wmepa_suitable, crs = nad27_12N)
  wmepa_suitable_nad83 <- st_transform(wmepa_suitable, crs = nad83)
  
  #'  State, highway, and waterbody shapefiles
  usa <- st_read("./Shapefiles/tl_2012_us_state/tl_2012_us_state.shp")
  hwys <- st_read("./Shapefiles/GEE/PrimaryRoads_AZ_NM.shp") %>% st_transform(wgs84)
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
  ggplot(wmepa_wgs84) + geom_sf() + geom_sf(data = wmz1) #'  Recovery zone 1
  ggplot(wmepa_wgs84) + geom_sf() + geom_sf(data = wmz2) #'  Recovery zone 2
  ggplot(wmepa_wgs84) + geom_sf() + geom_sf(data = wmz3) #'  Recovery zone 3
  ggplot(wmepa_wgs84) + geom_sf() + geom_sf(data = wmepa_suitable, fill = "blue") + geom_sf(data = wmz1, fill = "orange")
  
  #'  Recovery zone 1 bounding box
  st_bbox(wmz1)
  wmz1_bbox <- st_as_sfc(st_bbox(wmz1))
  wmz1_extent_wgs84 <- st_transform(wmz1_bbox, wgs84)
  
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
  ggplot(wmz1_bbox) + geom_sf() + geom_sf(data = wmz1) + 
    geom_sf(data = homesites_wgs84[homesites_wgs84$Site_Type == "Den",], aes(color = Year), shape = 16, size = 3) 
  ggplot(wmz1_bbox) + geom_sf() + geom_sf(data = wmz1) + 
    geom_sf(data = homesites_wgs84[homesites_wgs84$Site_Type == "Rendezvous",], aes(color = Year), shape = 16, size = 3) #+
    #' #'  Constrain plot to bbox of the experimental population area
    #' coord_sf(xlim = c(-114.53588, -103.04233), ylim = c(31.96210, 35.53438), expand = FALSE)
    
  #'  Which den site falls way outside experimental population area?
  home_exp_intersection <- st_intersection(homesites_wgs84, wmepa_wgs84)
  ggplot(wmepa_wgs84) + geom_sf() + geom_sf(data = hwys) + geom_sf(data = home_exp_intersection, aes(color = Year), shape = 16, size = 3)
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
  ggplot(wmepa_wgs84) + geom_sf() + 
    geom_sf(data = dens, aes(color = Year), shape = 16) +
    geom_sf(data = homesite_buffers[homesite_buffers$Site_Type == "Den",], colour = "black", fill = NA) +
    coord_sf(xlim = c(-109.8417, -107.3039), ylim = c(33.0215, 34.2875), expand = FALSE) # wmepa_sbb_defined
    #coord_sf(xlim = c(606823.5, 840562.3), ylim = c(3656509.8, 3798894.4), expand = FALSE) # wmepa
  ggplot(wmepa_wgs84) + geom_sf() + 
    geom_sf(data = rnds, aes(color = Year), shape = 16) +
    geom_sf(data = homesite_buffers[homesite_buffers$Site_Type == "Rendezvous",], colour = "black", fill = NA) +
    coord_sf(xlim = c(-109.8417, -107.3039), ylim = c(33.0215, 34.2875), expand = FALSE) # wmepa_sbb_defined
    #coord_sf(xlim = c(606823.5, 840562.3), ylim = c(3656509.8, 3798894.4), expand = FALSE) # wmepa
    
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
  den_sample <- reduce_sites(den_weights) %>%
    dplyr::select(-NewDen_SameYear)
  rnd_sample <- reduce_sites(rnd_weights)

  #'  Merge sampled den and rendezvous sites into single df
  used_homesites <- den_sample %>% bind_rows(rnd_sample)
  used_homesites_nad27 <- st_transform(used_homesites, nad27_12N)
  used_homesites_nad83 <- st_transform(used_homesites, nad83)

  #'  --------------------------------
  ####  Generate available locations  ####
  #'  --------------------------------
  #####  Create buffered sampling extent  #####
  #'  ------------------------------------
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
  
  #'  Mask out large water bodies so not available when drawing random locations
  homesite_mcp_buff_watermask <- st_difference(homesite_mcp_buff, st_union(bigwater_nad27))
  plot(homesite_mcp_buff_watermask[1])
  
  #'  Mask out unsuitable habitat so not available
  #'  First generate polygon of unsuitable habitat
  homesite_mcp_buff_UNsuitablemask <- st_difference(homesite_mcp_buff, st_union(wmepa_suitable_nad27))
  #'  Mask out unsuitable habitat from buffer with large water bodies already masked out
  homesite_mcp_buff_suitablemask <- st_difference(homesite_mcp_buff_watermask, st_union(homesite_mcp_buff_UNsuitablemask))
  plot(homesite_mcp_buff_suitablemask[1])
  
  #'  Visualize (note the coordinate system!)
  #'  100% MCP and buffered MCP
  ggplot(homesite_mcp_buff) + geom_sf() + geom_sf(data = homesite_mcp_sf)
  #'  Den/rendezvous sites within maksed & buffered MCP and Zone 1 for context within Exp. Pop. Area
  ggplot(st_transform(wmepa, nad27_12N)) + geom_sf() + geom_sf(data = homesite_mcp_buff_suitablemask, color = "red") + 
    geom_sf(data = homesite_mcp_sf, color = "blue") + geom_sf(data = st_transform(wmz1, nad27_12N), fill = "gray25", alpha = 0.30) +
    geom_sf(data = homesites_nad27_usa[homesites_nad27_usa$Site_Type == "Den",], aes(color = Year), shape = 16, size = 1.5) +
    ggtitle("Den sites, 100% MCP, and buffered MCP (excluding unsuitable habitat)")
  ggplot(st_transform(wmepa, nad27_12N)) + geom_sf() + geom_sf(data = homesite_mcp_buff_suitablemask, color = "red") + 
    geom_sf(data = homesite_mcp_sf, color = "blue") + geom_sf(data = st_transform(wmz1, nad27_12N), fill = "gray25", alpha = 0.30) +
    geom_sf(data = homesites_nad27_usa[homesites_nad27_usa$Site_Type == "Rendezvous",], aes(color = Year), shape = 16, size = 1.5) +
    ggtitle("Rendezvous sites, 100% MCP, and suitable buffered MCP \n(excluding unsuitable habitat)")
  
  #'  Save as shapefiles
  homesite_mcp_buff_wgs84 <- st_transform(homesite_mcp_buff, wgs84); st_bbox(homesite_mcp_buff_wgs84)
  homesite_mcp_buff_suitablemask_wgs84 <- st_transform(homesite_mcp_buff_suitablemask, wgs84); st_bbox(homesite_mcp_buff_suitablemask_wgs84)
  # st_write(homesite_mcp_buff_wgs84, "./Shapefiles/Homesites/Homesite_buffered_MCP.shp")
  # st_write(homesite_mcp_buff_suitablemask, "./Shapefiles/Homesites/Homesite_buffered_MCP_suitableHabitat.shp")
  
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
    avail_locs_plot <- ggplot(wmepa) + geom_sf() + 
      geom_sf(data = homesite_mcp_buff, color = "red", size = 1.2) + 
      geom_sf(data = wmz1, fill = "gray25", alpha = 0.30) +
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
  
  #'  Save all locations
  st_write(den_locs_wgs84, "./Shapefiles/Homesites/Used_Available_locations_den.shp")
  st_write(rnd_locs_wgs84, "./Shapefiles/Homesites/Used_Available_locations_rnd.shp")
  
  #'  ---------------------
  ####  Gather covariates  ####
  #'  ---------------------
  #'  Read in covariates extracted from Google Earth Engine
  #'  Apply 2000 data to sites older than 2000 and 2022 data to sites from 2023
  #'  due to temporal extent of MODIS and Hansen data
  canopy_den <- read_csv("./Data/GEE extracted data/percent_canopy_denSeason.csv") %>%
    dplyr::select(-site_year)
  canopy_rnd <- read_csv("./Data/GEE extracted data/percent_canopy_rndSeason.csv") %>%
    dplyr::select(-site_year)
  ndvi_den <- read_csv("./Data/GEE extracted data/mean_NDVI_denSeason.csv")
  ndvi_rnd <- read_csv("./Data/GEE extracted data/mean_NDVI_rndSeason.csv")
  
  #'  Create raster stack of all terrain covariates (must be same grid & res)
  terrain_stack <- c(elev, slope, rough, curve) 
  
  #'  Extract covariate values at each used and available location
  get_covs <- function(locs_nad83, locs_wgs84, canopy, ndvi) {
    #'  Extract covariate values at each location (use correct projections!)
    terrain <- terra::extract(terrain_stack, locs_nad83) 
    names(terrain) <- c("ID", "Elevation_m", "Slope_degrees", "Roughness_VRM", "Gaussian_curvature")
    rds <- terra::extract(roads, locs_nad83) %>% rename("Nearest_road_m" = "mosaic_dist2road")
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
  write_csv(all_data_den, "./Data/all_data_den.csv")
  write_csv(all_data_rnd, "./Data/all_data_rnd.csv")
  # write_csv(all_data_den, "./Data/all_data_den_1to1000ratio.csv")
  # write_csv(all_data_rnd, "./Data/all_data_rnd_1to1000ratio.csv")
  
  #'  -----------------------
  #####  Explore covaraites  #####
  #'  -----------------------
  #'  Compare spread of covaraite values between use and available locations
  all_data_den <- mutate(all_data_den, used = ifelse(used == 0, "available", "used"))
  all_data_rnd <- mutate(all_data_rnd, used = ifelse(used == 0, "available", "used"))
  
  ######  Den site histograms  ######
  ggplot(all_data_den, aes(x = Elevation_m, color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Elevation m", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Den \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Den \nLocations") +
    geom_vline(xintercept = mean(all_data_den$Elevation_m[all_data_den$used == "used"]), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(all_data_den$Elevation_m[all_data_den$used == "available"]), linetype = "dashed", color = "#7570b3")
  ggplot(all_data_den, aes(x = Slope_degrees, color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Slope", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Den \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Den \nLocations") +
    geom_vline(xintercept = mean(all_data_den$Slope_degrees[all_data_den$used == "used"]), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(all_data_den$Slope_degrees[all_data_den$used == "available"]), linetype = "dashed", color = "#7570b3")
  ggplot(all_data_den, aes(x = Roughness_VRM, color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Roughness (VRM)", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Den \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Den \nLocations") +
    geom_vline(xintercept = mean(all_data_den$Roughness_VRM[all_data_den$used == "used"]), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(all_data_den$Roughness_VRM[all_data_den$used == "available"]), linetype = "dashed", color = "#7570b3")
  ggplot(all_data_den, aes(x = Gaussian_curvature, color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Curvature", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Den \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Den \nLocations") +
    geom_vline(xintercept = mean(all_data_den$Gaussian_curvature[all_data_den$used == "used"]), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(all_data_den$Gaussian_curvature[all_data_den$used == "available"]), linetype = "dashed", color = "#7570b3")
  ggplot(all_data_den, aes(x = Nearest_road_m, color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Distance to nearest road (m)", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Den \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Den \nLocations") +
    geom_vline(xintercept = mean(all_data_den$Nearest_road_m[all_data_den$used == "used"]), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(all_data_den$Nearest_road_m[all_data_den$used == "available"]), linetype = "dashed", color = "#7570b3")
  ggplot(all_data_den, aes(x = log(Nearest_road_m), color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Log of distance to nearest road (m)", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Den \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Den \nLocations") +
    geom_vline(xintercept = mean(log(all_data_den$Nearest_road_m[all_data_den$used == "used"])), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(log(all_data_den$Nearest_road_m[all_data_den$used == "available"])), linetype = "dashed", color = "#7570b3")
  ggplot(all_data_den, aes(x = Human_mod_index, color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Human modification", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Den \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Den \nLocations") +
    geom_vline(xintercept = mean(all_data_den$Human_mod_index[all_data_den$used == "used"]), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(all_data_den$Human_mod_index[all_data_den$used == "available"]), linetype = "dashed", color = "#7570b3")
  ggplot(all_data_den, aes(x = Nearest_water_m, color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Distance to nearest water (m)", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Den \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Den \nLocations") +
    geom_vline(xintercept = mean(all_data_den$Nearest_water_m[all_data_den$used == "used"]), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(all_data_den$Nearest_water_m[all_data_den$used == "available"]), linetype = "dashed", color = "#7570b3")
  ggplot(all_data_den, aes(x = log(Nearest_water_m), color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Log of distance to nearest water (m)", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Den \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Den \nLocations") +
    geom_vline(xintercept = mean(log(all_data_den$Nearest_water_m[all_data_den$used == "used"])), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(log(all_data_den$Nearest_water_m[all_data_den$used == "available"])), linetype = "dashed", color = "#7570b3")
  ggplot(all_data_den, aes(x = Mean_percent_canopy, color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Avg. percent canopy cover (250 m radius)", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Den \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Den \nLocations") +
    geom_vline(xintercept = mean(all_data_den$Mean_percent_canopy[all_data_den$used == "used"]), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(all_data_den$Mean_percent_canopy[all_data_den$used == "available"]), linetype = "dashed", color = "#7570b3")
  ggplot(all_data_den, aes(x = meanNDVI, color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Avg. seasonal NDVI (250 m radius)", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Den \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Den \nLocations") +
    geom_vline(xintercept = mean(all_data_den$meanNDVI[all_data_den$used == "used"]), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(all_data_den$meanNDVI[all_data_den$used == "available"]), linetype = "dashed", color = "#7570b3")
  
  ######  Rendezvous site histograms  ######
  ggplot(all_data_rnd, aes(x = Elevation_m, color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Elevation m", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Rnd \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Rnd \nLocations") +
    geom_vline(xintercept = mean(all_data_rnd$Elevation_m[all_data_rnd$used == "used"]), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(all_data_rnd$Elevation_m[all_data_rnd$used == "available"]), linetype = "dashed", color = "#7570b3")
  ggplot(all_data_rnd, aes(x = Slope_degrees, color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Slope", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Rnd \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Rnd \nLocations") +
    geom_vline(xintercept = mean(all_data_rnd$Slope_degrees[all_data_rnd$used == "used"]), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(all_data_rnd$Slope_degrees[all_data_rnd$used == "available"]), linetype = "dashed", color = "#7570b3")
  ggplot(all_data_rnd, aes(x = Roughness_VRM, color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Roughness (VRM)", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Rnd \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Rnd \nLocations") +
    geom_vline(xintercept = mean(all_data_rnd$Roughness_VRM[all_data_rnd$used == "used"]), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(all_data_rnd$Roughness_VRM[all_data_rnd$used == "available"]), linetype = "dashed", color = "#7570b3")
  ggplot(all_data_rnd, aes(x = Gaussian_curvature, color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Curvature", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Rnd \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Rnd \nLocations") +
    geom_vline(xintercept = mean(all_data_rnd$Gaussian_curvature[all_data_rnd$used == "used"]), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(all_data_rnd$Gaussian_curvature[all_data_rnd$used == "available"]), linetype = "dashed", color = "#7570b3")
  ggplot(all_data_rnd, aes(x = log(Nearest_road_m), color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Log of distance to nearest road (m)", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Rnd \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Rnd \nLocations") +
    geom_vline(xintercept = mean(log(all_data_rnd$Nearest_road_m[all_data_rnd$used == "used"])), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(log(all_data_rnd$Nearest_road_m[all_data_rnd$used == "available"])), linetype = "dashed", color = "#7570b3")
  ggplot(all_data_rnd, aes(x = Human_mod_index, color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Human modification", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Rnd \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Rnd \nLocations") +
    geom_vline(xintercept = mean(all_data_rnd$Human_mod_index[all_data_rnd$used == "used"]), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(all_data_rnd$Human_mod_index[all_data_rnd$used == "available"]), linetype = "dashed", color = "#7570b3")
  ggplot(all_data_rnd, aes(x = log(Nearest_water_m), color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Log of distance to nearest water (m)", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Rnd \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Rnd \nLocations") +
    geom_vline(xintercept = mean(log(all_data_rnd$Nearest_water_m[all_data_rnd$used == "used"])), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(log(all_data_rnd$Nearest_water_m[all_data_rnd$used == "available"])), linetype = "dashed", color = "#7570b3")
  ggplot(all_data_rnd, aes(x = Mean_percent_canopy, color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Avg. percent canopy cover (250 m radius)", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Rnd \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Rnd \nLocations") +
    geom_vline(xintercept = mean(all_data_rnd$Mean_percent_canopy[all_data_rnd$used == "used"]), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(all_data_rnd$Mean_percent_canopy[all_data_rnd$used == "available"]), linetype = "dashed", color = "#7570b3")
  ggplot(all_data_rnd, aes(x = meanNDVI, color = used, fill = used)) + 
    geom_histogram(alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Used vs Available Locations", x = "Avg. seasonal NDVI (250 m radius)", y = "Frequency") + 
    scale_fill_manual(values = c("available" = "#7570b3","used" = "#d95f02"), name = "Rnd \nLocations") +
    scale_color_manual(values = c("available" = "#7570b3", "used" = "#d95f02"), name = "Rnd \nLocations") +
    geom_vline(xintercept = mean(all_data_rnd$meanNDVI[all_data_rnd$used == "used"]), linetype = "dashed", color = "#d95f02") +
    geom_vline(xintercept = mean(all_data_rnd$meanNDVI[all_data_rnd$used == "available"]), linetype = "dashed", color = "#7570b3")
  
  #'  -----------------------------------------
  ####  Extract covariates for reference grid  ####
  #'  -----------------------------------------
  #' #'  Mask out waterbodies
  #' masked_water <- mask(empty_rast, bigwater_nad83, inverse = TRUE)
  #' plot(masked_water)
  #' #'  Mask out unsuitable areas (so areas NOT in the wmep_suitable polygon)
  #' wmepa_grid <- mask(masked_water, wmepa_suitable_nad83)
  #' plot(wmepa_grid); plot(homesites_nad83, add = TRUE)
  #' #'  Save
  #' writeRaster(wmepa_grid, "./Shapefiles/WMEPA_masked_grid.tif", overwrite = TRUE)
  
  #'  Convert masked wmepa_grid to a polygon
  wmepa_grid <- terra::rast("./Shapefiles/WMEPA_masked_grid.tif")
  wmepa_poly <- as.polygons(wmepa_grid)
  wmepa_poly_sf <- st_as_sf(wmepa_poly)
  st_write(wmepa_poly_sf, "./Shapefiles/WMEPA_masked_polygon.shp")
  
  
  #'  Crop to just buffered MCP area
  
  