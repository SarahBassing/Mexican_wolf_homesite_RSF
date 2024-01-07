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
  library(spdep)
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
  
  #'  Identify pairwise combinations of sites that are within 250m of each other
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
  
  #'  Assign a confidence weight to sites based on evidence/comments 
  den_weights <- den_clusters %>%
    #'  Highest weight goes to sites that were confirmed or have detailed information about the site (visited)
    mutate(wgts = ifelse(grepl("Confirmed", Comments), 4, NA),
           wgts = ifelse(grepl("natal", Comments), 4, wgts),
           wgts = ifelse(grepl("Natal", Comments), 4, wgts),
           wgts = ifelse(grepl("Den Description", Comments), 4, wgts),
           wgts = ifelse(grepl("Failed", Comments), 4, wgts),
           wgts = ifelse(grepl("failed", Comments), 4, wgts), 
           #'  Middle weight goes to sites with strong evidence, often identified via GPS clustering
           wgts = ifelse(grepl("Strong", Comments), 3, wgts),
           wgts = ifelse(grepl("based on GPS", Comments), 3, wgts),
           wgts = ifelse(grepl("cluster", Comments), 3, wgts),
           wgts = ifelse(Pack_year == "Willow Creek_2022", 3, wgts), # assuming ID'd by GPS locations
           #'  Lowest weight goes to sites with weak evidence or based on approximate triangulation
           wgts = ifelse(grepl("Weak", Comments), 1, wgts),
           wgts = ifelse(Pack_year == "Luna_2019" | Pack_year == "Luna_2020", 1, wgts),
           #'  Assign middle-low weight value for sites with no comments (assuming 
           #'  these sites were identified by GPS)
           wgts = ifelse(is.na(wgts), 2, wgts))
  
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
           #'  Assign middle-low weight value for sites with no comments (assuming 
           #'  these sites were identified by GPS)
           wgts = ifelse(is.na(wgts), 2, wgts))
  
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

  #'  Create single MCP that includes den and rendezvous sites
  used_homesites_sp <- as(used_homesites, "Spatial")
  homesite_mcp <- mcp(used_homesites_sp, percent = 100) 
  homesite_mcp_sf <- st_as_sf(homesite_mcp)
  #'  Hacky way to estimate "radius" of polygon (pretending it's a perfect square)
  (mcp_radius <- as.numeric(sqrt(st_area(homesite_mcp_sf))/2)) # as.numeric(sqrt(st_area(homesite_mcp_sf)/pi)/2) if it was a perfect circle
  #'  Buffer MCP by 50% available habitat extends beyond known used sites
  library(rgeos)
  homesite_mcp_buff <- gBuffer(homesite_mcp, width = mcp_radius)
  # homesite_mcp_buff <- st_buffer(homesite_mcp_sf, mcp_radius)
  # st_area(homesite_mcp_buff)/st_area(homesite_mcp_sf)

  plot(homesite_mcp_buff[1])
  plot(homesite_mcp[1], add = T, col = "red")
  
  
  ##### BUFFERING NOT WORKING  ####
  