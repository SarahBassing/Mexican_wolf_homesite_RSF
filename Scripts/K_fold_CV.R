  #'  -------------------
  #'  K-fold Cross-Validation for RSFs
  #'  Mexican wolf project
  #'  Sarah B. Bassing
  #'  March 2024
  #'  -------------------
  #'  Script to run k-fold cross-validation for most supported RSFs following
  #'  methods described by Boyce et al. (2002). Script pulls in used/available 
  #'  location data, partitions data into K fold trainign and testing data sets,
  #'  trains models, predicts trained model results and then tests predictions
  #'  using the training data.
  #'  -------------------
  
  #'  Clear memory
  rm(list = ls())
  
  #'  Load packages
  library(car)
  library(groupdata2)
  library(sf)
  library(terra)
  library(stars)
  library(ggplot2)
  library(tidyverse)
  
  #'  Load used/available location data and covariates
  all_data_den <- read_csv("./Data/all_data_den.csv")
  all_data_rnd <- read_csv("./Data/all_data_rnd.csv")
  
  #'  Load MWEPA masked grid and covariate data
  grid_covs <- read_csv("./Data/MWEPA_suitable_grid_covs.csv")
  
  #'  --------------------------------
  ####  K-fold training/testing data  ####
  #'  --------------------------------
  #'  Partition each data set into K-folds
  fold_in_the_cheese <- function(dat, mu.sd, K) {
    set.seed(2024)
    #'  Use groupdata2 package with cat_col = "used" to balance folds proportional
    #'  to 0's and 1's in each data set
    fold_df <- fold(dat, k = K, cat_col = "used")
    #'  Turn into a dataframe
    fold_df <- as.data.frame(fold_df)
    
    #'  Create K training and testing data sets
    training_sets <- list()
    testing_sets <- list()
    for(i in 1:K) {
      training_sets[[i]] <- fold_df[fold_df$.folds != i,]
      testing_sets[[i]] <- fold_df[fold_df$.folds == i,]
    }
    
    #'  Double check these training sets split correctly
    print("Training data: all folds but #1?")
    print(unique(training_sets[[1]][".folds"])) # should have folds 2+
    print("Training data: all folds but #2?")
    print(unique(training_sets[[2]][".folds"])) # should have all folds BUT 2
    
    print("Testing data: only fold #1?")
    print(unique(testing_sets[[1]][".folds"])) # should have fold 1 only
    print("Testing data: only fold #2?")
    print(unique(testing_sets[[2]][".folds"])) # should have fold 2 only
    
    #'  List testing and training lists together
    k_folded_dat <- list(training_sets, testing_sets)
    return(k_folded_dat)
  }
  data_den_k <- fold_in_the_cheese(all_data_den, K = 10)
  data_rnd_k <- fold_in_the_cheese(all_data_rnd, K = 10)
  
  #'  Standardize covariates for each training data set using mean & SD of original data
  standardize_training_covs <- function(dat) {
    #'  Center and scale covariates
    dat_z <- dat %>%
      transmute(Pack_year = Pack_year,
                Homesite = Homesite,
                used = as.numeric(used),
                wgts = as.numeric(wgts),
                Elev = scale(as.numeric(Elevation_m)),
                Slope = scale(as.numeric(Slope_degrees)),
                Rough = scale(as.numeric(Roughness_VRM)),
                Curve = scale(as.numeric(Gaussian_curvature)),
                Dist2Water = scale(as.numeric(Nearest_water_m)),
                logDist2Water = scale(log(as.numeric(Nearest_water_m))),
                HumanMod = scale(as.numeric(Human_mod_index)),
                Dist2Road = scale(as.numeric(Nearest_road_m)),
                logDist2Road = scale(log(as.numeric(Nearest_road_m))),
                CanopyCov = scale(as.numeric(Mean_percent_canopy)),
                AvgCanopyCov = scale(as.numeric(avg_MCP_canopycover)),
                SeasonalNDVI = scale(as.numeric(meanNDVI)),
                AvgSeasonalNDVI = scale(as.numeric(avg_MCP_meanNDVI)))
    return(dat_z)
  } 
  training_den <- lapply(data_den_k[[1]], standardize_training_covs)  #'  List #1 contains lists of training data
  training_rnd <- lapply(data_rnd_k[[1]], standardize_training_covs)  #'  List #1 contains lists of training data
  
  #'  -------------------------------
  ####  Train RSFs on training data  ####
  #'  -------------------------------
  #'  Top den model
  h4.den <- "used ~ Elev + Slope + Rough + Dist2Water + CanopyCov + CanopyCov:AvgCanopyCov + HumanMod + Dist2Road"
  #'  Top rendezvous site model
  h2.rnd <- "used ~ Elev + Rough + Curve + Dist2Water + SeasonalNDVI + SeasonalNDVI:AvgSeasonalNDVI"
  
  #'  Refit model K-fold times
  train_rsf <- function(dat, mod) {
    #'  Fit the model
    mod <- glm(formula = mod, data = dat, weight = wgts, family = binomial) 
    #'  Inspect the model
    print(summary(mod))
    print(car::vif(mod))
    return(mod)
  }
  trained_den_k <- lapply(training_den, train_rsf, mod = h4.den)
  trained_rnd_k <- lapply(training_rnd, train_rsf, mod = h2.rnd)
  
  #'  -------------------------------------
  ####  Predict trained data across MWEPA  ####
  #'  -------------------------------------
  #'  Grab mean and SD of each covariate from original data set
  cov_summary_stats <- function(dat) {
    #'  Drop site identifying info
    dat <- dplyr::select(dat, -c("ID", "Pack_year", "Homesite", "used", "wgts", "geometry"))
    #'  Calculate mean and SD of each covariate
    cov_means <- dat %>%
      summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))
    cov_sd <- dat %>%
      summarise(across(where(is.numeric), \(x) sd(x, na.rm = TRUE)))
    #'  Bind and return
    cov_stats <- bind_rows(cov_means, cov_sd) 
    cov_stats <- as.data.frame(cov_stats)
    row.names(cov_stats) <- c("Mean", "SD")
    return(cov_stats)
  }
  data_den_stats <- cov_summary_stats(all_data_den)
  data_rnd_stats <- cov_summary_stats(all_data_rnd)
  
  #'  Standardize based on mean and SD of covariates in original model
  standardize_mwepa_covs <- function(dat, mu.sd) {
    zcovs <- dat %>%
      transmute(cellID = cellID,
                ID = ID,
                Elev = (Elevation_m - mu.sd$Elevation_m[1])/mu.sd$Elevation_m[2],
                Slope = (Slope_degrees - mu.sd$Slope_degrees[1])/mu.sd$Slope_degrees[2],
                Rough = (Roughness_VRM - mu.sd$Roughness_VRM[1])/mu.sd$Roughness_VRM[2],
                Curve = (Gaussian_curvature - mu.sd$Gaussian_curvature[1])/mu.sd$Gaussian_curvature[2],
                Dist2Water = (Nearest_water_m - mu.sd$Nearest_water_m[1])/mu.sd$Nearest_water_m[2],
                # logDist2Water = scale(log(as.numeric(Nearest_water_m))),
                HumanMod = (Human_mod_index - mu.sd$Human_mod_index[1])/mu.sd$Human_mod_index[2],
                Dist2Road = (Nearest_road_m - mu.sd$Nearest_road_m[1])/mu.sd$Nearest_road_m[2],
                # logDist2Road = scale(log(as.numeric(Nearest_road_m))),
                CanopyCov = (Mean_percent_canopy - mu.sd$Mean_percent_canopy[1])/mu.sd$Mean_percent_canopy[2],
                AvgCanopyCov = (avg_MWEPA_canopycover - mu.sd$avg_MCP_canopycover[1])/mu.sd$avg_MCP_canopycover[2],
                ModifiedCanopyCov = CanopyCov*AvgCanopyCov,
                SeasonalNDVI = (meanNDVI - mu.sd$meanNDVI[1])/mu.sd$meanNDVI[2],
                AvgSeasonalNDVI = (avg_MCP_meanNDVI - mu.sd$avg_MCP_meanNDVI[1])/mu.sd$avg_MCP_meanNDVI[2],
                ModifiedNDVI = SeasonalNDVI * AvgSeasonalNDVI,
                x = as.numeric(X),
                y = as.numeric(Y))
    return(zcovs)
  }
  zcovs_den_mwepa <- standardize_mwepa_covs(grid_covs, mu.sd = data_den_stats)
  zcovs_rnd_mwepa <- standardize_mwepa_covs(grid_covs, mu.sd = data_rnd_stats)
  
  #'  Function to save parameter estimates from each trained model
  #'  Use coef(mod) to look at random effects estimates
  rsf_out <- function(mod){
    beta <- mod$coefficients
    out <- as.data.frame(beta) %>%
      mutate(Parameter = row.names(.),
             Parameter = ifelse(Parameter == "(Intercept)", "alpha", Parameter),
             Parameter = ifelse(Parameter == "Elev", "b.elev", Parameter),
             Parameter = ifelse(Parameter == "Slope", "b.slope", Parameter),
             Parameter = ifelse(Parameter == "Rough", "b.rough", Parameter),
             Parameter = ifelse(Parameter == "Curve", "b.curve", Parameter),
             Parameter = ifelse(Parameter == "Dist2Water", "b.water", Parameter),
             Parameter = ifelse(Parameter == "HumanMod", "b.hm", Parameter),
             Parameter = ifelse(Parameter == "Dist2Road", "b.road", Parameter),
             Parameter = ifelse(Parameter == "CanopyCov", "b.canopy", Parameter),
             Parameter = ifelse(Parameter == "CanopyCov:AvgCanopyCov", "b.canopyXavgcanopy", Parameter),
             Parameter = ifelse(Parameter == "SeasonalNDVI", "b.ndvi", Parameter),
             Parameter = ifelse(Parameter == "SeasonalNDVI:AvgSeasonalNDVI", "b.ndviXavgndvi", Parameter)) %>%
      relocate(Parameter, .before = beta) %>%
      #'  Spread data so each coefficient is it's own column
      pivot_wider(names_from = Parameter, values_from = beta) 
  
    return(out)
  }
  #'  Extract coefficient estimates for each trained model
  trained_den_k_coefs <- lapply(trained_den_k, rsf_out)
  trained_rnd_k_coefs <- lapply(trained_rnd_k, rsf_out)
  
  #'  Functions to predict across all grid cells based on RSF results
  #'  Should end up with 1 predicted value per grid cell
  #'  NOTE: want the predict relative probability of selection from RSF so not 
  #'  using a logit transformation. Drop intercept from the model and exponentiate 
  #'  coefs*covs (Fieberg et al. 2020)
  predict_den_rsf <- function(coef, cov) {
    #'  Generate modified canopy cover variable
    
    predict_rsf <- c()
    #'  Predict across each grid cell
    for(i in 1:nrow(cov)) {
      predict_rsf[i] <- exp(coef$b.elev*cov$Elev[i] + coef$b.slope*cov$Slope[i] + 
                              coef$b.rough*cov$Rough[i] + coef$b.water*cov$Dist2Water[i] + 
                              coef$b.canopy*cov$CanopyCov[i] + 
                              coef$b.canopyXavgcanopy*cov$ModifiedCanopyCov[i] + 
                              coef$b.hm*cov$HumanMod[i] + coef$b.road*cov$Dist2Road[i])}  
    predict_rsf <- as.data.frame(predict_rsf)
    predict_rsf <- cbind(cov$ID, cov$x, cov$y, predict_rsf)
    colnames(predict_rsf) <- c("ID", "x", "y", "predict_rsf")
    
    return(predict_rsf)
  }
  #'  Predict relative probability of selection for den habitat across MWEPA for k training models 
  den_Kpredict <- lapply(trained_den_k_coefs, predict_den_rsf, cov = zcovs_den_mwepa)
  head(den_Kpredict[[1]]); head(den_Kpredict[[7]])
  
  save(den_Kpredict, file = "./Outputs/kfold_predicted_den.RData")
  
  predict_rnd_rsf <- function(coef, cov) {
    predict_rsf <- c()
    #'  Predict across each grid cell
    for(i in 1:nrow(cov)) {
      predict_rsf[i] <- exp(coef$b.elev*cov$Elev[i] + coef$b.rough*cov$Rough[i] + 
                              coef$b.curve*cov$Curve[i] + coef$b.water*cov$Dist2Water[i] + 
                              coef$b.ndvi*cov$SeasonalNDVI[i] + coef$b.ndviXavgndvi*cov$ModifiedNDVI[i])}
    predict_rsf <- as.data.frame(predict_rsf)
    predict_rsf <- cbind(cov$ID, cov$x, cov$y, predict_rsf)
    colnames(predict_rsf) <- c("ID", "x", "y", "predict_rsf")
    
    return(predict_rsf)
  }
  #'  Predict relative probability of selection for rendezvous habitat across MWEPA for k training models 
  rnd_Kpredict <- lapply(trained_rnd_k_coefs, predict_rnd_rsf, cov = zcovs_rnd_mwepa)
  head(rnd_Kpredict[[1]]); head(rnd_Kpredict[[10]])
  
  save(rnd_Kpredict, file = "./Outputs/kfold_predicted_rnd.RData")
  
  #'  ---------------------------------------
  ####  Equal area bins for RSF predictions  ####
  #'  ---------------------------------------
  #' #'  Identify outliers in the predictions
  #' outliers <- function(predicted) { 
  #'   #'  Summarize predicted values
  #'   hist(predicted$predict_rsf, breaks = 100)
  #'   boxplot(predicted$predict_rsf)
  #'   #'  What value represents the 99th percentile in the predicted RSF values
  #'   quant <- quantile(predicted$predict_rsf, c(0.99), na.rm = TRUE)
  #'   #'  Print that value and maximum prediction
  #'   print(quant); print(max(predicted$predict_rsf, na.rm = TRUE))
  #'   #'  Identify the 1% most extreme values and set to 99th percentile value
  #'   predicted <- predicted %>%
  #'     mutate(outlier = ifelse(predict_rsf > quant, "outlier", "not_outlier"),
  #'            adjusted_rsf = ifelse(outlier == "outlier", quant, predict_rsf))
  #'   #'  How many predicted values are above the 99th percentile?
  #'   outlier <- predicted[predicted$outlier == "outlier",]
  #'   outlier <- filter(outlier, !is.na(outlier))
  #'   print(nrow(outlier))
  #'   
  #'   return(predicted)
  #' }
  #' den_Kpredict_outliers <- lapply(den_Kpredict, outliers)
  #' rnd_Kpredict_outliers <- lapply(rnd_Kpredict, outliers)
  
  #'  Load reference raster and identify coordinate system
  ref_raster <- terra::rast("./Shapefiles/WMEPA_masked_grid.tif"); res(ref_raster); crs(ref_raster)
  grid_poly <- st_read("./Shapefiles/WMEPA_masked_polygon.shp") # THIS TAKES AWHILE
  wmepa_grid_pts <- read_csv("./Data/WMEPA_suitable_grid_points.csv")
  nad83 <- st_crs(ref_raster)
  wmepa_grid_pts <- st_as_sf(wmepa_grid_pts, coords = c("X", "Y"), crs = nad83) 
  
  #'  Reclassify RSF predictions into 10 equal area bins (Boyce et al. 2002) and rasterize
  #'  Robust to extreme outlier on either end of distribution
  reclassify_RSF <- function(dat, covs) {
    #'  Create 10 breaks for bins of equal area
    bin_rsf <- quantile(dat$predict_rsf, seq(0, 1, by = 0.1))
    print(bin_rsf)
    
    #'  Cut predicted RSF values into bins based on break points
    dat$equal_area_bins <- cut(dat$predict_rsf, breaks = bin_rsf, include.lowest = TRUE)
    #'  Same thing but relable them to something more useful
    dat$bins <- cut(dat$predict_rsf, breaks = bin_rsf, include.lowest = TRUE,
                    labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
   
    #'  Double check bins are of equal size
    print(table(dat$bins))
    
    #'  Grab cellID from covs and add to dat
    cellID <- dplyr::select(covs, c("cellID", "ID"))
    dat <- full_join(dat, cellID, by = "ID")
    
    #'  Convert to sf object
    dat_poly <- full_join(wmepa_grid_pts, dat, by = "cellID") %>%  
      # full_join(grid_poly, dat, by = c("CellID" = "cellID"))
      # dplyr::select(c("CellID", "bins")) %>% rename("value" = "bins")
      dplyr::select(c("cellID", "bins")) %>% rename("x" = "bins"); crs(dat_poly)
    # dat_poly <- full_join(grid_poly, dat, by = c("CellID" = "cellID"))
    # dplyr::select(c("CellID", "bins")) %>% rename("value" = "bins")
    
    dat_vect <- vect(dat_poly)
    predicted_raster <- rasterize(dat_vect, ref_rast, field = "x")
    
    # dat <- st_as_sf(dat, coords = c("x", "y"), crs = nad83)
    # names(dat) <- c("ID", "predict_rsf", "equal_area_bins", "value", "geometry")
    
    #' #'  Rasterize predictions
    #' predicted_raster <- dat_poly %>% 
    #'   dplyr::select(c(bins, geometry)) %>%
    #'   st_rasterize(st_as_stars(st_bbox(ref_raster), nx = nrow(ref_raster), ny = ncol(ref_raster)))
    #' #template = read_stars("./Shapefiles/WMEPA_masked_grid.tif"), #template = st_as_stars(ref_raster)
    #'                #align = TRUE)
    
    return(dat)
  }
  den_Kpredict_binned <- lapply(den_Kpredict, reclassify_RSF, covs = zcovs_den_mwepa); head(den_Kpredict_binned[[1]])
  rnd_Kpredict_binned <- lapply(rnd_Kpredict, reclassify_RSF, covs = zcovs_rnd_mwepa)
  
  #'  -----------------------
  ####  Cross-validate RSFs  ####
  #'  -----------------------
  #'  Set projection 
  wgs84 <- "+proj=longlat +datum=WGS84 +no_defs"
  
  #'  Retain used locations from each training data set
  testing_pts <- function(test_dat) {
    used_locs <- test_dat %>%
      filter(used == 1) %>%
      dplyr::select(c(Pack_year, used, .folds, geometry)) %>%
      mutate(geometry = gsub("[()]", "", geometry)) %>%
      #'  Split geometry into two columns
      separate_wider_delim(cols = geometry, delim = ',', names = c("x", "y")) %>%
      mutate(x = gsub(".*c", "", x)) %>%
    #'  Convert to sf object
    st_as_sf(coords = c("x", "y"), crs = wgs84)
    return(used_locs)
  }
  testing_pts_den_k <- lapply(data_den_k[[2]], testing_pts); head(testing_pts_den_k[[1]])
  testing_rnd_den_k <- lapply(data_rnd_k[[2]], testing_pts)
  
  
  
  