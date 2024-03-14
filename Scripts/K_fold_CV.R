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
  # library(knitr)
  library(sf)
  library(terra)
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
  h4.den <- "used ~ Elev + Slope + Rough + Dist2Water + CanopyCov:AvgCanopyCov + HumanMod + Dist2Road"
  #'  Top rendezvous site model
  h2.rnd <- "used ~ Elev + Rough + Curve + Dist2Water + SeasonalNDVI:AvgSeasonalNDVI"
  
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
      summarise(across(where(is.numeric), sd, na.rm = TRUE))
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
      transmute(ID = ID,
                Elev = (Elevation_m - mu.sd$Elevation_m[1])/mu.sd$Elevation_m[2],
                Slope = (Slope_degrees - mu.sd$Slope_degrees[1])/mu.sd$Slope_degrees[2],
                Rough = (Roughness_VRM - mu.sd$Roughness_VRM[1])/mu.sd$Roughness_VRM[2],
                Curve = (Gaussian_curvature - mu.sd$Gaussian_curvature[1])/mu.sd$Gaussian_curvature[2],
                Dist2Water = (Nearest_water_m - mu.sd$Nearest_water_m[1])/mu.sd$Nearest_water_m[2],
                # logDist2Water = scale(log(as.numeric(Nearest_water_m))),
                HumanMod = (Human_mod_index - mu.sd$Human_mod_index[1])/mu.sd$Human_mod_index[2],
                Dist2Road = (Nearest_road_m - mu.sd$Nearest_road_m[1])/mu.sd$Nearest_road_m[2],
                # logDist2Road = scale(log(as.numeric(Nearest_road_m))),
                # CanopyCov = (Mean_percent_canopy - mu.sd$Mean_percent_canopy[1])/mu.sd$Mean_percent_canopy[2],
                # AvgCanopyCov = (avg_MCP_canopycover - mu.sd$avg_MCP_canopycover[1])/mu.sd$avg_MCP_canopycover[2],
                SeasonalNDVI = (meanNDVI - mu.sd$meanNDVI[1])/mu.sd$meanNDVI[2],
                AvgSeasonalNDVI = (avg_MCP_meanNDVI - mu.sd$avg_MCP_meanNDVI[1])/mu.sd$avg_MCP_meanNDVI[2],
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
             Parameter = ifelse(Parameter == "CanopyCov:AvgCanopyCov", "b.canopyXavgcanopy", Parameter),
             Parameter = ifelse(Parameter == "SeasonalNDVI:AvgSeasonalNDVI", "b.ndviXavgndvi", Parameter)) %>%
             # Parameter = ifelse(Parameter == "CanopyCov", "b.canopy", Parameter),
             # Parameter = ifelse(Parameter == "AvgCanopyCov", "b.avgcanopy", Parameter),
             # Parameter = ifelse(Parameter == "SeasonalNDVI", "b.ndvi", Parameter),
             # Parameter = ifelse(Parameter == "AvgSeasonalNDVI", "b.avgndvi", Parameter)) %>%
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
    predict_rsf <- c()
    #'  Predict across each grid cell
    for(i in 1:nrow(cov)) {
      predict_rsf[i] <- exp(coef$b.elev*cov$Elev[i] + coef$b.slope*cov$Slope[i] + 
                              coef$b.rough*cov$Rough[i] + coef$b.water*cov$Dist2Water[i] + 
                              coef$b.canopyXavgcanopy*cov$CanopyCov[i]*cov$AvgCanopyCov[i] +   ### REALLY NOT SURE IF THIS IS THE RIGHT WAY TO DO IT...
                              # coef$b.canopy*cov$CanopyCov[i]:coef$b.avgcanopy*cov$AvgCanopyCov[i] + 
                              coef$b.hm*cov$HumanMod[i] + coef$b.road*cov$Dist2Road[i])}  
    predict_rsf <- as.data.frame(predict_rsf)
    predict_rsf <- cbind(cov$ID, cov$x, cov$y, predict_rsf)
    colnames(predict_rsf) <- c("ID", "x", "y", "predict_rsf")
    
    return(predict_rsf)
  }
  #'  Predict relative probability of selection for den habitat across MWEPA for k training models 
  den_Kpredict <- lapply(trained_den_k_coefs, predict_den_rsf, cov = zcovs_den_mwepa)
  
  save(den_Kpredict, "./Outputs/kfold_predicted_den.RData")
  
  predict_rnd_rsf <- function(coef, cov) {
    predict_rsf <- c()
    #'  Predict across each grid cell
    for(i in 1:nrow(cov)) {
      predict_rsf[i] <- exp(coef$b.elev*cov$Elev[i] + coef$b.rough*cov$Rough[i] + 
                              coef$b.curve*cov$Curve[i] + coef$b.water*cov$Dist2Water[i] + 
                              coef$b.ndviXavgndvi*cov$SeasonalNDVI[i]*cov$AvgSeasonalNDVI[i])} #coef$b.ndvi*cov$SeasonalNDVI[i]:coef$b.avgndvi*cov$AvgSeasonalNDVI[i]
    predict_rsf <- as.data.frame(predict_rsf)
    predict_rsf <- cbind(cov$ID, cov$x, cov$y, predict_rsf)
    colnames(predict_rsf) <- c("ID", "x", "y", "predict_rsf")
    
    return(predict_rsf)
  }
  #'  Predict relative probability of selection for den habitat across MWEPA for k training models 
  rnd_Kpredict <- lapply(trained_rnd_k_coefs, predict_rnd_rsf, cov = zcovs_rnd_mwepa)
  
  save(rnd_Kpredict, "./Outputs/kfold_predicted_rnd.RData")
  