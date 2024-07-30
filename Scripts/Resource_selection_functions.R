  #'  ------------------------
  #'  Homesite resource selection functions(RSF)
  #'  Mexican wolf project
  #'  Sarah B. Bassing
  #'  February 2024
  #'  ------------------------
  #'  Script to run competing RSFs to estimate habitat selection for Mexican wolf
  #'  den and rendezvous sites and then identify the most supported models with
  #'  model selection.
  #'  
  #'  Testing competing hypotheses about the influence of resource availability
  #'  (e.g., water), protection from predators (e.g., slope, forest cover), and
  #'  human disturbance (e.g., distance to nearest road) on den and rendezvous
  #'  site seleciton.
  #'  ------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(car)
  library(MuMIn)
  library(tidyverse)
  
  #'  Load used/available location data and covariates
  all_data_den <- read_csv("./Data/all_data_den.csv")
  all_data_rnd <- read_csv("./Data/all_data_rnd.csv")
  
  #'  Standardize covariate data
  ztransform <- function(dat) {
    dat <- as.data.frame(dat) %>%
      dplyr::select(-geometry)
    #'  Center and scale covariates
    dat_z <- dat %>%
      transmute(Pack_year = Pack_year,
                Homesite = Homesite,
                used = as.numeric(used),
                wgts = as.numeric(wgts),
                Elev = as.numeric(scale(Elevation_m)),
                Slope = as.numeric(scale(Slope_degrees)),
                Rough = as.numeric(scale(Roughness_VRM)),
                Curve = as.numeric(scale(Gaussian_curvature)),
                Dist2Water = as.numeric(scale(Nearest_water_m)),
                logDist2Water = as.numeric(scale(log(Nearest_water_m))),
                HumanMod = as.numeric(scale(Human_mod_index)),
                Dist2Road = as.numeric(scale(Nearest_road_m)),
                logDist2Road = as.numeric(scale(log(Nearest_road_m))),
                CanopyCov = as.numeric(scale(Mean_percent_canopy)),
                AvgCanopyCov = as.numeric(scale(avg_MCP_canopycover)),
                SeasonalNDVI = as.numeric(scale(meanNDVI)),
                AvgSeasonalNDVI = as.numeric(scale(avg_MCP_meanNDVI)))  
    
    #'  Assess correlation among all scaled variables
    covs <- dat_z %>%
      dplyr::select(-c(Pack_year, Homesite, used, wgts))
    cor_matrix <- cor(covs, use = "complete.obs")
    print(cor_matrix)
    
    return(dat_z)
  }
  den_dataz <- ztransform(all_data_den)
  rnd_dataz <- ztransform(all_data_rnd)
  
  #'  ------------------------
  ####  Denning habitat RSFs  ####
  #'  ------------------------
  #####  H0: null model  #####
  h0.den <- glm(used ~ 1, data = den_dataz, weights = wgts, family = binomial)
  summary(h0.den)
  
  #####  H1: physical protection  #####
  h1.den <- glm(used ~ Elev + Slope + Rough + CanopyCov + CanopyCov:AvgCanopyCov, data = den_dataz, weight = wgts, family = binomial) 
  summary(h1.den)
  car::vif(h1.den)
  
  h1.den.v2 <- glm(used ~ Elev + I(Elev^2) + Slope + Rough + CanopyCov + CanopyCov:AvgCanopyCov, data = den_dataz, weight = wgts, family = binomial) 
  summary(h1.den.v2)
  
  model.sel(h1.den, h1.den.v2)
  #'  h1.den.v2 is 22.88 deltaAICc better than h1.den
  
  #####  H2: physical protection and water availability  ####
  h2.den <- glm(used ~ Elev + Slope + Rough + Dist2Water + CanopyCov + CanopyCov:AvgCanopyCov, data = den_dataz, weight = wgts, family = binomial) 
  summary(h2.den)
  car::vif(h2.den)
  #'  What about using the log of the distance to nearest waster?
  h2.den.v2 <- glm(used ~ Elev + Slope + Rough + logDist2Water + CanopyCov + CanopyCov:AvgCanopyCov, data = den_dataz, weight = wgts, family = binomial) 
  summary(h2.den.v2)
  
  h2.den.v3 <- glm(used ~ Elev + I(Elev^2) + Slope + Rough + Dist2Water + CanopyCov + CanopyCov:AvgCanopyCov, data = den_dataz, weight = wgts, family = binomial) 
  summary(h2.den.v3)
  
  #'  Which version of the model is most supported
  model.sel(h2.den, h2.den.v2, h2.den.v3)
  #'  h3.den.v3 is 17.41 deltaAICc better than the next best model (h2.den)
  
  #####  H3: human disturbance  #####
  h3.den <- glm(used ~ HumanMod + Dist2Road, data = den_dataz, weight = wgts, family = binomial)
  summary(h3.den)
  car::vif(h3.den)
  #'  What about using the log of the distance to nearest road?
  h3.den.v2 <- glm(used ~ HumanMod + logDist2Road, data = den_dataz, weight = wgts, family = binomial)
  summary(h3.den.v2)
  
  #'  Which version of the model is most supported
  model.sel(h3.den, h3.den.v2)
  #'  DeltaAICc < 2
  #'  Dist2Road pval >0.1, logDist2Road pval <0.1 but >0.05 
  
  #####  H4: global model  #####
  h4.den <- glm(used ~ Elev + Slope + Rough + Dist2Water + CanopyCov + CanopyCov:AvgCanopyCov + HumanMod + Dist2Road, 
                data = den_dataz, weight = wgts, family = binomial) 
  summary(h4.den)
  car::vif(h4.den)
  #'  What about using the log of the distance to nearest water and road?
  h4.den.v2 <- glm(used ~ Elev + Slope + Rough + logDist2Water + CanopyCov + CanopyCov:AvgCanopyCov + HumanMod + logDist2Road, 
                data = den_dataz, weight = wgts, family = binomial)
  summary(h4.den.v2)
  #'  What about mixing and matching based on which variable was most supported for H2 & H3?
  h4.den.v3 <- glm(used ~ Elev + Slope + Rough + Dist2Water + CanopyCov + CanopyCov:AvgCanopyCov + HumanMod + logDist2Road, 
                   data = den_dataz, weight = wgts, family = binomial)
  summary(h4.den.v3)
  
  h4.den.v4 <- glm(used ~ Elev + I(Elev^2) + Slope + Rough + Dist2Water + CanopyCov + CanopyCov:AvgCanopyCov + HumanMod + Dist2Road, 
                data = den_dataz, weight = wgts, family = binomial) 
  summary(h4.den.v4)
  
  model.sel(h4.den, h4.den.v2, h4.den.v3, h4.den.v4)
  #'  h4.den.v4 is 15.52 deltaAICc better than next best (h4.den)
  
  #####  Den RSF model selection using AICc  #####
  (den_ModSelect <- model.sel(h0.den, h1.den.v2, h2.den.v3, h3.den.v2, h4.den.v4))
  #'  h4.den.v4 lowest AICc, next closest model is h2.den.v3 9.59 deltaAICc away
  
  
  
  #'  ---------------------------
  ####  Rendezvous habitat RSFs  ####    
  #'  ---------------------------
  ##### H0: null model  #####
  h0.rnd <- glm(used ~ 1, data = rnd_dataz, weights = wgts, family = binomial)
  summary(h0.rnd)
  
  #####  H1: Ausband model (wet meadows)  #####
  h1.rnd <- glm(used ~ Rough + Curve + SeasonalNDVI + SeasonalNDVI:AvgSeasonalNDVI, data = rnd_dataz, weight = wgts, family = binomial) 
  summary(h1.rnd)
  car::vif(h1.rnd)
  
  #####  H2: water availability  #####
  h2.rnd <- glm(used ~ Elev + Rough + Curve + SeasonalNDVI + SeasonalNDVI:AvgSeasonalNDVI + Dist2Water, data = rnd_dataz, weight = wgts, family = binomial) 
  summary(h2.rnd)
  car::vif(h2.rnd)
  h2.rnd.v2 <- glm(used ~ Elev + Rough + Curve + SeasonalNDVI + SeasonalNDVI:AvgSeasonalNDVI + logDist2Water, data = rnd_dataz, weight = wgts, family = binomial) 
  summary(h2.rnd.v2)
  h2.rnd.v3 <- glm(used ~ Elev + I(Elev^2) + Rough + Curve + SeasonalNDVI + SeasonalNDVI:AvgSeasonalNDVI + Dist2Water, data = rnd_dataz, weight = wgts, family = binomial) 
  summary(h2.rnd.v3)
  model.sel(h2.rnd, h2.rnd.v2, h2.rnd.v3)
  #'  h2.rnd.v3 DeltaAICc = 24.88 better than h2.rnd
    
  #####  H3: human disturbance  #####
  h3.rnd <- glm(used ~ HumanMod + Dist2Road, data = rnd_dataz, weight = wgts, family = binomial)
  summary(h3.rnd)
  car::vif(h3.rnd)
  h3.rnd.v2 <- glm(used ~ HumanMod + logDist2Road, data = rnd_dataz, weight = wgts, family = binomial)
  summary(h3.rnd.v2)

  model.sel(h3.rnd, h3.rnd.v2)
  #'  DeltaAICc < 2
  
  #####  H4: global  #####
  h4.rnd <- glm(used ~ Elev + Rough + Curve + SeasonalNDVI + SeasonalNDVI:AvgSeasonalNDVI + Dist2Water + HumanMod + Dist2Road, 
                data = rnd_dataz, weight = wgts, family = binomial)
  summary(h4.rnd)
  car::vif(h4.rnd)
  h4.rnd.v2 <- glm(used ~ Elev + Rough + Curve + SeasonalNDVI + SeasonalNDVI:AvgSeasonalNDVI + logDist2Water + HumanMod + logDist2Road, 
                data = rnd_dataz, weight = wgts, family = binomial)
  summary(h4.rnd.v2)
  h4.rnd.v3 <- glm(used ~ Elev + I(Elev^2) + Rough + Curve + SeasonalNDVI + SeasonalNDVI:AvgSeasonalNDVI + Dist2Water + HumanMod + Dist2Road, 
                data = rnd_dataz, weight = wgts, family = binomial)
  summary(h4.rnd.v3)
  
  model.sel(h4.rnd, h4.rnd.v2, h4.rnd.v3)
  #'  h4.rnd.v3 deltaAICc = 23.51 better than h4.rnd
  
  #####  Rendezvous site RSF model selection using AICc  #####
  (rnd_ModSelect <- model.sel(h0.rnd, h1.rnd, h2.rnd.v3, h3.rnd, h4.rnd.v3))
  #'  h2.rnd.v3 lowest AICc, h4.rnd.v3 2.5 delta AICc away
  
  
  #'  ------------------
  ####  Results tables  ####
  #'  ------------------
  #####  Coefficient estimates  #####
  #'  Great result table of top model outputs
  top_rsf_out <- function(mod, sitetype){
    #'  Grab coefficient estimate info
    beta <- summary(mod)$coefficients[,1]
    se <- summary(mod)$coefficients[,2]
    pval <- summary(mod)$coefficients[,4]
    #'  Combine into single data frame and organize
    out <- as.data.frame(beta) %>%
      bind_cols(se, pval) %>% 
      #'  Grab parameter names and rename
      mutate(Parameter = row.names(.),
             Parameter = ifelse(Parameter == "(Intercept)", "Intercept", Parameter),
             Parameter = ifelse(Parameter == "Elev", "Elevation", Parameter),
             Parameter = ifelse(Parameter == "I(Elev^2)", "Elevation^2", Parameter),
             Parameter = ifelse(Parameter == "Slope", "Slope", Parameter),
             Parameter = ifelse(Parameter == "Rough", "Roughness", Parameter),
             Parameter = ifelse(Parameter == "Curve", "Curvature", Parameter),
             Parameter = ifelse(Parameter == "Dist2Water", "Distance to water", Parameter),
             Parameter = ifelse(Parameter == "HumanMod", "Human modification", Parameter),
             Parameter = ifelse(Parameter == "Dist2Road", "Distance to road", Parameter),
             Parameter = ifelse(Parameter == "CanopyCov", "Mean canopy cover", Parameter),
             Parameter = ifelse(Parameter == "CanopyCov:AvgCanopyCov", "Modified mean canopy cover", Parameter),
             Parameter = ifelse(Parameter == "SeasonalNDVI", "Mean NDVI", Parameter),
             Parameter = ifelse(Parameter == "SeasonalNDVI:AvgSeasonalNDVI", "Modified mean NDVI", Parameter),
             Site = sitetype) %>%
      #'  Reorganize and rename table
      relocate(Parameter, .before = beta) %>%
      relocate(Site, .before = Parameter) 
    row.names(out) <- NULL
    names(out) <- c("Site type", "Parameter", "Estimate", "SE", "p-value")
    #'  Round estiamtes to more manageable values
    out <- mutate(out, Estimate = round(Estimate, 2),
                  SE = round(SE, 2),
                  `p-value` = round(`p-value`, 3))
    return(out)
  }
  #'  Extract coefficient estimates for each trained model
  topmod_den_coefs <- top_rsf_out(h4.den.v4, sitetype = "Den")          
  topmod_rnd_coefs <- top_rsf_out(h2.rnd.v3, sitetype = "Rendezvous")   
  
  #'  Coefficient results table
  topmod_coefs <- bind_rows(topmod_den_coefs, topmod_rnd_coefs)
  
  #'  Write file for publication
  write_csv(topmod_coefs, file = "./Outputs/Tables/Top_Model_Coefficients.csv")
  
  #####  AIC model ranks  #####
  #'  Generate modelselection objects
  den_ModSelect <- model.sel(h0.den, h1.den.v2, h2.den.v3, h3.den.v2, h4.den.v4)
  rnd_ModSelect <- model.sel(h0.rnd, h1.rnd, h2.rnd.v3, h3.rnd, h4.rnd.v3)
  
  AIC_tbl <- function(aictable, sitetype) {
    #'  Grab relevant information from modelselection object
    modname <- row.names(aictable)
    aicc <- aictable$AICc
    deltaaic <- aictable$delta
    modwgt <- aictable$weight
    #'  Merge and organize into a single table
    modaicselect <- bind_cols(modname, aicc, deltaaic, modwgt)
    names(modaicselect) <- c("Model", "AICc", "deltaAICc", "Model weight")
    modaicselect <- modaicselect %>%
      #'  Add site type to dataframe
      mutate(Site = sitetype,
             #'  Rename models to something more meaningful
             Model = ifelse(Model == "h0.den", "Null", Model),
             Model = ifelse(Model == "h0.rnd", "Null", Model),
             Model = ifelse(Model == "h1.den.v2", "Physical protection", Model),
             Model = ifelse(Model == "h1.rnd", "Wet meadows", Model),
             Model = ifelse(Model == "h2.den.v3", "Physical protection and water", Model),
             Model = ifelse(Model == "h2.rnd.v3", "Water availability", Model),
             Model = ifelse(Model == "h3.den.v2", "Human disturbance", Model),
             Model = ifelse(Model == "h3.rnd", "Human disturbance", Model),
             Model = ifelse(Model == "h4.den.v4", "Global", Model),
             Model = ifelse(Model == "h4.rnd.v3", "Global", Model),
             #'  Round outputs to more manageable values
             AICc = round(AICc, 2),
             deltaAICc = round(deltaAICc, 2),
             `Model weight` = round(`Model weight`, 2)) %>%
      #'  Reorganize
      relocate(Site, .before = Model) %>%
      rename("Site type" = "Site")
    return(modaicselect)
  }
  AIC_table_den <- AIC_tbl(den_ModSelect, sitetype = "Den")
  AIC_table_rnd <- AIC_tbl(rnd_ModSelect, sitetype = "Rendezvous")
  
  #'  Combine into single AIC model ranking table
  AIC_table <- bind_rows(AIC_table_den, AIC_table_rnd)
  
  #'  Write table for publication
  write_csv(AIC_table, file = "./Outputs/Tables/AICc_Model_Rank_Table.csv")
  
  #'  ---------------------
  ####  Visualize results  ####
  #'  ---------------------
  #'  Load MWEPA masked grid and covariate data
  grid_covs <- read_csv("./Data/MWEPA_suitable_grid_covs.csv")
  
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
  
  #'  Standardize grid covs based on mean and SD of covariates in original model
  standardize_mwepa_covs <- function(dat, mu.sd) {
    zcovs <- dat %>%
      transmute(#cellID = cellID,
                ID = ID,
                Elev = (Elevation_m - mu.sd$Elevation_m[1])/mu.sd$Elevation_m[2],
                Elev2 = (Elev ^2), # square the standardize covariate
                Slope = (Slope_degrees - mu.sd$Slope_degrees[1])/mu.sd$Slope_degrees[2],
                Rough = (Roughness_VRM - mu.sd$Roughness_VRM[1])/mu.sd$Roughness_VRM[2],
                Curve = (Gaussian_curvature - mu.sd$Gaussian_curvature[1])/mu.sd$Gaussian_curvature[2],
                Dist2Water = (Nearest_water_m - mu.sd$Nearest_water_m[1])/mu.sd$Nearest_water_m[2],
                HumanMod = (Human_mod_index - mu.sd$Human_mod_index[1])/mu.sd$Human_mod_index[2],
                Dist2Road = (Nearest_road_m - mu.sd$Nearest_road_m[1])/mu.sd$Nearest_road_m[2],
                CanopyCov = (Mean_percent_canopy - mu.sd$Mean_percent_canopy[1])/mu.sd$Mean_percent_canopy[2],
                AvgCanopyCov = (avg_MWEPA_canopycover - mu.sd$avg_MCP_canopycover[1])/mu.sd$avg_MCP_canopycover[2],
                ModifiedCanopyCov = CanopyCov*AvgCanopyCov,
                SeasonalNDVI = (meanNDVI - mu.sd$meanNDVI[1])/mu.sd$meanNDVI[2],
                AvgSeasonalNDVI = (avg_MWEPA_meanNDVI - mu.sd$avg_MCP_meanNDVI[1])/mu.sd$avg_MCP_meanNDVI[2],
                ModifiedNDVI = SeasonalNDVI * AvgSeasonalNDVI,
                x = as.numeric(X),
                y = as.numeric(Y))
               
    return(zcovs)
  }
  zcovs_den_mwepa <- standardize_mwepa_covs(grid_covs, mu.sd = data_den_stats)
  zcovs_rnd_mwepa <- standardize_mwepa_covs(grid_covs, mu.sd = data_rnd_stats)
  
  #'  Function to save parameter estimates from each model
  rsf_out <- function(mod){
    beta <- mod$coefficients
    out <- as.data.frame(beta) %>%
      mutate(Parameter = row.names(.),
             Parameter = ifelse(Parameter == "(Intercept)", "alpha", Parameter),
             Parameter = ifelse(Parameter == "Elev", "b.elev", Parameter),
             Parameter = ifelse(Parameter == "I(Elev^2)", "b.elev2", Parameter),
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
  coefs_h4.den <- rsf_out(h4.den.v4)   
  coefs_h2.rnd <- rsf_out(h2.rnd.v3)   
  
  ####  Predict across suitable habitat ####
  #'  ----------------------------------
  #'  Functions to predict across all grid cells based on RSF results
  #'  NOTE: want the predict relative probability of selection from RSF so not 
  #'  using a logit transformation. Drop intercept from the model and exponentiate 
  #'  coefs*covs (Fieberg et al. 2020)
  predict_den_rsf <- function(coef, cov) {
    #'  Generate modified canopy cover variable
    
    predict_rsf <- c()
    #'  Predict across each grid cell
    for(i in 1:nrow(cov)) {
      predict_rsf[i] <- exp(coef$b.elev*cov$Elev[i] + coef$b.elev2*cov$Elev2[i] + 
                              coef$b.slope*cov$Slope[i] + coef$b.rough*cov$Rough[i] + 
                              coef$b.water*cov$Dist2Water[i] + coef$b.canopy*cov$CanopyCov[i] +
                              coef$b.canopyXavgcanopy*cov$ModifiedCanopyCov[i] +
                              coef$b.hm*cov$HumanMod[i] + coef$b.road*cov$Dist2Road[i])}
    predict_rsf <- as.data.frame(predict_rsf)
    predict_rsf <- cbind(cov$ID, cov$x, cov$y, predict_rsf)
    colnames(predict_rsf) <- c("ID", "x", "y", "predict_rsf")
    
    return(predict_rsf)
  }
  #'  Predict relative probability of selection for den habitat across MWEPA for k training models 
  den_h4.predict <- predict_den_rsf(coefs_h4.den, cov = zcovs_den_mwepa)
  head(den_h4.predict); tail(den_h4.predict)
  
  save(den_h4.predict, file = "./Outputs/den_h4.predict_elev2.RData")
  
  predict_rnd_rsf <- function(coef, cov) {
    predict_rsf <- c()
    #'  Predict across each grid cell
    for(i in 1:nrow(cov)) {
      predict_rsf[i] <- exp(coef$b.elev*cov$Elev[i] + coef$b.elev2*cov$Elev2[i] + 
                              coef$b.rough*cov$Rough[i] + coef$b.curve*cov$Curve[i] + 
                              coef$b.water*cov$Dist2Water[i] + coef$b.ndvi*cov$SeasonalNDVI[i] + 
                              coef$b.ndviXavgndvi*cov$ModifiedNDVI[i])}
    predict_rsf <- as.data.frame(predict_rsf)
    predict_rsf <- cbind(cov$ID, cov$x, cov$y, predict_rsf)
    colnames(predict_rsf) <- c("ID", "x", "y", "predict_rsf")
    
    return(predict_rsf)
  }
  #'  Predict relative probability of selection for rendezvous habitat across MWEPA for k training models 
  rnd_h2.predict <- predict_rnd_rsf(coefs_h2.rnd, cov = zcovs_rnd_mwepa)
  head(rnd_h2.predict); tail(rnd_h2.predict)
  
  save(rnd_h2.predict, file = "./Outputs/rnd_h2.predict_elev2.RData")
  
  ######  Map predicted habitat selection  ######
  #'  -------------------------------------
  library(sf)
  library(terra)
  library(stars)
  library(tidyterra)
  
  #'  Load reference raster and reference polygons
  grid_pts <- st_read("./Shapefiles/WMEPA_grid_clip_pts.shp"); crs(grid_pts)
  xy <- st_coordinates(grid_pts)
   
  #'  Reclassify RSF predictions into 10 equal area bins (Boyce et al. 2002) and rasterize
  #'  Robust to extreme outlier on either end of distribution
  reclassify_RSF <- function(dat, covs, listcount, sitetype) {
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
    
    #' #'  Grab cellID from covs and add to dat
    #' cellID <- dplyr::select(covs, c("cellID", "ID"))
    #' dat <- full_join(dat, cellID, by = "ID") %>%
    #'   rename("CellID" = "cellID") %>%
    #'   mutate(bins = as.numeric(bins))
    #' 
    #' #'  Append data to polygon sf object
    #' predict_poly <- full_join(grid_poly, dat, by = "CellID"); crs(predict_poly)
    #' #'  Rename for st_rasterize
    #' names(predict_poly) <- c("cellID", "newID", "x", "y", "predictions", "equal_area", "value", "geometry")
    
    #'  Append data to polygon sf object
    predict_pts <- dat %>%
      mutate(bins = as.numeric(bins)) %>%
      relocate("ID", .before = "predict_rsf") %>%
      dplyr::select(-c(x, y)) %>%
      cbind(xy) %>%
      relocate(X, .before = "ID") %>%
      relocate(Y, .after = "X")
    #'  Rename for st_rasterize
    names(predict_pts) <- c("x", "y", "ID", "predictions", "equal_area", "value")  
    
    return(predict_pts)
  }
  den_predict_binned <- reclassify_RSF(den_h4.predict, covs = zcovs_den_mwepa)
  rnd_predict_binned <- reclassify_RSF(rnd_h2.predict, covs = zcovs_rnd_mwepa)
  
  #'  Save binned classifications per fold
  save(den_predict_binned, file = "./Outputs/den_predict_binned_elev2.RData")
  save(rnd_predict_binned, file = "./Outputs/rnd_predict_binned_elev2.RData")
  
  #'  Grab the coordinate system
  ref_grid <- terra::rast("./Shapefiles/WMEPA_buffer_grid_clip.tif")
  nad83 <- crs(ref_grid)
  
  #'  Function to rasterize binned data
  rasterize_rsf <- function(predicted_rast) {
    #' #'  Use MWEPA masked grid as the template for rasterizing so the resolution, 
    #' #'  extent, and coordinate system are correct
    #' prediction_raster <- st_rasterize(predicted_rast %>% dplyr::select(value, geometry), 
    #'                                   template = read_stars("./Shapefiles/WMEPA_masked_grid.tif"), 
    #'                                   align = TRUE)
    #' plot(prediction_raster)
    #' 
    #' #'  Convert to a terra raster object
    #' prediction_raster_terra <- rast(prediction_raster)
    #' names(prediction_raster_terra) <- "RSF_bin"
    
    predicted_rast <- dplyr::select(predicted_rast, c(x, y, value))
    prediction_raster <- terra::rast(predicted_rast, type = "xyz", crs = nad83, digits = 6, extent = NULL)
    names(prediction_raster) <- "RSF_bin"
    plot(prediction_raster)
    
    return(prediction_raster)
  }
  den_predict_rast <- rasterize_rsf(den_predict_binned)
  rnd_predict_rast <- rasterize_rsf(rnd_predict_binned)
  
  #'  Crop raster to buffer
  buff <- st_read("./Shapefiles/Homesites/Homesite_buffered_MCP_suitableHabitat.shp") %>%
    st_transform(nad83)
  plot(den_predict_rast)
  plot(buff[1], add = T)
  crop_rast <- crop(den_predict_rast, buff)
  cropped_raster <- mask(crop_rast, buff)
  plot(cropped_raster)
  
  #'  Calculate area of each bin
  #'  Using full raster because this is the scale at which binning occurred 
  #'  (slightly higher proportion of higher bins within the cropped raster based
  #'  on how the full raster was binned and where the buffer overlays the raster)
  bin_area <- expanse(den_predict_rast, unit="km", transform=TRUE, byValue=TRUE, wide=FALSE)
  bin_area$percent = round(100 * bin_area$area / sum(bin_area$area), 3)
  bin_area
  print("Raster area"); sum(bin_area$area)
  
  #'  Save rasterized binned RSFs
  writeRaster(den_predict_rast, filename = "./Shapefiles/Predicted RSFs/den_predict_raster_elev2.tif", overwrite = TRUE)
  writeRaster(rnd_predict_rast, filename = "./Shapefiles/Predicted RSFs/rnd_predict_raster_elev2.tif", overwrite = TRUE)

  #'  Map predicted RSFs
  #'  Define color palette (only bins 7-10 really identifiable)
  mycolors <- c("ivory1", "beige", "lightyellow", "lemonchiffon1", "lemonchiffon2", "khaki1", "gold", "darkgoldenrod1", "darkorange", "firebrick1")
  den_rsf_plot <- ggplot() +
    geom_spatraster(data = den_predict_rast, aes(fill = RSF_bin)) +
    coord_sf(crs = nad83) +
    scale_fill_gradientn(colours = mycolors, na.value = NA, 
                         breaks = seq(0, 10, 2)) +
    labs(fill = "RSF bin") + xlab("Longitude") + ylab("Latitude") +
    ggtitle("Predicted selection bin for Mexican wolf den sites") + 
    theme_bw() + 
    theme(text = element_text(size = 18))
  den_rsf_plot
  
  rnd_rsf_plot <- ggplot() +
    geom_spatraster(data = rnd_predict_rast, aes(fill = RSF_bin)) +
    coord_sf(crs = nad83) +
    scale_fill_gradientn(colours = mycolors, na.value = NA, 
                         breaks = seq(0, 10, 2)) +
    labs(fill = "RSF bin") + xlab("Longitude") + ylab("Latitude") +
    ggtitle("Predicted selection bin for Mexican wolf rendezvous sites") + 
    theme_bw() + 
    theme(text = element_text(size = 18))
  rnd_rsf_plot
  
  ggsave("./Outputs/Figures/RSF_binned_den_plot.tiff", den_rsf_plot, units = "in", 
         height = 6, width = 10, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/Figures/RSF_binned_rnd_plot.tiff", rnd_rsf_plot, units = "in", 
         height = 6, width = 10, dpi = 600, device = 'tiff', compression = 'lzw')
  
   
  #####  Non-linear response to elevation  #####
  #'  -------------------------------------
  #'  Create range of elevation values to predict across and standardize
  min_elev_den <- min(all_data_den$Elevation_m); max_elev_den <- max(all_data_den$Elevation_m)
  min_elev_rnd <- min(all_data_rnd$Elevation_m); max_elev_rnd <- max(all_data_rnd$Elevation_m)
  range_elev_den <- seq(min_elev_den, max_elev_den, length.out = 500)
  range_elev_rnd <- seq(min_elev_rnd, max_elev_rnd, length.out = 500)
  zrange_elev_den <- (range_elev_den - mean(all_data_den$Elevation_m))/sd(all_data_den$Elevation_m)
  zrange_elev_rnd <- (range_elev_rnd - mean(all_data_rnd$Elevation_m))/sd(all_data_rnd$Elevation_m)
  #'  Merge into single data frame and square z-transformed elevation
  elev_den_cov <- cbind(range_elev_den, zrange_elev_den) %>% as.data.frame() %>% mutate(zrange_elev2_den = zrange_elev_den^2)
  elev_rnd_cov <- cbind(range_elev_rnd, zrange_elev_rnd) %>% as.data.frame() %>% mutate(zrange_elev2_rnd = zrange_elev_rnd^2)
  
  #'  DEN RSF PREDICTIONS
  #'  Predict relative prob. of den habitat selection across available elevational gradient
  predict_cov_effect_den <- function(coef, cov) {
    #'  Predict across range of elevation values while holding all other covs at their mean
    #'  (which is 0 since working with standardized covariates)
    predict_rsf <- exp(coef$b.elev*cov$zrange_elev_den + coef$b.elev2*cov$zrange_elev2_den + 
                         coef$b.slope*0 + coef$b.rough*0 + coef$b.water*0 + coef$b.canopy*0 +
                         coef$b.canopyXavgcanopy*0 + coef$b.hm*0 + coef$b.road*0)
    #'  Scale predicted rsf btween 0 and 1
    scaled_rsf <- (predict_rsf - min(predict_rsf))/(max(predict_rsf) - min(predict_rsf))
    #'  Create a data frame with raw elevation, standardized elevation, and predicted RSF
    predict_rsf <- as.data.frame(predict_rsf)
    predict_rsf <- bind_cols(cov$range_elev_den, cov$zrange_elev_den, cov$zrange_elev2_den, predict_rsf, scaled_rsf)
    colnames(predict_rsf) <- c("Elevation", "Elevation_zscale", "Elevation_zscale_sqr", "predicted_rsf", "normalized_rsf")
    return(predict_rsf)
  }
  #'  Predict relative probability of selection for rendezvous habitat across MWEPA for k training models
  predicted_elev_effect_den <- predict_cov_effect_den(coefs_h4.den, cov = elev_den_cov) %>%
    mutate(Site_type = "Den")
  #'  Elevation at highest relative prob. of selection during den season
  max(predicted_elev_effect_den$Elevation)
  max(predicted_elev_effect_den$predicted_rsf)
  (max_selected_elev_den <- with(predicted_elev_effect_den, Elevation[which.max(predicted_rsf)]))

  #'  RENDEZVOUS SITE PREDICTIONS
  #'  Predict relative prob. of rendezvous site habitat selection across available elevational gradient
  predict_cov_effect_rnd <- function(coef, cov) {
    #'  Predict across range of elevation values while holding all other covs at their mean
    #'  (which is 0 since working with standardized covariates)
    predict_rsf <- exp(coef$b.elev*cov$zrange_elev_rnd + coef$b.elev2*cov$zrange_elev2_rnd + 
                         coef$b.rough*0 + coef$b.curve*0 + coef$b.water*0 +
                         coef$b.ndvi*0 + coef$b.ndviXavgndvi*0)
    #'  Scale predicted rsf btween 0 and 1
    scaled_rsf <- (predict_rsf - min(predict_rsf))/(max(predict_rsf) - min(predict_rsf))
    #'  Create a data frame with raw elevation, standardized elevation, and predicted RSF
    predict_rsf <- as.data.frame(predict_rsf)
    predict_rsf <- bind_cols(cov$range_elev_rnd, cov$zrange_elev_rnd, cov$zrange_elev2_rnd, predict_rsf, scaled_rsf)
    colnames(predict_rsf) <- c("Elevation", "Elevation_zscale", "Elevation_zscale_sqr", "predicted_rsf", "normalized_rsf")
    return(predict_rsf)
  }
  #'  Predict relative probability of selection for rendezvous habitat across MWEPA for k training models
  predicted_elev_effect_rnd <- predict_cov_effect_rnd(coefs_h2.rnd, cov = elev_rnd_cov) %>%
    mutate(Site_type = "Rendezvous site")
  #'  Elevation at highest relative prob. of selection during den season
  max(predicted_elev_effect_rnd$Elevation)
  max(predicted_elev_effect_rnd$predicted_rsf)
  (max_selected_elev_rnd <- with(predicted_elev_effect_rnd, Elevation[which.max(predicted_rsf)]))
  
  #'  Merge into one data set
  predicted_elev_effect <- bind_rows(predicted_elev_effect_den, predicted_elev_effect_rnd)
  
  #'  Visualize
  elev_effect_plot <- ggplot(predicted_elev_effect, 
                                          aes(x = Elevation, y = normalized_rsf, color = Site_type)) +
    geom_line() +
    #' #'  Add intervals
    #' geom_ribbon(aes(ymin = lci, ymax = uci, color = Site_type), alpha = 0.3) +
    geom_vline(xintercept = max_selected_elev_den, linetype = 'dashed', col = 'salmon')+
    geom_vline(xintercept = max_selected_elev_rnd, linetype = 'dashed', col = 'darkcyan')+
    theme_bw() +
    theme(text = element_text(size = 16)) +
    xlim(c(1000, max(predicted_elev_effect$Elevation))) + 
    xlab("Elevation (m)") +
    ylab("Relative probability of seleciton") + 
    labs(color = "Pup-rearing site") +
    ggtitle("Predicted effect of elevation habitat selection")
  elev_effect_plot
  
  
  ######  BOOTSTRAP predictions  ######
  #'  ---------------------------
  #'  Bootstrap prediction intervals for predicted effect of elevation on relative probability of selection
  #'  Helpful code: https://cran.r-project.org/web/packages/glm.predict/vignettes/glm.predict.html
  #'  Number of bootstraps
  nboot <- 100
  
    
  #'  Function to bootstrap logistic regression 
  boot <- function(x, model){
    #'  Grab input data from model
    data <- model.frame(model)
    #'  Generate new data set by resampling data using same sample size as original data with replacement
    data_sample <- data[sample(seq_len(nrow(data)), replace = TRUE), ]
    #'  Rename the weights column so it matches original input data
    data_sample$wgts <- data_sample$`(weights)`
    #'  Refit model with new data set
    coef(update(model, data = data_sample))
  }
  
  #'  Den RSF
  rsf.den <- glm(used ~ Elev + I(Elev^2) + Slope + Rough + Dist2Water + CanopyCov + CanopyCov:AvgCanopyCov + HumanMod + Dist2Road, 
                 data = den_dataz, weight = wgts, family = binomial)
  
  #'  Call function nboot times and append bootstrapped coefficients into a single data frame
  betas_boot_den <- do.call(rbind, lapply(1:nboot, boot, rsf.den))
  
  #'  Create new data for each covariate to predict across
  #'  Note: holding all covs but elevation at their mean = 0 & setting intercept to 0 because only interested in linear predictors for RSF
  x.intercept <- x.slope <- x.rough <- x.water <- x.cover <- x.hum <- x.road <- x.avgcover <- rep(0, 100)
  #'  Creating sequence of standardized elevation values based on min and max of standardized data
  x.elev <- seq(min(den_dataz$Elev), max(den_dataz$Elev), length.out = 100)
  #'  Create standardized elevation^2 values based on sequence above
  x.elev2 <- x.elev^2
  #'  Create massive data frame of new data
  newdat <- cbind(x.intercept, x.elev, x.elev2, x.slope, x.rough, x.water, x.cover, x.hum, x.road, x.avgcover) 
  #'  Predict y based on the bootstrapped coefficients and new data
  yhat <- matrix(ncol = nboot, nrow = nrow(newdat))
  predProb <- matrix(ncol = nboot, nrow = nrow(newdat))
  muPredProb <- c()
  lowerp <- c()
  upperp <- c()
  for(i in 1:nboot){
    dat <- as.vector(newdat[i,])
    yhat[,i] <- betas_boot %*% dat
    predProb[,i] <- exp(yhat[,i]) 
    muPredProb[i] <- mean(predProb[,i])
    lowerp[i] <- quantile(predProb[,i], probs = c(0.025))
    upperp[i] <- quantile(predProb[,i], probs = c(0.975))
  }
  predicted_den_rsf_elev <- cbind(x.elev, muPredProb, lowerp, upperp) %>% as.data.frame()
  
  ggplot(predicted_den_rsf_elev, 
         aes(x = x.elev, y = muPredProb)) +
    geom_line() +
    #'  Add intervals
    geom_ribbon(aes(ymin = lowerp, ymax = upperp), alpha = 0.3) 
  
  #'  For a very specific value (highest elevation a max relative prob. of selection)
  newdat2 <- c(0, 2.9789171, 8.87394683, 0, 0, 0, 0, 0, 0, 0)
  yhat2 = betas_boot %*% newdat2
  predProb2 = exp(yhat2) 
  mean(predProb2)
  quantile(predProb2, probs = c(0.025, 0.975))
  
  
  
  
  
  #####  Functional response to available NDVI  ####
  #'  ------------------------------------------
  #'  Create range of meanNDVI values to predict across and standardize
  min_meanNDVI <- min(all_data_rnd$meanNDVI)
  max_meanNDVI <- max(all_data_rnd$meanNDVI)
  range_meanNDVI <- seq(min_meanNDVI, max_meanNDVI, length.out = 500)
  range_meanNDVIz <- (range_meanNDVI - mean(all_data_rnd$meanNDVI))/sd(all_data_rnd$meanNDVI)
  #'  Identify min, mean & max available meanNDVI across MCP during rnd season (unstandardized)
  min_avg_MCP_meanNDVI <- min(all_data_rnd$avg_MCP_meanNDVI)
  mean_avg_MCP_meanNDVI <- mean(all_data_rnd$avg_MCP_meanNDVI)
  max_avg_MCP_meanNDVI <- max(all_data_rnd$avg_MCP_meanNDVI)
  #'  Center and scale
  min_avg_MCP_meanNDVIz <- (min_avg_MCP_meanNDVI - data_rnd_stats$avg_MCP_meanNDVI[1])/data_rnd_stats$avg_MCP_meanNDVI[2]
  mean_avg_MCP_meanNDVIz <- (mean_avg_MCP_meanNDVI - data_rnd_stats$avg_MCP_meanNDVI[1])/data_rnd_stats$avg_MCP_meanNDVI[2]
  max_avg_MCP_meanNDVIz <- (max_avg_MCP_meanNDVI - data_rnd_stats$avg_MCP_meanNDVI[1])/data_rnd_stats$avg_MCP_meanNDVI[2]
  #'  Create modified covariate (meanNDVI*avg_MCP_meanNDIV)
  min_ModifiedNDVIz <- range_meanNDVIz*min_avg_MCP_meanNDVIz
  mean_ModifiedNDVIz <- range_meanNDVIz*mean_avg_MCP_meanNDVIz
  max_ModifiedNDVIz <- range_meanNDVIz*max_avg_MCP_meanNDVIz
  
  #' #'  Create input data set for predictions
  #' func_rsp_df <- function(zcov, modNDVI) {
  #'   dat <- dplyr::select(zcov, c("Elev", "Rough", "Curve", "SeasonalNDVI", "AvgSeasonalNDVI", "Dist2Water"))
  #'   dat <- dat[1:500,]
  #'   dat <- mutate(dat, Elev = 0)
  #'   dat$Elev <- c(0)
  #'   dat$Rough <- c(0)
  #'   dat$Curve <- c(0)
  #'   dat$SeasonalNDVI <- c(range_meanNDVIz)
  #'   dat$AvgSeasonalNDVI <- c(modNDVI)
  #'   dat$Dist2Water <- c(0)
  #'   dat$meanNDVI <- range_meanNDVI
  #'   return(dat)
  #' }
  #' func_rsp_AvailNDIV_low <- func_rsp_df(rnd_dataz, min_avg_MCP_meanNDVIz)
  #' func_rsp_AvailNDIV_mid <- func_rsp_df(rnd_dataz, mean_avg_MCP_meanNDVIz)
  #' func_rsp_AvailNDIV_high <- func_rsp_df(rnd_dataz, max_avg_MCP_meanNDVIz)
  #' 
  #' functional_response_rnd_rsf <- function(mod, cov, availndvi) { 
  #'   #'  Predict across range of NDVI values while holding all other covs at their mean
  #'   #'  (which is 0 since working with standardized covariates) 
  #'   predict_rsf <- predict(mod, newdata = func_rsp_AvailNDIV_low, type = "terms", se.fit = TRUE) 
  #'   #'  Grab fitted values and SEs
  #'   predicted_response <- as.data.frame(predict_rsf[[1]])
  #'   predicted_response <- predicted_response$SeasonalNDVI
  #'   predicted_se <- as.data.frame(predict_rsf[[2]])
  #'   predicted_se <- predicted_se$SeasonalNDVI
  #'   predict_intercept <- attributes(predict(mod, newdata = func_rsp_AvailNDIV_low, type = 'terms'))$constant 
  #'   #'  Create data frame of NDVI effect and prediction intervals
  #'   predicted_effect <- cbind(predicted_response, predicted_se) 
  #'   names(predicted_effect) <- c("predicted_response", "predicted_se")
  #'   predicted_effect <- as.data.frame(predicted_effect) %>% 
  #'     mutate(NDVI = predicted_response,
  #'            lci = NDVI - (predicted_se*1.96),
  #'            uci = NDVI + (predicted_se*1.96),
  #'            Average_AvailNDVI = availndvi,
  #'            meanNDVI = cov$meanNDVI)
  #'   return(predicted_effect)
  #' }
  #' #'  Predict relative probability of selection for rendezvous habitat across MWEPA for k training models 
  #' rnd_h2.predict_AvgNDVI_low <- functional_response_rnd_rsf(h2.rnd, cov = func_rsp_AvailNDIV_low, availndvi = "Low"); head(rnd_h2.predict_AvgNDVI_low); tail(rnd_h2.predict_AvgNDVI_low)
  #' rnd_h2.predict_AvgNDVI_mid <- functional_response_rnd_rsf(h2.rnd, cov = func_rsp_AvailNDIV_mid, availndvi = "Medium")
  #' rnd_h2.predict_AvgNDVI_high <- functional_response_rnd_rsf(h2.rnd, cov = func_rsp_AvailNDIV_high, availndvi = "High")
  #' 
  #' #'  Create one data frame with all predictions
  #' rnd_h2.predict_AvailAvgNDVI <- bind_rows(rnd_h2.predict_AvgNDVI_low, rnd_h2.predict_AvgNDVI_mid, rnd_h2.predict_AvgNDVI_high) %>%
  #'   mutate(Average_AvailNDVI = factor(Average_AvailNDVI, levels = c("Low", "Medium", "High")))
  
  #'  Create input data set for predictions
  func_rsp_df <- function(modNDVI) {
    AvailNDIV <- bind_cols(range_meanNDVI, range_meanNDVIz, modNDVI) %>%
      as.data.frame(.)
    names(AvailNDIV) <- c("meanNDVI", "meanNDVIz", "ModifiedNDVIz")
    return(AvailNDIV)
  }
  func_rsp_AvailNDIV_low <- func_rsp_df(min_ModifiedNDVIz); head(func_rsp_AvailNDIV_low); tail(func_rsp_AvailNDIV_low)
  func_rsp_AvailNDIV_mid <- func_rsp_df(mean_ModifiedNDVIz); head(func_rsp_AvailNDIV_mid); tail(func_rsp_AvailNDIV_mid)
  func_rsp_AvailNDIV_high <- func_rsp_df(max_ModifiedNDVIz); head(func_rsp_AvailNDIV_high); tail(func_rsp_AvailNDIV_high)

  #'  Predict across range of meanNDVI and different levels of average meanNDVI
  #'  across MCP during rnd season
  functional_response_rnd_rsf <- function(coef, cov) {
    #'  Predict across range of NDVI values while holding all other covs at their mean
    #'  (which is 0 since working with standardized covariates)
    predict_rsf <- exp(coef$b.elev*0 + coef$b.rough*0 + coef$b.curve*0 + coef$b.water*0 +
                         coef$b.ndvi*cov$meanNDVIz + coef$b.ndviXavgndvi*cov$ModifiedNDVIz)
    #'  Create a data frame with raw NDVI, standardized NDVI, and predicted RSF
    predict_rsf <- as.data.frame(predict_rsf)
    predict_rsf <- bind_cols(cov$meanNDVI, cov$meanNDVIz, cov$ModifiedNDVIz, predict_rsf)
    colnames(predict_rsf) <- c("meanNDVI", "meanNDVIz", "ModifiedNDVIz", "predict_rsf")
    return(predict_rsf)
  }
  #'  Predict relative probability of selection for rendezvous habitat across MWEPA for k training models
  rnd_h2.predict_AvgNDVI_low <- functional_response_rnd_rsf(coefs_h2.rnd, cov = func_rsp_AvailNDIV_low) %>%
    mutate(Average_AvailNDVI = "Low"); head(rnd_h2.predict_AvgNDVI_low); tail(rnd_h2.predict_AvgNDVI_low)
  rnd_h2.predict_AvgNDVI_mid <- functional_response_rnd_rsf(coefs_h2.rnd, cov = func_rsp_AvailNDIV_mid) %>%
    mutate(Average_AvailNDVI = "Medium"); head(rnd_h2.predict_AvgNDVI_mid); tail(rnd_h2.predict_AvgNDVI_mid)
  rnd_h2.predict_AvgNDVI_high <- functional_response_rnd_rsf(coefs_h2.rnd, cov = func_rsp_AvailNDIV_high) %>%
    mutate(Average_AvailNDVI = "High"); head(rnd_h2.predict_AvgNDVI_high); tail(rnd_h2.predict_AvgNDVI_high)

  #'  Create one data frame with all predictions
  rnd_h2.predict_AvailAvgNDVI <- bind_rows(rnd_h2.predict_AvgNDVI_low, rnd_h2.predict_AvgNDVI_mid, rnd_h2.predict_AvgNDVI_high) %>%
    mutate(Average_AvailNDVI = factor(Average_AvailNDVI, levels = c("Low", "Medium", "High")))
  
  #'  Plot with three slope (min, mean, max) and prediction intervals (NEED TO FIGURE OUT HOW TO CALCULATE VARIANCE!)
  NDVI_functional_response_plot <- ggplot(rnd_h2.predict_AvailAvgNDVI, 
                                          aes(x = meanNDVI, y = predict_rsf, color = Average_AvailNDVI)) +
    geom_line() +
    #' #'  Add intervals
    #' geom_ribbon(aes(ymin = lci, ymax = uci, color = Average_AvailNDVI), alpha = 0.3) +
    geom_hline(yintercept = 1.0, linetype = 'dashed', col = 'gray15')+
    theme_bw() +
    theme(text = element_text(size = 18)) +
    xlab("Mean NDVI within 250m radius of site") +
    ylab("Relative selection strength") + 
    labs(color = "Average available \ngreenleaf biomass") +
    ggtitle("Predicted functional response to NDVI as availability changes")
  NDVI_functional_response_plot
  
  ggsave("./Outputs/Figures/NDVI_functional_response_plot.tiff", NDVI_functional_response_plot, 
         units = "in", height = 8, width = 10, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  
  
  
  