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
                logDist2Road = scale(log(as.numeric(Nearest_road_m))))
    
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
  #####  H1: physical protection  #####
  h1.den <- glm(used ~ Elev + Slope + Rough, data = den_dataz, weight = wgts, family = binomial) #+Forest
  summary(h1.den)
  car::vif(h1.den)
  
  #####  H2: physical protection and water availability  ####
  h2.den <- glm(used ~ Elev + Slope + Rough + Dist2Water, data = den_dataz, weight = wgts, family = binomial) #+Forest
  summary(h2.den)
  car::vif(h2.den)
  #'  What about using the log of the distance to nearest waster?
  h2.den.v2 <- glm(used ~ Elev + Slope + Rough + logDist2Water, data = den_dataz, weight = wgts, family = binomial) #+Forest
  summary(h2.den.v2)
  
  #'  Which version of the model is most supported
  model.sel(h2.den, h2.den.v2)
  #'  DeltaAICc 31.33 (Dist2Water model weight = 1)
  #'  Dist2Water pval <<<0.01, logDist2Water pval <<0.01
  
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
  h4.den <- glm(used ~ Elev + Slope + Rough + Dist2Water + HumanMod + Dist2Road, #+Forest
                data = den_dataz, weight = wgts, family = binomial) 
  summary(h4.den)
  car::vif(h4.den)
  #'  What about using the log of the distance to nearest water and road?
  h4.den.v2 <- glm(used ~ Elev + Slope + Rough + logDist2Water + HumanMod + logDist2Road, #+Forest
                data = den_dataz, weight = wgts, family = binomial)
  summary(h4.den.v2)
  #'  What about mixing and matching based on which variable was most supported for H2 & H3?
  h4.den.v3 <- glm(used ~ Elev + Slope + Rough + Dist2Water + HumanMod + logDist2Road, #+Forest
                   data = den_dataz, weight = wgts, family = binomial)
  summary(h4.den.v3)
  
  model.sel(h4.den, h4.den.v2, h4.den.v3)
  #'  DeltaAICc h4.den and h4.den.v3 <2, h4.den.v2 delta = 41.43 from h4.den
  
  #####  Den RSF model selection using AICc  #####
  (den_ModSelect <- model.sel(h1.den, h2.den, h3.den, h4.den))
  
  #'  ---------------------------
  ####  Rendezvous habitat RSFs  ####    
  #'  ---------------------------
  #####  H1: Ausband model (wet meadows)  #####
  h1.rnd <- glm(used ~ Rough + Curve, data = rnd_dataz, weight = wgts, family = binomial) #+ NDVI
  summary(h1.rnd)
  car::vif(h1.rnd)
  
  #####  H2: water availability  #####
  h2.rnd <- glm(used ~ Elev + Rough + Curve + Dist2Water, data = rnd_dataz, weight = wgts, family = binomial) #+ NDVI
  summary(h2.rnd)
  car::vif(h2.rnd)
  h2.rnd.v2 <- glm(used ~ Rough + Curve + logDist2Water, data = rnd_dataz, weight = wgts, family = binomial) #+ NDVI
  summary(h2.rnd.v2)
  
  model.sel(h2.rnd, h2.rnd.v2)
  #'  DeltaAICc = 6.48 (h2.rnd model weight ~1)
    
  #####  H3: human disturbance  #####
  h3.rnd <- glm(used ~ HumanMod + Dist2Road, data = rnd_dataz, weight = wgts, family = binomial)
  summary(h3.rnd)
  car::vif(h3.rnd)
  h3.rnd.v2 <- glm(used ~ HumanMod + logDist2Road, data = rnd_dataz, weight = wgts, family = binomial)
  summary(h3.rnd.v2)

  model.sel(h3.rnd, h3.rnd.v2)
  #'  DeltaAICc < 2
  
  #####  H4: global  #####
  h4.rnd <- glm(used ~ Elev + Rough + Curve + Dist2Water + HumanMod + Dist2Road, #+ NDVI
                data = rnd_dataz, weight = wgts, family = binomial)
  summary(h4.rnd)
  car::vif(h4.rnd)
  h4.rnd.v2 <- glm(used ~ Elev + Rough + Curve + logDist2Water + HumanMod + logDist2Road, #+ NDVI
                data = rnd_dataz, weight = wgts, family = binomial)
  summary(h4.rnd.v2)
  
  model.sel(h4.rnd, h4.rnd.v2)
  #'  DeltaAICc = 7.17 (h4.rnd model weight ~1)
  
  #####  Rendezvous site RSF model selection using AICc  #####
  (rnd_ModSelect <- model.sel(h1.rnd, h2.rnd, h3.rnd, h4.rnd))
  
  