  #'  ----------------------
  #'  Sensitivity analysis for used:available ratio
  #'  Mexican wolf project
  #'  Sarah B. Bassing
  #'  Februrary 2024
  #'  ----------------------
  #'  Randomly sample available points into different use:available ratios and run
  #'  model to identify the ratio at which slope coefficients asymptope
  #'  ----------------------
  
  #'  Clean work station
  rm(list =ls())
  
  #'  Load packages
  library(purrr)
  library(ggplot2)
  library(tidyverse)
  
  #'  Load data
  data_den <- read_csv("./Data/all_data_den_1to1000ratio.csv")
  data_rnd <- read_csv("./Data/all_data_rnd_1to1000ratio.csv")
  
  #'  Sub-sample data sets to fewer available locations (2, 10, 20, 50, 100, 200, 500) 
  #'  to evaluate the point at which estimated slope coefficients stabilize
  subsample_avial <- function(dat, navail) {
    #'  Pull out available locations
    avail_locs <- filter(dat, used == 0)
    
    #'  Sub-sample number of available locations based on the number of used locations
    #'  per individual/season/year and a desired ratio of used:available locations
    sub_sample <- avail_locs %>%
      group_by(Pack_year) %>%
      slice_sample(n = navail) %>%
      ungroup()
    
    #'  Append sub-sampled available points to used points and re-organize
    full_dat <- dat %>%
      filter(used == 1) %>%
      bind_rows(sub_sample) %>%
      arrange(Pack_year)
    
    #'  Center and scale covariates
    full_dat_z <- full_dat %>%
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
    covs <- full_dat_z %>%
      dplyr::select(-c(Pack_year, Homesite, used, wgts))
    cor_matrix <- cor(covs, use = "complete.obs")
    print(cor_matrix)
    
    return(full_dat_z)
  }
  #'  Sub-sample available locations 100 times with different numbers of samples
  den_1_2_ratio <- rerun(50, subsample_avial(data_den, navail = 2))
  den_1_10_ratio <- rerun(50, subsample_avial(data_den, navail = 10))
  den_1_20_ratio <- rerun(50, subsample_avial(data_den, navail = 20))
  den_1_50_ratio <- rerun(50, subsample_avial(data_den, navail = 50))
  den_1_100_ratio <- rerun(50, subsample_avial(data_den, navail = 100))
  den_1_200_ratio <- rerun(50, subsample_avial(data_den, navail = 200))
  den_1_500_ratio <- rerun(50, subsample_avial(data_den, navail = 500))
  str(den_1_2_ratio[[1]]); head(den_1_2_ratio[[1]])
  
  rnd_1_2_ratio <- rerun(50, subsample_avial(data_rnd, navail = 2))
  rnd_1_10_ratio <- rerun(50, subsample_avial(data_rnd, navail = 10))
  rnd_1_20_ratio <- rerun(50, subsample_avial(data_rnd, navail = 20))
  rnd_1_50_ratio <- rerun(50, subsample_avial(data_rnd, navail = 50))
  rnd_1_100_ratio <- rerun(50, subsample_avial(data_rnd, navail = 100))
  rnd_1_200_ratio <- rerun(50, subsample_avial(data_rnd, navail = 200))
  rnd_1_500_ratio <- rerun(50, subsample_avial(data_rnd, navail = 500))
  head(rnd_1_2_ratio[[1]])
  
  #'  Fit basic den model
  fit_den_rsf <- function(dat){
    mod <- glm(used ~ Elev + Slope + Dist2Water, weight = wgts, data = dat, family = binomial)
    return(mod)
  }
  den_1_2_mod <- lapply(den_1_2_ratio, fit_den_rsf); summary(den_1_2_mod[[1]]); car::vif(den_1_2_mod[[1]])
  den_1_10_mod <- lapply(den_1_10_ratio, fit_den_rsf); summary(den_1_10_mod[[1]]); car::vif(den_1_10_mod[[1]])
  den_1_20_mod <- lapply(den_1_20_ratio, fit_den_rsf); summary(den_1_20_mod[[1]]); car::vif(den_1_20_mod[[1]])
  den_1_50_mod <- lapply(den_1_50_ratio, fit_den_rsf); summary(den_1_50_mod[[1]]); car::vif(den_1_50_mod[[1]])
  den_1_100_mod <- lapply(den_1_100_ratio, fit_den_rsf); summary(den_1_100_mod[[1]]); car::vif(den_1_100_mod[[1]])
  den_1_200_mod <- lapply(den_1_200_ratio, fit_den_rsf); summary(den_1_200_mod[[1]]); car::vif(den_1_200_mod[[1]])
  den_1_500_mod <- lapply(den_1_500_ratio, fit_den_rsf); summary(den_1_500_mod[[1]]); car::vif(den_1_500_mod[[1]])
  
  #'  Fit simplified version of Ausband et al. (2010) RSF  
  fit_rnd_rsf <- function(dat){
    mod <- lm(used ~ Rough + Curve, weight = wgts, data = dat, family = bimomial(link = "logit"))
    return(mod)
  }
  rnd_1_2_mod <- lapply(rnd_1_2_ratio, fit_rnd_rsf); summary(rnd_1_2_mod[[1]]); car::vif(rnd_1_2_mod[[1]])
  rnd_1_10_mod <- lapply(rnd_1_10_ratio, fit_rnd_rsf); summary(rnd_1_10_mod[[1]]); car::vif(rnd_1_10_mod[[1]])
  rnd_1_20_mod <- lapply(rnd_1_20_ratio, fit_rnd_rsf); summary(rnd_1_20_mod[[1]]); car::vif(rnd_1_20_mod[[1]])
  rnd_1_50_mod <- lapply(rnd_1_50_ratio, fit_rnd_rsf); summary(rnd_1_50_mod[[1]]); car::vif(rnd_1_50_mod[[1]])
  rnd_1_100_mod <- lapply(rnd_1_100_ratio, fit_rnd_rsf); summary(rnd_1_100_mod[[1]]); car::vif(rnd_1_100_mod[[1]])
  rnd_1_200_mod <- lapply(rnd_1_200_ratio, fit_rnd_rsf); summary(rnd_1_200_mod[[1]]); car::vif(rnd_1_200_mod[[1]])
  rnd_1_500_mod <- lapply(rnd_1_500_ratio, fit_rnd_rsf); summary(rnd_1_500_mod[[1]]); car::vif(rnd_1_500_mod[[1]])
  
  #'  Save parameters and combine into single data frame per homesite and ratio
  # rounddig <- 2
  rsf_out <- function(mod, navail) {
    est <- as.data.frame(coef(summary(mod)))
    mod_out <- est %>%
      transmute(Parameter = row.names(.),
                Estimate = Estimate) %>%
      mutate(Parameter = ifelse(Parameter == "(Intercept)", "Intercept", Parameter),
             nAvailable = navail)
    rownames(mod_out) <- NULL
    return(mod_out)
  }
  den_1_2_mod_out <- lapply(den_1_2_mod, rsf_out, navail = 2); den_1_2_mod_out <- do.call(rbind.data.frame, den_1_2_mod_out)
  den_1_10_mod_out <- lapply(den_1_10_mod, rsf_out, navail = 10); den_1_10_mod_out <- do.call(rbind.data.frame, den_1_10_mod_out)
  den_1_20_mod_out <- lapply(den_1_20_mod, rsf_out, navail = 20); den_1_20_mod_out <- do.call(rbind.data.frame, den_1_20_mod_out)
  den_1_50_mod_out <- lapply(den_1_50_mod, rsf_out, navail = 50); den_1_50_mod_out <- do.call(rbind.data.frame, den_1_50_mod_out)
  den_1_100_mod_out <- lapply(den_1_100_mod, rsf_out, navail = 100); den_1_100_mod_out <- do.call(rbind.data.frame, den_1_100_mod_out)
  den_1_200_mod_out <- lapply(den_1_200_mod, rsf_out, navail = 200); den_1_200_mod_out <- do.call(rbind.data.frame, den_1_200_mod_out)
  den_1_500_mod_out <- lapply(den_1_500_mod, rsf_out, navail = 500); den_1_500_mod_out <- do.call(rbind.data.frame, den_1_500_mod_out)
  
  rnd_1_2_mod_out <- lapply(rnd_1_2_mod, rsf_out, navail = 2); rnd_1_2_mod_out <- do.call(rbind.data.frame, rnd_1_2_mod_out)
  rnd_1_10_mod_out <- lapply(rnd_1_10_mod, rsf_out, navail = 10); rnd_1_10_mod_out <- do.call(rbind.data.frame, rnd_1_10_mod_out)
  rnd_1_20_mod_out <- lapply(rnd_1_20_mod, rsf_out, navail = 20); rnd_1_20_mod_out <- do.call(rbind.data.frame, rnd_1_20_mod_out)
  rnd_1_50_mod_out <- lapply(rnd_1_50_mod, rsf_out, navail = 50); rnd_1_50_mod_out <- do.call(rbind.data.frame, rnd_1_50_mod_out)
  rnd_1_100_mod_out <- lapply(rnd_1_100_mod, rsf_out, navail = 100); rnd_1_100_mod_out <- do.call(rbind.data.frame, rnd_1_100_mod_out)
  rnd_1_200_mod_out <- lapply(rnd_1_200_mod, rsf_out, navail = 200); rnd_1_200_mod_out <- do.call(rbind.data.frame, rnd_1_200_mod_out)
  rnd_1_500_mod_out <- lapply(rnd_1_500_mod, rsf_out, navail = 500); rnd_1_500_mod_out <- do.call(rbind.data.frame, rnd_1_500_mod_out)
  
  #'  Format data for plotting
  all_coefs <- function(ratio1, ratio2, ratio3, ratio4, ratio5, ratio6, ratio7) {
    coefs <- bind_rows(ratio1, ratio2, ratio3, ratio4, ratio5, ratio6, ratio7) %>%
      mutate(nAvailable = factor(nAvailable, levels = c("2", "10", "20", "50", "100", "200", "500")))
    return(coefs)
  }
  den_coefs <- all_coefs(den_1_2_mod_out, den_1_10_mod_out, den_1_20_mod_out, den_1_50_mod_out, 
                         den_1_100_mod_out, den_1_200_mod_out, den_1_500_mod_out)
  rnd_coefs <- all_coefs(rnd_1_2_mod_out, rnd_1_10_mod_out, rnd_1_20_mod_out, rnd_1_50_mod_out, 
                         rnd_1_100_mod_out, rnd_1_200_mod_out, rnd_1_500_mod_out)
  
  #'  Plot estimated coefficients as number of available locations increases
  boxplot_den_dat <- function(df) {
    #'  Intercept
    plot_intercept <- ggplot(df[df$Parameter == "Intercept",], aes(factor(nAvailable), y = Estimate)) + 
      geom_boxplot() + ggtitle("Intercept") + theme_light() +
      labs(x = "Number of available locations", y = "Average coefficient estimate") 
    plot(plot_intercept)
    #'  Slope parameter #1
    plot_elevation <- ggplot(df[df$Parameter == "Elev",], aes(factor(nAvailable), y = Estimate)) +
      geom_boxplot() + ggtitle("Elevation") + theme_light() +
      labs(x = "Number of available locations", y = "Average coefficient estimate") 
    plot(plot_elevation)
    #'  Slope parameter #2
    plot_slope <- ggplot(df[df$Parameter == "Slope",], aes(factor(nAvailable), y = Estimate)) +
      geom_boxplot() + ggtitle("Slope") + theme_light() +
      labs(x = "Number of available locations", y = "Average coefficient estimate")
    plot(plot_slope)
    #'  Slope parameter #3
    plot_water <- ggplot(df[df$Parameter == "Dist2Water",], aes(factor(nAvailable), y = Estimate)) +
      geom_boxplot() + ggtitle("Distance to nearest waterbody") + theme_light() + 
      labs(x = "Number of available locations", y = "Average coefficient estimate") 
    plot(plot_water)
  }
  den_coef_boxplots <- boxplot_den_dat(den_coefs)
  
  boxplot_rnd_dat <- function(df) {
    #'  Intercept
    plot_intercept <- ggplot(df[df$Parameter == "Intercept",], aes(factor(nAvailable), y = Estimate)) +
      geom_boxplot() + ggtitle("Intercept") + theme_light() +
      labs(x = "Number of available locations", y = "Average coefficient estimate")
    plot(plot_intercept)
    #'  Slope parameter #1
    plot_rough <- ggplot(df[df$Parameter == "Rough",], aes(factor(nAvailable), y = Estimate)) +
      geom_boxplot() + ggtitle("Ruggedness (VRM)") + theme_light() +
      labs(x = "Number of available locations", y = "Average coefficient estimate") 
    plot(plot_rough)
    #'  Slope parameter #3
    plot_curve <- ggplot(df[df$Parameter == "Curve",], aes(factor(nAvailable), y = Estimate)) +
      geom_boxplot() + ggtitle("Gaussian curvature") + theme_light() +
      labs(x = "Number of available locations", y = "Average coefficient estimate") 
    plot(plot_curve)
  }
  rnd_coef_boxplots <- boxplot_rnd_dat(rnd_coefs)
  
  #'  Stand along ggplots
  plot_elev <- ggplot(den_coefs[den_coefs$Parameter == "Elev",], aes(factor(nAvailable), y = Estimate)) +
    geom_boxplot() + ggtitle("Elevation") + theme_light() + theme(axis.title.x=element_blank()) +
    labs(x = "Number of available locations", y = "Average coefficient estimate") 
  plot_slope <- ggplot(den_coefs[den_coefs$Parameter == "Slope",], aes(factor(nAvailable), y = Estimate)) +
    geom_boxplot() + ggtitle("Slope") + theme_light() + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    labs(x = "Number of available locations", y = "Average coefficient estimate")
  plot_rough <- ggplot(rnd_coefs[rnd_coefs$Parameter == "Rough",], aes(factor(nAvailable), y = Estimate)) +
    geom_boxplot() + ggtitle("Ruggedness (VRM)") + theme_light() +
    labs(x = "Number of available locations", y = "Average coefficient estimate") 
  plot_curve <- ggplot(rnd_coefs[rnd_coefs$Parameter == "Curve",], aes(factor(nAvailable), y = Estimate)) +
    geom_boxplot() + ggtitle("Gaussian curvature") + theme_light() + theme(axis.title.y=element_blank()) + 
    labs(x = "Number of available locations", y = "Average coefficient estimate") 
  
  
  library(patchwork) 
  den_sensitivity <- plot_elev + plot_slope + plot_annotation(title = "Effect of sample rate on den RSF coefficents")
  rnd_sensitivity <- plot_rough + plot_curve + plot_annotation(title = "Effect of sample rate on rendezvous site RSF coefficients")
  sensitivity_plots <- den_sensitivity / rnd_sensitivity + plot_annotation(tag_levels = 'a')
  
  #'  Alternative layout
  thm <- theme(plot.title = element_text(face = 2, size = 16))
  top_plot <- wrap_elements((plot_elev + plot_slope) + 
                                   plot_annotation(title = "Effect of sample rate on den RSF coefficents", theme = thm))
  bottom_plot <- wrap_elements(plot_rough + plot_curve + 
                                 plot_annotation(title = "Effect of sample rate on rendezvous site RSF coefficients", theme = thm))
  sensitivity_plots <- top_plot / bottom_plot 
  
  #'  Save for publication
  ggsave("./Outputs/Figures/Sample_rate_sensitivity_plots.tiff", sensitivity_plots, 
         units = "in", height = 7, width = 7, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  