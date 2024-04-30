# compare hypotheses for species with SoE trends
library(here)
library(dplyr)
library(nlme)
library(data.table)
library(AICcmodavg)
library(MuMIn)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(ggthemes)
library(hrbrthemes)
library(colorspace)
# get SoE for all species
# read in results files
# get SoE for all species
# read in results files
soe_results <- list.files(path=here("Results/SoE_BBS/BCRs/Jan2024"), pattern = '*.Rdata', full.names = T) %>%
  purrr::map_df(~ get(load(file = .x)))

# remove models that have not converged
soe_results <- soe_results %>% 
  filter(converged == 1)

# weird duplicates
inds <- soe_results %>%
  group_indices(AOU, BCR, radius)
soe_results$index <- inds
soe_results <- soe_results[!duplicated(soe_results$index),]
soe_results$index <- NULL
rm(inds)

# get radius in KM
soe_results <- soe_results %>%
  mutate(rad_km = radius / 1000) 

# Create a summary table of group counts
inds <- soe_results %>%
  group_indices(AOU, BCR)
soe_results$group_index <- inds
rm(inds)

group_counts <- soe_results %>% 
  group_by(group_index) %>% 
  summarize(count = n())

# Filter the original data frame based on the count condition
soe_results <- soe_results %>% 
  filter(group_index %in% group_counts$group_index[group_counts$count >= 10])


# identify SoE for all species and calculate weights
soe_results <- soe_results %>% 
  group_by(AOU, BCR) %>%
  mutate(minAIC = min(AIC)) %>%
  mutate(delAIC = AIC - minAIC) %>%
  mutate(relLik = exp(-0.5 * delAIC)) %>%
  mutate(aicweight = relLik / sum(relLik)) %>%
  ungroup()

# Sample and summarize within groups
soe_results <- soe_results %>%
  group_by(AOU, BCR) %>%
  summarise(
    y_med = median(replicate(1000, sample(rad_km, replace = TRUE, prob = aicweight))),
    y_sd = sd(replicate(1000, sample(rad_km, replace = TRUE, prob = aicweight))),
    y_min = min(replicate(1000, sample(rad_km, replace = TRUE, prob = aicweight))),
    y_max = max(replicate(1000, sample(rad_km, replace = TRUE, prob = aicweight)))
  )

# bring in complete species list and focal species list 
focal_species <- fread(here('Data/BBS/focal_species.csv'))

# get AOU codes from species list to filter observations
focal_species <- merge(focal_species, species_list[,c('genus_species', "AOU")], 
                       by.x = 'Latin name', by.y = 'genus_species')
rm(species_list)

# Bring in BCR codes
bcr_codes <- fread(here('Data/TheNorthAmerica/bcr_codes.csv'))


# Test hypotheses
# 1. Change in species-specific rel. abundance over time
# 2. Increased fragmentation over time
# 3. Change in species richness over time
# 4. Change in overall rel. abundance over time
# 5. Reduced resources (decrease forest cover)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# 1. Differences in species-specific rel. abundance over sapce ####

# Bring in index values
# load in index data
abun_index <- fread(here('Data/TheNorthAmerica/Index_Best_1966-2019_core.csv')) %>%
  filter(AOU %in% focal_species$AOU) %>%
  filter(Year == 2018)

# just get BCRS
abun_index <- abun_index %>%
  filter(grepl("^[A-Za-z]+[0-9]+$", Region))

# grab full code names
abun_index <- merge(abun_index, bcr_codes, by.x = 'Region', by.y = 'BCR_code')


# grab latin, common name, and class from focal_species
abun_index <- merge(abun_index, focal_species, 'AOU')


# average index and then keep subset of columns
abun_index <- abun_index %>%
  group_by(AOU, Region) %>%
  mutate(abun_index = mean(Index, na.rm = T)) %>%
  ungroup() %>%
  select(AOU, Region, Year, abun_index, BCR_name, BCR_num, `Common name`, Class)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# 2. Change in resources across space (fragmentation) ####

# bring in BBS observations to get StateNum, route, and BCR
source(here('Scripts/01 - Data Prep - Get Species Observations.R'))
bbs_observations <- bbs_observations %>%
  select(StateNum, Route, BCR)


# values were calculated with this script:
# source(here('Scripts/06 - Data Analysis - patch size over time.R'))

# bring in output from analysis script
frag_list <- list()
counter = 0
for(birds in focal_species$AOU){
  counter = counter + 1
  load(paste0(here('Results/SoE_BBS/Mean_Forest'), '/SOE_BBS_mean_forest_species_AOU_',birds, '.Rdata'))
  load(paste0(here('Results/SoE_BBS/Mean_Forest'), '/SOE_BBS_mean_forest_species_AOU_',birds, '_BCRs.Rdata'))
  
  # keep most recnet year only
  holder <- holder %>% 
    filter(Year ==2018)
  holder_bcr <- holder_bcr %>%
    filter(Year == 2018)
  
  # add in BCR
  holder <- merge(holder, holder_bcr)
  
  # mean_forest_area by BCR
  holder <- holder %>%
    group_by(BCR) %>%
    summarise(mean_forest_area = mean(mean_forest_area, na.rm = T)) %>%
    mutate(AOU = birds)
  frag_list[[counter]] <- holder
}
fragmentation <- bind_rows(frag_list)
rm(frag_list, holder, counter, birds)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# 3. Change in species richness over time ####
# 4. Change in overall rel. abundance over time (stored in same output) ####

if(!file.exists(here('Results/SoE_BBS/Species_Richness/observed_spp_richness_BCRs.RDS'))){
  source(here('Scripts/06 - Data Analysis - Observed Species Richness - BCRs.R'))
}else{
  spp_rich <- readRDS(here('Results/SoE_BBS/Species_Richness/observed_spp_richness_BCRs.RDS')) %>%
    filter(AOU %in% focal_species$AOU)
}

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# 5. Reduced habitat (Change in Forest Cover over time) ####

# values were calculated with this script:
# source(here('Scripts/06 - Data Analysis - prop forest over time.R'))

# bring in output from analysis script
prop_list <- list()
counter = 0
for(birds in focal_species$AOU){
  counter = counter + 1
  load(paste0(here('Results/SoE_BBS/Prop_Forest'), '/SOE_BBS_prop_forest_species_AOU_',birds, '.Rdata'))
  load(paste0(here('Results/SoE_BBS/Prop_Forest'), '/SOE_BBS_prop_forest_species_AOU_',birds, '_BCRs.Rdata'))
  
  # keep most recnet year only
  holder <- holder %>% 
    filter(Year == 2018)
  holder_bcr <- holder_bcr %>%
    filter(Year == 2018) %>%
    select(!forest)
  
  # add in BCR
  holder <- merge(holder, holder_bcr)
  
  # mean_forest_area by BCR
  holder <- holder %>%
    group_by(BCR) %>%
    summarise(prop_forest = mean(forest, na.rm = T)) %>%
    mutate(AOU = birds)
  
  prop_list[[counter]] <- holder
}
forest_cover <- bind_rows(prop_list)
rm(prop_list, holder, counter, birds)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
## 6. Change in Developed Areas over time ####

# values were calculated with this script:
# source(here('Scripts/06 - Data Analysis - prop developed over time.R'))

# bring in output from analysis script
# files are larger for some reason so they take longer to read in. A couple of minutes
prop_list <- list()
counter = 0
for(birds in focal_species$AOU){
  counter = counter + 1
  load(paste0(here('Results/SoE_BBS/prop_devel'), '/SOE_BBS_prop_devel_species_AOU_',birds, '.Rdata'))
  holder <- holder %>% 
    filter(Year != 0) %>% 
    group_by(BCR) %>%
    summarise(prop_devel = mean(developed, na.rm = T)) %>%
    mutate(AOU = birds)
  prop_list[[counter]] <- holder
}
devel <- bind_rows(prop_list)
rm(prop_list, holder, counter, birds)


# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# # -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# # CHECKING RESIDUALS AND CORRELATION
# check_holder <- data.frame(AOU = 0, mod.ols = 0, fixed_var = 0,pow_var = 0, exp_var = 0, constp_var = 0, best = 0)
# counter = 0
# for(birds in focal_species$AOU){
#   counter = counter + 1
#   
#   
#   model_species <- soe_results %>% filter(AOU == birds) %>% select(rad_km, BCR, AOU)
#   model_species <- merge(model_species, abun_index[abun_index$AOU == birds, c( 'abun_index', 'BCR_num')], by.x = 'BCR', by.y = 'BCR_num')
#   model_species <- merge(model_species, fragmentation[fragmentation$AOU == birds, c('mean_forest_area', 'BCR')])
#   model_species <- merge(model_species, spp_rich %>% filter(AOU == birds))
#   model_species <- merge(model_species, forest_cover[forest_cover$AOU == birds, c( 'prop_forest', 'BCR')])
#   
# 
#   
#   mod.ols <- lm(rad_km ~ scale(mean_forest_area), data = model_species)
#   fixedvar <- gls(rad_km ~ scale(mean_forest_area), weights = varFixed(~ BCR),data = model_species)
#   varpow <- gls(rad_km ~ scale(mean_forest_area),weights = varPower(form = ~ BCR), data = model_species)
#   varexp <- gls(rad_km ~ scale(mean_forest_area),weights = varExp(form = ~ BCR), data = model_species) 
#   #varconstp <- gls(rad_km ~ scale(mean_forest_area),weights = varConstPower(form = ~ BCR), data = model_species)
#   
#   ms <- model.sel(mod.ols, fixedvar, varpow, varexp, varconstp)
#   
#   
#   check_holder[counter,'AOU'] <- birds
#   check_holder[counter,'mod.ols'] <- ms['mod.ols', 'delta']
#   check_holder[counter,'fixed_var'] <- ms['fixedvar', 'delta']
#   check_holder[counter,'pow_var'] <- ms['varpow', 'delta']
#   check_holder[counter,'exp_var'] <- ms['varexp', 'delta']
#   check_holder[counter,'constp_var'] <- ms['varconstp', 'delta']
#   check_holder[counter,'best'] <- rownames(ms)[1]
# }
# 
# check_holder
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# FIT MODELS REFLECTING HYPOTHESES ####
hypo_holder <- data.frame(AOU = 0, mod0 = 0, mod1 = 0, mod2 = 0, mod3 = 0, mod4 = 0, mod5 = 0, mod6 = 0)

counter = 0
for(birds in focal_species$AOU){
  
  #if(birds == 7290) next
  counter = counter + 1
  
  
  model_species <- soe_results %>% filter(AOU == birds) %>% select(y_med, y_sd, BCR, AOU)
  model_species <- merge(model_species, abun_index[abun_index$AOU == birds, c( 'abun_index', 'BCR_num')], by.x = 'BCR', by.y = 'BCR_num')
  model_species <- merge(model_species, fragmentation[fragmentation$AOU == birds, c('mean_forest_area', 'BCR')])
  model_species <- merge(model_species, spp_rich %>% filter(AOU == birds))
  model_species <- merge(model_species, forest_cover[forest_cover$AOU == birds, c( 'prop_forest', 'BCR')])
  model_species <- merge(model_species, devel[, c('AOU', 'prop_devel')])
  

  # calculate model weights based on y_sd
  model_species$spp_weights <- 1/model_species$y_sd^2
  
  # null
  
  mod0 <- lm(y_med ~ 1,
             data = model_species, na.action = na.omit, weights = spp_weights)
  
  
  # 1. Change in species-specific rel. abundance over time
  
  mod1 <- lm(y_med ~ scale(abun_index),
             data = model_species, na.action = na.omit,
             weights = spp_weights)
  
  # 2. Change in resources over time (fragmentation)
  
  mod2 <- lm(y_med ~ scale(mean_forest_area),
             data = model_species, na.action = na.omit,
             weights = spp_weights)
  
  # 3. Change in species richness over time
  
  mod3 <- lm(y_med ~ scale(observed_richness),
             data = model_species, na.action = na.omit,
             weights = spp_weights)
  
  # 4. Change in overall rel. abundance over time (stored in same output)
  
  mod4 <- lm(y_med ~ scale(total_indv_counted),
             data = model_species, na.action = na.omit,
             weights = spp_weights)
  
  # 5. Change in forest cover over time
  
  mod5 <- lm(y_med ~ scale(prop_forest),
             data = model_species, na.action = na.omit,
             weights = spp_weights)
  
  # 6. Change in Devleoped Areas over time
  
  mod6 <- lm(y_med ~ scale(prop_devel),
             data = model_species, na.action = na.omit,
             weights = spp_weights)
  
  ms <- model.sel(mod0,mod1, mod2, mod3, mod4, mod5, mod6)
  
  hypo_holder[counter,'AOU'] <- birds
  hypo_holder[counter,'mod0'] <- ms['mod0', 'delta']
  hypo_holder[counter,'mod1'] <- ms['mod1', 'delta']
  hypo_holder[counter,'mod2'] <- ms['mod2', 'delta']
  hypo_holder[counter,'mod3'] <- ms['mod3', 'delta']
  hypo_holder[counter,'mod4'] <- ms['mod4', 'delta']
  hypo_holder[counter,'mod5'] <- ms['mod5', 'delta']
  hypo_holder[counter,'mod6'] <- ms['mod6', 'delta']
  hypo_holder[counter,'best'] <- rownames(ms)[1]
  
}

# get plots
hypo_holder$intercept <- 0
hypo_holder$slope <- 0
counter = 0
mod_holder <- list()
for(birds in hypo_holder$AOU){
  
  counter = counter + 1
  
  
  model_species <- soe_results %>% filter(AOU == birds) %>% select(y_med, y_sd, BCR, AOU)
  model_species <- merge(model_species, abun_index[abun_index$AOU == birds, c( 'abun_index', 'BCR_num')], by.x = 'BCR', by.y = 'BCR_num')
  model_species <- merge(model_species, fragmentation[fragmentation$AOU == birds, c('mean_forest_area', 'BCR')])
  model_species <- merge(model_species, spp_rich %>% filter(AOU == birds))
  model_species <- merge(model_species, forest_cover[forest_cover$AOU == birds, c( 'prop_forest', 'BCR')])
  model_species <- merge(model_species, devel[, c('AOU', 'prop_devel')])
  
  # calculate model weights based on y_sd
  model_species$spp_weights <- 1/model_species$y_sd^2
  best_model <- hypo_holder$best[counter]
  
  if(best_model == 'mod0'){
    next
  } else{
    # Fit the best model
    best_mod <- switch(best_model,
                       mod1 = lm(y_med ~ scale(abun_index), data = model_species, na.action = na.omit, weights = spp_weights),
                       mod2 = lm(y_med ~ scale(mean_forest_area), data = model_species, na.action = na.omit, weights = spp_weights),
                       mod3 = lm(y_med ~ scale(observed_richness), data = model_species, na.action = na.omit, weights = spp_weights),
                       mod4 = lm(y_med ~ scale(total_indv_counted), data = model_species, na.action = na.omit, weights = spp_weights),
                       mod5 = lm(y_med ~ scale(prop_forest), data = model_species, na.action = na.omit, weights = spp_weights),
                       mod6 = lm(y_med ~ scale(prop_devel), data = model_species, na.action = na.omit, weights = spp_weights)
    )
    mod_pred <- as.data.frame(predict(best_mod, interval = 'confidence'))
    mod_pred$model <- best_model
    mod_pred$AOU <- birds
    if(best_model == 'mod1'){
      mod_pred$cov_name <- 'Species-specific Abundance Index'
      mod_pred$covariate <- scale(model_species$abun_index)
    } else{
      if(best_model == 'mod2'){
        mod_pred$cov_name <- 'Forest Patch Size (Fragmentation)'
        mod_pred$covariate <- scale(model_species$mean_forest_area)
      } else{
        if(best_model == 'mod3'){
          mod_pred$cov_name <- 'Species Richness'
          mod_pred$covariate <- scale(model_species$observed_richness)
        } else{
          if(best_model == 'mod4'){
            mod_pred$cov_name <- 'Community Relative Abundance'
            mod_pred$covariate <- scale(model_species$total_indv_counted)
          } else{
            if(best_model == 'mod5'){
              mod_pred$cov_name <- 'Resource Availability (Forest Cover)'
              mod_pred$covariate <- scale(model_species$prop_forest)
            } else{
              mod_pred$cov_name <- 'Human Development'
              mod_pred$covariate <- scale(model_species$prop_devel)
            }
          }
        }
      }
    }
    mod_holder[[counter]] <- mod_pred
    hypo_holder$intercept[counter] <- best_mod$coefficients[1]
    hypo_holder$slope[counter] <- best_mod$coefficients[2]
    
  }
  
  
}

plot_mod <- bind_rows(mod_holder)
plot_mod <- merge(plot_mod, focal_species, by = 'AOU')

(ggg <- ggplot(plot_mod, aes(x=covariate,y=fit, group = factor(AOU), color = Class)) +
    geom_line(linewidth = 1.5, alpha = 0.5) +
    scale_color_discrete_qualitative(palette = "Harmonic") +
    labs(x = 'Scaled Covariate', y = 'Scale of Effect (km)') +
    geom_line(aes(y=upr), linetype="dotted", linewidth = 1, alpha = 0.35) + 
    geom_line(aes(y=lwr), linetype="dotted", linewidth = 1, alpha = 0.35) +
    theme_clean() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 13),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 13)) +
    facet_wrap(~cov_name))


ggsave(filename = here('Results/Figures/spatial_hypotheses_results_Dec2023.png'), plot = last_plot(), width = 12, height = 8, device = 'png',
       dpi = 300)


plot_mod <- merge(plot_mod, hypo_holder[, c('AOU', 'best', 'intercept', 'slope')], by = 'AOU')
plot_table <- plot_mod %>%
  select(`Common name`, `Latin name`, Class, cov_name, intercept, slope) %>%
  filter(!duplicated(`Common name`)) %>%
  mutate(intercept = round(intercept, 2),
         slope = round(slope, 2)) %>%
  arrange(cov_name, Class, `Common name`)
write.csv(plot_table, file = here('Results/spatial_hypotheses_results_Dec2023.csv'))

# summarize information for figure 5 in main manuscript
hypo_holder %>%
  count(best)

hypo_holder <- merge(hypo_holder, focal_species, by= 'AOU')

hypo_holder %>%
  group_by(Class) %>%
  count(best)

good <- hypo_holder %>% filter(best != "mod0")
good <- merge(good, focal_species, by = 'AOU')

not_good <-  hypo_holder %>% filter(best == "mod0")
not_good <- merge(not_good, focal_species, by = 'AOU')

# for(i in 1:nrow(good)){
#   
#   birds <- good[i, 'AOU']
#   
#   model_species <- soe_results %>% filter(AOU == birds) %>% select(rad_km, BCR, AOU)
#   model_species <- merge(model_species, abun_index[abun_index$AOU == birds, c( 'abun_index', 'BCR_num')], by.x = 'BCR', by.y = 'BCR_num')
#   model_species <- merge(model_species, fragmentation[fragmentation$AOU == birds, c('mean_forest_area', 'BCR')])
#   model_species <- merge(model_species, spp_rich %>% filter(AOU == birds))
#   model_species <- merge(model_species, forest_cover[forest_cover$AOU == birds, c( 'prop_forest', 'BCR')])
#   
#   best <- good[i, 'best']
#   
#   if(best == 'mod1'){
#     
#     model_species$abun_index_sc <- as.vector(scale(model_species$abun_index))
#     mod <-   gls(rad_km ~ abun_index_sc, 
#                  weights = varFixed(~ BCR),
#                  data = model_species, na.action = na.omit)
#     
#     MyData <- data.frame(abun_index_sc = seq(min(model_species$abun_index_sc), max(model_species$abun_index_sc),, 100))
#     mod_pred <- predictSE.gls(mod, MyData, se.fit = T)
#     
#     plot(x = model_species$abun_index_sc,
#          y = model_species$rad_km,
#          xlab = "Abundance Index",
#          ylab = "Scale of Effect (km)",
#          pch = 16,
#          type = "n",
#          ylim = c(min(mod_pred[[1]])-1.96*min(mod_pred[[2]]) - 0.5, max(mod_pred[[1]])+1.96*max(mod_pred[[2]]) + 0.5),
#          main = focal_species[AOU == birds, "Common name"])
#     
#     lines(x = MyData$abun_index_sc,y = mod_pred[[1]],lwd = 5)
#     lines(MyData$abun_index_sc, mod_pred[[1]]+1.96*mod_pred[[2]], lty=2, col = "red", lwd =5)
#     lines(MyData$abun_index_sc, mod_pred[[1]]-1.96*mod_pred[[2]], lty=2, col = "red", lwd = 5)
#     
#   } else{
#     if(best == 'mod2'){
#       # 2. Change in resources over time (fragmentation)
#       
#       mod <-   gls(rad_km ~ scale(mean_forest), 
#                    weights = varFixed(~ BCR),
#                    data = model_species, na.action = na.omit)
#       
#     } else{
#       if(best == 'mod3'){
#         model_species$richness_sc <- as.vector(scale(model_species$observed_richness))
#         
#         # 3. Change in species richness over time
#         mod <-   gls(rad_km ~ richness_sc, 
#                      weights = varFixed(~ BCR),
#                      data = model_species, na.action = na.omit)
#         
#         
#         MyData <- data.frame(range_index = seq(min(model_species$rad_km), max(model_species$rad_km),, 100),
#                              richness_sc = seq(min(model_species$richness_sc), max(model_species$richness_sc),, 100))
#         
#         mod_pred <- predictSE.gls(mod, MyData, se.fit = T)
#         
#         plot(x = model_species$richness_sc,
#              y = model_species$rad_km,
#              xlab = "Observed Richness",
#              ylab = "Scale of Effect (km)",
#              pch = 16,
#              type = "n",
#              ylim = c(min(mod_pred[[1]])-1.96*min(mod_pred[[2]]) - 0.5, max(mod_pred[[1]])+1.96*max(mod_pred[[2]]) + 0.5),
#              main = focal_species[AOU == birds, "Common name"])
#         
#         lines(x = MyData$richness_sc,y = mod_pred[[1]],lwd = 5)
#         lines(MyData$richness_sc, mod_pred[[1]]+1.96*mod_pred[[2]], lty=2, col = "red", lwd =5)
#         lines(MyData$richness_sc, mod_pred[[1]]-1.96*mod_pred[[2]], lty=2, col = "red", lwd = 5)
#         
#       } else{
#         if(best == 'mod4'){
#           # 4. Change in overall rel. abundance over time
#           
#           model_species$total_indv_counted_sc <- as.vector(scale(model_species$total_indv_counted))
#           
#           mod <-   gls(rad_km ~ total_indv_counted_sc, 
#                        weights = varFixed(~ BCR),
#                        data = model_species, na.action = na.omit)
#           
#           MyData <- data.frame(range_index = seq(min(model_species$rad_km), max(model_species$rad_km),, 100),
#                                total_indv_counted_sc = seq(min(model_species$total_indv_counted_sc), max(model_species$total_indv_counted_sc),, 100))
#           
#           mod_pred <- predictSE.gls(mod, MyData, se.fit = T)
#           
#           plot(x = model_species$total_indv_counted_sc,
#                y = model_species$rad_km,
#                xlab = "Community Abundance Index",
#                ylab = "Scale of Effect (km)",
#                pch = 16,
#                type = "n",         
#                ylim = c(min(mod_pred[[1]])-1.96*min(mod_pred[[2]]) - 0.5, max(mod_pred[[1]])+1.96*max(mod_pred[[2]]) + 0.5),
#                main = focal_species[AOU == birds, "Common name"])
#           
#           lines(x = MyData$total_indv_counted_sc,y = mod_pred[[1]],lwd = 5)
#           lines(MyData$total_indv_counted_sc, mod_pred[[1]]+1.96*mod_pred[[2]], lty=2, col = "red", lwd =5)
#           lines(MyData$total_indv_counted_sc, mod_pred[[1]]-1.96*mod_pred[[2]], lty=2, col = "red", lwd = 5)
#         } else{
#           # 5. Change in forest cover over time
#           
#           model_species$prop_forest_sc <- as.vector(scale(model_species$prop_forest))
#           
#           
#           mod <-   gls(rad_km ~ prop_forest_sc, 
#                        weights = varFixed(~ BCR),
#                        data = model_species, na.action = na.omit)
#           
#           MyData <- data.frame(range_index = seq(min(model_species$rad_km), max(model_species$rad_km),, 100),
#                                prop_forest_sc = seq(min(model_species$prop_forest_sc), max(model_species$prop_forest_sc),, 100))
#           P0 <- predict(mod, newdata = MyData)
#           
#           plot(x = model_species$prop_forest_sc,
#                y = model_species$rad_km,
#                xlab = "Proportion of Forest Cover",
#                ylab = "SoE",
#                pch = 16,
#                type = "n",
#                ylim = c(min(mod_pred[[1]])-1.96*min(mod_pred[[2]]) - 0.5, max(mod_pred[[1]])+1.96*max(mod_pred[[2]]) + 0.5),
#                main = focal_species[AOU == birds, "Common name"])
# 
#           lines(x = MyData$prop_forest_sc,y = mod_pred[[1]],lwd = 5)
#           lines(MyData$prop_forest_sc, mod_pred[[1]]+1.96*mod_pred[[2]], lty=2, col = "red", lwd =5)
#           lines(MyData$prop_forest_sc, mod_pred[[1]]-1.96*mod_pred[[2]], lty=2, col = "red", lwd = 5)
#         }
#       }
#     }
#   }
#   
# }
# m0 <- gls(rad_km ~ range_index, data = model_species, na.action = na.omit) 
# MyData <- data.frame(range_index = seq(4.92, 6.26, length = 25))
# P0 <- predict(m0, newdata = MyData)
# par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
# plot(x = model_species$range_index,
#      y = model_species$rad_km,
#      xlab = "Range_index",
#      ylab = "SoE",
#      pch = 16,
#      type = "n")
# lines(x = MyData$range_index,
#       y = P0,
#       lwd = 5)





plot_holder <- data.frame()
for(i in 1:nrow(good)){
  
  birds <- good[i, 'AOU']
  
  model_species <- soe_results %>% filter(AOU == birds) %>% select(y_med, BCR, AOU)
  model_species <- merge(model_species, abun_index %>% 
                           filter(AOU == birds) %>% 
                           select(abun_index, BCR_num), by.x = 'BCR', by.y = 'BCR_num')
  model_species <- merge(model_species, fragmentation[, c('AOU', 'mean_forest_area', 'BCR')])
  model_species <- merge(model_species, spp_rich)
  model_species <- merge(model_species, forest_cover[, c('AOU', 'prop_forest', 'BCR')])
  model_species <- merge(model_species, devel[, c('AOU', 'prop_devel', 'BCR')])
  model_species$best <- good[i, 'best']
  plot_holder <- rbind(plot_holder, model_species)
}

plot_holder <- plot_holder %>%
  group_by(AOU) %>%
  mutate(abun_index_sc = as.vector(scale(abun_index)),
         mean_forest_sc = as.vector(scale(mean_forest_area)),
         richness_sc = as.vector(scale(observed_richness)),
         total_indv_counted_sc = as.vector(scale(total_indv_counted)),
         prop_forest_sc = as.vector(scale(prop_forest)),
         prop_devel_sc = as.vector(scale(prop_devel))) %>%
  ungroup()

# get AOU codes from species list to filter observations
plot_holder <- merge(plot_holder, focal_species[,c('Common name', 'Class', "AOU", "Family")], 
                     by = 'AOU')

my_dark = brewer.pal(n = 3, "Dark2")[c(1,3)] #there are 9, I exluded the two lighter hues

#test <- plot_holder %>% filter(best == 'mod1', AOU == 4090)

m1 <- ggplot(data = plot_holder %>% filter(best == 'mod1'), mapping = aes(y = y_med, x = abun_index_sc, color = factor(Class), group = factor(Class))) +
  geom_smooth(method = 'lm')+
  scale_y_continuous(limits = c(0,6)) +
  scale_color_brewer(palette = 'Dark2', name = 'Species Group')+
  labs(x = 'Species-specific Abundance Index', y = 'Scale of Effect (km)')+
  theme_minimal()
m2 <- ggplot(data = plot_holder %>% filter(best == 'mod2'), mapping = aes(y = y_med, x = mean_forest_sc, color = factor(Class), group = factor(Class))) +
  geom_smooth(method = 'lm')+
  scale_color_manual(values = my_dark, name = 'Species Group')+
  labs(x = 'Forest Patch Size (Fragmentation)', y = 'Scale of Effect (km)')+
  scale_y_continuous(limits = c(0,6))+
  theme_minimal()
m3 <- ggplot(data = plot_holder %>% filter(best == 'mod3'), mapping = aes(y = y_med, x = richness_sc, color = factor(Class), group = factor(Class))) +
  geom_smooth(method = 'lm')+
  scale_color_brewer(palette = 'Dark2', name = 'Species Group')+
  labs(x = 'Species Richness', y = 'Scale of Effect (km)')+
  scale_y_continuous(limits = c(0,6))+
  theme_minimal()
m4 <- ggplot(data = plot_holder %>% filter(best == 'mod4'), mapping = aes(y = y_med, x = prop_forest_sc, color = factor(Class), group = factor(Class))) +
  geom_smooth(method = 'lm')+
  scale_color_brewer(palette = 'Dark2', name = 'Species Group')+
  labs(x = 'Community Relative Abundance', y = 'Scale of Effect (km)')+
  scale_y_continuous(limits = c(0,6))+
  theme_minimal()
m5 <- ggplot(data = plot_holder %>% filter(best == 'mod5'), mapping = aes(y = y_med, x = prop_forest_sc, color = factor(Class), group = factor(Class))) +
  geom_smooth(method = 'lm')+
  scale_color_brewer(palette = 'Dark2', name = 'Species Group')+
  labs(x = 'Resource Availability (Forest Cover)', y = 'Scale of Effect (km)')+
  scale_y_continuous(limits = c(0,6))+
  theme_minimal()
m6 <- ggplot(data = plot_holder %>% filter(best == 'mod6'), mapping = aes(y = y_med, x = prop_devel_sc, color = factor(Class), group = factor(Class))) +
  geom_smooth(method = 'lm')+
  scale_color_brewer(palette = 'Dark2', name = 'Species Group')+
  labs(x = 'Human Development', y = 'Scale of Effect (km)')+
  scale_y_continuous(limits = c(0,6))+
  theme_minimal()

(m1 + m2 + m3) / (m4 + m5)
ggsave(filename = here('Results/Figures/spatial_hypotheses_results_dec2023.png'),plot = last_plot(), width = 11, height = 8, device = 'png')


## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# families ####
m1 <- ggplot(data = plot_holder %>% filter(best == 'mod1'), mapping = aes(y = rad_km, x = abun_index_sc, color = factor(Family), group = factor(Family))) +
  geom_smooth(method = 'lm')+
  scale_y_continuous(limits = c(0,6)) +
  scale_color_discrete(name = 'Species Group')+
  labs(x = 'Species-specific Abundance Index', y = 'Scale of Effect (km)')+
  theme_minimal()
m2 <- ggplot(data = plot_holder %>% filter(best == 'mod2'), mapping = aes(y = rad_km, x = mean_forest_sc, color = factor(Class), group = factor(Class))) +
  geom_smooth(method = 'lm')+
  scale_color_discrete(name = 'Species Group')+
  labs(x = 'Forest Patch Size (Fragmentation)', y = 'Scale of Effect (km)')+
  scale_y_continuous(limits = c(0,6))+
  theme_minimal()
m3 <- ggplot(data = plot_holder %>% filter(best == 'mod3'), mapping = aes(y = rad_km, x = richness_sc, color = factor(Class), group = factor(Class))) +
  geom_smooth(method = 'lm')+
  scale_color_discrete(name = 'Species Group')+
  labs(x = 'Species Richness', y = 'Scale of Effect (km)')+
  scale_y_continuous(limits = c(0,6))+
  theme_minimal()
m4 <- ggplot(data = plot_holder %>% filter(best == 'mod4'), mapping = aes(y = rad_km, x = prop_forest_sc, color = factor(Class), group = factor(Class))) +
  geom_smooth(method = 'lm')+
  scale_color_discrete(name = 'Species Group')+
  labs(x = 'Resource Availability (Forest Cover)', y = 'Scale of Effect (km)')+
  scale_y_continuous(limits = c(0,6))+
  theme_minimal()
m5 <- ggplot(data = plot_holder %>% filter(best == 'mod5'), mapping = aes(y = rad_km, x = prop_devel_sc, color = factor(Class), group = factor(Class))) +
  geom_smooth(method = 'lm')+
  scale_color_discrete(name = 'Species Group')+
  labs(x = 'Human Development', y = 'Scale of Effect (km)')+
  scale_y_continuous(limits = c(0,6))+
  theme_minimal()

(m1 + m2 + m3) / (m4 + m5)
ggsave(filename = here('Results/Figures/spatial_hypotheses_results_dec2023.png'),plot = last_plot(), width = 11, height = 8, device = 'png')
