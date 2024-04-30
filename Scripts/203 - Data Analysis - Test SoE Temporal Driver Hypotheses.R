# compare hypotheses for species with SoE trends
library(here)
library(dplyr)
library(nlme)
library(data.table)
library(AICcmodavg)
library(tseries)
library(MuMIn)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(purrr)

library(ggthemes)
library(hrbrthemes)
library(colorspace)
# get SoE for all species
# read in results files
soe_results <- list.files(path=here("Results/SoE_BBS/Year"), pattern = '*.Rdata', full.names = T) %>%
  purrr::map_df(~ get(load(file = .x)))

# weird duplicates
inds <- soe_results %>%
  group_indices(AOU, year, radius)
soe_results$index <- inds
soe_results <- soe_results[!duplicated(soe_results$index),]
soe_results$index <- NULL

# get radius in KM
soe_results <- soe_results %>%
  mutate(rad_km = radius / 1000)

# identify SoE for all species and calculate weights
soe_results <- soe_results %>% 
  group_by(AOU, year) %>%
  mutate(minAIC = min(AIC)) %>%
  mutate(delAIC = AIC - minAIC) %>%
  mutate(relLik = exp(-0.5 * delAIC)) %>%
  mutate(aicweight = relLik / sum(relLik)) %>%
  ungroup()

# Sample and summarize within groups
results <- soe_results %>%
  group_by(AOU, year) %>%
  summarise(
    y_med = median(replicate(1000, sample(rad_km, replace = TRUE, prob = aicweight))),
    y_sd = sd(replicate(1000, sample(rad_km, replace = TRUE, prob = aicweight)))
  )


# bring in complete species list and focal species list 
species_list <- fread(here('Data/BBS/SpeciesList.csv'))
focal_species <- fread(here('Data/BBS/focal_species.csv'))

# get AOU codes from species list to filter observations
focal_species <- merge(focal_species, species_list[,c('genus_species', "AOU")], 
                       by.x = 'Latin name', by.y = 'genus_species')
rm(species_list)
# Test hypotheses
# 1. Change in species-specific rel. abundance over time
# 2. Increased fragmentation over time
# 3. Change in species richness over time
# 4. Change in overall rel. abundance over time
# 5. Reduced resources (decrease forest cover)
# 6. Change in human development/modification over time

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# 1. Change in species-specific rel. abundance over time ####

# Bring in index values
# load in index data
abun_index <- fread(here('Data/TheNorthAmerica/Index_Best_1966-2019_core.csv')) %>%
  filter(AOU %in% focal_species$AOU)

abun_index <- abun_index %>%
  group_by(AOU, Year) %>%
  summarize(abun_index = mean(Index, na.rm = T))

abun_index <- abun_index %>%
  mutate(Year = as.Date(paste(Year, "01", "01",sep="-"), format="%Y-%m-%d"))
# mutate(Year = format(as.Date(paste(Year, "01", "01",sep="-"), format="%Y-%m-%d"), "%Y" ))

# grab latin, common name, and class from focal_species
abun_index <- merge(abun_index, focal_species, 'AOU')


# subset to time period of interest
abun_index <- abun_index %>% 
  filter(Year > '2000-01-01' & Year < '2019-01-01')

abun_index <- abun_index %>%
  mutate(year = year(Year)) %>%
  select(-Year)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# 2. Change in resources over time (fragmentation) ####

# values were calculated with this script:
# source(here('Scripts/06 - Data Analysis - Fragmentation Hypothesis.R'))

# bring in output from analysis script
frag_list <- list()
counter = 0
for(birds in focal_species$AOU){
  counter = counter + 1
  load(paste0(here('Results/SoE_BBS/Mean_Forest'), '/SOE_BBS_mean_forest_species_AOU_',birds, '.Rdata'))
  holder <- holder %>% 
    filter(Year != 0) %>% 
    group_by(Year) %>%
    summarise(mean_forest_area = mean(mean_forest_area, na.rm = T)) %>%
    mutate(AOU = birds)
  frag_list[[counter]] <- holder
}
fragmentation <- bind_rows(frag_list)
rm(frag_list, holder, counter, birds)

fragmentation <- fragmentation %>%
  rename(year = Year)
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# 3. Change in species richness over time ####
# 4. Change in overall rel. abundance over time (stored in same output) ####

if(!file.exists(here('Results/SoE_BBS/Species_Richness/observed_spp_richness.RDS'))){
  source(here('Scripts/206 - Data Analysis - Spp Rich Hypothesis.R'))
}else{
  spp_rich <- readRDS(here('Results/SoE_BBS/Species_Richness/observed_spp_richness.RDS')) %>%
    filter(AOU %in% focal_species$AOU)
}

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# 5. Reduced habitat (Change in Forest Cover over time) ####

# values were calculated with this script:
# source(here('Scripts/207 - Data Analysis - Forest Hypothesis.R'))

# bring in output from analysis script
prop_list <- list()
counter = 0
for(birds in focal_species$AOU){
  counter = counter + 1
  load(paste0(here('Results/SoE_BBS/Prop_Forest'), '/SOE_BBS_prop_forest_species_AOU_',birds, '.Rdata'))
  holder <- holder %>% 
    filter(Year != 0) %>% 
    group_by(Year) %>%
    summarise(prop_forest = mean(forest, na.rm = T)) %>%
    mutate(AOU = birds)
  prop_list[[counter]] <- holder
}
forest_cover <- bind_rows(prop_list)
rm(prop_list, holder, counter, birds)

forest_cover <- forest_cover %>%
  rename(year = Year)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# 6. Change in Developed Areas over time ####

# values were calculated with this script:
# source(here('Scripts/208 - Data Analysis - Development Hypothesis.R'))

# bring in output from analysis script
# files are larger for some reason so they take longer to read in. A couple of minutes
prop_list <- list()
counter = 0
for(birds in focal_species$AOU){
  counter = counter + 1
  load(paste0(here('Results/SoE_BBS/prop_devel'), '/SOE_BBS_prop_devel_species_AOU_',birds, '.Rdata'))
  holder <- holder %>% 
    filter(Year != 0) %>% 
    group_by(Year) %>%
    summarise(prop_devel = mean(developed, na.rm = T)) %>%
    mutate(AOU = birds)
  prop_list[[counter]] <- holder
}
devel <- bind_rows(prop_list)
rm(prop_list, holder, counter, birds)

devel <- devel %>%
  rename(year = Year)



# FIT MODELS REFLECTING HYPOTHESES ####
hypo_holder <- data.frame(AOU = 0, mod0 = 0, mod1 = 0, mod2 = 0, mod3 = 0, mod4 = 0, mod5 = 0, mod6 = 0)

counter = 0
for(birds in focal_species$AOU){
  counter = counter + 1
  
  
  model_species <- results %>% filter(AOU == birds) %>% select(year, y_med, y_sd, AOU)
  model_species <- merge(model_species, abun_index[, c('AOU', 'abun_index', 'year')])
  model_species <- merge(model_species, fragmentation[, c('AOU', 'mean_forest_area', 'year')])
  model_species <- merge(model_species, spp_rich)
  model_species <- merge(model_species, forest_cover[, c('AOU', 'prop_forest', 'year')])
  model_species <- merge(model_species, devel[, c('AOU', 'prop_devel', 'year')])
  
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
  
  mod3 <- lm(y_med ~ scale(species_observed),
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
  
  
  model_species <- results %>% filter(AOU == birds) %>% select(year, y_med, y_sd, AOU)
  model_species <- merge(model_species, abun_index[, c('AOU', 'abun_index', 'year')])
  model_species <- merge(model_species, fragmentation[, c('AOU', 'mean_forest_area', 'year')])
  model_species <- merge(model_species, spp_rich)
  model_species <- merge(model_species, forest_cover[, c('AOU', 'prop_forest', 'year')])
  model_species <- merge(model_species, devel[, c('AOU', 'prop_devel', 'year')])
  
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
                       mod3 = lm(y_med ~ scale(species_observed), data = model_species, na.action = na.omit, weights = spp_weights),
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
          mod_pred$covariate <- scale(model_species$species_observed)
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
  scale_color_discrete_qualitative(palette = "Dynamic") +
  labs(x = 'Scaled Covariate', y = 'Scale of Effect (km)') +
  geom_line(aes(y=upr), linetype="dotted", linewidth = 1, alpha = 0.35) + 
  geom_line(aes(y=lwr), linetype="dotted", linewidth = 1, alpha = 0.35) +
  theme_clean() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 13)) +
  facet_wrap(~cov_name))
  

ggsave(filename = here('Results/Figures/temporal_hypotheses_results_Dec2023.png'), plot = last_plot(), width = 12, height = 8, device = 'png',
       dpi = 300)


plot_mod <- merge(plot_mod, hypo_holder[, c('AOU', 'best', 'intercept', 'slope')], by = 'AOU')
plot_table <- plot_mod %>%
  select(`Common name`, `Latin name`, Class, cov_name, intercept, slope) %>%
  filter(!duplicated(`Common name`)) %>%
  mutate(intercept = round(intercept, 2),
         slope = round(slope, 2)) %>%
  arrange(cov_name, Class, `Common name`)
write.csv(plot_table, file = here('Results/temporal_hypotheses_results_Dec2023.csv'))

# # summarize information for figure 5 in main manuscript
hypo_holder <- merge(hypo_holder, focal_species, by= 'AOU')
hypo_holder %>%
   group_by(Class) %>%
   count(best)
