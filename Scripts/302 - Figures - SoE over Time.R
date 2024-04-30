# compare hypotheses for species with SoE trends
library(here)
library(dplyr)
library(nlme)
library(data.table)
library(ggplot2)
library(ggthemes)
library(hrbrthemes)
library(broom)
library(viridis)
library(colorspace)
library(patchwork)

# EDIT: DEC 2023: Dr. Rob Fletcher recommended Hand-wing index in place of wing length
hwi <- read.csv(here('Data/Global-HWI-v1.1/catherinesheard-Global-HWI-1533081/Dataset HWI 2020-04-10.csv'))


# AVONET morphological traits
traits <- read.csv(here('Data/AVONET/ELEData/ELEData/TraitData/AVONET_Raw_Data.csv'))

# some times a given species is reported more than once, so I am using average values here
traits <- traits %>%
  group_by(Species1_BirdLife) %>%
  summarize(mean.wing = mean(Wing.Length, na.rm = T),
            mean.tarsus = mean(Tarsus.Length, na.rm = T),
            mean.kipps = mean(Kipps.Distance, na.rm = T))

# bring in complete species list and focal species list (clark rushing 2019-2020 papers)
species_list <- fread(here('Data/BBS/SpeciesList.csv'))
focal_species <- fread(here('Data/BBS/focal_species.csv'))

# get AOU codes from species list to filter observations
focal_species <- merge(focal_species, species_list[,c('genus_species', "AOU")], 
                       by.x = 'Latin name', by.y = 'genus_species')
rm(species_list)

# add in info from AviBase, Bird et al 2020
focal_species$body.mass <- c(21.6, 2159, 285, 12.6, 14, 14.16, 50.1, 19.44, 18.99, 69.5, #last is M. carolinus
                             71.6, 48.5, 19.9, 14.69, 40.03, 29.13, 9.99, 5.8, 14.3, 9.04,  # last is S. cerulea
                             10.54, 7.64, 9.69, 27.5, 10.2, 26.18, 12.5,18.96, 68.8, # last is T rufum
                             8.74, 18, 11.4)
focal_species$clutch <- c(6, 2, 4.5, 2.5, 4.5, 4.5, 3.5, 4, 3, 5, #last is M. carolinus
                          5.5, 4, 5, 3.5, 3.5, 4, 5.5, 3, 5, 4, # last is S. cerulea
                          4, 4, 4, 4.5, 5, 4, 4.5, 4.5,4, # last is T rufum
                          4.5,4, 4)

focal_species <- merge(focal_species, y = traits[, c('mean.wing', 'Species1_BirdLife')], by.x = 'Latin name', by.y = 'Species1_BirdLife')
focal_species <- merge(focal_species, y = hwi, by.x = 'Latin name', by.y = 'IUCN.name')


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
rm(inds)
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
soe_results <- soe_results %>%
  group_by(AOU, year) %>%
  summarise(
    y_med = median(replicate(1000, sample(rad_km, replace = TRUE, prob = aicweight))),
    y_sd = sd(replicate(1000, sample(rad_km, replace = TRUE, prob = aicweight))),
    y_min = min(replicate(1000, sample(rad_km, replace = TRUE, prob = aicweight))),
    y_max = max(replicate(1000, sample(rad_km, replace = TRUE, prob = aicweight)))
  )


# Grab species information - add it to Soe_results
soe_results <- merge(soe_results, focal_species)

# deal with column and bird names
soe_results <- rename(soe_results, common_name = `Common name`, latin_name = `Latin name`)
soe_results$common_name[soe_results$common_name == "Swainson\x92s warbler"] <- 'Swainsons warbler'


# save summaries for supplementary info

# hold <- soe_results %>% 
#   filter(!duplicated(AOU)) %>% 
#   arrange(y_med) %>% 
#   mutate(y_sd = round(y_sd, 1)) %>%
#   select(common_name, y_med, y_sd, y_min, y_max)
# 
# write.csv(hold, row.names = F, file = here('soe_over_time_dec2023.csv'))


# make factors 
soe_results$year_factor <- factor(soe_results$year - 2001)
soe_results <- soe_results %>%
  mutate(class_fac = factor(Class))

soe_results <- soe_results %>% 
  mutate(hwi_scaled = scale(HWI))

soe_results$year_0 <- soe_results$year - 2001

# estimated scales of effect over time
test <- soe_results %>% filter(AOU == 3260)

out <- list()
counter = 0
for(a in unique(soe_results$AOU)){
  counter = counter + 1
  spp <- soe_results %>% filter(AOU == a)
  fit <- lm(y_med ~ year, data = spp, weights = 1/y_sd^2)
  fitsum <- summary(fit)
  fitsum <- pf(fitsum$fstatistic[1L], fitsum$fstatistic[2L], fitsum$fstatistic[3L], lower.tail = FALSE)
  these <- fit$coefficients[2]
  hold <- augment(x = fit, data = spp, se_fit = T)
  hold$effect <- these
  hold$pvalue <- fitsum
  out[[counter]] <- hold
}

out <- rbindlist(out)

test <- out %>% 
  group_by(common_name) %>%
  dplyr::mutate(
    first = dplyr::first(.fitted),
    last = dplyr::last(.fitted)) %>%
  mutate(direction = ifelse(first > last, 'Decrease', 'Increase')) %>%
  mutate(facet = paste0(direction, '_', Class)) %>%
  mutate(significant = ifelse(pvalue < 0.05, "Significant", "Not_Significant"))

labels <- test %>%
  group_by(common_name) %>%
  mutate(mean_year = mean(year_0),
         mean_rad = mean(y_med, na.rm = T)) %>%
  ungroup() %>%
  filter(!duplicated(common_name))

labels$mean_year[8] <- 2.5 # towhee
labels$mean_year[32] <- 12 #bluebird
labels$mean_year[6] <- 12.5 #orchard
labels$mean_year[23] <- 13 # hooded
labels$mean_year[22] <- 3 #kentucky
labels$mean_year[10] <- 13 #dick
labels$mean_year[28] <- 13 #titmouse
labels$mean_year[12] <- 5.5 #swainsons warbler
labels$mean_year[19] <- 9 #yellow throated warbler
labels$mean_year[13] <- 3 #white eyed vireo
labels$mean_year[18] <- 13 #cewa
labels$mean_year[31] <- 15 #woth
labels$mean_year[21] <- 4 #luwa
labels$mean_year[15] <- 3 #swwa
labels$mean_year[27] <- 13 #bhnu
labels$mean_year[10] <- 2 #dick
labels$mean_year[17] <- 12 #gwwa
labels$mean_year[4] <- 3 #acfl
labels$mean_year[16] <- 4.5 #wewa
labels$mean_year[20] <- 12.5 #prwa

labels$mean_rad[8] <- 2 # towhee
labels$mean_rad[5] <- 1 # fish crow
labels$mean_rad[27] <- 2.4 # nuthatch
labels$mean_rad[28] <- .9 # titmouse
labels$mean_rad[26] <- 1.25 # wren
labels$mean_rad[13] <- 1.75 #white eyed vireo
labels$mean_rad[18] <- 3.5 #cewa
labels$mean_rad[31] <- 1.6 #woth
labels$mean_rad[21] <- 1.0 #luwa
labels$mean_rad[22] <- 2.75 #kewa
labels$mean_rad[15] <- 2.25 #kewa
labels$mean_rad[32] <- 2.45 #eabl
labels$mean_rad[1] <- 0.45 #blvu
labels$mean_rad[25] <- 3.75 #brth
labels$mean_rad[27] <- 3 #bhnu
labels$mean_rad[24] <- 3.75 #nomo
labels$mean_rad[3] <- 3.5 #rbwo
labels$mean_rad[7] <- 4.25 #fisp
labels$mean_rad[10] <- 5.1 #dick
labels$mean_rad[4] <- 3.75 #acfl
labels$mean_rad[29] <- 4.61 #acfl
labels$mean_rad[15] <- 2.5 #swwa
labels$mean_rad[12] <- 2 #swwa
labels$mean_rad[23] <- 1 #hooded
labels$mean_rad[16] <- .45 #hooded
labels$mean_rad[19] <- .65 #ytwa
labels$mean_rad[6] <- 1.2 #ytwa
labels$mean_rad[20] <- 3 #prwa
#labels %>% filter(common_name == 'Northern mockingbird') %>% pull(mean_rad)
#which(labels$common_name == 'Brown thrasher')

(ggj <- ggplot(test, aes(x=year_0,y=.fitted, color = common_name)) +
  geom_line(linewidth = 1.5, alpha = 0.5) +
  scale_y_continuous(limits = c(-0.5,6)) +
  scale_color_discrete_sequential(palette = "my_pal") +
  labs(x = 'Year', y = 'Scale of Effect (km)') +
  geom_line(aes(y=.fitted + 1.96*.se.fit), linetype="dotted", linewidth = 1, alpha = 0.35) + 
  geom_line(aes(y=.fitted - 1.96*.se.fit), linetype="dotted", linewidth = 1, alpha = 0.35) +
  geom_text(data = labels, aes(x = mean_year, y = mean_rad + 0.05, label = common_name, vjust = -0.25, angle = effect*100), size = 5) + 
  theme_clean() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 13),
        legend.position="none") +
  facet_wrap(~direction + Class))

ggsave(filename = here('Results/Figures/scale_over_time_v02_dec2023.png'), plot = last_plot(), width = 12, height = 8, device = 'png',
       dpi = 300)



## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

(ggplot(data = soe_results, mapping = aes(x = year, y = y_med, color = Class)) +
    scale_color_brewer(palette = "Dark2") +
    stat_smooth(stat="smooth",method = "lm",
              size = 1.5,
              se = F,
              linetype ="solid",
              alpha = 0.55))



soe_results <- soe_results %>% 
  group_by(AOU) %>%
  arrange(year) %>%
  dplyr::mutate(
    first = dplyr::first(y_med),
    last = dplyr::last(y_med),
    mean_rad = mean(y_med, na.rm = T),
    min_rad = min(y_med, na.rm = T),
    max_rad = max(y_med, na.rm = T),
    low_rad = mean(y_med - 1.96*y_sd,na.rm = T),
    up_rad = mean(y_med + 1.96*y_sd,na.rm = T),
    common_name = gsub(" ", "_", common_name)) %>%
  mutate(direction = ifelse(first > last, "decreasing", 'increasing')) %>%
  ungroup()


soe_results <- soe_results %>% 
  mutate(diff = abs(first - last)) %>% 
  mutate(diff = ifelse(diff >= first*0.25, 1, 0)) %>%
  mutate(AOU = factor(AOU))


hold <- ggplot(data = soe_results, mapping = aes(x = year, y = y_med, group = factor(AOU), color = factor(diff))) +
  scale_color_brewer(palette = "Dark2") +
  geom_line(stat="smooth",method = "lm",
            size = 1.5,
            linetype ="solid",
            alpha = 0.55) +
  facet_wrap(~Class, scales = 'fixed') 

(hold <- ggplot(data = soe_results, mapping = aes(x = year, y = y_med, color = common_name)) +
 # scale_color_brewer(palette = "Dark2") +
  geom_line(stat="smooth",method = "lm",
            size = 1.5,
            linetype ="solid",
            alpha = 0.55) +
  facet_wrap(~Class, scales = 'fixed') )
# reorder(common_name,mean_rad)



ggplot() +
  scale_color_brewer(palette = "Dark2") +
  geom_point(data = soe_results, mapping = aes(x = reorder(common_name, mean_rad), y = y_med, group = factor(AOU), color = Class)) +
  geom_point(data = soe_results, mapping = aes(x =reorder(common_name,mean_rad), y = mean_rad, group = factor(AOU)), size = 3) +
  coord_flip() +
  facet_wrap(~Class, scales = 'fixed')

test <- soe_results %>%
  arrange(Class, mean_rad) %>%
  mutate(index = 1:n())

test <- test %>%
  mutate(common_name = gsub("_", " ", common_name))

(v01 <- ggplot() +
  scale_color_brewer(palette = "Dark2") +
  geom_linerange(data = test, mapping = aes(ymin = low_rad, ymax = up_rad, x =reorder(common_name, index), y = mean_rad, color = Class), linewidth = 0.75) +#
 # geom_point(data = test, mapping = aes(x = reorder(common_name, index), y = y_med,  color = Class)) +
  geom_point(data = test, mapping = aes(x =reorder(common_name, index), y = mean_rad, color = Class), size = 2) +
  coord_flip()+
  labs(x = 'Species', y = 'Scale of Effect (km)') +
  theme(axis.text.x = element_text(hjust = 2)) +
  theme_clean() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 13),
        legend.position="none") +
  facet_wrap(~Class, scales = 'fixed'))

ggsave(filename = here('Results/Figures/scale_over_time_dec2023.png'), plot = last_plot(), device = 'png',width = 8, height = 4, units = 'in',dpi = 300)

v01 / ggj
ggsave(filename = here('Results/Figures/scale_over_time_combo_dec2023.png'), plot = last_plot(), device = 'png',width = 16, height = 16, units = 'in',dpi = 300)


save_this <- test %>% group_by(common_name) %>% summarize(mean_soe = mean(rad_km), sd_soe = sd(rad_km), min_soe = min(rad_km), max_soe = max(rad_km)) %>% arrange(mean_soe)
write.csv(save_this, file = here('Results/soe_over_time.csv'))
