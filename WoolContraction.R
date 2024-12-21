## Wool Contraction
## Last updated: 12/20/2024


library(readxl)
library(tidyverse)
library(prioritizr) #https://prioritizr.net/articles/prioritizr.html
library(rgdal)
library(rcbc)
library(terra)
library(highs)
library(rnaturalearth)
library(extrafont)
library(showtext)
library(reshape2)

loadfonts(device = "win")

wd <- #Set working directory to where repo is stored

setwd(wd)

#Suffix for Wool Contraction: wc

#Need:
# 1. Land Footprint of Existing Wool (goats, GLWv4)
# 2. Target Layers: biodiversity (species richness), carbon potential, NCP (local)

#1.1 Existing Wool
existing_sheep_wc <- rast('Wool/Sheep_DA_heads_ek4_snap.tif')

#1.2 Land template
land_template_wc <- rast('template_raster_10k_ek4.tif')

#2.1 Carbon
animal_potential_carbon_wc <- rast('Animal_COC_ek4_10k_area_snap_B1.tif')
names(animal_potential_carbon_wc) <- "Carbon"

#2.2 Biodiversity
biodiversity_wc <- rast('vert_rsr_richness.tif') %>%
  project(crs(land_template_wc)) # CB
names(biodiversity_wc) <- "Biodiversity"

#2.3 NCP
ncp_wc <- rast('local_NCP_all_targets_ek4_10k_snap.tif')
names(ncp_wc) <- "NCP"

#####################
## Additional Prep For Prioritizr
#####################

#Cropping template to least common extent
template_wc <- crop(biodiversity_wc, animal_potential_carbon_wc)
template_wc <- crop(template_wc, biodiversity_wc)
template_wc <- crop(template_wc, existing_sheep_wc)
template_wc <- crop(template_wc, ncp_wc)

print(template_wc)

#Cropping each layer to template
existing_sheep_wc <- crop(existing_sheep_wc, template_wc)
biodiversity_wc <- crop(biodiversity_wc, template_wc)
animal_potential_carbon_wc <- crop(animal_potential_carbon_wc, template_wc)
land_wc <- crop(land_template_wc, template_wc)
ncp_wc <- crop(ncp_wc, template_wc)

#Scaling Because Prioritizr Throws Errors
existing_sheep_wc <- existing_sheep_wc * 0.01

animal_potential_carbon_wc <- animal_potential_carbon_wc * 0.000000001
animal_potential_carbon_wc[animal_potential_carbon_wc < 0.001] <- NA

biodiversity_wc[biodiversity_wc < 0.001] <- NA

#Masking Biodiversity & Carbon To Footprint
animal_potential_carbon_wc <- mask(animal_potential_carbon_wc, existing_sheep_wc)
biodiversity_wc <- mask(biodiversity_wc, existing_sheep_wc)
land_wc <- mask(land_wc, existing_sheep_wc)
ncp_wc <- mask(ncp_wc, existing_sheep_wc)

#Convert to "raster" because Prioritizr handles it better
existing_sheep_wc_r <- raster(existing_sheep_wc)
biodiversity_wc_r <- raster(biodiversity_wc)
animal_potential_carbon_wc_r <- raster(animal_potential_carbon_wc)
land_wc_r <- raster(land_wc)
ncp_wc_r <- raster(ncp_wc)

#Stack Target Features
feats_wc_r <- stack(biodiversity_wc_r, animal_potential_carbon_wc_r, ncp_wc_r)

#Plot Check
plot(feats_wc_r)

###########################################
###########################################
###########################################

#####################
## LAND BUDGET - PRIORITIZING NATURE
#####################

# 1 pct of land footprint
land_budget_wc <- 0.01 * cellStats(land_wc_r, "sum", na.rm = TRUE)

# Set up prioritizr problem
land_p1_wc <- problem(x = land_wc_r, features = feats_wc_r) %>%
  add_min_shortfall_objective(land_budget_wc) %>% 
  add_relative_targets(c(1,1,1)) %>% 
  add_gurobi_solver(gap = 0.01) ## AQUIRE GUROBI SOLVER -- IF IN ACADEMIA, THIS IS FREE
  ## If no access to gurobi solver, use CBC solver, although this solver is significantly
  ## if not prohibitively slower https://prioritizr.net/reference/index.html#solvers 

# Solve prioritizr problem
land_s1_wc <- solve(land_p1_wc)

# Plot solution
plot(land_s1_wc)

# Confirm pixel counts
freq(land_s1_wc)

# Create summary eval table
eval_p05_land_wc <- eval_target_coverage_summary(land_p1_wc,land_s1_wc)

writeRaster(
  land_s1_wc, "Wool/land_s1_wc_1pct.tif")


#####################
## LAND BUDGET LOOP
#####################

# Add the percent column to the eval table
eval_p05_land_wc$budgetpct <- c(1)

# Step wise, from 5% to 100% by 5
for (i in seq(from = 5, to = 100, by = 5)) {
  land_budget_wc <- (i * 0.01) * cellStats(land_wc_r, "sum", na.rm = TRUE)
  
  land_p1_wc_i <- problem(x = land_wc_r, features = feats_wc_r) %>%
    add_min_shortfall_objective(land_budget_wc) %>% 
    add_relative_targets(c(1,1,1)) %>% 
    add_gurobi_solver(gap = 0.05)
  
  print(paste0("Just checked the problem set for ", i))
  
  land_s1_wc_i <- solve(land_p1_wc_i)
  
  eval_p05_land_wc_i <- eval_target_coverage_summary(land_p1_wc_i,land_s1_wc_i)
  
  eval_p05_land_wc_i$budgetpct <- c(i)
  
  eval_p05_land_wc <- rbind(eval_p05_land_wc,eval_p05_land_wc_i)
  
  filename <- paste0("Wool/land_s1_wc", i, "pct.tif", overwrite = TRUE)
  
   writeRaster(
     land_s1_wc_i, filename)
  
}

# Plot step wise scenario results

ggplot(eval_p05_land_wc, aes(x = budgetpct, y = relative_held, color = feature)) +
  #geom_point() +
  geom_abline(intercept = 0, slope = 0.01, linetype = 2) +
  scale_colour_manual(values = c("#16c6dd", "#8ec92f", "#9391ff"), name = "Feature",
                      labels = c("Biodiversity", "Carbon", "NCP")) +
  theme_minimal() + 
  theme(text = element_text(family = "Corbel", size = 16)) +
  labs(x = "Land Budget Percentage", title = "Wool Scenario 3 (Contraction + Nature)" , y = "Ratio of Feature Retained", ) +
  geom_line(size = 1)


ggsave("Wool/feature_curves_WC.png",
       width = 6, height = 3, units = "in",
       bg = "white")

# Save table for any future plots
write.csv(eval_p05_land_wc, 'Wool/eval_p05_land_wc.csv')


#####################
## LAND BUDGET -- NO BIO / NO CARBON
#####################

existing_sheep_wc_r_inv <- 1/existing_sheep_wc_r

# Run once if Prioritizr is throwing scale errors
existing_cotton_cc_r_inv <- existing_cotton_cc_r_inv * 0.1

#1 pct of land
land_budget_wc <- 0.01 * cellStats(land_wc_r, "sum", na.rm = TRUE)

land_p1_wc_nb <- problem(x = land_wc_r, features = existing_sheep_wc_r_inv) %>%
  add_min_shortfall_objective(land_budget_wc) %>% 
  add_relative_targets(c(1)) %>% 
  add_gurobi_solver(gap = 0.01)

land_s1_wc_nb <- solve(land_p1_wc_nb)

plot(land_s1_wc_nb)

freq(land_s1_wc_nb)

writeRaster(
  land_s1_wc_nb, "Wool/land_s1_wc_1pct_noBioNoC.tif", overwrite = TRUE)

eval_p05_land_wc_nb <- eval_target_coverage_summary(land_p1_wc_nb,land_s1_wc_nb)

#####################
## LAND BUDGET LOOP _ NO BIO / NO CARBON
#####################

# Add the percent column to the eval table
eval_p05_land_wc_nb$budgetpct <- c(1)

land_s1_wc_nb[land_s1_wc_nb < 1] <- NA

# Add variable ratios to eval table 


#Biodiversity
biodiversity_wc_r_i <- mask(biodiversity_wc_r, land_s1_wc_nb)
bio_ratio_i <- cellStats(biodiversity_wc_r_i, "sum", na.rm = TRUE) / cellStats(biodiversity_wc_r, "sum", na.rm = TRUE)

#Carbon
animal_potential_carbon_wc_r_i <- mask(animal_potential_carbon_wc_r, land_s1_wc_nb)
carbon_ratio_i <- cellStats(animal_potential_carbon_wc_r_i, "sum", na.rm = TRUE) / cellStats(animal_potential_carbon_wc_r, "sum", na.rm = TRUE)

#NCP
ncp_wc_r_i <- mask(ncp_wc_r, land_s1_wc_nb)
ncp_ratio_i <- cellStats(ncp_wc_r_i, "sum", na.rm = TRUE) / cellStats(ncp_wc_r, "sum", na.rm = TRUE)

eval_p05_land_wc_nb$bio_held <- c(bio_ratio_i)
eval_p05_land_wc_nb$carbon_held <- c(carbon_ratio_i)
eval_p05_land_wc_nb$ncp_held <- c(ncp_ratio_i)


for (i in seq(from = 5, to = 100, by = 5)) {
  land_budget_nobio <- (i * 0.01) * cellStats(land_wc_r, "sum", na.rm = TRUE)
  
  land_p1_wc_i <- problem(x = land_wc_r, features = existing_sheep_wc_r_inv) %>%
    add_min_shortfall_objective(land_budget_nobio) %>% 
    add_relative_targets(1) %>% 
    add_gurobi_solver(gap = 0.01)
  
  print(paste0("Just checked the problem set for ", i))
  
  land_s1_wc_i <- solve(land_p1_wc_i)
  
  eval_p05_land_wc_nb_i <- eval_target_coverage_summary(land_p1_wc_i,land_s1_wc_i)
  
  land_s1_wc_i[land_s1_wc_i < 1] <- NA
  
  #Bio
  biodiversity_wc_r_i <- mask(biodiversity_wc_r, land_s1_wc_i)
  bio_ratio_i <- cellStats(biodiversity_wc_r_i, "sum", na.rm = TRUE) / cellStats(biodiversity_wc_r, "sum", na.rm = TRUE)
  
  #Carbon
  animal_potential_carbon_wc_r_i <- mask(animal_potential_carbon_wc_r, land_s1_wc_i)
  carbon_ratio_i <- cellStats(animal_potential_carbon_wc_r_i, "sum", na.rm = TRUE) / cellStats(animal_potential_carbon_wc_r, "sum", na.rm = TRUE)
  
  #NCP
  ncp_wc_r_i <- mask(ncp_wc_r, land_s1_wc_i)
  ncp_ratio_i <- cellStats(ncp_wc_r_i, "sum", na.rm = TRUE) / cellStats(ncp_wc_r, "sum", na.rm = TRUE)
  
  eval_p05_land_wc_nb_i$budgetpct <- c(i)
  eval_p05_land_wc_nb_i$bio_held <- c(bio_ratio_i)
  eval_p05_land_wc_nb_i$carbon_held <- c(carbon_ratio_i)
  eval_p05_land_wc_nb_i$ncp_held <- c(ncp_ratio_i)
  
  eval_p05_land_wc_nb <- rbind(eval_p05_land_wc_nb,eval_p05_land_wc_nb_i)
  
  filename <- paste0("Wool/land_s1_wc_", i, "pct_NoBioNoC.tif")
  
  writeRaster(
   land_s1_wc_i, filename)
  
}


melt()
eval_p05_land_wc_nb_melt <- melt(eval_p05_land_wc_nb[10:13], id.vars="budgetpct")

ggplot(eval_p05_land_wc_nb_melt, aes(x = budgetpct, y = value, color = variable)) +
  geom_abline(intercept = 0, slope = 0.01, linetype = 2) +
  scale_colour_manual(values = c("#16c6dd", "#8ec92f", "#9391ff"), name = "Feature",
                      labels = c("Biodiversity", "Carbon", "NCP")) +
  theme_minimal() + 
  theme(text = element_text(family = "Corbel", size = 16)) +
  labs(x = "Land Budget Percentage", title = "Wool Scenario 2 (Contraction)", y = "Ratio of Feature Retained", ) +
  geom_line(size = 1)


ggsave("Wool/feature_curves_WoolWC_noBioC_inverse.png",
       width = 6, height = 3, units = "in",
       bg = "white")

ggsave("Wool/feature_curves_WC_noBioC_inverse.png",
       width = 6, height = 3, units = "in",
       bg = "white")


write.csv(eval_p05_land_wc_nb_melt, 'Wool/eval_p05_land_wc_nb_inverse.csv')

write.csv(eval_p05_land_wc_nb, 'Wool/eval_p05_land_wc_nb_inverse.csv')


