## Cashmere Contraction
## Last updated: 12/19/2023

library(readxl)
library(tidyverse)
library(prioritizr)
library(rgdal)
library(rcbc)
library(terra)
library(highs)
library(rnaturalearth)
library(extrafont)
library(showtext)
library(reshape2)


wd <- #Set working directory to where repo is stored
  
setwd(wd)

#Code for Cashmere Contraction: cshc

#Need:
# 1. Land Footprint of Existing Cashmere (goats, GLWv4)
# 2. Target Layers: biodiversity (species richness), carbon potential, NCP (local)

#1.1 Existing Cashmere
existing_goats_cshc <- rast('Cashmere/Goats_DA_heads_ek4_10ksnap.tif')

#1.2 Land template
land_template_cshc <- rast('template_raster_10k_ek4.tif')


#2.1 Potential Carbon
animal_potential_carbon_cshc <- rast('Animal_COC_ek4_10k_area_snap_B1.tif')
names(animal_potential_carbon_cshc) <- "Carbon"

#2.2 Biodiversity
biodiversity_cshc <- rast('vert_rsr_richness.tif') %>% 
  project(crs(land_template_cshc)) # CB
names(biodiversity_cshc) <- "Biodiversity"

#2.3 NCP
ncp_cshc <- rast('local_NCP_all_targets_ek4_10k_snap.tif')
names(ncp_cshc) <- "NCP"


#####################
## Additional Prep For Prioritizr
#####################

#Cropping template to least common extent
template_cshc <- crop(biodiversity_cshc, animal_potential_carbon_cshc)
template_cshc <- crop(template_cshc, biodiversity_cshc)
template_cshc <- crop(template_cshc, existing_goats_cshc)
template_cshc <- crop(template_cshc, ncp_cshc)

#Cropping each layer to template
existing_goats_cshc <- crop(existing_goats_cshc, template_cshc)
biodiversity_cshc <- crop(biodiversity_cshc, template_cshc)
animal_potential_carbon_cshc <- crop(animal_potential_carbon_cshc, template_cshc)
land_cshc <- crop(land_template_cshc, template_cshc)
ncp_cshc <- crop(ncp_cshc, template_cshc)

#Scaling Because Prioritizr Throws Errors

existing_goats_cshc <- existing_goats_cshc * 0.01

animal_potential_carbon_cshc <- animal_potential_carbon_cshc * 0.000000001
animal_potential_carbon_cshc[animal_potential_carbon_cshc < 0.001] <- 0

biodiversity_cshc[biodiversity_cshc < 0.001] <- 0

# Masking Biodiversity & Carbon to Goat Footprint
animal_potential_carbon_cshc <- mask(animal_potential_carbon_cshc, existing_goats_cshc)
biodiversity_cshc <- mask(biodiversity_cshc, existing_goats_cshc)
land_cshc <- mask(land_cshc, existing_goats_cshc)
ncp_cshc <- mask(ncp_cshc, existing_goats_cshc)

#Convert to "raster" because Prioritizr handles it better 
existing_goats_cshc_r <- raster(existing_goats_cshc)
biodiversity_cshc_r <- raster(biodiversity_cshc)
animal_potential_carbon_cshc_r <- raster(animal_potential_carbon_cshc)
land_cshc_r <- raster(land_cshc)
ncp_cshc_r <- raster(ncp_cshc)

#Stack Target Features
feats_cshc_r <- stack(biodiversity_cshc_r, animal_potential_carbon_cshc_r, ncp_cshc_r)

#Plot Check
plot(feats_cshc_r)

###########################################
###########################################
###########################################

#####################
## LAND BUDGET - PRIORITIZING NATURE
#####################

#1 pct of land
land_budget_cshc <- 0.01 * cellStats(land_cshc_r, "sum", na.rm = TRUE)

# Set up prioritizr problem
land_p1_cshc <- problem(x = land_cshc_r, features = feats_cshc_r) %>%
  add_min_shortfall_objective(land_budget_cshc) %>% 
  add_relative_targets(c(1, 1, 1)) %>% 
  add_gurobi_solver(gap = 0.01) ## AQUIRE GUROBI SOLVER -- IF IN ACADEMIA, THIS IS FREE
## If no access to gurobi solver, use CBC solver, although this solver is significantly
## if not prohibitively slower https://prioritizr.net/reference/index.html#solvers 

# Solve prioritizr problem
land_s1_cshc <- solve(land_p1_cshc)

# Plot solution
plot(land_s1_cshc)

# Confirm pixel counts
freq(land_s1_cshc)

# Create summary eval table
eval_p05_land_cshc <- eval_target_coverage_summary(land_p1_cshc,land_s1_cshc)

writeRaster(
  land_s1_cshc, "Cashmere/land_s1_cshc_1pct.tif", overwrite = TRUE)


#####################
## LAND BUDGET LOOP - CASHMERE
#####################

# Add the percent column to the eval table
eval_p05_land_cshc$budgetpct <- c(1)

# Step wise, from 5% to 100% by 5
for (i in seq(from = 5, to = 100, by = 5)) {
  land_budget_cshc <- (i * 0.01) * cellStats(land_cshc_r, "sum", na.rm = TRUE)
  
  land_p1_cshc_i <- problem(x = land_cshc_r, features = feats_cshc_r) %>%
    add_min_shortfall_objective(land_budget_cshc) %>% 
    add_relative_targets(c(1,1,1)) %>% 
    add_gurobi_solver(gap = 0.01)
  
  print(paste0("Just checked the problem set for ", i))
  
  land_s1_cshc_i <- solve(land_p1_cshc_i)
  
  eval_p05_land_cshc_i <- eval_target_coverage_summary(land_p1_cshc_i,land_s1_cshc_i)
  
  eval_p05_land_cshc_i$budgetpct <- c(i)
  
  eval_p05_land_cshc <- rbind(eval_p05_land_cshc,eval_p05_land_cshc_i)
  
  filename <- paste0("Cashmere/land_s1_cshc", i, "pct.tif", overwrite = TRUE)
  
  writeRaster(
    land_s1_cshc_i, filename)
  
}

ggplot(eval_p05_land_cshc, aes(x = budgetpct, y = relative_held, color = feature)) +
  geom_abline(intercept = 0, slope = 0.01, linetype = 2) +
  scale_colour_manual(values = c("#16c6dd", "#8ec92f", "#9391ff"), name = "Feature",
                      labels = c("Biodiversity", "Carbon", "NCP")) +
  theme_minimal() + 
  theme(text = element_text(family = "Corbel", size = 16)) +
  labs(x = "Land Budget Percentage", title = "Cashmere Scenario 3 (Contraction + Nature)", y = "Ratio of Feature Retained", ) +
  geom_line(size = 1)

ggsave("Cashmere/feature_curves_cshc.png",
       width = 6, height = 3, units = "in",
       bg = "white")

write.csv(eval_p05_land_cshc, 'Cashmere/eval_p05_land_cshc_1219.csv')


#####################
## LAND BUDGET -- NO BIO / NO CARBON
#####################

existing_goats_cshc_r_inv <- 1/existing_goats_cshc_r

#1 pct of land
land_budget_cshc <- 0.01 * cellStats(land_cshc_r, "sum", na.rm = TRUE)

land_p1_cshc_nb <- problem(x = land_cshc_r, features = existing_goats_cshc_r_inv) %>%
  add_min_shortfall_objective(land_budget_cshc) %>% 
  add_relative_targets(c(1)) %>% 
  add_gurobi_solver(gap = 0.01)

land_s1_cshc_nb <- solve(land_p1_cshc_nb)

plot(land_s1_cshc_nb)
freq(land_s1_cshc_nb)

eval_p05_land_cshc_nb <- eval_target_coverage_summary(land_p1_cshc_nb,land_s1_cshc_nb)

writeRaster(
  land_s1_cshc_nb, "Cashmere/land_s1_cshc_1pct_NoBioNoC.tif", overwrite = TRUE)

#####################
## LAND BUDGET LOOP - NO BIO / NO CARBON - CASHMERE
#####################

# Add the percent column to the eval table
eval_p05_land_cshc_nb$budgetpct <- c(1)


land_s1_cshc_nb[land_s1_cshc_nb < 1] <- NA

# Add variable ratios to eval table 

#Bio
biodiversity_cshc_r_i <- mask(biodiversity_cshc_r, land_s1_cshc_nb)
bio_ratio_i <- cellStats(biodiversity_cshc_r_i, "sum", na.rm = TRUE) / cellStats(biodiversity_cshc_r, "sum", na.rm = TRUE)

#Carbon
animal_potential_carbon_cshc_r_i <- mask(animal_potential_carbon_cshc_r, land_s1_cshc_nb)
carbon_ratio_i <- cellStats(animal_potential_carbon_cshc_r_i, "sum", na.rm = TRUE) / cellStats(animal_potential_carbon_cshc_r, "sum", na.rm = TRUE)

#NCP
ncp_cshc_r_i <- mask(ncp_cshc_r, land_s1_cshc_nb)
ncp_ratio_i <- cellStats(ncp_cshc_r_i, "sum", na.rm = TRUE) / cellStats(ncp_cshc_r, "sum", na.rm = TRUE)

eval_p05_land_cshc_nb$bio_held <- c(bio_ratio_i)
eval_p05_land_cshc_nb$carbon_held <- c(carbon_ratio_i)
eval_p05_land_cshc_nb$ncp_held <- c(ncp_ratio_i)


for (i in seq(from = 5, to = 100, by = 5)) {
  land_budget_nobio <- (i * 0.01) * cellStats(land_cshc_r, "sum", na.rm = TRUE)

  land_p1_cshc_i <- problem(x = land_cshc_r, features = existing_goats_cshc_r_inv) %>%
    add_min_shortfall_objective(land_budget_nobio) %>% 
    add_relative_targets(1) %>% 
    add_gurobi_solver(gap = 0.01)
  
  print(paste0("Just checked the problem set for ", i))
  
  land_s1_cshc_i <- solve(land_p1_cshc_i)
  
  eval_p05_land_cshc_nb_i <- eval_target_coverage_summary(land_p1_cshc_i,land_s1_cshc_i)
  
  land_s1_cshc_i[land_s1_cshc_i < 1] <- NA
  
  #Bio
  biodiversity_cshc_r_i <- mask(biodiversity_cshc_r, land_s1_cshc_i)
  bio_ratio_i <- cellStats(biodiversity_cshc_r_i, "sum", na.rm = TRUE) / cellStats(biodiversity_cshc_r, "sum", na.rm = TRUE)
  
  #Carbon
  animal_potential_carbon_cshc_r_i <- mask(animal_potential_carbon_cshc_r, land_s1_cshc_i)
  carbon_ratio_i <- cellStats(animal_potential_carbon_cshc_r_i, "sum", na.rm = TRUE) / cellStats(animal_potential_carbon_cshc_r, "sum", na.rm = TRUE)
  
  #NCP
  ncp_cshc_r_i <- mask(ncp_cshc_r, land_s1_cshc_i)
  ncp_ratio_i <- cellStats(ncp_cshc_r_i, "sum", na.rm = TRUE) / cellStats(ncp_cshc_r, "sum", na.rm = TRUE)
  
  eval_p05_land_cshc_nb_i$budgetpct <- c(i)
  eval_p05_land_cshc_nb_i$bio_held <- c(bio_ratio_i)
  eval_p05_land_cshc_nb_i$carbon_held <- c(carbon_ratio_i)
  eval_p05_land_cshc_nb_i$ncp_held <- c(ncp_ratio_i)
  
  eval_p05_land_cshc_nb <- rbind(eval_p05_land_cshc_nb,eval_p05_land_cshc_nb_i)
  
  filename <- paste0("Cashmere/land_s1_cshc", i, "pct_NoBioNoC.tif")
  
  writeRaster(
    land_s1_cshc_i, filename)
  
}

melt()
eval_p05_land_cshc_nb_melt <- melt(eval_p05_land_cshc_nb[10:13], id.vars="budgetpct")

ggplot(eval_p05_land_cshc_nb_melt, aes(x = budgetpct, y = value, color = variable)) +
  #geom_point() +
  geom_abline(intercept = 0, slope = 0.01, linetype = 2) +
  scale_colour_manual(values = c("#16c6dd", "#8ec92f", "#9391ff"), name = "Feature",
                      labels = c("Biodiversity", "Carbon", "NCP")) +
  theme_minimal() + 
  theme(text = element_text(family = "Corbel", size = 16)) +
  labs(x = "Land Budget Percentage", title = "Cashmere Scenario 2 (Contraction)", y = "Ratio of Feature Retained", ) +
  geom_line(size = 1)


ggsave("Cashmere/feature_curves_cshc_noBioC.png",
       width = 6, height = 3, units = "in",
       bg = "white")

write.csv(eval_p05_land_cshc_nb_melt, 'Cashmere/eval_p05_land_cshc_nb_melt.csv')

write.csv(eval_p05_land_cshc_nb, 'Cashmere/eval_p05_land_cshc_nb.csv')

