## Cotton Reduction
## Last updated: 12/20/2024

library(readxl)
library(tidyverse)
library(prioritizr) #https://prioritizr.net/articles/prioritizr.html
library(rgdal)
library(rcbc)
library(terra)
library(highs)
library(rnaturalearth)
library(ggplot2)
library(gurobi)
library(extrafont)
library(showtext)
library(reshape2)

loadfonts(device = "win")

wd <- #Set working directory to where repo is stored
  
setwd(wd)

#Suffix for Cotton Reduction: cc

#Need:
# 1. Land Footprint of Existing Cotton (SPAM)
# 2. Target Layers: biodiversity (species richness), carbon potential, NCP (local)

#1.1 Existing Cotton
existing_cotton_cc <- rast('Cotton/CottonSpamGT0_P_ek410ksnap.tif')

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
template_cc <- crop(biodiversity_cc, animal_potential_carbon_cc)
template_cc <- crop(template_cc, biodiversity_cc)
template_cc <- crop(template_cc, ncp_cc)

#Cropping each layer to template
existing_cotton_cc <- crop(existing_cotton_cc, template_cc)
biodiversity_cc <- crop(biodiversity_cc, template_cc)
animal_potential_carbon_cc <- crop(animal_potential_carbon_cc, template_cc)
land_cc <- crop(land_template, template_cc)
ncp_cc <- crop(ncp_cc, template_cc)

#Scaling Because Prioritizr Throws Errors

#1.1 Cotton
existing_cotton_cc <- existing_cotton_cc * 0.01

#2.1 Carbon
animal_potential_carbon_cc <- animal_potential_carbon_cc * 0.000000001
animal_potential_carbon_cc[animal_potential_carbon_cc < 0.001] <- 0

#2.2 Biodiversity
biodiversity_cc[biodiversity_cc < 0.001] <- 0

# Masking Biodiversity & Carbon To the Cotton Footprint
animal_potential_carbon_cc <- mask(animal_potential_carbon_cc, existing_cotton_cc)
biodiversity_cc <- mask(biodiversity_cc, existing_cotton_cc)
land_cc <- mask(land_cc, existing_cotton_cc)
ncp_cc <- mask(ncp_cc, existing_cotton_cc)

#Convert to "raster" because Prioritizr handles it better 
existing_cotton_cc_r <- raster(existing_cotton_cc)
biodiversity_cc_r <- raster(biodiversity_cc)
animal_potential_carbon_cc_r <- raster(animal_potential_carbon_cc)
land_cc_r <- raster(land_cc)
ncp_cc_r <- raster(ncp_cc)

#Stack Target Features
feats_cc_r <- stack(biodiversity_cc_r, animal_potential_carbon_cc_r, ncp_cc_r)

#Plot Check
plot(feats_cc_r)


###########################################
###########################################
###########################################


#####################
## LAND BUDGET - PRIORITIZING NATURE
#####################

# 1 pct of land footprint
land_budget <- 0.01 * cellStats(land_cc_r, "sum", na.rm = TRUE)

# Set up prioritizr problem
land_p1_cc <- problem(x = land_cc_r, features = feats_cc_r) %>%
  add_min_shortfall_objective(land_budget) %>% 
  add_relative_targets(c(1,1,1)) %>% 
  add_gurobi_solver(gap = 0.01) ## AQUIRE GUROBI SOLVER -- IF IN ACADEMIA, THIS IS FREE
  ## If no access to gurobi solver, use CBC solver, although this solver is significantly
  ## if not prohibitively slower https://prioritizr.net/reference/index.html#solvers 

# Solve prioritizr problem
land_s1_cc <- solve(land_p1_cc)

# Plot solution
plot(land_s1_cc)

# Confirm pixel counts
freq(land_s1_cc)

# Create summary eval table
eval_p05_land_cc <- eval_target_coverage_summary(land_p1_cc,land_s1_cc)

writeRaster(
  land_s1_cc, "Cotton/land_s1_cc_1pct.tif")

#####################
## LAND BUDGET LOOP -- PRIORITIZING NATURE
#####################

# Add the percent column to the eval table
eval_p05_land_cc$budgetpct <- c(1)

# Step wise, from 5% to 100% by 5
for (i in seq(from = 5, to = 100, by = 5)) {
  land_budget <- (i * 0.01) * cellStats(land_cc_r, "sum", na.rm = TRUE)
  
  land_p1_cc_i <- problem(x = land_cc_r, features = feats_cc_r) %>%
    add_min_shortfall_objective(land_budget) %>% 
    add_relative_targets(c(1,1,1)) %>% 
    add_gurobi_solver(gap = 0.01)
  
  print(paste0("Just checked the problem set for ", i))
  
  land_s1_cc_i <- solve(land_p1_cc_i)
  
  eval_p05_land_i <- eval_target_coverage_summary(land_p1_cc_i,land_s1_cc_i)
  
  eval_p05_land_i$budgetpct <- c(i)
  
  eval_p05_land_cc <- rbind(eval_p05_land_cc,eval_p05_land_i)
  
  filename <- paste0("Cotton/land_s1_cc_", i, "pct.tif")
  
  writeRaster(
    land_s1_cc_i, filename)

  }
  
ggplot(eval_p05_land_cc, aes(x = budgetpct, y = relative_held, color = feature)) +
  geom_abline(intercept = 0, slope = 0.01, linetype = 2) +
  scale_colour_manual(values = c("#16c6dd", "#8ec92f", "#9391ff"), name = "Feature",
                      labels = c("Biodiversity", "Carbon", "NCP")) +
  theme_minimal() + 
  theme(text = element_text(family = "Corbel", size = 16)) +
  labs(x = "Land Budget Percentage", title = "Cotton Scenario 3 (Contraction + Nature)" , y = "Ratio of Feature Retained", ) +
  geom_line(size = 1)


ggsave("Cotton/feature_curves_CottonCC.png",
       width = 6, height = 3, units = "in",
       bg = "white")


write.csv(eval_p05_land_cc, 'Cotton/eval_p05_land_cotton_cc.csv')


#####################
## LAND BUDGET -- NO BIO / NO CARBON
#####################

land_budget <- 0.01 * cellStats(land_cc_r, "sum", na.rm = TRUE)

existing_cotton_cc_r_inv <- 1/existing_cotton_cc_r

existing_cotton_cc_r_inv <- existing_cotton_cc_r_inv * 0.01


land_p1_cc_noBioC <- problem(x = land_cc_r , features = existing_cotton_cc_r_inv) %>%
  add_min_shortfall_objective(land_budget) %>%
  add_relative_targets(1) %>%
  add_gurobi_solver(gap = 0.01)

land_s1_cc_noBioC <- solve(land_p1_cc_noBioC)

plot(land_s1_cc_noBioC)

freq(land_s1_cc_noBioC)

eval_p05_land_cc_noBio <- eval_target_coverage_summary(land_p1_cc_noBioC,land_s1_cc_noBioC)

writeRaster(
  land_s1_cc_noBioC, "Cotton/land_s1_1pct_NoBioNoC.tif")

#####################
## LAND BUDGET LOOP - NO BIO / NO CARBON
#####################

# Add the percent column to the eval table
eval_p05_land_cc_noBio$budgetpct <- c(1)

land_s1_cc_noBioC[land_s1_cc_noBioC < 1] <- NA

# Add variable ratios to eval table 

#Bio
biodiversity_cc_r_i <- mask(biodiversity_cc_r, land_s1_cc_noBioC)
bio_ratio_i <- cellStats(biodiversity_cc_r_i, "sum", na.rm = TRUE) / cellStats(biodiversity_cc_r, "sum", na.rm = TRUE)

#Carbon
animal_potential_carbon_cc_r_i <- mask(animal_potential_carbon_cc_r, land_s1_cc_noBioC)
carbon_ratio_i <- cellStats(animal_potential_carbon_cc_r_i, "sum", na.rm = TRUE) / cellStats(animal_potential_carbon_cc_r, "sum", na.rm = TRUE)

#NCP
ncp_cc_r_i <- mask(ncp_cc_r, land_s1_cc_noBioC)
ncp_ratio_i <- cellStats(ncp_cc_r_i, "sum", na.rm = TRUE) / cellStats(ncp_cc_r, "sum", na.rm = TRUE)

eval_p05_land_cc_noBio$bio_held <- c(bio_ratio_i)
eval_p05_land_cc_noBio$carbon_held <- c(carbon_ratio_i)
eval_p05_land_cc_noBio$ncp_held <- c(ncp_ratio_i)


for (i in seq(from = 5, to = 100, by = 5)) {
  land_budget_nobio <- (i * 0.01) * cellStats(land_cc_r, "sum", na.rm = TRUE)
  
  #not prioritizing bio/c
  
  land_p1_cc_i <- problem(x = land_cc_r, features = existing_cotton_cc_r_inv) %>%
    add_min_shortfall_objective(land_budget_nobio) %>% 
    add_relative_targets(1) %>% 
    add_gurobi_solver(gap = 0.05)
  
  print(paste0("Just checked the problem set for ", i))
  
  land_s1_cc_i <- solve(land_p1_cc_i)
  
  eval_p05_land_i <- eval_target_coverage_summary(land_p1_cc_i,land_s1_cc_i)

  land_s1_cc_i[land_s1_cc_i < 1] <- NA
  
  #Bio
  biodiversity_cc_r_i <- mask(biodiversity_cc_r, land_s1_cc_i)
  bio_ratio_i <- cellStats(biodiversity_cc_r_i, "sum", na.rm = TRUE) / cellStats(biodiversity_cc_r, "sum", na.rm = TRUE)
  
  #Carbon
  animal_potential_carbon_cc_r_i <- mask(animal_potential_carbon_cc_r, land_s1_cc_i)
  carbon_ratio_i <- cellStats(animal_potential_carbon_cc_r_i, "sum", na.rm = TRUE) / cellStats(animal_potential_carbon_cc_r, "sum", na.rm = TRUE)
  
  #NCP
  ncp_cc_r_i <- mask(ncp_cc_r, land_s1_cc_i)
  ncp_ratio_i <- cellStats(ncp_cc_r_i, "sum", na.rm = TRUE) / cellStats(ncp_cc_r, "sum", na.rm = TRUE)
    
  
  eval_p05_land_i$budgetpct <- c(i)
  eval_p05_land_i$bio_held <- c(bio_ratio_i)
  eval_p05_land_i$carbon_held <- c(carbon_ratio_i)
  eval_p05_land_i$ncp_held <- c(ncp_ratio_i)
  
  eval_p05_land_cc_noBio <- rbind(eval_p05_land_cc_noBio,eval_p05_land_i)
  
  filename <- paste0("Cotton/land_s1_cc_", i, "pct_NoBioNoC.tif")
  
  writeRaster(
   land_s1_cc_i, filename)
  
}

melt()
eval_p05_land_cc_noBio_melt <- melt(eval_p05_land_cc_noBio[10:13], id.vars="budgetpct")

# Potentially remove erroneous scenarios (Prioritizr issue)
#eval_p05_land_cc_noBio_melt <- eval_p05_land_cc_noBio_melt[-c(17, 18, 38, 39, 59, 60),]

ggplot(eval_p05_land_cc_noBio_melt, aes(x = budgetpct, y = value, color = variable)) +

  geom_abline(intercept = 0, slope = 0.01, linetype = 2) +
  scale_colour_manual(values = c("#16c6dd", "#8ec92f", "#9391ff"), name = "Feature",
                      labels = c("Biodiversity", "Carbon", "NCP")) +
  theme_minimal() + 
  theme(text = element_text(family = "Corbel", size = 16)) +
  labs(x = "Land Budget Percentage", title = "Cotton Scenario 2 (Contraction)", y = "Ratio of Feature Retained", ) +
  geom_line(size = 1)

write.csv(eval_p05_land_cc_noBio_melt, 'Cotton/eval_p05_land_cotton_cc_NB.csv')

write.csv(eval_p05_land_cc_noBio_melt, 'Cotton/eval_p05_land_cotton_cc_NB.csv')

ggsave("Cotton/feature_curves_CottonCC_noBioC.png",
       width = 6, height = 3, units = "in",
       bg = "white")


###########################################
###########################################
###########################################

#####################
## PRODUCTION BUDGET BASED - PRIORITIZING NATURE
#####################

#Making the budget as a production of just 1%
prod_budget_cc_r <- 0.01 * cellStats(existing_cotton_cc_r, "sum", na.rm = TRUE)

prod_p1_cc <- problem(x = existing_cotton_cc_r, features = feats_cc_r) %>%
  add_min_shortfall_objective(prod_budget_cc_r) %>%
  add_relative_targets(c(1, 1, 1)) %>%
  add_gurobi_solver(gap = 0.01)

prod_s1_cc <- solve(prod_p1_cc)

plot(prod_s1_cc)

eval_p05_prod_cc <- eval_target_coverage_summary(prod_p1_cc,prod_s1_cc)

writeRaster(
  prod_s1_cc, "Cotton/prod_s1_1pct_cc.tif", overwrite = TRUE)

eval_p05_prod_cc$budgetpct <- c(1)


#Land
prod_s1_cc_l <- mask(prod_s1_cc, land_cc_r)

land_ratio <- cellStats(prod_s1_cc_l, "sum", na.rm = TRUE) / cellStats(land_cc_r, "sum", na.rm = TRUE)

eval_p05_prod_cc$land_held <- c(land_ratio)

#####################
## PRODUCTION BUDGET LOOP - PRIORITIZING NATURE
#####################

for (i in seq(from = 5, to = 100, by = 5)) {
  prod_budget <- (i * 0.01) * cellStats(existing_cotton_cc_r, "sum", na.rm = TRUE)
  
  prod_p1_cc_i <- problem(x = existing_cotton_cc_r, features = feats_cc_r) %>%
    add_min_shortfall_objective(prod_budget) %>%
    add_relative_targets(c(1,1,1)) %>%
    add_gurobi_solver(gap = 0.01)
  
  print(paste0("Just checked the problem set for ", i))
  
  prod_s1_cc_i <- solve(prod_p1_cc_i)
  
  eval_p05_prod_i <- eval_target_coverage_summary(prod_p1_cc_i,prod_s1_cc_i)
  
  eval_p05_prod_i$budgetpct <- c(i)
  
  #Land
  prod_s1_cc_i <- mask(prod_s1_cc_i, land_cc_r)
  land_ratio_i <- cellStats(prod_s1_cc_i, "sum", na.rm = TRUE) / cellStats(land_cc_r, "sum", na.rm = TRUE)
  
  eval_p05_prod_i$land_held <- c(land_ratio_i)
  
  eval_p05_prod_cc <- rbind(eval_p05_prod_cc,eval_p05_prod_i)
  
  filename <- paste0("Cotton/prod_s1_cc", i, "pct.tif")
  
  writeRaster(
    prod_s1_cc_i, filename, overwrite = TRUE)
}

ggplot(eval_p05_prod_cc, aes(x = budgetpct, y = relative_held, color = feature)) +
  #geom_point() +
  geom_abline(intercept = 0, slope = 0.01, linetype = 2) +
  scale_colour_manual(values = c("#16c6dd", "#8ec92f", "#9391ff"), name = "Feature",
                      labels = c("Biodiversity", "Carbon", "NCP")) +
  theme_minimal() + 
  theme(text = element_text(family = "Corbel", size = 16)) +
  labs(x = "Production Budget Percentage", title = "Cotton Production Scenario 2 (Contraction + Nature)" , y = "Ratio of Feature Retained", ) +
  geom_line(size = 1)


ggsave("Cotton/feature_curves_CottonProdCC.png",
       width = 6, height = 3, units = "in",
       bg = "white")

write.csv(eval_p05_prod_cc, 'Cotton/eval_p05_prod_cotton_cc.csv')


#####################
## PRODUCTION BUDGET BASED - NO BIO / NO C
#####################

#Making the budget as a production of just 1%
prod_budget_cc_r <- 0.01 * cellStats(existing_cotton_cc_r, "sum", na.rm = TRUE)

existing_cotton_cc_r_inv <- 1/existing_cotton_cc_r

existing_cotton_cc_r_inv <- existing_cotton_cc_r_inv * 0.01

prod_p1_cc_nb <- problem(x = existing_cotton_cc_r, features = existing_cotton_cc_r_inv) %>%
  add_min_shortfall_objective(prod_budget_cc_r) %>%
  add_relative_targets(1) %>%
  add_gurobi_solver(gap = 0.01)

prod_s1_cc_nb <- solve(prod_p1_cc_nb)

eval_p05_prod_cc_nb <- eval_target_coverage_summary(prod_p1_cc_nb,prod_s1_cc_nb)

writeRaster(prod_s1_cc_nb, "Cotton/prod_s1_1pct_cc_NoBioNoC.tif", overwrite = TRUE)

#Bio
biodiversity_cc_r_i <- mask(biodiversity_cc_r, prod_s1_cc_nb)
bio_ratio_i <- cellStats(biodiversity_cc_r_i, "sum", na.rm = TRUE) / cellStats(biodiversity_cc_r, "sum", na.rm = TRUE)

#Carbon
animal_potential_carbon_cc_r_i <- mask(animal_potential_carbon_cc_r, prod_s1_cc_nb)
carbon_ratio_i <- cellStats(animal_potential_carbon_cc_r_i, "sum", na.rm = TRUE) / cellStats(animal_potential_carbon_cc_r, "sum", na.rm = TRUE)

#NCP
ncp_cc_r_i <- mask(ncp_cc_r, prod_s1_cc_nb)
ncp_ratio_i <- cellStats(ncp_cc_r_i, "sum", na.rm = TRUE) / cellStats(ncp_cc_r, "sum", na.rm = TRUE)

#Land
prod_s1_cc_nb_i <- mask(land_cc_r, prod_s1_cc_nb)
land_ratio_i <- cellStats(prod_s1_cc_nb_i, "sum", na.rm = TRUE) / cellStats(land_cc_r, "sum", na.rm = TRUE)


eval_p05_prod_cc_nb$bio_held <- c(bio_ratio_i)
eval_p05_prod_cc_nb$carbon_held <- c(carbon_ratio_i)
eval_p05_prod_cc_nb$ncp_held <- c(ncp_ratio_i)
eval_p05_prod_cc_nb$land_held <- c(land_ratio_i)

#####################
## PRODUCTION BUDGET LOOP - NO BIO / NO C
#####################
eval_p05_prod_cc_nb$budgetpct <- c(1)

prod_budget_cc_r <- 0.01 * cellStats(existing_cotton_cc_r, "sum", na.rm = TRUE)


for (i in seq(from = 5, to = 100, by = 5)) {
  prod_budget <- (i * 0.01) * cellStats(existing_cotton_cc_r, "sum", na.rm = TRUE)
  
  prod_p1_cc_i <- problem(x = existing_cotton_cc_r, features = existing_cotton_cc_r_inv) %>%
    add_min_shortfall_objective(prod_budget) %>%
    add_relative_targets(1) %>%
    add_gurobi_solver(gap = 0.01)
  
  print(paste0("Just checked the problem set for ", i))
  
  prod_s1_cc_i <- solve(prod_p1_cc_i)
  
  eval_p05_prod_i <- eval_target_coverage_summary(prod_p1_cc_i,prod_s1_cc_i)
  
  prod_s1_cc_i[prod_s1_cc_i < 1] <- NA
  
  #Bio
  biodiversity_cc_r_i <- mask(biodiversity_cc_r, prod_s1_cc_i)
  bio_ratio_i <- cellStats(biodiversity_cc_r_i, "sum", na.rm = TRUE) / cellStats(biodiversity_cc_r, "sum", na.rm = TRUE)
  
  #Carbon
  animal_potential_carbon_cc_r_i <- mask(animal_potential_carbon_cc_r, prod_s1_cc_i)
  carbon_ratio_i <- cellStats(animal_potential_carbon_cc_r_i, "sum", na.rm = TRUE) / cellStats(animal_potential_carbon_cc_r, "sum", na.rm = TRUE)
  
  #NCP
  ncp_cc_r_i <- mask(ncp_cc_r, prod_s1_cc_i)
  ncp_ratio_i <- cellStats(ncp_cc_r_i, "sum", na.rm = TRUE) / cellStats(ncp_cc_r, "sum", na.rm = TRUE)
  
  #Land
  prod_s1_cc_nb_i <- mask(land_cc_r, prod_s1_cc_i)
  land_ratio_i <- cellStats(prod_s1_cc_nb_i, "sum", na.rm = TRUE) / cellStats(land_cc_r, "sum", na.rm = TRUE)
  
  eval_p05_prod_i$budgetpct <- c(i)
  
  eval_p05_prod_i$bio_held <- c(bio_ratio_i)
  eval_p05_prod_i$carbon_held <- c(carbon_ratio_i)
  eval_p05_prod_i$ncp_held <- c(ncp_ratio_i)
  eval_p05_prod_i$land_held <- c(land_ratio_i)
  
  eval_p05_prod_cc_nb <- rbind(eval_p05_prod_cc_nb,eval_p05_prod_i)
  
  filename <- paste0("Cotton/prod_s1_cc_", i, "pct_NoBioNoC.tif")
  
  writeRaster(
     prod_s1_cc_i, filename, overwrite = TRUE)
}

melt()
eval_p05_prod_cc_nb_melt <- melt(eval_p05_prod_cc_nb[10:14], id.vars="budgetpct")

eval_p05_prod_cc_nb_melt <- eval_p05_prod_cc_nb_melt[-c(4, 25, 46, 67),]

ggplot(eval_p05_prod_cc_nb_melt, aes(x = budgetpct, y = value, color = variable)) +
  #geom_point() +
  geom_abline(intercept = 0, slope = 0.01, linetype = 2) +
  scale_colour_manual(values = c("#16c6dd", "#8ec92f", "#9391ff", "#000000"), name = "Feature",
                      labels = c("Biodiversity", "Carbon", "NCP", "Land")) +
  theme_minimal() + 
  theme(text = element_text(family = "Corbel", size = 16)) +
  labs(x = "Production Budget Percentage", title = "Cotton Production Scenario 2 (Contraction)", y = "Ratio of Feature Retained", ) +
  geom_line(size = 1)


ggsave("Cotton/feature_curves_Cotton_noBioC.png",
       width = 6, height = 3, units = "in",
       bg = "white")

write.csv(eval_p05_prod_cc_nb_melt, 'Cotton/eval_p05_prod_cc_nb_melt.csv')

write.csv(eval_p05_prod_cc_nb, 'Cotton/eval_p05_prod_cc_nb.csv')
