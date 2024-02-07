##Carbon Impacts
## Last updated: 12/20/2023

library(readxl)
library(tidyverse)
library(prioritizr)
library(rgdal)
library(rcbc)
library(terra)
library(highs)
library(rnaturalearth)
library(ggplot2)


wd <- #Set working directory to where repo is stored
  
setwd(wd)


# Methods
# 1) Take Animal Carbon COC (both median & minimum), multiply the carbon by the area per pixel of each commodity
## This has been completed and is provided in the "Carbon" folder
# 2) Resample scenarios to Animal Carbon COC resolution, mask animal carbon w/ scenarios 
## This was completed in Arc due to processing speed; this is not within this code but should be done for any relevant scenarios (can be batch processed)
# 3) Use this script to summarize carbon values


solutions <- list.files(workdir, pattern="*.tif$", full.names = FALSE)

# Return file names
solutions.names <- c( unlist( lapply(str_sub(solutions,1, -5), FUN=function(x) { x[1] })))   
print(solutions.names)


# Create table to populate

### TABLE WITH ANIMAL C - MIN

min_Carbon_COC <- data.frame(
  scenario = character(), 
  carbon_Gt = double(), 
  CO2 = double())


# Make rasters for each file
for (i in 1:length(solutions)){
  rasterforsc <- assign(solutions.names[i], raster(solutions[i]))
  scenario_name = solutions.names[i]
  carbonGt = cellStats(rasterforsc, "sum") * 6136.98325321 * 1e-9
  carbondioxide = carbonGt * 3.67
  fp_df <- data.frame(scenario = scenario_name, 
                      carbon_Gt = carbonGt,
                      CO2 = carbondioxide)
  min_Carbon_COC <- bind_rows(min_Carbon_COC, fp_df)
}

write.csv(min_Carbon_COC, "Carbon/min_Carbon_COC.csv")

### TABLE WITH ANIMAL C - MEDIAN

med_Carbon_COC <- data.frame(
  scenario = character(), 
  carbon_Gt = double(), 
  CO2 = double())


# Making rasters for each file
for (i in 1:length(solutions)){
  rasterforsc <- assign(solutions.names[i], raster(solutions[i]))
  scenario_name = solutions.names[i]
  carbonGt = cellStats(rasterforsc, "sum") * 6136.98325321 * 1e-9
  carbondioxide = carbonGt * 3.67
  fp_df <- data.frame(scenario = scenario_name, 
                      carbon_Gt = carbonGt,
                      CO2 = carbondioxide)
  med_Carbon_COC <- bind_rows(med_Carbon_COC, fp_df)
}

write.csv(med_Carbon_COC, "Carbon/med_Carbon_COC.csv")




