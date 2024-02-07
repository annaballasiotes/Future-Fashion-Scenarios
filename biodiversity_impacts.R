## Extract biodiversity data from scenario footprints
## Last updated: 01/30/2023

library(tidyverse)
library(Matrix)
library(terra)
library(sf)
library(rnaturalearth)
library(extrafont)
library(showtext)
library(reshape2)
library(tools)
loadfonts(device = "win")

wd <- #Set working directory to where repo is stored
  
setwd(wd)

template <- rast('template_raster_10k_ek4.tif')


# NOTE: This biodiversity data folder is where IUCN data is located
# This data cannot be freely distributed, and is excluded from this repository

biodiv_folder <- "Biodiversity/biodiversity_data/"

# Read in biodiversity matrices
taxas <- c("amphibians", "birds", "mammals", "reptiles")

# Read in IUCN data
iucn <- read_csv(paste0(biodiv_folder, "all_verts_aoh_10km_iucn_info_list.csv")) %>%
  dplyr::select(species, common_name, taxonomic_group, iucn_category)

for (i in seq_along(taxas)){
  taxa <- taxas[i]
  mat <- readRDS(paste0(biodiv_folder, taxa, "_redlist_range_10km.rds"))

  if(i == 1){mats <- mat} else{mats <- cbind(mats, mat)}
}


# Create Footprints from All Solutions
## This ensures that only the "1" value (solution) is stored
## to eliminate any data errors

cash_files <- list.files("Cashmere", pattern="*.tif$", full.names = FALSE)
cott_files <- list.files("Cotton", pattern="*.tif$", full.names = FALSE)
wool_files <- list.files("Wool", pattern="*.tif$", full.names = FALSE)

for (i in cash_files){
  cash_name <- str_sub(i, start = 0, end = -5)
  filename1 <- paste0("Cashmere/", cash_name,".tif")
  cash_i <- rast(filename1)
  cash_i[cash_i<1] <- NA
  filename2 <- paste0("Biodiversity/", cash_name, ".tif")
  writeRaster(cash_i, filename2)
}

for (i in cott_files){
  cott_name <- str_sub(i, start = 0, end = -5)
  filename1 <- paste0("Cotton/", cott_name, ".tif")
  cott_i <- rast(filename1)
  cott_i[cott_i<1] <- NA
  filename2 <- paste0("Biodiversity/", cott_name,".tif")
  writeRaster(cott_i, filename2)
}

for (i in wool_files){
  wool_name <- str_sub(i, start = 0, end = -5)
  filename1 <- paste0("Cotton/", wool_name, ".tif")
  wool_i <- rast(filename1)
  wool_i[wool_i<1] <- NA
  filename2 <- paste0("Biodiversity/", wool_name, ".tif")
  writeRaster(wool_i, filename2)
}

fp_files <- list.files("Biodiversity", pattern="*.tif$", full.names = TRUE)

scenarios <- str_sub(fp_files, start = 13, end = -5)


# Extract biodiversity data for all footprints
# This requires the IUCN data
for (i in seq_along(fp_files)){
  fp_file <- fp_files[i]
  scenario <- str_sub(fp_file, start = 18, end = -5) #May need to change start depending on wd
  print(paste0("Scenario: ", scenario))
  
  rast <- rast(fp_file) %>% 
    project(crs(template),
            method = "near") %>% 
    extend(template)
  
  fp_cells <- cells(rast)
  
  fp_spec_ids <- which(colSums(mats[fp_cells, , drop = FALSE]) > 0)
  
  # restrain to those w/ more than 20% of their range exposed
  fp_within_ids <- which(colSums(mats[fp_cells, fp_spec_ids, drop = FALSE])/
                           colSums(mats[ , fp_spec_ids, drop = FALSE]) > 0.2)
  
  gc()
  
  fp_spec <- colnames(mats[, fp_spec_ids])[fp_within_ids]
  
  fp_spec_exp <- (
    colSums(mats[fp_cells, fp_spec_ids, drop = FALSE])/
      colSums(mats[ , fp_spec_ids, drop = FALSE])
  )[fp_within_ids]
  
  if (length(fp_spec) > 0){
    
    fp_df <- data.frame(species = fp_spec, 
                        # country = country,
                        scenario = scenario) %>% 
      left_join(iucn, by = "species") %>% 
      mutate(pct_exposed = fp_spec_exp,
             range_km2 = colSums(mats[ , fp_spec_ids, drop = FALSE])[fp_within_ids])
    
    gc()
    
    fp_df_sum <- fp_df %>% 
      group_by(scenario, 
               # country, 
               taxonomic_group, iucn_category) %>% 
      summarize(n_species = n(), 
                .groups = "drop") 
    
    df_ALL <- bind_rows(df_ALL, fp_df)
    df_ALL_sum <- bind_rows(df_ALL_sum, fp_df_sum)
    
  }
}


write_csv(df_ALL, "species_by_scenario.csv")
write_csv(df_ALL_sum, "species_summary_by_scenario.csv")



summary_results <- df_ALL %>% 
  group_by(scenario, iucn_category) %>% 
  summarize(n_species = sum(n_species),
            .groups = "drop") %>% 
  dplyr::filter(!iucn_category %in% c("LR/cd", "LR/lc", "NA")) %>% 
  drop_na() %>% 
  mutate(iucn_category = factor(
    iucn_category,
    levels = c("DD", "LC", "NT", "VU", "EN", "CR")
  ))



contraction_summary <- write_csv("Biodiversity/summary_results.csv")


# Function to extract rows based on a string in a specific column
extractRowsByString <- function(data, column, substring) {
  # Check if the column exists in the data frame
  if (!(column %in% names(data))) {
    stop("Column not found in the data frame.")
  }
  
  # Extract rows where the specified string is present in the specified column
  selectedRows <- data[grepl(substring, data[[column]]), ]
  
  return(selectedRows)
}

cotton_contraction = extractRowsByString(summary_results, "scenario", "_cc")

wool_contraction = extractRowsByString(summary_results, "scenario", "wc")

cashmere_contraction = extractRowsByString(summary_results, "scenario", "cshc")

contraction_summary_land <- contraction_summary %>% filter(type == "land") # May need to add "land" column for relevant scenarios,
                                                                           # if production AND land scenarios were run

contraction_summary_land$commodity <- factor(contraction_summary_land$commodity, levels = c("cashmere", "cotton", "wool"), 
                  labels = c("Cashmere", "Cotton", "Wool"))

ggplot(contraction_summary_land, aes(x = pct, y = total_species, color = scenario)) +
  #geom_point() +
  #geom_abline(intercept = 0, slope = 0.01, linetype = 2) +
  scale_colour_manual(values = c("#8ec92f", "#bebada"), name = "Scenario",
                      labels = str_wrap(c("Prioritizing Nature", "Ignoring Nature"), width = 20)) +
  theme_minimal() + 
  facet_wrap(vars(commodity), scales = "free") + 
  theme(text = element_text(family = "Corbel", size = 16)) +
  labs(x = "% Land Footprint Removed", title = "Red List Species (Biodiversity)", y = "Species Count", ) +
  geom_line(size = 1)



