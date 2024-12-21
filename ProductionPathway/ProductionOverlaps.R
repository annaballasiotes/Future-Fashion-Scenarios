## Production Overlaps
## Last updated: 12/20/2024

library(tidyverse)
library(Matrix)
library(terra)
library(sf)
library(rnaturalearth)
library(tools)
library(scales)

# Cotton

Cotton_ProductionOverlap <- read_csv("Cotton_ProductionOverlap.csv")

Cotton_ProductionOverlap$cotton_kg_pct <- 100*(Cotton_ProductionOverlap$cotton_kg / 5.148972e+07)

ggplot(Cotton_ProductionOverlap, aes(x = pct, y = cotton_kg_pct, color = nature)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_colour_manual(values = c("#35D3E7", "#9A33BF"), name = "Pathway",
                      labels = c("Nature Focused", "Production Focused")) +
  theme_minimal() + 
  theme(text = element_text(family = "Corbel", size = 14)) +
  labs(x = "% Land Footprint Removed", title = "Cotton Production in Removed Land Footprint", y = "% Cotton Production Removed") +
  scale_y_continuous(labels = scales::comma) +
  geom_line(size = 1)

# Wool Sheep

Wool_SheepOverlap <- read_csv("Wool_SheepOverlap.csv")
Wool_SheepOverlap$sheep_heads_pct <- 100* (Wool_SheepOverlap$sheep_hds /5044846.72)

ggplot(Wool_SheepOverlap, aes(x = pct, y = sheep_heads_pct, color = nature)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_colour_manual(values = c("#35D3E7", "#9A33BF"), name = "Pathway",
                      labels = c("Nature Focused", "Production Focused")) +
  theme_minimal() + 
  theme(text = element_text(family = "Corbel", size = 14)) +
  labs(x = "% Land Footprint Removed", title = "Sheep in Removed Land Footprint", y = "% Sheep Removed") +
  scale_y_continuous(labels = scales::comma) +
  geom_line(size = 1)

# Cashmere Goats

Cshmr_GoatOverlap <- read_csv("Cashmere_GoatOverlap.csv")
Cshmr_GoatOverlap$goat_hds_pct <- 100 * (Cshmr_GoatOverlap$goat_hds / 1343837.52)


ggplot(Cshmr_GoatOverlap, aes(x = pct, y = goat_hds_pct, color = nature)) +
  #geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_colour_manual(values = c("#35D3E7", "#9A33BF"), name = "Pathway",
                      labels = c("Nature Focused", "Production Focused", "Identity")) +
  scale_linetype_manual(values = c("Identity Line" = "dotted")) + 
  labs(linetype = "Line Type") + 
  theme_minimal() + 
  theme(text = element_text(family = "Corbel", size = 14)) +
  labs(x = "% Land Footprint Removed", title = "Goats in Removed Land Footprint", y = "% Goats Removed") +
  scale_y_continuous(labels = scales::comma) +
  geom_line(size = 1)
