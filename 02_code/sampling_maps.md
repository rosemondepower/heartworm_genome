# Sampling maps

```R
# Map samples

## Load packages
library(ggplot2)
library(dplyr)
library(purrr)
library(stats)
library(graphics)
library(grDevices)
library(utils)
library(datasets)
library(methods)
library(base)
require(maps)
library(mapdata)
library(tidyverse)
library(readxl)
library(ozmaps) 
library(grid)
library(RColorBrewer)
library(gridExtra)
library(ggrepel)



#############################################################################
# Current map
## Samples I already have data for
#############################################################################

# Set working directory
setwd("C:/Users/rpow2134/OneDrive - The University of Sydney (Staff)/Documents/HW_WGS/R_analysis/batch4/map")

#-- metadata that describes information about the samples, such as country of origin, and GPS coordinates
location_file <- "location_map.csv"

# Read the actual data into R
location <- read.csv(location_file, header = TRUE)

## Make world map data
world_map <- map_data("world")

# Set colors for the points
red_palette1 <- brewer.pal(n = 9, name = "YlOrRd")
red_palette2 <- brewer.pal(n=9, name = "YlOrBr")
blue_palette <- brewer.pal(n = 9, name = "Blues")
green_palette1 <- brewer.pal(n = 9, name = "BuGn")
green_palette2 <- brewer.pal(n=9, name = "YlGn")


scale_colour_pop <-
      c('USA: Wisconsin' = red_palette2[8],
        'USA: Illinois' = red_palette1[8],
        'USA: Missouri' = red_palette2[7],
        'USA: Texas' = red_palette1[7],
        'USA: Louisiana' = red_palette1[6],
        'USA: Georgia' = red_palette2[5],
        'USA: Florida' = red_palette1[5],
        'CRI: San Jose' = "purple4",
        'PAN: Boca Chica' = "darkviolet",
        'PAN: Puerto Armuelles' = "blueviolet",
        'PAN: San Lorenzo' = "mediumpurple1",
        'ITA: Northeast' = green_palette1[9], 
        'ITA: Pavia' =  green_palette2[8],
        'ROU: Bucharest' = green_palette1[7],
        'ROU: Comana' =  green_palette2[7],
        'ROU: Giurgiu' = green_palette2[5],
        'GRC: Thess/Xanthi' = green_palette1[5],
        'VNM: Hai Phong' = "deeppink3",
        'VNM: Ha Noi' = "deeppink1",
        'THA: Bangkok' = "violetred1",
        'MYS: Selangor' = "hotpink",
        'AUS: Lockhart River' = blue_palette[9],
        'AUS: Cairns' = blue_palette[8],
        'AUS: Townsville' = blue_palette[7],
        'AUS: Rockhampton' = blue_palette[6],
        'AUS: Brisbane' = blue_palette[5],
        'AUS: Sydney' = blue_palette[4])



# Create a ggplot for the world map
world_plot <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey90") +
  geom_point(data = location, aes(x = longitude, y = latitude, color = city, shape = collection), size = 3) + # Map "collected" to shape
  geom_text_repel(data = location, aes(x = longitude, y = latitude, label = label), size = 5, fontface = "bold", max.overlaps = 20, force = 4) +  # Add labels
  scale_color_manual(values = scale_colour_pop, limits = c("USA: Wisconsin", "USA: Illinois", "USA: Missouri", "USA: Texas", "USA: Louisiana", "USA: Georgia", "USA: Florida",
"CRI: San Jose", "PAN: Boca Chica", "PAN: Puerto Armuelles", "PAN: San Lorenzo",
"ITA: Northeast", "ITA: Pavia", "ROU: Bucharest", "ROU: Comana", "ROU: Giurgiu", "GRC: Thess/Xanthi",
"VNM: Hai Phong", "VNM: Ha Noi", "THA: Bangkok", "MYS: Selangor",
"AUS: Lockhart River", "AUS: Cairns", "AUS: Townsville", "AUS: Rockhampton", "AUS: Brisbane", "AUS: Sydney")) +
  scale_shape_manual(values = c("This study" = 19, "Public data" = 17), breaks = c("This study", "Public data"), labels = c("This study", "Public data")) +  # Define shape values for samples we did/did not collect
  theme_void() +
  guides(
    color = guide_legend(override.aes = list(size = 5), title = "Location", keywidth = unit(0.5, "lines"), keyheight = unit(0.5, "lines"), order = -1, ncol = 1),  # Add title to the Location legend and set order to -1
    shape = guide_legend(override.aes = list(size = 5), title = "Collection", order = 1)  # Set order for "Collected" legend to 1
  ) +
  theme(
    legend.position = "right",  # Move the legend to the right
    legend.direction = "vertical",  # Vertical layout for the legend
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),  # Decrease the font size of the legend text
    legend.background = element_rect(color = NA, fill = "white"),  # Add a box around the legend
    legend.box.spacing = unit(0.2, "lines"),  # Adjust the spacing around the legend box
    plot.margin = unit(c(1, 1, 1, 4), "lines")  # Adjust the spacing around the plot (top, right, bottom, left)
  )



# Display the world plot
print(world_plot)

# Save plot
ggsave("map.png", height=8, width=16, bg = "white", dpi = 300)
ggsave("map.tif", height=8, width=16, bg = "white", dpi = 300)
ggsave("map.pdf", height=8, width=16, bg = "white", dpi = 300)






#################################
# AUSTRALIA MAP
#################################

# Add map inset to zoom in on Australian samples
# Manually specify the coordinates for the area of the world map to show in the inset
aus_xmin <- 141
aus_xmax <- 155
aus_ymin <- -60
aus_ymax <- -9

# Filter the world map data for the inset area
aus_data <- subset(world_map, long >= aus_xmin & long <= aus_xmax & lat >= aus_ymin & lat <= aus_ymax)

# Aus location data
location_aus <- read.csv("location_map_aus.csv", header = TRUE)


# Create a ggplot for the inset
aus_plot <- ggplot() +
  geom_polygon(data = aus_data, aes(x = long, y = lat, group = group), fill = "grey85") +
  geom_point(data = location_aus, aes(x = longitude, y = latitude, color = city, shape = Collection), size = 6) +
  geom_text_repel(data = location_aus, aes(x = longitude, y = latitude, label = label), size = 5, fontface = "bold",  vjust = 0, hjust = -0.1) +  # Add labels
  scale_color_manual(values = scale_colour_pop) +
  scale_shape_manual(values = c("This study" = 19, "Public data" = 17), breaks = c("This study", "Public data"), labels = c("This study", "Public data")) +  # Define shape values for samples we did/did not collect
  coord_fixed(ratio=1) +
  theme_void() +
  guides(color = "none") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.25, 0.13), 
    legend.title = element_text(size = 18),
    legend.direction = "vertical",  # Vertical layout for the legend
    legend.text = element_text(size = 16),  # Decrease the font size of the legend text
    legend.background = element_rect(color = NA, fill = NA),  # Add a box around the legend
    legend.box.spacing = unit(0.2, "lines"),  # Adjust the spacing around the legend box
  )

# Add title "Australia" above the inset_plot
#aus_plot <- aus_plot +
  #annotate("text",
         #  x = min(aus_data$long) + 0.5,   # Adjust x-coordinate for inside placement
        #  y = max(aus_data$lat) - 0.5,    # Adjust y-coordinate for inside placement
        #  label = "AUSTRALIA",
        #  hjust = 0, vjust = 1,           # Align the text appropriately
        #  size = 8, fontface = "bold")

print(aus_plot)

# Save plot
ggsave("map_aus.png", height=10, width=5, dpi = 300)
ggsave("map_aus.tif", height=10, width=5, dpi = 300)
ggsave("map_aus.pdf", height=10, width=5, dpi = 300)

```