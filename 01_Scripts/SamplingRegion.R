# Load required libraries
library(ggplot2)
library(sf)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(osmdata)
library(tidyverse)
library(rnaturalearthhires)
library(ggmapinset)
library(here)

# Load your data (assuming data has lat and lot columns)

data <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/Combined.csv")

data_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326)

# Get Quebec boundary
quebec <- ne_states(country = "Canada", returnclass = "sf") %>%
  filter(name == "Québec")
canada <- ne_countries(country = "Canada", returnclass = "sf")
canada_us <- ne_countries(country = c("Canada", "United States of America"), returnclass = "sf")

ne_countries(returnclass = "sf")

#test
roads <- ne_download(scale = 10, type = "roads", category = "cultural")

cities <- ne_download(scale = 10, type = "populated_places", category = "cultural")

states <- ne_download(scale = 50, type = "states", category = "cultural")

c <- ne_download(scale = 50, type = "countries", category = "cultural")

test <- st_crop(roads, quebec)
testc <- st_crop(test, bbox)
plot(testc)
CLASSI_ECO_QC.gdb
# zones:
zones <- st_layers(here("zones-Raw", "CLASSI_ECO_QC_GDB", "CLASSI_ECO_QC.gdb"))
reg <- st_read(here("zones-Raw", "CLASSI_ECO_QC_GDB", "CLASSI_ECO_QC.gdb"),
               layer = "N3_DOM_BIO")

plot(reg$SREG_ECO)
r2 <- reg[2,]
r3 <- reg[3,]
r4 <- reg[4,]
r5 <- reg[5,]

r2t <- st_transform(r2, st_crs(bbox))
r3t <- st_transform(r3, st_crs(bbox))
r4t <- st_transform(r4, st_crs(bbox))
r5t <- st_transform(r5, st_crs(bbox))


# Define bounding box for the region of interest
bbox <- st_bbox(c(xmin = -78, xmax = -74.5, ymin = 45.5, ymax = 48.2), crs = 4326)

study_bbox <- st_as_sfc(st_bbox(c(xmin = -78, xmax = -74.5, ymin = 45.5, ymax = 48.2), crs = 4326))

# Crop the data to the bounding box
quebec_cropped <- st_crop(quebec, study_bbox)
roads_cropped <- st_crop(roads, study_bbox)
cities_cropped <- st_crop(cities, study_bbox)
states_cropped <- states[states$admin == "United States of America" | states$admin == "Canada",]
data_cropped <- st_crop(data_sf, bbox)
c_cropped <- c[c$NAME_EN == "United States of America" | c$NAME_EN == "Canada",]

r2c <- st_crop(r2t, bbox)
r3c <- st_crop(r3t, bbox)
r4c <- st_crop(r4t, bbox)
r5c <- st_crop(r5t, bbox)

states$admin[states$admin == "Canada"]
states$admin[states$admin == "United States of America"]
states$admin[states$admin == "United States of America" | states$admin == "Canada"]

# Plot map
main_map <- ggplot() +
  geom_sf(data = quebec_cropped, fill = "lightgray", color = "black")+
  geom_sf(data = r3c, fill = "orange", alpha = 0.3) + 
  geom_sf(data = r4c, fill = "yellow", alpha = 0.3) + 
  geom_sf(data = r5c, fill = "darkgreen", alpha = 0.3) +
  geom_sf(data = roads_cropped, color = "darkgrey", size = 0.5, alpha = 0.8) +
  geom_sf(data = data_cropped, color = "maroon", size = 1.5)  + 
  geom_sf(data = cities_cropped, shape = 15, color = "black", size = 2) +
  geom_text(data = cities_cropped, aes(x = LONGITUDE, y = LATITUDE, label = NAME_EN), 
          size = 3, hjust = c(0.45, 0.5, 0.45), vjust = c(-1.2, -1.2, 2), fontface = "bold") +
  #annotation_scale(location = "bl", pad_x = unit(0.55, "cm"), pad_y = unit(0.65, "cm")) +
  #annotation_north_arrow(location = "bl", style = north_arrow_fancy_orienteering, pad_x = unit(0.55, "cm"), pad_y = unit(1.15, "cm")) +
  labs(title = "Study Region",
       x = "", y = "") +
  theme_minimal()
 
# Inset map (Canada with study region highlighted)
inset_map <- ggplot() +
  geom_sf(data = canada, fill = "white", color = "black") +
  geom_sf(data = study_bbox, fill = "red", alpha = 0.5) +  # Highlight study region in red
  theme_void() +
  labs(title = "")

final_plot <- cowplot::ggdraw() +
  cowplot::draw_plot(main_map) + 
  cowplot::draw_plot(inset_map, x = 0.54, y = 0.59, width = 0.35, height = 0.35) +
  cowplot::draw_line(x = c(0.58, 0.58), y = c(0.59, 0.895), color = "black", size = 0.2) +
  cowplot::draw_line(x = c(0.58, 0.842), y = c(0.59, 0.59), color = "black", size = 0.2)

inset_map <- ggplot() +
  geom_sf(data = c_cropped, fill = "white", color = "black") +
  #geom_sf(data = states_cropped, fill = "white", color = "black") +
  geom_sf(data = study_bbox, color = "red", alpha = 0.5) +  # Highlight study region
  coord_sf(xlim = c(-90, -55), ylim = c(40, 57)) +  # Zoom in closer
  geom_rect(aes(xmin = -120, xmax = -60, ymin = 30, ymax = 60), 
            fill = NA, color = "black", linewidth = 0.1) +  # Box around inset
  theme_void() +
  labs(title = "")

combined_plot <- cowplot::plot_grid(main_map, inset_map, ncol = 2, rel_widths = c(3, 1))



final_plot <- cowplot::ggdraw(combined_plot) +
  cowplot::draw_plot(annotation_scale(location = "br"), x = 0.85, y = 0.1, width = 0.1, height = 0.1) +
  cowplot::draw_plot(annotation_north_arrow(location = "br", which_north = "true", style = north_arrow_fancy_orienteering()), 
            x = 0.85, y = 0.05, width = 0.1, height = 0.1)

final_plot


cowplot::ggdraw() +
  cowplot::draw_plot(main_map, x = 0, y = 0, width = 0.7, height = 1) +  # Main figure takes most space
  cowplot::draw_plot(inset_map, x = 0.7, y = 0.4, width = 0.3, height = 0.9)

