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
library(cowplot)

# Load your data (assuming data has lat and lot columns)

data <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/Combined.csv")

betall <- st_read("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Other/BetAllegh/data/shapefile/betualle.shp",
        layer = "betualle")
plot(ba_sf)
ba_sf <- st_as_sf(betall, coords = c("Longitude", "Latitude"), crs = 4326)

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

states <- ne_download(scale = 10, type = "states", category = "cultural")

c <- ne_download(scale = 50, type = "countries", category = "cultural")

test <- st_crop(roads, quebec)
testc <- st_crop(test, bbox)
plot(testc)

# zones:
zones <- st_layers(here("zones-Raw", "CLASSI_ECO_QC_GDB", "CLASSI_ECO_QC.gdb"))
reg <- st_read(here("zones-Raw", "CLASSI_ECO_QC_GDB", "CLASSI_ECO_QC.gdb"),
               layer = "N3_DOM_BIO")

quebec <- states[states$name == "Québec",]

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
bbox <- st_bbox(c(xmin = -78, xmax = -74.5, ymin = 45.8, ymax = 48.2), crs = 4326)

study_bbox <- st_as_sfc(st_bbox(c(xmin = -78, xmax = -74.5, ymin = 45.8, ymax = 48.2), crs = 4326))

# Crop the data to the bounding box
quebec_cropped <- st_crop(quebec, study_bbox)
roads_cropped <- st_intersection(roads, quebec_cropped)
cities_cropped <- st_crop(cities, study_bbox)
states_cropped <- states[states$admin == "United States of America" | states$admin == "Canada",]
data_cropped <- st_crop(data_sf, bbox)
c_cropped <- c[c$NAME_EN == "United States of America" | c$NAME_EN == "Canada",]

r2c <- st_crop(r2t, bbox)
r3c <- st_crop(r3t, bbox)
r4c <- st_crop(r4t, bbox)
r5c <- st_crop(r5t, bbox)
?st_crop
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

inset_map

combined_plot <- cowplot::plot_grid(main_map, inset_map, ncol = 2, rel_widths = c(3, 1))



final_plot <- cowplot::ggdraw(combined_plot) +
  cowplot::draw_plot(annotation_scale(location = "br"), x = 0.85, y = 0.1, width = 0.1, height = 0.1) +
  cowplot::draw_plot(annotation_north_arrow(location = "br", which_north = "true", style = north_arrow_fancy_orienteering()), 
                     x = 0.85, y = 0.05, width = 0.1, height = 0.1)

final_plot


cowplot::ggdraw() +
  cowplot::draw_plot(main_map, x = 0, y = 0, width = 0.7, height = 1) +  # Main figure takes most space
  cowplot::draw_plot(inset_map, x = 0.7, y = 0.4, width = 0.3, height = 0.9)

#### gginset
inset_cities <- data.frame(name = c("Montreal"), lat = c(45.5019), long = c(-73.5674)) %>% 
  mutate(LONG = long, LAT = lat) %>% 
  st_as_sf(coords = c("LONG", "LAT"), crs = 4326)
inset_regions <- data.frame(name=c("ONTARIO", "QUEBEC", "USA"), lat=c(49.441, 50.391, 34.586), long=c(-86.582, -74.016, -87.103))

main_map <- ggplot() +
  geom_sf(data = quebec_cropped, fill = "lightgray", color = "black", linewidth=0.3, alpha = 0.5) +
  geom_sf(data = r3c, fill = "orange", alpha = 0.3, show.legend = T) + 
  geom_sf(data = r4c, fill = "yellow", alpha = 0.3) + 
  geom_sf(data = r5c, fill = "lightgreen", alpha = 0.3) +
  geom_sf(data = roads_cropped, color = "darkgrey", size = 0.5, alpha = 2) +
  geom_sf(data = data_cropped, color = "maroon", size = 1.5)  + 
  geom_sf(data = cities_cropped, shape = 15, color = "black", size = 2) +
  geom_text(data = cities_cropped, aes(x = LONGITUDE, y = LATITUDE, label = NAME_EN), 
            size = 3, hjust = c(0.45, 0.3, 0.4), vjust = c(-1.1, -1.4, 2), fontface = "bold") +
  annotation_scale(location = "bl", pad_x = unit(0.55, "cm"), pad_y = unit(0.8, "cm"), style = "bar", width_hint=0.15) +
  annotation_north_arrow(location = "bl", style = north_arrow_fancy_orienteering, pad_x = unit(0.75, "cm"), pad_y = unit(1.4, "cm"),
                         height=unit(1, "cm"), width = unit(1, "cm")) +
  labs(title = "",
       x = "", y = "") +
  theme_minimal() 

inset_map <- ggplot() +
  #geom_sf(data = c_cropped, fill = "white", color = "black") +
  geom_sf(data = states_cropped, fill = "white", color = "black", alpha=0, linewidth = 0.1) +
  geom_sf(data = ba_sf, aes(fill="HR"), color="grey", show.legend = T) +
  scale_fill_manual(name=" ", values=c("HR" = "grey"), labels = c("HR" = "Historical range")) + #\n
  geom_sf(data = study_bbox, color = "red", alpha = 0.1, linewidth=0.5) +
  geom_sf(data = inset_cities, color = "black", size = 1.2, shape = 15) +
  geom_text(data = inset_regions, aes(x = long, y = lat, label =name), 
            size = 2.5, hjust = c(0.45, 0.4, 0.1), vjust = c(-1.1, -1.4, 2), fontface = "bold") +
  geom_text(data = inset_cities, aes(x = long, y = lat, label =name), 
            size = 1.8, hjust = c(0.0), vjust = c(-1.1), fontface = "bold") +
  coord_sf(xlim = c(-90, -58), ylim = c(30, 60)) +
  #annotation_scale(location = "br", pad_x = unit(0.3, "cm"), pad_y = unit(0.2, "cm"), style = "bar") +  # Zoom in closer
  #geom_rect(aes(xmin = -125, xmax = -65, ymin = 20, ymax = 70), 
  #          fill = NA, color = "black", linewidth = 0.1) +  # Box around inset
  labs(title = "") + 
  theme_void() +
  theme(legend.position = c(0.95, 0.02),
        legend.justification = c(1,0),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.4,"cm")) + # 0.3. with \n
  theme(panel.background = element_rect(fill = "white", color = "white")) 

inset_map
final_plot <- cowplot::ggdraw() +
  cowplot::draw_plot(main_map) + 
  cowplot::draw_plot(inset_map, x = 0.579, y = 0.5437, width = 0.4, height = 0.4) +
  cowplot::draw_line(x = c(0.656, 0.656), y = c(0.544, 0.902), color = "black", size = 0.3) +
  cowplot::draw_line(x = c(0.656, 0.902), y = c(0.544, 0.544), color = "black", size = 0.3)

final_plot










extra_elements <- ggplot() +
  geom_sf(data = r3c, fill = "orange", alpha = 0.3) +
annotation_scale(location = "bl", width_hint = 0.2) +  # Scale bar
  annotation_north_arrow(location = "b", which_north = "true", 
                         style = north_arrow_fancy_orienteering()) +  # North arrow
  theme_void()

cowplot::plot_grid(
  main_map,  # Main map on the left (70%)
  plot_grid(inset_map, extra_elements, ncol = 1, rel_heights = c(0.7, 0.3)),  # Inset on top-right, extras below
  ncol = 2, rel_widths = c(0.7, 0.3)  # Main = 70%, Inset + Extras = 30%
)
legend <- get_legend(main_map)

plot_grid(
  main_map,  
  plot_grid(inset_map, legend, extra_elements, ncol = 1, rel_heights = c(0.6, 0.2, 0.2)),  
  ncol = 2, rel_widths = c(0.7, 0.3)  # Left = Main map (70%), Right = Inset + Legend + Scale (30%)
)

ggdraw(main_map, xlim = c(0.5,1), clip = "off")
