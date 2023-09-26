# Date: 2023-09-25
# Author: Lukas Van Riel

### Load packages ----
library(tidyverse)
library(plantecophys)
library(here)
library(rticles)

### Load data ----
data <- read.csv(file = here("00_rawdata", "LI-COR6400_leaf_measurements.csv"))

data1 <- data %>% 
  filter(Curve==unique(Curve)[4])
  
### Analysis ----

p <- fitaci(data1)

plot(p)

