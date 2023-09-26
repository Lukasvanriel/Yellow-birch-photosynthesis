# Date: 2023-09-25
# Author: Lukas Van Riel

### Load packages ----
library(tidyverse)
library(plantecophys)
library(here)
library(rticles)

# Get current working directory:
here()

### Load data ----
data <- read.csv(file = here("00_rawdata", "LI-COR6400_leaf_measurements.csv"))

data.single <- data %>% 
  filter(Curve==unique(Curve)[1])

### Analysis ----

# Run analysis for one single tree: 
plot.single <- fitaci(data=data.single)

plot(plot.single)

# Now run for all trees at the same time:
fits <- fitacis(data, "Curve", fitmethod="bilinear") 
# From the plantecophys package: The bilinear method is much faster than the default method.

with(coef(fits), plot(Vcmax, Jmax))

# combine all resulting coeficients into one dataframe
coef.mult <- coef(fits)

# Add the life stage again and convert the adult column into a factor for the plot and ANOVA
coef.comb <- coef.mult %>% 
  left_join(., select(data, Curve, adult)) %>% 
  mutate(adult=as.factor(adult))

with(coef.comb, plot(Vcmax, Jmax))

# Plot the values per life stage

# Plot the Vcmax
ggplot(data=coef.comb, aes(x = adult, y = Vcmax)) +
  geom_boxplot() +
  ggtitle("Vcmax distributions of juveniles vs. adults") +
  labs(y=bquote('Vcmax (\u03BCmol'~m^-2 ~s^-1 ))

    
# Plot the Jmax
ggplot(data=coef.comb, aes(x = adult, y = Jmax)) +
  geom_boxplot() +
  ggtitle("Jmax distributions of juveniles vs. adults") +
  labs(y=bquote('Jmax (\u03BCmol'~m^-2 ~s^-1 ))

# Anova for each
aov.Vcmax <- aov(Vcmax ~ adult, data = coef.comb)
summary(aov.Vcmax)

aov.Jmax <- aov(Jmax ~ adult, data = coef.comb)
summary(aov.Jmax)

