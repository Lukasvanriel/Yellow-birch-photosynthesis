# Date: 2023-09-25
# Author: Lukas Van Riel

### Restore the correct R environment and package versions -----------------------
require(renv)

renv::status() # See which packages need installing
renv::restore() # Disclaimer: The renv.lock was created on Mac, other OS might alter package versions

### Load packages ----------------------------------------------------------------
library(tidyverse) 
library(plantecophys) # For the fitting of the ACi curves and parameter extraction
library(here) # For convenient pathing


### Create raw data --------------------------------------------------------------
# For this project I will use data similar to the manyacidat dataset of the plantecophys package
raw.data <- plantecophys::manyacidat

# For the sake of this mock project, I will randomly identify trees as adult/juvenile
raw.data <- raw.data %>% 
  mutate(adult = ifelse(treatment %in% c(1000, 5, 10, 15), 1, 0)) %>% 
  select(-treatment)

# Write out the resulting 'raw' data file
write.csv(raw.data, file=here("00_rawdata", "VanRiel_LI-COR6400_20230920.csv"),
          row.names = FALSE)


### Load data ----------------------------------------------------------------
# Get current working directory:
here() # Make sure this is the root folder of the project

data <- read.csv(file = here("00_rawdata", "VanRiel_LI-COR6400_20230920.csv"))

# Extract one single individual to illustrate the fitting process later
data.single <- data %>% 
  filter(Curve==unique(Curve)[1])

### Analysis ---------------------------------------------------------------------

# Run analysis for one single tree: 
fit.single <- fitaci(data=data.single)
# The fitaci function fits the ACi curve and returns parameter values
      
# Plot the ACi curve (black dots) together with the fit of the fitaci function 

plot.single <- plot(fit.single)

# Now fit the curves for all trees at the same time:

fits <- fitacis(data, group="Curve", fitmethod="bilinear")
# From the plantecophys package documentation: The bilinear method is much faster
# than the default method.

# Visually inspect the relationship between Vcmax and Jmax. ideally this is linear
with(coef(fits), plot(Vcmax, Jmax))


# combine all resulting coeficients into one dataframe and select the relevant coefficients
coef.mult <- coef(fits) %>% 
  select(Curve, Vcmax, Jmax, Rd)

# Add the life stage again and convert the adult column into a factor for the boxplot and ANOVA.
# First only select the adult column since that is the only one we want to add. 
# Then filter out all the duplicate rows that occur since we dropped the other columns.
adult.df <- data %>% 
  select(Curve, adult) %>% 
  distinct(Curve, .keep_all = TRUE)
  
# Now join the two dataframes
coef.comb <- coef.mult %>% 
  left_join(., adult.df) %>% 
  mutate(adult=as.factor(adult))


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

# Carry out an Anova for both Jmax and Vcmax
aov.Vcmax <- aov(Vcmax ~ adult, data = coef.comb)
summary(aov.Vcmax)

aov.Jmax <- aov(Jmax ~ adult, data = coef.comb)
summary(aov.Jmax)

# Write out the data that was used for the analysis and figures:
write.csv(coef.comb, file=here("02_outdata", "VanRiel_20230926_coefficients.csv"),
          row.names = FALSE)
