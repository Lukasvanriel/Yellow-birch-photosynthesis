## File to conduct photosynthetis analysis
## Lukas Van Riel       13/11/2025
##### Load packages #####
library(tidyverse)
library(conflicted)
library(vegan)
#library(car)
#library(MVN)
#library(sf)

conflicts_prefer(dplyr::filter)

##### Load data #####
data.all <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/Combined.csv")
elevation <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Elevation/el.csv")
soil <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP-SoilLayers.csv") 
plot <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP-PlotInfo.csv") 
comm.tree <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP_TreeCommunity.csv")
comm.under <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP-UnderstoryCommunity.csv")

### Wrangle:
data.all <- data.all %>% 
  mutate(
    Bin = factor(format(round(as.numeric(substring(data.all$Bin, 2, 4)) / 10, 1), nsmall = 1))
  ) %>% 
  rename(Tree = tree) %>% 
  filter(Tree != "P475.1-e") %>%  # The extra measurements
  filter(LeafArea != 40.5634)

data.elevation <- elevation %>% 
  mutate(Tree = data.all$Tree) %>% 
  rename(Elevation = Value) %>% 
  select(-ID) %>% 
  relocate(Tree, .before = Elevation)

data.soil <- soil %>% 
  filter(Horizon == "B1") %>% 
  select(Plot, Horizon, soil_class, pH) %>% 
  rbind(c("P470.3", "B1", "LoSa", 4.91)) %>% 
  mutate(pH = as.numeric(pH))

# Remove some extra measuremnets (for now?)
data.all <- data.all %>% 
  filter(Tree != "P475.1-e") %>%  # The extra measurements
  filter(LeafArea != 40.5634) # Double measurements for P480.1-9

data.cca <- data.all %>% 
  left_join(data.elevation, by = "Tree") %>% 
  select(Tree, Plot, Bin,  Stage, V_cmax, J_max, diameter, LeafArea, Longitude, Latitude, age, Elevation)  %>% 
  left_join(data.soil, by = "Plot") %>% 
  mutate(Stage = factor(Stage),
         soil_class = factor(soil_class)) %>% 
  select(-Horizon) %>% 
  left_join(plot %>% select(Plot, Exposure, Slope), by = "Plot") %>% 
  mutate(Exposure = factor(Exposure),
         Slope = factor(Slope)) 


comm.tree.summ <- comm.tree %>% 
  select(-Plot) %>% 
  decostand(method = "hellinger") %>% 
  bind_cols(select(comm.tree, Plot)) %>% 
  relocate(Plot, .before = everything()) %>% 
  mutate(
    Bin = factor(format(round(as.numeric(substring(comm.tree$Plot, 2, 4)) / 10, 1), nsmall = 1))
  )

comm.under.summ <- comm.under %>% 
  select(-Plot) %>% 
  decostand(method = "hellinger") %>% 
  bind_cols(select(comm.under, Plot)) %>% 
  mutate(
    Bin = factor(format(round(as.numeric(substring(comm.tree$Plot, 2, 4)) / 10, 1), nsmall = 1))
  )


##### Functions #####


##### Body #####

#### RDA exploratory analysis ####

# Tree community PCA (exclude Plot and Bin columns)
tree_pca <- rda(comm.tree.summ[, -c(1, 14)])  # Exclude Plot and Bin

# Check variance explained
summary(tree_pca)$cont$importance[2, 1:2]  # % variance by PC1 and PC2

# Extract scores
tree_scores <- scores(tree_pca, display = "sites", choices = 1:2)
tree_axis_data <- data.frame(
  Plot = comm.tree.summ$Plot,
  Tree_PC1 = tree_scores[, 1],
  Tree_PC2 = tree_scores[, 2]
)

# Understory community PCA
under_pca <- rda(comm.under.summ[, -c(58, 59)])  # Exclude Plot and Bin

summary(under_pca)$cont$importance[2, 1:2]

under_scores <- scores(under_pca, display = "sites", choices = 1:2)
under_axis_data <- data.frame(
  Plot = comm.under.summ$Plot,
  Under_PC1 = under_scores[, 1],
  Under_PC2 = under_scores[, 2]
)

# Merge with your main dataset
data.all <- data.cca %>%
  left_join(tree_axis_data, by = "Plot") %>%
  left_join(under_axis_data, by = "Plot")

# ---- Tree Community PCA Plot ----
tree_sites <- as.data.frame(scores(tree_pca, display = "sites", choices = 1:2))
tree_species <- as.data.frame(scores(tree_pca, display = "species", choices = 1:2))

tree_sites$Plot <- comm.tree.summ$Plot
tree_sites$Bin <- as.numeric(as.character(comm.tree.summ$Bin))  # Convert to numeric
tree_species$Species <- rownames(tree_species)

tree_var <- summary(tree_pca)$cont$importance[2, 1:2] * 100

ggplot() +
  geom_point(data = tree_sites, aes(x = PC1, y = PC2, color = Bin), size = 3) +
  geom_text(data = tree_sites, aes(x = PC1, y = PC2, label = Plot), vjust = -1, size = 3) +
  geom_segment(data = tree_species, aes(x = 0, y = 0, xend = PC1*1.5, yend = PC2*1.5), 
               arrow = arrow(length = unit(0.2, "cm")), color = "gray30", alpha = 0.6) +
  geom_text(data = tree_species, aes(x = PC1*1.5, y = PC2*1.5, label = Species), 
            color = "gray30", size = 3, fontface = "bold") +
  scale_color_gradient(low = "red", high = "blue", name = "Latitude") +
  labs(x = paste0("PC1 [", round(tree_var[1], 1), "%]"),
       y = paste0("PC2 [", round(tree_var[2], 1), "%]"),
       title = "Tree Community Composition (Hellinger PCA)") +
  theme_minimal() +
  theme(legend.position = "right")

# ---- Understory Community PCA Plot ----
under_sites <- as.data.frame(scores(under_pca, display = "sites", choices = 1:2))
under_species <- as.data.frame(scores(under_pca, display = "species", choices = 1:2))

under_sites$Plot <- comm.under.summ$Plot
under_sites$Bin <- as.numeric(as.character(comm.under.summ$Bin))
under_species$Species <- rownames(under_species)

under_var <- summary(under_pca)$cont$importance[2, 1:2] * 100

ggplot() +
  geom_point(data = under_sites, aes(x = PC1, y = PC2, color = Bin), size = 3) +
  geom_text(data = under_sites, aes(x = PC1, y = PC2, label = Plot), vjust = -1, size = 2.5) +
  geom_segment(data = under_species, aes(x = 0, y = 0, xend = PC1*2, yend = PC2*2), 
               arrow = arrow(length = unit(0.15, "cm")), color = "gray30", alpha = 0.4) +
  #geom_text(data = under_species, aes(x = PC1*2, y = PC2*2, label = Species), 
  #          color = "gray30", size = 2, alpha = 0.7) +
  scale_color_gradient(low = "red", high = "blue", name = "Latitude") +
  labs(x = paste0("PC1 [", round(under_var[1], 1), "%]"),
       y = paste0("PC2 [", round(under_var[2], 1), "%]"),
       title = "Understory Community Composition (Hellinger PCA)") +
  theme_minimal() +
  theme(legend.position = "right")


######RDA:

# Prepare data - complete cases only
rda_data <- data.all %>%
  select(V_cmax, J_max, Latitude, Elevation, Slope, Exposure, 
         LeafArea, diameter, Stage, Tree_PC1, Tree_PC2, Under_PC1, Under_PC2) %>%
  drop_na()

# Response matrix
photo_matrix <- rda_data[, c("V_cmax", "J_max")]

# RDA model
rda_model <- rda(photo_matrix ~ Latitude + Elevation + Slope + Exposure + 
                   LeafArea + diameter + Stage + 
                   Tree_PC1 + Tree_PC2 + Under_PC1 + Under_PC2,
                 data = rda_data)

# Summary
summary(rda_model)

# Variance explained
RsquareAdj(rda_model)

# Significance tests
anova(rda_model, permutations = 999)  # Overall
anova(rda_model, by = "terms", permutations = 999)  # By variable
anova(rda_model, by = "axis", permutations = 999)  # By RDA axis

# Biplot
plot(rda_model, scaling = 2)



#######LMM
library(nlme)

# Simple correlation
cor.test(data.all$V_cmax, data.all$J_max)

# Visualize
plot(data.all$V_cmax, data.all$J_max, 
     xlab = "Vcmax", ylab = "Jmax",
     main = "Correlation between Vcmax and Jmax")
abline(lm(J_max ~ V_cmax, data = data.all), col = "red")

# With ggplot
library(ggplot2)
ggplot(data.all, aes(V_cmax, J_max)) +
  geom_point(aes(color = Stage), size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Vcmax (µmol/m²/s)", 
       y = "Jmax (µmol/m²/s)",
       title = "Vcmax-Jmax Coordination") +
  theme_minimal()

# By group if you want
data.all %>%
  group_by(Stage) %>%
  summarize(Correlation = cor(V_cmax, J_max),
            N = n())

## OK so absolutely necessary to use multivariate technique

# Create long format
data_long <- data.all %>%
  select(Tree, Plot, Bin, Stage, Latitude, Longitude, Elevation, Slope, Exposure,
         LeafArea, diameter, Tree_PC1, Tree_PC2, Under_PC1, Under_PC2,
         V_cmax, J_max) %>%
  pivot_longer(cols = c(V_cmax, J_max),
               names_to = "Parameter",
               values_to = "Value") %>%
  drop_na()

# Make sure factors are factors
data_long$Parameter <- factor(data_long$Parameter)
data_long$Stage <- factor(data_long$Stage)
data_long$Exposure <- factor(data_long$Exposure)
data_long$Bin <- factor(data_long$Bin)
data_long$Plot <- factor(data_long$Plot)
data_long$Tree <- factor(data_long$Tree)


# Start with main predictors from RDA

# Step 1: Simplify random effects - just random intercepts
mv_model <- lme(
  Value ~ Parameter * (Latitude + Slope + Exposure + Stage + diameter),
  random = ~ 1 | Bin/Plot,  # Random intercepts only
  weights = varIdent(form = ~1|Parameter),
  data = data_long,
  method = "REML"
)

summary(mv_model)
anova(mv_model)




mv_model_simple <- lme(
  Value ~ Parameter * (Latitude + Slope + Exposure + Stage + diameter),
  random = ~1 | Bin,
  weights = varIdent(form = ~1|Parameter),
  data = data_long,
  method = "REML"
)

summary(mv_model_simple)
anova(mv_model_simple)



library(lme4)

mv_model_fixed <- lmer(
  Value ~ Parameter * (Latitude + Slope + Exposure + Stage + diameter + Bin) +
    (1|Plot),
  data = data_long
)

summary(mv_model_fixed)
anova(mv_model_fixed)

# 1. Basic residual plots
plot(mv_model_simple)  # Residuals vs fitted
qqnorm(residuals(mv_model_simple, type = "normalized"))
qqline(residuals(mv_model_simple, type = "normalized"))

# 2. Check residuals by groups
boxplot(residuals(mv_model_simple, type = "normalized") ~ data_long$Bin)
boxplot(residuals(mv_model_simple, type = "normalized") ~ data_long$Parameter)

# 3. Check for spatial autocorrelation
plot(ACF(mv_model_simple, resType = "normalized"), alpha = 0.05)

# 4. Variance homogeneity check (are the different variances for Vcmax/Jmax working?)
plot(mv_model_simple, resid(., type = "normalized") ~ fitted(.) | Parameter)

# 5. Check influential observations
plot(mv_model_simple, id = 0.05)  # Labels influential points




#### Spatial autocorrelation is a problem!

# Compound symmetry within Bin
mv_model_spatial_simple <- lme(
  Value ~ Parameter * (Latitude + Slope + Exposure + Stage + diameter),
  random = ~1 | Bin,
  correlation = corCompSymm(form = ~ 1 | Bin),
  weights = varIdent(form = ~1|Parameter),
  data = data_long,
  method = "REML"
)

summary(mv_model_spatial_simple)
anova(mv_model_spatial_simple)

# Check if ACF is now clean
plot(ACF(mv_model_spatial_simple, resType = "normalized"), alpha = 0.05)

# Compare models
anova(mv_model_simple, mv_model_spatial_simple)

# Get unique plot coordinates
plot_coords <- data_long %>%
  group_by(Plot, Bin) %>%
  summarize(Latitude = unique(Latitude),
            Longitude = unique(Longitude),
            .groups = "drop")

# Check - should have 15 rows
nrow(plot_coords)  # Should be 15

# Merge back to ensure coordinates are consistent
data_long <- data_long %>%
  select(-Latitude, -Longitude) %>%  # Remove old ones
  left_join(plot_coords, by = c("Plot", "Bin"))

# Now fit with 2D spatial correlation at PLOT level
mv_spatial_2D <- lme(
  Value ~ Parameter * (Latitude + Slope + Exposure + Stage + diameter),
  random = ~1 | Bin/Plot,  # Plot nested in Bin
  correlation = corExp(form = ~ Latitude + Longitude | Bin/Plot),  # Spatial correlation between plots
  weights = varIdent(form = ~1|Parameter),
  data = data_long,
  method = "REML"
)

summary(mv_spatial_2D)
anova(mv_spatial_2D)


# Option 1: Spatial correlation across all observations
mv_spatial_unnest <- lme(
  Value ~ Parameter * (Latitude + Slope + Exposure + Stage + diameter),
  random = ~1 | Bin/Plot,
  correlation = corExp(form = ~ Latitude + Longitude),  # No grouping factor
  weights = varIdent(form = ~1|Parameter),
  data = data_long,
  method = "REML"
)

set.seed(123)
data_long <- data_long %>%
  mutate(Lat_jitter = Latitude + rnorm(n(), 0, 0.00001),
         Lon_jitter = Longitude + rnorm(n(), 0, 0.00001))

# Now try with jittered coordinates
mv_spatial_jitter <- lme(
  Value ~ Parameter * (Latitude + Slope + Exposure + Stage + diameter),
  random = ~1 | Bin/Plot,
  correlation = corExp(form = ~ Lat_jitter + Lon_jitter),
  weights = varIdent(form = ~1|Parameter),
  data = data_long,
  method = "REML"
)

summary(mv_spatial_jitter)
anova(mv_spatial_jitter)

# Check ACF
plot(ACF(mv_spatial_jitter, resType = "normalized"), alpha = 0.05)

# Compare models
anova(mv_model_simple, mv_spatial_jitter)


###### 
# 1. Look at residuals by group - are outliers driving results?
boxplot(residuals(mv_model_simple, type = "normalized") ~ data_long$Exposure)
boxplot(residuals(mv_model_simple, type = "normalized") ~ data_long$Slope)

# 2. Check group means - are single-plot groups extreme?
data_long %>%
  group_by(Exposure, Parameter) %>%
  summarize(Mean = mean(Value), N_plots = n_distinct(Plot))

data_long %>%
  group_by(Slope, Parameter) %>%
  summarize(Mean = mean(Value), N_plots = n_distinct(Plot))

# 3. Sensitivity: Drop the single-plot exposures
mv_model_balanced <- lme(
  Value ~ Parameter * (Latitude + Slope + Stage + diameter),  # Drop Exposure
  random = ~1 | Bin,
  weights = varIdent(form = ~1|Parameter),
  data = data_long,
  method = "REML"
)

anova(mv_model_balanced)  # Is Slope still significant?




### SLope effect seems robust:
library(emmeans)
emmeans(mv_model_simple, pairwise ~ Slope, at = list(Parameter = "V_cmax"))
emmeans(mv_model_simple, pairwise ~ Exposure, at = list(Parameter = "V_cmax"))

# Extract predictions
library(effects)
plot(Effect(c("Slope", "Parameter"), mv_model_simple))





# Calculate ratio
data.all$JV_ratio <- data.all$J_max / data.all$V_cmax

# Check distribution
hist(data.all$JV_ratio)
summary(data.all$JV_ratio)

# Model it
library(lme4)
ratio_model <- lmer(JV_ratio ~ Latitude + Slope + Exposure + Stage + diameter +
                      (1|Bin) + (1|Bin:Plot),
                    data = data.all)

summary(ratio_model)
anova(ratio_model)

# Visualize
boxplot(JV_ratio ~ Slope, data = data.all)
boxplot(JV_ratio ~ Stage, data = data.all)



library(emmeans)
emmeans(ratio_model, pairwise ~ Slope)




# 1. Full model with community
mv_model_community <- lme(
  Value ~ Parameter * (Latitude + Slope + Exposure + Stage + diameter + 
                         Tree_PC1 + Tree_PC2 + Under_PC1 + Under_PC2),
  random = ~1 | Bin,
  weights = varIdent(form = ~1|Parameter),
  data = data_long,
  method = "REML"
)

anova(mv_model_community)

# 2. Compare to model without community
anova(mv_model_simple, mv_model_community)

# 3. Check if community correlates with topography
cor.test(data.all$Tree_PC1, as.numeric(data.all$Slope))
cor.test(data.all$Under_PC1, as.numeric(data.all$Slope))

# Also visual
boxplot(Tree_PC1 ~ Slope, data = data.all)
boxplot(Under_PC1 ~ Slope, data = data.all)

# 4. Ratio model with community
ratio_community <- lmer(JV_ratio ~ Slope + Exposure + Stage + 
                          Tree_PC1 + Tree_PC2 + Under_PC1 + Under_PC2 +
                          (1|Bin) + (1|Bin:Plot),
                        data = data.all)
anova(ratio_community)




# Refit both with ML (not REML)
mv_simple_ML <- update(mv_model_simple, method = "ML")
mv_community_ML <- update(mv_model_community, method = "ML")

anova(mv_simple_ML, mv_community_ML)

# Check if community correlates with latitude
cor.test(data.all$Tree_PC1, data.all$Latitude)
cor.test(data.all$Under_PC1, data.all$Latitude)




####USing Bin as 
# Option 1: Drop random Bin, keep Plot random
mv_model_bin <- lme(
  Value ~ Parameter * (Bin + Slope + Exposure + Stage + diameter),
  random = ~1 | Plot,
  weights = varIdent(form = ~1|Parameter),
  data = data_long,
  method = "REML"
)

anova(mv_model_bin)

# Post-hoc comparisons between Bins for Vcmax
library(emmeans)
emmeans(mv_model_bin, pairwise ~ Bin, at = list(Parameter = "V_cmax"))

# Visualize pattern across bins
emmeans_result <- emmeans(mv_model_bin, ~ Bin | Parameter)
plot(emmeans_result)

# Or simpler visualization
data_long %>%
  group_by(Bin, Parameter) %>%
  summarize(Mean = mean(Value), SE = sd(Value)/sqrt(n())) %>%
  ggplot(aes(Bin, Mean, color = Parameter, group = Parameter)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = Mean-SE, ymax = Mean+SE), width = 0.2) +
  theme_minimal()


# Check pH range and distribution
summary(data.cca$pH)
hist(data.cca$pH)

# Option 1: Recreate data_long from data.cca (cleanest)
data_long <- data.cca %>%
  select(Tree, Plot, Bin, Stage, Latitude, Longitude, Elevation, Slope, Exposure,
         LeafArea, diameter, age, pH, soil_class,  # Added age, pH, soil_class
         V_cmax, J_max) %>%
  pivot_longer(cols = c(V_cmax, J_max),
               names_to = "Parameter",
               values_to = "Value") %>%
  drop_na(Value) %>%  # Keep rows with valid Vcmax/Jmax
  mutate(
    Parameter = factor(Parameter),
    Stage = factor(Stage),
    Exposure = factor(Exposure),
    Bin = factor(Bin),
    Plot = factor(Plot),
    Tree = factor(Tree),
    Slope = factor(Slope),
    soil_class = factor(soil_class)
  )

# Now test
mv_model_pH <- lme(
  Value ~ Parameter * (Slope + pH),
  random = ~1 | Bin,
  weights = varIdent(form = ~1|Parameter),
  data = data_long
)
anova(mv_model_pH)

# Test 1: Main effect
mv_model_pH <- lme(
  Value ~ Parameter * (Slope + pH),
  random = ~1 | Bin,
  weights = varIdent(form = ~1|Parameter),
  data = data_long
)

# Test 2: Does pH modify slope effect?
mv_model_pH_int <- lme(
  Value ~ Parameter * (Slope * pH),  # Interaction
  random = ~1 | Bin,
  weights = varIdent(form = ~1|Parameter),
  data = data_long
)

anova(mv_model_pH_int)

mv_model_full <- lme(
  Value ~ Parameter * (Slope + pH + age),
  random = ~1 | Bin,
  weights = varIdent(form = ~1|Parameter),
  data = data_long %>% filter(Stage == "adult")  # age only for adults
)

anova(mv_model_full)

## Decline per decade?
# Quick check
data.all %>% 
  filter(Stage == "adult") %>%
  summarize(cor_age_vcmax = cor(age, V_cmax, use="complete.obs"))

data.all %>% filter(Stage == "adult") %>%
  group_by(Slope) %>%
  summarize(
    mean_age = mean(age, na.rm=T),
    sd_age = sd(age, na.rm=T),
    n = n()
  )

# Visual check
boxplot(age ~ Slope, data = data.all %>% filter(Stage == "adult"))

# Stats
anova(lm(age ~ Slope, data = data.all %>% filter(Stage == "adult")))


##### Soil texture
# Check texture distribution by slope
table(data.all$Slope, data.all$soil_class)

mv_model_texture <- lme(
  Value ~ Parameter * (Slope + soil_class + age),
  random = ~1 | Bin,
  weights = varIdent(form = ~1|Parameter),
  data = data_long %>% filter(Stage == "adult")
)
anova(mv_model_texture)


# 1. Post-hoc: Which textures differ?
library(emmeans)
emmeans(mv_model_texture, pairwise ~ soil_class, at = list(Parameter = "V_cmax"))

# 2. Test slope×texture interaction
mv_model_texture_int <- lme(
  Value ~ Parameter * (Slope * soil_class + age),
  random = ~1 | Bin,
  weights = varIdent(form = ~1|Parameter),
  data = data_long %>% filter(Stage == "adult")
)
anova(mv_model_texture_int)

# 3. Visual check - does texture matter WITHIN slope B (the only testable slope)?
data.all %>% 
  filter(Stage == "adult", Slope == "B") %>%
  group_by(soil_class) %>%
  summarize(
    mean_Vcmax = mean(V_cmax, na.rm=T),
    n = n()
  )

boxplot(V_cmax ~ soil_class, 
        data = data.all %>% filter(Stage == "adult", Slope == "B"),
        main = "Vcmax by texture (Slope B only)")


#### Is there confounding:
# Slope B is the only one with texture variation
data.all %>% 
  filter(Stage == "adult", Slope == "B") %>%
  group_by(soil_class) %>%
  summarize(
    mean_Vcmax = mean(V_cmax, na.rm=T),
    sd_Vcmax = sd(V_cmax, na.rm=T),
    n = n()
  )

# Stats within Slope B only
slope_b_test <- lm(V_cmax ~ soil_class + age, 
                   data = data.all %>% filter(Stage == "adult", Slope == "B"))
summary(slope_b_test)
anova(slope_b_test)



##### Final model

# Your robust model
mv_model_final <- lme(
  Value ~ Parameter * (Slope + Latitude + age),
  random = ~1 | Bin,
  weights = varIdent(form = ~1|Parameter),
  data = data_long %>% filter(Stage == "adult")
)

