## File to conduct first exploratory analysis
## Lukas Van Riel       17/12/2024
##### Load packages #####
library(tidyverse)
library(conflicted)
library(vegan)
library(car)
library(MVN)

conflicts_prefer(dplyr::filter)

##### Load data #####
data.all <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/Combined.csv")
comm.tree <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP_TreeCommunity.csv")
comm.under <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP-UnderstoryCommunity.csv")

### Wrangle:
data.all <- data.all %>% 
  mutate(
    Bin = factor(format(round(as.numeric(substring(data.all$Bin, 2, 4)) / 10, 1), nsmall = 1))
  )

# Remove some extra measuremnets (for now?)
data.all <- data.all %>% 
  filter(tree != "P475.1-e") %>%  # The extra measurements
  filter(LeafArea != 40.5634) # Double measurements for P480.1-9

comm.tree.summ <- comm.tree %>% 
  select(-Plot) %>% 
  decostand(method = "hellinger") %>% 
  bind_cols(select(comm.tree, Plot)) %>% 
  relocate(Plot, .before = everything()) %>% 
  mutate(
    Bin = factor(format(round(as.numeric(substring(comm.tree$Plot, 2, 4)) / 10, 1), nsmall = 1))
  )

comm.tree.pa <- comm.tree %>% 
  select(-Plot) %>% 
  decostand(method = "pa") %>% 
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

### PERMANOVA ####
## Assumptions:
#Multivariate normality test
mvn(data = data.all[, c("J_max", "V_cmax")], mvnTest = "mardia") # signficant

## Adults/Saplings for all bins
data.adult <- data.all %>% 
  filter(Stage == "adult")
data.sapling <- data.all %>% 
  filter(Stage == "sapling")

dist.matrix.adult <- dist(data.adult[, c("V_cmax", "J_max")], method = "euclidean")
dist.matrix.sapling <- dist(data.sapling[, c("V_cmax", "J_max")], method = "euclidean")

permanova.adult <- adonis2(dist.matrix.adult ~ Bin, data = data.adult)
permanova.sapling <- adonis2(dist.matrix.sapling ~ Bin, data = data.sapling)

# Summary
permanova.adult # Not significant
permanova.sapling # Not significant


dispersion <- betadisper(dist.matrix.adult, data.adult %>% pull(Bin))
anova(dispersion)

dispersion <- betadisper(dist.matrix.sapling, data.sapling %>% pull(Bin))
anova(dispersion)


## Adults vs Saplings per bin
permanova.bins <- lapply(levels(data.all$Bin), FUN = function(x){
  data.bin <- data.all %>% 
    filter(Bin == x)
  
  dist.matrix.bin <- dist(data.bin[, c("V_cmax", "J_max")], method = "euclidean")
  
  permanova.bin <- adonis2(dist.matrix.bin ~ Stage, data = data.bin)
  
  permanova.bin
} )

permanova.bins[[5]]

####
# Convert data to long format for ggplot
df_long <- data.all %>%
  pivot_longer(cols = c(V_cmax, J_max), names_to = "Variable", values_to = "Value") %>% 
  mutate(Stage = ifelse(Stage=="adult", "Adults", "Saplings"),
         Variable = factor(Variable, levels = c("V_cmax", "J_max"),
                           labels = c(expression(paste(V[paste(c, ",", max, sep = "")], ~(mu*mol / m^2 / s))),
                                      expression(paste(J[max], ~(mu*mol / m^2 / s))) ) ))


# Create boxplot
ggplot(df_long, aes(x = Bin, y = Value, fill = Stage)) +
  geom_boxplot(alpha = 0.9, outlier.shape = NA) +  # Boxplot with transparent fill
  #geom_jitter(width = 0.2, alpha = 0.2, size = 1.5) +  # Adds jitter for visibility
  facet_grid(rows = vars(Variable), cols = vars(Stage), scales = "free_y", switch = "y", labeller = label_parsed) +  # Facet by Type (row) and Variable (column)
  scale_fill_brewer(palette = "Set3") +  # Nice color scheme
  labs(#title = "Distribution of Vcmax and Jmax by stage and Latitude",
       x = "Latitudinal bin (°N)", y = NULL) + #expression(paste(V[paste(c, ",", max, sep = "")], ~(mu*mol / m^2 / s)))) +
  theme_minimal(base_size = 12) +  # Clean theme +
  theme(
    legend.position = "none",  # Hide legend (optional),  # Left-aligned facet labels for "a" and "s"
    strip.placement = "outside",# Move facet labels outside the plot area
    strip.background = element_blank(),
    strip.text.x = element_text(size = 13, angle = 0, hjust = 0.5, face = "bold"),
    strip.text.y.left = element_text(size = 12, angle = 90, hjust = 0.5, face = "bold")
  ) # Hide legend (optional)


### Model relationships ####
model <- lm(V_cmax ~ J_max , data=data.all)
summary(model)


ggplot(data.all, aes(x = J_max, y = V_cmax)) +
  geom_point(alpha = 0.6) +  # Scatterplot points
  geom_smooth(method = "lm", color = "blue", se = TRUE) +  # Linear fit with confidence interval
  labs(title = "Relationship between Vcmax and Jmax",
       x = "Vcmax",
       y = "Jmax") +
  theme_minimal()


model <- lm(V_cmax ~ diameter + LeafArea + Width + Length + age, data=data.all)
summary(model)

p1 <- ggplot(data.all, aes(x = diameter, y = V_cmax)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  labs(title = "", x = "diameter", y = "V_cmax") +
  theme_minimal()

p2 <- ggplot(data.all, aes(x = LeafArea, y = V_cmax)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  labs(title = "", x = "LeafArea", y = "V_cmax") +
  theme_minimal()

p3 <- ggplot(data.all, aes(x = Width, y = V_cmax)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  labs(title = "", x = "Width", y = "V_cmax") +
  theme_minimal()

p4 <- ggplot(data.all, aes(x = age, y = V_cmax)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  labs(title = "", x = "Minimum age", y = "V_cmax") +
  theme_minimal()

model <- lm(V_cmax ~ age, data=data.all)
summary(model)

####

mvn(data = data.all[, c("age")], mvnTest = "mardia") # signficant
hist(data.all$age)

m1 <- aov(age ~ Bin, data = data.all)
summary(m1)

ggplot(data.all) +
  geom_boxplot(aes(x=Bin, y=age, fill = Bin)) +
  scale_fill_manual(values = c("#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac")) +
  labs(
    title = expression(paste("Boxplot of ", V[paste(c, ",", max, sep = "")], " per latitude")),
    x = "Latitude (°N)",
    y = "Minimum age (yr)"
  ) + 
  theme_minimal()

### community compositions ####
comm.tree.rel <- comm.tree %>% 
  select(-Plot) %>% 
  decostand(method = "total") %>% 
  bind_cols(select(comm.tree, Plot)) %>% 
  relocate(Plot, .before = everything())

dissimilarity_matrix <- vegdist(comm.tree.rel %>% select(-Plot), method = "bray")
dissimilarity_matrix
adonis2(dissimilarity_matrix ~ Bin, data = data.frame(Bin=rep(seq(46, 48, by=0.5), each=3)), permutations = 999)


#### effect of species presence/abundance on vcmax/jmax##
data.comb <- data.all %>% 
  left_join(comm.tree.rel, by="Plot")
data.comb.pa <- data.all %>% 
  left_join(comm.tree.pa, by="Plot")

cover.V <- lm(V_cmax ~ AS + BA + FG + PG + TA + TC + AB + OV + BP + TO + AR + FN, data = data.comb)
summary(cover.V)

cover.test <- lm(J_max ~ AB, data = data.comb)
summary(cover.test)


ggplot(data.comb.pa, aes(x = AR, y = J_max)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  labs(title = "", x = "Relative cover by Acer Sacch.", y = "V_cmax") +
  theme_minimal()


heatmap(round(cor(data.comb[,96:107]),2))

plot(data.comb$PG, data.comb$TA)


#### SOIL
## pH: Linear mixed model:
data.soil <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP-SoilLayers.csv")

data.soil <- data.soil %>%
  mutate(
    Bin = factor(format(round(as.numeric(substring(Plot, 2, 4)) / 10, 1), nsmall = 1))
  )
  

lmer_pH <- lmer(pH ~ Bin + (1 | Plot), data = data.soil)
summary(lmer_pH)

library(lmerTest)
anova(lmer_pH)
hist(data.soil$pH)

library(emmeans)
emm <- emmeans(lmer_pH, pairwise ~ Bin)

plot(emm, comparisons = TRUE)

## Soil texture class:
# By bin
chisq.test(table(data.soil$Bin, data.soil$soil_class))

table_soil <- table(data.soil$Bin, data.soil$soil_class)
prop.table(table_soil, margin = 1)  #

chisq <- chisq.test(table_soil)
chisq$residuals

library(corrplot)
corrplot(chisq$residuals, is.cor = FALSE)

# By plot
chisq.test(table(data.soil$Plot, data.soil$soil_class))

table_soil <- table(data.soil$Plot, data.soil$soil_class)
prop.table(table_soil, margin = 1)  #

chisq <- chisq.test(table_soil)
chisq$residuals

corrplot(chisq$residuals, is.cor = FALSE)

## Permanova for community compositions

dist.matrix.tree <- vegdist(comm.tree.summ%>% select(-c(Plot, Bin)), method = "bray")


adonis2(dist.matrix.tree ~ Bin, data = comm.tree.summ)

betadisper_result <- betadisper(dist.matrix.tree, comm.tree.summ %>% pull(Bin))
anova(betadisper_result)

dist.matrix.under <- vegdist(comm.under.summ%>% select(-c(Plot, Bin)), method = "bray")

adonis2(dist.matrix.under ~ Bin, data = comm.under.summ)

betadisper_result <- betadisper(dist.matrix.under, comm.under.summ %>% pull(Bin))
anova(betadisper_result)


##### Dendro ####
library(treeclim)

plot.data <- read.csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP-PlotInfo.csv")
coords <- plot.data %>% 
  select(Plot, Longitude, Latitude) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

data.rings.480.1 <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/TreeCores-Detrended/P480_1-Detrended-ModNegExp.csv")
clim.TN <- read_csv("/Users/lukas/Downloads/meanT/tg_mean.csv")
clim.PN <- read_csv("/Users/lukas/Downloads/Prectot/prcptot.csv")
clim.TS <- read_csv("/Users/lukas/Downloads/Tmean/tg_mean.csv")
clim.PS <- read_csv("/Users/lukas/Downloads/TotalPrec/prcptot.csv")

clim.T <- rbind(clim.TN, clim.TS)
clim.P <- rbind(clim.PN, clim.PS)

## Extract correct climate data:
clim.T.sf <- clim.T %>% 
  mutate(time = substring(time, 1,7)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

times <- unique(clim.T.sf$time)
results <- list()

for (t in times) {
  print(t)
  # Filter raster data for the current year
  climate_year <- clim.T.sf %>% dplyr::filter(time == t)
  
  # Perform a spatial join for the current year
  matched <- st_join(coords, climate_year, join = st_nearest_feature)
  
  # Add the year column explicitly (it will be consistent)
  matched <- matched %>%
    mutate(time = t) %>%
    select(Plot, time, 
           paste0("ssp245_", "tg_mean", "_p50"),
           paste0("ssp126_", "tg_mean", "_p50"),
           paste0("ssp585_", "tg_mean", "_p50"),
    ) %>% 
    separate(time, into = c("year", "month"), sep = "-") %>% 
    mutate(year = as.numeric(year), month = as.numeric(month)) %>% 
    st_drop_geometry()# Keep relevant columns
  
  # Append the results to the list
  results[[as.character(t)]] <- matched
}
results.comb.T <- bind_rows(results, .id = "column_label")
temp.sel <- results.comb.T  %>% 
  rename(temp = ssp245_tg_mean_p50) %>% 
  select(Plot, year, month, temp)
  
## Extract correct climate data:
clim.P.sf <- clim.P %>% 
  mutate(time = substring(time, 1,7)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

times <- unique(clim.P.sf$time)
results <- list()

for (t in times) {
  print(t)
  # Filter raster data for the current year
  climate_year <- clim.P.sf %>% dplyr::filter(time == t)
  
  # Perform a spatial join for the current year
  matched <- st_join(coords, climate_year, join = st_nearest_feature)
  
  # Add the year column explicitly (it will be consistent)
  matched <- matched %>%
    mutate(time = t) %>%
    select(Plot, time, 
           paste0("ssp245_", "prcptot", "_p50"),
           paste0("ssp126_", "prcptot", "_p50"),
           paste0("ssp585_", "prcptot", "_p50"),
    ) %>% 
    separate(time, into = c("year", "month"), sep = "-") %>% 
    mutate(year = as.numeric(year), month = as.numeric(month)) %>% 
    st_drop_geometry()# Keep relevant columns
  
  # Append the results to the list
  results[[as.character(t)]] <- matched
}
results.comb.P <- bind_rows(results, .id = "column_label")
prec.sel <- results.comb.P %>% 
  rename(prec = ssp245_prcptot_p50) %>% 
  select(Plot, year, month, prec)

head(temp.sel)
head(prec.sel)


## Try for one plot
plot_name <- "P460_1"
data.rings <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/TreeCores-Detrended/", plot_name, "-Detrended-ModNegExp.csv"))

plot.name <- "P460.1"
data.climate <- merge(temp.sel %>% filter(Plot == plot.name),
                   prec.sel %>% filter(Plot == plot.name)) %>% 
  select(-Plot) %>% 
  arrange(year, month)



## Now the treerings
tr <- data.rings %>% 
  column_to_rownames("year")

## Try the dcc
# dc_resp <- dcc(muc_spruce, muc_clim)
test <- dcc(dplR::chron(tr), data.48.1)
test <- dcc(tr, data.48.1, dynamic = "moving")

test$coef
plot(test) + scale_color_manual(values = c("steelblue", "tomato3"))











### OLD ####
plot(data.soil$Bin, data.soil$pH)

#### Community PCA's ------
# Create tree community composition per plot
comm.tree.summ <- comm.tree %>% 
  select(-Plot) %>% 
  decostand(method = "hellinger") %>% 
  bind_cols(select(comm.tree, Plot)) %>% 
  relocate(Plot, .before = everything())

comm.tree.pa <- comm.tree %>% 
  select(-Plot) %>% 
  decostand(method = "pa") %>% 
  bind_cols(select(comm.tree, Plot)) %>% 
  relocate(Plot, .before = everything())

# Perform PCA on the Hellinger-transformed data
pca_comm <- prcomp(comm.tree.summ %>% select(-Plot), scale. = TRUE)

# Summarize PCA results
summary(pca_comm)

# Extract scores and loadings
scores <- as.data.frame(pca_comm$x)  # PCA scores
loadings <- as.data.frame(pca_comm$rotation)  # PCA loadings

# Add row names for labeling
scores$Row <- comm.tree.summ$Plot
loadings$Var <- rownames(loadings)

# Plot PCA scores with original axes
ggplot() +
  # PCA scores (rows)
  geom_point(data = scores, aes(x = PC1, y = PC2), size = 2) +
  geom_text(data = scores, aes(x = PC1, y = PC2, label = Row), vjust = -1) +
  
  # Original axes (loadings)
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1*5, yend = PC2*5), 
               arrow = arrow(length = unit(0.2, "cm")), color = "red") +
  geom_text(data = loadings, aes(x = PC1*5, y = PC2*5, label = Var), color = "red", position=position_jitter(width=0.3,height=0.4)) +
  
  # Add labels and theme
  labs(x = paste0("PC1 [", as.character(round(pca_comm$sdev[1]**2/sum(pca_comm$sdev**2) * 100, 1)), "%]"),
       y = paste0("PC2 [", as.character(round(pca_comm$sdev[2]**2/sum(pca_comm$sdev**2) * 100, 1)), "%]"),
       title = "PCA with Hellinger Transformation") +
  theme_minimal()

## Now for the understory:
comm.under.summ <- comm.under %>% 
  select(-Plot) %>% 
  decostand(method = "hellinger") %>% 
  bind_cols(select(comm.under, Plot))

# Perform PCA on the Hellinger-transformed data
pca_under <- prcomp(comm.under.summ %>% select(-Plot), scale. = TRUE)

# Summarize PCA results
summary(pca_under)

# Extract scores and loadings
scores.under <- as.data.frame(pca_under$x)  # PCA scores
loadings.under <- as.data.frame(pca_under$rotation)  # PCA loadings

# Add row names for labeling
scores.under$Row <- comm.under.summ$Plot
loadings.under$Var <- rownames(loadings.under)

# Plot PCA scores with original axes
ggplot() +
  # PCA scores (rows)
  geom_point(data = scores.under, aes(x = PC1, y = PC2), size = 2) +
  geom_text(data = scores.under, aes(x = PC1, y = PC2, label = Row), vjust = -1) +
  
  # Original axes (loadings)
   geom_segment(data = loadings.under, aes(x = 0, y = 0, xend = PC1*20, yend = PC2*20), 
                arrow = arrow(length = unit(0.2, "cm")), color = "red") +
  # geom_text(data = loadings.under, aes(x = PC1*20, y = PC2*20, label = Var), color = "red", vjust = 1.9) +
  
  # Add labels and theme
  labs(x = paste0("PC1 [", as.character(round(pca_under$sdev[1]**2/sum(pca_under$sdev**2) * 100, 1)), "%]"),
       y = paste0("PC2 [", as.character(round(pca_under$sdev[2]**2/sum(pca_under$sdev**2) * 100, 1)), "%]"),
       title = "PCA with Hellinger Transformation") +
  theme_minimal()

#### Indicator species ----

### Check out indicator species framework
library(indicspecies)
comm.tree.summ$group <- rep(1:5, each=3)

indval <- multipatt(comm.tree.summ %>% select(-Plot, -group), comm.tree.summ$group, 
                    control = how(nperm=999))
summary(indval) 

comm.under.summ$group <- rep(1:5, each=3)
indval.under <- multipatt(comm.under.summ %>% select(-Plot, -group), comm.under.summ$group, 
                    control = how(nperm=999))

summary(indval.under) 

####  Run the ANOVA's:-----

## plot: 

#### Trying to make nice figure
ggplot(data.all) +
  geom_boxplot(aes(x=Bin, y=V_cmax, fill = Bin)) +
  scale_fill_manual(values = c("#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8", "#253494")) +
  labs(
    title = expression(paste("Boxplot of ", V[paste(c, ",", max, sep = "")], " per latitude")),
    x = "Latitude (°N)",
    y = expression(paste(V[paste(c, ",", max, sep = "")], ~(mu*mol / m^2 / s)))
  ) + 
  theme_minimal() + guides(fill="none")

ggplot(data.all) +
  geom_boxplot(aes(x=Bin, y=V_cmax, fill = Bin)) +
  scale_fill_manual(values = c("#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac")) +
  labs(
    title = expression(paste(V[paste(c, ",", max, sep = "")], " per latitude")),
    x = "Latitude (°N)",
    y = expression(paste(V[paste(c, ",", max, sep = "")], ~(mu*mol / m^2 / s)))
  ) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + guides(fill=guide_legend(title="Latitude [°N]"))

ggplot(data.all) +
  geom_boxplot(aes(x=Bin, y=V_cmax, fill = Bin)) +
  scale_fill_manual(values = c("#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac")) +
  labs(
    title = expression(paste(V[paste(c, ",", max, sep = "")], " per latitude")),
    x = "Latitude (°N)",
    y = expression(paste(V[paste(c, ",", max, sep = "")], ~(mu*mol / m^2 / s)))
  ) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
        , panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank()) + guides(fill=guide_legend(title="Latitude [°N]"))

ggplot(data.all) +
  geom_boxplot(aes(x=Bin, y=V_cmax, fill = Bin)) +
  scale_fill_manual(values = c("#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac")) +
  labs(
    title = expression(paste(V[paste(c, ",", max, sep = "")], " per latitude")),
    x = "Latitude (°N)",
    y = expression(paste(V[paste(c, ",", max, sep = "")], ~(mu*mol / m^2 / s)))
  ) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        , panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank()) + guides(fill=guide_legend(title="Latitude [°N]"))

ggplot(data.all) +
  geom_boxplot(aes(x=Bin, y=V_cmax, fill = Bin), outlier.shape = NA) +
  scale_fill_manual(values = c("#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac")) +
  labs(
    title = expression(paste(V[paste(c, ",", max, sep = "")], " per latitude")),
    x = "Latitude (°N)",
    y = expression(paste(V[paste(c, ",", max, sep = "")], ~(mu*mol / m^2 / s)))
  ) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        , panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank()) + guides(fill=guide_legend(title="Latitude [°N]"))

letters <- data.frame(x = levels(data.all$Bin),
                      y = data.all %>% group_by(Bin) %>% summarise(m = max(V_cmax)) %>% pull(m),
                      sign = rep("a", 5))

ggplot(data.all) +
  geom_boxplot(aes(x=Bin, y=V_cmax, fill = Bin), outlier.shape = NA) +
  scale_fill_manual(values = c("#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac")) +
  geom_text(letters, mapping = aes(x = x, y = y, label = sign)) +
  labs(
    title = expression(paste(V[paste(c, ",", max, sep = "")], " per latitude")),
    x = "Latitude (°N)",
    y = expression(paste(V[paste(c, ",", max, sep = "")], ~(mu*mol / m^2 / s)))
  ) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        , panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank()) + guides(fill=guide_legend(title="Latitude [°N]"))

ggplot(data.all) +
  geom_boxplot(aes(x=Bin, y=V_cmax, fill = Bin)) +
  scale_fill_manual(values = c("#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac")) +
  geom_text(letters, mapping = aes(x = x, y = y, label = sign), vjust = -0.8) +
  labs(
    title = expression(paste(V[paste(c, ",", max, sep = "")], " per latitude")),
    x = "Latitude (°N)",
    y = expression(paste(V[paste(c, ",", max, sep = "")], ~(mu*mol / m^2 / s)))
  ) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        , panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank()) + guides(fill=guide_legend(title="Latitude [°N]"))


letters <- data.frame(x = levels(data.all$Bin),
                      y = data.all %>% group_by(Bin) %>% summarise(m = min(V_cmax)) %>% pull(m),
                      sign = rep("a", 5))

ggplot(data.all) +
  geom_boxplot(aes(x=Bin, y=V_cmax, fill = Bin)) +
  scale_fill_manual(values = c("#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac")) +
  geom_text(letters, mapping = aes(x = x, y = y, label = sign), vjust = 1.3) +
  labs(
    title = expression(paste(V[paste(c, ",", max, sep = "")], " per latitude")),
    x = "Latitude (°N)",
    y = expression(paste(V[paste(c, ",", max, sep = "")], ~(mu*mol / m^2 / s)))
  ) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        , panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank()) + guides(fill=guide_legend(title="Latitude [°N]"))


letters <- data.frame(x = levels(data.all$Bin),
                      y = rep(4,5),
                      sign = rep("a", 5))

ggplot(data.all) +
  geom_boxplot(aes(x=Bin, y=V_cmax, fill = Bin)) +
  scale_fill_manual(values = c("#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac")) +
  geom_text(letters, mapping = aes(x = x, y = y, label = sign), vjust = 1.3) +
  labs(
    title = expression(paste(V[paste(c, ",", max, sep = "")], " per latitude")),
    x = "Latitude (°N)",
    y = expression(paste(V[paste(c, ",", max, sep = "")], ~(mu*mol / m^2 / s)))
  ) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        , panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank()) + guides(fill=guide_legend(title="Latitude [°N]"))

letters <- data.frame(x = levels(data.all$Bin),
                      y = rep(75, 5),
                      sign = rep("a", 5))

ggplot(data.all) +
  geom_boxplot(aes(x=Bin, y=V_cmax, fill = Bin)) +
  scale_fill_manual(values = c("#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac")) +
  geom_text(letters, mapping = aes(x = x, y = y, label = sign), vjust = 1.3) +
  labs(
    title = expression(paste(V[paste(c, ",", max, sep = "")], " per latitude")),
    x = "Latitude (°N)",
    y = expression(paste(V[paste(c, ",", max, sep = "")], ~(mu*mol / m^2 / s)))
  ) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        , panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank()) + guides(fill=guide_legend(title="Latitude [°N]"))

#####

ggplot(data.all %>% filter(Stage=="adult")) +
  geom_boxplot(aes(y=V_cmax, x=Bin)) +
  labs(
    title = paste("Boxplot of adult Vcmax per latitude"),
    x = "Latitude [°N]",
    y = "Vc,max [umol / m2 / s]"
  ) + 
  theme_minimal() 

ggplot(data.all %>% filter(Stage=="sapling")) +
  geom_boxplot(aes(y=V_cmax, x=Bin)) +
  labs(
    title = paste("Boxplot of sapling Vcmax per latitude"),
    x = "Latitude [°N]",
    y = "Vc,max [umol / m2 / s]"
  ) + 
  theme_minimal() 

ggplot(data.all) +
  geom_boxplot(aes(y=J_max, x=Bin)) +
  labs(
    title = paste("Boxplot of Jmax per latitude"),
    x = "Latitude [°N]",
    y = "J,max [umol / m2 / s]"
  ) + 
  theme_minimal() 

ggplot(data.all %>% filter(Stage=="adult")) +
  geom_boxplot(aes(y=J_max, x=Bin)) +
  labs(
    title = paste("Boxplot of adult Jmax per latitude"),
    x = "Latitude [°N]",
    y = "J,max [umol / m2 / s]"
  ) + 
  theme_minimal() 

ggplot(data.all %>% filter(Stage=="sapling")) +
  geom_boxplot(aes(y=J_max, x=Bin)) +
  labs(
    title = paste("Boxplot of sapling Jmax per latitude"),
    x = "Latitude [°N]",
    y = "J,max [umol / m2 / s]"
  ) + 
  theme_minimal() 

anov.V <- aov(V_cmax ~ Bin, data = data.all)
summary(anov.V)
anov.J <- aov(J_max ~ Bin, data = data.all)
summary(anov.J)

manov <- manova(cbind(V_cmax, J_max) ~ Bin, data = data.all)
summary(manov)

summary.aov(manov)
TukeyHSD(anov.J, conf.level=.95)

### Check out if conditions for ANOVA apply:
## Normality: 
hist(data.all$V_cmax, breaks = 10)
hist(data.all$J_max, breaks = 10)

# V does not look great, J could maybe be normal:

qqnorm(anov.V$residuals)
qqline(anov.V$residuals)

qqnorm(anov.J$residuals)
qqline(anov.J$residuals)

shapiro.test(data.all$V_cmax) # Not normal
shapiro.test(data.all$J_max) # Not normal


library(MVN)
result <- mvn(data = data.all[, c("V_cmax", "J_max")], mvnTest = "mardia")
print(result$multivariateNormality) # not normal

## Equal variance (Levene test is more robust than Bartlett if normality is not certain): 
leveneTest(V_cmax ~ Bin, data=data.all) # Good
leveneTest(J_max ~ Bin, data=data.all) # Good

leveneTest(V_cmax ~ Bin, data=data.all %>% filter(Stage == "adult")) # Good
leveneTest(J_max ~ Bin, data=data.all %>% filter(Stage == "adult")) # Good

leveneTest(V_cmax ~ Bin, data=data.all %>% filter(Stage == "sapling")) # Good
leveneTest(J_max ~ Bin, data=data.all %>% filter(Stage == "sapling")) # Good



# Kruskal-Wallis test should be robust in this case:
kruskal.test(V_cmax ~ Bin, data = data.all) # Not significant
kruskal.test(J_max ~ Bin, data = data.all) # Not significant

kruskal.test(V_cmax ~ Bin, data = data.all %>% filter(Stage == "adult")) # Not significant
kruskal.test(J_max ~ Bin, data = data.all %>% filter(Stage == "adult")) # Not significant

kruskal.test(V_cmax ~ Bin, data = data.all %>% filter(Stage == "sapling")) # Not significant
kruskal.test(J_max ~ Bin, data = data.all %>% filter(Stage == "sapling")) # Not significant

#TODO: Do I need to check similar distribution requirement? Don't think so since the expected distribtuions should be equal since it is of the same species

# Don't think this is the best! No it is very sensitive to normality (which is not the case)
install.packages("biotools")
library(biotools)

box_m_test <- boxM(data.all[, c("V_cmax", "J_max")], data.all$Bin)
print(box_m_test)


## Robust alternative to Manova:
library(robustbase)

manova.robust <- lmrob(cbind("V_cmax", "J_max") ~ Bin, data = data.all)



dist_matrix <- dist(data.all[, c("V_cmax", "J_max")])

# Perform PERMANOVA
permanova <- adonis2(dist_matrix ~ Bin, data = data.all) # Not significant

# Summary of results
print(permanova)

dist_matrix_adult <- dist(data.all[data.all$Stage == "adult", c("V_cmax", "J_max")])
dist_matrix_sapling <- dist(data.all[data.all$Stage == "sapling", c("V_cmax", "J_max")])
?dist()
# Perform PERMANOVA
permanova_adult <- adonis2(dist_matrix_adult ~ Bin, data = data.all %>% filter(Stage == "adult"))
permanova_sapling <- adonis2(dist_matrix_sapling ~ Bin, data = data.all %>% filter(Stage == "sapling"))
?adonis2
# Summary of results
print(permanova_adult) # Not significant
print(permanova_sapling) # Not significant

dispersion <- betadisper(dist_matrix_adult, data.all %>% filter(Stage == "adult") %>% pull(Bin))
anova(dispersion)

#### Regression of Vcmax/Jmax by various parameters  ------

ggplot(data.all, aes(x = age, y = V_cmax, colour = Bin)) +
  geom_point() + geom_smooth(method = "lm", se = T) +
  labs(
    title = paste(""),
    x = "Age [yr]",
    y = "Vcmax  [umol/m2/s]"
  ) + 
  theme_minimal() 

# now A with different independent variables: 

ggplot(data.all, aes(x = mean20, y = `A_sim420-25`, colour = Bin)) +
  geom_point() + geom_smooth(method = "lm", se = T) +
  labs(
    title = paste("Amax at CO2 = 420 ppm, T = 25°"),
    x = "Average growth index for last 20 years",
    y = "Amax [umol/m2/s]"
  ) + 
  theme_minimal() 

# now Amax at different temperatures
simsC <- data.all %>% 
  filter(Stage == "adult") %>% 
  select(c(tree, Bin, colnames(data.all)[grepl("sim", colnames(data.all))])) %>% 
  pivot_longer(cols = -c(Bin, tree), names_to = "Scenario") %>% 
  separate_wider_delim(Scenario, delim = "_sim", names = c("Parameter", "Values")) %>% 
  separate_wider_delim(Values, delim = "-", names = c("Cair", "Tair")) %>% 
  pivot_wider(names_from = Parameter, values_from = value) %>% 
  group_by(Bin, Cair, Tair) %>% 
  summarise(Amean = mean(A, na.rm = T)) %>% 
  mutate(Cair=as.numeric(Cair))

ggplot(simsC %>% filter(Cair == 840), aes(x = Tair, y = Amean, colour = Bin)) +
  geom_point() + #geom_smooth(method = "lm", se = T) +
  labs(
    title = paste("At 840 ppm"),
    x = "Air temperature [°C]",
    y = "Mean Amax [umol/m2/s]"
  ) + 
  theme_minimal() 

ggplot(simsC %>% filter(Tair == 30), aes(x = Cair, y = Amean, colour = Bin)) +
  geom_point() + geom_smooth(method = "loess", se = F) +
  labs(
    title = paste("At Tair 30 °C"),
    x = "Ambient CO2 concentration [ppm]",
    y = "Mean Amax [umol/m2/s]"
  ) + 
  theme_minimal() 





#### After meeting:

# Let's compare adults and saplings:
# All
ggplot(data.all) +
  geom_boxplot(aes(y = V_cmax, x = Stage)) +
  labs(
    title = paste("Boxplot of Vcmax by stage"),
    x = "Stage",
    y = "Vc,max [umol / m2 / s]"
  ) + 
  theme_minimal() 

ggplot(data.all) +
  geom_boxplot(aes(y = J_max, x = Stage)) +
  labs(
    title = paste("Boxplot of Jmax by stage"),
    x = "Stage",
    y = "Vc,max [umol / m2 / s]"
  ) + 
  theme_minimal() 

ggplot(data.all %>% filter(Bin=="46.0")) +
  geom_boxplot(aes(y = V_cmax, x = Stage)) +
  labs(
    title = paste("Boxplot of Vcmax by stage for 46.0 °N"),
    x = "Latitude [°N]",
    y = "Vc,max [umol / m2 / s]"
  ) + 
  theme_minimal() 

ggplot(data.all %>% filter(Bin=="46.5")) +
  geom_boxplot(aes(y = V_cmax, x = Stage)) +
  labs(
    title = paste("Boxplot of Vcmax by stage for 46.5 °N"),
    x = "Latitude [°N]",
    y = "Vc,max [umol / m2 / s]"
  ) + 
  theme_minimal() 

ggplot(data.all %>% filter(Bin=="47.0")) +
  geom_boxplot(aes(y = V_cmax, x = Stage)) +
  labs(
    title = paste("Boxplot of Vcmax by stage for 47.0 °N"),
    x = "Latitude [°N]",
    y = "Vc,max [umol / m2 / s]"
  ) + 
  theme_minimal() 

ggplot(data.all %>% filter(Bin=="47.5")) +
  geom_boxplot(aes(y = V_cmax, x = Stage)) +
  labs(
    title = paste("Boxplot of Vcmax by stage for 47.5 °N"),
    x = "Latitude [°N]",
    y = "Vc,max [umol / m2 / s]"
  ) + 
  theme_minimal() 

ggplot(data.all %>% filter(Bin=="48.0")) +
  geom_boxplot(aes(y = V_cmax, x = Stage)) +
  labs(
    title = paste("Boxplot of Vcmax by stage for 48.0 °N"),
    x = "Latitude [°N]",
    y = "Vc,max [umol / m2 / s]"
  ) + 
  theme_minimal() 



ggplot(data.all %>% filter(Bin=="46.0")) +
  geom_boxplot(aes(y = J_max, x = Stage)) +
  labs(
    title = paste("Boxplot of Jmax by stage for 46.0 °N"),
    x = "Latitude [°N]",
    y = "Vc,max [umol / m2 / s]"
  ) + 
  theme_minimal() 

ggplot(data.all %>% filter(Bin=="46.5")) +
  geom_boxplot(aes(y = J_max, x = Stage)) +
  labs(
    title = paste("Boxplot of Jmax by stage for 46.5 °N"),
    x = "Latitude [°N]",
    y = "Vc,max [umol / m2 / s]"
  ) + 
  theme_minimal() 

ggplot(data.all %>% filter(Bin=="47.0")) +
  geom_boxplot(aes(y = J_max, x = Stage)) +
  labs(
    title = paste("Boxplot of Jmax by stage for 47.0 °N"),
    x = "Latitude [°N]",
    y = "Vc,max [umol / m2 / s]"
  ) + 
  theme_minimal() 

ggplot(data.all %>% filter(Bin=="47.5")) +
  geom_boxplot(aes(y = J_max, x = Stage)) +
  labs(
    title = paste("Boxplot of Jmax by stage for 47.5 °N"),
    x = "Latitude [°N]",
    y = "Vc,max [umol / m2 / s]"
  ) + 
  theme_minimal() 

ggplot(data.all %>% filter(Bin=="48.0")) +
  geom_boxplot(aes(y = J_max, x = Stage)) +
  labs(
    title = paste("Boxplot of Jmax by stage for 48.0 °N"),
    x = "Latitude [°N]",
    y = "Vc,max [umol / m2 / s]"
  ) + 
  theme_minimal() 

# So which statistical test do we need here? Manova it seems, so let's check assumptions:
# Normality:

result <- mvn(data = data.all %>% filter(Bin == "47.0", Stage == "adult") %>% select(c(V_cmax, J_max)), mvnTest = "mardia")
print(result$multivariateNormality) # not normal

hist(data.all %>% filter(Bin == "47.0", Stage == "adult") %>% pull(V_cmax), breaks = 4)

#Results: sapling and adult Vcmax not normal, Jmax normal
# 46.0-48.0: all normal but is this correct or jsut because such low

# variances are all good!

leveneTest(V_cmax ~ Stage, data = data.all %>% filter(Bin == "47.5")) # Good

leveneTest(J_max ~ Stage, data = data.all %>% filter(Bin == "47.5")) # Good

#### So let's do a Permanova:

dist_matrix_460 <- dist(data.all[data.all$Bin == "46.0", c("V_cmax", "J_max")])
dist_matrix_465 <- dist(data.all[data.all$Bin == "46.5", c("V_cmax", "J_max")])
dist_matrix_470 <- dist(data.all[data.all$Bin == "47.0", c("V_cmax", "J_max")])
dist_matrix_475 <- dist(data.all[data.all$Bin == "47.5", c("V_cmax", "J_max")])
dist_matrix_480 <- dist(data.all[data.all$Bin == "48.0", c("V_cmax", "J_max")])

# Perform PERMANOVA
permanova_460 <- adonis2(dist_matrix_460 ~ Stage, data = data.all %>% filter(Bin == "46.0"))
permanova_465 <- adonis2(dist_matrix_465 ~ Stage, data = data.all %>% filter(Bin == "46.5"))
permanova_470 <- adonis2(dist_matrix_470 ~ Stage, data = data.all %>% filter(Bin == "47.0"))
permanova_475 <- adonis2(dist_matrix_475 ~ Stage, data = data.all %>% filter(Bin == "47.5"))
permanova_480 <- adonis2(dist_matrix_480 ~ Stage, data = data.all %>% filter(Bin == "48.0"))

# Summary of results
print(permanova_460) # Not significant
print(permanova_465) # Not significant
print(permanova_470) # Not significant
print(permanova_475) # Not significant
print(permanova_480) # Not significant





