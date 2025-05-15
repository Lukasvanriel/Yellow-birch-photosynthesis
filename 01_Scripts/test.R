# Load necessary library
library(vegan)
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
library(ggrepel)

# Simulate species abundance data (response matrix)
set.seed(42)
species_data <- data.frame(
  species1 = rpois(20, lambda = 10),
  species2 = rpois(20, lambda = 20),
  species3 = rpois(20, lambda = 5),
  species4 = rpois(20, lambda = 12)
)

# Simulate environmental variables (explanatory matrix)
environmental_data <- data.frame(
  pH = runif(20, 4.5, 7.5),
  moisture = runif(20, 10, 80),
  light = runif(20, 100, 1000)
)


# Perform CCA
cca_model <- cca(species_data ~ ., data = environmental_data)

# View summary
summary(cca_model)

# Plotting
plot(cca_model, scaling = 2, main = "CCA Biplot: Species vs Environment")


########
data.all <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/Combined.csv")
elevation <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Elevation/el.csv")
soil <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP-SoilLayers.csv") 
plot <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP-PlotInfo.csv") 

### Wrangle:
data.all <- data.all %>% 
  mutate(
    Bin = factor(format(round(as.numeric(substring(data.all$Bin, 2, 4)) / 10, 1), nsmall = 1))
  ) %>% 
  filter(tree != "P475.1-e") %>%  # The extra measurements
  filter(LeafArea != 40.5634)

data.elevation <- elevation %>% 
  mutate(tree = data.all$tree) %>% 
  rename(Elevation = Value) %>% 
  select(-ID) %>% 
  relocate(tree, .before = Elevation)

data.soil <- soil %>% 
  filter(Horizon == "B1") %>% 
  select(Plot, Horizon, soil_class, pH) %>% 
  rbind(c("P470.3", "B1", "LoSa", 4.91)) %>% 
  mutate(pH = as.numeric(pH))

# Remove some extra measuremnets (for now?)
data.all <- data.all %>% 
  filter(tree != "P475.1-e") %>%  # The extra measurements
  filter(LeafArea != 40.5634) # Double measurements for P480.1-9

data.cca <- data.all %>% 
  left_join(data.elevation, by = "tree") %>% 
  select(tree, Plot, Stage, V_cmax, J_max, diameter, LeafArea, Longitude, Latitude, age, Elevation)  %>% 
  left_join(data.soil, by = "Plot") %>% 
  mutate(Stage = factor(Stage),
         soil_class = factor(soil_class)) %>% 
  select(-Horizon) %>% 
  left_join(plot %>% select(Plot, Exposure, Slope), by = "Plot") %>% 
  mutate(Exposure = factor(Exposure),
         Slope = factor(Slope))
  
data.cca$Exposure

#### 
response.data <- data.cca %>% 
  select(V_cmax, J_max)

env.data <- data.cca %>% 
  select(-tree, -V_cmax, -J_max, -Plot,- age)

######### RDA ##############
response.data <- data.cca %>% 
  select(V_cmax, J_max)

env.data <- data.cca %>% 
  select(-tree, -Stage, -V_cmax, -J_max, -Plot,-age)

response.data.scaled <- response.data %>% 
  scale() %>% 
  as.data.frame()

env.data.scaled <- env.data %>% 
  select(-c(Exposure, Slope, soil_class)) %>% 
  scale() %>% 
  as.data.frame() %>% 
  mutate(Exposure = env.data$Exposure,
         Slope = env.data$Slope,
         soil_class = env.data$soil_class)

table(env.data.scaled$Exposure)/10

#### PCA on community compositions to add to RDA analysis:
comm.tree <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP_TreeCommunity.csv")
comm.under <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP-UnderstoryCommunity.csv")

comm.under.summ <- comm.under %>% 
  select(-Plot) %>% 
  decostand(method = "hellinger") %>% 
  bind_cols(select(comm.under, Plot))

comm.tree.summ <- comm.tree %>% 
  select(-Plot) %>% 
  decostand(method = "hellinger") %>% 
  bind_cols(select(comm.tree, Plot)) %>% 
  relocate(Plot, .before = everything())

##
pca.tree <- rda(comm.tree.summ %>% select(-Plot))
pca.under <- rda(comm.under.summ %>% select(-Plot))

summary(pca.tree)
scores(pca.tree)
summary(pca.under)
scores(pca.under)

plot(pca.tree, scaling=1)
plot(pca.tree, scaling=2)

ev <- pca.tree$CA$eig
ev[ev > mean(ev)]

ev <- pca.under$CA$eig
ev[ev > mean(ev)]

#Think only first 2 seems legit

head(bstick(pca.tree))
screeplot(pca.tree, 
          bstick = TRUE, type = "lines")
head(bstick(pca.under))
screeplot(pca.under, 
          bstick = TRUE, type = "lines")

### visualisation:
biplot(pca.tree)
biplot(pca.under)

biplot(pca.tree, scaling = 1)
biplot(pca.under, scaling = 1)




# Extract scores
site_scores <- scores(pca.tree, display = "sites")
species_scores <- scores(pca.tree, display = "species")

# Convert to data frames
site_df <- as.data.frame(site_scores)
site_df$Site <- comm.tree$Plot

species_df <- as.data.frame(species_scores)
species_df$Species <- rownames(species_df)

site_df$Latitude <- factor(rep(c("46.0", "46.5", "47.0", "47.5", "48.0"), each = 3),
                           levels = c("48.0", "47.5", "47.0", "46.5", "46.0"))

# Plot with ggplot2
main.plot <- ggplot() +
  geom_point(data = site_df, aes(x = PC1, y = PC2, colour = Latitude), size = 3) +
  scale_color_manual(values = c("46.0" = "#bd0026", "46.5" = "#f03b20", "47.0" = "#fdae61", "47.5" = "#abd9e9", "48.0" = "#2c7bb6")) +
  geom_text_repel(data = site_df, aes(x = PC1, y = PC2, label = Site), color = "blue") +
  geom_segment(data = species_df,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "darkgreen") +
  geom_text_repel(data = species_df, aes(x = PC1, y = PC2, label = Species), color = "darkgreen") +
  labs(x = "PC1 (43.4%)", y = "PC2 (21.9%)", title = "PCA biplot of tree community composition",color = "Latitude (°N)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position="inside",
    legend.position.inside = c(0.10, 0.45)
  )

library(cowplot)
species_legend <- "AB: Abies balsamea\nAR: Acer rubrum\nAS: Acer saccharum
BA: Betula alleghaniensis\nBP: Betula papyrifera\nFG: Fagus grandifolia
FN: Fraxinus nigra\nOV: Ostrya virginiana\nPG: Populus grandidentata
TA: Tilia americana\nTC: Tsuga canadensis\nTO: Thuja occidentalis"

ggdraw(main.plot) +
  draw_label(species_legend, x = 0.1, y = 0.94, hjust = 0, vjust = 1,
             size = 10, fontface = "plain", lineheight = 1.2)

###### Understory
# Extract scores
site_scores <- scores(pca.under, display = "sites")
species_scores <- scores(pca.under, display = "species")

# Convert to data frames
site_df <- as.data.frame(site_scores)
site_df$Site <- comm.tree$Plot

species_df <- as.data.frame(species_scores)
species_df$Species <- rownames(species_df)

# Plot with ggplot2
ggplot() +
  geom_point(data = site_df, aes(x = PC1, y = PC2), color = "blue", size = 3) +
  geom_text_repel(data = site_df, aes(x = PC1, y = PC2, label = Site), color = "blue") +
  geom_segment(data = species_df,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "darkred") +
  geom_text_repel(data = species_df, aes(x = PC1, y = PC2, label = Species), color = "darkred") +
  labs(x = "PC1 (28.3%)", y = "PC2 (16.3%)", title = "PCA biplot of tree community composition") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )


###### Add axes to dataset for RDA:
pca.ax.tree <- scores(pca.tree, display = "sites")[, 1:2]
pca.ax.under <- scores(pca.under, display = "sites")[, 1:2]

env.data.scaled <- env.data.scaled %>% 
  mutate(pca.tree1 = rep(pca.ax.tree[,1], each=10),
         pca.tree2 = rep(pca.ax.tree[,2], each=10),
         pca.under1 = rep(pca.ax.under[,1], each=10),
         pca.under2 = rep(pca.ax.under[,2], each=10))

summary(env.data.scaled)

### Run RDA ####
rda.test <- rda(response.data.scaled ~., env.data.scaled)

if(F){
rda.under <- rda(response.data.scaled ~., env.data.scaled %>% select(pca.tree1,
                                                                     pca.tree2, pca.under1, pca.under2))
summary(rda.under)
anova.cca(rda.under, by="axis", step = 100000)
anova.cca(rda.under, by="terms", step = 100000)
}

# collinearity?
heatmap(abs(cor(env.data.scaled %>% select(-Exposure, -Slope, -soil_class))), ## Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))

vif.cca(rda.test)

rda <- rda(response.data.scaled ~., env.data.scaled %>% select(-pca.under2, -pca.tree2, -pca.under1))
heatmap(abs(cor(env.data.scaled %>% select(-Exposure, -Slope, -soil_class, -pca.under2, -pca.tree2, -pca.under1))), ## Compute pearson correlation (note they are absolute values)
         col = rev(heat.colors(6)), 
         Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))

#########Interpret:
##########
summary(rda)
coef(rda)
scores(rda)

R_adj <- RsquareAdj(rda)

R_adj$r.squared
R_adj$adj.r.squared

plot(rda)

var.sc <- scores(rda, choices = 1:2,display = "sp")
plot(rda)
arrows(0,0,var.sc[,1], var.sc[,2], length = 0, lty = 1, col="red")

plot(rda)
arrows(0,0,var.sc[,1], var.sc[,2], length = 0, lty = 1, col="red")


plot(rda, display = c("sp", "lc", "cn"))
#Scaling2 
plot(rda, scaling=1, display = c("lc", "cn"))

plot(env.data$Longitude, response.data$J_max)
lmodel <- lm(response.data$J_max~env.data$Longitude)
summary(lmodel)
coef(lmodel)


# Significance: 
anova.cca(rda, step = 10000)
anova.cca(rda, by="axis", step = 10000)
anova.cca(rda, by="terms", step = 10000)

rda$CA$eig[rda$CA$eig > mean(rda$CA$eig)]


#### Try forward modelling:

fwd.sel <- ordiR2step(rda(response.data.scaled ~ 1, data = env.data.scaled %>% select(-pca.under2, -pca.tree2, -pca.under1)), # lower model limit (simple!)
                      scope = formula(rda), # upper model limit (the "full" model)
                      direction = "forward",
                      R2scope = TRUE, # can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = FALSE)

fwd.sel$call

spe.rda.signif <- rda(response.data.scaled ~ Exposure, env.data.scaled %>% select(-pca.under2, -pca.tree2, -pca.under1))
# check the adjusted R2
RsquareAdj(spe.rda.signif)

anova.cca(spe.rda.signif, step = 1000)


ordiplot(spe.rda.signif,
         scaling = 2)

ordiplot(spe.rda.signif,
         scaling = 2)

spe.rda.signif <- rda(response.data.scaled ~ Longitude, env.data.scaled)
# check the adjusted R2
RsquareAdj(spe.rda.signif)

anova.cca(spe.rda.signif, step = 1000)


##### Plot RDA:

# Extract site and variable scores
site_scores <- scores(rda, display = "sites")
response_scores <- scores(rda, display = "species")  # traits
biplot_scores <- scores(rda, display = "bp")  # explanatory variables

# Convert to data frames
site_df <- as.data.frame(site_scores)
site_df$Site <- rownames(site_df)
site_df$Stage <- rep(rep(c("adult", "sapling"), each=5), 15)

response_df <- as.data.frame(response_scores)
response_df$Trait <- rownames(response_df)

biplot_df <- as.data.frame(biplot_scores)
biplot_df$Variable <- rownames(biplot_df)

# Plot RDA biplot
ggplot() +
  # Plot sites
  geom_point(data = site_df, aes(x = RDA1, y = RDA2, col=Stage), size = 1.5) +
  #geom_text_repel(data = site_df, aes(x = RDA1, y = RDA2, label = Site), size = 3, color = "black") +
  
  # Plot trait arrows
  geom_segment(data = response_df,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "blue") +
  geom_text_repel(data = response_df, aes(x = RDA1, y = RDA2, label = Trait), color = "blue") +
  
  # Plot explanatory variable arrows
  geom_segment(data = biplot_df,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "darkred") +
  geom_text_repel(data = biplot_df, aes(x = RDA1, y = RDA2, label = Variable), color = "darkred") +
  
  labs(x = "RDA1 (24.3%)", y = "RDA2 (0.3%)", title = "RDA Biplot") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )



###### What if we use slightly other Exposure?


############ Now check on the 
env.data.scaled


boxplot(V_cmax ~ Exposure, data = cbind(env.data, response.data), title="Test")


# Convert data to long format for ggplot
df_comb <- cbind(response.data, env.data)
  
  
  
  data.all %>%
  pivot_longer(cols = c(V_cmax, J_max), names_to = "Variable", values_to = "Value") %>% 
  mutate(Stage = ifelse(Stage=="adult", "Adults", "Saplings"),
         Variable = factor(Variable, levels = c("V_cmax", "J_max"),
                           labels = c(expression(paste(V[paste(c, ",", max, sep = "")], ~(mu*mol *~ m^{-2} *~ s^{-1}))),
                                      expression(paste(J[max], ~(mu*mol *~ m^{-2} *~ s^{-1}))) ) ))


# Create boxplot
ggplot(df_comb, aes(x = Exposure, y = V_cmax)) +
  geom_boxplot(fill="grey", alpha = 0.9, outlier.shape = NA) + 
  labs(#title = "Distribution of Vcmax and Jmax by stage and Latitude",
    x = "Exposure", y = expression(paste(V[paste(c, ",", max, sep = "")], ~(mu*mol~m^{-2}~s^{-1})))) + #expression(paste(V[paste(c, ",", max, sep = "")], ~(mu*mol m^{-2} s^-1)))) +
  theme_minimal(base_size = 12)  # Clean theme
#######


a <- expression(paste(V[paste(c, ",", max, sep = "")], ~(mu*mol m^{-2} s^-1)))

data.adult <- data.all %>% 
  filter(Stage == "adult")
data.sapling <- data.all %>% 
  filter(Stage == "sapling")

dist.matrix.adult <- dist(response.data, method = "euclidean")
dist.matrix.sapling <- dist(, method = "euclidean")
dist.matrix.all <- dist(response.data)

permanova.adult <- adonis2(dist.matrix.adult ~ Bin, data = data.adult)
permanova.sapling <- adonis2(dist.matrix.sapling ~ Bin, data = data.sapling)
permanova.all <- adonis2(dist.matrix.all ~ Exposure, data = cbind(response.data, env.data))

dispersion <- betadisper(dist.matrix.all, cbind(response.data, env.data) %>% pull(Exposure))
anova(dispersion)


##### Variance partitioning:

colnames(env.data.scaled)
phys.data.scaled <- env.data.scaled %>% 
  select(Longitude, Latitude, Elevation, Exposure, Slope)
biol.data.scaled <- env.data.scaled %>% 
  select(diameter, LeafArea, pH, soil_class, pca.tree1)

rda.part <- varpart(response.data.scaled, phys.data.scaled, biol.data.scaled)
rda.part$part

plot(rda.part,
     Xnames = c("Geo/Topo", "Bio/Chem"), # name the partitions
     bg = c("orange4", "seagreen3"), alpha = 90, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.2)

## Test significance:
anova.cca(rda(response.data.scaled, phys.data.scaled))
anova.cca(rda(response.data.scaled, biol.data.scaled))
anova.cca(rda(response.data.scaled, phys.data.scaled, biol.data.scaled))
anova.cca(rda(response.data.scaled, biol.data.scaled, phys.data.scaled))







###
summary(rda.test)
coef(rda.test)

R_adj <- RsquareAdj(rda.test)

R_adj$r.squared
R_adj$adj.r.squared

plot(rda.test)

var.sc <- scores(rda.test, choices = 1:2,display = "sp")
plot(rda.test)
arrows(0,0,var.sc[,1], var.sc[,2], length = 0, lty = 1, col="red")

plot(rda.test)
arrows(0,0,var.sc[,1], var.sc[,2], length = 0, lty = 1, col="red")


plot(rda.test, display = c("sp", "lc", "cn"))
#Scaling2 
plot(rda.test, scaling=1, display = c("lc", "cn"))

plot(env.data$Longitude, response.data$J_max)
lmodel <- lm(response.data$J_max~env.data$Longitude)
summary(lmodel)
coef(lmodel)


# Significance: 
anova.cca(rda.test, by="axis", step = 100000)
anova.cca(rda.test, by="terms", step = 100000)

rda.test$CA$eig[rda.test$CA$eig > mean(rda.test$CA$eig)]















### PCoA:
library(ape)

comm.tree <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP_TreeCommunity.csv")
comm.under <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP-UnderstoryCommunity.csv")


comm.under.summ <- comm.under %>% 
  select(-Plot) %>% 
  decostand(method = "hellinger") %>% 
  bind_cols(select(comm.under, Plot))

comm.tree.summ <- comm.tree %>% 
  select(-Plot) %>% 
  decostand(method = "hellinger") %>% 
  bind_cols(select(comm.tree, Plot)) %>% 
  relocate(Plot, .before = everything())

community.data <- data.frame(tree=1:150)

pcoa.tree <- pcoa(dist(comm.tree.summ %>% select(-Plot)))
summary(pcoa.tree)
head(pcoa.tree$values)

biplot.pcoa(pcoa.tree, comm.tree.summ %>% select(-Plot))

community.data$tree.ax1 <- rep(pcoa.tree$vectors[,1], each=10)
community.data$tree.ax2 <- rep(pcoa.tree$vectors[,2], each=10)

pcoa.under <- pcoa(dist(comm.under.summ %>% select(-Plot)))
summary(pcoa.under)
head(pcoa.under$values)

biplot.pcoa(pcoa.under, comm.under.summ %>% select(-Plot))

community.data$under.ax1 <- rep(pcoa.under$vectors[,1], each=10)
community.data$under.ax2 <- rep(pcoa.under$vectors[,2], each=10)





########## PCA FOR CEF in French
biplot(pca.tree)
biplot(pca.under)

biplot(pca.tree, scaling = 1)
biplot(pca.under, scaling = 1)




# Extract scores
site_scores <- scores(pca.tree, display = "sites")
species_scores <- scores(pca.tree, display = "species")

# Convert to data frames
site_df <- as.data.frame(site_scores)
site_df$Site <- comm.tree$Plot

species_df <- as.data.frame(species_scores)
species_df$Species <- c("ERS", "BOJ", "HEG", "PEG", "TIL", "PRU", "SAB", "OSV", "BOP", "THO", "ERR", "FRN")
  

site_df$Latitude <- factor(rep(c("46.0", "46.5", "47.0", "47.5", "48.0"), each = 3),
                           levels = c("48.0", "47.5", "47.0", "46.5", "46.0"))

# Plot with ggplot2
main.plot <- ggplot() +
  geom_point(data = site_df, aes(x = PC1, y = PC2, colour = Latitude), size = 3) +
  scale_color_manual(values = c("46.0" = "#bd0026", "46.5" = "#f03b20", "47.0" = "#fdae61", "47.5" = "#abd9e9", "48.0" = "#2c7bb6")) +
  geom_text_repel(data = site_df, aes(x = PC1, y = PC2, label = Site), color = "blue") +
  geom_segment(data = species_df,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "darkgreen") +
  geom_text_repel(data = species_df, aes(x = PC1, y = PC2, label = Species), color = "darkgreen") +
  labs(x = "PC1 (43.4%)", y = "PC2 (21.9%)", title = "PCA biplot of tree community composition",color = "Latitude (°N)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position="inside",
    legend.position.inside = c(0.10, 0.45)
  )

library(cowplot)
species_legend <- "AB: Abies balsamea\nAR: Acer rubrum\nAS: Acer saccharum
BA: Betula alleghaniensis\nBP: Betula papyrifera\nFG: Fagus grandifolia
FN: Fraxinus nigra\nOV: Ostrya virginiana\nPG: Populus grandidentata
TA: Tilia americana\nTC: Tsuga canadensis\nTO: Thuja occidentalis"

sort(species_df$Species)
species_legend <-"BOJ: Betula alleghaniensis\nBOP: Betula papyrifera\nERR: Acer rubrum
ERS: Acer saccharum\nFRN: Fraxinus nigra\nHEG: Fagus grandifolia
OSV: Ostrya virginiana\nPEG: Populus grandidentata\nPRU: Tsuga canadensis
SAB: Abies balsamea\nTHO: Thuja occidentalis\nTIL: Tilia americana"

ggdraw(main.plot) +
  draw_label(species_legend, x = 0.1, y = 0.94, hjust = 0, vjust = 1,
             size = 10, fontface = "plain", lineheight = 1.2)








### Multivariate assessment
library(mvabund)

fit <- manyglm(response.data.scaled ~ tree.ax1 + tree.ax2 + under.ax1 + under.ax2, 
               data = community.data %>% select(-tree), family = "gamma")
summary(fit)
anova(fit, resamp = "pit.trap")

mvabund(response.data.scaled)


glm.commV <- glm(response.data.scaled$V_cmax ~ tree.ax1 + tree.ax2 + under.ax1 + under.ax2, 
    data = community.data %>% select(-tree))
summary(glm.commV)

glm.commJ <- glm(response.data.scaled$J_max ~ tree.ax1 + tree.ax2 + under.ax1 + under.ax2, 
                data = community.data %>% select(-tree))
summary(glm.commJ)

####Forward modelling:


fwd.sel <- ordiR2step(rda(response.data.scaled ~ 1, data = env.data.scaled), # lower model limit (simple!)
                      scope = formula(rda.test), # upper model limit (the "full" model)
                      direction = "forward",
                      R2scope = TRUE, # can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = FALSE)

fwd.sel$call

spe.rda.signif <- rda(response.data.scaled ~ Exposure, env.data.scaled)
# check the adjusted R2
RsquareAdj(spe.rda.signif)

anova.cca(spe.rda.signif, step = 1000)


ordiplot(spe.rda.signif,
         scaling = 2)

ordiplot(spe.rda.signif,
         scaling = 2)

spe.rda.signif <- rda(response.data.scaled ~ Longitude, env.data.scaled)
# check the adjusted R2
RsquareAdj(spe.rda.signif)

anova.cca(spe.rda.signif, step = 1000)



########OLD#####

library(yacca)
ccacor <- cancor(env.data, response.data)
#X-morph-env
#Y-phys-respon


cca <- cca(response.data, env.data)  # note the reversed order in yacca
summary(cca_model)

U <- as.matrix(env.data) %*% ccacor$xcoef  # canonical variates for X
V <- as.matrix(response.data) %*% ccacor$ycoef  # canonical variates for Y


plot(U[, 3], V[, 2],
     xlab = "Canonical Variate 1 (X set)",
     ylab = "Canonical Variate 1 (Y set)",
     main = "Canonical Correlation Scatterplot",
     pch = 19, col = "steelblue")
abline(lm(V[, 2] ~ U[, 3]), col = "red")


cor_X_U1 <- cor(env.data, U[, 1])
cor_Y_V1 <- cor(response.data, V[, 1])

barplot(cor_X_U1, main = "X vars vs Canonical Variate 1", ylim = c(-1, 1))
barplot(cor_Y_V1, main = "Y vars vs Canonical Variate 1", ylim = c(-1, 1))

library(ggplot2)

df_scores <- data.frame(U1 = U[, 1], V1 = V[, 1])
ggplot(df_scores, aes(x = U1, y = V1)) +
  geom_point(color = "darkgreen") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "Canonical Variate 1 (X set)", y = "Canonical Variate 1 (Y set)",
       title = "Canonical Correlation Plot") +
  theme_minimal()

##### Biplot?
structure_X <- cor(env.data, U)  # structure coefficients for X
structure_Y <- cor(response.data, V)

scores <- data.frame(U1 = U[, 1], V1 = V[, 1])
vars_X <- data.frame(structure_X[, 1, drop = FALSE])
vars_Y <- data.frame(structure_Y[, 1, drop = FALSE])

# Add variable names
vars_X$var <- rownames(vars_X)
vars_Y$var <- rownames(vars_Y)
colnames(vars_X)[1] <- "cor"
colnames(vars_Y)[1] <- "cor"

# Biplot: Canonical scores and variable vectors
ggplot(scores, aes(x = U1, y = V1)) +
  geom_point(color = "steelblue", size = 2) +
  geom_segment(data = vars_X, aes(x = 0, y = 0, xend = cor * 1, yend = 0), 
               arrow = arrow(length = unit(0.2, "cm")), color = "tomato") +
  #geom_text_repel(data = vars_X, aes(x = cor * 5, y = 0, label = var), 
  #                color = "tomato", direction = "y") +
  geom_segment(data = vars_Y, aes(x = 0, y = 0, xend = 0, yend = cor * 1), 
               arrow = arrow(length = unit(0.2, "cm")), color = "darkgreen") +
  #geom_text_repel(data = vars_Y, aes(x = 0, y = cor * 5, label = var), 
  #                color = "darkgreen", direction = "x") +
  labs(x = "Canonical Variate 1 (X set)", y = "Canonical Variate 1 (Y set)",
       title = "Canonical Correlation Biplot") +
  theme_minimal()




plot(cca)
structure.plot(cca)


scores <- comp.scores(cca)

plot(scores$cancor.scores$X.scores[, 1],
     scores$cancor.scores$Y.scores[, 1],
     xlab = "Canonical Variate 1 (X set)",
     ylab = "Canonical Variate 1 (Y set)",
     main = "Canonical Scores Plot",
     pch = 19, col = "dodgerblue")

abline(lm(scores$cancor.scores$Y.scores[, 1] ~ scores$cancor.scores$X.scores[, 1]),
       col = "red", lwd = 2)

