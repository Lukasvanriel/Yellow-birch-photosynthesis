## File to conduct first exploratory analysis
## Lukas Van Riel       17/12/2024
##### Load packages #####
library(tidyverse)
library(conflicted)
library(vegan)
library(car)

##### Load data #####
data.all <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/Combined.csv")
comm.tree <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP_TreeCommunity.csv")
comm.under <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP-UnderstoryCommunity.csv")

### Wrangle:
data.all <- data.all %>% 
  mutate(
    Bin = factor(format(round(as.numeric(substring(data.all$Bin, 2, 4)) / 10, 1), nsmall = 1))
  )

##### Functions #####


##### Body #####

#### Community PCA's ------
# Create tree community composition per plot
comm.tree.summ <- comm.tree %>% 
  select(-Plot) %>% 
  decostand(method = "hellinger") %>% 
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
ggplot(data.all) +
  geom_boxplot(aes(y=V_cmax, x=Bin)) +
  labs(
    title = paste("Boxplot of Vcmax per latitude"),
    x = "Latitude [°N]",
    y = "Vc,max [umol / m2 / s]"
  ) + 
  theme_minimal() 

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

# Perform PERMANOVA
permanova_adult <- adonis2(dist_matrix_adult ~ Bin, data = data.all %>% filter(Stage == "adult"))
permanova_sapling <- adonis2(dist_matrix_sapling ~ Bin, data = data.all %>% filter(Stage == "sapling"))

# Summary of results
print(permanova_adult) # Not significant
print(permanova_sapling) # Not significant



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

