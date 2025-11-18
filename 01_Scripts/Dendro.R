## File to conduct dendro analysis
## Lukas Van Riel       13/11/2025
##### Load packages #####
library(tidyverse)
library(conflicted)
library(dplR)
library(treeclim)

conflicts_prefer(dplyr::filter)

# Set path
rwi_path <- "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/TreeCores-Detrended"

##### Load data #####
# Get all detrended files
rwi_files <- list.files(rwi_path, pattern = "-Detrended-ModNegExp.csv", full.names = TRUE)

# Load and build chronologies
chronologies <- list()

# Load without setting row names first
for (file in rwi_files) {
  plot_name <- str_extract(basename(file), "P\\d+_\\d")
  
  # Read without row names
  rwi_raw <- read.csv(file)
  
  # Check structure
  if (plot_name == names(chronologies)[1] || length(chronologies) == 0) {
    print(head(rwi_raw))  # See what we're working with
    print(str(rwi_raw))
  }
  
  # Set years as row names manually
  rownames(rwi_raw) <- rwi_raw$year
  rwi <- rwi_raw[, -which(names(rwi_raw) == "year")]  # Remove year column
  
  # Build chronology
  crn <- chron(rwi, prefix = plot_name)
  chronologies[[plot_name]] <- crn
  
  cat(plot_name, ":", nrow(crn), "years,", ncol(rwi), "trees\n")
}


# Set climate path
clim_path <- "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/ClimateData/Monthly"

# Load monthly data
temp_monthly <- read.csv(file.path(clim_path, "Tmean.csv"))
precip_monthly <- read.csv(file.path(clim_path, "TotPrec.csv"))

# Filter historical period only (pre-2024)
temp_hist <- temp_monthly %>% filter(year <= 2024)
precip_hist <- precip_monthly %>% filter(year <= 2024)

# Check
range(temp_hist$year)  # Should be 1950-2023
unique(temp_hist$Plot)  # Should match your plot names


### Correlatio nanalysis#######
# Storage for results
climate_sensitivity <- data.frame()

for (plot_name in names(chronologies)) {
  crn <- chronologies[[plot_name]]
  
  # Prepare climate data
  temp_plot <- temp_hist %>% 
    filter(Plot == plot_name) %>%
    select(year, month, temp = ssp245_tg_mean_p50)
  
  precip_plot <- precip_hist %>%
    filter(Plot == plot_name) %>%
    select(year, month, precip = ssp245_prcptot_p50)
  
  # Run correlations
  dcc_temp <- dcc(crn, temp_plot, 
                  selection = -6:9, 
                  method = "correlation",
                  dynamic = "static",
                  var_names = "temp")
  
  dcc_precip <- dcc(crn, precip_plot,
                    selection = -6:9,
                    method = "correlation",
                    dynamic = "static",
                    var_names = "precip")
  
  # Extract using rownames
  summer_temp_cor <- mean(dcc_temp$coef[c("temp.curr.jun", "temp.curr.jul", "temp.curr.aug"), "coef"])
  spring_precip_cor <- mean(dcc_precip$coef[c("precip.curr.apr", "precip.curr.may"), "coef"])
  
  climate_sensitivity <- rbind(climate_sensitivity, 
                               data.frame(Plot = plot_name,
                                          Summer_temp_cor = summer_temp_cor,
                                          Spring_precip_cor = spring_precip_cor))
  
  cat(plot_name, "done\n")
}

# View results
print(climate_sensitivity)


#### Link to physiology

# Extract bin and latitude from plot names
climate_sensitivity <- climate_sensitivity %>%
  mutate(
    Bin = factor(as.numeric(str_extract(Plot, "\\d{3}")) / 10),
    Plot_num = as.numeric(str_extract(Plot, "\\d$"))
  )

# Create site physiology with correct columns
site_physiology <- data.all %>%
  filter(Stage == "adult") %>%  # lowercase!
  group_by(Plot) %>%
  summarize(
    Mean_Vcmax = mean(V_cmax, na.rm = TRUE),
    Mean_Jmax = mean(J_max, na.rm = TRUE),
    Mean_JV_ratio = mean(J_max / V_cmax, na.rm = TRUE),
    Latitude = first(Latitude),
    Slope = first(Slope),
    Exposure = first(Exposure),
    Bin = first(Bin),
    n_trees = n(),
    .groups = "drop"
  )
# Convert climate_sensitivity plot names to match data.all format
climate_sensitivity$Plot <- gsub("_", ".", climate_sensitivity$Plot)

# Convert both to numeric
climate_sensitivity$Bin <- as.numeric(as.character(climate_sensitivity$Bin))
site_physiology$Bin <- as.numeric(as.character(site_physiology$Bin))

# Now join
phys_growth <- left_join(climate_sensitivity, site_physiology, 
                         by = c("Plot", "Bin"))
print(phys_growth)



# Key plot 1: Vcmax vs Temperature Sensitivity by Slope
ggplot(phys_growth, aes(Mean_Vcmax, Summer_temp_cor, color = Slope)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Mean Vcmax", 
       y = "Summer Temperature Sensitivity (r)",
       title = "Does low Vcmax predict high climate sensitivity?")

# Key plot 2: Slope effect
ggplot(phys_growth, aes(Slope, Summer_temp_cor)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 2) +
  labs(y = "Summer Temperature Sensitivity")

# Quick stats
cor.test(phys_growth$Mean_Vcmax, phys_growth$Summer_temp_cor)

# Model
lm_sensitivity <- lm(Summer_temp_cor ~ Mean_Vcmax + Slope + Latitude, 
                     data = phys_growth)
summary(lm_sensitivity)


# The key figure for your paper
phys_growth$Vulnerability <- ifelse(phys_growth$Slope == "B", "Steep (Vulnerable)", "Moderate")

ggplot(phys_growth, aes(Mean_Vcmax, Summer_temp_cor)) +
  geom_point(aes(color = Slope, size = 3)) +
  geom_smooth(method = "lm") +
  labs(x = "Mean Vcmax (μmol m⁻² s⁻¹)", 
       y = "Summer Temperature Sensitivity",
       title = "Physiological capacity predicts climate vulnerability")





######## What other things to explore?
# JV ratio vs climate sensitivity

# Stats
cor.test(phys_growth$Mean_JV_ratio, phys_growth$Summer_temp_cor)

# Model with JV ratio
lm_jv <- lm(Summer_temp_cor ~ Mean_JV_ratio + Slope + Latitude, 
            data = phys_growth)
summary(lm_jv)

ggplot(phys_growth, aes(Mean_JV_ratio, Summer_temp_cor, color = Slope)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  labs(x = "Jmax:Vcmax Ratio", y = "Summer Temperature Sensitivity")

# Correlation between slope type and 
# JV ratio by slope
ggplot(phys_growth, aes(Slope, Mean_JV_ratio)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 2) +
  labs(y = "Jmax:Vcmax Ratio")

# Stats
summary(aov(Mean_JV_ratio ~ Slope, data = phys_growth))

# Or as table
phys_growth %>%
  group_by(Slope) %>%
  summarize(Mean_JV = mean(Mean_JV_ratio),
            SD_JV = sd(Mean_JV_ratio))

###### Moving windwo analysis:
## Two periods first:
# Compare early vs late periods
temporal_change <- data.frame()

for (plot_name in names(chronologies)) {
  crn <- chronologies[[plot_name]]
  
  # Get climate data (no gsub!)
  temp_plot <- temp_hist %>% 
    filter(Plot == plot_name) %>%
    select(year, month, temp = ssp245_tg_mean_p50)
  
  # Split into two periods
  temp_early <- temp_plot %>% filter(year >= 1950, year <= 1987)
  temp_late <- temp_plot %>% filter(year >= 1988, year <= 2024)
  
  # Early period correlation
  dcc_early <- dcc(crn, temp_early,
                   selection = -6:9,
                   method = "correlation",
                   dynamic = "static",
                   var_names = "temp")
  
  # Late period correlation
  dcc_late <- dcc(crn, temp_late,
                  selection = -6:9,
                  method = "correlation",
                  dynamic = "static",
                  var_names = "temp")
  
  # Extract summer temperature sensitivity
  summer_early <- mean(dcc_early$coef[c("temp.curr.jun", "temp.curr.jul", "temp.curr.aug"), "coef"])
  summer_late <- mean(dcc_late$coef[c("temp.curr.jun", "temp.curr.jul", "temp.curr.aug"), "coef"])
  
  temporal_change <- rbind(temporal_change,
                           data.frame(Plot = plot_name,
                                      Early_1950_1987 = summer_early,
                                      Late_1988_2024 = summer_late,
                                      Change = summer_late - summer_early))
  
  cat(plot_name, "done\n")
}

# Merge with plot characteristics
temporal_change$Plot_fixed <- gsub("_", ".", temporal_change$Plot)
temporal_change <- temporal_change %>%
  left_join(phys_growth %>% select(Plot, Bin, Slope, Mean_Vcmax), 
            by = c("Plot_fixed" = "Plot"))

# View results
print(temporal_change)

# Is sensitivity becoming more negative over time?
t.test(temporal_change$Change)

# By slope
temporal_change %>%
  group_by(Slope) %>%
  summarize(Mean_change = mean(Change),
            SD = sd(Change))

ggplot(temporal_change, aes(x = Early_1950_1987, y = Late_1988_2024, color = Slope)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(size = 4) +
  labs(x = "Early Period (1950-1987)", 
       y = "Late Period (1988-2024)",
       title = "Temporal change in climate sensitivity") +
  annotate("text", x = -Inf, y = Inf, 
           label = "Below line = MORE NEGATIVE (increasing stress)", 
           hjust = -0.1, vjust = 1.5)


# Overall trend test
t.test(temporal_change$Change)

# By slope
temporal_change %>%
  group_by(Slope) %>%
  summarize(Mean_change = mean(Change),
            SD = sd(Change),
            n = n())

# Does Vcmax predict temporal change?
cor.test(temporal_change$Mean_Vcmax, temporal_change$Change)

ggplot(temporal_change, aes(Mean_Vcmax, Change, color = Slope)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  labs(x = "Mean Vcmax", 
       y = "Change in Temperature Sensitivity (Late - Early)",
       title = "Do low Vcmax sites show increasing stress?")

## Real running window:
# Try treeclim's moving window (30-year windows, 5-year steps)
moving_window_results <- data.frame()

for (plot_name in names(chronologies)) {
  crn <- chronologies[[plot_name]]
  
  temp_plot <- temp_hist %>% 
    filter(Plot == plot_name) %>%
    select(year, month, temp = ssp245_tg_mean_p50)
  
  tryCatch({
    dcc_moving <- dcc(crn, temp_plot,
                      selection = -6:9,
                      dynamic = "moving",
                      win_size = 30,
                      win_offset = 5,
                      var_names = "temp")
    
    # Get coefficient dataframe and month info
    coef_df <- dcc_moving$coef$coef
    month_info <- dcc_moving$coef$pretty_names
    
    # Find summer months (current year = uppercase)
    summer_rows <- which(month_info$month_label %in% c("JUN", "JUL", "AUG"))
    
    # Extract for each window
    window_names <- colnames(coef_df)
    
    for (window in window_names) {
      year <- as.numeric(sub("-.*", "", window)) + 15  # Middle year
      summer_cor <- mean(coef_df[summer_rows, window])
      
      moving_window_results <- rbind(moving_window_results,
                                     data.frame(Plot = plot_name,
                                                Window = window,
                                                Year = year,
                                                Summer_cor = summer_cor))
    }
    
    cat(plot_name, "done\n")
    
  }, error = function(e) {
    cat(plot_name, "FAILED:", e$message, "\n")
  })
}

# Merge with plot info
moving_window_results$Plot_fixed <- gsub("_", ".", moving_window_results$Plot)
moving_window_results <- moving_window_results %>%
  left_join(phys_growth %>% select(Plot, Slope, Mean_Vcmax, Bin), 
            by = c("Plot_fixed" = "Plot"))

# Visualize
ggplot(moving_window_results, aes(x = Year, y = Summer_cor, group = Plot, color = Slope)) +
  geom_line(alpha = 0.6) +
  geom_smooth(aes(group = Slope), method = "lm", se = FALSE, linewidth = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year (window midpoint)", 
       y = "Summer Temperature Sensitivity",
       title = "Temporal trends by slope")



# Test for temporal trend
lm_trend <- lm(Summer_cor ~ Year * Slope, data = moving_window_results)
summary(lm_trend)

# By slope
moving_window_results %>%
  group_by(Slope) %>%
  summarize(
    Trend_slope = coef(lm(Summer_cor ~ Year))[2],
    p_value = summary(lm(Summer_cor ~ Year))$coefficients[2,4]
  )

# Visualize individual plots
ggplot(moving_window_results, aes(x = Year, y = Summer_cor, color = Slope)) +
  geom_line(aes(group = Plot), alpha = 0.4) +
  geom_smooth(aes(group = Slope), method = "loess", se = TRUE, linewidth = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~Slope) +
  labs(x = "Year", y = "Temperature Sensitivity",
       title = "Non-linear temporal trends")



# Temporal trend by latitude bin
moving_window_results %>%
  group_by(Bin) %>%
  summarize(
    Trend_slope = coef(lm(Summer_cor ~ Year))[2],
    p_value = summary(lm(Summer_cor ~ Year))$coefficients[2,4],
    n_plots = n_distinct(Plot)
  )

# Test for latitude × time interaction
lm_lat_time <- lm(Summer_cor ~ Year * Bin, data = moving_window_results)
summary(lm_lat_time)

# Visualize
ggplot(moving_window_results, aes(x = Year, y = Summer_cor, color = factor(Bin))) +
  geom_line(aes(group = Plot), alpha = 0.4) +
  geom_smooth(aes(group = Bin), method = "lm", se = TRUE, linewidth = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year", y = "Temperature Sensitivity",
       title = "Temporal trends by latitude",
       color = "Latitude Bin")

# Probably confounding:
table(phys_growth$Bin, phys_growth$Slope)


### BAI:

## 0:
# Track how much we're losing per tree
bark_check <- data.frame()

for (plot_name in names(chronologies)) {
  raw_file <- file.path(raw_path, paste0(plot_name, "c.rwl"))
  
  if (file.exists(raw_file)) {
    rwl <- read.rwl(raw_file)
    
    plot_trees <- data.all %>%
      filter(Stage == "adult", Plot == gsub("_", ".", plot_name))
    
    for (i in 1:nrow(plot_trees)) {
      tree_id <- plot_trees$Tree[i]
      
      if (tree_id %in% colnames(rwl)) {
        current_radius <- plot_trees$Inside_bark_radius_cm[i]
        dbh <- plot_trees$DBH_cm[i]
        
        # Get ring widths
        ring_data <- rwl[, tree_id]
        valid_rows <- !is.na(ring_data)
        ring_widths <- ring_data[valid_rows] / 10
        years <- as.numeric(rownames(rwl)[valid_rows])
        
        # Calculate radius backwards
        n_rings <- length(ring_widths)
        radius <- rep(NA, n_rings)
        radius[n_rings] <- current_radius
        
        for (j in (n_rings-1):1) {
          radius[j] <- radius[j+1] - ring_widths[j+1]
        }
        
        # Check losses
        n_negative <- sum(radius < 0)
        total_ring_width <- sum(ring_widths)
        
        bark_check <- rbind(bark_check,
                            data.frame(
                              Tree = tree_id,
                              DBH_cm = dbh,
                              Current_radius = current_radius,
                              Total_ring_sum = total_ring_width,
                              Difference = total_ring_width - current_radius,
                              Years_lost = n_negative,
                              Pct_lost = n_negative / n_rings * 100
                            ))
      }
    }
  }
}

# Summary
summary(bark_check$Years_lost)
summary(bark_check$Pct_lost)
summary(bark_check$Difference)

# How many trees losing >20% of rings?
sum(bark_check$Pct_lost > 20)

# Plot
hist(bark_check$Difference, breaks = 20,
     main = "Ring sum exceeds current radius",
     xlab = "Difference (cm)")
abline(v = 0, col = "red", lwd = 2)

# If many trees have large positive difference, bark correction is too aggressive
mean(bark_check$Difference)

# ============================================
# STEP 1: Bark Correction (adults only)
# ============================================

data.all <- data.all %>%
  mutate(
    DBH_cm = diameter / 10,  # mm to cm
    Bark_thickness_cm = 0.03 * DBH_cm,  # 3% rule
    Inside_bark_radius_cm = (DBH_cm - 2 * Bark_thickness_cm) / 2
  )

# Check
head(data.all %>% filter(Stage == "adult") %>% 
       select(Tree, diameter, DBH_cm, Inside_bark_radius_cm))

# ============================================
# BAI CALCULATION - FORWARD METHOD
# ============================================

# Set path to raw data
raw_path <- "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Raw/CDendro-output/Corrected"

# Initialize empty list (easier than dataframe)
bai_list <- list()

for (plot_name in names(chronologies)) {
  raw_file <- file.path(raw_path, paste0(plot_name, "c.rwl"))
  
  if (file.exists(raw_file)) {
    rwl <- read.rwl(raw_file)
    plot_trees <- data.all %>%
      filter(Stage == "adult", Plot == gsub("_", ".", plot_name))
    
    for (i in 1:nrow(plot_trees)) {
      tree_id <- plot_trees$Tree[i]
      
      if (tree_id %in% colnames(rwl)) {
        ring_data <- rwl[, tree_id]
        valid_rows <- !is.na(ring_data)
        ring_widths <- ring_data[valid_rows] / 10
        years <- as.numeric(rownames(rwl)[valid_rows])
        
        if (length(ring_widths) > 10) {
          radius <- cumsum(ring_widths)
          n_rings <- length(radius)
          bai <- pi * (radius^2 - c(0, radius[-n_rings])^2)
          recent_years <- years >= 2004
          
          if (sum(recent_years) > 5) {
            recent_bai <- bai[recent_years]
            recent_yrs <- years[recent_years]
            trend_lm <- lm(recent_bai ~ recent_yrs)
            
            # Add to list
            bai_list[[tree_id]] <- data.frame(
              Tree = tree_id,
              Plot = gsub("_", ".", plot_name),
              Mean_recent_BAI = mean(recent_bai),
              BAI_trend = coef(trend_lm)[2],
              n_years = sum(recent_years)
            )
          }
        }
      }
    }
    cat(plot_name, "done\n")
  }
}

# Convert list to dataframe
bai_results <- bind_rows(bai_list)


# ============================================
# STEP 3: Merge with Physiology
# ============================================

bai_results <- bind_rows(bai_list)

# Merge with physiology
bai_results <- bai_results %>%
  left_join(
    data.all %>% 
      filter(Stage == "adult") %>%
      select(Tree, V_cmax, J_max, Slope, Bin, Latitude),
    by = "Tree"
  )

# Check
print(paste("Total trees:", nrow(bai_results)))
head(bai_results)
nrow(bai_results)  # Should be ~75 trees

# ============================================
# STEP 4: Analysis
# ============================================

# Does Vcmax predict recent BAI?
cor.test(bai_results$V_cmax, bai_results$Mean_recent_BAI)

# Does Vcmax predict growth trend?
cor.test(bai_results$V_cmax, bai_results$BAI_trend)

# By slope
bai_results %>%
  group_by(Slope) %>%
  summarize(
    Mean_BAI = mean(Mean_recent_BAI, na.rm = TRUE),
    Mean_trend = mean(BAI_trend, na.rm = TRUE),
    n = n()
  )

# Models
lm_bai <- lm(Mean_recent_BAI ~ V_cmax + Slope + Latitude, data = bai_results)
summary(lm_bai)

lm_trend <- lm(BAI_trend ~ V_cmax + Slope + Latitude, data = bai_results)
summary(lm_trend)

# ============================================
# STEP 5: Visualizations
# ============================================

# BAI vs Vcmax
ggplot(bai_results, aes(V_cmax, Mean_recent_BAI, color = Slope)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  labs(x = "Vcmax (μmol m⁻² s⁻¹)", 
       y = "Mean Recent BAI (cm²/yr)",
       title = "Does high Vcmax predict better growth?")

# Growth trend vs Vcmax
ggplot(bai_results, aes(V_cmax, BAI_trend, color = Slope)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Vcmax (μmol m⁻² s⁻¹)", 
       y = "BAI Trend (cm²/yr²)",
       title = "Does high Vcmax predict growth trends?")


#### Precipitation sensitifivity 

# Correlation with Vcmax
cor.test(phys_growth$Mean_Vcmax, phys_growth$Spring_precip_cor)

# Precipitation sensitivity model
lm_precip <- lm(Spring_precip_cor ~ Mean_Vcmax + Slope + Latitude, 
                data = phys_growth)
summary(lm_precip)

# Visualization
ggplot(phys_growth, aes(Mean_Vcmax, Spring_precip_cor, color = Slope)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  labs(x = "Mean Vcmax", 
       y = "Spring Precipitation Sensitivity",
       title = "Does water limitation explain topographic effects?")

# Compare temp vs precip sensitivity by slope
phys_growth %>%
  group_by(Slope) %>%
  summarize(
    Mean_temp_sens = mean(Summer_temp_cor),
    Mean_precip_sens = mean(Spring_precip_cor)
  )

# Are temp and precip sensitivities correlated?
cor.test(phys_growth$Summer_temp_cor, phys_growth$Spring_precip_cor)


###Seasonal patterns

# === 1. Collect seasonal correlation data ===
seasonal_patterns <- data.frame()

for (plot_name in names(chronologies)) {
  crn <- chronologies[[plot_name]]
  
  temp_plot <- temp_hist %>% 
    filter(Plot == plot_name) %>%
    select(year, month, temp = ssp245_tg_mean_p50)
  
  dcc_temp <- dcc(crn, temp_plot, 
                  selection = -6:9, 
                  method = "correlation",
                  dynamic = "static",
                  var_names = "temp")
  
  temp_coefs <- dcc_temp$coef %>%
    mutate(Plot = plot_name)
  
  seasonal_patterns <- rbind(seasonal_patterns, temp_coefs)
  cat(plot_name, "done\n")
}

# === 2. Merge with plot characteristics (fix Plot name format) ===
seasonal_patterns$Plot_fixed <- gsub("_", ".", seasonal_patterns$Plot)

seasonal_patterns <- seasonal_patterns %>%
  left_join(phys_growth %>% select(Plot, Bin, Latitude, Slope, Mean_Vcmax), 
            by = c("Plot_fixed" = "Plot"))

# === 3. Order months properly ===
seasonal_patterns$month_order <- factor(seasonal_patterns$month,
                                        levels = c("Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
                                                   "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP"))

# === 4. Plot by SLOPE ===
seasonal_avg_slope <- seasonal_patterns %>%
  group_by(Slope, month_order) %>%
  summarize(Mean_coef = mean(coef), .groups = "drop")

ggplot(seasonal_avg_slope, aes(x = month_order, y = Mean_coef, color = Slope, group = Slope)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Month (lowercase=prev year, CAPS=current)", 
       y = "Mean Correlation",
       title = "Seasonal temperature sensitivity by slope") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# === 5. Plot by LATITUDE BIN ===
seasonal_avg_bin <- seasonal_patterns %>%
  group_by(Bin, month_order) %>%
  summarize(Mean_coef = mean(coef), .groups = "drop")

ggplot(seasonal_avg_bin, aes(x = month_order, y = Mean_coef, color = factor(Bin), group = Bin)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Month (lowercase=prev year, CAPS=current)", 
       y = "Mean Correlation",
       title = "Seasonal temperature sensitivity by latitude",
       color = "Latitude Bin") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



### Moving window and physiology trend:

# Calculate temporal trend for each plot
window_trends <- moving_window_results %>%
  group_by(Plot) %>%
  summarize(
    Temporal_trend = coef(lm(Summer_cor ~ Year))[2],
    Early_sensitivity = mean(Summer_cor[Year <= 1980], na.rm = TRUE),
    Late_sensitivity = mean(Summer_cor[Year >= 2000], na.rm = TRUE),
    Change = Late_sensitivity - Early_sensitivity,
    .groups = "drop"
  )

# Fix plot names and merge
window_trends$Plot_fixed <- gsub("_", ".", window_trends$Plot)

window_trends <- window_trends %>%
  left_join(phys_growth %>% select(Plot, Mean_Vcmax, Mean_Jmax, Slope, Bin),
            by = c("Plot_fixed" = "Plot"))

# Check
head(window_trends)
sum(is.na(window_trends$Slope))  # Should be 0

# Analyses
cor.test(window_trends$Mean_Vcmax, window_trends$Temporal_trend)

window_trends %>%
  group_by(Slope) %>%
  summarize(
    Mean_trend = mean(Temporal_trend),
    Mean_change = mean(Change),
    n = n()
  )

lm_temporal <- lm(Temporal_trend ~ Mean_Vcmax + Slope, data = window_trends)
summary(lm_temporal)

# Visualization
ggplot(window_trends, aes(Mean_Vcmax, Temporal_trend, color = Slope)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Mean Vcmax (μmol m⁻² s⁻¹)",
       y = "Temporal Trend (Δr/year)",
       title = "Does high Vcmax predict faster adaptation?")



####### Summary stat table:
# Create comprehensive summary table
summary_stats <- data.frame(
  Variable = c(
    "Sample Size",
    "Adults", "Saplings", "Plots", "Latitude bins",
    "",
    "Vcmax (μmol m⁻² s⁻¹)",
    "  Gentle slopes (A)", "  Steep slopes (B)", "  Very steep (C)",
    "Jmax (μmol m⁻² s⁻¹)",
    "  Gentle slopes (A)", "  Steep slopes (B)", "  Very steep (C)",
    "JV Ratio",
    "  Gentle slopes (A)", "  Steep slopes (B)", "  Very steep (C)",
    "",
    "Summer Temp Sensitivity (r)",
    "  Gentle slopes (A)", "  Steep slopes (B)", "  Very steep (C)",
    "Spring Precip Sensitivity (r)",
    "",
    "Temporal Trend (Δr/year)",
    "  Gentle slopes (A)", "  Steep slopes (B)", "  Very steep (C)"
  ),
  Value = c(
    "150 trees (75+75), 15 plots, 5 bins",
    "", "", "", "",
    "",
    paste0(round(mean(data.all$V_cmax[data.all$Stage=="adult"], na.rm=T), 1), " ± ", 
           round(sd(data.all$V_cmax[data.all$Stage=="adult"], na.rm=T), 1)),
    "", "", "",
    paste0(round(mean(data.all$J_max[data.all$Stage=="adult"], na.rm=T), 1), " ± ",
           round(sd(data.all$J_max[data.all$Stage=="adult"], na.rm=T), 1)),
    "", "", "",
    paste0(round(mean(data.all$JV_ratio[data.all$Stage=="adult"], na.rm=T), 2), " ± ",
           round(sd(data.all$JV_ratio[data.all$Stage=="adult"], na.rm=T), 2)),
    "", "", "",
    "",
    paste0(round(mean(phys_growth$Summer_temp_cor), 2), " ± ",
           round(sd(phys_growth$Summer_temp_cor), 2)),
    "", "", "",
    paste0(round(mean(phys_growth$Spring_precip_cor), 2), " ± ",
           round(sd(phys_growth$Spring_precip_cor), 2)),
    "",
    paste0("+", round(mean(window_trends$Temporal_trend)*1000, 2), " ± ",
           round(sd(window_trends$Temporal_trend)*1000, 2), " ×10⁻³"),
    "", "", ""
  )
)

# Fill in slope-specific values
# I'll create a cleaner version...



######## DCC for multiple variables
# ============================================
# PHASE 2 EXPANSION: Multi-variable Climate Analysis
# ============================================

library(tidyverse)
library(treeclim)
library(dplR)

# ============================================
# SETUP: Load data
# ============================================

clim.files <- list.files("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/ClimateData/Monthly", 
                         pattern = "\\.csv", full.names = TRUE)

monthly.clim <- lapply(clim.files, FUN = read_csv)
names(monthly.clim) <- substr(basename(clim.files), 1, nchar(basename(clim.files)) - 4)

tree.rings <- lapply(plot.data$Plot, FUN = function(p) { 
  read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/TreeCores-Detrended/", 
                  p, "-Detrended-ModNegExp.csv"))
})
names(tree.rings) <- plot.data$Plot

# Variable definitions
var_to_file <- c(
  Tmean = "Tmean",
  Tmin = "Tmin", 
  Tmax = "Tmax",
  GrowingDays = "CumulatDaysAbove0",
  ColdestDay = "ColdestDay",
  ColdDays15 = "Tmin-15",
  ColdDays25 = "Tmin-25",
  HottestDay = "HottestDay",
  HeatDays32 = "Tmax32",
  TotPrecip = "TotPrec",
  WetDays10 = "WetDays10mm",
  MaxPrecip1d = "MaxTotPrec1"
)

clim.data <- lapply(var_to_file, function(fname) {
  monthly.clim[[fname]]
})
names(clim.data) <- names(var_to_file)

# ============================================
# ANALYSIS 1: CORRELATION METHOD (All 12 variables)
# ============================================
# ============================================
# SIMPLIFIED APPROACH: Correlation only, then targeted response
# ============================================

# First, get correlation working with shorter window
cat("\n=== STEP 1: Correlation analysis ===\n")

# ============================================
# SIMPLIFIED: Test each variable separately
# ============================================

climate_sensitivity_individual <- data.frame()

for (var_name in names(var_to_file)) {
  cat("\n--- Testing", var_name, "---\n")
  
  for (plot_name in plot.data$Plot) {
    cat(plot_name, "...")
    
    tryCatch({
      # Prepare data
      rings <- tree.rings[[plot_name]] %>% column_to_rownames("year")
      
      # Get ONLY this climate variable
      clim_single <- clim.data[[var_name]] %>%
        rename(clim.var = matches("^ssp245.*_p50$")) %>%
        select(Plot, year, month, clim.var) %>%
        filter(Plot == plot_name) %>%
        select(-Plot) %>%
        arrange(year, month) %>%
        rename(!!var_name := clim.var) %>%
        as.data.frame()
      
      # Run dcc for this one variable
      dcc_single <- dcc(dplR::chron(rings), clim_single,
                        method = "correlation",
                        selection = -6:9,  # Full window works with 1 variable
                        var_names = var_name)
      
      # Extract summer correlation
      coef_df <- dcc_single$coef
      summer_rows <- grep("curr\\.(jun|jul|aug)", 
                          rownames(coef_df), ignore.case = TRUE)
      
      if (length(summer_rows) > 0) {
        summer_cor <- mean(coef_df[summer_rows, "coef"], na.rm = TRUE)
        
        climate_sensitivity_individual <- rbind(
          climate_sensitivity_individual,
          data.frame(
            Plot = plot_name,
            Variable = var_name,
            Summer_cor = summer_cor
          )
        )
      }
      
      cat(" ✓")
      
    }, error = function(e) {
      cat(" FAILED:", e$message)
    })
  }
  cat("\n")
}

# Fix plot names and reshape
climate_sensitivity_individual$Plot <- gsub("_", ".", climate_sensitivity_individual$Plot)

climate_wide <- climate_sensitivity_individual %>%
  pivot_wider(names_from = Variable, 
              values_from = Summer_cor,
              names_prefix = "Sens_")

# Merge with physiology
phys_climate_full <- phys_growth %>%
  left_join(climate_wide, by = "Plot")

# Compare which variable best predicts Vcmax
corr_results <- data.frame(
  Variable = names(var_to_file),
  Cor_with_Vcmax = sapply(names(var_to_file), function(v) {
    col_name <- paste0("Sens_", v)
    if (col_name %in% names(phys_climate_full)) {
      cor(phys_climate_full$Mean_Vcmax, 
          phys_climate_full[[col_name]], 
          use = "complete.obs")
    } else NA
  }),
  p_value = sapply(names(var_to_file), function(v) {
    col_name <- paste0("Sens_", v)
    if (col_name %in% names(phys_climate_full)) {
      test <- cor.test(phys_climate_full$Mean_Vcmax, 
                       phys_climate_full[[col_name]])
      test$p.value
    } else NA
  })
) %>%
  arrange(desc(abs(Cor_with_Vcmax)))

cat("\n=== RESULTS ===\n")
print(corr_results)

# Visualize
ggplot(corr_results, 
       aes(x = reorder(Variable, Cor_with_Vcmax), 
           y = Cor_with_Vcmax)) +
  geom_col(aes(fill = p_value < 0.05)) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("gray70", "steelblue"),
                    name = "p < 0.05") +
  labs(x = NULL, 
       y = "Correlation: Climate Sensitivity → Vcmax",
       title = "Which climate variable best predicts photosynthetic capacity?",
       subtitle = "Climate sensitivity (summer growth response) vs mean Vcmax by plot") +
  theme_minimal()














##########################################OLD
##### Read in data
plot.data <- read.csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP-PlotInfo.csv") %>% 
  mutate(Plot = gsub("\\.", "_", Plot))

coords <- plot.data %>% 
  select(Plot, Longitude, Latitude) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

tree.rings <- lapply(plot.data$Plot, FUN = function(p) { 
  read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/TreeCores-Detrended/", p, "-Detrended-ModNegExp.csv"))
} )
names(tree.rings) <- plot.data$Plot
tree.rings

tree.rings.raw <- lapply(plot.data$Plot, FUN = function(p) { 
  read.rwl(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Raw/CDendro-output/Corrected/" , p, "c.rwl"), format="tucson")
} )
tree.rings.raw

plot(tree.rings.raw[[4]])

##### Check crossdating:
rwl20.2 <- corr.rwl.seg(tree.rings.raw[[4]], seg.length = 30, bin.floor = 10)
?corr.rwl.seg




##### Functions


##### Analysis
t <- detrend.series(tree.rings.raw[[1]][,2], verbose = T)

cv_per_tree <- apply(tree.rings.raw[[1]], 2, function(x) {
  x <- na.omit(x)
  sd(x) / mean(x)
})
cv_per_tree <- apply(tree.rings.raw[[2]], 2, function(x) {
  x <- na.omit(x)
  sd(x) / mean(x)
})




# Load necessary packages
library(dplR)
library(tidyverse)

# Set path to folder containing .rwl files
data_dir <- "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Raw/CDendro-output/Corrected"  # <-- Replace with your folder path
rwl_files <- list.files(data_dir, pattern = "\\.rwl$", full.names = TRUE)

# Initialize output list
results_list <- list()

# Loop through each site file
for (file in rwl_files) {
  
  # Extract site name (e.g., "Site01")
  site_name <- tools::file_path_sans_ext(basename(file))
  cat("Processing", site_name, "...\n")
  
  # Load RWL data
  rwl <- read.rwl(file)
  
  # ---- 1. CV per tree ----
  cv_per_tree <- apply(rwl, 2, function(x) {
    x <- na.omit(x)
    if(length(x) > 1) sd(x) / mean(x) else NA
  })
  
  # ---- 3. Interseries correlation and EPS ----
  rwi <- detrend(rwl, method = "ModNegExp")
  crn <- chron(rwi, prewhiten = TRUE)
  
  # Get Rbar and EPS
  stats <- rwi.stats(rwi)
  eps_vals <- rwi.stats.running(rwi, window.length = 30)
  ?rwi.stats
  # Get average Rbar and EPS over entire period
  mean_rbar <- mean(stats$rbar.bt, na.rm = TRUE)
  mean_eps <- mean(eps_vals$eps, na.rm = TRUE)
  
  # ---- 4. Summary ring width stats ----
  rw_summary <- apply(rwl, 2, function(x) {
    x <- na.omit(x)
    c(mean = mean(x), sd = sd(x), min = min(x), max = max(x))
  }) %>% t() %>% as.data.frame()
  
  rw_summary$tree_id <- rownames(rw_summary)
  rw_summary$cv <- cv_per_tree[rownames(rw_summary)]
  rw_summary$site <- site_name
  
  # ---- 5. Save results ----
  site_summary <- data.frame(
    site = site_name,
    n_trees = ncol(rwl),
    n_years = nrow(rwl),
    mean_cv = mean(cv_per_tree, na.rm = TRUE),
    mean_rbar = mean_rbar,
    mean_eps = mean_eps
  )
  
  results_list[[site_name]] <- list(
    site_summary = site_summary,
    tree_summary = rw_summary,
    chron = crn
  )
  
  # ---- 6. Plot raw data and chronology ----
  pdf(paste0(data_dir, "/", site_name, "_plots.pdf"))
  plot(rwl, plot.type = "spag", main = paste("Raw Ring Widths -", site_name))
  plot(crn, main = paste("Chronology -", site_name))
  dev.off()
}

# ---- 7. Combine and export summary results ----
site_summaries <- bind_rows(lapply(results_list, `[[`, "site_summary"))
tree_summaries <- bind_rows(lapply(results_list, `[[`, "tree_summary"))

# Export as CSVs
#write.csv(site_summaries, file = file.path(data_dir, "site_diagnostics_summary.csv"), row.names = FALSE)
#write.csv(tree_summaries, file = file.path(data_dir, "tree_diagnostics_summary.csv"), row.names = FALSE)


#### Test differences between latitudes
tree_summaries <- tree_summaries %>% 
  mutate(bin = as.numeric(substr(site, 2, 4)) / 10)
colnames(tree_summaries)

anova_model <- kruskal.test(mean ~ as.factor(bin), data = tree_summaries)
anova_model

ggplot(tree_summaries, aes(x = factor(bin), y = mean)) +
  geom_boxplot() +
  labs(x = "Latitude", y = "Coefficient of Variation (CV)") +
  theme_minimal()



###What about the age effect? Do older trees show diff climate sens indep. of ther Vcmax?
# Calculate plot-level mean age
plot_age <- data.all %>%
  filter(Stage == "adult") %>%
  group_by(Plot) %>%
  summarize(Mean_age = mean(age, na.rm=T))

phys_growth <- phys_growth %>%
  left_join(plot_age, by = "Plot")

# Test
lm_age <- lm(Summer_temp_cor ~ Mean_Vcmax + Slope + Latitude + Mean_age, 
             data = phys_growth)
summary(lm_age)

