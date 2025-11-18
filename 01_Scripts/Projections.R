library(tidyverse)

# Load monthly mean temperature
temp_monthly <- read.csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/ClimateData/Monthly/Tmean.csv")

# Calculate climate changes
climate_change <- temp_monthly %>%
  # Calculate summer temp (Jun-Aug) and annual mean
  group_by(Plot, year) %>%
  summarize(
    # Summer (matches Phase 2 analysis)
    Summer_126 = mean(ssp126_tg_mean_p50[month %in% 6:8]),
    Summer_245 = mean(ssp245_tg_mean_p50[month %in% 6:8]),
    Summer_585 = mean(ssp585_tg_mean_p50[month %in% 6:8]),
    # Annual mean
    Annual_126 = mean(ssp126_tg_mean_p50),
    Annual_245 = mean(ssp245_tg_mean_p50),
    Annual_585 = mean(ssp585_tg_mean_p50),
    .groups = "drop"
  ) %>%
  # Define time periods
  mutate(
    Period = case_when(
      year >= 1950 & year <= 2024 ~ "Historical",
      year >= 2041 & year <= 2060 ~ "2050s",
      year >= 2061 & year <= 2080 ~ "2070s",
      year >= 2081 & year <= 2100 ~ "2090s",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Period))

# Calculate baseline and changes
climate_summary <- climate_change %>%
  pivot_longer(cols = contains("_"), 
               names_to = c("Season", "Scenario"), 
               names_sep = "_",
               values_to = "Temperature") %>%
  group_by(Plot, Period, Season, Scenario) %>%
  summarize(Mean_temp = mean(Temperature), .groups = "drop") %>%
  pivot_wider(names_from = Period, values_from = Mean_temp) %>%
  mutate(
    Delta_2050s = `2050s` - Historical,
    Delta_2070s = `2070s` - Historical,
    Delta_2090s = `2090s` - Historical
  )

# Check results
head(climate_summary)

# Quick summary
climate_summary %>%
  filter(Season == "Summer") %>%
  group_by(Scenario) %>%
  summarize(
    Mean_warming_2070s = mean(Delta_2070s),
    Range_2070s = paste(round(range(Delta_2070s), 2), collapse = " to ")
  )


######1: climate departure:
# Fix Plot names for merging
climate_summary <- climate_summary %>%
  mutate(Plot = gsub("_", ".", Plot))

# Calculate current climate envelope
current_envelope <- climate_summary %>%
  filter(Season == "Summer", Scenario == "245") %>%
  summarize(
    Current_min = min(Historical),
    Current_max = max(Historical),
    Current_range = Current_max - Current_min
  )

print(current_envelope)

# Identify sites exceeding envelope
climate_departure <- climate_summary %>%
  filter(Season == "Summer", Scenario == "245") %>%
  mutate(
    Exceeds_2070s = `2070s` > current_envelope$Current_max,
    Departure_2070s = `2070s` - current_envelope$Current_max
  ) %>%
  select(Plot, Historical, `2070s`, Delta_2070s, Exceeds_2070s, Departure_2070s)

# Summary
table(climate_departure$Exceeds_2070s)
climate_departure %>% arrange(desc(Departure_2070s))

########2: Combine with climate vulnerability
# Merge with vulnerability data
future_risk <- phys_growth %>%
  left_join(climate_departure, by = "Plot") %>%
  mutate(
    # Combined risk metric
    Climate_stress = Delta_2070s * -1 * Summer_temp_cor,  # Negative corr × warming = stress
    
    # Risk categories
    High_warming = Delta_2070s > median(Delta_2070s),
    Negative_response = Summer_temp_cor < 0,
    
    Double_jeopardy = High_warming & Negative_response
  ) %>%
  select(Plot, Bin, Latitude, Slope, 
         Mean_Vcmax, Summer_temp_cor, 
         Delta_2070s, Departure_2070s,
         Climate_stress, Double_jeopardy)

# Summary
table(future_risk$Double_jeopardy)

# Most at-risk sites
future_risk %>%
  arrange(desc(Climate_stress)) %>%
  select(Plot, Slope, Summer_temp_cor, Delta_2070s, Climate_stress, Double_jeopardy)

ggplot(future_risk, aes(Delta_2070s, Summer_temp_cor, 
                        color = Double_jeopardy, shape = Slope)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 4) +
  geom_text(aes(label = Plot), hjust = -0.2, size = 3) +
  scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "red")) +
  labs(x = "Projected Warming by 2070s (°C)", 
       y = "Current Temperature Sensitivity",
       title = "Future Climate Risk") +
  theme_bw()


# Plot 2: Stress by topography
ggplot(future_risk, aes(Slope, Climate_stress, fill = Slope)) +
  geom_boxplot() +
  geom_point(size = 3, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Climate Stress Index", 
       title = "Future Climate Stress by Topography") +
  theme_bw()


########3: Other time periods, more ssp
# Expand climate departure for ALL scenarios and periods
climate_departure_all <- climate_summary %>%
  filter(Season == "Summer") %>%
  pivot_longer(cols = c(`2050s`, `2070s`, `2090s`),
               names_to = "Future_period",
               values_to = "Future_temp") %>%
  mutate(
    Delta = Future_temp - Historical,
    Exceeds = Future_temp > current_envelope$Current_max,
    Departure = Future_temp - current_envelope$Current_max
  ) %>%
  select(Plot, Scenario, Future_period, Historical, Future_temp, Delta, Exceeds, Departure)

# Summary table: How many sites exceed envelope?
climate_departure_all %>%
  group_by(Scenario, Future_period) %>%
  summarize(
    n_exceeds = sum(Exceeds),
    Mean_warming = mean(Delta),
    Max_departure = max(Departure),
    .groups = "drop"
  )

# Integrate with vulnerability for all combinations
future_risk_all <- phys_growth %>%
  select(Plot, Slope, Latitude, Mean_Vcmax, Summer_temp_cor) %>%
  left_join(climate_departure_all, by = "Plot", relationship = "many-to-many") %>%
  mutate(
    Climate_stress = Delta * -1 * Summer_temp_cor,
    Double_jeopardy = Delta > median(Delta) & Summer_temp_cor < 0
  )

# Count double jeopardy sites by scenario/period
future_risk_all %>%
  group_by(Scenario, Future_period) %>%
  summarize(
    n_double_jeopardy = sum(Double_jeopardy),
    Mean_stress = mean(Climate_stress),
    .groups = "drop"
  )
