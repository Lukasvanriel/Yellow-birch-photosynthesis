### Script to funnel all the different raw data files into one coherent database
### Author: Lukas Van Riel
### Date: 25/09/2024

### Load packages
library(tidyverse)
library(conflicted)

#### 0 - Basic permutations #### 
## Diameters
data.diam <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Raw/Tree_diameters-wide.csv")
diam <- data.diam %>% 
  pivot_longer(cols = -Plot, values_to = "diameter") %>% 
  mutate(tree = paste(Plot, ifelse(substring(name, 1, 1) == "a", substring(name, nchar(name), nchar(name)),
                                   as.character(as.numeric(substring(name, nchar(name), nchar(name))) + 5)),
                      sep="-"),
         diameter = 10 * diameter) %>% 
  relocate(tree, .before = diameter) %>% 
  select(tree, diameter)
  
write_csv(diam, "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Tree_diameters.csv")

## Tree community
tree.comm.long <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Community-dbh.csv")

tree.comm <- tree.comm.long %>% 
  mutate(Diameter = 10 * Diameter) %>% 
  mutate(A = (Diameter/2)**2 * pi) %>% 
  group_by(Plot, Species) %>% 
  summarise(Atot = sum(A)) %>% 
  pivot_wider(names_from = Species, values_from = Atot, values_fill = 0.0) %>% 
  ungroup()

write_csv(tree.comm, "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Community-dbh.csv")

#### 1 - Leaf parameters #### 
#Load txt data file:
data.leaf.raw <- read_table("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Raw/Lukas-LeafData/Leaf-Scans.txt")

data.leaf <- data.leaf.raw %>% 
  mutate(LeafName = as.character(WinFOLIA),
         LeafArea = as.numeric(AvgLeafArea), 
         Width = as.numeric(AvgLeafHorizWidth), 
         Length = as.numeric(AvgLeafVertLength)) %>% 
  select(LeafName, LeafArea, Width, Length)

write_csv(data.leaf, "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Leaf-Dimensions.csv")

#### 2 - Split Racir Data into individual measurements #### 
library(readxl)
library(plantecophys)
library(racir)

#functions
extract_licor_data <- function(path){
  sheets <- excel_sheets(path)
  
  # Split into individual measurements ------
  
  data.input <- read_excel(path, sheets[1], skip = 14, col_types = "numeric") %>% 
    filter(!is.na(obs)) %>% 
    mutate(A=as.numeric(A), Ca=as.numeric(Ca), CO2_r=as.numeric(CO2_r),
           E=as.numeric(E), gtc=as.numeric(gtc), Tleaf=as.numeric(Tleaf)) %>% 
    mutate(t_diff=c(0, diff(as.numeric(elapsed)))) %>% 
    relocate(t_diff, .after = obs)
  
  # Split into the individual measurements
  
  data.breaks <- data.input %>% 
    filter(t_diff>20) %>% # There should easily be more than 20s between measurements 
    mutate(n_meas=c(obs[1], diff(as.numeric(obs)))) %>% 
    relocate(n_meas, .after = obs)
  
  breaks <- c(data.breaks$n_meas, nrow(data.input)-sum(data.breaks$n_meas))
  
  test <- lapply(1:length(breaks), function(i) rep(as.character(i), breaks[i])) %>% 
    unlist()
  
  data.split <- split(data.input, test)
  
  #Only keep real measurements, which take at least 5 minutes (150 datapoints)
  measurements <- data.split[which(sapply(data.split, function(x) nrow(x) > 150))]
  
  #Order the measurements correctly 
  measurements <-  measurements[order(as.numeric(names( measurements)))]
  
  #####Analysis
  # Manually inspect the number of measurements of each data frame
  print(sapply(measurements, nrow))
  
  measurements
}

racir_to_csv <- function(racir_file_path) {
  racirs <- extract_licor_data(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Raw/ACiResponse/",
                                      racir_file_path, ".xlsx"))
  
  for (i in 1:length(racirs)) {write_csv(racirs[[i]], paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas_measurements/",
                                                             substring(racir_file_path, nchar(racir_file_path)-5, nchar(racir_file_path)), "-", as.character(i), ".csv"))}
}

#Create list of plotnames
plotnames <- c("2024-07-29-0958_P465.1", "2024-07-30-0928_P465.2", "2024-07-30-1342_P465.3", "2024-08-05-1224_P475.1-Restored",
               "2024-08-05-1733_2", "2024-08-06-0928_P480.1", "2024-08-07-0933_P480.3", "2024-08-07-1316_P480.2",
               "2024-08-08-0922_P475.2", "2024-08-08-1346_P475.3", "2024-08-13-1157_P470.1", "2024-08-14-0912_P470.2",
               "2024-08-14-1311_P470.3", "2024-08-15-0928_P460.1", "2024-08-15-1323_P460.2", "2024-08-16-0941_P460.3")

sapply(plotnames, racir_to_csv)

#### 3 - Aci curve data analysis #### 
#functions

racir_bulk <- function(data, empty, 
                       default.range = as.logical(readline("Use CO2 range [50,1050]? (T/F) ")),
                       A.limits = as.logical(readline("Are A limits needed? (T/F) "))){
  # Fix some weird behaviour with tibble
  data <- as.data.frame(data)
  empty <- as.data.frame(empty)
  
  # Start function
  if(default.range) {
    cmin <- 50
    cmax <- 1050
    racircalcheck(data = data, mincut = cmin, maxcut = cmax)
  } else {
    
    proceed <- F
    while(! proceed) {
      racircalcheck(data = data)
      
      cmin <- as.numeric(readline("What is the low point cutoff? "))
      cmax <- as.numeric(readline("What is the high point cutoff? "))
      
      racircalcheck(data = data, mincut = cmin, maxcut = cmax)
      
      proceed <- as.logical(readline("Proceed with these cutoffs? (T/F) "))
    }
  }
  
  if(A.limits){
    Amax <- as.numeric(readline("Provide the upper limit for A: "))
    Amin <- as.numeric(readline("Provide the lower limit for A: "))
    data <- data %>% filter(A < Amax, A>Amin)
  }
  
  corrected <- racircal(data = data, caldata = empty, mincut = cmin, maxcut = cmax)
  corrected <- corrected %>% 
    filter(Cicor<1500, Cicor>0)
  
  #Now add fitaci
  fit <- tryCatch(
    expr = {
      fitaci(corrected, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                        PPFD = "PARi", Rd ="Rd"), Tcorrect=T)
    },
    error = function(e){
      message('Normal fitting is not working! Trying the bilinear fit method now. ')
      fitaci(corrected, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                        PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='bilinear')
    },
    finally = {
      message('All done, quitting.')
    }
  )   
  
  plot(fit)
  
  print(fit$pars)
  fit
}

which_tree <- function(string) {
  if (substr(string, 1,1) == "a") {
    as.numeric(substr(string, nchar(string), nchar(string)))
  } else(as.numeric(substr(string, nchar(string), nchar(string))) + 5)
}

Extract_curve_parameters <- function(plot, tree){
  data <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                          plot, "-", as.character(info[info$Plot==plot,tree]), ".csv"))
  empty <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                           plot, "-", as.character(info[info$Plot==plot,"empty"]), ".csv"))
  # Run the function
  fit <- racir_bulk(data, empty, default.range = T, A.limits = F)
  
  # Get rid of outliers (Q1/Q3 +- 1.5 x IQR)
  A_range <- c(quantile(fit$df$Ameas)[[2]] - IQR(fit$df$Ameas) *1.5, quantile(fit$df$Ameas)[[4]] + IQR(fit$df$Ameas) *1.5)
  
  fit$df <- fit$df[fit$df$Ameas >= A_range[1] & fit$df$Ameas <= A_range[2],]
  
  # Plot results
  pdf(file=paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Figures/",
                      paste(plot,which_tree(tree), sep = "-"), ".pdf"))
  plot(fit)
  text(x = par("usr")[1], y = par("usr")[4], 
       labels = paste0(paste("Vc,max =", round(fit$pars["Vcmax", "Estimate"],1)), "\n", paste(" Jmax   =", round(fit$pars["Jmax", "Estimate"], 1)),
                       "\n", paste("Method:", fit$fitmethod)), 
       adj = c(-0.18, 1.3), cex = 0.9)
  dev.off()
  # Return results
  return(c(paste(plot,which_tree(tree), sep = "-"), unlist(fit$pars), fit$fitmethod))
}

#
info <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Measurements-information.csv")

parameters_df <- data.frame(ID=rep("",150),
                            Vcmax=rep(0,150),
                            Jmax=rep(0,150),
                            Rd=rep(0,150),
                            SE_Vcmax=rep(0,150),
                            SE_Jmax=rep(0,150),
                            SE_Rd=rep(0,150),
                            Method=rep("",150)) 

sapply(info$Plot, FUN = function(p) {sapply(colnames(info)[3:ncol(info)], 
                                            FUN = function(t) {parameters_df[min(which(parameters_df$Vcmax == 0)), ] <<- Extract_curve_parameters(p, t)})})

write_csv(parameters_df, "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/parameters.csv")

hist(as.numeric(parameters_df$Vcmax))
hist(as.numeric(parameters_df$Jmax))

table(parameters_df$Method)

parameters_df[min(which(parameters_df$Vcmax == 0)), ] <- Extract_curve_parameters(p, t)


#### 4 - Check out data #### 
data_raw <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Parameters_fitted.csv")
data_plot <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Raw/Plot-info.csv")


data <- data_raw %>% 
  select(-Num) %>% 
  mutate(full=tree) %>% 
  separate_wider_delim(tree, "-", names=c("Plot", "Tree")) %>% 
  filter(Tree != "e") %>% # drop the extra tree for now
  left_join(data_plot, by="Plot") %>% 
  separate_wider_delim(Plot, ".", names=c("Bin", "Plot")) %>% 
  mutate(Bin=factor(paste0(substring(Bin, 2,3), ".", substring(Bin, 4,4)))) %>% 
  mutate(sapling=as.numeric(Tree)%/%6) %>% 
  mutate(Class=factor(ifelse(sapling==1, "sapling", "adult"))) %>% 
  relocate(Class, .before = V_cmax)

ggplot(data) +
  geom_point(aes(x=V_cmax, y=J_max, colour = Class))

ggplot(data) +
  geom_point(aes(x=Latitude, y=V_cmax, colour = Class))

ggplot(data) +
  geom_point(aes(x=Latitude, y=J_max, colour = Class))

ggplot(data) +
  geom_boxplot(aes(y=J_max, x=Bin, colour = Class))

ggplot(data) +
  geom_boxplot(aes(y=V_cmax, x=Bin, colour = Class))

# Run some anova's
aov_Va <- aov(V_cmax ~ Bin,
                         data = data %>% filter(Class=="adult"))
summary(aov_Va)
ggplot(data %>% filter(Class=="adult")) +
  geom_boxplot(aes(y=V_cmax, x=Bin))

aov_Vs <- aov(V_cmax ~ Bin,
              data = data %>% filter(Class=="sapling"))
summary(aov_Vs)
ggplot(data %>% filter(Class=="sapling")) +
  geom_boxplot(aes(y=Vcmax, x=Bin))

aov_Ja <- aov(J_max ~ Bin,
              data = data %>% filter(Class=="adult"))
summary(aov_Ja)
ggplot(data %>% filter(Class=="adult")) +
  geom_boxplot(aes(y=J_max, x=Bin))

aov_Js <- aov(J_max ~ Bin,
              data = data %>% filter(Class=="sapling"))
summary(aov_Js)
ggplot(data %>% filter(Class=="sapling")) +
  geom_boxplot(aes(y=J_max, x=Bin))


shapiro.test(aov_Va$residuals)
hist(aov_Va$residuals)
par(mfrow = c(1, 2)) # combine plots

# histogram
hist(aov_Va$residuals)

# QQ-plot
qqPlot(aov_Va$residuals,
       id = FALSE # id = FALSE to remove point identification
)


#### 5 - Granulometrie #### 
data_gran <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Raw/Granulometrie_YBP.csv",
                      col_types = c("i", "c", "c", rep("d", 15), "c", "c"))

# Drop unimportant columns, then average over three runs
data_split <- data_gran %>% 
  select(-c(Measurement, Date, `Ultrason-mode`, `Ultrasons-durée`)) %>% 
  filter(Name != "test1-test") %>% 
  mutate(clay=bin0to2, 
         silt= bin2to4 + bin4to8 + bin8to15 + bin15to31 + bin31to53,
         sand= bin53to106 + bin106to250 + bin250to500 + bin500to1000 + bin1000to2000) %>% 
  mutate(clay_n = clay/(clay+silt+sand),
         silt_n = silt/(clay+silt+sand),
         sand_n = sand/(clay+silt+sand)) %>%
  group_by(Name) %>% 
  summarise(across(c(clay_n, silt_n, sand_n), mean))

  
# Create new columns for sand, silt and clay
data_class <- data_split %>% 
  mutate(soil_class=ifelse(sand_n<0.50, "SiLo",
                      ifelse(sand_n<0.70, "SaLo",
                      ifelse(sand_n<0.90, "LoSa", "Sa")))) %>% 
  separate_wider_delim(Name, names = c("Plot", "Horizon"), delim = "-")

# Correct 1 measurement with double measurement: P475.1-A2
data_class$soil_class[36] <- data_class$soil_class[37]

write_csv(data_class, "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Soil_Granulometrie.csv")

#### 6 - Soil horizons description #### 
data_soil <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Raw/Soil-layers.csv",
                      col_types = c())

# Drop some of the columns (only for reference) and merge with the soil classes
soil_merge <- data_soil %>% 
  select(Plot, Horizon, pH, Munsell_color) %>% 
  left_join(select(data_class, c(Plot, Horizon, soil_class)), by=c("Plot", "Horizon")) %>% 
  relocate(soil_class, .after = Horizon) %>% 
  separate_wider_delim(Munsell_color, delim = " ", names = c("Munsell_hue", "Munsell_tosplit")) %>% 
  separate_wider_delim(Munsell_tosplit, delim = "/", names = c("Munsell_value", "Munsell_chroma"))

write_csv(soil_merge, "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Soil_Layers.csv")

#### 7 - CDendro output #### 
## Functions

basic_ring_analysis <- function(file) {
  #Create list to be returned by function
  results <- list()
  
  ## Read in the file
  if(substr(file, nchar(file)-2, nchar(file)) == "rwl") {
    ring.data.raw <- read.rwl(file, format="tucson") %>% 
      mutate(year=rownames(ring.data.raw))
  } else {if(substr(file, nchar(file)-2, nchar(file)) == "csv") {
        ring.data.raw <- read_csv(file)
      }
    }
  
  ## Calculate age of each tree
  results$age <- ring.data.raw %>%
    pivot_longer(cols = -year, names_to = "Plot", values_to = "Width") %>% 
    group_by(Plot) %>% 
    mutate(non_na_flag = cumsum(! is.na(Width)) > 0) %>% 
    filter(non_na_flag) %>%                             # Filter to include only rows after the first non-NA
    summarise(age = n()) %>% 
    pivot_wider(names_from = Plot, values_from = age)
  
  ## Calculate the average ring widths per decade
  # Add years and decades as columns
  ring_data <- ring.data.raw %>% 
    mutate(Decade = floor(year / 10) * 10)

  # Pivot longer and summarise by decade and plot
  results$dec_averages <- ring_data %>%
    pivot_longer(
      cols = -c(year, Decade), # Keep Year and Decade as identifiers
      names_to = "Plot",
      values_to = "Width") %>%
    group_by(Decade, Plot) %>%
    summarise(
      Decade_Avg = ifelse(all(is.na(Width)), NA, mean(Width, na.rm = TRUE)),
      Valid_Years = sum(!is.na(Width)),
      NA_Years = sum(is.na(Width)),
      .groups = "drop" ) # Drop grouping for cleaner output
  
  ## Get growth indices for last 10, 20 to 50 years:
  results$growth_indices <- ring.data.raw %>% 
    pivot_longer(cols = -year, names_to = "Plot", values_to = "Width") %>% 
    group_by(Plot) %>% 
    summarise(mean10 = mean(tail(Width, 10), na.rm = TRUE),
              mean20 = mean(tail(Width, 20), na.rm = TRUE),
              mean30 = mean(tail(Width, 30), na.rm = TRUE),
              mean40 = mean(tail(Width, 40), na.rm = TRUE),
              mean50 = mean(tail(Width, 50), na.rm = TRUE))
  
  ## Return ages and decadal averages as a list
  results
}

write.detrend <- function(file) {
  plotname <- substring(basename(file), 1, 6)
  
  #
  rwl <- read.rwl(file, format="tucson")
  
  p.rwl <- rwl %>% 
    mutate(year = as.numeric(rownames(rwl))) %>% 
    pivot_longer(cols = -c(year), names_to = "plot", values_to = "width") %>% 
    dplyr::filter(! is.na(width)) %>% 
    ggplot(aes(x = year, y = width, col = plot)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) + 
    labs(
      title = paste(plotname, "- Original"),
      x = "Year",
      y = "Width"
    ) + 
    theme_minimal() 
  
  pdf(paste0('/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/TreeCores-Detrended/Figures/', plotname, '-Original.pdf'))
  print(p.rwl)
  dev.off()

  # ModNegExp
  rwi <- detrend(rwl, method = "ModNegExp")
  
  rwi <- mutate(rwi, year = as.numeric(rownames(rwi)))
  
  p.rwi <- rwi %>% 
    pivot_longer(cols = -c(year), names_to = "plot", values_to = "width") %>% 
    dplyr::filter(! is.na(width)) %>% 
    ggplot(aes(x = year, y = width, col = plot)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) + 
    labs(
      title = paste(plotname, "- Detrended"),
      x = "Year",
      y = "Detrended and standardised Width"
    ) + 
    theme_minimal()
  
  pdf(paste0('/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/TreeCores-Detrended/Figures/', plotname, '-Trends.pdf'))
  print(p.rwi)
  dev.off()
  
  lapply(names(rwl), function(x){
    pdf(paste0('/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/TreeCores-Detrended/Figures/', x, '-Series.pdf'))
    plot(detrend.series(y = rwl[,x], method = "ModNegExp", verbose=TRUE))
    dev.off()
  })
  


  write_csv(rwi, paste0('/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/TreeCores-Detrended/', plotname, '-Detrended-ModNegExp.csv'))
  
  if(F){
  # Spline
  rwi <- detrend(rwl, method = "Spline")
  
  p.rwi <- rwi %>% 
    mutate(year=as.numeric(rownames(rwi))) %>% 
    pivot_longer(cols = -c(year), names_to = "plot", values_to = "width") %>% 
    dplyr::filter(! is.na(width)) %>% 
    ggplot(aes(x = year, y = width, col = plot)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) + 
    labs(
      title = paste(plotname, "- Detrended"),
      x = "Year",
      y = "Detrended and standardised Width"
    ) + 
    theme_minimal()
  
  pdf(paste0('/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/TreeCores-Detrended/Figures/', plotname, '-Detrended-Spline.pdf'))
  print(p.rwi)
  dev.off()
  
  write_csv(rwi, paste0('/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/TreeCores-Detrended/', plotname, '-Detrended-Spline.csv'))
  }
}

## Analysis
library(dplR)
test <- read.rwl(files[1], format="tucson")

files <- Sys.glob("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Raw/CDendro-output/Corrected/*.rwl")

lapply(files, write.detrend)

files.csv <- Sys.glob("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/TreeCores-Detrended/*-ModNegExp.csv")

ring_detrended <- setNames(lapply(files.csv, FUN = basic_ring_analysis),
                           substring(basename(files.csv), 1, 6))


# Combine all and write out 
age_combined <- bind_rows(lapply(ring_detrended, FUN = function(x) {
  ages <- x[["age"]]
  colnames(ages) <- c("adult1", "adult2", "adult3", "adult4", "adult5")
  return(ages)
} ) )
colnames(age_combined) <- c(colnames(age_combined)[1:5], "extra")

age_comb_long <- age_combined %>% 
  pivot_longer(cols = everything(), names_to = "Plot", values_to = "age") %>% 
  mutate(Plot = paste(rep(plot.info$Plot, each=6), rep(c(as.character(1:5), "e"), 15), sep = "-")) %>% 
  filter(! is.na(age))


# Combine the second elements of each list
detrended_combined <- bind_rows(lapply(ring_detrended, `[[`,"dec_averages"))


# Combine the second elements of each list
tails_combined <- bind_rows(lapply(ring_detrended, `[[`,"growth_indices"))


write_csv(age_comb_long, "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/ages.csv")
write_csv(detrended_combined, "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Decadal-indices.csv")
write_csv(tails_combined, "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Growth-indices-last-years.csv")

#### Combine into one file:

ages <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/ages.csv")
growth_ind <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Growth-indices-last-years.csv")

ring.an.comb <- left_join(ages, growth_ind, by = "Plot") %>% 
  rename(tree=Plot)
  
write_csv(ring.an.comb, "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/tree.ring.analysis.csv")

if(F){
## OLD

ring_basics <- setNames(lapply(files, FUN = basic_ring_analysis),
                        substring(files, nchar(files[1]) - 10, nchar(files[1]) - 5))

t <- ring_basics$P460_1$averages

# Combine all and write out 
age_combined <- bind_rows(lapply(ring_basics, FUN = function(x) {
  ages <- x[["age"]]
  colnames(ages) <- c("adult1", "adult2", "adult3", "adult4", "adult5")
  return(ages)
  } ) )
colnames(age_combined) <- c(colnames(age_combined)[1:5], "extra")
                             
# Combine the second elements of each list
averages_combined <- bind_rows(lapply(ring_basics, `[[`,"averages"))

### Check out some of the results
summary(na.omit(unlist(age_combined)))

## Group by bin
library(gtools)

averages_binned <- averages_combined %>% 
  separate_wider_delim(cols = Plot, delim = "-", names = c("Plot", "adult")) %>% 
  separate_wider_delim(cols = Plot, delim = ".", names = c("bin", "plot")) %>% 
  group_by(Decade, bin) %>% 
  summarise(m=mean(Decade_Avg, na.rm = T)) %>% 
  mutate(bin=factor(bin, levels = c("P460", "P465", "P470", "P475", "P480")))

ggplot(data=averages_binned) +
  geom_point(aes(x=Decade, y=m, colour = bin)) +
  theme_minimal()

## plot linear trends over time for each bin

ggplot(data=averages_binned, aes(x=Decade, y=m, colour = bin)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  labs(
    title = "Linear Regression Trends of m Over Time by Bin",
    x = "Decade",
    y = "Average ring width"
  ) +
  theme_minimal()

class(averages)

### Is there a way to choose which detrending method?
# The book recommends plotting the mean trends over 20y periods, overlapping by 10:

x <- ring_basics[[1]]
rd <- x$averages
names(ring_basics) <- gsub("_", ".", names(ring_basics))

all.20yrav <- lapply(1:15, FUN = function(x){
  rd <- ring_basics[[x]]$averages
  plot <- names(ring_basics)[x]
  
  dec2av.list <- lapply(1:5, FUN=function(y){
    rd1 <- filter(rd, Plot == paste0(plot,"-", y))
    print(rd1)
    tree <- rep(paste0(plot, "-",  y), nrow(rd1) - 1)
    central.year <- rd1$Decade[2:nrow(rd1)]
    y20av <- vector("numeric", length = nrow(rd1)-1)
    
    for(d in 2:nrow(rd1)) {
      y20av[d-1] <- (rd1$Decade_Avg[d-1] * rd1$Valid_Years[d-1] + rd1$Decade_Avg[d] * rd1$Valid_Years[d]) / 
        (rd1$Valid_Years[d-1] + rd1$Valid_Years[d])
    }
    
    rd20y <- data.frame(tree = tree, centr.year = central.year, rw.av = y20av)
  })
  dec2av <- bind_rows(dec2av.list)
  
  fig <- ggplot(dec2av, aes(x = centr.year, y = rw.av, col=tree)) +
    geom_point() +
    geom_smooth(method = "lm", se = F) + 
    labs(
      title = paste(plot, "- 20year ring width averages"),
      x = "Central Year",
      y = "Average ring width"
    ) + 
    theme_minimal()
    
  print(fig)
  
  dec2av
})



###
### Let's go trhough the dplr package documentation--f
## Try for 1 plot:
rwl <- read.rwl(files[15], format="tucson")

#corr.rwl.seg(rwl, seg.length = 50, bin.floor = 20)

## Get some basic information
plot(rwl, plot.type="spag")

rwl.report(rwl)
rwl.stats(rwl)
# Mean sensitivity used to be included here as well, but book said to not use it!
# At best it is as good as the stdev of time series in cases of high autocorrelation

#Can also plot some stuff if would be interesting, eg for autocorrelation:
ar1 <- data.frame(x="P480.3",y=rwl.stats(rwl)$ar1)
ggplot(ar1,aes(x,y)) + geom_boxplot(width=.2) +
  geom_jitter(width=0.1) + 
  labs(y=expression(phi[1]),x=element_blank()) +
  theme_minimal()

## Detrending:
# Pivot to longer format 
rwl_long <- rwl %>% 
  mutate(year=as.numeric(rownames(rwl))) %>% 
  pivot_longer(cols = -c(year), names_to = "plot", values_to = "width") %>% 
  filter(! is.na(width))

# Inspect the trends
ggplot(data=rwl_long, aes(x = year, y = width, col = plot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  theme_minimal() 

# Seems to be a trend indeed, so let's detrend:
rwi <- detrend(rwl, method = "ModNegExp")
class(rwi)

rwi %>% 
  mutate(year=as.numeric(rownames(rwi))) %>% 
  pivot_longer(cols = -c(year), names_to = "plot", values_to = "width") %>% 
  dplyr::filter(!is.na(width)) %>% 
  ggplot(aes(x = year, y = width, col = plot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  theme_minimal() 

!is.na(t$width)

# See what the options are
detrend.series(y = rwl[, "P480.3-2"], method = "ModNegExp", verbose=TRUE)

# There are other detrend methods, such as converting to basal area increases:
bai <- bai.out(rwl, diam = NULL)
sum(rwi$`P460.1-1`) # Not correct
#TODO find way to estimate bark width
diam = data.frame(plot = colnames(rwl),
                  d = (c(275,	333, 390, 392, 251) - 8)) #0.8 cm bark?
bai <- bai.out(rwl, diam = diam)


bai %>% 
  mutate(year=as.numeric(rownames(rwl))) %>% 
  pivot_longer(cols = -c(year), names_to = "plot", values_to = "width") %>% 
  dplyr::filter(! is.na()) %>% 
  ggplot(aes(x = year, y = width, col = plot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  theme_minimal() 

#Plot with all of them: (skip bai for now since seem to have negative BAI's unless I use crazy bark thickness...)

p.rwi <- rwi %>% 
  mutate(year=as.numeric(rownames(rwi))) %>% 
  pivot_longer(cols = -c(year), names_to = "plot", values_to = "width") %>% 
  filter(! is.na(width)) %>% 
  ggplot(aes(x = year, y = width, col = plot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  labs(
    title = paste(plotname, "- Detrended"),
    x = "Year",
    y = "Detrended and standardised Width"
  ) + 
  theme_minimal()


# Check stats of detrended: (correlation between series (both within and between tree correlations) as well as the expressed population signal, signal-to-noise ratio, and subsample signal strength for a data set)
rwi.stats(rwi, prewhiten=TRUE)

# Or other interesting one: average interseries correlation
interseries.cor(rwi, prewhiten=TRUE,
                method="spearman")


# Now let's create a mean sample chronology:
rwchron <- chron(rwi)
plot(rwchron, add.spline=TRUE, nyrs=20)

## Ok let's detrend them all, using the ModNegExp, Spline: 

### Let's see how some of the trends are now, only looking at :
files.csv <- Sys.glob("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/TreeCores-Detrended/*-ModNegExp.csv")

ring_detrended <- setNames(lapply(files.csv, FUN = basic_ring_analysis),
                        substring(files.csv, nchar(files.csv[1]) - 9, nchar(files.csv[1]) - 4))

# Combine all and write out 
age_combined <- bind_rows(lapply(ring_detrended, FUN = function(x) {
  ages <- x[["age"]]
  colnames(ages) <- c("adult1", "adult2", "adult3", "adult4", "adult5")
  return(ages)
} ) )
colnames(age_combined) <- c(colnames(age_combined)[1:5], "extra")

# Combine the second elements of each list
detrended_combined <- bind_rows(lapply(ring_detrended, `[[`,"averages"))

### Check out some of the results
summary(na.omit(unlist(detrended_combined)))

## Group by bin
library(gtools)

detrended_binned <- detrended_combined %>% 
  separate_wider_delim(cols = Plot, delim = "-", names = c("Plot", "adult")) %>% 
  separate_wider_delim(cols = Plot, delim = ".", names = c("bin", "plot")) %>% 
  group_by(Decade, bin) %>% 
  summarise(m=mean(Decade_Avg, na.rm = T)) %>% 
  mutate(bin=factor(bin, levels = c("P460", "P465", "P470", "P475", "P480")))

ggplot(data=detrended_binned) +
  geom_point(aes(x=Decade, y=m, colour = bin)) +
  theme_minimal()

## plot linear trends over time for each bin

ggplot(data=detrended_binned, aes(x=Decade, y=m, colour = bin)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  labs(
    title = "Linear Regression Trends of m Over Time by Bin",
    x = "Decade",
    y = "Average ring width"
  ) +
  theme_minimal()

# Look at standard deviation as proxy for climate importance?
testfile <- read_csv(files.csv[1])

stdev.list <- setNames(lapply(files.csv, function(x) rwl.stats(read_csv(x))$stdev),
         substring(files.csv, nchar(files[1]) - 9, nchar(files[1]) - 4))
stdev.list$P475_1 <- stdev.list$P475_1[1:5] # Drop extra measurement for now

stdev.df <- bind_rows(stdev.list) %>% 
  pivot_longer(cols = everything(), names_to = "plot", values_to = "stdev") %>% 
  separate_wider_delim(plot, delim = "_", names = c("bin", "plot")) 

ggplot(stdev.df, aes(x=as.numeric(substring(bin, 3,4))/10, y=stdev, colour = plot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  labs(
    title = "Linear Regression Trends of stdev by plot",
    x = "Bin",
    y = "St Dev averaged by plot"
  ) +
  theme_minimal()

stdev_av <- stdev.df %>% 
  group_by(bin) %>% 
  summarise(m=mean(stdev))
  
ggplot(stdev_av, aes(x=bin, y=m)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  labs(
    title = "Linear Regression Trends of stdev by plot",
    x = "Bin",
    y = "St Dev averaged by plot"
  ) +
  theme_minimal()

stdev_av_ind <- bind_rows(stdev.list) %>%
  pivot_longer(cols = everything(), names_to = "plot", values_to = "stdev") %>% 
  group_by(plot) %>% 
  summarise(m=mean(stdev)) %>% 
  separate_wider_delim(cols = plot, delim = "_", names = c("bin", "plot"))

ggplot(stdev_av_ind, aes(x=as.numeric(substring(bin, 3,4))/10, y=m)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  labs(
    title = "Linear Regression Trends of stdev by plot",
    x = "Bin",
    y = "St Dev averaged by plot"
  ) +
  theme_minimal()


###Chronologies:
rwchron <- chron(rwi)
plot(rwchron$std)


## Some things we can do with chronologies:
# Visualise with smoothing:
yrs <- time(rwchron)
dat <- rwchron$std

spl128 <- caps(dat,nyrs=128)
spl64 <- caps(dat,nyrs=64)
spl32 <- caps(dat,nyrs=32)

par(mar=rep(2.5,4),mgp=c(1.2,0.25,0),tcl=0.5, xaxs="i",yaxs="i")
my.cols <- c("#1B9E77", "#D95F02", "#7570B3")
plot(yrs,dat,type="n",xlab="Year",ylab="RWI",axes=FALSE)
grid(col="black",lwd=0.5)
abline(h=1)
lines(yrs,dat,col="grey",lwd=1)
lines(yrs,spl128,col=my.cols[1],lwd=2)
lines(yrs,spl64,col=my.cols[2],lwd=2)
lines(yrs,spl32,col=my.cols[3],lwd=2)
axis(1);axis(2);axis(3);axis(4)
legend("topright", c("dat", "128yrs", "64yrs", "32yrs"), 
       lwd = 2, col = c("grey",my.cols),bg = "white")
box()

n <- length(yrs)
f128 <- 128/n
f128.lo <- lowess(x = yrs, y = dat, f = f128)
f64 <- 64/n
f64.lo <- lowess(x = yrs, y = dat, f = f64)
f32 <- 32/n
f32.lo <- lowess(x = yrs, y = dat, f = f32)

par(mar=rep(2.5,4),mgp=c(1.2,0.25,0),tcl=0.5,xaxs="i",yaxs="i")
plot(yrs,dat,type="n",xlab="Year",ylab="RWI",axes=FALSE)
grid(col="black",lwd=0.5)
abline(h=1)
lines(yrs,dat,col="grey",lwd=1)
lines(yrs,f128.lo$y,col=my.cols[1],lwd=2)
lines(yrs,f64.lo$y,col=my.cols[2],lwd=2)
lines(yrs,f32.lo$y,col=my.cols[3],lwd=2)
axis(1);axis(2);axis(3);axis(4)
legend("topright", c("dat", "f128", "f64", "f32"), 
       lwd = 2, col = c("grey",my.cols),bg = "white")
box()

## Characterise the temporal strucute
# Test significance of autocorreclation
par(mfcol=c(1, 2))
acf(dat) #up to 4 years
pacf(dat) #Not really any

#OR
dat.ar <- ar(dat,order.max = 10)
dat.ar #Also seems to suggest AR1

plot(0:10,dat.ar$aic,type="b",xlab="AR Order",ylab="AIC",
     main="Difference in AIC between each model and the best-fitting model") # 1 or 5

#OR
ar1 <- arima(dat,order=c(1,0,0))
ar2 <- arima(dat,order=c(2,0,0))
ar3 <- arima(dat,order=c(3,0,0))
ar4 <- arima(dat,order=c(4,0,0))
# example output
ar1

BIC(ar1,ar2,ar3,ar4)

#OR
library(forecast)
dat.arima <- auto.arima(dat, ic="bic")
summary(dat.arima)
acf(residuals(dat.arima))

dev.off()
### Frequency domain
redf.dat <- redfit(dat, nsim = 1000)

#With Fourier transform:

par(tcl = 0.5, mar = rep(2.2, 4), mgp = c(1.1, 0.1, 0),xaxs="i")

plot(redf.dat[["freq"]], redf.dat[["gxxc"]],
     ylim = range(redf.dat[["ci99"]], redf.dat[["gxxc"]]),
     type = "n", ylab = "Spectrum", xlab = "Frequency (cycles per year)",
     axes = FALSE)
grid()
lines(redf.dat[["freq"]], redf.dat[["gxxc"]], col = "black",lwd=1.5)
lines(redf.dat[["freq"]], smooth.spline(redf.dat[["ci99"]],spar = 0.8)$y, col = "#D95F02")
lines(redf.dat[["freq"]], smooth.spline(redf.dat[["ci95"]],spar = 0.8)$y, col = "#7570B3")
lines(redf.dat[["freq"]], smooth.spline(redf.dat[["ci90"]],spar = 0.8)$y, col = "#E7298A")
freqs <- pretty(redf.dat[["freq"]])
pers <- round(1 / freqs, 2)
axis(1, at = freqs, labels = TRUE)
axis(3, at = freqs, labels = pers)
mtext(text = "Period (year)", side = 3, line = 1.1)
axis(2); axis(4)
legend("topright", c("dat", "CI99", "CI95", "CI90"), lwd = 2,
       col = c("black", "#D95F02", "#7570B3", "#E7298A"),
       bg = "white")
box()

#Other method
out.wave <- morlet(y1 = dat, x1 = yrs, p2 = 8, dj = 0.1,
                   siglvl = 0.99)
wavelet.plot(out.wave, useRaster=NA, reverse.y = TRUE)

# Extract the frequencies:
install.packages("waveslim")
library(waveslim)

nPwrs2 <- trunc(log(n)/log(2)) - 1
dat.mra <- mra(dat, wf = "la8", J = nPwrs2, method = "modwt",
               boundary = "periodic")
yrsLabels <- paste(2^(1:nPwrs2),"yrs",sep="")

par(mar=c(3,2,2,2),mgp=c(1.25,0.25,0),tcl=0.5,xaxs="i",yaxs="i")
plot(yrs,rep(1,n),type="n", axes=FALSE, ylab="",xlab="",
     ylim=c(-3,38))
title(main="Multiresolution decomposition",line=0.75)
axis(side=1)
mtext("Years",side=1,line = 1.25)
Offset <- 0
dat.mra2 <- scale(as.data.frame(dat.mra))
for(i in nPwrs2:1){
  x <- dat.mra2[,i] + Offset
  lines(yrs,x)
  abline(h=Offset,lty="dashed")
  mtext(names(dat.mra)[[i]],side=2,at=Offset,line = 0)
  mtext(yrsLabels[i],side=4,at=Offset,line = 0)
  Offset <- Offset+5
}
box()

## high pass low pass:
library(signal)
# here are two frequencies that we will use for demonstrating filtering
f1 <- 0.2 # period of 5
f2 <- 0.05 # period of 20

# Initialize the butterworth filter.
# Multiply freq by 2
f1Low <- butter(n=4, W=f1*2, type="low")
f2Low <- butter(n=4, W=f2*2, type="low")

# We will use filtfilt do run the filter. 
# But before that there is a wrinkle. We have to 
# mirror and extend the data near the edge to 
# avoid end effects. The matlab/octave versions do something 
# like this behind the scenes I think. Regardless, it
# is a simple form of padding.

n <- length(dat)
# create a pas that is adds a 1/4 of the length to the front and
# back of the data. This can be tuned according to the freqs in question
# and there is nothing special about 1/4
pad <- floor(n/4)
# pad the data
datPad <- c(dat[pad:1],dat,dat[n:(n-pad)]) 

# run the filter  
datF1Low <- filtfilt(f1Low, datPad)
datF2Low <- filtfilt(f2Low, datPad)
# unpad the filtered data
datF1Low <- datF1Low[(pad+1):(n+pad)]
datF2Low <- datF2Low[(pad+1):(n+pad)]

par(mar=rep(2.5,4),mgp=c(1.2,0.25,0),tcl=0.5,xaxs="i",yaxs="i")
plot(yrs,dat,type="n",xlab="Year",ylab="",axes=FALSE)
grid(col="black",lwd=0.5)
abline(h=1)
lines(yrs,dat,col="grey80",lwd=1)
lines(yrs,datF1Low,col=my.cols[1],lwd=1.5)
lines(yrs,datF2Low,col=my.cols[2],lwd=2)
axis(1);axis(2);axis(3);axis(4)
legend("topright", c("f < 0.2", "f < 0.05"), 
       lwd = 2, col = c(my.cols[1:2]),bg = "white")
box()

redf.dat <- redfit(datF1Low, mctest=FALSE)
par(tcl = 0.5, mar = rep(2.2, 4), mgp = c(1.1, 0.1, 0),xaxs="i")
plot(redf.dat[["freq"]], redf.dat[["gxxc"]],
     xlim=c(0,0.25),
     type = "n", ylab = "Spectrum", xlab = "Frequency (cycles per year)",
     axes = FALSE)
grid()
lines(redf.dat[["freq"]], redf.dat[["gxxc"]], col = "black",lwd=1.5)
freqs <- pretty(redf.dat[["freq"]])
pers <- round(1 / freqs, 2)
axis(1, at = freqs, labels = TRUE)
axis(3, at = freqs, labels = pers)
mtext(text = "Period (year)", side = 3, line = 1.1)
axis(2); axis(4)
box()

f1f2Pass <- butter(n=4, W=c(f2,f1)*2, type="pass")
datPass <- filtfilt(f1f2Pass, datPad)
datPass <- datPass[(pad+1):(n+pad)]

par(mar=rep(2.5,4),mgp=c(1.2,0.25,0),tcl=0.5,xaxs="i",yaxs="i")
plot(yrs,datPass,type="n",xlab="Year",ylab="",axes=FALSE)
grid(col="black",lwd=0.5)
abline(h=0)
lines(yrs,datPass,col=my.cols[3],lwd=1.5)
axis(1);axis(2);axis(3);axis(4)
legend("topright", "0.05 < f < 0.2", 
       lwd = 2, col = c(my.cols[3]),bg = "white")
box()


redf.dat <- redfit(datPass, mctest=FALSE)
par(tcl = 0.5, mar = rep(2.2, 4), mgp = c(1.1, 0.1, 0),xaxs="i")
plot(redf.dat[["freq"]], redf.dat[["gxxc"]],
     xlim=c(0,0.25),
     type = "n", ylab = "Spectrum", xlab = "Frequency (cycles per year)",
     axes = FALSE)
grid()
lines(redf.dat[["freq"]], redf.dat[["gxxc"]], col = "black",lwd=1.5)
freqs <- pretty(redf.dat[["freq"]])
pers <- round(1 / freqs, 2)
axis(1, at = freqs, labels = TRUE)
axis(3, at = freqs, labels = pers)
mtext(text = "Period (year)", side = 3, line = 1.1)
axis(2); axis(4)
box()
}
# So corrected file are the ones where the dpi was corrected, use this!
test1 <- read.rwl("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Raw/CDendro-output/P460_1.rwl")
test1c <- read.rwl("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Raw/CDendro-output/Corrected/P460_1c.rwl")

#### 8 - Climatic data -----
library(sf)
# options(repos=c(CRAN="https://cran.rstudio.com/"))
# renv::snapshot()
# renv::clean()  # for good measure
# renv::rebuild()

plot.data <- read.csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Raw/Plot-info.csv")
coords <- plot.data %>% 
  select(Plot, Longitude, Latitude) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

directories <- basename(Sys.glob("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Raw/ClimateData/*"))

extract_all_climatic_variables <- function(directory, coordinates) {
  climate_data <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Raw/ClimateData/", 
                                  directory, "/", directory,".csv")) %>% 
    mutate(time = as.numeric(substring(time, 1,4))) %>% 
    rename(year=time)
  
  climate_sf <- climate_data %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326)
  
  
  years <- unique(climate_sf$year)
  results <- list()
  # Loop over each year and extract data
  for (yr in years) {
    # Filter raster data for the current year
    climate_year <- climate_sf %>% dplyr::filter(year == yr)
    
    # Perform a spatial join for the current year
    matched <- st_join(coords, climate_year, join = st_nearest_feature)
    
    # Add the year column explicitly (it will be consistent)
    matched <- matched %>%
      mutate(year = yr) %>%
      select(Plot, year, 
             paste0("ssp245_", directory, "_p50"),
             paste0("ssp126_", directory, "_p50"),
             paste0("ssp585_", directory, "_p50"),
      ) %>% 
      st_drop_geometry()# Keep relevant columns
    
    # Append the results to the list
    results[[as.character(yr)]] <- matched
  }
  
  results.comb <- bind_rows(results, .id = "column_label")
  write_csv(results.comb, 
            paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/ClimateData/", directory, ".csv"))
}

lapply(directories, extract_all_climatic_variables, coords)

#### 9 - Create csv of fitted parameters -----
 
extract.RDS.params <- function(tree, directory){
  # Check if -cc exists since this is the updated version
  if(file.exists(paste0(directory, tree, "-cc.RDS"))){
    readRDS(paste0(directory, tree, "-cc.RDS"))$`Fitted Parameters`
  } else { readRDS(paste0(directory, tree, ".RDS"))$`Fitted Parameters`}
}

param.list <- lapply(info$Plot, FUN = function(x){
  df <- extract.RDS.params(paste0(x, "-1"), "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Fits/")
  for(i in 2:10){
    df <- rbind(df, extract.RDS.params(paste0(x, "-", as.character(i)), "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Fits/"))
  }
  rownames(df) <- paste0(x, "-", 1:10)
  df
  }
)

params_all <- bind_rows(param.list, .id = "column_label") %>% 
  rownames_to_column(var = "tree") %>% 
  select(-column_label)

# Add extra tree: 
param_extra <- as.data.frame(c("P475.1-e", readRDS("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Fits/P475.1-e.RDS")$`Fitted Parameters`))
colnames(param_extra) <- colnames(params_all)

params_all <- rbind(params_all, param_extra)

write_csv(params_all, "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Parameters_fitted.csv")

#### 10: Simulate A values for other T and atmospheric Ci ----
library(dplyr)
library(magrittr)
library(photosynthesis)
library(future)
plan("multisession")

## Check out some of the function:
a <- readRDS(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Fits/", "P460.2-2", ".RDS"))
a$Plot +
  geom_point(data=a$Data, aes(x=Ci, y=A_carbox), color="red") +
  geom_point(data=a$Data, aes(x=Ci, y=A_regen), color="blue") +
  geom_point(data=a$Data, aes(x=Ci, y=A_model), color="darkgreen")

plot(a$Data$Ci, a$Data$A_model)
aD <- a$Data

### Compare the data/fit with the projected values by plantecophys and photosynthesis packages
# photosynthesis:
#-#

bake_par   = make_bakepar(
  replace = list(
    # Ea_gammastar = set_units(a$`Fitted Parameters`$Ea_gamma_star, "J / mol") # no effect
  ))                      # temperature response parameters
constants  = make_constants(use_tealeaves = F) # physical constants
leaf_par   = make_leafpar(
  replace = list(
    # g_uc = set_units(0.1, "mol / m^2 / s"), # Seems to shift it in the right way?
    #g_sc = set_units(0.2, "mol / m^2 / s"), # Seems ok
    T_leaf = set_units(mean(a$Data$Tleaf) + 273.15, "K"),
    V_cmax25 = set_units(a$`Fitted Parameters`$V_cmax, "umol / m^2 / s"),
    J_max25 = set_units(a$`Fitted Parameters`$J_max, "umol / m^2 / s"),
    R_d25 = set_units(abs(a$`Fitted Parameters`$R_d), "umol / m^2 / s"),
    gamma_star25 = set_units(a$`Fitted Parameters`$gamma_star25, "umol / mol"),
    g_mc25 = set_units(a$`Fitted Parameters`$g_mc25 + 0.1, "mol / m^2 / s"),
    theta_J = set_units(a$`Fitted Parameters`$theta_J, "")
    ), use_tealeaves = F)   # leaf parameters
enviro_par = make_enviropar(
  replace = list(
    C_air = set_units(seq(90, 900, by=90), "umol/mol"),
    PPFD = set_units(a$Data$Qabs[1], "umol/m^2/s")
  ), use_tealeaves = F) # environmental parameters
#-#
sim.photos18 <- photosynthesis(leaf_par, enviro_par, bake_par, constants, quiet = TRUE,
                    use_tealeaves = FALSE)
#plot to compare
plot(a$Data$Cicor, a$Data$Acor, pch=20)
points(a$Data$Cicor, a$Data$A_model, col="blue", pch=20)
points(sim.photos$C_i, sim.photos$A, col="red", pch=19) #photosynthesis

points(sim.plante$Ci, sim.plante$ALEAF, col="orange", pch=19) #planteco

#to compare mesophyl conductance
plot(a$Data$Cicor, a$Data$Acor, pch=20)
points(a$Data$Cicor, a$Data$A_model, col="blue", pch=20)
points(sim.photos$C_i, sim.photos$A, col="red", pch=19) # fitting default
points(sim.photos18$C_i, sim.photos18$A, col="orange", pch=19) # mesophyl at 0.18
points(sim.photos38$C_i, sim.photos38$A, col="cyan", pch=19) # mesophyl at 0.38

# plantecophys
library(plantecophys)
sim.plante <- Photosyn(VPD = 1.5, 
                       Ca = seq(90,900, by=90),
                       PPFD = a$Data$Qabs[1],
                       Tleaf = 25,
                       Patm = 101,325,
                       Jmax = a$`Fitted Parameters`$J_max,
                       Vcmax = a$`Fitted Parameters`$V_cmax)

points(sim.plante$Ci, sim.plante$ALEAF, col="orange", pch=19)

######## OLD ##
#-#
bake_par   = make_bakepar(replace = list())                      # temperature response parameters
constants  = make_constants(replace = list(
), use_tealeaves = FALSE) # physical constants
leaf_par   = make_leafpar(
  replace = list(
  T_leaf = set_units(c(288.15, 293.15, 298.15, 303.15), "K"),
  V_cmax25 = set_units(43.8, "umol / m^2 / s"),
  J_max25 = set_units(81.7, "umol / m^2 / s")
), use_tealeaves = FALSE)   # leaf parameters
enviro_par = make_enviropar(
  replace = list(
    C_air = set_units(c(300, 420, 450, 500, 600, 700, 800), "umol/mol"),
    PPFD = set_units(1681, "umol/m^2/s")
), use_tealeaves = FALSE) # environmental parameters
#-#

b <- photosynthesis(leaf_par, enviro_par, bake_par, constants, quiet = TRUE,
      use_tealeaves = FALSE)


####### Automate: 
info <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Measurements-information.csv")

extract.RDS <- function(tree, directory){
  # Check if -cc exists since this is the updated version
  if(file.exists(paste0(directory, tree, "-cc.RDS"))){
    readRDS(paste0(directory, tree, "-cc.RDS"))
  } else { readRDS(paste0(directory, tree, ".RDS"))}
}

simulate.photo <- function(tree){
  aci.fit <- extract.RDS(tree, 
                                "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Fits/")
  
  #-#
  bake_par   = make_bakepar()                      # temperature response parameters
  constants  = make_constants(use_tealeaves = FALSE) # physical constants
  leaf_par   = make_leafpar(
    replace = list(
      T_leaf = set_units(c(288.15, 293.15, 298.15, 303.15), "K"),
      V_cmax25 = set_units(aci.fit$`Fitted Parameters`$V_cmax, "umol / m^2 / s"),
      J_max25 = set_units(aci.fit$`Fitted Parameters`$J_max, "umol / m^2 / s"),
      R_d25 = set_units(abs(aci.fit$`Fitted Parameters`$R_d), "umol / m^2 / s"),
      gamma_star25 = set_units(aci.fit$`Fitted Parameters`$gamma_star25, "umol / mol"),
      #g_mc25 = set_units(aci.fit$`Fitted Parameters`$g_mc25 + 0.1, "mol / m^2 / s"),
      theta_J = set_units(aci.fit$`Fitted Parameters`$theta_J, "")
    ), use_tealeaves = FALSE)   # leaf parameters
  enviro_par = make_enviropar(
    replace = list(
      C_air = set_units(c(280, 400, 420, 450, 600, 800, 840, 870, 1140), "umol/mol"),
      PPFD = set_units(aci.fit$Data$Qabs[1], "umol/m^2/s")
    ), use_tealeaves = FALSE) # environmental parameters
  #-#
  
  sims <- photosynthesis(leaf_par, enviro_par, bake_par, constants, quiet = TRUE,
                         use_tealeaves = FALSE, parallel = TRUE)
  write_csv(sims, paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Simulations/",
                         tree,
                         ".csv"))
}

lapply(info$Plot, FUN = function(p) {
  lapply(as.character(1:10), FUN = function(q){
    simulate.photo(paste(p, q, sep = "-"))
  })
})

sim.all.list <- lapply(info$Plot, FUN = function(p) {
  sim.plot.list <- lapply(as.character(1:10), FUN = function(q) {
    sim <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Simulations/", paste(p, q, sep = "-"), ".csv"), show_col_types = F)
    sim.wide <- sim %>% 
      mutate(simulation = paste(as.character(C_air), as.character(T_leaf - 273.15), sep = "-")) %>% 
      select(A, C_i, simulation) %>% 
      pivot_wider(names_from = simulation, values_from = c(C_i, A), names_prefix = "sim")
  } )
  print(sim.plot.list)
  sim.plot <- bind_rows(sim.plot.list) %>% 
    mutate(tree = paste(p, as.character(1:10), sep = "-"))
  print(sim.plot)
} )

write_csv(bind_rows(sim.all.list), "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Simulations.csv")


#TEST:

tr <- "P465.2-5"
sims <- simulations %>% filter(tree == tr)
meas <- read.csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split-corrected/", tr, ".csv"))
test.plante <- simulate.plante(tr)
plot(meas$Cicor, meas$Acor, pch=20)
points(sims[, 1:28], sims[, 29:56], col = "blue", pch=20)
points(sims$`C_i_sim420-25`, sims$`A_sim420-25`, col = "red", pch=20)
points(test.plante$Ci, test.plante$ALEAF, col = "orange", pch=20)
points(sims$`C_i_sim600-30`, sims$`A_sim600-30`, col = "red", pch=20)

aci.fit$`Fitted Parameters`$V_cmax
aci.fit$Data$Qabs[1]

test.plante <- simulate.plante(tr)


simulate.plante <- function(tree){
  aci.fit <- extract.RDS(tree, 
                         "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Fits/")
  
  fitaci1 <- Photosyn(
    VPD = 1.5,
    Ca = 420,
    PPFD = aci.fit$Data$Qabs[1],
    Tleaf = 25,
    Patm = 101,
    RH = NULL,
    gsmodel = c("BBOpti", "BBLeuning", "BallBerry", "BBdefine"), #####
    g1 = 4,
    g0 = 0,
    gk = 0.5,
    vpdmin = 0.5,
    D0 = 5,
    GS = NULL,
    BBmult = NULL,
    alpha = 0.24,
    theta = 0.85,
    Jmax = aci.fit$`Fitted Parameters`$J_max,
    Vcmax = aci.fit$`Fitted Parameters`$V_cmax,
    gmeso = NULL,
    TPU = 1000,
    alphag = 0,
    Rd0 = 2,
    Q10 = 1.92,
    Rd = NULL,
    TrefR = 25,
    Rdayfrac = 1,
    EaV = 58550,
    EdVC = 2e+05,
    delsC = 629.26,
    EaJ = 29680,
    EdVJ = 2e+05,
    delsJ = 631.88,
    GammaStar = NULL,
    Km = NULL,
    Ci = NULL,
    Tcorrect = TRUE,
    returnParsOnly = FALSE,
    whichA = c("Ah", "Amin", "Ac", "Aj") )
}

#### 11 Create one master file ----
## Load in data
gasex.param <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP_GasExchangeParameters.csv")
leaf.dim <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP_LeafDimensions.csv")
tree.diam <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP_TreeDbh.csv")
plot.info <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP-PlotInfo.csv") 
ring.analysis <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP_RingAnalysis.csv") 

#To finalise:
simulations <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Simulations.csv") 

## (For now?) these do not need recombinations:
tree.community <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP_TreeCommunity.csv") 
understory.community <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP-UnderstoryCommunity.csv") 
soil <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/YBP-SoilLayers.csv") 

# Rework slightly to be able to join later
leaf.dim <- rename(leaf.dim, tree = LeafName)

data.master <- gasex.param  %>% 
  select(tree, V_cmax, V_cmax_se, J_max, J, J_se, R_d, R_d_se) %>% 
  mutate(Plot = substring(tree, 1, 6),
         Bin = substring(tree, 1, 4), 
         Stage = ifelse((as.numeric(substring(tree, nchar(tree), nchar(tree))) %in% 1:5) | (substring(tree, nchar(tree), nchar(tree)) == "e"),
                        "adult", 
                        "sapling")) %>% 
  relocate(c(Plot, Bin, Stage), .after = tree) %>% 
  left_join(select(tree.diam, tree, diameter), by = "tree") %>% 
  left_join(leaf.dim, by = "tree") %>% 
  left_join(select(plot.info, Plot, Longitude, Latitude), by = "Plot") %>% 
  left_join(simulations %>% relocate(tree, .before = everything()), by = "tree") %>% 
  left_join(ring.analysis, by = "tree")

write_csv(data.master, "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Final_Database/Combined.csv")

