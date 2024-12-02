### Script to funnel all the different raw data files into one coherent database
### Author: Lukas Van Riel
### Date: 25/09/2024

### Load packages
library(tidyverse)
library(conflicted)

#### 1 - Leaf parameters #### 
#Load txt data file:
data.leaf.raw <- read_table("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Raw/Lukas-LeafData/Leaf-Scans-copy.txt")

data.leaf <- data.leaf.raw %>% 
  mutate(LeafName = as.character(WinFOLIA),
         LeafArea = as.numeric(NofLeaves), 
         Width = as.numeric(TotLeafArea), 
         Length = as.numeric(AvgLeafArea)) %>% 
  select(LeafName, LeafArea, Width, Length)

write_csv(data.leaf, "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Leaf-Dimensions.csv")

#### 2 -  Split Racir Data into individual measurements #### 
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
data_raw <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/parameters.csv")
data_plot <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Raw/Plot-info.csv")


data <- data_raw %>% 
  mutate(full=ID) %>% 
  separate_wider_delim(ID, "-", names=c("Plot", "Tree")) %>% 
  left_join(data_plot, by="Plot") %>% 
  separate_wider_delim(Plot, ".", names=c("Bin", "Plot")) %>% 
  mutate(Bin=factor(paste0(substring(Bin, 2,3), ".", substring(Bin, 4,4)))) %>% 
  mutate(sapling=as.numeric(Tree)%/%6) %>% 
  mutate(Class=factor(ifelse(sapling==1, "sapling", "adult"))) %>% 
  relocate(Class, .before = Vcmax)

data$Longitude

ggplot(data) +
  geom_point(aes(x=Vcmax, y=Jmax, colour = Class))

ggplot(data) +
  geom_point(aes(x=Latitude, y=Vcmax, colour = Class))

ggplot(data) +
  geom_point(aes(x=Latitude, y=Jmax, colour = Class))

ggplot(data) +
  geom_boxplot(aes(y=Jmax, x=Bin, colour = Class))

ggplot(data) +
  geom_boxplot(aes(y=Vcmax, x=Bin, colour = Class))

#!!!
ggplot(data) +
  geom_boxplot(aes(y=Vcmax, x=Bin, colour = Method))

# Run some anova's
install.packages("car")
library(car)

aov_Va <- aov(Vcmax ~ Bin,
                         data = data %>% filter(Class=="adult"))
summary(aov_Va)
ggplot(data %>% filter(Class=="adult")) +
  geom_boxplot(aes(y=Vcmax, x=Bin))

aov_Vs <- aov(Vcmax ~ Bin,
              data = data %>% filter(Class=="sapling"))
summary(aov_Vs)
ggplot(data %>% filter(Class=="sapling")) +
  geom_boxplot(aes(y=Vcmax, x=Bin))

aov_Ja <- aov(Jmax ~ Bin,
              data = data %>% filter(Class=="adult"))
summary(aov_Ja)
ggplot(data %>% filter(Class=="adult")) +
  geom_boxplot(aes(y=Jmax, x=Bin))

aov_Js <- aov(Jmax ~ Bin,
              data = data %>% filter(Class=="sapling"))
summary(aov_Js)
ggplot(data %>% filter(Class=="sapling")) +
  geom_boxplot(aes(y=Jmax, x=Bin))




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

# Split some columns and drop some unimportant or constant
data_split <- data_gran %>% 
  separate_wider_delim(Name, delim = "-", names = c("Plot", "Horizon")) %>% 
  separate_wider_delim(Date, delim = " ", names = c("Date", "Time"))
  
# Create new columns for sand, silt and clay
data_class <- data_split %>% 
  mutate(clay=bin0to2, 
         silt= bin2to4 + bin4to8 + bin8to15 + bin15to31 + bin31to53,
         sand= bin53to106 + bin106to250 + bin250to500 + bin500to1000 + bin1000to2000) %>% 
  mutate(clay_n = clay/(clay+silt+sand),
         silt_n = silt/(clay+silt+sand),
         sand_n = sand/(clay+silt+sand)) %>% 
  mutate(soil_class=ifelse(sand_n<0.50, "SiLo",
                      ifelse(sand_n<0.70, "SaLo",
                      ifelse(sand_n<0.90, "LoSa", "Sa"))))

write_csv(data_class[-c(1:3),], "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Soil_Granulometrie.csv")

table(data_class$soil_class)


#### 6 - Soil horizons description #### 
data_soil <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Raw/Soil-layers.csv",
                      col_types = c())

data_soil_extended <- data_soil %>% 
  separate_wider_delim(Munsell_color, delim = " ", names = c("Munsell_hue", "Munsell_tosplit")) %>% 
  separate_wider_delim(Munsell_tosplit, delim = "/", names = c("Munsell_value", "Munsell_chroma"))

write_csv(data_soil_extended, "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Soil_Layers.csv")

#### 7 - CDendro output #### 
library(dplR)

files <- Sys.glob("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Raw/CDendro-output/Corrected/*.rwl")

basic_ring_analysis <- function(file) {
  #Create list to be returned by function
  results <- list()
  
  ## Read in the file
  if(substr(file, nchar(file)-2, nchar(file)) == "rwl") {
    ring.data.raw <- read.rwl(file, format="tucson")
  } else {if(substr(file, nchar(file)-2, nchar(file)) == "csv") {
        ring.data.raw <- read_csv(file)
      }
    }
  
  #print(ncol(ring.data.raw))
  
  ## Calculate age of each tree
  results$age <- ring.data.raw %>%
    summarise(across(everything(), ~ sum(!is.na(.))))
  
  ## Calculate the average ring widths per decade
  # Add years and decades as columns
  ring_data <- ring.data.raw %>% 
    mutate(Year = as.numeric(rownames(ring.data.raw)),
           Decade = floor(Year / 10) * 10)
  
  # Pivot longer and summarise by decade and plot
  results$averages <- ring_data %>%
    pivot_longer(
      cols = -c(Year, Decade), # Keep Year and Decade as identifiers
      names_to = "Plot",
      values_to = "Width") %>%
    group_by(Decade, Plot) %>%
    summarise(
      Decade_Avg = ifelse(all(is.na(Width)), NA, mean(Width, na.rm = TRUE)),
      Valid_Years = sum(!is.na(Width)),
      NA_Years = sum(is.na(Width)),
      .groups = "drop" ) # Drop grouping for cleaner output
  
  
  ## Return ages and decadal averages as a list
  results
}

ring_basics <- setNames(lapply(files, FUN = basic_ring_analysis),
                        substring(files, nchar(files[1]) - 10, nchar(files[1]) - 5))

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



### Let's go trhough the dplr package documentation----
## Try for 1 plot:
rwl <- read.rwl(files[1], format="tucson")

#corr.rwl.seg(rwl, seg.length = 50, bin.floor = 20)

## Get some basic information
plot(rwl, plot.type="spag")

rwl.report(rwl)
rwl.stats(rwl)
# Mean sensitivity used to be included here as well, but book said to not use it!
# At best it is as good as the stdev of time series in cases of high autocorrelation

#Can also plot some stuff if would be interesting, eg for autocorrelation:
ar1 <- data.frame(x="P460.1",y=rwl.stats(rwl)$ar1)
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
detrend.series(y = rwl[, "P460.1-1"], verbose=TRUE)

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

write.detrend <- function(file) {
  plotname <- substring(file, nchar(file)-10, nchar(file)-5)
 
  #
  rwl <- read.rwl(file, format="tucson")
  
  p.rwl <- rwl %>% 
    mutate(year=as.numeric(rownames(rwl))) %>% 
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
  
  pdf(paste0('/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/TreeCores-Detrended/Figures/', plotname, '-Detrended-ModNegExp.pdf'))
  print(p.rwi)
  dev.off()
  
  write_csv(rwi, paste0('/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/TreeCores-Detrended/', plotname, '-Detrended-ModNegExp.csv'))
  
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

lapply(files, write.detrend)



## Let's see now how some of the trends are now:
files.csv <- Sys.glob("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/TreeCores-Detrended/*.csv")

ring_detrended <- setNames(lapply(files.csv, FUN = basic_ring_analysis),
                        substring(files.csv, nchar(files[1]) - 9, nchar(files[1]) - 4))

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
