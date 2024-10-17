### Script to funnel all the different raw data files into one coherent database
### Author: Lukas Van Riel
### Date: 25/09/2024

### Load packages
library(tidyverse)

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

