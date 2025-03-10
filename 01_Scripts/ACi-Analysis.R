### Script to funnel all the different raw data files into one coherent database
### Author: Lukas Van Riel
### Date: 25/09/2024


### Load packages -----
library(tidyverse)
library(plantecophys)
library(racir)
library(photosynthesis)
#library(ptinpoly)

### Load data -----
info <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Measurements-information.csv")

### Functions -----
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
  
  #Now analyse with bilinear function
  fit <- fitaci(corrected, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                           PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='bilinear')
  
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

### Body -----
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

write_csv(parameters_df, "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/parameters_bilinear.csv")

# compare the values

params <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/parameters.csv")
params_bilinear <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/parameters_bilinear.csv")

params_bil <- params_bilinear[params$Method == "default",]
params <- params[params$Method == "default",]

params_bil$ID == params$ID

params_diff <- params_bil$Vcmax-params$Vcmax
names(params_diff) <- params$ID

params_diff

# Seems the bilinear one is a better fit, and is able to always fit the data. 
# We'll continue with this for now

data_plot <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Raw/Plot-info.csv")

data <- params_bilinear %>% 
  mutate(full=ID) %>% 
  separate_wider_delim(ID, "-", names=c("Plot", "Tree")) %>% 
  left_join(data_plot, by="Plot") %>% 
  separate_wider_delim(Plot, ".", names=c("Bin", "Plot")) %>% 
  mutate(Bin=factor(paste0(substring(Bin, 2,3), ".", substring(Bin, 4,4)))) %>% 
  mutate(sapling=as.numeric(Tree)%/%6) %>% 
  mutate(Class=factor(ifelse(sapling==1, "sapling", "adult"))) %>% 
  relocate(Class, .before = Vcmax)

ggplot(data) +
  geom_point(aes(x=Vcmax, y=Jmax, colour = Class))

ggplot(data) +
  geom_point(aes(x=Latitude, y=Vcmax, colour = Class))

ggpalot(data) +
  geom_point(aes(x=Latitude, y=Jmax, colour = Class))

ggplot(data) +
  geom_boxplot(aes(y=Jmax, x=Bin, colour = Class))

ggplot(data) +
  geom_boxplot(aes(y=Vcmax, x=Bin, colour = Class))


# Run some anova's
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


### Test some things for fitaci:----
# What is the difference between Ca and Ci?

data_test <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                             "P460.1", "-", as.character(info[info$Plot=="P460.1",3]), ".csv"))
empty_test <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                              "P460.1", "-", as.character(info[info$Plot=="P460.1","empty"]), ".csv"))

plot(data_test$Ci, data_test$Ca, pch=19)
plot(data_test$Ci)
points(data_test$Ca, col="red")

### Get some actual measurements: ---- 
# meas1: perfect, meas2: good, meas3: good with cut at end
# meas4: late transition point, meas5: tricky, meas6: bad/unusable(?)

# Get the data
m1 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                             "P460.2", "-", as.character(info[info$Plot=="P460.2","adult 2"]), ".csv"))
e1 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P460.2", "-", as.character(info[info$Plot=="P460.2","empty"]), ".csv"))
m2 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P465.1", "-", as.character(info[info$Plot=="P465.1","adult 2"]), ".csv"))
e2 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P465.1", "-", as.character(info[info$Plot=="P465.1","empty"]), ".csv"))
m3 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P460.1", "-", as.character(info[info$Plot=="P460.1","adult 5"]), ".csv"))
e3 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P460.1", "-", as.character(info[info$Plot=="P460.1","empty"]), ".csv"))
m4 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P470.3", "-", as.character(info[info$Plot=="P470.3","adult 4"]), ".csv"))
e4 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P470.3", "-", as.character(info[info$Plot=="P470.3","empty"]), ".csv"))
m5 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P465.1", "-", as.character(info[info$Plot=="P465.1","sapling 5"]), ".csv"))
e5 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P465.1", "-", as.character(info[info$Plot=="P465.1","empty"]), ".csv"))
m6 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P470.3", "-", as.character(info[info$Plot=="P470.3","adult 1"]), ".csv"))
e6 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P470.3", "-", as.character(info[info$Plot=="P470.3","empty"]), ".csv"))


plot(m1$Ca, m1$CO2_r)
plot(m1$CO2_s, m1$Ca)

# M1: P460.2-2
m1_cor <- racircal(data = as.data.frame(m1), caldata = as.data.frame(e1), mincut = 10, maxcut = 1050)
m1_cor <- m1_cor %>% filter(Cicor<1400, Cicor>0)

plot(m1_cor$Cicor, m1_cor$Acor, pch=20, ylim=c(0,20), xlim=c(0,1100))
points(m1$Ci, m1$A, pch=20, col="red")
points(e1$Ci, e1$A, pch=20, col="blue")

fit1 <- fitaci(m1_cor, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                         PPFD = "Qabs", Rd ="Rd"), Tcorrect=T, fitmethod='bilinear')

fit1$pars
plot(fit1, ylim=c(0,22))

# Test effect of mesophyll conductance: 

fit1 <- fitaci(m1_cor, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                       PPFD = "Qabs", Rd ="Rd"), Tcorrect=T)
plot(fit1, ylim=c(0,22))


# M2: P465.1-2
m2_cor <- racircal(data = as.data.frame(m2), caldata = as.data.frame(e2), mincut = 10, maxcut = 1050)
m2_cor <- m2_cor %>% filter(Cicor<1400, Cicor>0)

plot(m2_cor$Cicor, m2_cor$Acor, pch=20, ylim=c(-2,11), xlim=c(0,1100))
points(m2$Ci, m2$A, pch=20, col="red")
points(e2$Ci, e2$A, pch=20, col="blue")

fit2 <- fitaci(m2_cor, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                        PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='bilinear')

fit2$pars
plot(fit2)


#M3 P460.1-5
m3_cor <- racircal(data = as.data.frame(m3), caldata = as.data.frame(e3), mincut = 10, maxcut = 1050)
m3_cor <- m3_cor %>% filter(Cicor<1400, Cicor>0)

plot(m3_cor$Cicor, m3_cor$Acor, pch=20, ylim=c(-2,11), xlim=c(0,1200))
points(m3$Ci, m3$A, pch=20, col="red")
points(e3$Ci, e3$A, pch=20, col="blue")

fit3.d <- fitaci(m3_cor, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                       PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='default')
fit3.b <- fitaci(m3_cor, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                       PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='bilinear')

fit3.d$RMSE #RMSE = 14.47578 
fit3.b$RMSE

fit3$pars
plot(fit3, ylim=c(0,10))
#Does it now make sense to play a bit with the transition point?
fit3b <- fitaci(m3_cor, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                       PPFD = "PARi", Rd ="Rd"),
                Tcorrect=T,
                citransition = 700,
                fitmethod='default')
plot(fit3b, ylim=c(0,10))
fit3b$RMSE

#M4 P470.3-4
m4_cor <- racircal(data = as.data.frame(m4), caldata = as.data.frame(e4), mincut = 10, maxcut = 1100)
m4_cor <- m4_cor %>% filter(Cicor<1400, Cicor>0)

plot(m4_cor$Cicor, m4_cor$Acor, pch=20, ylim=c(-2,20), xlim=c(0,1100))
points(m4$Ci, m4$A, pch=20, col="red")
points(e4$Ci, e4$A, pch=20, col="blue")

fit4b <- fitaci(m4_cor, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                       PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='bilinear')
fit4d <- fitaci(m4_cor, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                        PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='default')

fit4d$pars
fit4b$pars

plot(fit4d)
plot(fit4b)

fit4d$RMSE 
fit4b$RMSE


#M5 P465.1-10
m5_cor <- racircal(data = as.data.frame(m5), caldata = as.data.frame(e5), mincut = 10, maxcut = 1100)
m5_cor <- m5_cor %>% filter(Cicor<1400, Cicor>0)

plot(m5_cor$Cicor, m5_cor$Acor, pch=20)
plot(m5_cor$Cicor, m5_cor$Acor, pch=20, ylim=c(-2,20), xlim=c(0,1100))
points(m5$Ci, m5$A, pch=20, col="red")
points(e5$Ci, e5$A, pch=20, col="blue")

fit5b <- fitaci(m5_cor, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                        PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='bilinear')
fit5d <- fitaci(m5_cor, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                        PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='default')

fit5b$pars
fit5d$pars

plot(fit5d)
plot(fit5b)
plot(fit5d, ylim=c(0,21))
plot(fit5b, ylim=c(0,21))

fit5d$RMSE 
fit5b$RMSE

#There seems to be one outlier skewing the RMSE's:
m5_cor2 <- m5_cor %>% filter(Acor<30)
plot(m5_cor2$Cicor, m5_cor2$Acor, pch=20)

fit5b2 <- fitaci(m5_cor2, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                        PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='bilinear')
fit5d2 <- fitaci(m5_cor2, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                        PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='default')

fit5b2$pars
fit5d2$pars

plot(fit5d2)
plot(fit5b2)
poly.vertices <- data.frame(x=c(0,0,1000, 1000, 425), y=c(0,23,23, 22, 22))
polygon(poly.vertices$x, poly.vertices$y, border="maroon", col=NA)
outside <- (pip2d(as.matrix(poly.vertices), as.matrix(m5_cor2 %>% select(Cicor, Acor))) < 0)

plot(m5_cor2$Cicor, m5_cor2$Acor, pch=20, col="black")
points(m5_cor2$Cicor[outside==F], m5_cor2$Acor[outside==F], pch=20, col="green")
points(m5_cor2$Cicor[outside==T], m5_cor2$Acor[outside==T], pch=20, col="red")

fit5d2$RMSE 
fit5b2$RMSE

#Try it again but drop some points: 

plot(fit5b2)
poly.vertices <- data.frame(x=c(0,0,1000, 1000, 425), y=c(0,23,23, 22, 22))
polygon(poly.vertices$x, poly.vertices$y, border="maroon", col=NA)
outside <- (pip2d(as.matrix(poly.vertices), as.matrix(m5_cor2 %>% select(Cicor, Acor))) < 0)

plot(m5_cor2$Cicor, m5_cor2$Acor, pch=20, col="black")
points(m5_cor2$Cicor[outside==F], m5_cor2$Acor[outside==F], pch=20, col="green")
points(m5_cor2$Cicor[outside==T], m5_cor2$Acor[outside==T], pch=20, col="red")

# correct based on this: 
m5_cor3 <- m5_cor2 %>% filter(! outside)

fit5b3 <- fitaci(m5_cor3, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                          PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='bilinear')
fit5d3 <- fitaci(m5_cor3, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                          PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='default')

fit5d3$pars
fit5b3$pars

plot(fit5d3)
plot(fit5b3)

fit6d3$RMSE
fit5b3$RMSE

#M6 P470.3-1
m6_cor <- racircal(data = as.data.frame(m6), caldata = as.data.frame(e6), mincut = 10, maxcut = 1200)
m6_cor <- m6_cor %>% filter(Cicor<1400, Cicor>0)

plot(m6_cor$Cicor, m6_cor$Acor, pch=20)
plot(m6_cor$Cicor, m6_cor$Acor, pch=20, ylim=c(-2,20), xlim=c(0,1100))
points(m6$Ci, m6$A, pch=20, col="red")
points(e6$Ci, e6$A, pch=20, col="blue")

fit6b <- fitaci(m6_cor, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                        PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='bilinear')
#fit6d <- fitaci(m6_cor, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
#                                        PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='default')

#fit6d$pars
fit6b$pars

plot(fit6b, ylim=c(0,12))

#plot(fit5b, ylim=c(0,21))

#fit6d$RMSE 
fit6b$RMSE

# Try to correct with the polygon:
poly.vertices <- data.frame(x=c(-1,-1,1100, 1100, 900, 50, 0), y=c(-1,30,30,10.5, 10.5, 1.5, -1))
polygon(poly.vertices$x, poly.vertices$y, border="maroon", col=NA)
outside <- (pip2d(as.matrix(poly.vertices), as.matrix(m6_cor %>% select(Cicor, Acor))) < 0)

plot(m6_cor$Cicor, m6_cor$Acor, pch=20, col="black")
points(m6_cor$Cicor[outside==F], m6_cor$Acor[outside==F], pch=20, col="green")
points(m6_cor$Cicor[outside==T], m6_cor$Acor[outside==T], pch=20, col="red")

# correct based on this: 
m6_cor2 <- m6_cor %>% filter(! outside)

fit6b2 <- fitaci(m6_cor2, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                        PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='bilinear')
fit6d2 <- fitaci(m6_cor2, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor",
                                        PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='default')

fit6d2$pars
fit6b2$pars

plot(fit6d2)
plot(fit6b2)

fit6d2$RMSE
fit6b2$RMSE

#M7 P470.2-4
m7 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P470.2", "-", as.character(info[info$Plot=="P470.2","adult 4"]), ".csv"))
e7 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P470.2", "-", as.character(info[info$Plot=="P470.2","empty"]), ".csv"))

m7_cor <- racircal(data = as.data.frame(m7), caldata = as.data.frame(e7), mincut = 10, maxcut = 1120)
m7_cor <- m7_cor %>% filter(Cicor<1400, Cicor>0)

plot(m7_cor$Cicor, m7_cor$Acor, pch=20)
plot(m7_cor$Cicor, m7_cor$Acor, pch=20, ylim=c(-2,20), xlim=c(0,1100))
points(m7$Ci, m7$A, pch=20, col="red")
points(e7$Ci, e7$A, pch=20, col="blue")

fit7b <- fitaci(m7_cor, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                        PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='bilinear')
fit6d <- fitaci(m6_cor, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor",
                                        PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='default')

plot(fit7b, ylim=c(0,15))

# correct
poly.vertices <- data.frame(x=c(-1,-1,1100, 1100, 550, 50, 0), y=c(-1,25,25,15, 15, 1.5, -1))
polygon(poly.vertices$x, poly.vertices$y, border="maroon", col=NA)
outside <- (pip2d(as.matrix(poly.vertices), as.matrix(m7_cor %>% select(Cicor, Acor))) < 0)

plot(m7_cor$Cicor, m7_cor$Acor, pch=20, col="black")
points(m7_cor$Cicor[outside==F], m7_cor$Acor[outside==F], pch=20, col="green")
points(m7_cor$Cicor[outside==T], m7_cor$Acor[outside==T], pch=20, col="red")

m7_cor2 <- m7_cor %>% filter(! outside)

fit7b2 <- fitaci(m7_cor2, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                          PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='bilinear')
fit7d2 <- fitaci(m7_cor2, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor",
                                          PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='default')

fit7d2$pars
fit7b2$pars

plot(fit7d2)
plot(fit7b2)

fit7d2$RMSE
fit7b2$RMSE

#M8 P470.2-1
m8 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P470.2", "-", as.character(info[info$Plot=="P470.2","adult 1"]), ".csv"))
e8 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P470.2", "-", as.character(info[info$Plot=="P470.2","empty"]), ".csv"))

m8_cor <- racircal(data = as.data.frame(m8), caldata = as.data.frame(e8), mincut = 10, maxcut = 1161)
m8_cor <- m8_cor %>% filter(Cicor<1400, Cicor>0)

plot(m8_cor$Cicor, m8_cor$Acor, pch=20)
plot(m8_cor$Cicor, m8_cor$Acor, pch=20, ylim=c(-2,20), xlim=c(0,1100))
points(m8$Ci, m8$A, pch=20, col="red")
points(e8$Ci, e8$A, pch=20, col="blue")

fit8b <- fitaci(m8_cor, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                        PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='bilinear')
fit8d <- fitaci(m8_cor, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor",
                                        PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='default')

plot(fit8b, ylim=c(0,15))
fit8b$pars
fit8b$RMSE

# correct
poly.vertices <- data.frame(x=c(225, 225, -1, -1, 650, 650, 1000, 1000), y=c(-5, 22, 22,23, 23, -4, -4, -5))
polygon(poly.vertices$x, poly.vertices$y, border="maroon", col=NA)
outside <- (pip2d(as.matrix(poly.vertices), as.matrix(m8_cor %>% select(Cicor, Acor))) < 0)

plot(m8_cor$Cicor, m8_cor$Acor, pch=20, col="black")
points(m8_cor$Cicor[outside==F], m8_cor$Acor[outside==F], pch=20, col="green")
points(m8_cor$Cicor[outside==T], m8_cor$Acor[outside==T], pch=20, col="red")

m8_cor2 <- m8_cor %>% filter(! outside)

fit8b2 <- fitaci(m8_cor2, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                          PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='bilinear')
fit8d2 <- fitaci(m8_cor2, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor",
                                          PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='default')

fit8d2$pars
fit8b2$pars

plot(fit8d2)
plot(fit8b2)

fit8d2$RMSE
fit8b2$RMSE

#M9 P465.2-1
m9 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P465.2", "-", as.character(info[info$Plot=="P465.2","adult 1"]), ".csv"))
e9 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P465.2", "-", as.character(info[info$Plot=="P465.2","empty"]), ".csv"))

m9_cor <- racircal(data = as.data.frame(m9), caldata = as.data.frame(e9), mincut = 10, maxcut =950)
m9_cor <- m9_cor %>% filter(Cicor<1400, Cicor>0)

plot(m9_cor$Cicor, m9_cor$Acor, pch=20)
plot(m9_cor$Cicor, m9_cor$Acor, pch=20, ylim=c(-2,20), xlim=c(0,1100))
points(m9$Ci, m9$A, pch=20, col="red")
points(e9$Ci, e9$A, pch=20, col="blue")

fit9b <- fitaci(m9_cor, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                        PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='bilinear')
fit9d <- fitaci(m9_cor, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor",
                                        PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='default')

plot(fit9b, ylim=c(0,10))
plot(fit9d, ylim=c(0,10))

fit9d$pars
fit9b$pars

fit9d$RMSE
fit9b$RMSE

# correct
poly.vertices <- data.frame(x=c(800, 800, 1300, 1300, -1, -1, 1000, 1000), y=c(7, -1, -1, 20, 10, 20, 20, -5))
polygon(poly.vertices$x, poly.vertices$y, border="maroon", col=NA)
outside <- (pip2d(as.matrix(poly.vertices), as.matrix(m8_cor %>% select(Cicor, Acor))) < 0)

plot(m8_cor$Cicor, m8_cor$Acor, pch=20, col="black")
points(m8_cor$Cicor[outside==F], m8_cor$Acor[outside==F], pch=20, col="green")
points(m8_cor$Cicor[outside==T], m8_cor$Acor[outside==T], pch=20, col="red")

m8_cor2 <- m8_cor %>% filter(! outside)

fit8b2 <- fitaci(m8_cor2, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                          PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='bilinear')
fit8d2 <- fitaci(m8_cor2, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor",
                                          PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='default')

fit8d2$pars
fit8b2$pars

plot(fit8d2)
plot(fit8b2)

fit8d2$RMSE
fit8b2$RMSE

# Other packages?
library(photosynthesis)

m1_cor$PARi1 <- rep(1800, nrow(m1_cor))  
m1_cor$PARi2 <- rep(1000, nrow(m1_cor)) 
m1_cor$PARi3 <- rep(3500, nrow(m1_cor)) 
m1_cor$g_mc1 <- rep(10, nrow(m1_cor)) 
m1_cor$g_mc2 <- rep(100, nrow(m1_cor)) 
m1_cor$g_mc3 <- rep(1000000, nrow(m1_cor)) 
m1_cor$g_mc4 <- rep(0.001, nrow(m1_cor)) 

m1_cor$T_leaf <- m1_cor$Tleaf+273.15

a1 <- fit_acianova_V_ad_response(data = m1_cor, varnames = list(A_net = "Acor", T_leaf = "T_leaf", C_i = "Cicor", 
                                            PPFD = "PARi1", g_mc = "g_mc1"))
a2 <- fit_aci_response(data = m1_cor, varnames = list(A_net = "Acor", T_leaf = "T_leaf", C_i = "Cicor", 
                                                     PPFD = "PARi", g_mc = "g_mc2"))
a3 <- fit_aci_response(data = m1_cor, varnames = list(A_net = "Acor", T_leaf = "T_leaf", C_i = "Cicor", 
                                                      PPFD = "PARi", g_mc = "g_mc3"))
a4 <- fit_aci_response(data = m1_cor, varnames = list(A_net = "Acor", T_leaf = "T_leaf", C_i = "Cicor", 
                                                      PPFD = "PARi", g_mc = "g_mc4"))

a5 <- fit_aci_response(data = m1_cor, varnames = list(A_net = "Acor", T_leaf = "T_leaf", C_i = "Cicor", 
                                                      PPFD = "PARi2", g_mc = "g_mc2"))
a6 <- fit_aci_response(data = m1_cor, varnames = list(A_net = "Acor", T_leaf = "T_leaf", C_i = "Cicor", 
                                                      PPFD = "PARi3", g_mc = "g_mc2"))

a1$`Fitted Parameters`
a2$`Fitted Parameters`
a3$`Fitted Parameters`
a4$`Fitted Parameters`
a5$`Fitted Parameters`
a6$`Fitted Parameters`

# So: g_mc does not have any effect (probably around 0.15)
# PPFD does have an effect!

m3_cor$PARi1 <- rep(1800, nrow(m3_cor))  
m3_cor$g_mc2 <- rep(100, nrow(m3_cor)) 
m3_cor$T_leaf <- m3_cor$Tleaf+273.15

m3_t <- fit_aci_response(data = m3_cor, varnames = list(A_net = "Acor", T_leaf = "T_leaf", C_i = "Cicor", 
                                                        PPFD = "Qabs", g_mc = "g_mc2"))

m4_cor$PARi1 <- rep(1800, nrow(m4_cor))  
m4_cor$g_mc2 <- rep(100, nrow(m4_cor)) 
m4_cor$T_leaf <- m4_cor$Tleaf+273.15

m4_t <- fit_aci_response(data = m4_cor, varnames = list(A_net = "Acor", T_leaf = "T_leaf", C_i = "Cicor", 
                                                        PPFD = "Qabs", g_mc = "g_mc2"))

m3_t$`Fitted Parameters`
m4_t$`Fitted Parameters`

# Harder fits
m5_cor$PARi1 <- rep(1800, nrow(m5_cor))  
m5_cor$g_mc2 <- rep(100, nrow(m5_cor)) 
m5_cor$T_leaf <- m5_cor$Tleaf+273.15

m5_t <- fit_aci_response(data = m5_cor, varnames = list(A_net = "Acor", T_leaf = "T_leaf", C_i = "Cicor", 
                                                        PPFD = "PARi1", g_mc = "g_mc2"))

m6_cor$PARi1 <- rep(1800, nrow(m6_cor))  
m6_cor$g_mc2 <- rep(100, nrow(m6_cor)) 
m6_cor$T_leaf <- m6_cor$Tleaf+273.15

m6_t <- fit_aci_response(data = m6_cor, varnames = list(A_net = "Acor", T_leaf = "T_leaf", C_i = "Cicor", 
                                                        PPFD = "Qabs", g_mc = "g_mc2"))

m5_t$`Fitted Parameters`
m6_t$`Fitted Parameters`


m6_cor2$PARi1 <- rep(1800, nrow(m6_cor2))  
m6_cor2$g_mc2 <- rep(100, nrow(m6_cor2)) 
m6_cor2$T_leaf <- m6_cor2$Tleaf+273.15

m6_t2 <- fit_aci_response(data = m6_cor2, varnames = list(A_net = "Acor", T_leaf = "T_leaf", C_i = "Cicor", 
                                                        PPFD = "Qabs", g_mc = "g_mc2"))

m6_t2$`Fitted Parameters`

sensit <- analyze_sensitivity(data = data.example[data.example$Q_2 == 1500, ],
                            funct = fit_aci_response,
                            varnames = list(A_net = "A",
                                            T_leaf = "T_leaf",
                                            C_i = "Ci",
                                            PPFD = "Qin"),
                            useg_mct = TRUE,
                            test1 = "gamma_star25",
                            element_out = 1,
                            test2 = "g_mc25",
                            fitTPU = TRUE,
                            Ea_gamma_star = 0,
                            Ea_g_mc = 0,
                            values1 = seq(from = 20,
                                          to = 60,
                                          by = 4),
                            values2 = seq(from = 0.2,
                                          to = 2,
                                          by = 0.1))


data.example <- read.csv(system.file("extdata", "A_Ci_Q_data_1.csv",
                     package = "photosynthesis"))
data.example$Q_2 <- as.factor((round(data.example$Qin, digits = 0)))
# Convert leaf temperature to K
data.example$T_leaf <- data.example$Tleaf + 273.15
# Run a sensitivity analysis on gamma_


ggplot(sensit, aes(x = gamma_star25, y = g_mc25, z = V_cmax)) +
geom_tile(aes(fill = V_cmax)) +
labs(
x = expression(Gamma * "*"[25] ~ "(" * mu * mol ~ mol^
{
-1
} * ")"),
y = expression(g[m][25] ~ "(" * mu * mol ~ m^{
-2
} ~ s^{
-1
} ~ Pa^
{
-1
} * ")")
) +
scale_fill_distiller(palette = "Greys") +
geom_contour(colour = "Black", size = 1) +
theme_bw()

plot(a1$Data$Ca, a1$Data$Acor)
identify(a1$Data$Ca, a1$Data$Acor)

plot(a1$Data$Ca[-c(58,85)], a1$Data$Acor[-c(58,85)])



# Plot with fitaci
# Remove with identify
# 

plot(a1$Data$Cicor, a1$Data$Acor)
indices <- identify(a1$Data$Cicor, a1$Data$Acor)
plot(a1$Data$Cicor[-indices], a1$Data$Acor[-indices])

a1 <- fit_aci_response(data = m1_cor, varnames = list(A_net = "Acor", T_leaf = "T_leaf", C_i = "Cicor", 
                                                      PPFD = "PARi1", g_mc = "g_mc1"))

a1$`Fitted Parameters`

## TEST the procedure:

m5_t

plot(m5_cor$Cicor, m5_cor$Acor)
indices <- identify(m5_cor$Cicor, m5_cor$Acor)
plot(m5_cor$Cicor[-indices], m5_cor$Acor[-indices])

m5_cor2 <- m5_cor %>% 
  slice(-indices)

plot(m5_cor2$Cicor, m5_cor2$Acor)
indices <- identify(m5_cor2$Cicor, m5_cor2$Acor)
plot(m5_cor2$Cicor[-indices], m5_cor2$Acor[-indices])

m5_cor2 <- m5_cor2 %>% 
  slice(-indices)

m5_t_cut <- fit_aci_response(data = m5_cor2, varnames = list(A_net = "Acor", T_leaf = "T_leaf", C_i = "Cicor", 
                                                             PPFD = "PARi1", g_mc = "g_mc2"))                                      

m5_t_cut

m5_t$`Fitted Parameters`

m5_t_cut$`Fitted Parameters`

#m6

m6_t

plot(m6_cor$Cicor, m6_cor$Acor)
indices <- identify(m6_cor$Cicor, m6_cor$Acor)
plot(m6_cor$Cicor[-indices], m6_cor$Acor[-indices])

m6_cor2 <- m6_cor %>% 
  slice(-indices)

plot(m6_cor2$Cicor, m6_cor2$Acor)
indices <- identify(m6_cor2$Cicor, m6_cor2$Acor)
plot(m6_cor2$Cicor[-indices], m6_cor2$Acor[-indices])

m6_cor2 <- m6_cor2 %>% 
  slice(-indices)

m6_t_cut <- fit_aci_response(data = m6_cor2, varnames = list(A_net = "Acor", T_leaf = "T_leaf", C_i = "Cicor", 
                                                             PPFD = "PARi1", g_mc = "g_mc2"))                                      

m6_t_cut

m6_t$`Fitted Parameters`

m6_t_cut$`Fitted Parameters`


#### X Get back to photosynthesis package: -----
library(dplyr)
library(photosynthesis)

## Get data:
# meas1: perfect, meas2: good, meas3: good with cut at end
# meas4: late transition point, meas5: tricky, meas6: bad/unusable(?)

m1 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P460.2", "-", as.character(info[info$Plot=="P460.2","adult 2"]), ".csv"))
e1 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P460.2", "-", as.character(info[info$Plot=="P460.2","empty"]), ".csv"))
m2 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P465.1", "-", as.character(info[info$Plot=="P465.1","adult 2"]), ".csv"))
e2 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P465.1", "-", as.character(info[info$Plot=="P465.1","empty"]), ".csv"))
m3 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P460.1", "-", as.character(info[info$Plot=="P460.1","adult 5"]), ".csv"))
e3 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P460.1", "-", as.character(info[info$Plot=="P460.1","empty"]), ".csv"))
m4 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P470.3", "-", as.character(info[info$Plot=="P470.3","adult 4"]), ".csv"))
e4 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P470.3", "-", as.character(info[info$Plot=="P470.3","empty"]), ".csv"))
m5 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P465.1", "-", as.character(info[info$Plot=="P465.1","sapling 5"]), ".csv"))
e5 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P465.1", "-", as.character(info[info$Plot=="P465.1","empty"]), ".csv"))
m6 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P470.3", "-", as.character(info[info$Plot=="P470.3","adult 1"]), ".csv"))
e6 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                      "P470.3", "-", as.character(info[info$Plot=="P470.3","empty"]), ".csv"))

##Correct measurements based on empty measurements:
m1_cor <- racircal(data = as.data.frame(m1), caldata = as.data.frame(e1), mincut = 10, maxcut = 1080)
m1_cor <- m1_cor %>% filter(Cicor<1400, Cicor>0)

plot(m1_cor$Cicor, m1_cor$Acor, pch=20, ylim=c(0,20), xlim=c(0,1100))
points(m1$Ci, m1$A, pch=20, col="red")
points(e1$Ci, e1$A, pch=20, col="blue")

## Change temperature to K for photosynthesis package, and set Qin and mesophyl conductance levels:

# So, tests reveal that g_mc does not have any effect (probably around 0.15)
# PPFD does have an effect! PPFD should be the incident irradiance in umol m-2 s-1, 
# assumed to be the absorbed irradiance, so either Qin or Qabs


m1_final <- m1_cor %>% 
  mutate(Tleaf_K = Tleaf+273.15)
  #Likely need to input mesophys conductance as well
 


## Automate deleting outliers: 
# 
correct_aci <- function(plot, tree) {
  ## Read in data
  data.split <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                                plot, "-", as.character(info[info$Plot==plot, tree]), ".csv"))
  data.empty <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                                plot, "-", as.character(info[info$Plot==plot,"empty"]), ".csv"))
  
  # Get rid of some extreme outliers (Q1/Q3 +- 1.5 x IQR)
  A_range <- c(quantile(data.split$A)[[2]] - IQR(data.split$A) *1.5, quantile(data.split$A)[[4]] + IQR(data.split$A) *1.5)
  data.split <- data.split[data.split$A >= A_range[1] & data.split$A <= A_range[2],]
  
  data.split <- filter(data.split, Ci < 2000, Ci > -100)
  ## Iterate until correct range of Ci is used:
  Cmin <- 10
  
  # Iterate at low resolution
  data.corrected <- lapply(c(800, 900, 1000, 1100, 1200), FUN=function(x){
    racircal(data = as.data.frame(data.split), caldata = as.data.frame(data.empty), 
             mincut = Cmin, maxcut = x)
  })
  sapply(length(data.corrected):1, FUN=function(x, colors=c("orange", "darkblue", "green", "red", "black")){
    if(x == length(data.corrected)){plot(data.corrected[[x]]$Cicor, data.corrected[[x]]$Acor, pch=20)
      } else{points(data.corrected[[x]]$Cicor, data.corrected[[x]]$Acor, col=colors[x], pch=20)}
  })
  legend("topleft", legend = as.character(c(1200, 1100, 1000, 900, 800)), 
         col = c("black", "red", "green", "darkblue", "orange"), pch=20)
  
  c.high <- as.numeric(readline("What is the higher colour to explore? "))
  c.low <- as.numeric(readline("What is the lower colour to explore? "))
  
  # Iterate at higher resolution
  data.corrected2 <- lapply(seq(c.high, c.low, by=-25), FUN=function(x){
    racircal(data = as.data.frame(data.split), caldata = as.data.frame(data.empty), 
             mincut = Cmin, maxcut = x)
  })
  print(length(data.corrected2))
  sapply(1:length(data.corrected2), FUN=function(x, colors=c("black", "red", "green", "darkblue", "orange")){
    if(x == 1){plot(data.corrected2[[x]]$Cicor, data.corrected2[[x]]$Acor, pch=20)
    } else{points(data.corrected2[[x]]$Cicor, data.corrected2[[x]]$Acor, col=colors[x], pch=20)}
  })
  legend("topleft", legend = as.character(seq(c.high, c.low, by=-25)), 
         col = c("black", "red", "green", "darkblue", "orange"), pch=20)

  Cmax <- as.numeric(readline("What is the final upper limit to consider? "))
  
  #TODO add a title (=argument) to racircal to get nice titles in the figures
  ## Correct the measurements using the empty measurement on the correct range: 
  final.corrected <- racircal(data = as.data.frame(data.split), caldata = as.data.frame(data.empty), 
           mincut = Cmin, maxcut = Cmax)
  
  treenumber <- as.numeric(substring(tree, nchar(tree), nchar(tree))) + 
    ifelse(substring(tree, 1, 1) == "a", 0, 5)
  
  ## Write out the corrected measurement: 
  write_csv(final.corrected, paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split-corrected/",
                  plot, "-", as.character(treenumber), ".csv"))
}

lapply(colnames(info[3:12]), function(x){correct_aci("P480.3", x)})


print("test")
## Run the actual analyses to extract parameters

curves <- Sys.glob("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split-corrected/*.csv")

for (curve in curves) {
  data <- read_csv(curve) %>% 
    mutate(Tleaf_K = Tleaf+273.15)
  tree <- substring(basename(curve), 1, nchar(basename(curve)) - 4)
  
  aci.fit <- fit_aci_response(data = as.data.frame(data), varnames = list(A_net = "Acor", T_leaf = "Tleaf_K", C_i = "Cicor",
                                                                   PPFD = "Qabs"))
  saveRDS(aci.fit,
          paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Fits/",
                 tree,
                 ".RDS"))
}


results <- Sys.glob("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Fits/*.RDS")

lapply(results, FUN=function(r){
  result <- readRDS(r)
  
  tree <- substring(basename(r), 1, nchar(basename(r)) - 4)
  
  pdf(file=paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Fits/Figures/",
                  tree, ".pdf"))
  print(result$Plot)
  dev.off()
})

#### Clean up the curves to hopefully improve some fitting: 

clean_aci <- function(tree, cleaned.before=F) {
  if(cleaned.before){
    curve <- readRDS(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Fits/",
                            tree, "-cc.RDS"))
    data <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split-corrected-cleaned/",
                            tree, "-cc.csv"))
  } else {
  curve <- readRDS(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Fits/",
                  tree, ".RDS"))
  data <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split-corrected/",
                         tree, ".csv"))
  }
  print(curve$Plot)
  
  data.clean <- data
  continue <- as.logical(readline("Iterate? (T/F) "))
  if(! continue) {
    # do nothing (?)
  } else {
  
  while(continue){
    plot(data.clean$Cicor, data.clean$Acor)
    outliers <- identify(data.clean$Cicor, data.clean$Acor)
    data.clean <- data.clean %>% 
      slice(-outliers)
    plot(fitaci(data.clean, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                           PPFD = "Qabs", Rd ="Rd"), Tcorrect=T, fitmethod='bilinear'))
    continue <- as.logical(readline("Re-iterate? (T/F) "))
  }
  
  write_csv(data.clean, paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split-corrected-cleaned/",
                  tree, "-cc.csv"))
    }
}

clean_aci("P475.3-10", cleaned.before = T)




### Try fitting again with cleaned datasets: 

#curves <- Sys.glob("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split-corrected-cleaned/*.csv")

curves <- "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split-corrected-cleaned/P475.3-10-cc.csv"

for (curve in curves) {
  data <- read_csv(curve) %>% 
    mutate(Tleaf_K = Tleaf+273.15)
  tree <- substring(basename(curve), 1, nchar(basename(curve)) - 7)
  
  aci.fit <- fit_aci_response(data = as.data.frame(data), varnames = list(A_net = "Acor", T_leaf = "Tleaf_K", C_i = "Cicor",
                                                                          PPFD = "Qabs")) #This should be absorbed irradiance
  
  pdf(file=paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Fits/Figures/",
                  tree, "-cc.pdf"))
  print(aci.fit$Plot)
  dev.off()
  
  saveRDS(aci.fit,
          paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Fits/",
                 tree,
                 "-cc.RDS"))
}
  

#### Try to fix the really bad ones: 
tree <- "P475.3-10"
clean_aci("P475.3-10", cleaned.before = T)
data.tocut <- read_csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split-corrected-cleaned/P475.3-10-cc.csv")
plot(data.tocut$Cicor, data.tocut$Acor, pch=20)

data.cut <- data.t0cut %>% 
  filter(Cicor < 300 | Cicor > 600)

plot(data.cut$Cicor, data.cut$Acor, pch=20)
plot(fitaci(data.cut, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                        PPFD = "Qabs", Rd ="Rd"), Tcorrect=T, fitmethod='bilinear'))
write_csv(data.cut, paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split-corrected-cleaned/",
                             tree, "-cc.csv"))


#### Now let's look at the empty measurements:

info$Plot

e <- lapply(info$Plot, function(x) read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                                                   x, "-", as.character(info[info$Plot==x,"empty"]), ".csv")))

# All seem relatively good! P460-1 deviates most from other measurements and P480-1 seem a bit off. 
# There is a second empty measurement of P480.1, let's see if that helps
e10 <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                       "P480.1-12", ".csv"))

plot(e[[13]]$Ci, e[[13]]$A, ylim=c(-4, 4), xlim=c(-10, 1400), pch=20)
points(e10$Ci, e10$A, ylim=c(-4, 4), pch=20, col="orange")


t <- racircal(data = as.data.frame(read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                                              "P480.1", "-", as.character(info[info$Plot=="P480.1","sapling 3"]), ".csv"))),
         caldata = as.data.frame(e[[13]]), mincut = 10, maxcut = 1000)

t2 <- racircal(data = as.data.frame(read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                                                   "P480.1", "-", as.character(info[info$Plot=="P480.1","sapling 3"]), ".csv"))),
              caldata = as.data.frame(e10), mincut = 10, maxcut = 1000)

plot(t$Cicor, t$Acor, pch=20, xlim=c(0, 1100), ylim=c(-2, 15))
points(t2$Cicor, t2$Acor, pch=20, col="red")

fitaci(t2, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                            PPFD = "PARi", Rd ="Rd"), Tcorrect=T, fitmethod='bilinear')

# Does not seem to change much to the corrected curves. 

#### The extra measurement: -----

correct_aci_extra <- function(data.split, data.empty) {
  
  # Get rid of some extreme outliers (Q1/Q3 +- 1.5 x IQR)
  A_range <- c(quantile(data.split$A)[[2]] - IQR(data.split$A) *1.5, quantile(data.split$A)[[4]] + IQR(data.split$A) *1.5)
  data.split <- data.split[data.split$A >= A_range[1] & data.split$A <= A_range[2],]
  
  data.split <- filter(data.split, Ci < 2000, Ci > -100)
  ## Iterate until correct range of Ci is used:
  Cmin <- 10
  
  # Iterate at low resolution
  data.corrected <- lapply(c(800, 900, 1000, 1100, 1200), FUN=function(x){
    racircal(data = as.data.frame(data.split), caldata = as.data.frame(data.empty), 
             mincut = Cmin, maxcut = x)
  })
  sapply(length(data.corrected):1, FUN=function(x, colors=c("orange", "darkblue", "green", "red", "black")){
    if(x == length(data.corrected)){plot(data.corrected[[x]]$Cicor, data.corrected[[x]]$Acor, pch=20)
    } else{points(data.corrected[[x]]$Cicor, data.corrected[[x]]$Acor, col=colors[x], pch=20)}
  })
  legend("topleft", legend = as.character(c(1200, 1100, 1000, 900, 800)), 
         col = c("black", "red", "green", "darkblue", "orange"), pch=20)
  
  c.high <- as.numeric(readline("What is the higher colour to explore? "))
  c.low <- as.numeric(readline("What is the lower colour to explore? "))
  
  # Iterate at higher resolution
  data.corrected2 <- lapply(seq(c.high, c.low, by=-25), FUN=function(x){
    racircal(data = as.data.frame(data.split), caldata = as.data.frame(data.empty), 
             mincut = Cmin, maxcut = x)
  })
  print(length(data.corrected2))
  sapply(1:length(data.corrected2), FUN=function(x, colors=c("black", "red", "green", "darkblue", "orange")){
    if(x == 1){plot(data.corrected2[[x]]$Cicor, data.corrected2[[x]]$Acor, pch=20)
    } else{points(data.corrected2[[x]]$Cicor, data.corrected2[[x]]$Acor, col=colors[x], pch=20)}
  })
  legend("topleft", legend = as.character(seq(c.high, c.low, by=-25)), 
         col = c("black", "red", "green", "darkblue", "orange"), pch=20)
  
  Cmax <- as.numeric(readline("What is the final upper limit to consider? "))
  
  ## Correct the measurements using the empty measurement on the correct range: 
  final.corrected <- racircal(data = as.data.frame(data.split), caldata = as.data.frame(data.empty), 
                              mincut = Cmin, maxcut = Cmax)
  
  treenumber <- as.numeric(substring(tree, nchar(tree), nchar(tree))) + 
    ifelse(substring(tree, 1, 1) == "a", 0, 5)
  
  ## Write out the corrected measurement: 
  write_csv(final.corrected, paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split-corrected/",
                                    "P475.1", "-", "e", ".csv"))
}

extra <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                         "P475.1", "-", as.character(5), ".csv"))
extra.empty <- read_csv(paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split/",
                               "P475.1", "-", as.character(info[info$Plot=="P475.1","empty"]), ".csv"))

correct_aci_extra(extra, extra.empty)

curve <- "/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split-corrected/P475.1-e.csv"

data <- read_csv(curve) %>% 
  mutate(Tleaf_K = Tleaf+273.15)

aci.fit <- fit_aci_response(data = as.data.frame(data), varnames = list(A_net = "Acor", T_leaf = "Tleaf_K", C_i = "Cicor",
                                                                        PPFD = "Qabs"))
result <- saveRDS(aci.fit,
        paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Fits/",
               "P475.1-e",
               ".RDS"))

pdf(file=paste0("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Fits/Figures/",
                "P475.1-e", ".pdf"))
print(result$Plot)
dev.off()


clean_aci("P475.1-e", cleaned.before = F)


### Test the effect of mesophyll conductance: ----
library(tidyverse)
library(plantecophys)
library(racir)
library(photosynthesis)
# read in curve  P460.2-2
curve <- read.csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split-corrected/P460.2-2.csv")
#curve <- read.csv("/Users/lukas/Library/CloudStorage/OneDrive-UniversitedeMontreal/Data/Ch2/Processed/Gas-Measurements/Split-corrected-cleaned/P460.2-2-cc.csv")

plot(curve$Cicor, curve$Acor, pch = 20)

curve <- curve %>% 
  mutate(
    Tleaf_K = Tleaf + 273.15,
    g_mc1 = 0.4,
    g_mc3 = 1000000)


## fits
fita.d0 <- fitaci(curve, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                        PPFD = "Qabs", Rd ="Rd"), Tcorrect=T, fitmethod='default')
fita.b0 <- fitaci(curve, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                       PPFD = "Qabs", Rd ="Rd"), Tcorrect=T, fitmethod='bilinear')

fita.d1 <- fitaci(curve, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                        PPFD = "Qabs", Rd ="Rd"), gmeso = 0.4, Tcorrect=T, fitmethod='default')
fita.b1 <- fitaci(curve, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                        PPFD = "Qabs", Rd ="Rd"), gmeso = 0.4, Tcorrect=T, fitmethod='bilinear')

fita.d2 <- fitaci(curve, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                        PPFD = "Qabs", Rd ="Rd"), gmeso = 0.087, Tcorrect=T, fitmethod='default')
fita.b2 <- fitaci(curve, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                        PPFD = "Qabs", Rd ="Rd"), gmeso = 0.087, Tcorrect=T, fitmethod='bilinear')

fita.d3 <- fitaci(curve, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                        PPFD = "Qabs", Rd ="Rd"), gmeso = 1000000, Tcorrect=T, fitmethod='default')
fita.b3 <- fitaci(curve, varnames = list(ALEAF = "Acor", Tleaf = "Tleaf", Ci = "Cicor", 
                                        PPFD = "Qabs", Rd ="Rd"), gmeso = 1000000, Tcorrect=T, fitmethod='bilinear')

plot(fita.d0)
plot(fita.b0)
plot(fita.d1)
plot(fita.b1)
plot(fita.d2)
plot(fita.b2)
plot(fita.d3)
plot(fita.b3)

fita.d0$pars
fita.b0$pars
fita.d1$pars
fita.b1$pars
fita.b2$pars
fita.b3$pars



fitp0 <- fit_aci_response(data = as.data.frame(curve), varnames = list(A_net = "Acor", T_leaf = "Tleaf_K", C_i = "Cicor", 
                                                 PPFD = "Qabs")) 
fitp0$`Fitted Parameters`

fitp1t <- fit_aci_response(data = as.data.frame(curve), varnames = list(A_net = "Acor", T_leaf = "Tleaf_K", C_i = "Cicor", 
                                                                       PPFD = "Qabs", g_mc = "g_mc1"), useg_mct = T) 
fitp1t$`Fitted Parameters`

fitp2 <- fit_aci_response(data = as.data.frame(curve), varnames = list(A_net = "Acor", T_leaf = "Tleaf_K", C_i = "Cicor", 
                                                                       PPFD = "Qabs", g_mc = "g_mc2"), useg_mc = T)
fitp2$`Fitted Parameters`

fitp3 <- fit_aci_response(data = as.data.frame(curve), varnames = list(A_net = "Acor", T_leaf = "Tleaf_K", C_i = "Cicor", 
                                                                       PPFD = "Qabs", g_mc = "g_mc3"), useg_mc = T)
fitp3$`Fitted Parameters`

fitp0$Plot
fitp1$Plot
fitp1t$Plot
fitp3$Plot

fitp2$`Fitted Parameters`
fitp3$`Fitted Parameters`
fitp1$`Fitted Parameters`
