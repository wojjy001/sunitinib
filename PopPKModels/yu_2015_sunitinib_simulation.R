# Yu H, Steeghs N, Kloth JS, Wit D, Hasselt J, Erp NP, Beijnen JH, Schellens JH,
# Mathijssen RH, Huitema AD. Integrated semi‐physiological pharmacokinetic model
# for both sunitinib and its active metabolite SU12662. British journal of
# clinical pharmacology. 2015;79(5):809-19
# ------------------------------------------------------------------------------
# Remove objects in workspace
  rm(list = ls(all = TRUE))
# Load package libraries
  library(MASS)	# mvrnorm function
  library(ggplot2)  # Plotting
  library(grid)   #Plotting
  library(plyr)  # Split and rearrange data, ddply function
  library(dplyr)  # New plyr
  library(mrgsolve) # Metrum differential equation solver for pharmacometrics
# Source functions file
  source("/Volumes/Prosecutor/sunitinib/functions.R")
# Set the working directory
  dir <- "/Volumes/Prosecutor/sunitinib/PopPKModels/"
  setwd(dir)
# Source model code
  source("yu_2015_sunitinib_model.R")

# ------------------------------------------------------------------------------
# Set simulation options
# Number of individuals to simulate
  nsim <- 1000
  ID.seq <- 1:nsim
# Simulation times (hours)
  times <- seq(from = 0,to = 1000,by = 1)
# Dose
  dose <- 50	# mg
  freq <- 24	# frequency, hours
  dose.times <- seq(from = 0,to = 672,by = 24)
# Weight
  wt <- exp(log(rnorm(nsim,mean = 82.3,sd = 19.4)))
  wt[wt < 39] <- 39
  wt[wt > 157] <- 157
  plot(density(wt))
# Random effects
  ETA.matrix <- mvrnorm(nsim,mu = c(0,0,0,0),OMEGA) %>% as.data.frame

# ------------------------------------------------------------------------------
# Input data frame for simulation
  ID.data <- data.frame(
    ID = ID.seq,
    ETA.matrix,
    WT = wt
  )
  names(ID.data)[c(2,3,4,5)] <- c("ETACLP","ETAVCP","ETACLM","ETAVCM")
  input.data <- lapply(ID.data,rep.int,times = length(times)) %>% as.data.frame
  input.data <- input.data[with(input.data,order(input.data$ID)),]
  input.data$time <- times
  input.data$amt <- 0
  input.data$amt[input.data$time %in% dose.times] <- dose
  input.data$evid <- 0
  input.data$evid[input.data$time %in% dose.times] <- 1
  input.data$cmt <- 1

# ------------------------------------------------------------------------------
# Simulate
  sim.data <- mod %>% mrgsim(data = input.data) %>% as.data.frame

# ------------------------------------------------------------------------------
# Plot results
  plotobj1 <- NULL
  plotobj1 <- ggplot(sim.data)
  plotobj1 <- plotobj1 + stat_summary(aes(x = time,y = IPREP),
    geom = "line",fun.y = median, colour = "red",size = 1)	# Parent median
  plotobj1 <- plotobj1 + stat_summary(aes(x = time,y = IPREP),
    geom = "ribbon",fun.ymin = "CI50lo",fun.ymax = "CI50hi",
    fill = "red",alpha = 0.3)	# Parent confidence intervals
  plotobj1 <- plotobj1 + stat_summary(aes(x = time,y = IPREM),
    geom = "line",fun.y = median, colour = "blue",size = 1)	# Metabolite median
  plotobj1 <- plotobj1 + stat_summary(aes(x = time,y = IPREM),
    geom = "ribbon",fun.ymin = "CI50lo",fun.ymax = "CI50hi",
    fill = "blue",alpha = 0.3)	# Metabolite confidence intervals
  plotobj1 <- plotobj1 + scale_y_continuous("Concentration (ng/mL)\n",
    breaks = c(0,20,40,60,80),lim = c(0,80))
  plotobj1 <- plotobj1 + scale_x_continuous("\nTime (hours)")
  plotobj1
