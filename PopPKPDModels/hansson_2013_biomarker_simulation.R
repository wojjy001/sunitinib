# Hansson EK, Amantea MA, Westwood P, Milligan PA, Houk BE, French J, Karlsson MO,
# Friberg LE. PKPD Modeling of VEGF, sVEGFR-2, sVEGFR-3, and sKIT as Predictors
# of Tumor Dynamics and Overall Survival Following Sunitinib Treatment in GIST.
# CPT: pharmacometrics & systems pharmacology. 2013;2(11):1-9.
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
# Set the working directory
  dir <- "/Volumes/Prosecutor/sunitinib/PopPKPDModels/"
  setwd(dir)
# Source model code
  source("hansson_2013_biomarker_model.R")

# ------------------------------------------------------------------------------
# Functions
# 95% prediction interval functions
  CI95lo <- function(x) quantile(x,probs = 0.025)
  CI95hi <- function(x) quantile(x,probs = 0.975)
# Interquartile range
  CI50lo <- function(x) quantile(x,probs = 0.25)
  CI50hi <- function(x) quantile(x,probs = 0.75)
# Define a custom ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))

# ------------------------------------------------------------------------------
# Set simulation options
# Simulation seed for reproducible numbers
  set.seed(123456)
# Number of individuals to simulate
  nsim <- 10
  ID.seq <- 1:nsim
# Random effects
  ETA.matrix <- mvrnorm(nsim,mu = rep(0,times = dim(OMEGA)[1]),OMEGA) %>%
  as.data.frame
# Time sequence
  times <- seq(from = 0,to = 18*7*24,by = 24)
# Dosing and dosing times
  DOSE <- 50
  on.times <- c(seq(from = 24,to = 672,by = 24),
    seq(from = 1032,to = 1680,by = 24))

# ------------------------------------------------------------------------------
# Input data frame for simulation
  ID.data <- data.frame(
    ID = ID.seq,
    DOSE = 0,
    ETA.matrix
  )
  names(ID.data)[c(3:dim(ID.data)[2])] <- c("ETAMRT23","ETABM0","ETADPSLO","ETABM02","ETABM03","ETABM0S","ETAMRTS",
  "ETADPSLOS","ETAIC50","ETAIC502","ETAIC503","ETAIC50S")
  input.data <- lapply(ID.data,rep.int,times = length(times)) %>% as.data.frame
  input.data <- input.data[with(input.data,order(input.data$ID)),]
  input.data$time <- times
  input.data$DOSE[input.data$time %in% on.times] <- DOSE
  input.data$cmt <- 1

# ------------------------------------------------------------------------------
# Simulate
  sim.data <- mod %>% mrgsim(data = input.data) %>% as.data.frame

# ------------------------------------------------------------------------------
# Plot
  plotobj <- NULL
  plotobj <- ggplot(sim.data)
  plotobj <- plotobj + geom_line(aes(x = time/7/24,y = IPRE_VEGF))
  plotobj <- plotobj + scale_y_log10("\nVEGF (pg/mL)",
    breaks = c(0.1,1,10,100),labels = c(0.1,1,10,100))
  plotobj <- plotobj + scale_x_continuous("Time (weeks)\n",
    breaks = seq(from = 0,to = 18,by = 2))
  plotobj <- plotobj + facet_wrap(~ID)
  plotobj
