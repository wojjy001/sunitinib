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
  library(reshape)	# melt function
# Set the working directory
  dir <- "/Volumes/Prosecutor/sunitinib/PopPKPDModels/"
  setwd(dir)
# Source model code
  source("hansson_2013_biomarker_tumour_model.R")

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
  nsim <- 50
  ID.seq <- 1:nsim
# Time sequence
  times <- seq(from = 0,to = 2016,by = 24)
# Dosing and dosing times
  DOSE <- 50
  on.times <- c(seq(from = 24,to = 672,by = 24),
    seq(from = 1032,to = 1680,by = 24))
# Random effects
  ETA.matrix <- mvrnorm(nsim,mu = rep(0,times = dim(OMEGA)[1]),OMEGA) %>%
  as.data.frame
  EPS.matrix <- mvrnorm(nsim*length(times),mu = rep(0,times = dim(SIGMA)[1]),SIGMA) %>% as.data.frame

# ------------------------------------------------------------------------------
# Input data frame for simulation
  ID.data <- data.frame(
    ID = ID.seq,
    DOSE = 0,
    ETA.matrix,
    EPS.matrix
  )
  names(ID.data)[c(3:dim(ID.data)[2])] <- c("ETAMRT23","ETABM0","ETADPSLO","ETABM02","ETABM03","ETABM0S","ETAMRTS",
  "ETADPSLOS","ETAIC50","ETAIC502","ETAIC503","ETAIC50S","ETAKG","ETAKSKIT","ETAKDRUG","EPSBASE")
  input.data <- lapply(ID.data,rep.int,times = length(times)) %>% as.data.frame
  input.data <- input.data[with(input.data,order(input.data$ID)),]
  input.data$time <- times
  input.data$DOSE[input.data$time %in% on.times] <- DOSE
  input.data$cmt <- 0
  input.data <- input.data[with(input.data,order(input.data$ID,input.data$time)),]

# ------------------------------------------------------------------------------
# Simulate
  sim.data <- mod %>% mrgsim(data = input.data) %>% as.data.frame

# ------------------------------------------------------------------------------
# Rearrange data frame for plotting
  melt.data <- melt(sim.data,id = c("ID","time","DOSE","CL"),
    measure = c("IPRE_VEGF","IPRE_VEGFR2","IPRE_VEGFR3","IPRE_SKIT","TUMOUR"))
# Plot
  axis.breaks <- c(seq(from = 0.1,to = 0.9,by = 0.1),
                    seq(from = 1,to = 9,by = 1),
                    seq(from = 10,to = 90,by = 10),
                    seq(from = 100,to = 1000,by = 100))
  plotobj <- NULL
  plotobj <- ggplot(melt.data)
  plotobj <- plotobj + stat_summary(aes(x = time/7/24,y = value),
    geom = "line",fun.y = median,colour = "red")
  plotobj <- plotobj + stat_summary(aes(x = time/7/24,y = value),
    geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
    fill = "red",alpha = 0.3)
  plotobj <- plotobj + scale_y_log10("\nBiomarker (pg/mL)",
    breaks = axis.breaks,labels = axis.breaks,lim = c(1,NA))
  plotobj <- plotobj + scale_x_continuous("Time (weeks)\n",
    breaks = seq(from = 0,to = 18,by = 2))
  plotobj <- plotobj + facet_wrap(~variable)
  plotobj
