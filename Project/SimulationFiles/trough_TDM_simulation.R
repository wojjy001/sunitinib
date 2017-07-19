# Trough TDM Simulation
# Simulation exercise to explore the range of AUCs that are generated when
# each individual has their trough concentrations to be target to 50 ng/mL
# ------------------------------------------------------------------------------
# Remove objects in workspace
  rm(list = ls(all = TRUE))
# Define global directory for project
  global.dir <- "/Volumes/Prosecutor/sunitinib/Project/"
# Source required files
  source(paste0(global.dir,"functions.R"))	# Universal functions file
  source(paste0(global.dir,"ModelFiles/yu_2015_sunitinib_model.R"))	# PK model
# Define simulation file directory
  sim.dir <- paste0(global.dir,"SimulationFiles/")
  source(paste0(sim.dir,"trough_auc_simulation_options.R"))	# Simulation options
# Define output directory
  output.dir <- paste0(global.dir,"Output/TroughAUCSimulation/")
# Set the working directory
  setwd(output.dir)

# ------------------------------------------------------------------------------
# Generate random effect parameters
# Between subject variability
  if (ntotal > 1) {
    pk.ETA.matrix <- mvrnorm(ntotal,
      mu = rep(0,times = dim(pk.OMEGA)[1]),pk.OMEGA) %>%
      as.data.frame
    names(pk.ETA.matrix) <- c("ETACLP","ETAVCP","ETACLM","ETAVCM")
  } else {
    pk.ETA.matrix <- data.frame(ETACLP = 0,ETAVCP = 0,ETACLM = 0,ETAVCM = 0)
  }

# ------------------------------------------------------------------------------
# Perform PK simulations
# Input data frame for simulation
  pk.ID.data <- data.frame(SIM = SIM.seq,ID = ID.seq,WT,pk.ETA.matrix)
  input.pk.data <- lapply(pk.ID.data,rep.int,times = length(pk.times)) %>%
    as.data.frame
  input.pk.data <- input.pk.data[with(input.pk.data,
    order(input.pk.data$SIM,input.pk.data$ID)),]
  input.pk.data$time <- pk.times
  input.pk.data$cmt <- 1
  input.pk.data$amt <- dose
  input.pk.data$evid <- 0
  input.pk.data$evid[input.pk.data$time %in% on.times] <- 1
# Simulate
  pk.data <- pk.mod %>% mrgsim(data = input.pk.data,
    carry.out = c("SIM","amt")) %>%
    as.data.frame
# Calculate 24-hour AUCs (of the combined parent and metabolite)
  pk.data <- ddply(pk.data, .(SIM,ID), auc24.function, .progress = "text")

# ------------------------------------------------------------------------------
# Plot AUC24 versus trough concentration on the last day of 4 weeks treatment
  last.pk.data <- pk.data[pk.data$time == max(pk.times),]

  plotobj1 <- NULL
  plotobj1 <- ggplot(last.pk.data)
  plotobj1 <- plotobj1 + geom_point(aes(x = IPRE,y = AUC24),
    colour = "blue",shape = 1,size = 2,alpha = 0.3)
  plotobj1 <- plotobj1 + geom_vline(aes(xintercept = 0.05),
    linetype = "dashed")
  plotobj1 <- plotobj1 + geom_vline(aes(xintercept = 0.1),
    linetype = "dashed")
  plotobj1 <- plotobj1 + scale_y_continuous("24-hour AUC (mg*h/L)")
  plotobj1 <- plotobj1 + scale_x_continuous("Trough Concentration (mg/L)")
  print(plotobj1)

# Subset the data to only include trough concentrations within the therapeutic
# range (0.05 - 0.1 mg/L)
  thera.pk.data <- last.pk.data[last.pk.data$IPRE >= 0.05 &
    last.pk.data$IPRE <= 0.1,]
  length(unique(thera.pk.data$ID))
# What is the minimum AUC when trough is in the therapeutic range?
  min.AUC24 <- min(thera.pk.data$AUC24)
  print(min.AUC24)
# What is the maximum AUC when trough is in the therapeutic range?
  max.AUC24 <- max(thera.pk.data$AUC24)
  print(max.AUC24)
# Round IPRE to the nearest 0.01
  thera.pk.data$rIPRE <- round(thera.pk.data$IPRE,digits = 2)
# What are the range of AUCs when rIPRE == 0.05
  summary.AUC24.lowIPRE <- summary.function(
    thera.pk.data$AUC24[thera.pk.data$rIPRE == 0.05])
  print(summary.AUC24.lowIPRE)
# What are the range of AUCs when rIPRE == 0.1
  summary.AUC24.highIPRE <- summary.function(
    thera.pk.data$AUC24[thera.pk.data$rIPRE == 0.1])
  print(summary.AUC24.highIPRE)

# ------------------------------------------------------------------------------
# Save plot and summary statistics
# Plot
  ggsave(plot = plotobj1,
    filename = "week4_auc_versus_trough.png",
    dpi = 300)
# Summary statistics
  summary.list <- list("AUC24 when IPRE == 0.05" = summary.AUC24.lowIPRE,
    "AUC24 when IPRE == 0.1" = summary.AUC24.highIPRE)
  capture.output(summary.list,file = "week4_auc_versus_trough_summary.txt")
