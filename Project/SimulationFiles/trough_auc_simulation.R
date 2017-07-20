# Trough-AUC Simulation
# Simulation exercise to explore the range of AUCs that are generated when
# various individual trough concentrations are 50 ng/mL (i.e., target trough
# concentration for sunitinib)
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

# Plot matrix of random effects
  plotmatrix <- ggpairs(data = pk.ETA.matrix,
    columnLabels = c("Parent Clearance","Parent Central Volume",
    "Metabolite Clearance","Metabolite Central Volume"),
    lower = list(
      continuous = wrap("points",alpha = 0.3,colour = "blue",shape = 1)
    ),
    upper = list(
      continuous = wrap("cor",colour = "black")
    ),
    diag = list(
      continuous = wrap("densityDiag",colour = "red")
    )
  )
  print(plotmatrix)
  ggsave(plot = plotmatrix,
    filename = "pk_model_correlations.png",
    dpi = 300,height = 9,width = 9)

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
  last.pk.data <- pk.data[pk.data$time == 4*7*24,]

  plotobj1 <- NULL
  plotobj1 <- ggplot(last.pk.data)
  plotobj1 <- plotobj1 + geom_point(aes(x = IPRE,y = AUC24),
    colour = "blue",shape = 1,size = 2,alpha = 0.3)
  plotobj1 <- plotobj1 + geom_vline(aes(xintercept = 0.05),
    linetype = "dashed")
  plotobj1 <- plotobj1 + geom_vline(aes(xintercept = 0.1),
    linetype = "dashed")
  plotobj1 <- plotobj1 + scale_y_continuous("24-hour AUC (mg*h/L)",
    lim = c(0,6.5))
  plotobj1 <- plotobj1 + scale_x_continuous("Trough Concentration (mg/L)",
    lim = c(0,0.3))
  print(plotobj1)

# Plot IPRE in the last dosing interval
  last.dose.pk.data <- pk.data[pk.data$time >= 4*7*24-24 &
    pk.data$time <= 4*7*24,]
  interval.times <- unique(last.dose.pk.data$time)
  unique.ID <- rep(1:ntotal,times = length(interval.times)) %>% sort
  last.dose.pk.data$uID <- unique.ID

  last.dose.pk.data$WTf <- "< 70 kg"
  last.dose.pk.data$WTf[last.dose.pk.data$WT >= 70] <- ">= 70 kg"
  last.dose.pk.data$WTf <- as.factor(last.dose.pk.data$WTf)

  plotobj2 <- NULL
  plotobj2 <- ggplot(last.dose.pk.data)
  plotobj2 <- plotobj2 + geom_line(aes(x = time,y = IPRE,group = uID,
    colour = WTf),
    alpha = 0.3)
  plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 0.05),
    linetype = "dashed")
  plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 0.1),
    linetype = "dashed")
  plotobj2 <- plotobj2 + scale_y_continuous("Sunitinib Concentration (mg/L)",
    lim = c(0,0.3))
  plotobj2 <- plotobj2 + scale_x_continuous("Time Since First Dose (hours)",
    breaks = seq(from = 648,to = 672,by = 4))
  plotobj2 <- plotobj2 + theme(legend.position = "none")
  print(plotobj2)

# # Subset the data to only include trough concentrations within the therapeutic
# # range (0.05 - 0.1 mg/L)
#   thera.pk.data <- last.pk.data[last.pk.data$IPRE >= 0.05 &
#     last.pk.data$IPRE <= 0.1,]
#   length(unique(thera.pk.data$ID))
# # What is the minimum AUC when trough is in the therapeutic range?
#   min.AUC24 <- min(thera.pk.data$AUC24)
#   print(min.AUC24)
# # What is the maximum AUC when trough is in the therapeutic range?
#   max.AUC24 <- max(thera.pk.data$AUC24)
#   print(max.AUC24)
# # Round IPRE to the nearest 0.01
#   thera.pk.data$rIPRE <- round(thera.pk.data$IPRE,digits = 2)
# # What are the range of AUCs when rIPRE == 0.05
#   summary.AUC24.lowIPRE <- summary.function(
#     thera.pk.data$AUC24[thera.pk.data$rIPRE == 0.05])
#   print(summary.AUC24.lowIPRE)
# # What are the range of AUCs when rIPRE == 0.1
#   summary.AUC24.highIPRE <- summary.function(
#     thera.pk.data$AUC24[thera.pk.data$rIPRE == 0.1])
#   print(summary.AUC24.highIPRE)

# ------------------------------------------------------------------------------
# Save plots and summary statistics
# Plot
  ggsave(plot = plotobj1,
    filename = "week4_auc_vs_trough_50mg.png",
    dpi = 300,height = 5,width = 7)
  ggsave(plot = plotobj2,
    filename = "week4_IPRE_vs_time_50mg.png",
    dpi = 300,height = 5,width = 7)
# # Summary statistics
#   summary.list <- list("AUC24 when IPRE == 0.05" = summary.AUC24.lowIPRE,
#     "AUC24 when IPRE == 0.1" = summary.AUC24.highIPRE)
#   capture.output(summary.list,file = "week4_auc_versus_trough_summary.txt")
