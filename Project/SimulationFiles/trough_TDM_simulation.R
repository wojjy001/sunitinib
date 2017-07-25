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
  source(paste0(sim.dir,"trough_TDM_simulation_options.R"))	# Simulation options
# Define output directory
  output.dir <- paste0(global.dir,"Output/TroughTDMSimulation/")
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
# Perform initial PK simulations (first week)
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
  initial.pk.data <- pk.mod %>% mrgsim(data = input.pk.data,
    carry.out = c("SIM","amt")) %>%
    as.data.frame

# ------------------------------------------------------------------------------
# Simulate standard TDM scenario
  standard.tdm <- function(initial.pk.data) {
  # Initiate a vector of sample times
    sample.times <- c(168)
  # Sampling frequency
    sample.freq <- 24	# hours
  # Define therapeutic window
    trough.target <- 0.05	# mg/L
    # auc.target <- 1.38	# mg*h/L
  # Make all predicted concentration (IPRE) after the first week (t = 168 hours)
  # become NA.  These will be filled with simulated concentrations during each
  # loop.
    pk.data <- initial.pk.data
    pk.data <- auc24.function(pk.data)
    pk.data$IPRE[pk.data$time > max(sample.times)] <- NA

  # Simulate concentrations until the end of TDM phase
    repeat {
    # Call information from previous dosing interval
      last.sample <- max(sample.times)	# Time of most recent sample
      prev.dose <- pk.data$amt[pk.data$time == last.sample-24]	# Previous dose
    # Calculate next dose based on information from previous dosing interval
      prev.IPRE <- pk.data$IPRE[pk.data$time == last.sample]	# Previous sampled conc
      if (prev.IPRE != trough.target) {
        next.dose <- trough.target/prev.IPRE*prev.dose
      } else {
        next.dose <- prev.dose
      }
      # prev.AUC24 <- pk.data$AUC24[pk.data$time == last.sample]	# Previous interval's AUC
      # if (prev.AUC24 != auc.target) {
      #   next.dose <- auc.target/prev.AUC24*prev.dose
      # } else {
      #   next.dose <- prev.dose
      # }
    # Input next.dose for simulation
      input.pk.data <- pk.data
      input.pk.data$amt[input.pk.data$time == last.sample] <- next.dose
      # Re-add evid and rate columns
        input.pk.data$cmt <- 1
        input.pk.data$evid <- 0
        input.pk.data$evid[input.pk.data$time %in% on.times] <- 1
    # Simulate concentrations
      pk.data <- pk.mod %>% mrgsim(data = input.pk.data,
        carry.out = c("SIM","amt")) %>% as.data.frame
    # Add the next sample time to the list of sample times
      next.sample <- last.sample + sample.freq
      sample.times <- sort(c(unique(c(sample.times,next.sample))))
    # Make all IPRE after the next sample == NA
      pk.data <- auc24.function(pk.data)
      pk.data$IPRE[pk.data$time > max(sample.times)] <- NA
    # If the last IPRE in the data frame is NA, then continue with the loop
      if (is.na(pk.data$IPRE[pk.data$time == 4*7*24]) == FALSE) break
    }	# repeat

    # Simulate the next series of concentrations when the individual's final dose
    # from the TDM phase is doubled
    # Call information from previous dosing interval
      last.sample <- max(sample.times)	# Time of most recent sample
      prev.dose <- pk.data$amt[pk.data$time == last.sample-24]	# Previous dose
      next.dose <- prev.dose*2
      after.TDM.times <- on.times[on.times >= 4*7*24]
    # Input next.dose for simulation
      input.pk.data <- pk.data
      input.pk.data$amt[input.pk.data$time %in% after.TDM.times] <- next.dose
      # Re-add evid and rate columns
        input.pk.data$cmt <- 1
        input.pk.data$evid <- 0
        input.pk.data$evid[input.pk.data$time %in% on.times] <- 1
    # Simulate concentrations
      pk.data <- pk.mod %>% mrgsim(data = input.pk.data,
        carry.out = c("SIM","amt")) %>% as.data.frame
      pk.data <- auc24.function(pk.data)
      pk.data
  }	# standard.tdm

# Simulate concentrations arising from standard TDM
  pk.data <- ddply(initial.pk.data, .(SIM,ID), standard.tdm, .progress = "text")

# ------------------------------------------------------------------------------
# Plot AUC24 versus trough concentration on the last day of 4 weeks treatment
  endTDM.pk.data <- pk.data[pk.data$time == 4*7*24,]
  next.pk.data <- pk.data[pk.data$time == 4*7*24+24,]

  plotobj1 <- NULL
  plotobj1 <- ggplot()
  plotobj1 <- plotobj1 + geom_point(aes(x = IPRE,y = AUC24),
    data = endTDM.pk.data,
    colour = "blue",shape = 1,size = 2,alpha = 0.3)
  plotobj1 <- plotobj1 + geom_point(aes(x = IPRE,y = AUC24),
    data = next.pk.data,
    colour = "red",shape = 1,size = 2,alpha = 0.3)
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
  endTDM.dose.pk.data <- pk.data[pk.data$time >= 4*7*24-24 &
    pk.data$time <= 4*7*24,]
  interval.times <- unique(last.dose.pk.data$time)
  unique.ID <- rep(1:ntotal,times = length(interval.times)) %>% sort
  last.dose.pk.data$uID <- unique.ID

  plotobj2 <- NULL
  plotobj2 <- ggplot(last.dose.pk.data)
  plotobj2 <- plotobj2 + geom_line(aes(x = time,y = IPRE,group = uID),
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

# ------------------------------------------------------------------------------
# Randomly sample some individuals and plot their concentration profiles
  rand.ID <- sample(ID.seq,25)
  rand.pk.data <- pk.data[pk.data$ID %in% rand.ID,]

# IPRE versus time
  plotobj3 <- NULL
  plotobj3 <- ggplot(rand.pk.data)
  plotobj3 <- plotobj3 + geom_line(aes(x = time/24/7,y = IPRE),
    colour = "red")
  plotobj3 <- plotobj3 + geom_point(aes(x = time/24/7,y = IPRE),
    data = rand.pk.data[rand.pk.data$time == 4*7*24,])
  plotobj3 <- plotobj3 + geom_vline(aes(xintercept = 3),
    linetype = "dashed")
  plotobj3 <- plotobj3 + geom_hline(aes(yintercept = 0.05),
    linetype = "dashed")
  plotobj3 <- plotobj3 + scale_y_continuous("Total Sunitinib Concentration (mg/L)",
    breaks = seq(from = 0,to = 0.2,by = 0.05),
    labels = seq(from = 0,to = 0.2,by = 0.05),lim = c(0,0.2))
  plotobj3 <- plotobj3 + scale_x_continuous("Time (weeks)",
    breaks = 0:8)
  plotobj3 <- plotobj3 + facet_wrap(~ID,ncol = 5)
  print(plotobj3)

# AUC24 versus time
  plotobj4 <- NULL
  plotobj4 <- ggplot(rand.pk.data)
  plotobj4 <- plotobj4 + geom_line(aes(x = time/24/7,y = AUC24),
    colour = "blue")
  plotobj4 <- plotobj4 + geom_point(aes(x = time/24/7,y = AUC24),
    data = rand.pk.data[rand.pk.data$time == 4*7*24,])
  plotobj4 <- plotobj4 + geom_vline(aes(xintercept = 3),
    linetype = "dashed")
  plotobj4 <- plotobj4 + scale_y_continuous("Total Sunitinib 24-hour AUC (mg*h/L)",
    breaks = c(0,1,2,3,4),labels = c(0,1,2,3,4))
  plotobj4 <- plotobj4 + scale_x_continuous("Time (weeks)",
    breaks = 0:8)
  plotobj4 <- plotobj4 + facet_wrap(~ID,ncol = 5)
  print(plotobj4)

# Dose versus time
  plotobj5 <- NULL
  plotobj5 <- ggplot(rand.pk.data[rand.pk.data$amt != 0,])
  plotobj5 <- plotobj5 + geom_step(aes(x = time/24/7,y = amt),
    colour = "darkgreen")
  plotobj5 <- plotobj5 + geom_point(aes(x = time/24/7,y = amt),
    data = rand.pk.data[rand.pk.data$time == 4*7*24-24,])
  plotobj5 <- plotobj5 + geom_vline(aes(xintercept = 3),
    linetype = "dashed")
  plotobj5 <- plotobj5 + scale_y_continuous("Sunitinib Daily Dose (mg)")
  plotobj5 <- plotobj5 + scale_x_continuous("Time (weeks)",
    breaks = 0:8)
  plotobj5 <- plotobj5 + facet_wrap(~ID,ncol = 5)
  print(plotobj5)

# ------------------------------------------------------------------------------
# Summaries
  summary.AUC <- summary.function(last.pk.data$AUC24)
  summary.IPRE <- summary.function(last.pk.data$IPRE)
  summary.amt <- summary.function(pk.data$amt[pk.data$time == 4*7*24-24])

  AUC.trough.data <- data.frame(
    ID = ID.seq,
    IPRE4 = pk.data$IPRE[pk.data$time == 4*7*24],
    AUC244 = pk.data$AUC24[pk.data$time == 4*7*24],
    IPRE4.next = pk.data$IPRE[pk.data$time == 4*7*24+24],
    AUC244.next = pk.data$AUC24[pk.data$time == 4*7*24+24],
    IPRE8 = pk.data$IPRE[pk.data$time == 8*7*24],
    AUC248 = pk.data$AUC24[pk.data$time == 8*7*24]
  )
  AUC.trough.data$IPRE.ratio <- AUC.trough.data$IPRE8/AUC.trough.data$IPRE4
  AUC.trough.data$AUC24.ratio <- AUC.trough.data$AUC248/AUC.trough.data$AUC244
  AUC.trough.data$IPRE.next.ratio <- AUC.trough.data$IPRE4.next/AUC.trough.data$IPRE4
  AUC.trough.data$AUC24.next.ratio <- AUC.trough.data$AUC244.next/AUC.trough.data$AUC244

  summary.AUC24.next.ratio <- summary.function(AUC.trough.data$AUC24.next.ratio)
  summary.IPRE.next.ratio <- summary.function(AUC.trough.data$IPRE.next.ratio)

  summary.AUC24.ratio <- summary.function(AUC.trough.data$AUC24.ratio)
  summary.IPRE.ratio <- summary.function(AUC.trough.data$IPRE.ratio)

# ------------------------------------------------------------------------------
# Save plots and summary statistics
# Plot
  ggsave(plot = plotobj1,
    filename = "week4_auc_vs_trough_TDM.png",
    dpi = 300,height = 5,width = 7)
  ggsave(plot = plotobj2,
    filename = "week4_IPRE_vs_time_TDM.png",
    dpi = 300,height = 5,width = 7)
  ggsave(plot = plotobj3,
    filename = "IPRE_vs_time_TDM_25rand.png",
    dpi = 300,height = 7,width = 10)
  ggsave(plot = plotobj4,
    filename = "AUC24_vs_time_TDM_25rand.png",
    dpi = 300,height = 7,width = 10)
  ggsave(plot = plotobj5,
    filename = "Dose_vs_time_TDM_25rand.png",
    dpi = 300,height = 7,width = 10)
