# Sunitinib Simulation
# This script simulates sunitinib concentrations following a specified dosing
# regimen for a population of unique individuals
# Simulation options are saved in a text file for reference
# ------------------------------------------------------------------------------
# Define and create output directory
# Name of simulation type - which will be used to name the output directory
  output.dir.name <- "bayes_tdm06"
  output.dir <- paste0(global.dir,"Output/",output.dir.name,"/")
  dir.create(file.path(output.dir),showWarnings = FALSE)
# Set the working directory
  setwd(output.dir)

# ------------------------------------------------------------------------------
# Population characteristics
# Number of individuals to simulate
  nid <- 1000
  nsim <- 1
  ntotal <- nid*nsim
  SIM.seq <- rep(1:nsim,times = nid) %>% sort
  ID.seq <- rep(1:nid,times = nsim)

# Covariates
# Weight (kg)
  WT.mean <- 82.3
  WT.sd <- 0.2
  WT.min <- 39
  WT.max <- 157
  WT <- rlnorm(nid,meanlog = log(WT.mean),sd = WT.sd)
  WT[WT < WT.min] <- WT.min
  WT[WT > WT.max] <- WT.max

# ------------------------------------------------------------------------------
# Dosing specifications
  dose <- 50	# mg
# Create sequence of dosing times based on "on" and "off" periods in cycles
  ncycles <- 33
  cycles <- 1:ncycles
  cycle.duration <- 6	# weeks
  on.duration <- 4	# Duration drug is actively administered (weeks)
  on.times <- llply(cycles,function(x) {
    on.times <- (cycles[x]-1)*cycle.duration*7*24+
    seq(from = 0,to = on.duration*7*24-24,by = 24)
  })	#llply
  on.times <- unlist(on.times)

# ------------------------------------------------------------------------------
# Time sequences
  nweeks <- ncycles*cycle.duration
# Pharmacokinetics
# Simulation times (hours)
  end.pk.time <- nweeks*24*7	# end simulation time
  pk.increment <- 24 # concentration collection times
  pk.times <- seq(from = 0,to = end.pk.time,by = pk.increment) # hours
# Create sequence of cycles for each individual
# Always begin with cycle 1
# Length of pk.times will always be odd
  cyc <- c(rep(cycles,(length(pk.times)-1)/ncycles),ncycles) %>% sort

# ------------------------------------------------------------------------------
# Therapeutic drug monitoring specifications
  TDM <- 1
  # 0 = Do not perform therapeutic drug monitoring
  # 1 = Perform therapeutic drug monitoring
  if (TDM == 0) {
    trough.day <- NA
    dose.day <- NA
    target <- NA
    trough.target <- NA
    trough.lower <- 0.05	# mg/L
    trough.upper <- 0.1	# mg/L
    AUC.target <- NA
    AUC.lower <- 1.4
    AUC.upper <- 2.6
    optim.day1 <- NA
    optim.day2 <- NA
    optim.day3 <- NA
  }
  if (TDM == 1) {
  # Sample trough on day Z of the first cycle
    trough.day <- 14
    # trough.sample <- NULL
  # Code if you want to sample and dose adjustment at the same time in every cycle
    # trough.sample <- llply(cycles,function(x) {
    #   for (i in 1:length(trough.day)) {
    #     trough.sample[i] <- (cycles[x]-1)*cycle.duration*7*24+trough.day[i]*24
    #   }
    #   trough.sample
    # })	#llply
    # trough.sample <- unlist(trough.sample)
    trough.sample <- trough.day*24
  # Dose adjustment on day X of the first cycle
    dose.day <- 16
    # dose.adjust <- NULL
    # dose.adjust <- llply(cycles,function(x) {
    #   for (i in 1:length(dose.day)) {
    #     dose.adjust[i] <- (cycles[x]-1)*cycle.duration*7*24+dose.day[i]*24
    #   }
    #   dose.adjust
    # })	#llply
    # dose.adjust <- unlist(dose.adjust)
    dose.adjust <- dose.day*24
    target <- "auc"	# "trough" or "AUC"
    trough.target <- NA	# mg/L
    trough.lower <- 0.05	# mg/L
    trough.upper <- 0.1	# mg/L
    AUC.target <- 2	# mg*h/L
    AUC.lower <- 1.4	# mg*h/L
    AUC.upper <- 2.6	# mg*h/L
    optim.day1 <- 28*24	# hours
    optim.day2 <- NA*24	# hours
    optim.day3 <- NA*24	# hours
  }
  target.days <- c(optim.day1,optim.day2,optim.day3)

# ------------------------------------------------------------------------------
# Save simulation options to file
  sim.options <- data.frame(output.dir.name,
    nid,nsim,ntotal,WT.mean,WT.sd,WT.min,WT.max,
    nweeks,dose,ncycles,cycle.duration,on.duration,
    TDM,trough.day = paste(trough.day[1],trough.day[2],sep = "_"),
    dose.day,target,
    trough.target,trough.lower,trough.upper,
    AUC.target,AUC.lower,AUC.upper,
    target.days = paste(optim.day1,optim.day2,optim.day3,sep = "_")
  )
# Read in previous simulation options
  if (file.exists(paste0(global.dir,"Output/sim_options.csv")) == TRUE) {
    prev.options <- read.csv(file = paste0(global.dir,"Output/sim_options.csv"))
    prev.options <- prev.options[prev.options$output.dir.name != output.dir.name,]
    sim.options <- rbind(prev.options,sim.options)
  }
  write.csv(sim.options,file = paste0(global.dir,"Output/sim_options.csv"),
    quote = FALSE,row.names = FALSE)
