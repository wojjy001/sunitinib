# Sunitinib Simulation
# This script simulates sunitinib concentrations following an optimisied dose
# for each individual that will achieve a target AUC
# ------------------------------------------------------------------------------
# Remove objects in workspace
  rm(list = ls(all = TRUE))
# Define global directory for project
  global.dir <- "/Volumes/Prosecutor/sunitinib/Project/"
# Source required files
  source(paste0(global.dir,"functions.R"))	# Universal functions file
  source(paste0(global.dir,"ModelFiles/yu_2015_sunitinib_model.R"))	# PK model
# Set output directory (folder outside of git folder)
  WT <- 100	# kg, all individuals will have the same weight
  output.dir <- paste0("/Volumes/Prosecutor/sunitinib_nogit/SensitivityAnalysis/",
    WT,"kg/")
  dir.create(output.dir)

# ------------------------------------------------------------------------------
# Create population
  nid <- 500	# Number of individuals
  nsim <- 1	# Number of times individuals with x characteristics will be
  # simulated
  ntotal <- nid*nsim	# Total number of individuals to be simulated
  SIM.seq <- rep(1:nsim,times = nid) %>% sort
  ID.seq <- rep(1:nid,times = nsim)
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

# Dosing specifications
  dose <- 50	# mg, an initial dose that will be overwritten during dose optim
# Create sequence of dosing times based on "on" and "off" periods in cycles
  ncycles <- 30
  cycles <- 1:ncycles
  cycle.duration <- 6	# weeks
  on.duration <- 4	# Duration drug is actively administered (weeks)
  on.times <- llply(cycles,function(x) {
    on.times <- (cycles[x]-1)*cycle.duration*7*24+
    seq(from = 0,to = on.duration*7*24-24,by = 24)
  })	#llply
  on.times <- unlist(on.times)

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

# Input data frame for simulation
  pk.ID.data <- data.frame(SIM = SIM.seq,ID = ID.seq,WT,pk.ETA.matrix)
  input.pk.data <- lapply(pk.ID.data,rep.int,times = length(pk.times)) %>%
    as.data.frame
  input.pk.data <- input.pk.data[with(input.pk.data,
    order(input.pk.data$SIM,input.pk.data$ID)),]
  input.pk.data$time <- pk.times
  input.pk.data$cyc <- cyc
  input.pk.data$amt <- dose
  input.pk.data$cmt <- 1
  input.pk.data$evid <- 0
  input.pk.data$evid[input.pk.data$time %in% on.times] <- 1

# ------------------------------------------------------------------------------
# Define target auc concentrations
  target.auc <- c(1,2,3,4)
# Run simulations process for each value in "target.auc"
  for (i in 1:length(target.auc)) {
  # Create a folder for output to be saved
  # Each "target.auc" value will have its own directory
    target.dir.name <- paste0("target_auc_",target.auc[i],"mgL_",
      WT,"kg")
    target.dir <- paste0(output.dir,target.dir.name,"/")
    dir.create(target.dir)
    setwd(target.dir)
    print(paste0("Simulating ",target.dir.name))
  # Optimise the dose for each individual that will achieve the target auc
  # concentration by Day 28 of the first cycle of treatment
    pop.optimise.dose <- function(input.df) {
    # Subset the input data frame for only the first cycle
    # Optimisation takes too long if all cycles are included
      input.optim.data <- input.df[input.df$cyc == 1,]
    # Initial parameter estimates
      initial.dose <- dose	# standard 50 mg dose
      initial.err <- 0.01
      par <- c(initial.dose,initial.err)
    # Optimise dose function
      optimise.dose <- function(par) {
      # Assign estimable parameters to objects
        input.optim.data$amt[input.optim.data$evid == 1] <- par[1]
        err <- par[2]
      # Simulate concentration-time profile with each iteration of dose
        optim.data <- pk.mod %>% mrgsim(data = input.optim.data) %>%
          as.data.frame
        optim.data <- auc24.function(optim.data)
      # Pull out the predicted AUC at the target time
        yhat <- optim.data$AUC24[optim.data$time == 28*24]
      # Find the value of dose that maximises the likelihood of the predicted
      # concentration being the target AUC
        loglik <- dnorm(target.auc[i],yhat,yhat*err,log = TRUE)
      # Define the objective function value that will be optimised by "optim"
        objective <- -1*sum(loglik)
      }
    # Run the optim function to obtain a value for dose
      optimised.dose <- optim(par,
        optimise.dose,
        hessian = FALSE,
        method = "L-BFGS-B",
        lower = c(0.0001,0.0001),upper = c(Inf,Inf)
      )
    # Create output object for each individual with their individual dose
      dose.data <- data.frame(dose = optimised.dose$par[1])
    }
  # Run "pop.optimise.dose" for each individual in the population
    dose.data <- ddply(input.pk.data, .(SIM,ID), pop.optimise.dose,
      .progress = "text")
  # Merge dose.data into input.pk.data
    add.doses <- function(dose.data) {
      ID <- dose.data$ID[1]
      SIM <- dose.data$SIM[1]
      input.sim.data <- input.pk.data[input.pk.data$SIM == SIM &
        input.pk.data$ID == ID,]
      input.sim.data$amt[input.sim.data$evid == 1] <- dose.data$dose[1]
      input.sim.data
    }
    input.sim.data <- ddply(dose.data, .(SIM,ID), add.doses)
  # Simulate concentrations for each individual given their optimised doses
    sim.data <- pk.mod %>% mrgsim(data = input.sim.data,
      carry.out = c("SIM","cyc","amt")) %>% as.data.frame
    sim.data <- ddply(sim.data, .(SIM,ID), auc24.function)
  # Plot concentrations and 24-hour AUC for the first 2 cycles for the population
    early.data <- sim.data[sim.data$cyc <= 2,]
    summary.IPRE <- ddply(early.data, .(time),
      function(early.data) graded.summary(early.data$IPRE))
    summary.AUC24 <- ddply(early.data, .(time),
      function(early.data) graded.summary(early.data$AUC24))
    # Plot IPRE
      plotobj1 <- NULL
      plotobj1 <- ggplot(summary.IPRE)
      plotobj1 <- plotobj1 + geom_ribbon(aes(x = time/24,ymin = CI90lo,
        ymax = CI90hi),fill = "red",alpha = 0.1)
      plotobj1 <- plotobj1 + geom_ribbon(aes(x = time/24,ymin = CI80lo,
        ymax = CI80hi),fill = "red",alpha = 0.1)
      plotobj1 <- plotobj1 + geom_ribbon(aes(x = time/24,ymin = CI60lo,
        ymax = CI60hi),fill = "red",alpha = 0.1)
      plotobj1 <- plotobj1 + geom_ribbon(aes(x = time/24,ymin = CI40lo,
        ymax = CI40hi),fill = "red",alpha = 0.1)
      plotobj1 <- plotobj1 + geom_ribbon(aes(x = time/24,ymin = CI20lo,
        ymax = CI20hi),fill = "red",alpha = 0.1)
      plotobj1 <- plotobj1 + geom_line(aes(x = time/24,y = med),
        colour = "red")
      plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 0.05),
        linetype = "dashed")
      plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 0.1),
        linetype = "dashed")
      plotobj1 <- plotobj1 + scale_y_continuous("Total Sunitinib Concentration (mg/L)",lim = c(0,NA))
      plotobj1 <- plotobj1 + scale_x_continuous("Time (days)",
        breaks = seq(from = 0,to = 84,by = 14),
        labels = seq(from = 0,to = 84,by = 14))
      print(plotobj1)
      ggsave(plot = plotobj1,filename = paste0(target.dir.name,"_IPREvstime.png"),
        height = 15,width = 20,units = "cm",dpi = 300)
    # Plot AUC24
      plotobj2 <- NULL
      plotobj2 <- ggplot(summary.AUC24)
      plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24,ymin = CI90lo,
        ymax = CI90hi),fill = "blue",alpha = 0.1)
      plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24,ymin = CI80lo,
        ymax = CI80hi),fill = "blue",alpha = 0.1)
      plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24,ymin = CI60lo,
        ymax = CI60hi),fill = "blue",alpha = 0.1)
      plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24,ymin = CI40lo,
        ymax = CI40hi),fill = "blue",alpha = 0.1)
      plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24,ymin = CI20lo,
        ymax = CI20hi),fill = "blue",alpha = 0.1)
      plotobj2 <- plotobj2 + geom_line(aes(x = time/24,y = med),
        colour = "blue")
      plotobj2 <- plotobj2 + geom_hline(aes(yintercept = target.auc[i]),
        linetype = "dashed")
      plotobj2 <- plotobj2 + scale_y_continuous("Total Sunitinib 24-hour AUC (mg*h/L)",lim = c(0,target.auc[i]+0.5))
      plotobj2 <- plotobj2 + scale_x_continuous("Time (days)",
        breaks = seq(from = 0,to = 84,by = 14),
        labels = seq(from = 0,to = 84,by = 14))
      print(plotobj2)
      ggsave(plot = plotobj2,filename = paste0(target.dir.name,"_AUC24vstime.png"),
        height = 15,width = 20,units = "cm",dpi = 300)
  # Save output data frame
    output.sim.data <- sim.data
    write.csv(output.sim.data,file = paste0(target.dir.name,"_pk_data.csv"),
      quote = FALSE,row.names = FALSE)
    target.sim.data <- sim.data[sim.data$time == 28*24,]
    write.csv(target.sim.data,file = paste0(target.dir.name,"_day28_data.csv"),
      quote = FALSE,row.names = FALSE)
    write.csv(dose.data,file = paste0(target.dir.name,"_dose_data.csv"),
      quote = FALSE,row.names = FALSE)
  }
