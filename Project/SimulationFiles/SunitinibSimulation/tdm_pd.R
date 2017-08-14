# Sunitinib PD Simulation
# This script reads in previously simulated sunitinib concentrations and
# simulates biomarkers and adverse effects
# ------------------------------------------------------------------------------
# Remove objects in workspace
  rm(list = ls(all = TRUE))
################################################################################
# ONLY DIRECTORIES TO BE CHANGED DEPENDING ON USER
# Define global directory for project
  global.dir <- "/Volumes/Prosecutor/sunitinib/Project/"
# Create output directory (folder outside of git folder) if not already
  output.dir <- "/Volumes/Prosecutor/sunitinib_nogit/"
################################################################################
# Create output directory specifically for sensitivity analysis (if not already)
  tdm.output.dir <- paste0(output.dir,"TDMSimulations/")

# ------------------------------------------------------------------------------
# Source required files
  source(paste0(global.dir,"functions.R"))	# Universal functions file
  source(paste0(global.dir,"ModelFiles/hansson_2013_sunitinib_models.R"))	# PD model

# ------------------------------------------------------------------------------
# Define simulation file directory
# Name of simulation type - which will be used to name the output directory
  sim.dir.name <- "standard_dose03"
  sim.dir <- paste0(tdm.output.dir,sim.dir.name)
  setwd(sim.dir)
# Read in simulated data
  pk.data <- read.csv(file = paste0(sim.dir,"/",sim.dir.name,"_pk_data.csv"))
  nid <- length(unique(pk.data$ID))
  ID.seq <- unique(pk.data$ID)
  nsim <- length(unique(pk.data$SIM))
  SIM.seq <- unique(pk.data$SIM)
  ntotal <- nid*nsim

# ------------------------------------------------------------------------------
# Covariates
# Baseline tumour size
  first.OBASE <- exp(log(rnorm(nid,mean = 195,sd = 120)))
  repeat {
    new.OBASE <- first.OBASE[first.OBASE == "NaN"]
    new.OBASE <- exp(log(rnorm(length(new.OBASE),mean = 195,sd = 120)))
    if (length(new.OBASE[new.OBASE == "NaN"]) == 0) break
  }
  OBASE <- c(first.OBASE[first.OBASE != "NaN"],new.OBASE)
  # plot(hist(OBASE))
# Baseline hand-foot syndrome grade
  hfs.grade0.prob <- 0.95
  hfs.grade1.prob <- 0.02
  hfs.grade2.prob <- 0.02
  hfs.grade3.prob <- 0.01
  base.hfs.probs <- c(hfs.grade0.prob,hfs.grade1.prob,hfs.grade2.prob,
    hfs.grade3.prob)
  HFSBASE <- unlist(llply(seq_len(ntotal),
    function(x) sample(c(0,1,2,3),size = x/x,prob = base.hfs.probs)))
# Baseline fatigue grade
  fat.grade0.prob <- 0.85
  fat.grade1.prob <- 0.09
  fat.grade2.prob <- 0.05
  fat.grade3.prob <- 0.01
  base.fat.probs <- c(fat.grade0.prob,fat.grade1.prob,fat.grade2.prob,
    fat.grade3.prob)
  FATBASE <- unlist(llply(seq_len(ntotal),
    function(x) sample(c(0,1,2,3),size = x/x,prob = base.fat.probs)))

# ------------------------------------------------------------------------------
# Time sequences for simulation
  nweeks <- 200
# Pharmacodynamics
  end.pd.time <- nweeks*24*7	# end simulation time
  pd.increment <- 24*7	# measurement times
  pd.times <- seq(from = 0,to = end.pd.time,by = pd.increment)

# Time to event/overall survival
  end.surv.time <- nweeks	# end simulation time (200 weeks)
  surv.increment <- 2	# measurement times
  surv.times <- seq(from = 0,to = end.surv.time,by = surv.increment)

# ------------------------------------------------------------------------------
# Generate random effect parameters
# Between subject variability
  if (ntotal > 1) {
    pd.ETA.matrix <- mvrnorm(ntotal,
      mu = rep(0,times = dim(pd.OMEGA)[1]),pd.OMEGA) %>%
      as.data.frame
    names(pd.ETA.matrix) <- c("ETAVEGFR3BASE","ETAVEGFR3MRT","ETAVEGFR3I50",
      "ETASKITBASE","ETASKITMRT","ETASKITI50","ETASKITSLP","ETAKG","ETAKRSKIT",
      "ETAKRD","ETAOBASE","ETAANCBASE","ETAANCMTT","ETAANCEMAX","ETAANCE50",
      "ETABPBASE","ETABPSLP","ETABPMRT","ETAHFS0","ETAHFS1","ETAHFS2","ETAFAT0",
      "ETAFAT1","ETAFAT2","ETAFAT3")
  } else {
    pd.ETA.matrix <- data.frame(ETAVEGFR3BASE = 0,ETAVEGFR3MRT = 0,
      ETASKITBASE = 0,ETASKITMRT = 0,ETASKITSLP = 0,ETAVEGFR3I50 = 0,
      ETASKITI50 = 0,ETAKG = 0,ETAKRSKIT = 0,ETAKRD = 0,ETAOBASE = 0,
      ETAANCBASE = 0,ETAANCMTT = 0,ETAANCEMAX = 0,ETAANCE50 = 0,ETABPBASE = 0,
      ETABPSLP = 0,ETABPMRT = 0,ETAHFS0 = 0,ETAHFS1 = 0,ETAHFS2 = 0,ETAFAT0 = 0,
      ETAFAT1 = 0,ETAFAT2 = 0,ETAFAT3 = 0)
  }

# ------------------------------------------------------------------------------
# Input data frame for simulation
  pkpd.data <- pk.data[c("ID","time","amt","SIM","cyc","WT","IPREP","IPREM",
    "IPRE","AUC24")]
  pd.ID.data <- data.frame(SIM = SIM.seq,ID = ID.seq,OBASE,HFSBASE,FATBASE,
    pd.ETA.matrix)
  input.pd.data <- merge(pkpd.data,pd.ID.data,by = c("SIM","ID"),all = T)
  input.pd.data$cmt <- 0
  input.pd.data$evid <- 0
  input.pd.data$RSURV <- runif(dim(input.pd.data)[1],min = 0,max = 1)
  input.pd.data$RDROP <- runif(dim(input.pd.data)[1],min = 0,max = 1)

# ------------------------------------------------------------------------------
# Simulate biomarker and adverse effect profiles
  pd.data <- pd.mod %>% mrgsim(data = input.pd.data,
    carry.out = c("SIM","cyc","amt","IPREP","IPREM","IPRE","AUC24")) %>%
    as.data.frame
# Use hand-foot syndrome state probabilities to simulate HFS profile
  pd.data <- ddply(pd.data, .(ID), simulate.HFS.grade, .progress = "text")
# Use fatigue state probabilities to simulate FAT profile
  pd.data <- ddply(pd.data, .(ID), simulate.FAT.grade, .progress = "text")
# Determine overall survival for the population over the study period
  pd.data <- ddply(pd.data, .(ID), alive.function, .progress = "text")

# ------------------------------------------------------------------------------
# Clean up and save the PD simulated data
  output.pd.data <- pd.data[c("SIM","ID","time","cyc","amt","IPREP","IPREM",
    "IPRE","AUC24","TUMOUR","ANC","BP","WT","OBASE","HFSBASE","FATBASE",
    "IPRE_VEGFR3","IPRE_SKIT","HFS","FAT","status")]
  write.csv(output.pd.data,file = paste0(sim.dir.name,"_pd_data.csv"),
    quote = FALSE,row.names = FALSE)
