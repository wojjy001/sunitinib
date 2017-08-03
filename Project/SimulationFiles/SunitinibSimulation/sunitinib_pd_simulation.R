# Sunitinib PD Simulation
# This script reads in previously simulated sunitinib concentrations and
# simulates biomarkers and adverse effects
# ------------------------------------------------------------------------------
# Remove objects in workspace
  rm(list = ls(all = TRUE))
# Define global directory for project
  global.dir <- "/Volumes/Prosecutor/sunitinib/Project/"
# Source required files
  source(paste0(global.dir,"functions.R"))	# Universal functions file
  source(paste0(global.dir,"ModelFiles/hansson_2013_sunitinib_models.R"))	# PD model
# Define simulation file directory
  options.dir <- paste0(global.dir,"SimulationFiles/")
  source(paste0(options.dir,"sunitinib_pd_simulation_options.R"))	# Simulation options
  setwd(sim.dir)

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
