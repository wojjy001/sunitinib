# Hansson EK, Amantea MA, Westwood P, Milligan PA, Houk BE, French J,
# Karlsson MO, Friberg LE. PKPD Modeling of VEGF, sVEGFR-2, sVEGFR-3, and sKIT
# as Predictors of Tumor Dynamics and Overall Survival Following Sunitinib
# Treatment in GIST. CPT: pharmacometrics & systems pharmacology. 2013;2(11):1-9
# ------------------------------------------------------------------------------
# Remove objects in workspace
  rm(list = ls(all = TRUE))
# Source functions file
  source("/Volumes/Prosecutor/sunitinib/functions.R")
# Source simulation options file
  source("/Volumes/Prosecutor/sunitinib/sim_options.R")
# Set the working directory
  dir <- "/Volumes/Prosecutor/sunitinib/PopPKPDModels/"
  setwd(dir)
# Source previous tumour simulation
  source("hansson_2013_tumour_simulation.R")
# Source model code
  source("hansson_2013_bsurvival_model.R")

# ------------------------------------------------------------------------------
# Input data frame for simulation
  ID.data <- data.frame(
    ID = ID.seq,
    DOSE = 0,
    REVENT = runif(length(ID.seq),min = 0,max = 1),
    RDROP = runif(length(ID.seq),min = 0,max = 1)
  )
  input.data <- lapply(ID.data,rep.int,times = length(sur.times)) %>%
    as.data.frame
  input.data <- input.data[with(input.data,order(input.data$ID)),]
  input.data$time <- sur.times
  input.data$DOSE[input.data$time %in% on.times.weeks] <- DOSE
  input.data$CL <- PK.CL
  input.data$cmt <- 0

# Call on parameters from biomarker simulation for each individual
  bio.data.par <- ddply(bio.data, .(ID), headperID)
  input.data <- merge(input.data,
    bio.data.par[c("ID","BM03","IC503","MRT3")],
    by = c("ID"),all = TRUE)
  input.data <- input.data[with(input.data,order(input.data$ID,input.data$time)),]

# ------------------------------------------------------------------------------
# Simulate
  bsur.data <- bsur.mod %>%
    mrgsim(data = input.data,carry.out = c("REVENT","RDROP")) %>%
    as.data.frame

# ------------------------------------------------------------------------------
# Calculate proportion surviving (accounting for drop out)
  # Assign whether individual is dead or alive based on survival and drop out
  # probabilities
  bsur.data <- alive.function(bsur.data)
  # Remove individuals who died or dropped out within the first 4 weeks
  early.death <- bsur.data$ID[bsur.data$time == 5 & bsur.data$ALIVE == 0]
  # input.alive.data <- bsur.data[!bsur.data$ID %in% early.death,]
  input.alive.data <- bsur.data
  alive.data <- ddply(input.alive.data, .(time), pro.alive.function)

# ------------------------------------------------------------------------------
# Plot
  plotobj <- NULL
  plotobj <- ggplot(alive.data)
  plotobj <- plotobj + geom_step(aes(x = time,y = pro.alive*100),
    colour = "red")
  plotobj <- plotobj + scale_y_continuous("Survival (%)\n",lim = c(0,NA))
  plotobj <- plotobj + scale_x_continuous("\nTime (weeks)")
  plotobj
