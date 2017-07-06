# Hansson EK, Amantea MA, Westwood P, Milligan PA, Houk BE, French J,
# Karlsson MO, Friberg LE. PKPD Modeling of VEGF, sVEGFR-2, sVEGFR-3, and sKIT
# as Predictors of Tumor Dynamics and Overall Survival Following Sunitinib
# Treatment in GIST. CPT: pharmacometrics & systems pharmacology. 2013;2(11):1-9
# ------------------------------------------------------------------------------
# Remove objects in workspace
  rm(list = ls(all = TRUE))
# Source simulation options file
  source("/Volumes/Prosecutor/sunitinib/sim_options.R")
# Source functions file
  source("/Volumes/Prosecutor/sunitinib/functions.R")
# Set the working directory
  dir <- "/Volumes/Prosecutor/sunitinib/PopPKPDModels/"
  setwd(dir)
# Source previous biomarker simulation
  source("hansson_2013_biomarker_simulation.R")
# Source model code
  source("hansson_2013_tumour_model.R")

# ------------------------------------------------------------------------------
# Generate random effects
  tum.ETA.matrix <- mvrnorm(nsim,
    mu = rep(0,times = dim(tum.OMEGA)[1]),tum.OMEGA) %>%
    as.data.frame
  tum.EPS.matrix <- mvrnorm(nsim*length(tum.times),
    mu = rep(0,times = dim(tum.SIGMA)[1]),tum.SIGMA) %>%
    as.data.frame

# ------------------------------------------------------------------------------
# Input data frame for simulation
  ID.data <- data.frame(
    ID = ID.seq,
    DOSE = 0,
    tum.ETA.matrix
  )
  names(ID.data)[c(3:dim(ID.data)[2])] <- c("ETAKG","ETAKSKIT","ETAKDRUG")
  input.data <- lapply(ID.data,rep.int,times = length(tum.times)) %>%
    as.data.frame
  input.data <- input.data[with(input.data,order(input.data$ID)),]
  input.data$time <- tum.times
  input.data$EPSBASE <- tum.EPS.matrix[,1]
  input.data$DOSE[input.data$time %in% on.times] <- DOSE
  input.data$CL <- PK.CL
  input.data$cmt <- 0

# Call on parameters from biomarker simulation for each individual
  bio.data.par <- ddply(bio.data, .(ID), headperID)
  input.data <- merge(input.data,
    bio.data.par[c("ID","BM03","IC503","MRT3","BM0S","IC50S","MRTS","DPSLOS")],
    by = c("ID"),all = TRUE)
  input.data <- input.data[with(input.data,order(input.data$ID,input.data$time)),]

# ------------------------------------------------------------------------------
# Simulate
  tum.data <- tum.mod %>% mrgsim(data = input.data) %>% as.data.frame

# ------------------------------------------------------------------------------
# Plot
  plotobj <- NULL
  plotobj <- ggplot(tum.data)
  plotobj <- plotobj + stat_summary(aes(x = time/7/24,y = TUMOUR),
    geom = "line",fun.y = median,colour = "red")
  plotobj <- plotobj + stat_summary(aes(x = time/7/24,y = TUMOUR),
    geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
    fill = "red",alpha = 0.3)
  plotobj <- plotobj + scale_y_continuous("Sum of Longest Diameters (mm)\n",lim = c(0,NA))
  plotobj <- plotobj + scale_x_continuous("\nTime (weeks)",
    breaks = seq(from = 0,to = 18,by = 2))
  plotobj
