# Sunitinib PD Simulation
# This script reads in previously simulated sunitinib concentrations and
# simulates biomarkers and adverse effects
# ------------------------------------------------------------------------------
# Define and create output directory
# Name of simulation type - which will be used to name the output directory
  sim.dir.name <- "bayes_tdm06"
  sim.dir <- paste0(global.dir,"Output/",sim.dir.name,"/")
# Read in simulated data
  pk.data <- read.csv(file = paste0(sim.dir,sim.dir.name,"_pk_data.csv"))
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
