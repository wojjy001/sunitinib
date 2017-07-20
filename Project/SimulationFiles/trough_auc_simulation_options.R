# Trough-AUC Simulation
# Simulation exercise to explore the range of AUCs that are generated when
# various individual trough concentrations are 50 ng/mL (i.e., target trough
# concentration for sunitinib)
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
  WT <- rlnorm(nid,meanlog = log(82.3),sd = 0.2)
  WT[WT < 39] <- 39
  WT[WT > 157] <- 157

# ------------------------------------------------------------------------------
# Time sequences
  nweeks <- 4
# Pharmacokinetics
# Simulation times (hours)
  end.pk.time <- nweeks*24*7	# end simulation time
  pk.increment <- 1 # concentration collection times
  pk.times <- seq(from = 0,to = end.pk.time,by = pk.increment) # hours

# ------------------------------------------------------------------------------
# Dosing specifications
  dose <- 50	# mg
  ncycles <- 1
  cycles <- 1:ncycles	# number of cycles
  cycle.duration <- nweeks	# weeks
  on.duration <- nweeks	# Duration drug is actively administered (weeks)
  off.duration <- cycle.duration-on.duration
  on.times <- seq(from = 0,to = on.duration*7*24-24,by = 24)	# 4 weeks in hours
  if (ncycles > 1) {
    add.times <- llply(seq_len(max(cycles)),function(x) {
      add.times <- cycles[x]*cycle.duration*7*24+
      seq(from = 0,to = on.duration*7*24-24,by = 24)
    })
    add.times <- unlist(add.times)
    on.times <- c(on.times,add.times)
  }
