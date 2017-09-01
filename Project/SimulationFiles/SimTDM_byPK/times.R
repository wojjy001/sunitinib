# Define time sequences for simulation
# Required to create population
# ------------------------------------------------------------------------------
# Create sequence of dosing times based on "on" and "off" periods in cycles
  ncycles <- 2	# Number of treatment cycles
  cycles <- 1:ncycles	# Sequence of cycle numbers
  cycle.duration <- 6	# weeks
  on.duration <- 4	# Duration drug is actively administered (weeks)
  on.times <- llply(cycles,function(x) {
    on.times <- (cycles[x]-1)*cycle.duration*7*24+
    seq(from = 0,to = on.duration*7*24-24,by = 24)
  })	#llply
  on.times <- unlist(on.times)

# ------------------------------------------------------------------------------
# Pharmacokinetics
# Time sequences
  nweeks <- ncycles*cycle.duration
# Simulation times (hours)
  end.pk.time <- cycle.duration*ncycles*24*7	# end PK simulation time (hours)
  pk.increment <- 24 # concentration collection times (hours)
  pk.times <- seq(from = 0,to = end.pk.time,by = pk.increment) # hours
# Create sequence of cycles for each individual
# Always begin with cycle 1
# Length of pk.times will always be odd
  cyc <- c(rep(cycles,(length(pk.times)-1)/ncycles),ncycles) %>% sort

# ------------------------------------------------------------------------------
# Pharmacodynamics
  pd.ncycles <- 35	# Number of treatment cycles to be simulated
  pd.cycles <- 1:pd.ncycles	# Sequence of cycle numbers
  pd.nweeks <- pd.ncycles*cycle.duration	# Number of weeks
  pd.increment <- 24*7	# measurement times in hours
  pd.times <- seq(from = 0,to = pd.nweeks*pd.increment,by = pd.increment)
# Create sequence of cycles for each individual
# Always begin with cycle 1
  pd.cyc <- c(rep(pd.cycles,(length(pd.times)-1)/pd.ncycles),pd.ncycles) %>%
    sort
