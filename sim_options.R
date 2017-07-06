# Simulation options for sunitinib project
# ------------------------------------------------------------------------------
# Population characteristics
# Number of individuals to simulate
  nsim <- 500
  ID.seq <- 1:nsim

  PK.CL <- 34

# ------------------------------------------------------------------------------
# Time sequences
# Biomarkers (VEGF,VEGFR2,VEGFR3,SKIT)
  bio.times <- seq(from = 0,to = 200*7*24,by = 24)	# hours
# Tumour
  tum.times <- seq(from = 0,to = 200*7*24,by = 24)	# hours
# Survival
  sur.times <- seq(from = 0,to = 200,by = 1)	# weeks

# ------------------------------------------------------------------------------
# Dosing specifications
  DOSE <- 50	# mg
  cycles <- 1:30
  on.times <- seq(from = 24,to = 3*7*24,by = 24)	# 4 weeks in hours
  add.times <- llply(seq_len(30),function(x) {
    add.times <- cycles[x]*6*7*24+seq(from = 24,to = 3*7*24,by = 24)
  })
  add.times <- unlist(add.times)
  on.times <- c(on.times,add.times)
  on.times.weeks <- unique(round(on.times/7/24,digits = 0))
