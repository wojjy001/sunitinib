# Create population of individuals to be simulated
# ------------------------------------------------------------------------------
# Basic population
  nid <- 1000	# Number of unique individuals
  nsim <- 1	# Number of times individuals with x characteristics will be
  # simulated
  ntotal <- nid*nsim	# Total number of individuals to be simulated
  SIM.seq <- rep(1:nsim,times = nid) %>% sort	# Sequence of simulation numbers
  ID.seq <- rep(1:nid,times = nsim)	# ID sequence

# ------------------------------------------------------------------------------
# Pharmacokinetic covariate parameters
# Weight (kg)
  WT.mean <- 82.3
  WT.sd <- 0.2
  WT.min <- 39
  WT.max <- 157
  WT <- rlnorm(nid,meanlog = log(WT.mean),sd = WT.sd)
  WT[WT < WT.min] <- WT.min
  WT[WT > WT.max] <- WT.max

# ------------------------------------------------------------------------------
# Pharmacodynamic covariate parameters
# Baseline tumour size (sum of longest diameters, mm)
  first.OBASE <- exp(log(rnorm(nid,mean = 195,sd = 120)))
  repeat {
    new.OBASE <- first.OBASE[first.OBASE == "NaN"]
    new.OBASE <- exp(log(rnorm(length(new.OBASE),mean = 195,sd = 120)))
    if (length(new.OBASE[new.OBASE == "NaN"]) == 0) break
  }
  OBASE <- c(first.OBASE[first.OBASE != "NaN"],new.OBASE)
# Baseline hand-foot syndrome grade
  hfs.grade0.prob <- 0.95
  hfs.grade1.prob <- 0.02
  hfs.grade2.prob <- 0.02
  hfs.grade3.prob <- 0.01
  base.hfs.probs <- c(hfs.grade0.prob,hfs.grade1.prob,hfs.grade2.prob,
    hfs.grade3.prob)
  HFSBASE <- unlist(llply(seq_len(nid),
    function(x) sample(c(0,1,2,3),size = x/x,prob = base.hfs.probs)))
# Baseline fatigue grade
  fat.grade0.prob <- 0.85
  fat.grade1.prob <- 0.09
  fat.grade2.prob <- 0.05
  fat.grade3.prob <- 0.01
  base.fat.probs <- c(fat.grade0.prob,fat.grade1.prob,fat.grade2.prob,
    fat.grade3.prob)
  FATBASE <- unlist(llply(seq_len(nid),
    function(x) sample(c(0,1,2,3),size = x/x,prob = base.fat.probs)))

# ------------------------------------------------------------------------------
# Pharmacokinetic random effect parameters
# Covariance terms, therefore need to use multivariate random number generator
  if (ntotal > 1) {
    pk.ETA.matrix <- mvrnorm(ntotal,
      mu = rep(0,times = dim(pk.OMEGA)[1]),pk.OMEGA) %>%
      as.data.frame
    names(pk.ETA.matrix) <- c("ETACLP","ETAVCP","ETACLM","ETAVCM")
  } else {
    pk.ETA.matrix <- data.frame(ETACLP = 0,ETAVCP = 0,ETACLM = 0,ETAVCM = 0)
  }

# ------------------------------------------------------------------------------
# Pharmacodynamic random effect parameters
# Covariance terms, therefore need to use multivariate random number generator
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
# Input PK data frame for simulation
  pk.ID.data <- data.frame(SIM = SIM.seq,ID = ID.seq,WT,OBASE,HFSBASE,FATBASE,
    pk.ETA.matrix)
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
