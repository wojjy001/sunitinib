# Sunitinib Simulation
# This script simulates sunitinib concentrations following a specified dosing
# regimen for a population of unique individuals
# Dose adjustments can be made by standard therapeutic drug monitoring or
# based on Bayes predictions
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
  dir.create(output.dir)
# Create output directory specifically for sensitivity analysis (if not already)
  tdm.output.dir <- paste0(output.dir,"TDMSimulations/")
  dir.create(tdm.output.dir)

# ------------------------------------------------------------------------------
# Source required files
  source(paste0(global.dir,"functions.R"))	# Universal functions file
  source(paste0(global.dir,"ModelFiles/yu_2015_sunitinib_model.R"))	# PK model

# ------------------------------------------------------------------------------
# Population characteristics
# Number of individuals to simulate
  nid <- 10
  nsim <- 1
  ntotal <- nid*nsim
  SIM.seq <- rep(1:nsim,times = nid) %>% sort
  ID.seq <- rep(1:nid,times = nsim)

# Covariates
# Weight (kg)
  WT.mean <- 82.3
  WT.sd <- 0.2
  WT.min <- 39
  WT.max <- 157
  WT <- rlnorm(nid,meanlog = log(WT.mean),sd = WT.sd)
  WT[WT < WT.min] <- WT.min
  WT[WT > WT.max] <- WT.max

# ------------------------------------------------------------------------------
# Dosing specifications
  dose <- 50	# mg - Initial dose administered prior to therapeutic drug monitoring
# Create sequence of dosing times based on "on" and "off" periods in cycles
  ncycles <- 33
  cycles <- 1:ncycles
  cycle.duration <- 6	# weeks
  on.duration <- 4	# Duration drug is actively administered (weeks)
  on.times <- llply(cycles,function(x) {
    on.times <- (cycles[x]-1)*cycle.duration*7*24+
    seq(from = 0,to = on.duration*7*24-24,by = 24)
  })	#llply
  on.times <- unlist(on.times)

# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# Therapeutic drug monitoring specifications
################################################################################
  TDM <- 0
  # 0 = Do not perform therapeutic drug monitoring
  # 1 = Perform therapeutic drug monitoring
  output.dir.name <- "standard_dose03"
  # standard_doseXX = no TDM, TDM has to equal zero
  # standard_tdmXX = standard TDM, TDM == 1
  # pro_tdmXX = TDM with proportional dose adjustment, TDM == 1
  # bayes_tdmXX = TDM based on Bayes predictions, TDM == 1
  output.dir <- paste0(tdm.output.dir,output.dir.name,"/")
  dir.create(output.dir)
  setwd(output.dir)
################################################################################
  if (TDM == 0) {	# Values within these brackets do not need to be changed
    trough.day <- NA
    dose.day <- NA
    target <- NA
    trough.target <- NA
    trough.lower <- 0.05	# mg/L
    trough.upper <- 0.1	# mg/L
    AUC.target <- NA
    AUC.lower <- 1.4
    AUC.upper <- 2.6
    optim.day1 <- NA
    optim.day2 <- NA
    optim.day3 <- NA
  }
  if (TDM == 1) {
################################################################################
  # Sample trough on day X of the first cycle
    trough.day <- 14
    trough.sample <- trough.day*24	# Convert to hours
  # Dose adjustment on day Y of the first cycle
    dose.day <- 16
    dose.adjust <- dose.day*24	# Convert to hours
    target <- "auc"	# "trough" or "AUC"
    trough.target <- NA	# mg/L
    trough.lower <- 0.05	# mg/L
    trough.upper <- 0.1	# mg/L
    AUC.target <- 2	# mg*h/L
    AUC.lower <- 1.4	# mg*h/L
    AUC.upper <- 2.6	# mg*h/L
    optim.day1 <- 28*24	# hours
    optim.day2 <- NA*24	# hours
    optim.day3 <- NA*24	# hours
################################################################################
  }
  target.days <- c(optim.day1,optim.day2,optim.day3)

# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
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
# Create optim.data object to store the Bayes predicted profile and predicted profile following dose optimisation
  optim.data <- data.frame(ID = 0,time = 0,SIM = 0,cyc = 0,amt = 0,DEPOT = 0,
    CENTP = 0,CENTM = 0,PERIM = 0,AUC = 0,WT = 0,OBASE = 0,HFSBASE = 0,
    FATBASE = 0,IPREP = 0,IPREM = 0,IPRE = 0,DVP = 0,DVM = 0,DV = 0,CLP = 0,
    VCP = 0,CLM = 0,VCM = 0,QI = 0,VPM = 0,AUC24 = 0)

# Simulate concentrations
  simulate.pk <- function(input.pk.data) {
  # No therapeutic drug monitoring
    if (TDM == 0) {
    # Simulate as normal
      pk.data <- pk.mod %>% mrgsim(data = input.pk.data,
        carry.out = c("SIM","cyc","amt")) %>% as.data.frame
    }	# if TDM == 0
  # Therapeutic drug monitoring
    if (TDM == 1) {
    # Times where the new dose will be administered
      new.on.times <- on.times[on.times >= dose.adjust]
    # Simulate standard dosing
      initial.pk.data <- pk.mod %>%
        mrgsim(data = input.pk.data,carry.out = c("SIM","cyc","amt")) %>%
        as.data.frame
    # At day 14 sample a trough concentration
      sample.DV <- initial.pk.data$DV[initial.pk.data$time %in% trough.sample]
      prev.dose <- initial.pk.data$amt[initial.pk.data$time ==
        min(trough.sample)]
    # On day 16, make a single dose adjustment based on trough concentration
      if (contains("standard_tdm",vars = output.dir.name) == 1) {
        if (sample.DV < trough.lower) {
          new.dose <- prev.dose+12.5	# mg
        } else if (sample.DV > trough.upper) {
          new.dose <- prev.dose-12.5 # mg
        } else {
          new.dose <- prev.dose	# mg
        }
      }	#if contains "standard_tdm"
      if (contains("pro_tdm",vars = output.dir.name) == 1) {
        if (sample.DV < trough.lower | sample.DV > trough.upper) {
          new.dose <- trough.target/sample.DV*prev.dose
        } else {
          new.dose <- prev.dose
        }
      }	#if contains "pro_tdm"
      if (contains("bayes_tdm",vars = output.dir.name) == 1) {
      # Obtain empirical Bayes estimates of parameters for the individual
        input.bayes.data <- input.pk.data[input.pk.data$cyc == 1,]
      # Initial estimates for Bayes parameters
      # Use population typical the first time bayes parameters are estimated
        initial.bayes.par <- exp(c(0,0,0,0))
        run.once <- FALSE
        repeat {
        # Use previous parameter results as initial estimates mulitpled by a
        # random number so that they are not exactly the same
          if (run.once == TRUE) {
            initial.bayes.par <- initial.bayes.par*exp(runif(4,min = -0.01,
              max = 0.01))
          }
          par.bayes <- initial.bayes.par
        # Bayesian estimation using population PK model
        # Maximum likelihood estimation for Bayes parameters
          bayes.estimate <- function(par) {
          # Describe parameters to be optimised in the exponential domain
            ETACLPfit <- log(par[1])
            ETAVCPfit <- log(par[2])
            ETACLMfit <- log(par[3])
            ETAVCMfit <- log(par[4])
          # Add fitted parameters to the input data frame
            input.bayes.data$ETACLP <- ETACLPfit
            input.bayes.data$ETAVCP <- ETAVCPfit
            input.bayes.data$ETACLM <- ETACLMfit
            input.bayes.data$ETAVCM <- ETAVCMfit
          # Simulate concentration-time profiles with fitted parameters
            bayes.data <- pk.mod %>% mrgsim(data = input.bayes.data) %>%
              as.data.frame
            bayes.data$IPRE[is.finite(bayes.data$IPRE) == FALSE |
              bayes.data$IPRE < .Machine$double.eps] <- .Machine$double.eps
          # Pull out the predicted trough concentrations at sample time
            yhat <- bayes.data$IPRE[bayes.data$time %in% trough.sample]
          # Posterior log-likelihood
          # Error model: Y = IPRE*(1+ERRPRO), Y = IPRE + IPRE*ERRPRO
            loglikpost.sd <- ERRPROP1
            loglikpost <- dnorm(sample.DV,mean = yhat,sd = yhat*loglikpost.sd,
              log = TRUE)
          # Prior log-likelihood
            ETA <- c(ETACLPfit,ETAVCPfit,ETACLMfit,ETAVCMfit)
            ETABSV <- c(BSVCLP,BSVVCP,BSVCLM,BSVVCM)
            loglikprior <- dnorm(ETA,mean = 0,sd = ETABSV,log = TRUE)
          # Objective function value and minimise the value
            objective <- -1*sum(loglikpost,loglikprior)
          }
        # Run bayes.estimate function through optim
          bayes.result <- optim(par = par.bayes,
            bayes.estimate,
            hessian = FALSE,
            method = "L-BFGS-B",
            lower = c(0.001,0.001,0.001,0.001),upper = c(Inf,Inf,Inf,Inf)
          )
          run.once <- TRUE
          if (bayes.result$message ==
            "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH") break
        }	#repeat
      # Calculate an optimised dose for the individual given Bayes estimated
      # paremeters
      # Add empirical Bayes estimates to date frame
      # Estimates for parameters are in the exponential domain
        input.optim.data <- input.bayes.data
        input.optim.data$ETACLP <- log(bayes.result$par[1])
        input.optim.data$ETAVCP <- log(bayes.result$par[2])
        input.optim.data$ETACLM <- log(bayes.result$par[3])
        input.optim.data$ETAVCM <- log(bayes.result$par[4])
      # Dose optimisation
      # Initial estimates for dose and error
        initial.dose <- dose	# 50 mg - standard dose
        initial.error <- 0.01
        par.optim <- c(initial.dose,initial.error)
      # Find the doses that maximum the likelihood of trough
      # concentrations or AUC being the target
        optimise.dose <- function(par) {
        # Add fitted parameters to the input data frame
          input.optim.data$amt[input.optim.data$time %in%
            new.on.times] <- par[1]
          err <- par[2]
        # Simulate concentration-time profiles with fitted doses
          optim.data <- pk.mod %>% mrgsim(data = input.optim.data) %>%
            as.data.frame
          optim.data$IPRE[is.finite(optim.data$IPRE) == FALSE |
            optim.data$IPRE < .Machine$double.eps] <- .Machine$double.eps
          optim.data <- auc24.function(optim.data)
        # Minimise the error between target trough and predicted trough
        # concentrations
        # OR Minimise the error between target AUC and predicted AUC
          if (contains("trough",vars = target) == 1) {
            optim.target <- trough.target
          # Pull out the predicted trough concentrations with the fitted
          # doses for the interval
            yhat <- optim.data$IPRE[optim.data$time %in% target.days]
          }
          if (contains("AUC",vars = target) == 1) {
            optim.target <- AUC.target
          # Pull out the predicted AUC with the fitted doses for the interval
            yhat <- optim.data$AUC24[optim.data$time %in% target.days]
          }
          res <- dnorm(optim.target,yhat,yhat*err,log = TRUE)
        # Objective function value and minimise the value
          objective <- -1*sum(res)
        }
      # Optimise dose depending on trough or AUC target
        optimised.dose <- optim(par = par.optim,
          optimise.dose,
          hessian = FALSE,
          method = "L-BFGS-B",
          lower = c(0.0001,0.0001),upper = c(Inf,Inf)
          # control = list(parscale = par.optim,factr = 1e7)
        )
        new.dose <- optimised.dose$par[1]
      # Continue dosing given dose adjustment using Bayes parameters
        input.optim.data$amt[input.optim.data$time %in% new.on.times] <- new.dose
        ind.optim.data <- pk.mod %>% mrgsim(data = input.optim.data,
          carry.out = c("SIM","cyc","amt")) %>% as.data.frame
        ind.optim.data <- auc24.function(ind.optim.data)
        optim.data <<- rbind(optim.data,ind.optim.data)
      }	#if contains "bayes_tdm"
    # Continue dosing given dose adjustment using true parameters
      new.dose <- round(new.dose/12.5,digits = 0)*12.5
      new.dose[new.dose < 12.5] <- 12.5
      new.dose[new.dose > 87.5] <- 87.5
      input.pk.data$amt[input.pk.data$time %in% new.on.times] <- new.dose
      pk.data <- pk.mod %>% mrgsim(data = input.pk.data,
        carry.out = c("SIM","cyc","amt")) %>% as.data.frame
    }	#if TDM == 1
  # Calculate daily AUC (24 hours)
    pk.data <- auc24.function(pk.data)
    pk.data
  }	#simulate.pk

# Simulate concentrations arising from standard TDM
  pk.data <- ddply(input.pk.data, .(SIM,ID), simulate.pk, .progress = "text")

# ------------------------------------------------------------------------------
# Save output data frame
  output.pk.data <- pk.data
  write.csv(output.pk.data,file = paste0(output.dir.name,"_pk_data.csv"),
    quote = FALSE,row.names = FALSE)

# Save optim.data
  if (contains("bayes_tdm",vars = output.dir.name) == 1) {
    output.optim.data <- optim.data
    write.csv(output.optim.data,file = paste0(output.dir.name,"_optim_data.csv"),
      quote = FALSE,row.names = FALSE)
  }
