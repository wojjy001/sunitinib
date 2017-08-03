# Test efficiency of Bayesian estimation with trough concentrations
# ------------------------------------------------------------------------------
# Remove objects in workspace
  rm(list = ls(all = TRUE))
# Define global directory for project
  global.dir <- "/Volumes/Prosecutor/sunitinib/Project/"
# Source required files
  source(paste0(global.dir,"functions.R"))	# Universal functions file
  source(paste0(global.dir,"ModelFiles/yu_2015_sunitinib_model.R"))	# PK model
# Create and set working directory
  sample.time <- 1
  conc <- "combined"
  output.dir.name <- paste0("bayestest_sample",sample.time,"_conc",conc)
  output.dir <- paste0(global.dir,"Tests/",output.dir.name,"/")
  dir.create(file.path(output.dir),showWarnings = FALSE)
# Set the working directory
  setwd(output.dir)

# ------------------------------------------------------------------------------
# Create test population
  nid <- 500
  ID.seq <- 1:nid
# Assign covariate values
# Weight (kg)
  WT.mean <- 82.3
  WT.sd <- 0.2
  WT.min <- 39
  WT.max <- 157
  WT <- rlnorm(nid,meanlog = log(WT.mean),sd = WT.sd)
  WT[WT < WT.min] <- WT.min
  WT[WT > WT.max] <- WT.max
# Generate random effect parameters
  pk.ETA.matrix <- mvrnorm(nid,
    mu = rep(0,times = dim(pk.OMEGA)[1]),pk.OMEGA) %>%
    as.data.frame
  names(pk.ETA.matrix) <- c("ETACLP","ETAVCP","ETACLM","ETAVCM")

# ------------------------------------------------------------------------------
# Dosing specifications
  dose <- 50	# mg
# Create sequence of dosing times based on "on" and "off" periods in cycles
  ncycles <- 1
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
# Perform initial PK simulations (first week)
# Input data frame for simulation
  pk.ID.data <- data.frame(ID = ID.seq,WT = WT,pk.ETA.matrix)
  input.pk.data <- lapply(pk.ID.data,rep.int,times = length(pk.times)) %>%
    as.data.frame
  input.pk.data <- input.pk.data[with(input.pk.data,order(input.pk.data$ID)),]
  input.pk.data$time <- pk.times
  input.pk.data$cyc <- cyc
  input.pk.data$amt <- dose
  input.pk.data$cmt <- 1
  input.pk.data$evid <- 0
  input.pk.data$evid[input.pk.data$time %in% on.times] <- 1

# ------------------------------------------------------------------------------
# Simulate concentration-time profile for the population given standard dosing
  pop.pk.data <- pk.mod %>% mrgsim(data = input.pk.data,
    carry.out = c("amt")) %>% as.data.frame
  pop.pk.data <- auc24.function(pop.pk.data)
# Pull out each individual's true parameter values
  true.param.function <- function(pop.pk.data) {
    true.param <- data.frame(CLP = pop.pk.data$CLP[1],
      VCP = pop.pk.data$VCP[1],
      CLM = pop.pk.data$CLM[1],
      VCM = pop.pk.data$VCM[1])
  }
  true.param <- ddply(pop.pk.data, .(ID), true.param.function)

# ------------------------------------------------------------------------------
# Estimate empirical Bayes parameters and simulate resulting concentration-time
# profile.  Option for using sampled parent, sampled metabolite or combined.
  input.bayes.data <- input.pk.data
  bayesian.function <- function(input.bayes.data) {
    ID <- input.bayes.data$ID[1]
    pk.data <- pop.pk.data[pop.pk.data$ID == ID,]
  # Sample concentrations
    if (conc == "parent") {
      sample.DV <- pk.data$DVP[pk.data$time %in% c(sample.time*24)]
    }
    if (conc == "metabolite") {
      sample.DV <- pk.data$DVM[pk.data$time %in% c(sample.time*24)]
    }
    if (conc == "combined") {
      sample.DV <- pk.data$DV[pk.data$time %in% c(sample.time*24)]
    }
    if (conc == "separate") {
      sample.DVP <- pk.data$DVP[pk.data$time %in% c(sample.time*24)]
      sample.DVM <- pk.data$DVM[pk.data$time %in% c(sample.time*24)]
    }
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
        if (conc == "parent") {
        # Pull out the predicted trough concentrations at sample time
          yhat <- bayes.data$IPREP[bayes.data$time %in% c(sample.time*24)]
        # Posterior log-likelihood
        # Error model: Y = IPRE*(1+ERRPRO), Y = IPRE + IPRE*ERRPRO
          loglikpost.sd <- ERRPROP1
          loglikpost <- dnorm(sample.DV,mean = yhat,sd = yhat*loglikpost.sd,
            log = TRUE)
        }
        if (conc == "metabolite") {
        # Pull out the predicted trough concentrations at sample time
          yhat <- bayes.data$IPREM[bayes.data$time %in% c(sample.time*24)]
        # Posterior log-likelihood
        # Error model: Y = IPRE*(1+ERRPRO), Y = IPRE + IPRE*ERRPRO
          loglikpost.sd <- ERRPROM1
          loglikpost <- dnorm(sample.DV,mean = yhat,sd = yhat*loglikpost.sd,
            log = TRUE)
        }
        if (conc == "combined") {
        # Pull out the predicted trough concentrations at sample time
          yhat <- bayes.data$IPRE[bayes.data$time %in% c(sample.time*24)]
        # Posterior log-likelihood
        # Error model: Y = IPRE*(1+ERRPRO), Y = IPRE + IPRE*ERRPRO
          loglikpost.sd <- ERRPROP1
          loglikpost <- dnorm(sample.DV,mean = yhat,sd = yhat*loglikpost.sd,
            log = TRUE)
        }
        if (conc == "separate") {
        # Pull out the predicted trough concentrations at sample time
          yhatp <- bayes.data$IPREP[bayes.data$time %in% c(sample.time*24)]
          yhatm <- bayes.data$IPREM[bayes.data$time %in% c(sample.time*24)]
        # Posterior log-likelihood
        # Error model: Y = IPRE*(1+ERRPRO), Y = IPRE + IPRE*ERRPRO
          loglikpost.sdp <- ERRPROP1
          loglikpost.sdm <- ERRPROM1
          loglikpostp <- dnorm(sample.DVP,mean = yhatp,sd = yhatp*loglikpost.sdp,
            log = TRUE)
          loglikpostm <- dnorm(sample.DVM,mean = yhatm,sd = yhatm*loglikpost.sdm,
            log = TRUE)
          loglikpost <- c(loglikpostp,loglikpostm)
        }
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
    bayes.result <- data.frame(par = bayes.result$par)
  }

# Estimate Bayes parameters
  bayes.result <- ddply(input.bayes.data, .(ID), bayesian.function,
    .progress = "text")

# ------------------------------------------------------------------------------
# Calculate the parameter values for each individual (only ETA's were estimated)
  bayes.param.function <- function(bayes.result) {
    ID <- bayes.result$ID[1]
    WT <- head(input.bayes.data$WT[input.bayes.data$ID == ID],1)
    bayes.param <- data.frame(
      CLP = param(pk.mod)$POPCLP*bayes.result$par[1]*(WT/70)^0.75,
      VCP = param(pk.mod)$POPVCP*bayes.result$par[2]*(WT/70),
      CLM = param(pk.mod)$POPCLM*bayes.result$par[3]*(WT/70)^0.75,
      VCM = param(pk.mod)$POPVCM*bayes.result$par[4]*(WT/70)
    )
  }
  bayes.param <- ddply(bayes.result, .(ID), bayes.param.function)

# Calculate the error in Bayes parameters compared to the individual's true
# parameters
  param.err.function <- function(true.param) {
    ID <- true.param$ID[1]
    ind.bayes.param <- bayes.param[bayes.param$ID == ID,]
    abs.err <- true.param[-1]-ind.bayes.param[-1]
    rel.err <- true.param[-1]/ind.bayes.param[-1]*100
    names(abs.err) <- c("CLPabs","VCPabs","CLMabs","VCMabs")
    names(rel.err) <- c("CLPrel","VCPrel","CLMrel","VCMrel")
    param.err <- data.frame(abs.err,rel.err)
  }
  param.err <- ddply(true.param, .(ID), param.err.function)

# Summarise into median and confidence intervals
  melt.param.err <- melt(param.err,id.vars = "ID")
  melt.param.err$variable <- as.factor(melt.param.err$variable)
  summary.param.err <- ddply(melt.param.err, .(variable),
    function(melt.param.err) summary.function(melt.param.err$value))

# ------------------------------------------------------------------------------
# Plot box and whisker plots for error in parameter estimates
# Absolute error
  plotobj1 <- NULL
  plotobj1 <- ggplot(melt.param.err[
    contains("abs",vars = melt.param.err$variable),])
  plotobj1 <- plotobj1 + geom_boxplot(aes(x = variable,y = value,
    colour = variable),outlier.shape = NA)
  plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 0),linetype = "dashed")
  plotobj1 <- plotobj1 + scale_y_continuous("Absolute Error",
    lim = c(NA,1300))
  plotobj1 <- plotobj1 + scale_x_discrete("Pharmacokinetic Parameter")
  plotobj1 <- plotobj1 + theme(legend.position = "none")
  print(plotobj1)

  ggsave(plot = plotobj1,filename = paste0(output.dir.name,"_abserror.png"),
    width = 20,height = 15,unit = "cm",dpi = 300)

# Relative error
  plotobj2 <- NULL
  plotobj2 <- ggplot(melt.param.err[
    contains("rel",vars = melt.param.err$variable),])
  plotobj2 <- plotobj2 + geom_boxplot(aes(x = variable,y = value,
    colour = variable),outlier.shape = NA)
  plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 100),linetype = "dashed")
  plotobj2 <- plotobj2 + scale_y_continuous("Relative Error (%)",
    lim = c(0,250))
  plotobj2 <- plotobj2 + scale_x_discrete("Pharmacokinetic Parameter")
  plotobj2 <- plotobj2 + theme(legend.position = "none")
  print(plotobj2)

  ggsave(plot = plotobj2,filename = paste0(output.dir.name,"_relerror.png"),
    width = 20,height = 15,unit = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Simulate Bayes predicted profiles
  sim.bayes.function <- function(bayes.result) {
    ID <- bayes.result$ID[1]
    ind.input.bayes.data <- input.bayes.data[input.bayes.data$ID == ID,]
    ind.input.bayes.data$ETACLP <- log(bayes.result$par[1])
    ind.input.bayes.data$ETAVCP <- log(bayes.result$par[2])
    ind.input.bayes.data$ETACLM <- log(bayes.result$par[3])
    ind.input.bayes.data$ETAVCM <- log(bayes.result$par[4])
    bayes.data <- pk.mod %>% mrgsim(data = ind.input.bayes.data) %>%
      as.data.frame
  }
  bayes.data <- ddply(bayes.result, .(ID), sim.bayes.function)

# ------------------------------------------------------------------------------
# Plot concentration-time profile for 12 random individuals
  set.seed(123456)
  ID.rand <- sample(unique(ID.seq),12)
  rand.pop.pk.data <- pop.pk.data[pop.pk.data$ID %in% ID.rand,]
  rand.bayes.data <- bayes.data[bayes.data$ID %in% ID.rand,]

  plotobj3 <- NULL
  plotobj3 <- ggplot(rand.pop.pk.data)
  plotobj3 <- plotobj3 + geom_line(aes(x = time/24,y = IPRE),
    data = rand.pop.pk.data,colour = "red")
  plotobj3 <- plotobj3 + geom_line(aes(x = time/24,y = IPRE),
    data = rand.bayes.data,colour = "blue")
  plotobj3 <- plotobj3 + scale_y_continuous("Total Sunitinib Concentration (mg/L)")
  plotobj3 <- plotobj3 + scale_x_continuous("Time (days)")
  plotobj3 <- plotobj3 + facet_wrap(~ID,ncol = 4)
  print(plotobj3)

  ggsave(plot = plotobj3,filename = paste0(output.dir.name,"_randpred.png"),
    width = 20,height = 15,unit = "cm",dpi = 300)
