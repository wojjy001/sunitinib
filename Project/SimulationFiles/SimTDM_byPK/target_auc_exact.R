# Simulate PK following standard 50 mg dosing for weeks
# Then sample trough at the end of the first cycle
# Estimate individual PK parameters using Bayes' theorem
# Calculate new dose based on estimate of AUC, target AUC of 1.5 mg*h/L and
# previous dose and administer for the second (and subsequent) cycle
# Do NOT round doses to the nearest 12.5 mg (tablet sizes)
# ------------------------------------------------------------------------------
# Set working directory to project output directory
  setwd(output.dir)
  # Create a folder in the directory specific to this dosing regimen
    study.dir <- paste0(output.dir,"/",study.name)
    dir.create(study.dir)
    setwd(study.dir)
# For each individual:
  target.auc.exact <- function(input.pk.data) {
  # Simulate the first cycle given then 50 mg dose
    cyc1.data <- pk.mod %>%
      mrgsim(data = input.pk.data[input.pk.data$cyc == 1,],
      carry.out = c("SIM","amt","cyc")) %>% as.data.frame
  # Sample total measured sunitinib concentration at day 28
    sample.time <- 28*24	# hours
    sample.DV <- cyc1.data$DV[cyc1.data$time == sample.time]
  # Estimate individual PK parameters using Bayes' theorem
    # Initial estimates for Bayes' parameters
    # Use population typical values the first time Bayes' parameters are
    # estimated
    initial.bayes.par <- exp(c(0,0,0,0))	# exponential domain
    run.once <- FALSE	# Flag describing if the bayes.estimate function has been
    # run at least once for the individual
    repeat {	# Repeat Bayes' estimation until optim converges successfully
    # Use previous parameter results as initial estimates mulitpled by a
    # random number so that they are not exactly the same
      if (run.once == TRUE) {
        initial.bayes.par <- initial.bayes.par*exp(runif(4,min = -0.01,
          max = 0.01))
      }
      par.bayes <- initial.bayes.par
      input.bayes.data <- input.pk.data[input.pk.data$cyc == 1,]
    # Bayesian estimation using population PK model
    # Maximum likelihood estimation for Bayes' parameters
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
      # Pull out the predicted trough concentrations at sample time
        yhat <- bayes.data$IPRE[bayes.data$time == sample.time]
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
  # Simulate Bayes' estimated profile
    # Add Bayes' estimated parameters to input data frame for simulation
      input.bayes.data$ETACLP <- log(bayes.result$par[1])
      input.bayes.data$ETAVCP <- log(bayes.result$par[2])
      input.bayes.data$ETACLM <- log(bayes.result$par[3])
      input.bayes.data$ETAVCM <- log(bayes.result$par[4])
    # Simulate
      bayes.data <- pk.mod %>% mrgsim(data = input.bayes.data) %>%
        as.data.frame
  # Calculate Bayes' estimated 24-hour at time of sampling
    bayes.data <- auc24.function(bayes.data)
    sample.AUC24 <- bayes.data$AUC24[bayes.data$time == sample.time]
  # Calculate new dose for the individual
    new.dose <- 1.5/sample.AUC24*dose	# 1.5 is the target AUC24, dose = 50 mg
  # Add "new.dose" to the input data frame
    input.pk.data$amt[input.pk.data$cyc > 1] <- new.dose
  # Simulate
    pk.data <- pk.mod %>% mrgsim(data = input.pk.data,
      carry.out = c("SIM","amt","cyc")) %>% as.data.frame
    pk.data <- auc24.function(pk.data)	# Calculate AUC24 at each time-point
  }
# Simulate for the population
  pk.data <- ddply(input.pk.data, .(SIM,ID), target.auc.exact,
    .progress = "text")
  write.csv(pk.data,file = paste0(study.name,"_pk_data.csv"),
    quote = FALSE,row.names = FALSE)
