# Find the dose for each individual that results in a target 24-hour AUC
# Simulate the resulting profiles that will be provided as input for PD
# ------------------------------------------------------------------------------
# Set working directory to project output directory
  setwd(output.dir)
  # Create a folder in the directory specific to this dosing regimen
    study.dir <- paste0(output.dir,"/",study.name)
    dir.create(study.dir)
    setwd(study.dir)
# Optimise the dose for each individual that will achieve the target auc
# concentration by Day 28 of the first cycle of treatment
  pop.optimise.dose <- function(input.pk.data) {
  # Subset the input data frame for only the first cycle
  # Optimisation takes too long if all cycles are included
    input.optim.data <- input.pk.data[input.pk.data$cyc == 1,]
  # Initial parameter estimates
    initial.dose <- dose	# standard 50 mg dose
    initial.err <- 0.01
    par <- c(initial.dose,initial.err)
  # Optimise dose function
    optimise.dose <- function(par) {
    # Assign estimable parameters to objects
      input.optim.data$amt[input.optim.data$evid == 1] <- par[1]
      err <- par[2]
    # Simulate concentration-time profile with each iteration of dose
      optim.data <- pk.mod %>% mrgsim(data = input.optim.data) %>%
        as.data.frame
      optim.data <- auc24.function(optim.data)
    # Pull out the predicted AUC at the target time
      yhat <- optim.data$AUC24[optim.data$time == 28*24]
    # Find the value of dose that maximises the likelihood of the predicted
    # concentration being the target AUC
      loglik <- dnorm(target,yhat,yhat*err,log = TRUE)
    # Define the objective function value that will be optimised by "optim"
      objective <- -1*sum(loglik)
    }
  # Run the optim function to obtain a value for dose
    optimised.dose <- optim(par,
      optimise.dose,
      hessian = FALSE,
      method = "L-BFGS-B",
      lower = c(0.0001,0.0001),upper = c(Inf,Inf)
    )
  # Create output object for each individual with their individual dose
    dose.data <- data.frame(dose = optimised.dose$par[1])
  }	# pop.optimise.dose
# Run "pop.optimise.dose" for each individual in the population
  dose.data <- ddply(input.pk.data, .(SIM,ID), pop.optimise.dose,
    .progress = "text")
  write.csv(dose.data,file = paste0(study.name,"_dose_data.csv"),
    quote = FALSE,row.names = FALSE)
# Merge dose.data into input.pk.data
  add.doses <- function(dose.data) {
    ID <- dose.data$ID[1]
    SIM <- dose.data$SIM[1]
    input.pk.data <- input.pk.data[input.pk.data$SIM == SIM &
      input.pk.data$ID == ID,]
    input.pk.data$amt[input.pk.data$evid == 1] <- dose.data$dose[1]
    input.pk.data
  }
  input.pk.data <- ddply(dose.data, .(SIM,ID), add.doses)
# Simulate PK for the population
  pk.data <- pk.mod %>%
    mrgsim(data = input.pk.data,carry.out = c("SIM","amt","cyc")) %>%
    as.data.frame
# Calculate AUC24 at each time-point
  pk.data <- ddply(pk.data, .(SIM,ID), auc24.function)
  write.csv(pk.data,file = paste0(study.name,"_pk_data.csv"),
    quote = FALSE,row.names = FALSE)
