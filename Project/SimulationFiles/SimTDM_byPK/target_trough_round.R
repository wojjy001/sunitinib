# Simulate PK following standard 50 mg dosing for weeks
# Then sample trough at the end of the first cycle
# Calculate new dose based on sampled trough, target trough of 50 ng/mL and
# previous dose and administer for the second (and subsequent) cycle
# ROUND doses to the nearest 12.5 mg (tablet sizes)
# ------------------------------------------------------------------------------
# Set working directory to project output directory
  setwd(output.dir)
  # Create a folder in the directory specific to this dosing regimen
    study.name <- "target_trough_round"
    study.dir <- paste0(output.dir,"/",study.name)
    dir.create(study.dir)
    setwd(study.dir)

# ------------------------------------------------------------------------------
# For each individual:
  target.trough.round <- function(input.pk.data) {
  # Simulate the first cycle given then 50 mg dose
    cyc1.data <- pk.mod %>%
      mrgsim(data = input.pk.data[input.pk.data$cyc == 1,],
      carry.out = c("SIM","amt","cyc")) %>% as.data.frame
  # Sample total measured sunitinib concentration at day 28
    sample.time <- 28*24	# hours
    sample.DV <- cyc1.data$DV[cyc1.data$time == sample.time]
  # Calculate new dose for the individual
    new.dose <- 50/sample.DV*dose	# Where 50 is the target trough, dose = 50 mg
    new.dose <- round(new.dose/12.5,digits = 0)*12.5
  # Add "new.dose" to the input data frame
    input.pk.data$amt[input.pk.data$cyc > 1] <- new.dose
  # Simulate
    pk.data <- pk.mod %>% mrgsim(data = input.pk.data,
      carry.out = c("SIM","amt","cyc")) %>% as.data.frame
    pk.data <- auc24.function(pk.data)	# Calculate AUC24 at each time-point
  }
# Simulate for the population
  pk.data <- ddply(input.pk.data, .(SIM,ID), target.trough.round)
  write.csv(pk.data,file = paste0(study.name,"_pk_data.csv"),
    quote = FALSE,row.names = FALSE)
