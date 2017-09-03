# Simulate PK following standard 50 mg dosing for weeks
# Then sample AUC24 at the end of the first cycle
# Calculate new dose based on sampled AUC24, target AUC of 1.5 mg*h/L and
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
  target.auc.auc24 <- function(input.pk.data) {
  # Simulate the first cycle given then 50 mg dose
    cyc1.data <- pk.mod %>%
      mrgsim(data = input.pk.data[input.pk.data$cyc == 1,],
      carry.out = c("SIM","amt","cyc")) %>% as.data.frame
    cyc1.data <- auc24.function(cyc1.data)	# Calculate AUC24 at each time-point
  # Sample total measured sunitinib concentration at day 28
    sample.time <- 28*24	# hours
    sample.AUC24 <- cyc1.data$AUC24[cyc1.data$time == sample.time]
  # Calculate new dose for the individual
    new.dose <- 1.5/sample.AUC24*dose	# Where 1.5 is the target AUC, dose = 50 mg
  # Add "new.dose" to the input data frame
    input.pk.data$amt[input.pk.data$cyc > 1] <- new.dose
  # Simulate
    pk.data <- pk.mod %>% mrgsim(data = input.pk.data,
      carry.out = c("SIM","amt","cyc")) %>% as.data.frame
    pk.data <- auc24.function(pk.data)	# Calculate AUC24 at each time-point
  }
# Simulate for the population
  pk.data <- ddply(input.pk.data, .(SIM,ID), target.auc.auc24,
    .progress = "text")
  write.csv(pk.data,file = paste0(study.name,"_pk_data.csv"),
    quote = FALSE,row.names = FALSE)
