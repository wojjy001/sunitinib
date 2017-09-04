# Simulate PK following doses administered as 50 mg/70 kg * WT (kg)
# DO NOT round doses to the nearest 12.5 mg
# Cap doses between 12.5 and 87.5 mg
# Make all individuals have a weight of 70 kg
# ------------------------------------------------------------------------------
# Set working directory to project output directory
  setwd(output.dir)
  # Create a folder in the directory specific to this dosing regimen
    study.dir <- paste0(output.dir,"/",study.name)
    dir.create(study.dir)
    setwd(study.dir)
# Make all individuals have a weight of 70 kg
  input.pk.data$WT <- 70	# kg
# Calculate mg/kg dose for each individual
  input.pk.data$amt <- (50/70)*input.pk.data$WT
  input.pk.data$amt[input.pk.data$amt > 87.5] <- 87.5
  input.pk.data$amt[input.pk.data$amt < 12.5] <- 12.5
# Simulate PK for the population
  pk.data <- pk.mod %>%
    mrgsim(data = input.pk.data,carry.out = c("SIM","amt","cyc")) %>%
    as.data.frame
# Calculate AUC24 at each time-point
  pk.data <- ddply(pk.data, .(SIM,ID), auc24.function)
  write.csv(pk.data,file = paste0(study.name,"_pk_data.csv"),
    quote = FALSE,row.names = FALSE)
