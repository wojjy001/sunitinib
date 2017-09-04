# Simulate PK following doses administered as 50 mg/70 kg * WT (kg)
# ROUND doses to the nearest 12.5 mg
# ------------------------------------------------------------------------------
# Set working directory to project output directory
  setwd(output.dir)
  # Create a folder in the directory specific to this dosing regimen
    study.dir <- paste0(output.dir,"/",study.name)
    dir.create(study.dir)
    setwd(study.dir)
# Calculate mg/kg dose for each individual and round to nearest 12.5 mg
  input.pk.data$amt <- round((50/70)*input.pk.data$WT/12.5,digits = 0)*12.5
# Simulate PK for the population
  pk.data <- pk.mod %>%
    mrgsim(data = input.pk.data,carry.out = c("SIM","amt","cyc")) %>%
    as.data.frame
# Calculate AUC24 at each time-point
  pk.data <- ddply(pk.data, .(SIM,ID), auc24.function)
  write.csv(pk.data,file = paste0(study.name,"_pk_data.csv"),
    quote = FALSE,row.names = FALSE)
