# Simulate PK following standard 50 mg dosing and 4/2 schedule
# ------------------------------------------------------------------------------
# Set working directory to project output directory
  setwd(output.dir)
  # Create a folder in the directory specific to this dosing regimen
    study.name <- "standard"
    study.dir <- paste0(output.dir,"/",study.name)
    dir.create(study.dir)
    setwd(study.dir)

# ------------------------------------------------------------------------------
# Simulate PK for the population
  pk.data <- pk.mod %>%
    mrgsim(data = input.pk.data,carry.out = c("SIM","amt","cyc")) %>%
    as.data.frame
  pk.data <- auc24.function(pk.data)	# Calculate AUC24 at each time-point
  write.csv(pk.data,file = paste0(study.name,"_pk_data.csv"),
    quote = FALSE,row.names = FALSE)
