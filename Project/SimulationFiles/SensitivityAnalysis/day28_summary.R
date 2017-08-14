# Sensitivity Analysis
# Read in previously simulated data and summarise dose, IPRE and AUC24 at Day 28
# Day 28 is the end of the on-period of the first cycle
# This was the time for dose optimisation to achieve target AUC or trough
# concentration
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
# Set output directory for sensitivity analysis
  sens.output.dir <- paste0(output.dir,"SensitivityAnalysis/")
  setwd(sens.output.dir)

# ------------------------------------------------------------------------------
# Read in simulation files in the output directory
# List the folders in the output directory
  folder.list <- list.dirs(sens.output.dir)
  folder.list <- folder.list[contains("target",vars = folder.list)]
  study.list <- folder.list
  split.list <- folder.list
  type.list <- folder.list
  target.list <- folder.list
  weight.list <- folder.list
  for (i in 1:length(folder.list)) {
    folder.list[i] <- str_split(folder.list[i],pattern = "/")
    study.list[i] <- folder.list[[i]][8]
    split.list[i] <- str_split(study.list[i],pattern = "_")
    type.list[i] <- split.list[[i]][2]
    target.list[i] <- split.list[[i]][3]
    weight.list[i] <- split.list[[i]][4]
  }
  folder.data <- data.frame(study = study.list,type = type.list,
    target = target.list,weight = weight.list)
# Read in .csv files from each folder
# Day28 data
  read.day28.data <- function(folder.data) {
    study <- as.character(folder.data$study[1])
    type <- as.character(folder.data$type[1])
    target <- as.character(folder.data$target[1])
    weight <- as.character(folder.data$weight[1])
    day28.data <- read.csv(file = paste0(sens.output.dir,weight,"/",study,
      "/target_",type,"_",target,"_",weight,"_day28_data.csv"))
    day28.data$study <- study
    day28.data$type <- type
    day28.data$target <- target
    day28.data
  }
  day28.data <- ddply(folder.data, .(study), read.day28.data,
    .progress = "text")
# Dose data
  read.dose.data <- function(folder.data) {
    study <- as.character(folder.data$study[1])
    type <- as.character(folder.data$type[1])
    target <- as.character(folder.data$target[1])
    weight <- as.character(folder.data$weight[1])
    dose.data <- read.csv(file = paste0(sens.output.dir,weight,"/",study,
      "/target_",type,"_",target,"_",weight,"_dose_data.csv"))
    dose.data$study <- study
    dose.data$type <- type
    dose.data$target <- target
    dose.data
  }
  dose.data <- ddply(folder.data, .(study), read.dose.data,
    .progress = "text")
# Merge day28.data and dose.data
  sim.data <- merge(day28.data,dose.data,
    by = c("SIM","ID","study","type","target"),all = T)
  sim.data$targetf <- as.factor(sim.data$target)
  sim.data$WTf <- as.factor(sim.data$WT)

# ------------------------------------------------------------------------------
# Plot weight versus dose by IPRE target
  plotobj1 <- NULL
  plotobj1 <- ggplot(sim.data[sim.data$type == "auc",])
  plotobj1 <- plotobj1 + geom_boxplot(aes(x = WTf,y = dose,colour = WTf),
    outlier.shape = 1)
  plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 50),
    linetype = "dashed")
  plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 37.5),
    linetype = "dashed",colour = "grey")
  plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 62.5),
    linetype = "dashed",colour = "grey")
  plotobj1 <- plotobj1 + scale_y_continuous("Sunitinib Dose (mg)",
    breaks = seq(from = 0,to = max(sim.data$dose),by = 25),
    labels = seq(from = 0,to = max(sim.data$dose),by = 25))
  plotobj1 <- plotobj1 + scale_x_discrete("Weight (kg)")
  plotobj1 <- plotobj1 + facet_wrap(~targetf)
  plotobj1 <- plotobj1 + theme(legend.position = "none")
  print(plotobj1)
  ggsave(plot = plotobj1,filename = paste0(sim.dir,"/dosevsWTvsAUC24target.png"),
    height = 15,width = 20,units = "cm",dpi = 300)

# Plot dose versus IPRE
  plotobj2 <- NULL
  plotobj2 <- ggplot(sim.data[sim.data$type == "trough",])
  plotobj2 <- plotobj2 + geom_boxplot(aes(x = WTf,y = IPRE,colour = WTf),
    outlier.shape = 1)
  plotobj2 <- plotobj2 + geom_hline(aes(yintercept = min(sim.data$IPRE)),
    linetype = "dashed",colour = "grey")
  plotobj2 <- plotobj2 + geom_hline(aes(yintercept = max(sim.data$IPRE)),
    linetype = "dashed",colour = "grey")
  plotobj2 <- plotobj2 + geom_hline(aes(yintercept = median(sim.data$IPRE)),
    linetype = "dashed")
  plotobj2 <- plotobj2 + scale_y_continuous("Total Trough Sunitinib Concentration (mg/L)",
    lim = c(0,NA))
  plotobj2 <- plotobj2 + scale_x_discrete("Weight (kg)")
  plotobj2 <- plotobj2 + facet_wrap(~targetf)
  plotobj2 <- plotobj2 + theme(legend.position = "none")
  print(plotobj2)
  ggsave(plot = plotobj2,filename = paste0(sim.dir,"/WTvsIPREvsAUC24target.png"),
    height = 15,width = 20,units = "cm",dpi = 300)
