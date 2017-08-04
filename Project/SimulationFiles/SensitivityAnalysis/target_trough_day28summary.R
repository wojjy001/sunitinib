# Sunitinib Simulation
# Read in previously simulated data and summarise
# ------------------------------------------------------------------------------
# Remove objects in workspace
  rm(list = ls(all = TRUE))
# Define global directory for project
  global.dir <- "/Volumes/Prosecutor/sunitinib/Project/"
# Source required files
  source(paste0(global.dir,"functions.R"))	# Universal functions file
# Set output directory (folder outside of git folder)
  sim.dir <- paste0("/Volumes/Prosecutor/sunitinib_nogit/SensitivityAnalysis")

# Read in simulation files in the output directory
# List the folders in the output directory
  folder.list <- list.dirs(sim.dir)
  folder.list <- folder.list[contains("target_trough",vars = folder.list)]
  weight.list <- folder.list
  study.list <- folder.list
  for (i in 1:length(folder.list)) {
    folder.list[i] <- str_split(folder.list[i],pattern = "/")
    weight.list[i] <- folder.list[[i]][6]
    study.list[i] <- folder.list[[i]][7]
  }
  folder.data <- data.frame(weight = weight.list,
    study = study.list)
# Read in .csv files from each folder
  # Day28 data
    read.day28.data <- function(folder.data) {
      weight <- as.character(folder.data$weight[1])
      study <- as.character(folder.data$study[1])
      day28.data <- read.csv(file = paste0(sim.dir,"/",weight,"/",study,
        "/",study,"_day28_data.csv"))
      day28.data$study <- study
      day28.data
    }
    day28.data <- ddply(folder.data, .(study), read.day28.data)
  # Dose data
    read.dose.data <- function(folder.data) {
      weight <- as.character(folder.data$weight[1])
      study <- as.character(folder.data$study[1])
      dose.data <- read.csv(file = paste0(sim.dir,"/",weight,"/",study,
        "/",study,"_dose_data.csv"))
      dose.data$study <- study
      dose.data
    }
    dose.data <- ddply(folder.data, .(study), read.dose.data)
# Merge day28.data and dose.data
  sim.data <- merge(day28.data,dose.data,by = c("SIM","ID","study"),all = T)

# ------------------------------------------------------------------------------
# Plot weight versus dose by IPRE target
  sim.data$IPREf <- as.factor(round(sim.data$IPRE/0.025)*0.025)
  levels(sim.data$IPREf) <- c("Trough = 50 ng/mL",
    "Trough = 75 ng/mL",
    "Trough = 100 ng/mL")
  sim.data$WTf <- as.factor(sim.data$WT)
  plotobj1 <- NULL
  plotobj1 <- ggplot(sim.data)
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
  plotobj1 <- plotobj1 + facet_wrap(~IPREf)
  plotobj1 <- plotobj1 + theme(legend.position = "none")
  print(plotobj1)
  ggsave(plot = plotobj1,filename = paste0(sim.dir,"/dosevsWTvsIPREtarget.png"),
    height = 15,width = 20,units = "cm",dpi = 300)

# Plot dose versus IPRE
  plotobj2 <- NULL
  plotobj2 <- ggplot(sim.data)
  plotobj2 <- plotobj2 + geom_boxplot(aes(x = WTf,y = AUC24,colour = WTf),
    outlier.shape = 1)
  plotobj2 <- plotobj2 + geom_hline(aes(yintercept = min(sim.data$AUC24)),
    linetype = "dashed",colour = "grey")
  plotobj2 <- plotobj2 + geom_hline(aes(yintercept = max(sim.data$AUC24)),
    linetype = "dashed",colour = "grey")
  plotobj2 <- plotobj2 + geom_hline(aes(yintercept = median(sim.data$AUC24)),
    linetype = "dashed")
  plotobj2 <- plotobj2 + scale_y_continuous("Total Sunitinib 24-hour AUC (mg*h/L)",
    lim = c(0,NA))
  plotobj2 <- plotobj2 + scale_x_discrete("Weight (kg)")
  plotobj2 <- plotobj2 + facet_wrap(~IPREf)
  plotobj2 <- plotobj2 + theme(legend.position = "none")
  print(plotobj2)
  ggsave(plot = plotobj2,filename = paste0(sim.dir,"/WTvsAUC24vsIPREtarget.png"),
    height = 15,width = 20,units = "cm",dpi = 300)
