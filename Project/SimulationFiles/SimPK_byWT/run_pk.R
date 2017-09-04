# Load packages required to run all scripts
# Assign directories (simulate files and simulation output directories)
# Run PK simulation files
# ------------------------------------------------------------------------------
# Remove objects in workspace
  rm(list = ls(all = TRUE))
# Define simulation file directory for the project
  file.dir <- "/Volumes/Prosecutor/sunitinib/Project/SimulationFiles/"
  project.name <- "SimPK_byWT"
  project.dir <- paste0(file.dir,project.name)
  setwd(project.dir)
# Define output directory
# Output is saved into a totally different folder from scripts because my code
# is tracked by git and stored online.  Keeping all output (i.e., .csv and .png)
# files in the git repository slows down synchronisation
  output.dir <- paste0("/Volumes/Prosecutor/sunitinib_nogit/",project.name)
  dir.create(output.dir)	# Create the folder if not done so

# ------------------------------------------------------------------------------
# Source and run files required prior to PK simulation
# Functions file (package libraries, summary functions, plotting)
  source("functions.R")
# Population pharmacokinetic model of sunitinib parent and metabolite
  source("sunitinib_pk_model.R")
# Define PK and PD simulation times
  source("times.R")
# Create population of unique individuals
  dose <- 50	# mg, initial dose to be administered to all individuals
  source("population.R") # Resulting data frame is called "input.pk.data"
  # Ignore the warnings associated with:
  # In log(rnorm(length(new.OBASE), mean = 195, sd = 120)) : NaNs produced

# ------------------------------------------------------------------------------
# PK simulation
# Source the dosing regimen file
  study.name <- "standard_110kg"
  source(paste0(study.name,".R"))	# Resulting data frame is called "pk.data"

# ------------------------------------------------------------------------------
# Plot results
# Summarise median and confidence intervals for IPRE (total sunitinib) and AUC24
  summary.IPRE <- ddply(pk.data, .(time),
    function(pk.data) graded.summary(pk.data$IPRE))
  summary.AUC24 <- ddply(pk.data, .(time),
    function(pk.data) graded.summary(pk.data$AUC24))
# Plot IPRE
  plotobj1 <- NULL
  plotobj1 <- ggplot(summary.IPRE)
  plotobj1 <- plotobj1 + geom_ribbon(aes(x = time/24,ymin = CI90lo,
    ymax = CI90hi),fill = "brown3",alpha = 0.2)
  plotobj1 <- plotobj1 + geom_ribbon(aes(x = time/24,ymin = CI80lo,
    ymax = CI80hi),fill = "brown3",alpha = 0.2)
  plotobj1 <- plotobj1 + geom_ribbon(aes(x = time/24,ymin = CI60lo,
    ymax = CI60hi),fill = "brown3",alpha = 0.2)
  plotobj1 <- plotobj1 + geom_ribbon(aes(x = time/24,ymin = CI40lo,
    ymax = CI40hi),fill = "brown3",alpha = 0.2)
  plotobj1 <- plotobj1 + geom_ribbon(aes(x = time/24,ymin = CI20lo,
    ymax = CI20hi),fill = "brown3",alpha = 0.2)
  plotobj1 <- plotobj1 + geom_line(aes(x = time/24,y = med),
    colour = "brown3")
  plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 50),
    linetype = "dashed")
  plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 100),
    linetype = "dashed")
  plotobj1 <- plotobj1 + scale_y_continuous("Total Sunitinib Concentration (ng/mL)",
    lim = c(0,200),
    breaks = c(0,25,50,75,100,125,150,175,200),
    labels = c(0,25,50,75,100,125,150,175,200))
  plotobj1 <- plotobj1 + scale_x_continuous("Time (days)",
    breaks = seq(from = 0,to = max(pk.times),by = 14),
    labels = seq(from = 0,to = max(pk.times),by = 14))
  print(plotobj1)
  ggsave(plot = plotobj1,filename = paste0(study.name,"_IPREvstime.png"),
    height = 15,width = 20,units = "cm",dpi = 300)
# Plot AUC24
  plotobj2 <- NULL
  plotobj2 <- ggplot(summary.AUC24)
  plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24,ymin = CI90lo,
    ymax = CI90hi),fill = "skyblue4",alpha = 0.2)
  plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24,ymin = CI80lo,
    ymax = CI80hi),fill = "skyblue4",alpha = 0.2)
  plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24,ymin = CI60lo,
    ymax = CI60hi),fill = "skyblue4",alpha = 0.2)
  plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24,ymin = CI40lo,
    ymax = CI40hi),fill = "skyblue4",alpha = 0.2)
  plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24,ymin = CI20lo,
    ymax = CI20hi),fill = "skyblue4",alpha = 0.2)
  plotobj2 <- plotobj2 + geom_line(aes(x = time/24,y = med),
    colour = "skyblue4")
  plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 1.5),
    linetype = "dashed")
  plotobj2 <- plotobj2 + scale_y_continuous("Total Sunitinib 24-hour AUC (mg*h/L)",
    lim = c(0,NA))
  plotobj2 <- plotobj2 + scale_x_continuous("Time (days)",
    breaks = seq(from = 0,to = max(pk.times),by = 14),
    labels = seq(from = 0,to = max(pk.times),by = 14))
  print(plotobj2)
  ggsave(plot = plotobj2,filename = paste0(study.name,"_AUC24vstime.png"),
    height = 15,width = 20,units = "cm",dpi = 300)
# Plot IPRE versus AUC24
  plotobj3 <- NULL
  plotobj3 <- ggplot(pk.data[pk.data$time == max(on.times),])
  plotobj3 <- plotobj3 + geom_point(aes(x = IPRE,y = AUC24),
    colour = "skyblue4",shape = 1,size = 2)
  plotobj3 <- plotobj3 + scale_y_continuous("Total Sunitinib 24-hour AUC (mg*h/L)",
    lim = c(0,NA))
  plotobj3 <- plotobj3 + scale_x_continuous("Total Sunitinib Concentration (ng/mL)",
    lim = c(0,NA))
  print(plotobj3)
  ggsave(plot = plotobj3,filename = paste0(study.name,"_AUC24vsIPRE.png"),
    height = 15,width = 20,units = "cm",dpi = 300)
