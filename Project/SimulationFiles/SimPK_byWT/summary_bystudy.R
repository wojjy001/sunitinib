# Read in simulated pharmacodynamic data for studies
# Summarise simulation results against weight categories
# ------------------------------------------------------------------------------
# Remove objects in workspace
  rm(list = ls(all = TRUE))
# Define simulation file directory for the project
  file.dir <- "/Volumes/Prosecutor/sunitinib/Project/SimulationFiles/"
  project.name <- "SimPK_byWT"
  project.dir <- paste0(file.dir,project.name)
  setwd(project.dir)
# Functions file (package libraries, summary functions, plotting)
  source("functions.R")

# ------------------------------------------------------------------------------
# Define output directory and read simulated pharmacodynamic data
# Output is saved into a totally different folder from scripts because my code
# is tracked by git and stored online.  Keeping all output (i.e., .csv and .png)
# files in the git repository slows down synchronisation
# PD output is saved in the same folder as PK output
  output.dir <- paste0("/Volumes/Prosecutor/sunitinib_nogit/",project.name)
# List studies where mg/kg doses were exactly as calculated
  # study.type <- "mgkg_Xkg_exact"
  # study.names <- c("mgkg_40kg_exact","mgkg_50kg_exact","mgkg_60kg_exact",
  #   "mgkg_70kg_exact","mgkg_80kg_exact","mgkg_90kg_exact","mgkg_100kg_exact",
  #   "mgkg_110kg_exact")
# List studies where mg/kg doses were rounded to the nearest 12.5 mg
  # study.type <- "mgkg_Xkg_round"
  # study.names <- c("mgkg_40kg_round","mgkg_50kg_round","mgkg_60kg_round",
  #   "mgkg_70kg_round","mgkg_80kg_round","mgkg_90kg_round","mgkg_100kg_round",
  #   "mgkg_110kg_round")
# List studies that administered the standard 50 mg dose
  study.type <- "standard_Xkg"
  study.names <- c("standard_40kg","standard_50kg","standard_60kg",
    "standard_70kg","standard_80kg","standard_90kg","standard_100kg",
    "standard_110kg")
# Read in PD for each study listed in "study.names"
# Bind all simulation results together in one data frame
  read.study.data <- function(study.names) {
    study.dir <- paste0(output.dir,"/",study.names)	# Name study's directory
    setwd(study.dir)	# Set working directory to study directory
    pd.data <- read.csv(file = paste0(study.names,"_pd_data.csv"))	# Read data
  }
  pd.data <- ldply(study.names, read.study.data, .progress = "text")
  pd.data$WTf <- as.factor(pd.data$WT)	# Make weight a factor for plotting

# ------------------------------------------------------------------------------
# Create a folder to save summary results
  save.name <- "summary_bystudy"
  save.dir <- paste0(output.dir,"/",save.name)
  dir.create(save.dir)	# Create the folder if not done so
  setwd(save.dir)	# Set to working directory

# ------------------------------------------------------------------------------
# Summarise pharmacokinetic relationships
# Subset data to only include day 28
  day28.data <- pd.data[pd.data$time == 28*24,]
# Day 28 trough concentrations versus weight
  plotobj1 <- NULL
  plotobj1 <- ggplot(day28.data)
  plotobj1 <- plotobj1 + geom_boxplot(aes(x = WTf,y = IPRE,colour = WTf,
    fill = WTf),alpha = 0.3)
  plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 50),
    linetype = "dashed",colour = "brown3")
  plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 100),
    linetype = "dashed",colour = "brown3")
  plotobj1 <- plotobj1 + scale_y_continuous(
    "Total Sunitinib Concentration (ng/mL) at Day 28",lim = c(0,225),
    breaks = c(0,25,50,75,100,125,150,175,200,225),
    labels = c(0,25,50,75,100,125,150,175,200,225))
  plotobj1 <- plotobj1 + scale_x_discrete("Total Body Weight (kg)")
  plotobj1 <- plotobj1 + theme(legend.position = "none")
  # print(plotobj1)
  ggsave(plot = plotobj1,filename = paste0(study.type,"_IPREvsWT.png"),
    height = 15,width = 20,units = "cm",dpi = 300)
# Day 28 24-hour AUC versus weight
  plotobj2 <- NULL
  plotobj2 <- ggplot(day28.data)
  plotobj2 <- plotobj2 + geom_boxplot(aes(x = WTf,y = AUC24,colour = WTf,
    fill = WTf),alpha = 0.3)
  plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 1.5),
    linetype = "dashed",colour = "brown3")
  plotobj2 <- plotobj2 + scale_y_continuous(
    "Total Sunitinib 24-hour AUC (mg*h/L) at Day 28",lim = c(0,NA))
  plotobj2 <- plotobj2 + scale_x_discrete("Total Body Weight (kg)")
  plotobj2 <- plotobj2 + theme(legend.position = "none")
  # print(plotobj2)
  ggsave(plot = plotobj2,filename = paste0(study.type,"_AUC24vsWT.png"),
    height = 15,width = 20,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Summarise pharmacodynamic relationships
# Subset data to only include the last cycle of dosing
  max.cyc <- max(pd.data$cyc)
  lastcyc.data <- pd.data[pd.data$cyc == max.cyc,]
# For each individual, find the ANC nadir
  anc.nadir.data <- ddply(lastcyc.data, .(SIM,ID,WTf),
    function(lastcyc.data) int.min.max(lastcyc.data$ANC,lastcyc.data$time))
  # Plot ANC nadir in the last interval versus weight
    plotobj3 <- NULL
    plotobj3 <- ggplot(anc.nadir.data)
    plotobj3 <- plotobj3 + geom_boxplot(aes(x = WTf,y = min.value,colour = WTf,
      fill = WTf),alpha = 0.3)
    plotobj3 <- plotobj3 + geom_hline(aes(yintercept = 0.5),
      linetype = "dashed",colour = "brown3")
    plotobj3 <- plotobj3 + scale_y_log10(
      "Absolute Neutrophil Count (x10^9)\n(Nadir in Last Cycle)",
      breaks = c(0.1,0.5,1,5,10),labels = c(0.1,0.5,1,5,10))
    plotobj3 <- plotobj3 + scale_x_discrete("Total Body Weight (kg)")
    plotobj3 <- plotobj3 + theme(legend.position = "none")
    # print(plotobj3)
    ggsave(plot = plotobj3,filename = paste0(study.type,"_ANCvsWT.png"),
      height = 15,width = 20,units = "cm",dpi = 300)
# For each individual, find the maximum dBP
  max.dbp.data <- ddply(lastcyc.data, .(SIM,ID,WTf),
    function(lastcyc.data) int.min.max(lastcyc.data$BP,lastcyc.data$time))
  # Plot maximum diastolic blood pressure in the last interval versus weight
    plotobj4 <- NULL
    plotobj4 <- ggplot(max.dbp.data)
    plotobj4 <- plotobj4 + geom_boxplot(aes(x = WTf,y = max.value,colour = WTf,
      fill = WTf),alpha = 0.3)
    plotobj4 <- plotobj4 + scale_y_continuous(
      "Diastolic Blood Pressure (mmHg)\n(Maximum in Last Cycle)",
      breaks = seq(60,180,20),labels = seq(60,180,20))
    plotobj4 <- plotobj4 + scale_x_discrete("Total Body Weight (kg)")
    plotobj4 <- plotobj4 + theme(legend.position = "none")
    # print(plotobj4)
    ggsave(plot = plotobj4,filename = paste0(study.type,"_BPvsWT.png"),
      height = 15,width = 20,units = "cm",dpi = 300)
# For each individual, find the minimum sKIT
  min.skit.data <- ddply(lastcyc.data, .(SIM,ID,WTf),
    function(lastcyc.data) int.min.max(lastcyc.data$IPRE_SKIT,
    lastcyc.data$time))
  # Plot minimum sKIT in the last interval versus weight
    plotobj5 <- NULL
    plotobj5 <- ggplot(min.skit.data)
    plotobj5 <- plotobj5 + geom_boxplot(aes(x = WTf,y = min.value,colour = WTf,
      fill = WTf),alpha = 0.3)
    plotobj5 <- plotobj5 + scale_y_continuous(
      "sKIT Concentration (pg/mL)\n(Minimum in Last Cycle)",
      breaks = seq(8,20,2),labels = seq(8,20,2))
    plotobj5 <- plotobj5 + scale_x_discrete("Total Body Weight (kg)")
    plotobj5 <- plotobj5 + theme(legend.position = "none")
    # print(plotobj5)
    ggsave(plot = plotobj5,filename = paste0(study.type,"_sKITvsWT.png"),
      height = 15,width = 20,units = "cm",dpi = 300)
# For each individual, find the minimum s-VEGFR3
  min.vegfr3.data <- ddply(lastcyc.data, .(SIM,ID,WTf),
    function(lastcyc.data) int.min.max(lastcyc.data$IPRE_VEGFR3,
    lastcyc.data$time))
  # Plot minimum s-VEGFR3 in the last interval versus weight
    plotobj6 <- NULL
    plotobj6 <- ggplot(min.vegfr3.data)
    plotobj6 <- plotobj6 + geom_boxplot(aes(x = WTf,y = min.value,colour = WTf,
      fill = WTf),alpha = 0.3)
    plotobj6 <- plotobj6 + scale_y_continuous(
      "s-VEGFR3 Concentration (pg/mL)\n(Minimum in Last Cycle)",
      breaks = seq(8,14,1),labels = seq(8,14,1))
    plotobj6 <- plotobj6 + scale_x_discrete("Total Body Weight (kg)")
    plotobj6 <- plotobj6 + theme(legend.position = "none")
    # print(plotobj6)
    ggsave(plot = plotobj6,filename = paste0(study.type,"_sVEGFR3vsWT.png"),
      height = 15,width = 20,units = "cm",dpi = 300)
# For each individual, find the minimum tumour size and its corresponding time
  min.tumour.data <- ddply(pd.data, .(SIM,ID,WTf),
    function(pd.data) int.min.max(pd.data$TUMOUR,pd.data$time))
  # Plot minimum tumour size versus weight
    plotobj7 <- NULL
    plotobj7 <- ggplot(min.tumour.data)
    plotobj7 <- plotobj7 + geom_boxplot(aes(x = WTf,y = min.value,colour = WTf,
      fill = WTf),alpha = 0.3)
    plotobj7 <- plotobj7 + scale_y_continuous(
      "Tumour Size (Sum of Longest Diameters, mm)",
      lim = c(0,800),breaks = seq(0,800,100),labels = seq(0,800,100))
    plotobj7 <- plotobj7 + scale_x_discrete("Total Body Weight (kg)")
    plotobj7 <- plotobj7 + theme(legend.position = "none")
    # print(plotobj7)
    ggsave(plot = plotobj7,filename = paste0(study.type,"_minTumourSizevsWT.png"),
      height = 15,width = 20,units = "cm",dpi = 300)
  # Plot time at which minimum tumour size occurs versus weight
    plotobj8 <- NULL
    plotobj8 <- ggplot(min.tumour.data)
    plotobj8 <- plotobj8 + geom_boxplot(aes(x = WTf,y = min.time/24/7,
      colour = WTf,fill = WTf),alpha = 0.3)
    plotobj8 <- plotobj8 + scale_y_continuous(
      "Time Minimum Tumour Size Occurs (weeks)",
      lim = c(0,250),breaks = seq(0,250,50),labels = seq(0,250,50))
    plotobj8 <- plotobj8 + scale_x_discrete("Total Body Weight (kg)")
    plotobj8 <- plotobj8 + theme(legend.position = "none")
    # print(plotobj8)
    ggsave(plot = plotobj8,filename = paste0(study.type,"_minTumourTimevsWT.png"),
      height = 15,width = 20,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Subset "pd.data" for very last time-point
  max.time <- max(pd.data$time)
  last.data <- pd.data[pd.data$time == max.time,]
# Plot proportion of individuals in each fatigue grade at study conclusion
  # Make FAT a factor
    last.data$FATf <- as.factor(last.data$FAT)
  # Plot
    plotobj9 <- NULL
    plotobj9 <- ggplot(last.data)
    plotobj9 <- plotobj9 + geom_bar(aes(x = WTf,fill = FATf),stat = "count")
    plotobj9 <- plotobj9 + scale_y_continuous("Number of Individuals")
    plotobj9 <- plotobj9 + scale_x_discrete("Total Body Weight (kg)")
    plotobj9 <- plotobj9 + labs(fill = "Fatigue Grade")
    # print(plotobj9)
    ggsave(plot = plotobj9,filename = paste0(study.type,"_FATvsWT.png"),
      height = 15,width = 20,units = "cm",dpi = 300)
# Plot proportion of individuals in each hand-foot syndrome grade at study end
  # Make HFS a factor
    last.data$HFSf <- as.factor(last.data$HFS)
  # Plot
    plotobj10 <- NULL
    plotobj10 <- ggplot(last.data)
    plotobj10 <- plotobj10 + geom_bar(aes(x = WTf,fill = HFSf),stat = "count")
    plotobj10 <- plotobj10 + scale_y_continuous("Number of Individuals")
    plotobj10 <- plotobj10 + scale_x_discrete("Total Body Weight (kg)")
    plotobj10 <- plotobj10 + labs(fill = "Hand-Foot\nSyndrome Grade")
    # print(plotobj10)
    ggsave(plot = plotobj10,filename = paste0(study.type,"_HFSvsWT.png"),
      height = 15,width = 20,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Overall survival
