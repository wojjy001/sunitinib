# Read in simulated pharmacodynamic data for weight categories
# Summarise simulation results against studies
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
# List studies involving 40 kg individuals
  # study.type <- "40kg"
  # study.names <- c("standard_40kg","mgkg_40kg_exact","mgkg_40kg_round")
# List studies involving 50 kg individuals
  # study.type <- "50kg"
  # study.names <- c("standard_50kg","mgkg_50kg_exact","mgkg_50kg_round")
# List studies involving 60 kg individuals
  # study.type <- "60kg"
  # study.names <- c("standard_60kg","mgkg_60kg_exact","mgkg_60kg_round")
# List studies involving 70 kg individuals
  # study.type <- "70kg"
  # study.names <- c("standard_70kg","mgkg_70kg_exact","mgkg_70kg_round")
# List studies involving 80 kg individuals
  # study.type <- "80kg"
  # study.names <- c("standard_80kg","mgkg_80kg_exact","mgkg_80kg_round")
# List studies involving 90 kg individuals
  # study.type <- "90kg"
  # study.names <- c("standard_90kg","mgkg_90kg_exact","mgkg_90kg_round")
# List studies involving 100 kg individuals
  # study.type <- "100kg"
  # study.names <- c("standard_100kg","mgkg_100kg_exact","mgkg_100kg_round")
# List studies involving 110 kg individuals
  study.type <- "110kg"
  study.names <- c("standard_110kg","mgkg_110kg_exact","mgkg_110kg_round")
# Read in PD for each study listed in "study.names"
# Bind all simulation results together in one data frame
  read.study.data <- function(study.names) {
    study.dir <- paste0(output.dir,"/",study.names)	# Name study's directory
    setwd(study.dir)	# Set working directory to study directory
    pd.data <- read.csv(file = paste0(study.names,"_pd_data.csv"))	# Read data
    pd.data$study <- study.names	# Add labels for the study name
    pd.data
  }
  pd.data <- ldply(study.names, read.study.data, .progress = "text")
  mgkg.ratio <- 50/70
  exact.dose <- pd.data$WT[1]*mgkg.ratio
  round.dose <- round(pd.data$WT[1]/12.5*mgkg.ratio)*12.5
  pd.data$studyf <- as.factor(pd.data$study)	# Make study a factor for plotting
  levels(pd.data$studyf) <- c(
    paste0(round(exact.dose,digits = 2)," mg (",round(mgkg.ratio,digits = 2),
    " mg/kg)"),
    paste0(round.dose," mg (",round(round.dose/pd.data$WT[1],digits = 2),
    " mg/kg)"),
    "50 mg"
    )

# ------------------------------------------------------------------------------
# Create a folder to save summary results
  save.name <- "summary_byweight"
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
  plotobj1 <- plotobj1 + geom_boxplot(aes(x = studyf,y = IPRE,colour = studyf,
    fill = studyf),alpha = 0.3)
  plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 50),
    linetype = "dashed",colour = "brown3")
  plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 100),
    linetype = "dashed",colour = "brown3")
  plotobj1 <- plotobj1 + scale_y_continuous(
    "Total Sunitinib Concentration (ng/mL) at Day 28",lim = c(0,NA),
    breaks = seq(0,300,25),labels = seq(0,300,25))
  plotobj1 <- plotobj1 + scale_x_discrete("Dosing Strategy")
  plotobj1 <- plotobj1 + theme(legend.position = "none")
  # print(plotobj1)
  ggsave(plot = plotobj1,filename = paste0(study.type,"_IPREvsstudy.png"),
    height = 15,width = 20,units = "cm",dpi = 300)
# Day 28 24-hour AUC versus weight
  plotobj2 <- NULL
  plotobj2 <- ggplot(day28.data)
  plotobj2 <- plotobj2 + geom_boxplot(aes(x = studyf,y = AUC24,colour = studyf,
    fill = studyf),alpha = 0.3)
  plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 1.5),
    linetype = "dashed",colour = "brown3")
  plotobj2 <- plotobj2 + scale_y_continuous(
    "Total Sunitinib 24-hour AUC (mg*h/L) at Day 28",lim = c(0,NA))
  plotobj2 <- plotobj2 + scale_x_discrete("Dosing Strategy")
  plotobj2 <- plotobj2 + theme(legend.position = "none")
  # print(plotobj2)
  ggsave(plot = plotobj2,filename = paste0(study.type,"_AUC24vsstudy.png"),
    height = 15,width = 20,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Summarise pharmacodynamic relationships
# Subset data to only include the last cycle of dosing
  max.cyc <- max(pd.data$cyc)
  lastcyc.data <- pd.data[pd.data$cyc == max.cyc,]
# For each individual, find the ANC nadir
  anc.nadir.data <- ddply(lastcyc.data, .(SIM,ID,studyf),
    function(lastcyc.data) int.min.max(lastcyc.data$ANC,lastcyc.data$time))
  # Plot ANC nadir in the last interval versus weight
    plotobj3 <- NULL
    plotobj3 <- ggplot(anc.nadir.data)
    plotobj3 <- plotobj3 + geom_boxplot(aes(x = studyf,y = min.value,
      colour = studyf,fill = studyf),alpha = 0.3)
    plotobj3 <- plotobj3 + geom_hline(aes(yintercept = 0.5),
      linetype = "dashed",colour = "brown3")
    plotobj3 <- plotobj3 + scale_y_log10(
      "Absolute Neutrophil Count (x10^9)\n(Nadir in Last Cycle)",
      breaks = c(0.1,0.5,1,5,10),labels = c(0.1,0.5,1,5,10))
    plotobj3 <- plotobj3 + scale_x_discrete("Dosing Strategy")
    plotobj3 <- plotobj3 + theme(legend.position = "none")
    # print(plotobj3)
    ggsave(plot = plotobj3,filename = paste0(study.type,"_ANCvsstudy.png"),
      height = 15,width = 20,units = "cm",dpi = 300)
# For each individual, find the maximum dBP
  max.dbp.data <- ddply(lastcyc.data, .(SIM,ID,studyf),
    function(lastcyc.data) int.min.max(lastcyc.data$BP,lastcyc.data$time))
  # Plot maximum diastolic blood pressure in the last interval versus weight
    plotobj4 <- NULL
    plotobj4 <- ggplot(max.dbp.data)
    plotobj4 <- plotobj4 + geom_boxplot(aes(x = studyf,y = max.value,
      colour = studyf,fill = studyf),alpha = 0.3)
    plotobj4 <- plotobj4 + scale_y_continuous(
      "Diastolic Blood Pressure (mmHg)\n(Maximum in Last Cycle)",
      breaks = seq(60,260,20),labels = seq(60,260,20))
    plotobj4 <- plotobj4 + scale_x_discrete("Dosing Strategy")
    plotobj4 <- plotobj4 + theme(legend.position = "none")
    # print(plotobj4)
    ggsave(plot = plotobj4,filename = paste0(study.type,"_BPvsstudy.png"),
      height = 15,width = 20,units = "cm",dpi = 300)
# For each individual, find the minimum sKIT
  min.skit.data <- ddply(lastcyc.data, .(SIM,ID,studyf),
    function(lastcyc.data) int.min.max(lastcyc.data$IPRE_SKIT,
    lastcyc.data$time))
  # Plot minimum sKIT in the last interval versus weight
    plotobj5 <- NULL
    plotobj5 <- ggplot(min.skit.data)
    plotobj5 <- plotobj5 + geom_boxplot(aes(x = studyf,y = min.value,
      colour = studyf,fill = studyf),alpha = 0.3)
    plotobj5 <- plotobj5 + scale_y_continuous(
      "sKIT Concentration (pg/mL)\n(Minimum in Last Cycle)",
      breaks = seq(8,20,2),labels = seq(8,20,2))
    plotobj5 <- plotobj5 + scale_x_discrete("Dosing Strategy")
    plotobj5 <- plotobj5 + theme(legend.position = "none")
    # print(plotobj5)
    ggsave(plot = plotobj5,filename = paste0(study.type,"_sKITvsstudy.png"),
      height = 15,width = 20,units = "cm",dpi = 300)
# For each individual, find the minimum s-VEGFR3
  min.vegfr3.data <- ddply(lastcyc.data, .(SIM,ID,studyf),
    function(lastcyc.data) int.min.max(lastcyc.data$IPRE_VEGFR3,
    lastcyc.data$time))
  # Plot minimum s-VEGFR3 in the last interval versus weight
    plotobj6 <- NULL
    plotobj6 <- ggplot(min.vegfr3.data)
    plotobj6 <- plotobj6 + geom_boxplot(aes(x = studyf,y = min.value,
      colour = studyf,fill = studyf),alpha = 0.3)
    plotobj6 <- plotobj6 + scale_y_continuous(
      "s-VEGFR3 Concentration (pg/mL)\n(Minimum in Last Cycle)",
      breaks = seq(6,14,1),labels = seq(6,14,1))
    plotobj6 <- plotobj6 + scale_x_discrete("Dosing Strategy")
    plotobj6 <- plotobj6 + theme(legend.position = "none")
    # print(plotobj6)
    ggsave(plot = plotobj6,filename = paste0(study.type,"_sVEGFR3vsstudy.png"),
      height = 15,width = 20,units = "cm",dpi = 300)
# For each individual, find the minimum tumour size and its corresponding time
  min.tumour.data <- ddply(pd.data, .(SIM,ID,studyf),
    function(pd.data) int.min.max(pd.data$TUMOUR,pd.data$time))
  # Plot minimum tumour size versus weight
    plotobj7 <- NULL
    plotobj7 <- ggplot(min.tumour.data)
    plotobj7 <- plotobj7 + geom_boxplot(aes(x = studyf,y = min.value,
      colour = studyf,fill = studyf),alpha = 0.3)
    plotobj7 <- plotobj7 + scale_y_continuous(
      "Tumour Size (Sum of Longest Diameters, mm)",
      lim = c(0,NA),breaks = seq(0,800,100),labels = seq(0,800,100))
    plotobj7 <- plotobj7 + scale_x_discrete("Dosing Strategy")
    plotobj7 <- plotobj7 + theme(legend.position = "none")
    # print(plotobj7)
    ggsave(plot = plotobj7,
      filename = paste0(study.type,"_minTumourSizevsstudy.png"),
      height = 15,width = 20,units = "cm",dpi = 300)
  # Plot time at which minimum tumour size occurs versus weight
    plotobj8 <- NULL
    plotobj8 <- ggplot(min.tumour.data)
    plotobj8 <- plotobj8 + geom_boxplot(aes(x = studyf,y = min.time/24/7,
      colour = studyf,fill = studyf),alpha = 0.3)
    plotobj8 <- plotobj8 + scale_y_continuous(
      "Time Minimum Tumour Size Occurs (weeks)",
      lim = c(0,250),breaks = seq(0,250,50),labels = seq(0,250,50))
    plotobj8 <- plotobj8 + scale_x_discrete("Dosing Strategy")
    plotobj8 <- plotobj8 + theme(legend.position = "none")
    # print(plotobj8)
    ggsave(plot = plotobj8,
      filename = paste0(study.type,"_minTumourTimevsstudy.png"),
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
    plotobj9 <- plotobj9 + geom_bar(aes(x = studyf,fill = FATf),stat = "count")
    plotobj9 <- plotobj9 + scale_y_continuous("Number of Individuals")
    plotobj9 <- plotobj9 + scale_x_discrete("Dosing Strategy")
    plotobj9 <- plotobj9 + labs(fill = "Fatigue Grade")
    # print(plotobj9)
    ggsave(plot = plotobj9,filename = paste0(study.type,"_FATvsstudy.png"),
      height = 15,width = 20,units = "cm",dpi = 300)
# Plot proportion of individuals in each hand-foot syndrome grade at study end
  # Make HFS a factor
    last.data$HFSf <- as.factor(last.data$HFS)
  # Plot
    plotobj10 <- NULL
    plotobj10 <- ggplot(last.data)
    plotobj10 <- plotobj10 + geom_bar(aes(x = studyf,fill = HFSf),stat = "count")
    plotobj10 <- plotobj10 + scale_y_continuous("Number of Individuals")
    plotobj10 <- plotobj10 + scale_x_discrete("Dosing Strategy")
    plotobj10 <- plotobj10 + labs(fill = "Hand-Foot\nSyndrome Grade")
    # print(plotobj10)
    ggsave(plot = plotobj10,filename = paste0(study.type,"_HFSvsstudy.png"),
      height = 15,width = 20,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Overall survival
  # For overall survival plots, calculate the proportion of individuals alive
  # at each time-point
  # Kaplan Meier Plot for survival confidence intervals
  # All individuals have the same start time, i.e., time == 0
    pd.data$start <- 0
    pd.data <- ddply(pd.data, .(SIM,ID,studyf), stop.time.function) # For each
    # individual calculate their stop time
    km.data <- ddply(pd.data, .(SIM,ID,studyf), headperID)
    km.data <- km.data[c("SIM","ID","studyf","start","stop","event")]
    S <- Surv(time = km.data$start,time2 = km.data$stop,event = km.data$event)
    result <- survfit(formula = S ~ studyf,data = km.data)
    cols <- lapply(2:12, function(x) summary(result)[x])
    surv.data <- do.call(data.frame, cols)
    surv.data$strata <- as.factor(surv.data$strata)
    levels(surv.data$strata) <- levels(pd.data$studyf)
  # Plot overall survival
    plotobj11 <- NULL
    plotobj11 <- ggplot()
    plotobj11 <- plotobj11 + geom_ribbon(aes(x = time/24/7,ymin = lower,
      ymax = upper,fill = strata),data = surv.data,alpha = 0.2)
    plotobj11 <- plotobj11 + geom_line(aes(x = time/24/7,y = surv,
      colour = strata),data = surv.data)
    plotobj11 <- plotobj11 + scale_y_continuous("Probability of Survival",
      lim = c(0,1),
      breaks = seq(from = 0,to = 1,by = 0.2),
      labels = seq(from = 0,to = 1,by = 0.2))
    plotobj11 <- plotobj11 + scale_x_continuous("Time (weeks)",
      breaks = seq(from = 0,to = max(surv.data$time)/24/7,by = 30),
      lim = c(0,max(surv.data$time/24/7)))
    plotobj11 <- plotobj11 + labs(colour = "Dosing Strategy",
      fill = "Dosing Strategy")
    # print(plotobj11)
    ggsave(plot = plotobj11,
      filename = paste0(study.type,"_overallsurvival.png"),
      width = 20,height = 15,unit = "cm",dpi = 300)
