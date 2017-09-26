# Read in simulated pharmacodynamic data for each target AUC
# Summarise AUC24 versus PD outcome relationships
# ------------------------------------------------------------------------------
# Remove objects in workspace
# Define simulation file directory for the project
  file.dir <- "/Volumes/Prosecutor/sunitinib/Project/SimulationFiles/"
  project.name <- "ExposureResponse"
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
# List studies
  study.names <- c("target_1.0","target_1.2","target_1.4","target_1.5",
    "target_1.6","target_1.8","target_2.0","target_3.0")
# Bind all simulation results together in one data frame
  read.study.data <- function(study.names) {
    study.dir <- paste0(output.dir,"/",study.names)	# Name study's directory
    setwd(study.dir)	# Set working directory to study directory
    pd.data <- read.csv(file = paste0(study.names,"_pd_data.csv"))	# Read data
    pd.data$study <- study.names	# Add labels for the study name
    pd.data
  }
  pd.data <- ldply(study.names, read.study.data, .progress = "text")
# Make "study" a factor for plotting
  pd.data$study <- as.factor(pd.data$study)
# Round AUC24 to 2 significant figures
  pd.data$AUC24r <- round(pd.data$AUC24,digits = 1)
  pd.data$AUC24r[pd.data$AUC24r > 3] <- 3

# ------------------------------------------------------------------------------
# Create a folder to save summary results
  save.name <- "summary"
  save.dir <- paste0(output.dir,"/",save.name)
  dir.create(save.dir)	# Create the folder if not done so
  setwd(save.dir)	# Set to working directory

# ------------------------------------------------------------------------------
# Summarise outcomes at week 50
# Tumour model, for example, doesn't extrapolate well beyond this time-point
  pd52.data <- pd.data[pd.data$time == 52*7*24,]

# ------------------------------------------------------------------------------
# Population Summaries
# Tumour versus 24-hour AUC
  plotobj1 <- NULL
  plotobj1 <- ggplot(pd52.data)
  plotobj1 <- plotobj1 + geom_boxplot(aes(x = AUC24r,y = TUMOUR,
    colour = study,fill = study),alpha = 0.3)
  plotobj1 <- plotobj1 + scale_y_continuous(
    "Tumour (Sum of Longest Diameters, mm)",
    breaks = seq(0,5000,100),labels = seq(0,5000,100))
  plotobj1 <- plotobj1 + scale_x_continuous(
    "Total Sunitinib 24-hour AUC (mg*h/L)",
    breaks = c(1,1.2,1.4,1.5,1.6,1.8,2,3),
    labels = c(1,1.2,1.4,1.5,1.6,1.8,2,3))
  plotobj1 <- plotobj1 + theme(legend.position = "none")
  print(plotobj1)
  ggsave(plot = plotobj1,filename = paste0("TUMOURvsAUC24.png"),
    height = 15,width = 20,units = "cm",dpi = 300)
# ANC versus 24-hour AUC
  plotobj2 <- NULL
  plotobj2 <- ggplot(pd52.data)
  plotobj2 <- plotobj2 + geom_boxplot(aes(x = AUC24r,y = ANC,
    colour = study,fill = study),alpha = 0.3)
  plotobj2 <- plotobj2 + scale_y_log10("ANC (x10^9)",
    breaks = c(0.1,0.3,1,3,10),labels = c(0.1,0.3,1,3,10))
  plotobj2 <- plotobj2 + scale_x_continuous(
    "Total Sunitinib 24-hour AUC (mg*h/L)",
    breaks = c(1,1.2,1.4,1.5,1.6,1.8,2,3),
    labels = c(1,1.2,1.4,1.5,1.6,1.8,2,3))
  plotobj2 <- plotobj2 + theme(legend.position = "none")
  print(plotobj2)
  ggsave(plot = plotobj2,filename = paste0("ANCvsAUC24.png"),
    height = 15,width = 20,units = "cm",dpi = 300)
# dBP versus 24-hour AUC
  plotobj3 <- NULL
  plotobj3 <- ggplot(pd52.data)
  plotobj3 <- plotobj3 + geom_boxplot(aes(x = AUC24r,y = BP,
    colour = study,fill = study),alpha = 0.3)
  plotobj3 <- plotobj3 + scale_y_continuous("dBP (mmHg)",
    breaks = seq(40,200,20),labels = seq(40,200,20))
  plotobj3 <- plotobj3 + scale_x_continuous(
    "Total Sunitinib 24-hour AUC (mg*h/L)",
    breaks = c(1,1.2,1.4,1.5,1.6,1.8,2,3),
    labels = c(1,1.2,1.4,1.5,1.6,1.8,2,3))
  plotobj3 <- plotobj3 + theme(legend.position = "none")
  print(plotobj3)
  ggsave(plot = plotobj3,filename = paste0("BPvsAUC24.png"),
    height = 15,width = 20,units = "cm",dpi = 300)
# sVEGFR-3 versus 24-hour AUC
  plotobj4 <- NULL
  plotobj4 <- ggplot(pd52.data)
  plotobj4 <- plotobj4 + geom_boxplot(aes(x = AUC24r,y = IPRE_VEGFR3,
    colour = study,fill = study),alpha = 0.3)
  plotobj4 <- plotobj4 + scale_y_continuous("sVEGFR-3 (pg/mL)")
  plotobj4 <- plotobj4 + scale_x_continuous(
    "Total Sunitinib 24-hour AUC (mg*h/L)",
    breaks = c(1,1.2,1.4,1.5,1.6,1.8,2,3),
    labels = c(1,1.2,1.4,1.5,1.6,1.8,2,3))
  plotobj4 <- plotobj4 + theme(legend.position = "none")
  print(plotobj4)
  ggsave(plot = plotobj4,filename = paste0("VEGFR3vsAUC24.png"),
    height = 15,width = 20,units = "cm",dpi = 300)
# sKIT versus 24-hour AUC
  plotobj5 <- NULL
  plotobj5 <- ggplot(pd52.data)
  plotobj5 <- plotobj5 + geom_boxplot(aes(x = AUC24r,y = IPRE_SKIT,
    colour = study,fill = study),alpha = 0.3)
  plotobj5 <- plotobj5 + scale_y_continuous("sKIT (pg/mL)")
  plotobj5 <- plotobj5 + scale_x_continuous(
    "Total Sunitinib 24-hour AUC (mg*h/L)",
    breaks = c(1,1.2,1.4,1.5,1.6,1.8,2,3),
    labels = c(1,1.2,1.4,1.5,1.6,1.8,2,3))
  plotobj5 <- plotobj5 + theme(legend.position = "none")
  print(plotobj5)
  ggsave(plot = plotobj5,filename = paste0("SKITvsAUC24.png"),
    height = 15,width = 20,units = "cm",dpi = 300)
# Hand-Foot Syndrome Grade versus 24-hour AUC
  pd52.data$HFSf <- as.factor(pd52.data$HFS)
  plotobj6 <- NULL
  plotobj6 <- ggplot(pd52.data)
  plotobj6 <- plotobj6 + geom_bar(aes(x = AUC24r,fill = HFSf),stat = "count")
  plotobj6 <- plotobj6 + scale_y_continuous("Number of Individuals")
  plotobj6 <- plotobj6 + scale_x_continuous(
    "Total Sunitinib 24-hour AUC (mg*h/L)",
    breaks = c(1,1.2,1.4,1.5,1.6,1.8,2,3),
    labels = c(1,1.2,1.4,1.5,1.6,1.8,2,3))
  plotobj6 <- plotobj6 + labs(fill = "Hand-Foot\nSyndrome\nGrade")
  print(plotobj6)
  ggsave(plot = plotobj6,filename = paste0("HFSvsAUC24.png"),
    height = 15,width = 20,units = "cm",dpi = 300)
# Fatigue Grade versus 24-hour AUC
  pd52.data$FATf <- as.factor(pd52.data$FAT)
  plotobj7 <- NULL
  plotobj7 <- ggplot(pd52.data)
  plotobj7 <- plotobj7 + geom_bar(aes(x = AUC24r,fill = FATf),stat = "count")
  plotobj7 <- plotobj7 + scale_y_continuous("Number of Individuals")
  plotobj7 <- plotobj7 + scale_x_continuous(
    "Total Sunitinib 24-hour AUC (mg*h/L)",
    breaks = c(1,1.2,1.4,1.5,1.6,1.8,2,3),
    labels = c(1,1.2,1.4,1.5,1.6,1.8,2,3))
  plotobj7 <- plotobj7 + labs(fill = "Fatigue\nGrade")
  print(plotobj7)
  ggsave(plot = plotobj7,filename = paste0("FATvsAUC24.png"),
    height = 15,width = 20,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Overall survival
# For overall survival plots, calculate the proportion of individuals alive
# at each time-point
# Kaplan Meier Plot for survival confidence intervals
# All individuals have the same start time, i.e., time == 0
  pd.data$start <- 0
  pd.data <- ddply(pd.data, .(SIM,ID,study), stop.time.function) # For each
  # individual calculate their stop time
  km.data <- ddply(pd.data, .(SIM,ID,study), headperID)
  km.data <- km.data[c("SIM","ID","study","start","stop","event")]
  S <- Surv(time = km.data$start,time2 = km.data$stop,event = km.data$event)
  result <- survfit(formula = S ~ study,data = km.data)
  cols <- lapply(2:12, function(x) summary(result)[x])
  surv.data <- do.call(data.frame, cols)
  surv.data$strata <- as.factor(surv.data$strata)
  levels(surv.data$strata) <- levels(as.factor(pd52.data$AUC24r))
# Plot overall survival
  plotobj8 <- NULL
  plotobj8 <- ggplot()
  plotobj8 <- plotobj8 + geom_ribbon(aes(x = time/24/7,ymin = lower,
    ymax = upper,fill = strata),data = surv.data,alpha = 0.2)
  plotobj8 <- plotobj8 + geom_line(aes(x = time/24/7,y = surv,
    colour = strata),data = surv.data)
  plotobj8 <- plotobj8 + scale_y_continuous("Probability of Survival",
    lim = c(0,1),
    breaks = seq(from = 0,to = 1,by = 0.2),
    labels = seq(from = 0,to = 1,by = 0.2))
  plotobj8 <- plotobj8 + scale_x_continuous("Time (weeks)",
    breaks = seq(from = 0,to = max(surv.data$time)/24/7,by = 30),
    lim = c(0,max(surv.data$time/24/7)))
  plotobj8 <- plotobj8 + labs(fill = "Target 24-hour\nAUC (mg*h/L)",
    colour = "Target 24-hour\nAUC (mg*h/L)")
  print(plotobj8)
  ggsave(plot = plotobj8,filename = paste0("OSvsAUC24.png"),
    height = 15,width = 20,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Individual Summaries
# Find 12 random individuals
  rand.ind <- sample(unique(pd.data$ID),12)
  rand.data <- pd52.data[pd52.data$ID %in% rand.ind,]

# Tumour versus 24-hour AUC
  plotobj9 <- NULL
  plotobj9 <- ggplot(rand.data)
  plotobj9 <- plotobj9 + geom_smooth(aes(x = AUC24r,y = TUMOUR),
    method = "lm",se = F,colour = "brown3")
  plotobj9 <- plotobj9 + geom_point(aes(x = AUC24r,y = TUMOUR),
    colour = "skyblue4",size = 2)
  plotobj9 <- plotobj9 + scale_y_continuous(
    "Tumour (Sum of Longest Diameters, mm)",
    breaks = seq(0,5000,100),labels = seq(0,5000,100))
  plotobj9 <- plotobj9 + scale_x_continuous(
    "Total Sunitinib 24-hour AUC (mg*h/L)",
    breaks = c(1,1.5,2,3),
    labels = c(1,1.5,2,3))
  plotobj9 <- plotobj9 + facet_wrap(~ID)
  print(plotobj9)
  ggsave(plot = plotobj9,filename = paste0("12rand_TUMOURvsAUC24.png"),
    height = 15,width = 20,units = "cm",dpi = 300)
# ANC versus 24-hour AUC
  plotobj10 <- NULL
  plotobj10 <- ggplot(rand.data)
  plotobj10 <- plotobj10 + geom_smooth(aes(x = AUC24r,y = ANC),
    method = "lm",se = F,colour = "brown3")
  plotobj10 <- plotobj10 + geom_point(aes(x = AUC24r,y = ANC),
    colour = "skyblue4",size = 2)
  plotobj10 <- plotobj10 + scale_y_continuous(
    "ANC (x10^9)",
    breaks = seq(1,10,1),labels = seq(1,10,1))
  plotobj10 <- plotobj10 + scale_x_continuous(
    "Total Sunitinib 24-hour AUC (mg*h/L)",
    breaks = c(1,1.5,2,3),
    labels = c(1,1.5,2,3))
  plotobj10 <- plotobj10 + facet_wrap(~ID)
  print(plotobj10)
  ggsave(plot = plotobj10,filename = paste0("12rand_ANCvsAUC24.png"),
    height = 15,width = 20,units = "cm",dpi = 300)
# dBP versus 24-hour AUC
  plotobj11 <- NULL
  plotobj11 <- ggplot(rand.data)
  plotobj11 <- plotobj11 + geom_smooth(aes(x = AUC24r,y = BP),
    method = "lm",se = F,colour = "brown3")
  plotobj11 <- plotobj11 + geom_point(aes(x = AUC24r,y = BP),
    colour = "skyblue4",size = 2)
  plotobj11 <- plotobj11 + scale_y_continuous("dBP (mmHg)",
    breaks = seq(40,200,10),labels = seq(40,200,10))
  plotobj11 <- plotobj11 + scale_x_continuous(
    "Total Sunitinib 24-hour AUC (mg*h/L)",
    breaks = c(1,1.5,2,3),
    labels = c(1,1.5,2,3))
  plotobj11 <- plotobj11 + facet_wrap(~ID)
  print(plotobj11)
  ggsave(plot = plotobj11,filename = paste0("12rand_BPvsAUC24.png"),
    height = 15,width = 20,units = "cm",dpi = 300)
# sVEGFR-3 versus 24-hour AUC
  plotobj12 <- NULL
  plotobj12 <- ggplot(rand.data)
  plotobj12 <- plotobj12 + geom_smooth(aes(x = AUC24r,y = IPRE_VEGFR3),
    method = "lm",se = F,colour = "brown3")
  plotobj12 <- plotobj12 + geom_point(aes(x = AUC24r,y = IPRE_VEGFR3),
    colour = "skyblue4",size = 2)
  plotobj12 <- plotobj12 + scale_y_continuous("sVEGFR-3 (pg/mL)")
  plotobj12 <- plotobj12 + scale_x_continuous(
    "Total Sunitinib 24-hour AUC (mg*h/L)",
    breaks = c(1,1.5,2,3),
    labels = c(1,1.5,2,3))
  plotobj12 <- plotobj12 + facet_wrap(~ID)
  print(plotobj12)
  ggsave(plot = plotobj12,filename = paste0("12rand_VEGFR3vsAUC24.png"),
    height = 15,width = 20,units = "cm",dpi = 300)
# sKIT versus 24-hour AUC
  plotobj13 <- NULL
  plotobj13 <- ggplot(rand.data)
  plotobj13 <- plotobj13 + geom_smooth(aes(x = AUC24r,y = IPRE_SKIT),
    method = "lm",se = F,colour = "brown3")
  plotobj13 <- plotobj13 + geom_point(aes(x = AUC24r,y = IPRE_SKIT),
    colour = "skyblue4",size = 2)
  plotobj13 <- plotobj13 + scale_y_continuous("sKIT (pg/mL)")
  plotobj13 <- plotobj13 + scale_x_continuous(
    "Total Sunitinib 24-hour AUC (mg*h/L)",
    breaks = c(1,1.5,2,3),
    labels = c(1,1.5,2,3))
  plotobj13 <- plotobj13 + facet_wrap(~ID)
  print(plotobj13)
  ggsave(plot = plotobj13,filename = paste0("12rand_SKITvsAUC24.png"),
    height = 15,width = 20,units = "cm",dpi = 300)
# Hand-foot syndrome grade versus 24-hour AUC
  plotobj16 <- NULL
  plotobj16 <- ggplot(rand.data)
  plotobj16 <- plotobj16 + geom_step(aes(x = AUC24r,y = HFS),
    colour = "skyblue4",size = 2)
  plotobj16 <- plotobj16 + scale_y_continuous("Hand-Foot Syndrome Grade",
    lim = c(0,4),breaks = c(0,1,2,3,4),labels = c(0,1,2,3,4))
  plotobj16 <- plotobj16 + scale_x_continuous(
    "Total Sunitinib 24-hour AUC (mg*h/L)",
    breaks = c(1,1.5,2,3),
    labels = c(1,1.5,2,3))
  plotobj16 <- plotobj16 + facet_wrap(~ID)
  print(plotobj16)
  ggsave(plot = plotobj16,filename = paste0("12rand_HFSvsAUC24.png"),
    height = 15,width = 20,units = "cm",dpi = 300)
# Fatigue grade versus 24-hour AUC
  plotobj15 <- NULL
  plotobj15 <- ggplot(rand.data)
  plotobj15 <- plotobj15 + geom_step(aes(x = AUC24r,y = FAT),
    colour = "skyblue4",size = 2)
  plotobj15 <- plotobj15 + scale_y_continuous("Fatigue Grade",
    lim = c(0,4),breaks = c(0,1,2,3,4),labels = c(0,1,2,3,4))
  plotobj15 <- plotobj15 + scale_x_continuous(
    "Total Sunitinib 24-hour AUC (mg*h/L)",
    breaks = c(1,1.5,2,3),
    labels = c(1,1.5,2,3))
  plotobj15 <- plotobj15 + facet_wrap(~ID)
  print(plotobj15)
  ggsave(plot = plotobj15,filename = paste0("12rand_FATvsAUC24.png"),
    height = 15,width = 20,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Increasing dose for efficacy
  auc1 <- 1	# Original AUC
  auc2 <- 1.5	# Target AUC
  subset.data <- pd52.data[pd52.data$AUC24r == auc1 | pd52.data$AUC24r == auc2,]
  pd.ratio <- function(subset.data) {
    TUMOUR.ratio <- subset.data$TUMOUR[subset.data$AUC24r == auc2]/
      subset.data$TUMOUR[subset.data$AUC24r == auc1]
    TUMOUR.ratio <- round(TUMOUR.ratio/0.005)*0.005-1
    ANC.ratio <- subset.data$ANC[subset.data$AUC24r == auc2]/
      subset.data$ANC[subset.data$AUC24r == auc1]
    ANC.ratio <- round(ANC.ratio/0.005)*0.005-1
    BP.ratio <- subset.data$BP[subset.data$AUC24r == auc2]/
      subset.data$BP[subset.data$AUC24r == auc1]
    BP.ratio <- round(BP.ratio/0.005)*0.005-1
    IPRE_VEGFR3.ratio <- subset.data$IPRE_VEGFR3[subset.data$AUC24r == auc2]/
      subset.data$IPRE_VEGFR3[subset.data$AUC24r == auc1]
    IPRE_VEGFR3.ratio <- round(IPRE_VEGFR3.ratio/0.005)*0.005-1
    IPRE_SKIT.ratio <- subset.data$IPRE_SKIT[subset.data$AUC24r == auc2]/
      subset.data$IPRE_SKIT[subset.data$AUC24r == auc1]
    IPRE_SKIT.ratio <- round(IPRE_SKIT.ratio/0.005)*0.005-1
    result <- data.frame(TUMOUR.ratio,ANC.ratio,BP.ratio,IPRE_VEGFR3.ratio,
      IPRE_SKIT.ratio)
  }
  ratio.data <- ddply(subset.data, .(SIM,ID), pd.ratio)

# Count the number of individuals who had a 20% reduction in pd outcomes
  sig <- seq(-1,1,0.005)
  success.pro <- function(sig) {
    tumour.success <- length(ratio.data$TUMOUR.ratio[ratio.data$TUMOUR.ratio <=
      sig])/1000
    anc.success <- length(ratio.data$ANC.ratio[ratio.data$ANC.ratio <=
      sig])/1000
    bp.success <- length(ratio.data$BP.ratio[ratio.data$BP.ratio >= sig])/1000
    vegfr3.success <-
      length(ratio.data$IPRE_VEGFR3.ratio[ratio.data$IPRE_VEGFR3.ratio <=
      sig])/1000
    skit.success <-
      length(ratio.data$IPRE_SKIT.ratio[ratio.data$IPRE_SKIT.ratio <=
      sig])/1000
    result <- data.frame(sig,tumour.success,anc.success,bp.success,
      vegfr3.success,skit.success)
  }
  success.data <- ldply(sig, success.pro)
  # success.title <- paste0(
    # "50% Increased Exposure:\nTotal Sunitinib 24-hour AUC of ",
    # auc1," mg*h/L (reference) versus ",auc2," mg*h/L (target)")
  success.title <- "(a)"

  plotobj16 <- NULL
  plotobj16 <- ggplot()
  plotobj16 <- plotobj16 + ggtitle(success.title)
  plotobj16 <- plotobj16 + geom_line(aes(x = sig*100,
    y = tumour.success*100,colour = "Tumour"),data = success.data[success.data$sig <= 0,])
  plotobj16 <- plotobj16 + geom_line(aes(x = sig*100,y = anc.success*100,
    colour = "ANC"),
    data = success.data[success.data$sig <= 0,])
  plotobj16 <- plotobj16 + geom_line(aes(x = sig*100,y = bp.success*100,
    colour = "dBP"),
    data = success.data[success.data$sig >= 0,])
  plotobj16 <- plotobj16 + geom_line(aes(x = sig*100,y = vegfr3.success*100,
    colour = "sVEGFR-3"),
    data = success.data[success.data$sig <= 0,])
  plotobj16 <- plotobj16 + geom_line(aes(x = sig*100,y = skit.success*100,
    colour = "sKIT"),data = success.data[success.data$sig <= 0,])
  plotobj16 <- plotobj16 + geom_vline(aes(xintercept = 0),linetype = "dashed")
  plotobj16 <- plotobj16 + scale_y_continuous(
    "Percentage of Population",
    breaks = seq(0,100,10),labels = seq(0,100,10))
  plotobj16 <- plotobj16 + scale_x_continuous(
    "Relative Difference between Exposures (%)",
    breaks = seq(-100,100,20),
    labels = c(-100,-80,-60,-40,-20,0,"+20","+40","+60","+80","+100"))
  plotobj16 <- plotobj16 + labs(colour = "")
  plotobj16 <- plotobj16 + theme(legend.position = "bottom")
  plotobj16 <- plotobj16 + guides(col = guide_legend(ncol = 2,byrow = TRUE))
  print(plotobj16)
  # ggsave(plot = plotobj16,
  #   filename = paste0("50perIncrease_ExposureResponse.png"),
  #   height = 15,width = 25,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Decreasing dose for toxicity
  auc1 <- 2	# Original AUC
  auc2 <- 1.5	# Target AUC
  subset.data <- pd52.data[pd52.data$AUC24r == auc1 | pd52.data$AUC24r == auc2,]
  pd.ratio <- function(subset.data) {
    TUMOUR.ratio <- subset.data$TUMOUR[subset.data$AUC24r == auc2]/
      subset.data$TUMOUR[subset.data$AUC24r == auc1]
    TUMOUR.ratio <- round(TUMOUR.ratio/0.005)*0.005-1
    ANC.ratio <- subset.data$ANC[subset.data$AUC24r == auc2]/
      subset.data$ANC[subset.data$AUC24r == auc1]
    ANC.ratio <- round(ANC.ratio/0.005)*0.005-1
    BP.ratio <- subset.data$BP[subset.data$AUC24r == auc2]/
      subset.data$BP[subset.data$AUC24r == auc1]
    BP.ratio <- round(BP.ratio/0.005)*0.005-1
    IPRE_VEGFR3.ratio <- subset.data$IPRE_VEGFR3[subset.data$AUC24r == auc2]/
      subset.data$IPRE_VEGFR3[subset.data$AUC24r == auc1]
    IPRE_VEGFR3.ratio <- round(IPRE_VEGFR3.ratio/0.005)*0.005-1
    IPRE_SKIT.ratio <- subset.data$IPRE_SKIT[subset.data$AUC24r == auc2]/
      subset.data$IPRE_SKIT[subset.data$AUC24r == auc1]
    IPRE_SKIT.ratio <- round(IPRE_SKIT.ratio/0.005)*0.005-1
    result <- data.frame(TUMOUR.ratio,ANC.ratio,BP.ratio,IPRE_VEGFR3.ratio,
      IPRE_SKIT.ratio)
  }
  ratio.data <- ddply(subset.data, .(SIM,ID), pd.ratio)

# Count the number of individuals who had a 20% reduction in pd outcomes
  sig <- seq(-1,1,0.005)
  success.pro <- function(sig) {
    tumour.success <- length(ratio.data$TUMOUR.ratio[ratio.data$TUMOUR.ratio >=
      sig])/1000
    anc.success <- length(ratio.data$ANC.ratio[ratio.data$ANC.ratio >=
      sig])/1000
    bp.success <- length(ratio.data$BP.ratio[ratio.data$BP.ratio <= sig])/1000
    vegfr3.success <-
      length(ratio.data$IPRE_VEGFR3.ratio[ratio.data$IPRE_VEGFR3.ratio >=
      sig])/1000
    skit.success <-
      length(ratio.data$IPRE_SKIT.ratio[ratio.data$IPRE_SKIT.ratio >=
      sig])/1000
    result <- data.frame(sig,tumour.success,anc.success,bp.success,
      vegfr3.success,skit.success)
  }
  success.data <- ldply(sig, success.pro)
  # success.title <- paste0(
  # 	"25% Decreased Exposure:\nTotal Sunitinib 24-hour AUC of ",
  # 	auc1," mg*h/L (reference) versus ",auc2," mg*h/L (target)")
  success.title <- "(c)"

  plotobj17 <- NULL
  plotobj17 <- ggplot()
  plotobj17 <- plotobj17 + ggtitle(success.title)
  plotobj17 <- plotobj17 + geom_line(aes(x = sig*100,
    y = tumour.success*100,colour = "Tumour"),data = success.data[success.data$sig >= 0,])
  plotobj17 <- plotobj17 + geom_line(aes(x = sig*100,y = anc.success*100,
    colour = "ANC"),
    data = success.data[success.data$sig >= 0,])
  plotobj17 <- plotobj17 + geom_line(aes(x = sig*100,y = bp.success*100,
    colour = "dBP"),
    data = success.data[success.data$sig <= 0,])
  plotobj17 <- plotobj17 + geom_line(aes(x = sig*100,y = vegfr3.success*100,
    colour = "sVEGFR-3"),
    data = success.data[success.data$sig >= 0,])
  plotobj17 <- plotobj17 + geom_line(aes(x = sig*100,y = skit.success*100,
    colour = "sKIT"),data = success.data[success.data$sig >= 0,])
  plotobj17 <- plotobj17 + geom_vline(aes(xintercept = 0),linetype = "dashed")
  plotobj17 <- plotobj17 + scale_y_continuous(
    "Percentage of Population",
    breaks = seq(0,100,10),labels = seq(0,100,10))
  plotobj17 <- plotobj17 + scale_x_continuous(
    "Relative Difference between Exposures (%)",
    breaks = seq(-100,100,20),
    labels = c(-100,-80,-60,-40,-20,0,"+20","+40","+60","+80","+100"))
  plotobj17 <- plotobj17 + labs(colour = "")
  plotobj17 <- plotobj17 + theme(legend.position = "bottom")
  plotobj17 <- plotobj17 + guides(col = guide_legend(ncol = 2,byrow = TRUE))
  print(plotobj17)
  # ggsave(plot = plotobj17,
  #   filename = paste0("25perDecrease_ExposureResponse.png"),
  #   height = 15,width = 25,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Decreasing dose for toxicity
  auc1 <- 3	# Original AUC
  auc2 <- 1.5	# Target AUC
  subset.data <- pd52.data[pd52.data$AUC24r == auc1 | pd52.data$AUC24r == auc2,]
  pd.ratio <- function(subset.data) {
    TUMOUR.ratio <- subset.data$TUMOUR[subset.data$AUC24r == auc2]/
      subset.data$TUMOUR[subset.data$AUC24r == auc1]
    TUMOUR.ratio <- round(TUMOUR.ratio/0.005)*0.005-1
    ANC.ratio <- subset.data$ANC[subset.data$AUC24r == auc2]/
      subset.data$ANC[subset.data$AUC24r == auc1]
    ANC.ratio <- round(ANC.ratio/0.005)*0.005-1
    BP.ratio <- subset.data$BP[subset.data$AUC24r == auc2]/
      subset.data$BP[subset.data$AUC24r == auc1]
    BP.ratio <- round(BP.ratio/0.005)*0.005-1
    IPRE_VEGFR3.ratio <- subset.data$IPRE_VEGFR3[subset.data$AUC24r == auc2]/
      subset.data$IPRE_VEGFR3[subset.data$AUC24r == auc1]
    IPRE_VEGFR3.ratio <- round(IPRE_VEGFR3.ratio/0.005)*0.005-1
    IPRE_SKIT.ratio <- subset.data$IPRE_SKIT[subset.data$AUC24r == auc2]/
      subset.data$IPRE_SKIT[subset.data$AUC24r == auc1]
    IPRE_SKIT.ratio <- round(IPRE_SKIT.ratio/0.005)*0.005-1
    result <- data.frame(TUMOUR.ratio,ANC.ratio,BP.ratio,IPRE_VEGFR3.ratio,
      IPRE_SKIT.ratio)
  }
  ratio.data <- ddply(subset.data, .(SIM,ID), pd.ratio)

# Count the number of individuals who had a 20% reduction in pd outcomes
  sig <- seq(-1,1,0.005)
  success.pro <- function(sig) {
    tumour.success <- length(ratio.data$TUMOUR.ratio[ratio.data$TUMOUR.ratio >=
      sig])/1000
    anc.success <- length(ratio.data$ANC.ratio[ratio.data$ANC.ratio >=
      sig])/1000
    bp.success <- length(ratio.data$BP.ratio[ratio.data$BP.ratio <= sig])/1000
    vegfr3.success <-
      length(ratio.data$IPRE_VEGFR3.ratio[ratio.data$IPRE_VEGFR3.ratio >=
      sig])/1000
    skit.success <-
      length(ratio.data$IPRE_SKIT.ratio[ratio.data$IPRE_SKIT.ratio >=
      sig])/1000
    result <- data.frame(sig,tumour.success,anc.success,bp.success,
      vegfr3.success,skit.success)
  }
  success.data <- ldply(sig, success.pro)
  # success.title <- paste0(
  # 	"50% Decreased Exposure:\nTotal Sunitinib 24-hour AUC of ",
  # 	auc1," mg*h/L (reference) versus ",auc2," mg*h/L (target)")
  success.title <- "(e)"

  plotobj18 <- NULL
  plotobj18 <- ggplot()
  plotobj18 <- plotobj18 + ggtitle(success.title)
  plotobj18 <- plotobj18 + geom_line(aes(x = sig*100,
    y = tumour.success*100,colour = "Tumour"),data = success.data[success.data$sig >= 0,])
  plotobj18 <- plotobj18 + geom_line(aes(x = sig*100,y = anc.success*100,
    colour = "ANC"),
    data = success.data[success.data$sig >= 0,])
  plotobj18 <- plotobj18 + geom_line(aes(x = sig*100,y = bp.success*100,
    colour = "dBP"),
    data = success.data[success.data$sig <= 0,])
  plotobj18 <- plotobj18 + geom_line(aes(x = sig*100,y = vegfr3.success*100,
    colour = "sVEGFR-3"),
    data = success.data[success.data$sig >= 0,])
  plotobj18 <- plotobj18 + geom_line(aes(x = sig*100,y = skit.success*100,
    colour = "sKIT"),data = success.data[success.data$sig >= 0,])
  plotobj18 <- plotobj18 + geom_vline(aes(xintercept = 0),linetype = "dashed")
  plotobj18 <- plotobj18 + scale_y_continuous(
    "Percentage of Population",
    breaks = seq(0,100,10),labels = seq(0,100,10))
  plotobj18 <- plotobj18 + scale_x_continuous(
    "Relative Difference between Exposures (%)",
    breaks = seq(-100,100,20),
    labels = c(-100,-80,-60,-40,-20,0,"+20","+40","+60","+80","+100"))
  plotobj18 <- plotobj18 + labs(colour = "")
  plotobj18 <- plotobj18 + theme(legend.position = "bottom",
    legend.key.size = unit(0.2,units = "cm"))
  plotobj18 <- plotobj18 + guides(col = guide_legend(ncol = 5,byrow = TRUE))
  print(plotobj18)
  # ggsave(plot = plotobj18,
  #   filename = paste0("50perDecrease_ExposureResponse.png"),
  #   height = 15,width = 25,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Combine plots together
  plotobj19 <- plotobj16 + theme(legend.position = "none")
  plotobj20 <- plotobj17 + theme(legend.position = "none")
  plotobj21 <- plotobj18 + guides(col = guide_legend(ncol = 5,byrow = TRUE))
  plotobj22 <- grid.arrange(grobs = list(plotobj19,plotobj20,plotobj21),
    heights = c(5,5,6))
  print(plotobj22)
  # ggsave(plot = plotobj22,filename = paste0("ContExposureResponse.png"),
  #   width = 11.1,height = 12.6,units = "in",dpi = 300)

# ------------------------------------------------------------------------------
# Calculate the proportion of individuals who change Hand-foot syndrome and
# fatigue grades with a different level of sunitinib exposure
  auc1 <- 1	# Original AUC
  auc2 <- 1.5	# Target AUC
  subset.data <- pd52.data[pd52.data$AUC24r == auc1 | pd52.data$AUC24r == auc2,]
  cat.change <- function(subset.data) {
    HFS.change <- subset.data$HFS[subset.data$AUC24r == auc2]-
      subset.data$HFS[subset.data$AUC24r == auc1]
    FAT.change <- subset.data$FAT[subset.data$AUC24r == auc2]-
      subset.data$FAT[subset.data$AUC24r == auc1]
    result <- data.frame(HFS.change,FAT.change)
  }
  change.data <- ddply(subset.data, .(SIM,ID), cat.change)
  # change.title <- paste0(
  # 	"50% Increased Exposure:\nTotal Sunitinib 24-hour AUC of ",
  # 	auc1," mg*h/L (reference) versus ",auc2," mg*h/L (target)")
  change.title <- "(b)"

  plotobj23 <- NULL
  plotobj23 <- ggplot(change.data)
  plotobj23 <- plotobj23 + ggtitle(change.title)
  plotobj23 <- plotobj23 + geom_histogram(aes(x = HFS.change,
    y = ..count../1000*100,fill = "Hand-Foot Syndrome",
    colour = "Hand-Foot Syndrome"),alpha = 0.5,binwidth = 1)
  plotobj23 <- plotobj23 + geom_histogram(aes(x = FAT.change,
    y = ..count../1000*100,fill = "Fatigue",colour = "Fatigue"),
    alpha = 0.3,binwidth = 1)
  plotobj23 <- plotobj23 + scale_y_continuous(
    "Percentage of Population",
    lim = c(0,100),breaks = seq(0,100,10),labels = seq(0,100,10))
  plotobj23 <- plotobj23 + scale_x_continuous(
    "Difference in Grades between Exposures",
    breaks = seq(-4,4,1),labels = c(-4,-3,-2,-1,0,"+1","+2","+3","+4"))
  plotobj23 <- plotobj23 + labs(fill = "",
    colour = "")
  plotobj23 <- plotobj23 + theme(legend.position = "bottom",
    legend.key.size = unit(0.2,units = "cm"))
  print(plotobj23)
  # ggsave(plot = plotobj23,
  #   filename = paste0("50perIncrease_HFSandFAT.png"),
  #   height = 15,width = 25,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Calculate the proportion of individuals who change Hand-foot syndrome and
# fatigue grades with a different level of sunitinib exposure
  auc1 <- 2	# Original AUC
  auc2 <- 1.5	# Target AUC
  subset.data <- pd52.data[pd52.data$AUC24r == auc1 | pd52.data$AUC24r == auc2,]
  cat.change <- function(subset.data) {
    HFS.change <- subset.data$HFS[subset.data$AUC24r == auc2]-
      subset.data$HFS[subset.data$AUC24r == auc1]
    FAT.change <- subset.data$FAT[subset.data$AUC24r == auc2]-
      subset.data$FAT[subset.data$AUC24r == auc1]
    result <- data.frame(HFS.change,FAT.change)
  }
  change.data <- ddply(subset.data, .(SIM,ID), cat.change)
  # change.title <- paste0(
  # 	"25% Decreased Exposure:\nTotal Sunitinib 24-hour AUC of ",
  # 	auc1," mg*h/L (reference) versus ",auc2," mg*h/L (target)")
  change.title <- "(d)"

  plotobj24 <- NULL
  plotobj24 <- ggplot(change.data)
  plotobj24 <- plotobj24 + ggtitle(change.title)
  plotobj24 <- plotobj24 + geom_histogram(aes(x = HFS.change,
    y = ..count../1000*100,fill = "Hand-Foot Syndrome",
    colour = "Hand-Foot Syndrome"),alpha = 0.5,binwidth = 1)
  plotobj24 <- plotobj24 + geom_histogram(aes(x = FAT.change,
    y = ..count../1000*100,fill = "Fatigue",colour = "Fatigue"),
    alpha = 0.3,binwidth = 1)
  plotobj24 <- plotobj24 + scale_y_continuous(
    "Percentage of Population",
    lim = c(0,100),breaks = seq(0,100,10),labels = seq(0,100,10))
  plotobj24 <- plotobj24 + scale_x_continuous(
    "Difference in Grades between Exposures",
    breaks = seq(-4,4,1),labels = c(-4,-3,-2,-1,0,"+1","+2","+3","+4"))
  plotobj24 <- plotobj24 + labs(fill = "",
    colour = "")
  plotobj24 <- plotobj24 + theme(legend.position = "bottom",
    legend.key.size = unit(0.2,units = "cm"))
  print(plotobj24)
  # ggsave(plot = plotobj24,
  #   filename = paste0("25perDecrease_HFSandFAT.png"),
  #   height = 15,width = 25,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Calculate the proportion of individuals who change Hand-foot syndrome and
# fatigue grades with a different level of sunitinib exposure
  auc1 <- 3	# Original AUC
  auc2 <- 1.5	# Target AUC
  subset.data <- pd52.data[pd52.data$AUC24r == auc1 | pd52.data$AUC24r == auc2,]
  cat.change <- function(subset.data) {
    HFS.change <- subset.data$HFS[subset.data$AUC24r == auc2]-
      subset.data$HFS[subset.data$AUC24r == auc1]
    FAT.change <- subset.data$FAT[subset.data$AUC24r == auc2]-
      subset.data$FAT[subset.data$AUC24r == auc1]
    result <- data.frame(HFS.change,FAT.change)
  }
  change.data <- ddply(subset.data, .(SIM,ID), cat.change)
  # change.title <- paste0(
  # 	"50% Decreased Exposure:\nTotal Sunitinib 24-hour AUC of ",
  # 	auc1," mg*h/L (reference) versus ",auc2," mg*h/L (target)")
  change.title <- "(f)"

  plotobj25 <- NULL
  plotobj25 <- ggplot(change.data)
  plotobj25 <- plotobj25 + ggtitle(change.title)
  plotobj25 <- plotobj25 + geom_histogram(aes(x = HFS.change,
    y = ..count../1000*100,fill = "Hand-Foot Syndrome",
    colour = "Hand-Foot Syndrome"),alpha = 0.5,binwidth = 1)
  plotobj25 <- plotobj25 + geom_histogram(aes(x = FAT.change,
    y = ..count../1000*100,fill = "Fatigue",colour = "Fatigue"),
    alpha = 0.3,binwidth = 1)
  plotobj25 <- plotobj25 + scale_y_continuous(
    "Percentage of Population",
    lim = c(0,100),breaks = seq(0,100,10),labels = seq(0,100,10))
  plotobj25 <- plotobj25 + scale_x_continuous(
    "Difference in Grades between Exposures",
    breaks = seq(-4,4,1),labels = c(-4,-3,-2,-1,0,"+1","+2","+3","+4"))
  plotobj25 <- plotobj25 + labs(fill = "",
    colour = "")
  plotobj25 <- plotobj25 + theme(legend.position = "bottom",
    legend.key.size = unit(0.2,units = "cm"))
  print(plotobj25)
  # ggsave(plot = plotobj25,
  #   filename = paste0("50perDecrease_HFSandFAT.png"),
  #   height = 15,width = 25,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Combine plots together
  plotobj26 <- plotobj23 + theme(legend.position = "none")
  plotobj27 <- plotobj24 + theme(legend.position = "none")
  plotobj29 <- grid.arrange(grobs = list(plotobj26,plotobj27,plotobj25),
    heights = c(5,5,6))
  print(plotobj29)
  ggsave(plot = plotobj29,filename = paste0("CatExposureResponse.png"),
    width = 11.1,height = 12.6,units = "in",dpi = 300)

# ------------------------------------------------------------------------------
# Combine plots together
  plotobj30 <- grid.arrange(grobs = list(plotobj22,plotobj29),ncol = 2)
  ggsave(plot = plotobj30,filename = paste0("ExposureResponse.png"),
    width = 17.4,height = 21,units = "cm",dpi = 300)
  ggsave(plot = plotobj30,filename = paste0("ExposureResponse.pdf"),
    width = 17.4,height = 21,units = "cm",dpi = 300)
