# Sensitivity Analysis
# This script reads in previously simulated sunitinib concentrations and
# simulates biomarkers, adverse effects and overall survival
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
# Create output directory specifically for sensitivity analysis (if not already)
  sens.output.dir <- paste0(output.dir,"SensitivityAnalysis/70kg")
  dir.create(sens.output.dir)

# ------------------------------------------------------------------------------
# Source required files
  source(paste0(global.dir,"functions.R"))	# Universal functions file
  source(paste0(global.dir,"ModelFiles/hansson_2013_sunitinib_models.R")) # PD

# ------------------------------------------------------------------------------
# Read in simulation files in the output directory
# List the folders in the output directory
  folder.list <- list.dirs(sens.output.dir)
  folder.list <- folder.list[contains("target",vars = folder.list)]
  study.list <- folder.list
  for (i in 1:length(folder.list)) {
    folder.list[i] <- str_split(folder.list[i],pattern = "/")
    study.list[i] <- folder.list[[i]][7]
  }
  # study.list <- head(study.list,2)
# OR just list the one folder you want to process
  # study.list <- "target_standard_NA_70kg"

# ------------------------------------------------------------------------------
# Simulate and process pd output for each pk study design
  for (i in 1:length(study.list)) {
    study.dir <- paste0(sens.output.dir,"/",study.list[i],"/")
    setwd(study.dir)
  # Read in pk .csv file from the study folder
    pk.data <- read.csv(file = paste0(study.list[i],"_pk_data.csv"))
    pk.data$study <- study.list[i]
  # Read in population characteristics from pk.data
    nid <- length(unique(pk.data$ID))
    ID.seq <- 1:nid
    nsim <- length(unique(pk.data$SIM))
    SIM.seq <- 1:nsim
    ntotal <- nid*nsim

  # Covariates
  # Baseline tumour size
    OBASE <- exp(log(rnorm(nid,mean = 195,sd = 120)))
    repeat {
      new.OBASE <- OBASE[OBASE == "NaN"]
      new.OBASE <- exp(log(rnorm(length(new.OBASE),mean = 195,sd = 120)))
      OBASE <- c(OBASE[OBASE != "NaN"],new.OBASE)
      if (length(OBASE[OBASE == "NaN"]) == 0) break
    }
  # Baseline hand-foot syndrome grade
    hfs.grade0.prob <- 0.95
    hfs.grade1.prob <- 0.02
    hfs.grade2.prob <- 0.02
    hfs.grade3.prob <- 0.01
    base.hfs.probs <- c(hfs.grade0.prob,hfs.grade1.prob,hfs.grade2.prob,
      hfs.grade3.prob)
    HFSBASE <- unlist(llply(seq_len(ntotal),
      function(x) sample(c(0,1,2,3),size = x/x,prob = base.hfs.probs)))
  # Baseline fatigue grade
    fat.grade0.prob <- 0.85
    fat.grade1.prob <- 0.09
    fat.grade2.prob <- 0.05
    fat.grade3.prob <- 0.01
    base.fat.probs <- c(fat.grade0.prob,fat.grade1.prob,fat.grade2.prob,
      fat.grade3.prob)
    FATBASE <- unlist(llply(seq_len(ntotal),
      function(x) sample(c(0,1,2,3),size = x/x,prob = base.fat.probs)))

  # Generate random effect parameters
  # Between subject variability
    if (ntotal > 1) {
      pd.ETA.matrix <- mvrnorm(ntotal,
        mu = rep(0,times = dim(pd.OMEGA)[1]),pd.OMEGA) %>%
        as.data.frame
      names(pd.ETA.matrix) <- c("ETAVEGFR3BASE","ETAVEGFR3MRT","ETAVEGFR3I50",
        "ETASKITBASE","ETASKITMRT","ETASKITI50","ETASKITSLP","ETAKG","ETAKRSKIT",
        "ETAKRD","ETAOBASE","ETAANCBASE","ETAANCMTT","ETAANCEMAX","ETAANCE50",
        "ETABPBASE","ETABPSLP","ETABPMRT","ETAHFS0","ETAHFS1","ETAHFS2","ETAFAT0",
        "ETAFAT1","ETAFAT2","ETAFAT3")
    } else {
      pd.ETA.matrix <- data.frame(ETAVEGFR3BASE = 0,ETAVEGFR3MRT = 0,
        ETASKITBASE = 0,ETASKITMRT = 0,ETASKITSLP = 0,ETAVEGFR3I50 = 0,
        ETASKITI50 = 0,ETAKG = 0,ETAKRSKIT = 0,ETAKRD = 0,ETAOBASE = 0,
        ETAANCBASE = 0,ETAANCMTT = 0,ETAANCEMAX = 0,ETAANCE50 = 0,ETABPBASE = 0,
        ETABPSLP = 0,ETABPMRT = 0,ETAHFS0 = 0,ETAHFS1 = 0,ETAHFS2 = 0,ETAFAT0 = 0,
        ETAFAT1 = 0,ETAFAT2 = 0,ETAFAT3 = 0)
    }

  # Input data frame for simulation
    pkpd.data <- pk.data[c("SIM","ID","time","cyc","amt","WT","IPREP",
      "IPREM","IPRE","AUC24","study")]
    pd.ID.data <- data.frame(SIM = SIM.seq,ID = ID.seq,OBASE,HFSBASE,FATBASE,
      pd.ETA.matrix)
    input.pd.data <- merge(pkpd.data,pd.ID.data,by = c("SIM","ID"),all = T)
    input.pd.data$cmt <- 0
    input.pd.data$evid <- 0
    input.pd.data$RSURV <- runif(nid*length(unique(pk.data$time)),min = 0,max = 1)
    input.pd.data$RDROP <- runif(nid*length(unique(pk.data$time)),min = 0,max = 1)
    input.pd.data <- input.pd.data[with(input.pd.data,order(input.pd.data$SIM,
      input.pd.data$ID,input.pd.data$time)),]

  # Simulate biomarker and adverse effect profiles
    simulate.pd <- function(input.pd.data) {
      pd.data <- pd.mod %>% mrgsim(data = input.pd.data,
        carry.out = c("SIM","cyc","amt","IPREP","IPREM","IPRE","AUC24")) %>%
        as.data.frame
    }
    pd.data <- ddply(input.pd.data, .(study), simulate.pd,
      .progress = "text")
  # Use hand-foot syndrome state probabilities to simulate HFS profile
    pd.data <- ddply(pd.data, .(study,SIM,ID), simulate.HFS.grade,
      .progress = "text")
  # Use fatigue state probabilities to simulate FAT profile
    pd.data <- ddply(pd.data, .(study,SIM,ID), simulate.FAT.grade,
      .progress = "text")
  # Determine overall survival for the population over the study period
    pd.data <- ddply(pd.data, .(study,SIM,ID), alive.function,
      .progress = "text")
    AUCTarget <- str_split(unique(pd.data$study),pattern = "_")
    for (i in 1:length(AUCTarget)) {
      AUCTarget[i] <- AUCTarget[[i]][3]
    }
    AUCTarget <- unlist(AUCTarget)
    pd.data$AUCTarget <- as.factor(pd.data$study)
    levels(pd.data$AUCTarget) <- AUCTarget
  # For biomarker and adverse effect simulations, only plot the first 50 weeks
    early.data <- pd.data[pd.data$time <= 50*24*7,]

  # For overall survival plots, calculate the proportion of individuals alive
  # at each time-point
  # Kaplan Meier Plot for survival confidence intervals
  # All individuals have the same start time, i.e., time == 0
    pd.data$start <- 0
    pd.data <- ddply(pd.data, .(AUCTarget,SIM,ID), stop.time.function)  # For each individual calculate their stop time
    km.data <- ddply(pd.data, .(AUCTarget,SIM,ID), headperID)
    km.data <- km.data[c("AUCTarget","SIM","ID","start","stop","event")]
    S <- Surv(time = km.data$start,time2 = km.data$stop,event = km.data$event)
    result <- survfit(formula = S ~ AUCTarget,data = km.data)
    cols <- lapply(2:12, function(x) summary(result)[x])
    surv.data <- do.call(data.frame, cols)
    # surv.data$AUCTarget <- as.factor(surv.data$strata)
    # levels(surv.data$AUCTarget) <- AUCTarget
  # Plot overall survival over time
    plotobj1a <- NULL
    plotobj1a <- ggplot()
    plotobj1a <- plotobj1a + geom_line(aes(x = time/24/7,y = surv),
      colour = "red",data = surv.data)
    plotobj1a <- plotobj1a + geom_ribbon(aes(x = time/24/7,ymin = lower,
      ymax = upper),data = surv.data,fill = "red",alpha = 0.3)
    plotobj1a <- plotobj1a + scale_y_continuous("Probability of Survival",
      lim = c(0,1),
      breaks = seq(from = 0,to = 1,by = 0.2),
      labels = seq(from = 0,to = 1,by = 0.2))
    plotobj1a <- plotobj1a + scale_x_continuous("Time (weeks)",
      breaks = seq(from = 0,to = max(surv.data$time)/24/7,by = 12),
      lim = c(0,max(surv.data$time/24/7)))
    plotobj1a <- plotobj1a + theme(legend.position = "none")
    print(plotobj1a)

    ggsave(plot = plotobj1a,
      filename = paste0(study.list[i],"_overallsurvival.png"),
      width = 20,height = 15,unit = "cm",dpi = 300)

  # Summarise and plot sVEGFR-3 over time
    summary.VEGFR3 <- ddply(early.data, .(AUCTarget,time),
      function(early.data) graded.summary(early.data$IPRE_VEGFR3))
  # Plot median and confidence intervals and facet
    plotobj2a <- NULL
    plotobj2a <- ggplot(summary.VEGFR3)
    plotobj2a <- plotobj2a + geom_ribbon(aes(x = time/24/7,ymin = CI90lo,
      ymax = CI90hi,fill = AUCTarget),alpha = 0.1)
    plotobj2a <- plotobj2a + geom_ribbon(aes(x = time/24/7,ymin = CI80lo,
      ymax = CI80hi,fill = AUCTarget),alpha = 0.1)
    plotobj2a <- plotobj2a + geom_ribbon(aes(x = time/24/7,ymin = CI60lo,
      ymax = CI60hi,fill = AUCTarget),alpha = 0.1)
    plotobj2a <- plotobj2a + geom_ribbon(aes(x = time/24/7,ymin = CI40lo,
      ymax = CI40hi,fill = AUCTarget),alpha = 0.1)
    plotobj2a <- plotobj2a + geom_ribbon(aes(x = time/24/7,ymin = CI20lo,
      ymax = CI20hi,fill = AUCTarget),alpha = 0.1)
    plotobj2a <- plotobj2a + geom_line(aes(x = time/24/7,y = med,
      colour = AUCTarget))
    plotobj2a <- plotobj2a + scale_y_continuous("sVEGFR-3 Concentration (pg/mL)",
      lim = c(8,12))
    plotobj2a <- plotobj2a + scale_x_continuous("Time (weeks)")
    # plotobj2a <- plotobj2a + facet_wrap(~AUCTarget,ncol = 3)
    plotobj2a <- plotobj2a + theme(legend.position = "none")
    print(plotobj2a)

    ggsave(plot = plotobj2a,
      filename = paste0(study.list[i],"_sVEGFR3.png"),
      width = 30,height = 15,unit = "cm",dpi = 300)

  # Summarise and plot sKIT over time
    summary.SKIT <- ddply(early.data, .(AUCTarget,time),
      function(early.data) graded.summary(early.data$IPRE_SKIT))
  # Plot median and confidence intervals and facet
    plotobj3a <- NULL
    plotobj3a <- ggplot(summary.SKIT)
    plotobj3a <- plotobj3a + geom_ribbon(aes(x = time/24/7,ymin = CI90lo,
      ymax = CI90hi,fill = AUCTarget),alpha = 0.1)
    plotobj3a <- plotobj3a + geom_ribbon(aes(x = time/24/7,ymin = CI80lo,
      ymax = CI80hi,fill = AUCTarget),alpha = 0.1)
    plotobj3a <- plotobj3a + geom_ribbon(aes(x = time/24/7,ymin = CI60lo,
      ymax = CI60hi,fill = AUCTarget),alpha = 0.1)
    plotobj3a <- plotobj3a + geom_ribbon(aes(x = time/24/7,ymin = CI40lo,
      ymax = CI40hi,fill = AUCTarget),alpha = 0.1)
    plotobj3a <- plotobj3a + geom_ribbon(aes(x = time/24/7,ymin = CI20lo,
      ymax = CI20hi,fill = AUCTarget),alpha = 0.1)
    plotobj3a <- plotobj3a + geom_line(aes(x = time/24/7,y = med,
      colour = AUCTarget))
    plotobj3a <- plotobj3a + scale_y_continuous("sKIT Concentration (pg/mL)",
      lim = c(8,14))
    plotobj3a <- plotobj3a + scale_x_continuous("Time (weeks)")
    # plotobj3a <- plotobj3a + facet_wrap(~AUCTarget,ncol = 3)
    plotobj3a <- plotobj3a + theme(legend.position = "none")
    print(plotobj3a)

    ggsave(plot = plotobj3a,
      filename = paste0(study.list[i],"_sKIT.png"),
      width = 30,height = 15,unit = "cm",dpi = 300)

  # Summarise and plot tumour size over time
    summary.TUMOUR <- ddply(early.data, .(AUCTarget,time),
      function(early.data) graded.summary(early.data$TUMOUR))
  # Plot median and confidence intervals and facet
    plotobj4a <- NULL
    plotobj4a <- ggplot(summary.TUMOUR)
    plotobj4a <- plotobj4a + geom_ribbon(aes(x = time/24/7,ymin = CI90lo,
      ymax = CI90hi,fill = AUCTarget),alpha = 0.1)
    plotobj4a <- plotobj4a + geom_ribbon(aes(x = time/24/7,ymin = CI80lo,
      ymax = CI80hi,fill = AUCTarget),alpha = 0.1)
    plotobj4a <- plotobj4a + geom_ribbon(aes(x = time/24/7,ymin = CI60lo,
      ymax = CI60hi,fill = AUCTarget),alpha = 0.1)
    plotobj4a <- plotobj4a + geom_ribbon(aes(x = time/24/7,ymin = CI40lo,
      ymax = CI40hi,fill = AUCTarget),alpha = 0.1)
    plotobj4a <- plotobj4a + geom_ribbon(aes(x = time/24/7,ymin = CI20lo,
      ymax = CI20hi,fill = AUCTarget),alpha = 0.1)
    plotobj4a <- plotobj4a + geom_line(aes(x = time/24/7,y = med,
      colour = AUCTarget))
    plotobj4a <- plotobj4a + scale_y_continuous("Sum of Longest Diameters (mm)")
    plotobj4a <- plotobj4a + scale_x_continuous("Time (weeks)")
    # plotobj4a <- plotobj4a + facet_wrap(~AUCTarget,ncol = 3)
    plotobj4a <- plotobj4a + theme(legend.position = "none")
    print(plotobj4a)

    ggsave(plot = plotobj4a,filename = paste0(study.list[i],"_tumour_facet.png"),
      width = 30,height = 15,unit = "cm",dpi = 300)

  # Summarise and plot absolute neutrophil count over time
    summary.ANC <- ddply(early.data, .(AUCTarget,time),
      function(early.data) graded.summary(early.data$ANC))
  # Plot median and confidence intervals and facet
    plotobj5a <- NULL
    plotobj5a <- ggplot(summary.ANC)
    plotobj5a <- plotobj5a + geom_ribbon(aes(x = time/24/7,ymin = CI90lo,
      ymax = CI90hi,fill = AUCTarget),alpha = 0.1)
    plotobj5a <- plotobj5a + geom_ribbon(aes(x = time/24/7,ymin = CI80lo,
      ymax = CI80hi,fill = AUCTarget),alpha = 0.1)
    plotobj5a <- plotobj5a + geom_ribbon(aes(x = time/24/7,ymin = CI60lo,
      ymax = CI60hi,fill = AUCTarget),alpha = 0.1)
    plotobj5a <- plotobj5a + geom_ribbon(aes(x = time/24/7,ymin = CI40lo,
      ymax = CI40hi,fill = AUCTarget),alpha = 0.1)
    plotobj5a <- plotobj5a + geom_ribbon(aes(x = time/24/7,ymin = CI20lo,
      ymax = CI20hi,fill = AUCTarget),alpha = 0.1)
    plotobj5a <- plotobj5a + geom_line(aes(x = time/24/7,y = med,
      colour = AUCTarget))
    plotobj5a <- plotobj5a + scale_y_log10(
      "Absolute Neutrophil Count (x10^9)")
    plotobj5a <- plotobj5a + scale_x_continuous("Time (weeks)")
    # plotobj5a <- plotobj5a + facet_wrap(~AUCTarget,ncol = 3)
    plotobj5a <- plotobj5a + theme(legend.position = "none")
    print(plotobj5a)

    ggsave(plot = plotobj5a,filename = paste0(study.list[i],"_anc.png"),
      width = 30,height = 15,unit = "cm",dpi = 300)

  # Summarise and plot diastolic blood pressure over time
    summary.BP <- ddply(early.data, .(AUCTarget,time),
      function(early.data) graded.summary(early.data$BP))
  # Plot median and confidence intervals and facet
    plotobj6a <- NULL
    plotobj6a <- ggplot(summary.BP)
    plotobj6a <- plotobj6a + geom_ribbon(aes(x = time/24/7,ymin = CI90lo,
      ymax = CI90hi,fill = AUCTarget),alpha = 0.1)
    plotobj6a <- plotobj6a + geom_ribbon(aes(x = time/24/7,ymin = CI80lo,
      ymax = CI80hi,fill = AUCTarget),alpha = 0.1)
    plotobj6a <- plotobj6a + geom_ribbon(aes(x = time/24/7,ymin = CI60lo,
      ymax = CI60hi,fill = AUCTarget),alpha = 0.1)
    plotobj6a <- plotobj6a + geom_ribbon(aes(x = time/24/7,ymin = CI40lo,
      ymax = CI40hi,fill = AUCTarget),alpha = 0.1)
    plotobj6a <- plotobj6a + geom_ribbon(aes(x = time/24/7,ymin = CI20lo,
      ymax = CI20hi,fill = AUCTarget),alpha = 0.1)
    plotobj6a <- plotobj6a + geom_line(aes(x = time/24/7,y = med,
      colour = AUCTarget))
    plotobj6a <- plotobj6a + scale_y_continuous(
      "Diastolic Blood Pressure (mmHg)")
    plotobj6a <- plotobj6a + scale_x_continuous("Time (weeks)")
    # plotobj6a <- plotobj6a + facet_wrap(~AUCTarget,ncol = 3)
    plotobj6a <- plotobj6a + theme(legend.position = "none")
    print(plotobj6a)

    ggsave(plot = plotobj6a,filename = paste0(study.list[i],"_bp.png"),
      width = 30,height = 15,unit = "cm",dpi = 300)

  # Summarise and plot probability of hand-foot syndrome grade over time
  # Summary count function
    summary.count.function <- function(x) {
      total.n <- ntotal
      n <- length(x)
      result <- n/total.n
      names(result) <- "pro"
      result
    }
    summary.HFS <- ddply(early.data, .(AUCTarget,time,HFS),
      function(early.data) summary.count.function(early.data$HFS))
    summary.HFS$Grade <- as.factor(summary.HFS$HFS)
    plotobj7 <- NULL
    plotobj7 <- ggplot(summary.HFS)
    plotobj7 <- plotobj7 + geom_line(aes(x = time/24/7,y = pro,colour = Grade))
    plotobj7 <- plotobj7 + scale_y_continuous(
      "Probability of Hand-Foot Syndrome Grade",
      lim = c(0,1),
      breaks = seq(from = 0,to = 1,by = 0.2),
      labels = seq(from = 0,to = 1,by = 0.2))
    plotobj7 <- plotobj7 + scale_x_continuous("Time (weeks)")
    # plotobj7 <- plotobj7 + facet_wrap(~AUCTarget,ncol = 3)
    print(plotobj7)

    ggsave(plot = plotobj7,filename = paste0(study.list[i],"_HFS.png"),
      width = 30,height = 15,unit = "cm",dpi = 300)

  # Summarise and plot probability of fatigue grade over time
    summary.FAT <- ddply(early.data, .(AUCTarget,time,FAT),
      function(early.data) summary.count.function(early.data$FAT))
    summary.FAT$Grade <- as.factor(summary.FAT$FAT)
    plotobj8 <- NULL
    plotobj8 <- ggplot(summary.FAT)
    plotobj8 <- plotobj8 + geom_line(aes(x = time/24/7,y = pro,colour = Grade))
    plotobj8 <- plotobj8 + scale_y_continuous("Probability of Fatigue Grade",
      lim = c(0,1),
      breaks = seq(from = 0,to = 1,by = 0.2),
      labels = seq(from = 0,to = 1,by = 0.2))
    plotobj8 <- plotobj8 + scale_x_continuous("Time (weeks)")
    # plotobj8 <- plotobj8 + facet_wrap(~AUCTarget,ncol = 3)
    print(plotobj8)

    ggsave(plot = plotobj8,filename = paste0(study.list[i],"_FAT.png"),
      width = 30,height = 15,unit = "cm",dpi = 300)

  # Clean up and save the PD simulated data
    output.pd.data <- pd.data[c("study","SIM","ID","time","cyc","amt","IPREP",
      "IPREM","IPRE","AUC24","TUMOUR","ANC","BP","WT","OBASE","HFSBASE",
      "FATBASE","IPRE_VEGFR3","IPRE_SKIT","HFS","FAT","status")]
    write.csv(output.pd.data,file = paste0(study.list[i],"_pd_data.csv"),
      quote = FALSE,row.names = FALSE)
  }
