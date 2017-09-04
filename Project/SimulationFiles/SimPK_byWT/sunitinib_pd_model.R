# Hansson EK, Amantea MA, Westwood P, Milligan PA, Houk BE, French J,
# Karlsson MO, Friberg LE. PKPD Modeling of VEGF, sVEGFR-2, sVEGFR-3, and sKIT
# as Predictors of Tumor Dynamics and Overall Survival Following Sunitinib
# Treatment in GIST. CPT: pharmacometrics & systems pharmacology. 2013;2(11):1-9

# Hansson EK, Ma G, Amantea MA, French J, Milligan PA, Friberg LE, Karlsson MO.
# PKPD Modeling of Predictors for Adverse Effects and Overall Survival in
# Sunitinib-Treated Patients With GIST. CPT: pharmacometrics & systems
# pharmacology. 2013;2(12):1-9.
# ------------------------------------------------------------------------------
# Define model
pd.code <- '
$SET			atol = 1e-8, rtol = 1e-8
          maxsteps = 100000000

$CMT			// List of compartments
          // Initial conditions specified in $MAIN
          VEGFR3,
          SKIT,
          SKITP,
          TUMOUR,
          ANC,
          STEM,
          TRANSIT1,
          TRANSIT2,
          TRANSIT3,
          BP,
          CHZSURV,
          CHZDROP

$PARAM		// Population parameters
          // sVEGFR-3
          POPVEGFR3BASE = 63900,	// Baseline sVEGFR-3, pg/mL
          POPVEGFR3MRT = 16.7*24,	// Mean residence time, days (converted to hours)

          // sKIT
          POPSKITBASE = 39200,	// Baseline sKIT, pg/mL
          POPSKITMRT = 101*24,	// Mean residence time, days (converted to hours)
          POPSKITSLP = 0.035/1000,	// Slope for linear sKIT progression, hours^-1

          // Tumour
          POPKG = 0.0118/24/7,	// Tumour growth rate constant, h^-1
          LAMBDA = 0.0217/24/7,	// Tumour resistance development/regrowth rate constant, h^-1
          POPKRD = 0.00503/24/7, // Tumour size reduction rate constant due to drug effect, h^-1
          POPKRSKIT = 0.00282/24/7,	// Tumour size reduction rate constant due to sKIT, h^-1
          KRVEGFR3 = 0.0371/24/7,	// Tumour size reduction rate constant due to sVEGFR-3, h^-1

          // Absolute neutrophil count
          POPANCBASE = 4.94,	// Baseline ANC
          POPANCMTT = 248,	// ANC mean transit time, units
          ANCPOWER = 0.362,	// Feedback factor for mimicking the effect of endogenous growth factors, i.e., G-CSF (ANC model value)
          POPANCEMAX = 0.520,	// Maximum effect of sVEGFR-3 on ANC
          POPANCE50 = 0.552,	// sVEGFR-3 resulting in half of maximal effect on ANC
          ANCKE = log(2)/7,	// Elimination half-life of circulating neutrophils

          // Diastolic blood pressure
          POPBPBASE = 71.8,	// Baseline diastolic blood pressure
          POPBPSLP = 0.119,	// Slope for linear diastolic blood pressure progression due to drug exposure, hours^-1
          POPBPMRT = 361,	// Mean residence time, hours

          // Hand-foot syndrome
          HFSB01 = -10.4,	// Transition from grade 0 to grade 1
          HFSB02 = -0.974,	// Transition from grade 0 to grade 2
          HFSB03 = -1.59,	// Transition from grade 0 to grade 3
          HFSB11 = 2.29,	// Staying in grade 1
          HFSB12 = -9.53,	// Transition from grade 1 to grade 2
          HFSB13 = -1.33,	// Transition from grade 1 to grade 3
          HFSB21 = 3.04,	// Transition from grade 2 to grade 1
          HFSB22 = -0.747,	// Staying in grade 2
          HFSB23 = -9.09,	// Transition from grade 2 to grade 3
          HFSB31 = 3.4,	// Transition from grade 3 to grade 1
          HFSB32 = -1.65,	// Transition from grade 3 to grade 2
          HFSB33 = -0.453,	// Staying in grade 3
          HFS0VEGFR3 = -8,	// sVEGFR-3 effect on transitions from grade 0 HFS
          HFS1VEGFR3 = -6.05,	// sVEGFR-3 effect on transitions from grade 1 HFS
          HFS2VEGFR3 = -3.23,	// sVEGFR-3 effect on transitions from grade 2 HFS
          HFS3VEGFR3 = -4.75,	// sVEGFR-3 effect on transitions from grade 3 HFS

          // Fatigue
          FATB01 = -5.85,	// Transition from grade 0 to grade 1
          FATB02 = -1.14,	// Transition from grade 0 to grade 2
          FATB03 = -1.6,	// Transition from grade 0 to grade 3
          FATB11 = 2.63,	// Staying in grade 1
          FATB12 = -10.7,	// Transition from grade 1 to grade 2
          FATB13 = -1.77,	// Transition from grade 1 to grade 3
          FATB21 = 2.86,	// Transition from grade 2 to grade 1
          FATB22 = -0.427,	// Staying in grade 2
          FATB23 = -11.6,	// Transition from grade 2 to grade 3
          FATB31 = 3.06,	// Transition from grade 3 to grade 1
          FATB32 = -0.0903,	// Transition from grade 3 to grade 2
          FATB33 = -0.636,	// Staying in grade 3
          FAT0VEGFR3 = -1.93,	// sVEGFR-3 effect on transitions from grade 0 FAT
          FAT1VEGFR3 = -4.62,	// sVEGFR-3 effect on transitions from grade 1 FAT
          FAT2VEGFR3 = -4.64,	// sVEGFR-3 effect on transitions from grade 2 FAT
          FAT3VEGFR3 = -3.32,	// sVEGFR-3 effect on transitions from grade 3 FAT

          // Biomarker overall survival model
             // Survival
             LAMBH = 0.00596/7/24,	// Scale parameter in the Weibull density function for the survival model
             ALPHH = 1.226,	// Shape parameter in the Weibull density function for the survival model
             VEGFR3HAZ = -3.77,	// Parameter relating sVEGFR-3 to the hazard
             TBASEHAZ = -0.002371,	// Parameter relating observed baseline tumour size to the hazard

             // Drop out
             LAMBD = 0.00196/7/24,	// Scale parameter in the Weibull density function for the drop out model
             ALPHD = 1.273,	// Shape parameter in the Weibull density function for the drop out model

          // Pre-defined between subject variability
          ETAVEGFR3BASE = 0,
          ETAVEGFR3MRT = 0,
          ETAVEGFR3I50 = 0,
          ETASKITBASE = 0,
          ETASKITMRT = 0,
          ETASKITI50 = 0,
          ETASKITSLP = 0,
          ETAKG = 0,
          ETAKRSKIT = 0,
          ETAKRD = 0,
          ETAOBASE = 0,
          ETAANCBASE = 0,
          ETAANCMTT = 0,
          ETAANCEMAX = 0,
          ETAANCE50 = 0,
          ETABPBASE = 0,
          ETABPSLP = 0,
          ETABPMRT = 0,
          ETAHFS0 = 0,
          ETAHFS1 = 0,
          ETAHFS2 = 0,
          ETAFAT0 = 0,
          ETAFAT1 = 0,
          ETAFAT2 = 0,
          ETAFAT3 = 0,

          // Pre-defined random numbers for overall survival
          RSURV = 0,	// Probability of event
          RDROP = 0,	// Probability of drop out

          // Covariates
          AUC24 = 2,	// 24-hour area under the curve for combined parent and metabolite, mg*h/L
          WT = 70,	// Body weight, kg
          OBASE = 200,	// Observed tumour size at baseline (length of longest diameter, mm)
          HFSBASE = 0,	// Baseline hand-foot syndrome grade
          FATBASE = 0	// Baseline fatigue grade

$OMEGA		// Between-subject variability
          // Biomarkers
          @annotated
          BSVVEGFR3BASE: 0.186 : ETA on baseline sVEGFR-3
          BSVVEGFR3MRT : 0.06 : ETA on mean residence time for sVEGFR-3
          BSVSKITBASE: 0.254 : ETA on baseline sKIT
          BSVSKITMRT: 0.0753: ETA on mean residence time for sKIT
          BSVSKITSLP: 3.01 : ETA on slope for the sKIT progression model

$OMEGA		// Between-subject variability
          // Biomarkers
          @annotated @block
          BSVVEGFR3I50: 0.398 : ETA on daily sunitinib AUC resulting in half of its maximum effect on sVEGFR-3
          BSVSKITI50: 0.936 5.77 : ETA on daily sunitinib AUC resulting in half of its maximum effect on sKIT

$OMEGA		// Between-subject variability
          // Tumour
          @annotated
          BSVKG: 0.29	: ETA on tumour growth rate constant
          BSVKRSKIT: 5.91	: ETA on tumour size reduction rate constant due to sKIT
          BSVKRD: 1.41	: ETA on tumour size reduction rate constant due to drug effect
          BSVOBASE: 0.015625 : Proportional error in observed baseline tumour size

$OMEGA		// Between-subject variability
          // Absolute neutrophil count
          @annotated @ block
          BSVANCBASE: 0.172 : ETA for baseline ANC
          BSVANCMTT: 0.0477 0.0155	: ETA for ANC mean transit time
          BSVANCEMAX: 0 0.000914 0.0279	: ETA for ANC Emax

$OMEGA		// Between-subject variability
          // Absolute neutrophil count
          @annotated
          BSVANCE50: 0.209	: ETA for ANC EC50

$OMEGA		// Between-subject variability
          // Diastolic blood pressure
          @annotated @block
          BSVBPBASE: 0.0151 : ETA on baseline diastolic blood pressure
          BSVBPSLP: -0.0515 0.416 : ETA on linear progression of BP due to drug exposure

$OMEGA		// Between-subject variability
          // Diastolic blood pressure
          @annotated
          BSVBPMRT: 0.694 : ETA on mean residence time for BP

$OMEGA		// Between-subject variability
          // Hand-foot syndrome
          @annotated
          BSVHFS0: 9.42 : ETA on transition from grade 0
          BSVHFS1: 0.813 : ETA on transition from grade 1
          BSVHFS2: 0.0823 : ETA on transition from grade 2

$OMEGA		// Between-subject variability
          // Fatigue
          @annotated
          BSVFAT0: 1.12 : ETA on transition from grade 0
          BSVFAT1: 1.57 : ETA on transition from grade 1
          BSVFAT2: 1.68 : ETA on transition from grade 2
          BSVFAT3: 0.707281	: ETA on transition from grade 3

$MAIN			// Individual parameter values
          // sVEGFR-3
          double VEGFR3BASE = POPVEGFR3BASE*exp(ETAVEGFR3BASE);
          double VEGFR3MRT = POPVEGFR3MRT*exp(ETAVEGFR3MRT);
          double VEGFR3I50 = exp(ETAVEGFR3I50);

          // sKIT
          double SKITBASE = POPSKITBASE*exp(ETASKITBASE);
          double SKITMRT = POPSKITMRT*exp(ETASKITMRT);
          double SKITI50 = exp(ETASKITI50);
          double SKITSLP = POPSKITSLP*exp(ETASKITSLP);

          // Tumour
          double TBASE = OBASE*(1+ETAOBASE);
          double KG = POPKG*exp(ETAKG);
          double KRSKIT = POPKRSKIT*exp(ETAKRSKIT);
          double KRD = POPKRD*exp(ETAKRD);

          // Absolute neutrophil count
          double ANCBASE = POPANCBASE*exp(ETAANCBASE);
          double ANCMTT = POPANCMTT*exp(ETAANCMTT);
          double ANCEMAX = POPANCEMAX*exp(ETAANCEMAX);
          double ANCE50 = POPANCE50*exp(ETAANCE50);
          double ANCKTR = 4/ANCMTT;

          // Diastolic blood pressure
          double BPBASE = POPBPBASE*exp(ETABPBASE);
          double BPSLP = POPBPSLP*exp(ETABPSLP);
          double BPMRT = POPBPMRT*exp(ETABPMRT);

          // Natural linear progression of biomarker over time
          double VEGFR3DP = VEGFR3BASE;
          double SKITDP = SKITBASE*(1+SKITSLP*TIME);

          // Compartment initial conditions
          VEGFR3_0 = VEGFR3BASE;
          SKIT_0 = SKITBASE;
          SKITP_0 = SKITBASE;	// sKIT time-course in untreated patients
          TUMOUR_0 = TBASE;
          ANC_0 = ANCBASE;
          STEM_0 = (ANCKE*ANCBASE)/ANCKTR;
          TRANSIT1_0 = (ANCKE*ANCBASE)/ANCKTR;
          TRANSIT2_0 = (ANCKE*ANCBASE)/ANCKTR;
          TRANSIT3_0 = (ANCKE*ANCBASE)/ANCKTR;
          BP_0 = BPBASE;
          CHZSURV_0 = 0;
          CHZDROP_0 = 0;

$ODE			// Rate constants
          double VEGFR3KIN = VEGFR3DP*(1/VEGFR3MRT)*(1-((1*AUC24)/(VEGFR3I50+AUC24)));
          double VEGFR3KOUT = 1/VEGFR3MRT;
          double SKITKIN = SKITDP*(1/SKITMRT)*(1-((1*AUC24)/(SKITI50+AUC24)));
          double SKITKOUT = 1/SKITMRT;
          double BPKIN = BPBASE*(1/BPMRT)*(1+AUC24*BPSLP);
          double BPKOUT = 1/BPMRT;

          // Biomarkers
          dxdt_VEGFR3 = VEGFR3KIN -VEGFR3KOUT*VEGFR3;
          dxdt_SKIT = SKITKIN -SKITKOUT*SKIT;
          dxdt_SKITP = SKITKIN/(1-((1*AUC24)/(SKITI50+AUC24))) -SKITKOUT*SKITP;

          // Relative changes in biomarkers and their effect on tumour growth
          double rSKIT = (SKIT-SKITP)/SKITP;
          double rVEGFR3 = (VEGFR3-VEGFR3BASE)/VEGFR3BASE;

          // Tumour
          dxdt_TUMOUR = KG*TUMOUR -(KRD*AUC24 -KRSKIT*rSKIT -KRVEGFR3*rVEGFR3)*exp(-LAMBDA*SOLVERTIME)*TUMOUR;

          // sVEGFR-3 effect on absolute neutrophil count
          // Sunitinib AUC is indirectly related
          double eVEGFR3 = ANCEMAX*-rVEGFR3/(ANCE50-rVEGFR3);	// ANC model

          // Absolute neutrophil count
          dxdt_ANC = ANCKTR*TRANSIT3 -ANCKE*ANC;
          dxdt_STEM = ANCKTR*STEM*(1-eVEGFR3)*pow(ANCBASE/ANC,ANCPOWER) -ANCKTR*STEM;
          dxdt_TRANSIT1 = ANCKTR*STEM -ANCKTR*TRANSIT1;
          dxdt_TRANSIT2 = ANCKTR*TRANSIT1 -ANCKTR*TRANSIT2;
          dxdt_TRANSIT3 = ANCKTR*TRANSIT2 -ANCKTR*TRANSIT3;

          // Diastolic blood pressure
          dxdt_BP = BPKIN -BPKOUT*BP;

          // Convert time-scale to weeks
          double TIMEW = SOLVERTIME/24/7;

          // Cumulative hazard for survival
          dxdt_CHZSURV = 0;
          if (TIMEW > 4) dxdt_CHZSURV = LAMBH*ALPHH*pow(TIMEW,ALPHH-1)*exp(-VEGFR3HAZ*rVEGFR3-TBASEHAZ*TBASE);

          // Cumulative hazard for survival from drop out
          dxdt_CHZDROP = 0;
          if (TIMEW > 4) dxdt_CHZDROP = LAMBD*ALPHD*pow(TIMEW,ALPHD-1);

$TABLE		// Biomarkers
          double IPRE_VEGFR3 = log(VEGFR3);
          double IPRE_SKIT = log(SKIT);

          // Hand-foot syndrome
            // sVEGFR-3 effects on transition probabilities
            double HFS0VEGFR3EFF = HFS0VEGFR3*rVEGFR3;
            double HFS1VEGFR3EFF = HFS1VEGFR3*rVEGFR3;
            double HFS2VEGFR3EFF = HFS2VEGFR3*rVEGFR3;
            double HFS3VEGFR3EFF = HFS3VEGFR3*rVEGFR3;

            // Probability of transition from grade 0
            double HFSLG01 = HFSB01+HFS0VEGFR3EFF+ETAHFS0;
            double HFSLG02 = HFSLG01+HFSB02;
            double HFSLG03 = HFSLG02+HFSB03;
            double HFSP01 = exp(HFSLG01)/(1+exp(HFSLG01));
            double HFSP02 = exp(HFSLG02)/(1+exp(HFSLG02));
            double HFSP03 = exp(HFSLG03)/(1+exp(HFSLG03));
            double HFSPROB00 = 1-HFSP01;
            double HFSPROB01 = HFSP01-HFSP02;
            double HFSPROB02 = HFSP02-HFSP03;
            double HFSPROB03 = HFSP03;

            // Probability of transition from grade 1
            double HFSLG11 = HFSB11+HFS1VEGFR3EFF+ETAHFS1;
            double HFSLG12 = HFSLG11+HFSB12;
            double HFSLG13 = HFSLG12+HFSB13;
            double HFSP11 = exp(HFSLG11)/(1+exp(HFSLG11));
            double HFSP12 = exp(HFSLG12)/(1+exp(HFSLG12));
            double HFSP13 = exp(HFSLG13)/(1+exp(HFSLG13));
            double HFSPROB10 = 1-HFSP11;
            double HFSPROB11 = HFSP11-HFSP12;
            double HFSPROB12 = HFSP12-HFSP13;
            double HFSPROB13 = HFSP13;

            // Probability of transition from grade 2
            double HFSLG21 = HFSB21+HFS2VEGFR3EFF+ETAHFS2;
            double HFSLG22 = HFSLG21+HFSB22;
            double HFSLG23 = HFSLG22+HFSB23;
            double HFSP21 = exp(HFSLG21)/(1+exp(HFSLG21));
            double HFSP22 = exp(HFSLG22)/(1+exp(HFSLG22));
            double HFSP23 = exp(HFSLG23)/(1+exp(HFSLG23));
            double HFSPROB20 = 1-HFSP21;
            double HFSPROB21 = HFSP21-HFSP22;
            double HFSPROB22 = HFSP22-HFSP23;
            double HFSPROB23 = HFSP23;

            // Probability of transition from grade 3
            double HFSLG31 = HFSB31+HFS3VEGFR3EFF;
            double HFSLG32 = HFSLG31+HFSB32;
            double HFSLG33 = HFSLG32+HFSB33;
            double HFSP31 = exp(HFSLG31)/(1+exp(HFSLG31));
            double HFSP32 = exp(HFSLG32)/(1+exp(HFSLG32));
            double HFSP33 = exp(HFSLG33)/(1+exp(HFSLG33));
            double HFSPROB30 = 1-HFSP31;
            double HFSPROB31 = HFSP31-HFSP32;
            double HFSPROB32 = HFSP32-HFSP33;
            double HFSPROB33 = HFSP33;

          // Fatigue
            // sVEGFR-3 effects on transition probabilities
            double FAT0VEGFR3EFF = FAT0VEGFR3*rVEGFR3;
            double FAT1VEGFR3EFF = FAT1VEGFR3*rVEGFR3;
            double FAT2VEGFR3EFF = FAT2VEGFR3*rVEGFR3;
            double FAT3VEGFR3EFF = FAT3VEGFR3*rVEGFR3;

            // Probability of transition from grade 0
            double FATLG01 = FATB01+FAT0VEGFR3EFF+ETAFAT0;
            double FATLG02 = FATLG01+FATB02;
            double FATLG03 = FATLG02+FATB03;
            double FATP01 = exp(FATLG01)/(1+exp(FATLG01));
            double FATP02 = exp(FATLG02)/(1+exp(FATLG02));
            double FATP03 = exp(FATLG03)/(1+exp(FATLG03));
            double FATPROB00 = 1-FATP01;
            double FATPROB01 = FATP01-FATP02;
            double FATPROB02 = FATP02-FATP03;
            double FATPROB03 = FATP03;

            // Probability of transition from grade 1
            double FATLG11 = FATB11+FAT1VEGFR3EFF+ETAFAT1;
            double FATLG12 = FATLG11+FATB12;
            double FATLG13 = FATLG12+FATB13;
            double FATP11 = exp(FATLG11)/(1+exp(FATLG11));
            double FATP12 = exp(FATLG12)/(1+exp(FATLG12));
            double FATP13 = exp(FATLG13)/(1+exp(FATLG13));
            double FATPROB10 = 1-FATP11;
            double FATPROB11 = FATP11-FATP12;
            double FATPROB12 = FATP12-FATP13;
            double FATPROB13 = FATP13;

            // Probability of transition from grade 2
            double FATLG21 = FATB21+FAT2VEGFR3EFF+ETAFAT2;
            double FATLG22 = FATLG21+FATB22;
            double FATLG23 = FATLG22+FATB23;
            double FATP21 = exp(FATLG21)/(1+exp(FATLG21));
            double FATP22 = exp(FATLG22)/(1+exp(FATLG22));
            double FATP23 = exp(FATLG23)/(1+exp(FATLG23));
            double FATPROB20 = 1-FATP21;
            double FATPROB21 = FATP21-FATP22;
            double FATPROB22 = FATP22-FATP23;
            double FATPROB23 = FATP23;

            // Probability of transition from grade 3
            double FATLG31 = FATB31+FAT3VEGFR3EFF+ETAFAT3;
            double FATLG32 = FATLG31+FATB32;
            double FATLG33 = FATLG32+FATB33;
            double FATP31 = exp(FATLG31)/(1+exp(FATLG31));
            double FATP32 = exp(FATLG32)/(1+exp(FATLG32));
            double FATP33 = exp(FATLG33)/(1+exp(FATLG33));
            double FATPROB30 = 1-FATP31;
            double FATPROB31 = FATP31-FATP32;
            double FATPROB32 = FATP32-FATP33;
            double FATPROB33 = FATP33;

          // Survival functions
          double SURV = exp(-CHZSURV);	// Survival probability
          double DROP = exp(-CHZDROP);	// Drop out probability

$CAPTURE	AUC24 WT OBASE HFSBASE FATBASE
          IPRE_VEGFR3 IPRE_SKIT
          TIMEW SURV DROP RSURV RDROP
          HFSPROB00 HFSPROB01 HFSPROB02 HFSPROB03
          HFSPROB10 HFSPROB11 HFSPROB12 HFSPROB13
          HFSPROB20 HFSPROB21 HFSPROB22 HFSPROB23
          HFSPROB30 HFSPROB31 HFSPROB32 HFSPROB33
          FATPROB00 FATPROB01 FATPROB02 FATPROB03
          FATPROB10 FATPROB11 FATPROB12 FATPROB13
          FATPROB20 FATPROB21 FATPROB22 FATPROB23
          FATPROB30 FATPROB31 FATPROB32 FATPROB33
'
# Compile the model code
  pd.mod <- mcode("popPD",pd.code)

# ------------------------------------------------------------------------------
# Pull out variance parameters from model file
  pd.OMEGA <- as.matrix(omat(pd.mod))
