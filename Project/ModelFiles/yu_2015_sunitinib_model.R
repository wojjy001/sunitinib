# Yu H, Steeghs N, Kloth JS, Wit D, Hasselt J, Erp NP, Beijnen JH, Schellens JH,
# Mathijssen RH, Huitema AD. Integrated semi‐physiological pharmacokinetic model
# for both sunitinib and its active metabolite SU12662. British journal of
# clinical pharmacology. 2015;79(5):809-19
# ------------------------------------------------------------------------------
# Convert model's standard deviations and correlations to a variance-covariance
# matrix
# Define standard deviations for between-subject variability
  BSVCLP <- 0.339	# Clearance of parent, ETACLP (1)
  BSVVCP <- 0.324	# Central volume of parent, ETAVCP (2)
  BSVCLM <- 0.421	# Clearance of metabolite, ETACLM (3)
  BSVVCM <- 0.579	# Central volume of metabolite, ETAVCM (4)
  SD <- c(BSVCLP,BSVVCP,BSVCLM,BSVVCM)
# Define correlation matrix for between-subject variability
  R12 <- 0	# Clearance parent - central volume parent
  R13 <- 0.53	# Clearance parent - clearance metabolite
  R14 <- 0	# Clearance parent - central volume metabolite
  R23 <- 0	# Central volume parent - clearance metabolite
  R24 <- 0.48	# Central volume parent - central volume metabolite
  R34 <- 0.45	# Clearance metabolite - central volume metabolite
  CORR <- matrix(c(
    1,R12,R13,R14,
    R12,1,R23,R24,
    R13,R14,1,R34,
    R14,R24,R34,1
  ),4,4)
# Convert to variance-covariance matrix
  OMEGA <- cor2cov(CORR,SD)
# Define proportional residual error for measured concentrations
  ERRPROP1 <- sqrt(0.06)
  ERRPROM1 <- sqrt(0.03)

# ------------------------------------------------------------------------------
# Define model
pk.code <- '
$SET			atol = 1e-8, rtol = 1e-8
          maxsteps = 100000

$INIT			// Initial conditions for PK compartments
          // Parent, sunitinib
          DEPOT = 0,	// Depot compartment for dose
          CENTP = 0,	// Central compartment

          // Metabolite, SU12662
          CENTM = 0,	// Central compartment
          PERIM = 0,	// Peripheral compartment

          // Area under the curve
          AUC = 0	// Combined parent and metabolite

$PARAM		// Population parameters
          // Parent, sunitinib
          POPKA = 0.34,	// Absorption rate constant, h^-1
          POPCLP = 35.7,	// Clearance, L/h
          POPVCP = 1360,	// Central volume, L
          POPQH = 80,	// Hepatic blood flow, L/h

          // Metabolite, SU12662
          POPCLM = 17.1, // Clearance, L/h
          POPVCM = 635,	// Central volume, L
          POPQI = 20.1,	// Inter-compartmental clearance, L/h
          POPVPM = 388,	// Peripheral volume, L

          // Pre-defined between subject variability
          ETACLP = 0,	// ETA on clearance of parent
          ETAVCP = 0, // ETA on central volume of parent
          ETACLM = 0, // ETA on clearance of metabolite
          ETAVCM = 0, // ETA on central volume of metabolite

          // Covariates
          WT = 70,	// Body weight, kg
          OBASE = 200,	// Observed tumour size at baseline (length of longest diameter, mm)
          HFSBASE = 0,	// Baseline hand-foot syndrome grade
          FATBASE = 0	// Baseline fatigue grade

$OMEGA		// Between-subject variability
          @annotated @block
          BSVCLP:	0.114921 : ETA on clearance of parent
          BSVVCP:	0 0.104976 : ETA on central volume of parent
          BSVCLM:	0.07564107 0 0.1777241 : ETA on clearance of metabolite
          BSVVCM:	0 0.09004608 0.10969155 0.335241 : ETA on central volume of metabolite

$SIGMA		// Random unexplained variability
          @annotated
          ERRPROP1:	0.06 : Proportional error for parent in study 1
          ERRPROP2: 0.01 : Proportional error for parent in study 2 and study 3
          ERRPROM1: 0.03 : Proportional error for metabolite in study 1
          ERRPROM2: 0.01 : Proportional error for metabolite in study 2 and study 3

$MAIN			// Individual PK parameter values
          // Parent, sunitinib
          double CLP = POPCLP*pow(WT/70,0.75)*exp(ETACLP);
          double VCP = POPVCP*pow(WT/70,1)*exp(ETAVCP);
          double QH = POPQH*pow(WT/70,0.75);

          // Metabolite, SU12662
          double CLM = POPCLM*pow(WT/70,0.75)*exp(ETACLM);
          double VCM = POPVCM*pow(WT/70,1)*exp(ETAVCM);
          double QI = POPQI*pow(WT/70,0.75);
          double VPM = POPVPM*pow(WT/70,1);

          // Hepatic extraction
          double FH = QH/(QH+CLP);	// Parent fraction not metabolised by liver
          double EH = 1-FH;	// Parent fraction metabolised by liver
          double FM = 0.21;	// Parent fraction metabolised by liver to metabolite

          // Micro-rate constants
          double K12 = POPKA;
          double K20 = CLP/VCP;
          double K23 = QH/VCP;
          double K30 = CLM/VCM;
          double K34 = QI/VCM;
          double K43 = QI/VPM;

$ODE			// Parent, sunitinib
          dxdt_DEPOT = -K12*DEPOT;
          dxdt_CENTP = FH*(K12*DEPOT +K23*CENTP) -K23*CENTP;

          // Metabolite, SU12662
          dxdt_CENTM = FM*EH*(K12*DEPOT +K23*CENTP) -K34*CENTM +K43*PERIM -K30*CENTM;
          dxdt_PERIM = K34*CENTM -K43*PERIM;

          // Area under the curve
          dxdt_AUC = CENTP/VCP+CENTM/VCM;	// Combined parent and metabolite

$TABLE		// Calculate parent and metabolite concentrations
          double IPREP = CENTP/VCP;	// mg/L
          double IPREM = CENTM/VCM;	// mg/L
          double IPRE = IPREP+IPREM;	// mg/L
          double DVP = IPREP*(1+ERRPROP1);	// mg/L
          double DVM = IPREM*(1+ERRPROM1);	// mg/L
          double DV = IPRE*(1+ERRPROP1);	// mg/L

$CAPTURE	WT OBASE HFSBASE FATBASE
          IPREP IPREM IPRE DVP DVM DV
          CLP VCP CLM VCM QI VPM
'
# Compile the model code
  pk.mod <- mcode("popPK",pk.code)

# ------------------------------------------------------------------------------
# Pull out variance parameters from model file
  pk.OMEGA <- as.matrix(omat(pk.mod))
