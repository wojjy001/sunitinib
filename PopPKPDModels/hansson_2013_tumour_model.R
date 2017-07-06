# Hansson EK, Amantea MA, Westwood P, Milligan PA, Houk BE, French J,
# Karlsson MO, Friberg LE. PKPD Modeling of VEGF, sVEGFR-2, sVEGFR-3, and sKIT
# as Predictors of Tumor Dynamics and Overall Survival Following Sunitinib
# Treatment in GIST. CPT: pharmacometrics & systems pharmacology. 2013;2(11):1-9
# ------------------------------------------------------------------------------
# Define model
tum.code <- '
$SET			atol = 1e-8, rtol = 1e-8
          maxsteps = 100000000

$CMT			// List of compartments
          // Initial conditions specified in $MAIN
          VEGFR3,
          SKIT,
          SKITP,
          TUMOUR

$PARAM		// Population parameters
          // TUMOUR
          WTUM = 0.125,	// Residual error of baseline tumour sized
          POPKG = 0.0118/24/7,	// Tumour growth rate constant, h^-1
          POPKSKIT = 0.00282/24/7,	// Tumour size reduction rate constant, h^-1
          POPLAMBDA = 0.0217/24/7,	// Resistance appearance rate constant, h^-1
          POPKDRUG = 0.00503/24/7,	// Exposure driven drug effect, h^-1
          POPKVEGFR3 = 0.0371/24/7,	// Tumour size reduction rate constant, h^-1

          // Individual biomarker parameters
          // Common to some biomarkers
          IMAX = 1,

          // VEGFR3
          BM03 = 44085,
          IC503 = 0.94008,
          MRT3 = 392.95,

          // SKIT
          BM0S = 25731,
          IC50S = 0.27426,
          MRTS = 2353,
          DPSLOS = 0.267/1000,

          // Covariates
          OBASE = 70,	// Observed tumour size at baseline

          // Pre-defined between subject variability
          ETAKG = 0,
          ETAKSKIT = 0,
          ETAKDRUG = 0,

          // Pre-defined random unexplained variability
          EPSBASE = 0,

          // Pharmacokinetic parameters
          DOSE = 50,
          CL = 34.5

$OMEGA		// Between-subject variability
          @name tumour @annotated
          BSVKG:	0.29	: ETA on tumour growth rate constant
          BSVKSKIT:	5.91	: ETA on SKIT tumour size reduction rate constant
          BSVKDRUG:	1.41	: ETA on exposure driven drug effect

$SIGMA		// Random unexplained variability
          @annotated
          ERRBASE:	1	: Proportional error in observed baseline tumour size

$MAIN			// Individual parameter values
          // TUMOUR
          double W1 = WTUM*OBASE;
          double IBASE = OBASE+EPSBASE*W1;
          double KG = POPKG*exp(ETAKG);
          double KSKIT = POPKSKIT*exp(ETAKSKIT);
          double LAMBDA = POPLAMBDA;
          double KDRUG = POPKDRUG*exp(ETAKDRUG);
          double KVEGFR3 = POPKVEGFR3;

          // Compartment initial conditions
          VEGFR3_0 = BM03;
          SKIT_0 = BM0S;
          SKITP_0 = BM0S;
          TUMOUR_0 = IBASE;

          // Rate constants
          double KOUT3 = 1/MRT3;
          double KOUTS = 1/MRTS;

          // Sunitinib exposure
          double AUC = DOSE/CL;

$ODE			// Pharmacodynamics - sunitinib effect on biomarkers
          // VEGFR3
          double EFF3 = IMAX*AUC/(IC503+AUC);
          double KIN3 = BM03*KOUT3;

          // SKIT
          double EFFS = IMAX*AUC/(IC50S+AUC);
          double DPS = BM0S*(1+DPSLOS*SOLVERTIME);
          double KINS = DPS*KOUTS;

          // Pharmacodynamics - biomarker and sunitinib effects on tumour size
          double EFFSKIT = ((SKIT-SKITP)/SKITP)*KSKIT;
          double EFFVEGFR3 = ((VEGFR3-BM03)/BM03)*KVEGFR3;
          double AUC1 = AUC*KDRUG;

          // Differential equations
          dxdt_VEGFR3 = KIN3*(1-EFF3)-KOUT3*VEGFR3;
          dxdt_SKIT = KINS*(1-EFFS)-KOUTS*SKIT;
          dxdt_SKITP = KINS-KOUTS*SKITP;
          dxdt_TUMOUR = KG*TUMOUR -(AUC1-EFFSKIT-EFFVEGFR3)*exp(-LAMBDA*SOLVERTIME)*TUMOUR;

$CAPTURE	DOSE CL
'
# Compile the model code
  tum.mod <- mcode("popTUM",tum.code)

# ------------------------------------------------------------------------------
# Pull out variance parameters from model file
  tum.OMEGA <- as.matrix(omat(tum.mod))
  tum.SIGMA <- as.matrix(smat(tum.mod))
