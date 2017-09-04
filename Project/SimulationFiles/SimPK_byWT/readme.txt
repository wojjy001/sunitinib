#################################
### Directories to be changed ###
#################################
1. Open "run_pk.R"
2. Line 8; set "file.dir" to the directory where the “SimPK_byWT” folder is
saved (not including "SimPK_byWT" in the name)
3. Line 16; set "output.dir" to the same directory.  All simulated output will
be saved into the folder where the .R files are saved.
4. Line 37; "study.name" can only be one of the studies listed in the directory
Each of these "studies" have their own .R script in "SimPK_byWT".  Changing
"study.name" allows "run_pk.R" to source the correct study script
5. Open "run_pd.R"
6. Line 12; change "output.dir" to the same as Step 3
7. Line 13; change "study.name" to the same as Step 4
8. Line 21; change "file.dir" to the same as Step 2

###########################
### Simulation Sequence ###
###########################
1. Need to run "run_pk.R" first
2. Run "run_pd.R" second. You can run the PD simulation at a separate instance
if you had previously run and saved the PK simulation
