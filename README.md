# CnuBproject
Feasibility study: hydrophones and CnuB

##################

run "run_flux... .sh" file
  - in which type (BU/TD/OPT) is defined, 
  - in which energy resolution (Eres0, Eres30, Eres80) is defined,
  - in which flux (F), or flux*E^alpha (FE) is defined.

  1) It will create the configuration file, based on template "config.m4"
  2) It will create a .txt file with the flux (or flux*E^alpha) using script "src\NeutrinoFlux.py"
  3) It will apply the defined Energy resolution & Will perform Monte Carlo using script "src\Statistics.py"
	p-value, and CL-value for corresponding N are saved in .txt file
	plots are saved in .root file

run "src\N_p_CL_test.py" to analyze the p-values/CL-values

##################

For 201712_results:

m = 0.1eV is used for which the dip is within the analyzed range: 5e20 - 5e22 eV

the whole procedure is performed 3 times:
  - RUN15: N from 0 to 300 in steps of 50
  - RUN16: N from 301 to 650 in steps of 50
  - RUN17: N from 651 to 1400 in steps of 100

##################

Future improvements:

run_fluxes....sh
  - ideally you want to loop over different options instead of creating all different run_fluxes files, but we couldn't get it working somehow...

src\Statistics.py
  - The energy range to analyze is specified in src\Statistics.py by lower bound (Elb), and upperbound (Eub).
	It has to be changed if you want to analyze other masses (not very handy I know, sorry)
  - The range should be specified in src\Statistics.py, in the first for-loop in the main (not very handy I know, sorry)
  - Now the RUN number should be changed by hand (2x) in src\Statistics.py (not very handy I know, sorry)
