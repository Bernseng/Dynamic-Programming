********************************************
* 1. DOWNLOAD VARIABLES AND INSTALL TAXSIM * 
********************************************

a. SCF data for SAS CPORT
 Source: http://www.federalreserve.gov/econresdata/scf/scfindex.htm
 Destination: data_sas_in\

b. SCF data with TAXSIM variable for STATA
 Source: http://users.nber.org/~taxsim/to-taxsim/scf/dta/ (all swt*.dta)
 Destination: data_stata\

c. Install TAXSIM for STATA
 Source: http://users.nber.org/~taxsim/to-taxsim/scf/
 Install by:
  net from "http://www.nber.org/stata"
  net describe taxsim9
  net install taxsim9


*******************
* 2. RUN PROGRAMS * 
*******************

NOTE 1: Paths must be adjusted in all programs.
NOTE 2: 1_run_SCFaggregations.sas is based on SAS-file at:
        http://www.federalreserve.gov/econresdata/scf/scfindex.htm

0_run_CIMPORT.sas:          convert CPORT files to SAS datasets
1_run_SCFaggregations.sas:  calculate official SCF aggregates
2_run_exportSTATA.sas:      convert SAS datasets to STATA datasets
3_run_select_aggregates.do: select the needed aggregates
4_run_TAXSIM.do:            calculate tax variables for all years (see VARLISTS.do)
5_run_createSCF.do:         create combined SCF across all years -> scf_all.dta
6_run_calculations.do:      perform various calculations on scf_all.dta
7_run_tables_figures.py:    produce tables and figures (Python 2.7)