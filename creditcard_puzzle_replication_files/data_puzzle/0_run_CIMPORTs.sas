LIBNAME IN 'C:\Users\okoJD\Documents\SCF\data_sas_in\';

proc cimport data=IN.P13I6 infile='C:\Users\okoJD\Documents\SCF\data_sas_raw\scf2013';
run;

proc cimport data=IN.P10I6 infile='C:\Users\okoJD\Documents\SCF\data_sas_raw\scf2010';
run;

proc cimport data=IN.P07I6 infile='C:\Users\okoJD\Documents\SCF\data_sas_raw\scf2007';
run;

proc cimport data=IN.P04I6 infile='C:\Users\okoJD\Documents\SCF\data_sas_raw\scf2004';
run;

proc cimport data=IN.P01I6 infile='C:\Users\okoJD\Documents\SCF\data_sas_raw\scf2001';
run;

proc cimport data=IN.P98I6 infile='C:\Users\okoJD\Documents\SCF\data_sas_raw\scf98';
run;

proc cimport data=IN.P95I6 infile='C:\Users\okoJD\Documents\SCF\data_sas_raw\scf95';
run;

proc cimport data=IN.P92I4 infile='C:\Users\okoJD\Documents\SCF\data_sas_raw\scf92';
run;

proc cimport data=IN.P89I6 infile='C:\Users\okoJD\Documents\SCF\data_sas_raw\scf89';
run;
