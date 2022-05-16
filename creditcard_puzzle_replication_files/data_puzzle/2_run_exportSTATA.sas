LIBNAME WORKINGN 'C:\Users\okoJD\Documents\SCF\data_out_nom\';

proc export data=WORKINGN.Scfp1989 outfile="C:\Users\okoJD\Documents\SCF\data_stata\scf1989_aggregates.dta" DBMS=Stata REPLACE;
run;

proc export data=WORKINGN.Scfp1992 outfile="C:\Users\okoJD\Documents\SCF\data_stata\scf1992_aggregates.dta" DBMS=Stata REPLACE;
run;

proc export data=WORKINGN.Scfp1995 outfile="C:\Users\okoJD\Documents\SCF\data_stata\scf1995_aggregates.dta" DBMS=Stata REPLACE;
run;

proc export data=WORKINGN.Scfp1998 outfile="C:\Users\okoJD\Documents\SCF\data_stata\scf1998_aggregates.dta" DBMS=Stata REPLACE;
run;

proc export data=WORKINGN.Scfp2001 outfile="C:\Users\okoJD\Documents\SCF\data_stata\scf2001_aggregates.dta" DBMS=Stata REPLACE;
run;

proc export data=WORKINGN.Scfp2004 outfile="C:\Users\okoJD\Documents\SCF\data_stata\scf2004_aggregates.dta" DBMS=Stata REPLACE;
run;

proc export data=WORKINGN.Scfp2007 outfile="C:\Users\okoJD\Documents\SCF\data_stata\scf2007_aggregates.dta" DBMS=Stata REPLACE;
run;

proc export data=WORKINGN.Scfp2010 outfile="C:\Users\okoJD\Documents\SCF\data_stata\scf2010_aggregates.dta" DBMS=Stata REPLACE;
run;

proc export data=WORKINGN.Scfp2013 outfile="C:\Users\okoJD\Documents\SCF\data_stata\scf2013_aggregates.dta" DBMS=Stata REPLACE;
run;
