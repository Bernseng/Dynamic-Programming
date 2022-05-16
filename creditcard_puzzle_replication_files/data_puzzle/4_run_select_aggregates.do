clear
set more off
cd "C:\Users\okoJD\Dropbox\Projects\2015_CreditCardPuzzle\Data\puzzle"

****************
* 1. VARIABLES *
****************

run VARLISTS


***********
* 2. MAIN *
***********

foreach year of global years {

	if `year' == 1989 {
		use * using "data_stata\scf`year'_aggregates", clear
		rename xx1 YY1
		rename x1 Y1
	}
	else {
		use * using "data_stata\scf`year'_aggregates", clear
		rename yy1 YY1
		rename y1 Y1
	}
	keep YY1 Y1 networth homeeq houses
	save "data_stata\scf`year'_aggregates_selected", replace

}
