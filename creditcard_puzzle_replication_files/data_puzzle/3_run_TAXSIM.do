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
		use XX1 X1 $taxvars using "data_stata\swt`year'", clear
		rename XX1 YY1
		rename X1 Y1
	}
	else {
		use YY1 Y1 $taxvars using "data_stata\swt`year'", clear
	}
	taxsim9, replace	
	save "data_stata\scf`year'_tax_all", replace
	keep YY1 Y1 fiitax
	save "data_stata\scf`year'_tax_all", replace
	
}
