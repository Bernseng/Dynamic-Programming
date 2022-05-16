clear
set more off
cd "C:\Users\okoJD\Dropbox\Projects\2015_CreditCardPuzzle\Data\puzzle"

****************
* 1. VARIABLES *
****************

run VARLISTS

****************
* 2. LOAD DATA *
****************

foreach year of global years {
	
	local loadvars $incvars $checkvars $brokactsvars $ccdebtvars $otherccvars $demovars $edn_instvars_i
	
	* a. determine load variables
	if `year' == 1989 {
		local loadvars `loadvars' $savvars_2001pre $storebalvars_1989 $veh_instvars_1992pre $edn_instvars_1989
	}
	else if `year' == 1992 {
		local loadvars `loadvars' $savvars_2001pre $storebalvars_2007pre $veh_instvars_1992pre $edn_instvars_1992
	}
	else if `year' == 1995 {
		local loadvars `loadvars' $savvars_2001pre $storebalvars_2007pre $veh_instvars_1995post $edn_instvars_1995post
	}
	else if `year' == 1998 {
		local loadvars `loadvars' $savvars_2001pre $storebalvars_2007pre $veh_instvars_1995post $edn_instvars_1995post
	}
	else if `year' == 2001 {
		local loadvars `loadvars' $savvars_2001pre $eduactsvars_2001 $storebalvars_2007pre $veh_instvars_1995post $edn_instvars_1995post
	}
	else if `year' == 2004 {
		local loadvars `loadvars' $savvars_2004post $eduactsvars_2001post $storebalvars_2007pre $veh_instvars_1995post $edn_instvars_1995post
	}
	else if `year' == 2007 {
		local loadvars `loadvars' $savvars_2004post $eduactsvars_2001post $storebalvars_2007pre $veh_instvars_1995post $edn_instvars_1995post 
	}
	else if `year' == 2010 {
		local loadvars `loadvars' $savvars_2004post $eduactsvars_2001post $storebalvars_2010post $veh_instvars_1995post $edn_instvars_1995post
	}
	else if `year' == 2013 {
		local loadvars `loadvars' $savvars_2004post $eduactsvars_2001post $storebalvars_2010post $veh_instvars_1995post $edn_instvars_1995post
	}
	
	* b. load variables and merge in tax and aggregates
	if `year' == 1989 {
		use XX1 X1 X42001 `loadvars' using data_stata\swt`year', clear
		rename XX1 YY1
		rename X1 Y1
	}
	else {
		use YY1 Y1 X42001 `loadvars' using data_stata\swt`year', clear
	}
	merge 1:1 YY1 Y1 using data_stata\scf`year'_tax
	drop _merge
	merge 1:1 YY1 Y1 using data_stata\scf`year'_aggregates_selected
	drop _merge
	gen year =`year'
	
	* c. create variables not loaded
	if `year' == 1989 {
		local misvars $savvars_2004post $eduactsvars_2001 $eduactsvars_2001post	X7575 X7169 X7179 X7824 X7847 X7870 X7924 X7947 X7970
	}
	else if `year' == 1992 {
		local misvars $savvars_2004post $eduactsvars_2001 $eduactsvars_2001post	X7169 X7179
	}
	else if `year' == 1995 {
		local misvars $savvars_2004post $eduactsvars_2001 $eduactsvars_2001post	
	}
	else if `year' == 1998 {
		local misvars $savvars_2004post $eduactsvars_2001 $eduactsvars_2001post	
	}
	else if `year' == 2001 {
		local misvars $savvars_2004post $eduactsvars_2001post
	}
	else if `year' == 2004 {
		local misvars $savvars_2001pre $eduactsvars_2001
	}
	else if `year' == 2007 {
		local misvars $savvars_2001pre $eduactsvars_2001
	}
	else if `year' == 2010 {
		local misvars $savvars_2001pre $eduactsvars_2001 X424
	}
	else if `year' == 2013 {
		local misvars $savvars_2001pre $eduactsvars_2001 X424
	}
	foreach var of local misvars {
		gen `var' = .
	}
		
	* d. save yearly file
	save "data_stata\scf`year'", replace
	
}


***********************************
* 3. APPEND YEARLY FILES TOGETHER *
***********************************

clear all
foreach year of global years {
	append using "data_stata\scf`year'"
}
rename X42001 weight
save scf_all.dta, replace
