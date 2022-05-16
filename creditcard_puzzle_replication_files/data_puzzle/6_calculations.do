clear
set more off
cd "C:\Users\okoJD\Dropbox\Projects\2015_CreditCardPuzzle\Data\puzzle"

local puzzlecutval_telyukova = 0.5
local puzzlecutvals 0.5 // cuttofs for puzzle and borrower groups, re-scaled to mean income in other years

run VARLISTS
use scf_all.dta, clear 

*************
* 1. INCOME *
*************

// a. all monetary variables: missing, = 0 and divide by 1000
global installvars $storebalvars $veh_instvars $edn_instvars
global externalvars fiitax networth houses homeeq
global monvars $incvars $liqassvars $ccdebtvars $installvars $externalvars
foreach var of global monvars  {
	replace `var' = 0 if missing(`var')
	replace `var' = (`var'>0)*`var'/1000
}
rename X5729 income
rename X5702 wincome


// b. disp = after tax
gen income_disp  = income-fiitax
gen wincome_disp = wincome-fiitax

// c. mean income in each year
gen mean_income_disp = .
gen mean_income_disp_2001 = .

foreach year of global years {

	quietly tabstat income_disp if year == `year' [aw=weight], stat(mean) save
	matrix list r(StatTotal)
	matrix stats=r(StatTotal)
	local mean_income_now = stats[1,1]
	
	local mean_income_qrt = `mean_income_now'/4.0
	display "`year': `mean_income_qrt'"
	
	replace mean_income_disp = `mean_income_now' if year == `year'
	if `year' == 2001 {
		replace mean_income_disp_2001 = `mean_income_now'
	}
		
}
stop

// d. conversion factor for cut-offs
gen fac = mean_income_disp/mean_income_disp_2001


********************
* 2. LIQUID ASSETS *
********************

* a. checking
egen check = rowtotal($checkvars)

* b. saving
foreach var of global educacts {
	replace `var'= int(`var')
}

* remove saving accounts related to education, money market accounts or health
replace X3804 = 0 if year == 2001 & X6456 == 3 // not availible earlier
replace X3807 = 0 if year == 2001 & X6457 == 3 // not availible earlier
replace X3810 = 0 if year == 2001 & X6458 == 3 // not availible earlier
replace X3813 = 0 if year == 2001 & X6459 == 3 // not availible earlier
replace X3816 = 0 if year == 2001 & X6460 == 3 // not availible earlier
replace X3730 = 0 if inlist(year,2004,2007,2010,2013) & !inlist(X3732,1,3,5,6,12,20,-7)
replace X3736 = 0 if inlist(year,2004,2007,2010,2013) & !inlist(X3738,1,3,5,6,12,20,-7)
replace X3742 = 0 if inlist(year,2004,2007,2010,2013) & !inlist(X3744,1,3,5,6,12,20,-7)
replace X3748 = 0 if inlist(year,2004,2007,2010,2013) & !inlist(X3750,1,3,5,6,12,20,-7)
replace X3754 = 0 if inlist(year,2004,2007,2010,2013) & !inlist(X3756,1,3,5,6,12,20,-7)
replace X3760 = 0 if inlist(year,2004,2007,2010,2013) & !inlist(X3762,1,3,5,6,12,20,-7)
egen sav = rowtotal($savvars)

* c. brokage account
egen brokacts = rowtotal($brokactsvars)

* d. total
egen liqass = rowtotal(check sav brokacts)


***********************
* 3. CREDIT CARD DEBT *
***********************

rename X432 revolve
rename X410 gotCC

egen ccdebt = rowtotal($ccdebtvars)


***********
* 4. MORE *
***********

* LIQUID NET WORTH
gen liqnetworth = liqass - ccdebt

* INSTALLMENT LOANS
replace X2723 = 0 if X2710 != 83
replace X2740 = 0 if X2727 != 83
replace X2823 = 0 if X2810 != 83
replace X2840 = 0 if X2827 != 83
replace X2923 = 0 if X2910 != 83
replace X2940 = 0 if X2927 != 83
egen install = rowtotal($installvars)


****************************************
* 5. DEMOGRAPHICS AND SAMPLE SELECTION *
****************************************

rename X8022 age 

gen pop = 1
//replace pop = 0 if wincome_disp <= 0 | missing(wincome_disp)
replace pop = 0 if income_disp  <= 0 | missing(income_disp)
replace pop = 0 if !inrange(age,25,64)

drop X*
drop if pop != 1


**********************
* 6. AS IN TELYUKOVA *
**********************

// a. puzzle definition
gen int puzzle_telyukova =.
replace puzzle_telyukova = 1 if ccdebt > `puzzlecutval_telyukova'*fac & liqass > `puzzlecutval_telyukova'*fac & revolve > 1
replace puzzle_telyukova = 2 if ccdebt > `puzzlecutval_telyukova'*fac & puzzle_telyukova == . & revolve > 1
replace puzzle_telyukova = 3 if puzzle_telyukova == .

cap label define puzzle_telyukova 1 "Puzzle" 2 "Borrower" 3 "Saver"
label value puzzle_telyukova puzzle_telyukova

* b. shares
gen puzzle_share_telyukova = (puzzle_telyukova == 1)*100
gen borrower_share_telyukova = (puzzle_telyukova == 2)*100
gen saver_share_telyukova = (puzzle_telyukova == 3)*100

tabstat puzzle_share_telyukova borrower_share_telyukova saver_share_telyukova [aw=weight], stat(mean) by(year) format(%5.1f) 

* c. level comparisons with Telyukova
tabstat income_disp ccdebt liqass liqnetworth [aw=weight] if year == 2001, by(puzzle_telyukova) stat(mean p50) format(%7.3f) 


********************
* 7. AS IN FULFORD *
********************

// a. puzzle definition
gen sav_fulford = sav
replace sav_fulford = sav_fulford + (check-income/12) if (check-income/12) > 0

gen int puzzle_fulford =.
replace puzzle_fulford = 1 if ccdebt > 0.0001*income & sav_fulford > 0.0001*income & gotCC == 1
replace puzzle_fulford = 2 if ccdebt > 0.0001*income & puzzle_fulford == . & gotCC == 1
replace puzzle_fulford = 3 if sav_fulford > 0.0001*income & puzzle_fulford == . & gotCC == 1
replace puzzle_fulford = 4 if puzzle_fulford == . & gotCC == 1

cap label define puzzle_fulford 1 "Puzzle" 2 "Borrower" 3 "Saver" 4 "Corner"
label value puzzle_fulford puzzle_fulford

* b. shares (old definition)	
gen puzzle_share_fulford = (puzzle_fulford == 1)*100
gen borrower_share_fulford = (puzzle_fulford == 2)*100
gen saver_share_fulford = (puzzle_fulford == 3)*100
gen corner_share_fulford = (puzzle_fulford == 4)*100

tabstat puzzle_share_fulford borrower_share_fulford saver_share_fulford corner_share_fulford [aw=weight] if gotCC == 1, stat(mean) by(year) format(%5.1f) 


**********
* 8. NEW *
**********

foreach puzzlecutval of local puzzlecutvals {
	
	local cutval = `puzzlecutval'*1000
	
	* a. puzzle definition
	cap drop puzzle
	gen int puzzle =.
	replace puzzle = 1 if ccdebt > `puzzlecutval'*fac & liqass > `puzzlecutval'*fac & revolve > 1
	replace puzzle = 2 if ccdebt - liqass > `puzzlecutval'*fac & puzzle == . & revolve > 1
	replace puzzle = 3 if liqass - ccdebt > `puzzlecutval'*fac & puzzle == .
	replace puzzle = 4 if puzzle == .

	cap drop puzzle_alt
	gen int puzzle_alt =.
	replace puzzle_alt = 1 if ccdebt > `puzzlecutval'*fac & liqass > `puzzlecutval'*fac
	replace puzzle_alt = 2 if ccdebt - liqass > `puzzlecutval'*fac & puzzle_alt == .
	replace puzzle_alt = 3 if liqass - ccdebt > `puzzlecutval'*fac & puzzle_alt == .
	replace puzzle_alt = 4 if puzzle_alt == .

	cap label define puzzle 1 "Puzzle" 2 "Borrower" 3 "Saver" 4 "Corner"
	label value puzzle puzzle
	label value puzzle_alt puzzle
	
	cap drop liqass_above_income
	gen liqass_above_income = liqass > income_disp/12
	tab liqass_above_income puzzle [aw=weight], nofreq col

	cap drop liqass_above_income
	gen liqass_above_income = liqass > 0.5*income_disp/12
	tab liqass_above_income puzzle [aw=weight], nofreq col
	
	cap drop homeown
	gen homeown = houses > 0
	tab homeown puzzle [aw=weight], nofreq col
	
	gen puzzle_share = (puzzle == 1)*100
	gen borrower_share = (puzzle == 2)*100
	gen saver_share = (puzzle == 3)*100
	gen corner_share = (puzzle == 4)*100

	tabstat puzzle_share borrower_share saver_share corner_share [aw=weight], stat(mean) by(year) format(%5.1f) 

	gen puzzle_alt_share = (puzzle_alt == 1)*100
	gen borrower_alt_share = (puzzle_alt == 2)*100
	gen saver_alt_share = (puzzle_alt == 3)*100
	gen corner_alt_share = (puzzle_alt == 4)*100
	
	tabstat puzzle_alt_share borrower_alt_share saver_alt_share corner_alt_share [aw=weight], stat(mean) by(year) format(%5.1f) 
	
	* b. relative to quarterly income after tax
	gen ccdebt_rel = ccdebt/mean_income_disp*100*4
	gen liqass_rel = liqass/mean_income_disp*100*4
	gen liqnetworth_rel = liqnetworth/mean_income_disp*100*4

	gen income_disp_rel = income_disp/mean_income_disp*100
	gen houses_rel = houses/mean_income_disp*100*4
	gen homeeq_rel = homeeq/mean_income_disp*100*4
	gen networth_rel = networth/mean_income_disp*100*4
	gen install_rel = install/mean_income_disp*100*4
		
	* c. central tables
	local outvars ccdebt liqass liqnetworth income_disp houses homeeq networth install
	foreach var of local outvars {
		
		tabstat `var'_rel [aw=weight], stat(mean med) by(year) format(%5.1f) 
		tabstat `var'_rel [aw=weight], stat(mean med) by(puzzle) format(%5.1f) 
		
		preserve
		collapse (mean) mean = `var'_rel ///
						(p5) p5 = `var'_rel (p15) p15 = `var'_rel (p25) p25 = `var'_rel (p50) p50 = `var'_rel ///
						(p75) p75 = `var'_rel (p85) p85 = `var'_rel (p95) p95 = `var'_rel (p99) p99 = `var'_rel ///
						[aw=weight], by(puzzle)
		export excel using "output_`cutval'.xls", sheet("`var'") sheetreplace firstrow(variables)
		restore 

		preserve
		collapse (mean) mean =`var'_rel ///
						(p5) p5 = `var'_rel (p15) p15 = `var'_rel (p25) p25 = `var'_rel (p50) p50 = `var'_rel ///
						(p75) p75 = `var'_rel (p85) p85 = `var'_rel (p95) p95 = `var'_rel (p99) p99 = `var'_rel ///
						[aw=weight]
		export excel using "output_`cutval'.xls", sheet("`var'_tot") sheetreplace firstrow(variables)
		restore 
		
		preserve
		collapse (mean) mean = `var'_rel ///
						(p5) p5 = `var'_rel (p15) p15 = `var'_rel (p25) p25 = `var'_rel (p50) p50 = `var'_rel ///
						(p75) p75 = `var'_rel (p85) p85 = `var'_rel (p95) p95 = `var'_rel (p99) p99 = `var'_rel ///
						[aw=weight], by(year)
		export excel using "output_`cutval'.xls", sheet("`var'_year") sheetreplace firstrow(variables)
		restore 
		
	}


	// d. shares

		preserve
		collapse (mean) puzzle_share borrower_share saver_share corner_share [aw=weight], by(year)
		export excel using "output_`cutval'.xls", sheet("shares_year") sheetreplace firstrow(variables)
		restore 
		
		preserve
		collapse (mean) puzzle_share borrower_share saver_share corner_share [aw=weight]
		export excel using "output_`cutval'.xls", sheet("shares") sheetreplace firstrow(variables)
		restore 
		
	tabstat liqass_rel ccdebt_rel puzzle_share borrower_share saver_share corner_share [aw=weight], stat(mean) by(year) format(%5.1f) 
	
	drop *_rel *_share
}
