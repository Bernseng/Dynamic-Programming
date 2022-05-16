clear
set more off
cd "C:\Users\okoJD\Dropbox\Projects\2015_CreditCardPuzzle\Data\access\"

****************
* 1. LOAD DATA *
****************

local demovars P14 P19 X7020 P7020 P9002 P9003 X6809 X7004 X101 X5902 X6102 X5904 X6104 
local incvars X5729 X7650 X7362 P6781 P6785 X4106 X4706 P4106 P4706
local homevars X508 X601 X701 X714
local liquidvars X091001 X091008 X3929
local CCvars P410 X410 X09205 P09205 X09206 P09206 X414 P09207
use YY1 Y1 P42001 `liquidvars' `homevars' `CCvars' `demovars' `incvars' using p09pi6, clear

* renaming
rename YY1 id
rename Y1 imp
rename P42001 weight
rename X414 X09207

* specify mnumber
gen mnumber = imp-id*10


***********************
* 2. CREATE VARIABLES *
***********************

* a. demograhics
gen age = P14
gen age_sq = age*age
gen partner_age = P19

gen couple = X7020 == 2
gen stable_couple = inlist(P9002,1,2) & inlist(P9003,1,2)
gen stable_single = X7020 == 1 & P7020 == 1

gen minority = inlist(X6809,2,3,-7) | X7004 == 1
gen size = X101

gen educ = 0
replace educ = 1 if inlist(X5902,1,2) | inlist(X6102,1,2)
replace educ = 2 if X5904 == 1 | X6104 == 1
		
drop `demovars'


* b. income, unemployment and self-employed

gen income = .
replace income = X5729 if X7650 == 3
replace income = X7362 if X7650 != 3
gen log_income = log(income)

gen unemp = P6781+P6785 >= 1
gen some_unemp = P6781+P6785 >= 4
gen deep_unemp = P6781+P6785 >= 12

gen selfemp = X4106 == 2 | X4706 == 2 | P4106 ==2 | P4706 == 2

drop `incvars'


* c. homeownership
gen own = 0
replace own = 1 if inlist(X508,1)	
replace own = 1 if inlist(X601,1)
replace own = 1 if inlist(X701,1,3,4,5)

drop `homevars'


* d. liquid assets
gen checking = X091001
gen savings  = X091008
gen brokacts = X3929

gen liquid= checking + savings + brokacts 

drop `liquidvars'


*********************************
* 3. CREDIT CARD ACCESS AND USE *
*********************************

* 2007
gen noCC_07  = X410 == 5
gen gotCC_07 = X410 == 1 
gen useCC_07 = X09205 > 0
gen balCC_07 = X09206 > 0

* 2009
gen noCC_09  = P410 == 5
gen gotCC_09 = P410 == 1
gen useCC_09 = P09205 > 0
gen balCC_09 = P09206 > 0

* we focus on active credit card users with an outstanding balance
gen pop = gotCC_07 == 1 & useCC_07 == 1 & balCC_07 == 1
gen lostaccess = noCC_09 == 1


***********************
* 4. RESTRICT DATASET *
***********************

gen subpop = 1
replace subpop = 0 if income <= 0
replace subpop = 0 if !inrange(age,25,59) | (stable_couple == 1 & !inrange(partner_age,25,59))
replace subpop = 0 if stable_single == 0 & stable_couple == 0
replace subpop = 0 if pop != 1


*********
* 5. MI	*
*********

* full data set
save SCFpanel_mod_no_m_0, replace
local varlist lostaccess selfemp unemp some_unemp deep_unemp age age_sq size minority educ own log_income liquid weight stable_single subpop
keep `varlist' mnumber id

* m = 0 data set
foreach var in `varlist' {
	bysort id: egen i_`var' = sd(`var')
	replace i_`var' = 1 if i_`var' > 0
}
  
foreach var in `varlist' {
	replace `var' = . if i_`var' == 1
}
drop i_*
keep if mnumber == 1
replace mnumber = 0

save SCFpanel_mod_m_0, replace
append using SCFpanel_mod_no_m_0

mi import flong, m(mnumber) id(id) clear ///
	imputed(`varlist')	// 
	
mi passive: gen weight2 = weight
mi unregister weight2
mi svyset id [pweight = weight2]
mi register passive weight2
mi varying
mi describe

* save dataset for regressions
save SCFpanel_mod_MI, replace


****************
* 6. TABLE 5.2 *
****************

use SCFpanel_mod_no_m_0, clear

gen weight2 = weight/5
svyset id [pweight=weight2]

* a. unemployment
svy: tabulate unemp lostaccess if subpop==1, row
svy: tabulate unemp lostaccess if subpop==1, col
stop
* b. some unemployment
svy: tabulate some_unemp lostaccess if subpop==1, row
svy: tabulate some_unemp lostaccess if subpop==1, col

* c. deep unemployment
svy: tabulate deep_unemp lostaccess if subpop==1, row
svy: tabulate deep_unemp lostaccess if subpop==1, col


****************
* 7. TABLE 5.3 *
****************

use SCFpanel_mod_MI, clear

* a. unemployment
mi est, mcerr eform(or) post: svy, subpop(subpop): logit lostaccess i.unemp, difficult iter(200) technique(bhhh)
mi est, mcerr eform(or) post: svy, subpop(subpop): logit lostaccess i.unemp age age_sq size i.minority i.stable_single, difficult iter(200) technique(bhhh)
mi est, mcerr eform(or) post: svy, subpop(subpop): logit lostaccess i.unemp age age_sq size i.minority i.stable_single i.educ i.own log_income liquid i.selfemp, difficult iter(200) technique(bhhh)


** b. some unemployment
*mi est, mcerr eform(or) post: svy, subpop(subpop): logit lostaccess i.some_unemp, difficult iter(200) technique(bhhh)
*mi est, mcerr eform(or) post: svy, subpop(subpop): logit lostaccess i.some_unemp age age_sq size i.minority i.stable_single, difficult iter(200) technique(bhhh)
*mi est, mcerr eform(or) post: svy, subpop(subpop): logit lostaccess i.some_unemp age age_sq size i.minority i.stable_single i.educ i.own log_income liquid i.selfemp, difficult iter(200) technique(bhhh)


* c. deep unemployment
mi est, mcerr eform(or) post: svy, subpop(subpop): logit lostaccess i.deep_unemp, difficult iter(200) technique(bhhh)
mi est, mcerr eform(or) post: svy, subpop(subpop): logit lostaccess i.deep_unemp age age_sq size i.minority i.stable_single, difficult iter(200) technique(bhhh)
mi est, mcerr eform(or) post: svy, subpop(subpop): logit lostaccess i.deep_unemp age age_sq size i.minority i.stable_single i.educ i.own log_income liquid i.selfemp, difficult iter(200) technique(bhhh)
