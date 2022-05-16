* overall
global years 1989 1992 1995 1998 2001 2004 2007 2010 2013 

* a. income
global incvars X5729 X5702 // official definition
global taxvars year mstat depx agex pwages swages dividends otherprop pensions gssi transfers rentpaid otheritem childcare ui depchild mortgage stcg ltcg

* b. liquid assets
global checkvars X3506 X3510 X3514 X3518 X3522 X3526 X3529 // same as Telyuokova
global savvars_2001pre X3804 X3807 X3810 X3813 X3816 X3818 // same af Telyuokova
global savvars_2004post X3730 X3736 X3742 X3748 X3754 X3760 // close to Telyuokova (new variable names)
global eduactsvars_2001 X6456 X6457 X6458 X6459 X6460 // same af Telyuokova
global eduactsvars_2001post X3732 X3738 X3744 X3750 X3756 X3762 // close to Telyuokova (new variable names)
global brokactsvars X3930 X3932 // same af Telyuokova
	
	// composite (all variables needed to calculate row total)
	global savvars $savvars_2001pre $savvars_2004post
	global liqassvars $checkvars $savvars $brokactsvars
	global educacts $eduactsvars_2001 $eduactsvars_2001post
	
* c. credit card debt
global ccdebtvars X413 X421 // same as Telyuokova
global otherccvars X410 X432 // same as Telyuokova

* d. demographics
global demovars X8022 // same as Telyuokova

* e. installment loans
global storebalvars_1989 X424 X427 X430 // official defintion, not in ccdebtvars
global storebalvars_2007pre X424 X427 X430 X7575 // official defintion, not in ccdebtvars
global storebalvars_2010post X427 X430 X7575 // official defintion, not in ccdebtvars

	// composite (all variables needed to calculate row total)
	global storebalvars $storebalvars_2007pre

global veh_instvars_1992pre X2218 X2318 X2418 X2424 X2519 X2619 X2625 // official definition
global veh_instvars_1995post X2218 X2318 X2418 X2424 X2519 X2619 X2625 X7169 // official definition

	// composite (all variables needed to calculate row total)
	global veh_instvars $veh_instvars_1995post

global edn_instvars_1989 X2723 X2740 X2823 X2840 X2923 X2940 // official definition
global edn_instvars_1992 X7824 X7847 X7870 X7924 X7947 X7970 X2723 X2740 X2823 X2840 X2923 X2940 // official definition
global edn_instvars_1995post X7824 X7847 X7870 X7179 X7924 X7947 X7970 X2723 X2740 X2823 X2840 X2923 X2940 // official definition
global edn_instvars_i X2710 X2727 X2810 X2827 X2910 X2927 // official definition

	// composite (all variables needed to calculate row total)
	global edn_instvars $edn_instvars_1995post
