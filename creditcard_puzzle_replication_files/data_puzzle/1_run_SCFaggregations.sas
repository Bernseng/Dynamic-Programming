* this macro is designed to be used with SCF datasets from 1989
  forward to produce the core variables used in the current and
  historical Bulletin tables;

OPTIONS MPRINT;
OPTIONS SOURCE SOURCE2;
OPTIONS LINESIZE=78;
OPTIONS COMPRESS=YES;

LIBNAME IN 'C:\Users\okoJD\Documents\SCF\data_sas_in\';
LIBNAME WORKINGR 'C:\Users\okoJD\Documents\SCF\data_sas_out_real\';
LIBNAME WORKINGN 'C:\Users\okoJD\Documents\SCF\data_sas_out_nom\';

* IMPORTANT NOTE: This program assumes the SCF datasets stored in
  the IN library are named P89I6, P92I4, P95I6, P98I6, P01I6, P04I6,
  P07I6, P10I6, and P13I6.  If you wish to use a different naming convention, you will
  need to change the ELSE %LET SCFDS= statements in the section of the
  program where the CPI adjustments are done.;

***************************************************************************;
***************************************************************************;
* UTILITY MACROs HERE;
***************************************************************************;

* MACRO MCONV converts payments to MONTHLY basis;
* NOTE: this macro does not convert HOURLY (code 18) or DAILY (code 1);

  %MACRO MCONV(F);
    ( (&F=2)*52/12+(&F=3)*26/12+(&F=4)+(&F=5)/3+(&F=6)/12+(&F=11)/6+
    (&F=12)/2+(&F=31)*2+(&F=23)*13/12+(&F=24)*52/(6*12) )
  %MEND MCONV;

***************************************************************************;
***************************************************************************;
* MAIN MACRO HERE;
***************************************************************************;

* MACRO BULLIT:
  YEAR: SCF survey year
  REAL: set=YES (default) to compute dollar variables in real terms
  ADJINC: set=YES (default) to adjust lagged income to survey year dollars
  CPIBASE: give value of CPI-U-RS (times 100) for the period selected
  as the reference period for real dollars
  PUBLIC: set equal to NO (default) for work with an internal SCF
  dataset, set equal to YES for work with a disclosure-limited dataset;

%GLOBAL YEAR ID IID;
%MACRO BULLIT(YEAR=,REAL=YES,CPIBASE=,PUBLIC=NO,ADJINC=YES);

***************************************************************************;
***************************************************************************;
  %IF (&REAL EQ YES AND &ADJINC NE YES) %THEN %PUT WARNING! REAL=YES,
    BUT ADJINC^=YES!;
***************************************************************************;
***************************************************************************;

%*set CPI-U-RS factors for inflation adjustments and identify input
  datasets for various survey years;
%*CPI-U-RS is the current-methods version of the Consumer Price Index
  for all urban consumers;
%*for income, set CPILAG equal to the annual average of CPI-U-RS * 100
  for the year of the survey divided by that for the year preceding
  the survey;
%*for all other values, set CPICUR equal to CPIBASE divided by the
  value of CPI-U-RS * 100 for the September of the survey year;
%*NOTE: CPIBASE should be the index value for the base period
  (typically September for one of the survey years);
%* see http://www.bls.gov/cpi/cpirsai1978-2013.pdf
   NOTE: the file at the above link includes revised values of earlier
    estimates (accounting for a new method of imputing price change from
    rental vacancies);

%*set name of final input dataset;
  %IF (&YEAR=1989) %THEN %DO;
    %LET CPIADJ=&CPIBASE/1902;
    %LET CPILAG=1886/1808;
    %IF (&PUBLIC NE YES) %THEN %LET SCFDS=SCFR6;
    %ELSE %LET SCFDS=P89I6;
  %END;
  %ELSE %IF (&YEAR=1992) %THEN %DO;
    %LET CPIADJ=&CPIBASE/2116;
    %LET CPILAG=2103/2051;
    %IF (&PUBLIC NE YES) %THEN %LET SCFDS=SCF92I4;
    %ELSE %LET SCFDS=P92I4;
  %END;
  %ELSE %IF (&YEAR=1995) %THEN %DO;
    %LET CPIADJ=&CPIBASE/2265;
    %LET CPILAG=2254/2201;
    %IF (&PUBLIC NE YES) %THEN %LET SCFDS=SCF95I6;
    %ELSE %LET SCFDS=P95I6;
  %END;
  %ELSE %IF (&YEAR=1998) %THEN %DO;
    %LET CPIADJ=&CPIBASE/2405;
    %LET CPILAG=2397/2364;
    %IF (&PUBLIC NE YES) %THEN %LET SCFDS=SCF98I6;
    %ELSE %LET SCFDS=P98I6;
  %END;
  %ELSE %IF (&YEAR=2001) %THEN %DO;
    %LET CPIADJ=&CPIBASE/2618;
    %LET CPILAG=2600/2529;
    %IF (&PUBLIC NE YES) %THEN %LET SCFDS=SCF01I6;
    %ELSE %LET SCFDS=P01I6;
  %END;
  %ELSE %IF (&YEAR=2004) %THEN %DO;
    %LET CPIADJ=&CPIBASE/2788;
    %LET CPILAG=2774/2701;
    %IF (&PUBLIC NE YES) %THEN %LET SCFDS=SCF04I6;
    %ELSE %LET SCFDS=P04I6;
  %END;
  %ELSE %IF (&YEAR=2007) %THEN %DO;
    %LET CPIADJ=&CPIBASE/3062;
    %LET CPILAG=3045/2961;
    %IF (&PUBLIC NE YES) %THEN %LET SCFDS=SCF07I6;
    %ELSE %LET SCFDS=P07I6;
  %END;
  %ELSE %IF (&YEAR=2010) %THEN %DO;
    %LET CPIADJ=&CPIBASE/3208;
    %LET CPILAG=3202/3150;
    %IF (&PUBLIC NE YES) %THEN %LET SCFDS=SCF10I6;
    %ELSE %LET SCFDS=P10I6;
  %END;
  %ELSE %IF (&YEAR=2013) %THEN %DO;
    %LET CPIADJ=&CPIBASE/3438;
    %LET CPILAG=3421/3372;
    %IF (&PUBLIC NE YES) %THEN %LET SCFDS=SCF13I6;
    %ELSE %LET SCFDS=P13I6;
  %END;


***************************************************************************;
***************************************************************************;
* set names of the ID variables;
  %IF (&PUBLIC NE YES OR &YEAR EQ 1989) %THEN %DO;
    %LET ID=X1;
    %LET IID=XX1;
  %END;
  %ELSE %DO;
    %LET ID=Y1;
    %LET IID=YY1;
  %END;
  %IF (&PUBLIC NE YES OR &YEAR EQ 1989) %THEN %DO;
    %LET PID=X1;
    %LET PIID=XX1;
  %END;
  %ELSE %DO;
    %LET PID=Y1;
    %LET PIID=YY1;
  %END;
***************************************************************************;
***************************************************************************;
* specify list of variables to be kept in the dataset;
* to reduce the size of the public dataset, specify a more limited
  list of variables to be kept;
***************************************************************************;
%* specify here variables that are only to be kept for selected years;
  %LET ALSOKEEP=%STR();
  %IF (&YEAR GE 1992 AND &PUBLIC NE YES) %THEN %DO;
    %LET ALSOKEEP=%STR(Y1 YY1 WSAVED SAVED);
  %END;
  %ELSE %IF (&YEAR GE 1992) %THEN %DO;
    %LET ALSOKEEP=%STR(Y1 YY1 SAVED);
  %END;
  %ELSE %IF (&YEAR EQ 1989) %THEN %DO;
    %LET ALSOKEEP=%STR(X1 XX1);
  %END;
  %IF (&YEAR GE 1998 AND &PUBLIC NE YES) %THEN %DO;
    %LET ALSOKEEP=&ALSOKEEP SPENDMOR;
  %END;
  %IF (&YEAR GE 2004 AND &PUBLIC NE YES) %THEN %DO;
    %LET ALSOKEEP=&ALSOKEEP OMUTF;
  %END;
  %IF (&YEAR GE 1998) %THEN %DO;
    %LET ALSOKEEP=%STR(&ALSOKEEP BCALL BMAGZNEWS BMAILADTV BINTERNET
                       BFRIENDWORK BFINPRO BFINPLAN BSELF BDONT BOTHER 
                       ICALL IMAGZNEWS IMAILADTV IINTERNET IFRIENDWORK 
                       IFINPRO IFINPLAN ISELF IDONT IOTHER);                     
  %END;
  %IF (&YEAR GE 1995) %THEN %DO;
    %LET ALSOKEEP=%STR(&ALSOKEEP INTERNET BSHOPNONE BSHOPGRDL 
                       BSHOPMODR ISHOPNONE ISHOPGRDL ISHOPMODR
	           REFIN_EVER HEXTRACT_EVER PURCH1 );
  %END;
  %IF (&YEAR GE 2007 AND &PUBLIC NE YES) %THEN %DO;
    %LET ALSOKEEP=%STR(&ALSOKEEP ALT_XX1);
  %END;

* special keeplist for SDA, also added to public and internal keeplists;
  %LET SDAKEEP=%STR(PAYMORT1 PAYMORT2 PAYMORT3 PAYMORTO PAYLOC1
    PAYLOC2 PAYLOC3 PAYLOCO
    PAYHI1 PAYHI2 PAYLC1 PAYLC2 PAYLCO PAYORE1 PAYORE2 PAYORE3
    PAYOREV PAYVEH1 PAYVEH2 PAYVEH3 PAYVEH4 PAYVEHM PAYVEO1 PAYVEO2
    PAYVEOM PAYEDU1 PAYEDU2 PAYEDU3 PAYEDU4 PAYEDU5 PAYEDU6 PAYEDU7
    PAYILN1 PAYILN2 PAYILN3 PAYILN4 PAYILN5 PAYILN6 PAYILN7 PAYMARG
    PAYINS PAYPEN1 PAYPEN2 PAYPEN3 PAYPEN4 PAYPEN5 PAYPEN6 FARMBUS_KG MORT1
    MORT2 MORT3 REVPAY MORTPAY CONSPAY PIRMORT PIRCONS PIRREV
    PENACCTWD MMDA MMMF STMUTF TFBMUTF GBMUTF OBMUTF COMUTF 
    OMUTF HSTOCKS NSTOCKS WILSH NOTXBND MORTBND GOVTBND OBND 
    FUTPEN CURRPEN DBPLANCJ DBPLANT DCPLANCJ BPLANCJ ANNUIT TRUSTS 
    EQUITINC HBROK HTRAD NTRAD OWN NOWN LEASE NLEASE 
    NVEHIC NEWCAR1 NEWCAR2 FARMBUS NHNFIN LIFECL LF INDCAT 
    WSAVED SPENDMOR SPENDLESS EXPENSHILO EHCHKG IRAKH THRIFT HHSEX
    FOODHOME FOODAWAY FOODDELV PIR40 LEVRATIO DEBT2INC TURNDOWN
    FEARDENIAL TURNFEAR LATE HPAYDAY BNKRUPLAST5 NOCCBAL);

* remember to add dollar value variables to CPI adjustment array later
  in the program;
  %IF (&PUBLIC EQ YES) %THEN %LET KEEPLIST=&ALSOKEEP
    WGT MARRIED KIDS AGE INCOME NORMINC WAGEINC BUSSEFARMINC 
    INTDIVINC KGINC SSRETINC TRANSFOTHINC AGECL EDCL RACECL
    RACE HOUSECL ASSET FIN LIQ CDS SAVING NMMF STOCKS BOND RETQLIQ
    SAVBND CASHLI OTHMA OTHFIN CHECKING MMA CALL NFIN VEHIC HOUSES
    ORESRE NNRESRE BUS OTHNFIN DEBT PIRTOTAL MRTHEL RESDBT OTHLOC CCBAL
    INSTALL ODEBT NETWORTH HOMEEQ TPAY PLOAN1-PLOAN8 LLOAN1-LLOAN12
    EDUC OCCAT1 KGTOTAL SAVRES1-SAVRES9 ANYPEN BUSVEH
    NBUSVEH LATE60 EQUITY DEQ VLEASE RETEQ MERGEID HLIQ
    VEH_INST EDN_INST OTH_INST HELOC NH_MORT NOCHK WHYNOCKG HDEBT
    KGHOUSE KGORE KGBUS KGSTMF 
    HPRIM_MORT HSEC_MORT PURCH2 HMORT2 HELOC_YN ACTBUS NONACTBUS
    DONTWRIT MINBAL DONTLIKE SVCCHG NOMONEY CANTMANG CREDIT 
    DONTWANT OTHER CKLOCATION CKLOWFEEBAL CKMANYSVCS CKRECOMFRND
    CKPERSONAL CKCONNECTN CKLONGTIME CKSAFETY CKCONVPAYRL CKOTHCHOOSE
    FAMSTRUCT OCCAT2 HBUS &ALSOKEEP HCHECK HSAVING HMMA HCALL      
    HCDS HNMMF HBOND HRETQLIQ HSAVBND HCASHLI HOTHMA HOTHFIN HEQUITY    
    HFIN HVEHIC HHOUSES HORESRE HNNRESRE HOTHNFIN HNFIN HASSET HMRTHEL    
    HHELOC HNH_MORT HOTHLOC HRESDBT HCCBAL HVEH_INST HEDN_INST HOTH_INST  
    HINSTALL HODEBT HSTOCKS &ALSOKEEP &SDAKEEP    ; 

  %ELSE %LET KEEPLIST=&ALSOKEEP WGT0 WGT X1 XX1 WILSH MARRIED KIDS
    AGE INCCL2 INCOME NORMINC WAGEINC BUSSEFARMINC INTDIVINC KGINC 
    SSRETINC TRANSFOTHINC RENTINC X3 AGECL EDCL LIFECL FAMSTRUCT RACECL h_racecl
    RACE h_race HOUSECL ASSET FIN LIQ CDS
    NMMF STOCKS BOND RETQLIQ SAVBND CASHLI OTHMA OTHFIN CHECKING
    MMA CALL HFIN HLIQ HCDS SAVING HSAVING HNMMF HSTOCKS HBOND
    HRETQLIQ HSAVBND HCASHLI HOTHMA HOTHFIN HCHECK HMMA HCALL NFIN
    VEHIC HOUSES ORESRE NNRESRE BUS OTHNFIN HNFIN HVEHIC HHOUSES
    HORESRE HNNRESRE HBUS HOTHNFIN DEBT MRTHEL RESDBT OTHLOC CCBAL
    INSTALL ODEBT HDEBT HMRTHEL 
    HPRIM_MORT HSEC_MORT PURCH2 HMORT2
    HELOC_YN HRESDBT HOTHLOC HCCBAL HINSTALL
    HODEBT HASSET NETWORTH HOMEEQ TPAY PIRTOTAL TPLOAN PLOAN1-PLOAN8
    PLOANB1-PLOANB8
    TLLOAN LLOAN1-LLOAN12 EDUC X7401 X7402 OCCAT1 OCCAT2 INDCAT NHNFIN
    KGTOTAL KGHOUSE KGORE KGBUS KGSTMF SAVRES1-SAVRES9 ANYPEN DBPLANCJ
    DBPLANT DCPLANCJ BPLANCJ BUSVEH NBUSVEH LATE60 EQUITY EQUITINC
    HEQUITY DEQ OWN NOWN VOWN LEASE NLEASE VLEASE NVEHIC NEWCAR1 NEWCAR2
    NOCHK EHCHKG X3503 MERGEID REGION HBROK HTRAD NSTOCKS NTRAD
    RETEQ X3006 X30009 NOTXBND MORTBND GOVTBND OBND STMUTF TFBMUTF GBMUTF
    OBMUTF COMUTF IRAKH PENEQ THRIFT CURRPEN FUTPEN MORTPAY CONSPAY REVPAY 
    PIRMORT PIRCONS PIRREV WHYNOCKG VEH_INST HVEH_INST EDN_INST
    HEDN_INST OTH_INST HOTH_INST HELOC HHELOC NH_MORT
    HNH_MORT ANNUIT TRUSTS URBAN ACTBUS NONACTBUS 
    DONTWRIT MINBAL DONTLIKE SVCCHG NOMONEY CANTMANG CREDIT 
    DONTWANT OTHER CKLOCATION CKLOWFEEBAL CKMANYSVCS CKRECOMFRND
    CKPERSONAL CKCONNECTN CKLONGTIME CKSAFETY CKCONVPAYRL CKOTHCHOOSE 
    OUTMARG OUTPEN HHSEX SAVEQ &ALSOKEEP &SDAKEEP;

* NOTE: variable NWCAT is computed following this data step, and the
  value of this variable is merged back into the final dataset;
* MERGEID is dropped at the end;

***************************************************************************;

  DATA %UNQUOTE(SCFB&YEAR)(KEEP=&KEEPLIST);

    LENGTH VARNM $8;

    %IF (&PUBLIC EQ NO) %THEN %DO;
      SET %UNQUOTE(SCF&YEAR)%STR(.)&SCFDS;
    %END;
    %ELSE %IF (&PUBLIC EQ YES) %THEN %DO;
      SET IN%STR(.)&SCFDS;
    %END;


***************************************************************************;
*   sample and weight adjustments;

*   only keep observations with valid ID and weight;
    IF (&ID>0 & &IID>0 & X42001>0);

*   dummy merge variable for adding quantile calculations after this
    data step;
    MERGEID=1;

*   divide weight by 5 so totals estimated on the 5 implicates
    jointly are correct;

    WGT=X42001/5;

*   retain original weight: WGT0;
    WGT0=X42001;

*   include hardcoded adjustments to weights intended to dampen the
    influence of overly influential cases identified by graphical
    analysis of the distributions of data;
    %IF (&PUBLIC NE YES) %THEN %DO;
      %INCLUDE HC&YEAR; 
    %END;


***************************************************************************;
*   miscellaneous global data adjustments;

*   consider all FINANCE COMPANIES reported in mortgage 
    grids to be MORTGAGE COMPANIES; 
    ARRAY MVAR{*} X9083 X9084 X9085 X9099 X9100 X9101;
    DO I=1 TO DIM(MVAR);
      IF MVAR{I}=14 THEN MVAR{I}=18;
    END; 

***************************************************************************;
*   demographic variables;

*   sex of household head;
    HHSEX=X8021;

*   age of the household head, and categorical variable:
    1:<35, 2:35-44, 3:45-54, 4:55-64, 5:65-74, 6:>=75;
    AGE=X14;
    AGECL=1+(AGE GE 35)+(AGE GE 45)+(AGE GE 55)+(AGE GE 65)+(AGE GE 75);


*   education of the HH head, and categorical variable:
    1=no high school diploma/GED, 2=high school diploma or GED,
    3=some college, 4=college degree;
    EDUC=X5901;
    IF X5904 EQ 1 THEN EDCL=4;
    ELSE IF EDUC GE 13 THEN EDCL=3;
    ELSE IF (X5902 IN (1 2)) THEN EDCL=2;
    ELSE EDCL=1;

*   marital status of the HH head: 1=married/living with partner,
    2=neither married nor living with partner;
    IF (X8023 IN (1 2)) THEN MARRIED=1;
    ELSE MARRIED=2;

*   number of children (including natural children/step-children/
    foster children of head/spouse/partner);
*   NOTE: from 1995 forward, household listing information collected
    for one fewer HH member;
*   NOTE: code 36, foster child not available in the public data set;
    ARRAY REL{*} X108 X114 X120 X126 X132 X202 X208 X214 X220 X226;
    KIDS=0;
    DO I=1 TO DIM(REL);
      KIDS=KIDS+(REL{I}=4|REL{I}=13|REL{I}=36);
    END;

*   labor force participation: 1=working in some way, 0=not working
    at all;
    IF ((X4100 >=50 & X4100 <= 80)| X4100=97) THEN LF=0;
      ELSE LF=1;

*   life cycle variables: 1=head under 55 + not married/LWP + no
    children, 2=head under 55 + married/LWP + no children, 3=head
    under 55 + married/LWP + children, 4=head under 55 + not
    married/LWP + children, 5=head 55 or older and working,
    6=head 55 or older and not working;
    IF AGE<55 AND MARRIED NE 1 AND KIDS EQ 0 THEN LIFECL=1;
    ELSE IF AGE<55 AND MARRIED EQ 1 AND KIDS EQ 0 THEN LIFECL=2;
    ELSE IF AGE<55 AND MARRIED EQ 1 AND KIDS GT 0 THEN LIFECL=3;
    ELSE IF AGE<55 AND MARRIED NE 1 AND KIDS GT 0 THEN LIFECL=4;
    ELSE IF AGE>54 AND LF=1 THEN LIFECL=5;
    ELSE IF AGE>54 AND LF=0 THEN LIFECL=6;    

*   family structure: 1=not married/LWP + children, 
    2=not married/LWP + no children + head under 55, 
    3=not married/LWP + no children + head 55 or older, 
    4=married/LWP+ children, 5=married/LWP + no children;

    IF MARRIED NE 1 AND KIDS >= 1 THEN FAMSTRUCT=1;
    ELSE IF MARRIED NE 1 AND KIDS = 0 AND AGE<55 THEN FAMSTRUCT=2;
    ELSE IF MARRIED NE 1 AND KIDS = 0 AND AGE>54 THEN FAMSTRUCT=3;
    ELSE IF MARRIED EQ 1 AND KIDS >= 1 THEN FAMSTRUCT=4;
    ELSE IF MARRIED EQ 1 AND KIDS =0 THEN FAMSTRUCT=5;
    ELSE PUT "ERROR: UNCLASSIFIED FAMILY STRUCTURE " Y1= MARRIED= KIDS= AGE=;

*   race/ethnicity;
*   NOTE: prior to 1998, the SCF only asked for one response.
    In 1998, respondents were allowed to give multiple responses,
    but they were asked to give first the category they identified
    with most strongly.  Few people gave more than one response.
    For purposes of continuity with prior data, define the 1998+
    variable in terms of the strongest identification;
*   beginning in 2004, respondents were also asked a question to
    determine whether they were of Hispanic/Latino culture or origin;
    %IF (&YEAR GE 2004) %THEN %DO;
*     1=white non-Hispanic, 2=nonwhite or Hispanic;
*     For the public data, we only keep the first two race variables
      and only code x6810 with a 1 if there is any response or a 5 if
      there is no response;
      %IF (&PUBLIC EQ NO) %THEN %DO;
        RACECL=1+(X6809 ^= 1 | X6810 NOT IN (1 0) | X6811 NOT IN (1 0) |
          X6812 NOT IN (1 0) | X6813 NOT IN (1 0) | X6814 NOT IN (1 0));
        H_RACECL=1+(X6809 ^= 1 | X6810 NOT IN (1 0) | X6811 NOT IN (1 0) |
          X6812 NOT IN (1 0) | X6813 NOT IN (1 0) | X6814 NOT IN (1 0) |
          X7004=1);
      %END;
      %ELSE %IF (&PUBLIC EQ YES) %THEN %DO;
        RACECL=1+(X6809 ^= 1 | X6810 ^=5);
        H_RACECL=1+(X6809 ^= 1 | X6810 ^=5 | X7004=1);
      %END;
*     1=white non-Hispanic, 2=black/African-American, 3=Hispanic,
      4=Asian (only available in internal data set, see codebook), 
      5=other;
      IF X6809=1 THEN RACE=1;
      ELSE IF X6809=2 THEN RACE=2;
      ELSE IF X6809=3 THEN RACE=3;
      ELSE IF X6809=4 THEN RACE=4;
      ELSE RACE=5;

      IF X6809=1 & X7004^=1 THEN H_RACE=1;
      ELSE IF X6809=2 & X7004^=1 THEN H_RACE=2;
      ELSE IF X6809=3 | X7004=1 THEN H_RACE=3;
      ELSE IF X6809=4 THEN H_RACE=4;
      ELSE H_RACE=5;
    %END;
    %ELSE %IF (&YEAR GE 1998) %THEN %DO;
*     1=white non-Hispanic, 2=nonwhite or Hispanic;
*     For the public data, we only keep the first two race variables
      and only code x6810 with a 1 if there is any response or a 5 if
      there is no response;
      %IF (&PUBLIC EQ NO) %THEN %DO;
        RACECL=1+(X6809 ^= 1 | X6810 NOT IN (1 0) | X6811 NOT IN (1 0) |
          X6812 NOT IN (1 0) | X6813 NOT IN (1 0) | X6814 NOT IN (1 0));
      %END;
      %ELSE %IF (&PUBLIC EQ YES) %THEN %DO;
        RACECL=1+(X6809 ^= 1 | X6810 ^=5);
      %END;
*     1=white non-Hispanic, 2=black/African-American, 3=Hispanic,
      4=Asian (only available in internal data set, see codebook), 
      5=other;
      IF X6809=1 THEN RACE=1;
      ELSE IF X6809=2 THEN RACE=2;
      ELSE IF X6809=3 THEN RACE=3;
      ELSE IF X6809=4 THEN RACE=4;
      ELSE RACE=5;
      H_RACECL=RACECL;
      H_RACE=RACE;
    %END;
    %ELSE %DO;
      RACECL=1+(X5909 ^= 5);
      IF X5909=5 THEN RACE=1;
      ELSE IF X5909=4 THEN RACE=2;
      ELSE IF X5909=3 THEN RACE=3;
      ELSE IF X5909=2 THEN RACE=4;
      ELSE RACE=5;
      H_RACECL=RACECL;
      H_RACE=RACE;
    %END;
    
*   work status categories for head:
    1=work for someone else, 2=self-employed/partnership,
    3=retired/disabled + (student/homemaker/misc. not working and
    age 65 or older), 4=other groups not working (mainly those under
    65 and out of the labor force);
    IF X4106=1 THEN OCCAT1=1;
    ELSE IF (X4106 IN (2 3 4)) THEN OCCAT1=2;
    ELSE IF ((X4100 IN (50 52))|(X4100 IN (21 23 30 70 80 97 85 -7) &
      AGE>=65)) THEN OCCAT1=3;
    ELSE IF (X14<65) THEN OCCAT1=4;
    ELSE PUT "ERROR: UNCLASSIFIED WORK STATUS " Y1= X14= X4100= X4106=;
    
*   occupation classification for head:
    1=managerial/professional 2=technical/sales/services,
    3=other (incl. production/craft/repair workers, operators, 
    laborers, farmers, foresters, fishers) 4=not working;
    %IF (&PUBLIC EQ NO) %THEN %DO;
      %IF (&YEAR LT 2004) %THEN %DO;
        IF (3 <= X7401 <=37) THEN OCCAT2=1;
        ELSE IF (43 <= X7401 <=199) THEN OCCAT2=1;
        ELSE IF (203 <= X7401 <=235) THEN OCCAT2=2;
        ELSE IF (243 <= X7401 <=285) THEN OCCAT2=2;
        ELSE IF (303 <= X7401 <=389) THEN OCCAT2=2;
        ELSE IF (403 <= X7401 <=407) THEN OCCAT2=2;
        ELSE IF (413 <= X7401 <=427) THEN OCCAT2=2;
        ELSE IF (433 <= X7401 <=469) THEN OCCAT2=2;
        ELSE IF (473 <= X7401 <=499) THEN OCCAT2=3;
        ELSE IF (503 <= X7401 <=699) THEN OCCAT2=3;
        ELSE IF (703 <= X7401 <=799) THEN OCCAT2=3;
        ELSE IF (803 <= X7401 <=859) THEN OCCAT2=3;
        ELSE IF (863 <= X7401 <=889) THEN OCCAT2=3;
        ELSE IF (903 <= X7401 <=905) THEN OCCAT2=3;
        ELSE IF (X7401=0) THEN OCCAT2=4;
        ELSE PUT "ERROR: UNCLASSIFIED OCCUPATION STATUS " Y1= X14= X4100= X4106= X7401=;
      %END;
      %ELSE %IF (&YEAR GE 2004 AND &YEAR LE 2010) %THEN %DO;
        IF (10 <= X7401 <=200) THEN OCCAT2=1;
        ELSE IF (220 <= X7401 <=1530) THEN OCCAT2=1;
        ELSE IF (1600 <= X7401 <=1860) THEN OCCAT2=1;
        ELSE IF (2000 <= X7401 <=3650) THEN OCCAT2=1;
        ELSE IF (1540 <= X7401 <=1560) THEN OCCAT2=2;
        ELSE IF (4700 <= X7401 <=5930) THEN OCCAT2=2;
        ELSE IF (1900 <= X7401 <=1960) THEN OCCAT2=2;
        ELSE IF (7900 <= X7401 <=7900) THEN OCCAT2=2;
        ELSE IF (3700 <= X7401 <=4320) THEN OCCAT2=2;
        ELSE IF (4400 <= X7401 <=4400) THEN OCCAT2=2;
        ELSE IF (4420 <= X7401 <=4650) THEN OCCAT2=2;
        ELSE IF (9830 <= X7401 <=9840) THEN OCCAT2=2;
        ELSE IF (6200 <= X7401 <=7850) THEN OCCAT2=3;
        ELSE IF (8330 <= X7401 <=8330) THEN OCCAT2=3;
        ELSE IF (8350 <= X7401 <=8350) THEN OCCAT2=3;
        ELSE IF (8440 <= X7401 <=8630) THEN OCCAT2=3;
        ELSE IF (8740 <= X7401 <=8760) THEN OCCAT2=3;
        ELSE IF (8810 <= X7401 <=8810) THEN OCCAT2=3;
        ELSE IF (4410 <= X7401 <=4410) THEN OCCAT2=3;
        ELSE IF (7920 <= X7401 <=8320) THEN OCCAT2=3;
        ELSE IF (8340 <= X7401 <=8340) THEN OCCAT2=3;
        ELSE IF (8360 <= X7401 <=8430) THEN OCCAT2=3;
        ELSE IF (8640 <= X7401 <=8730) THEN OCCAT2=3;
        ELSE IF (8800 <= X7401 <=8800) THEN OCCAT2=3;
        ELSE IF (8830 <= X7401 <=9750) THEN OCCAT2=3;
        ELSE IF (210 <= X7401 <=210) THEN OCCAT2=3;
        ELSE IF (4340 <= X7401 <=4350) THEN OCCAT2=3;
        ELSE IF (6000 <= X7401 <=6130) THEN OCCAT2=3;
        ELSE IF (X7401=0) THEN OCCAT2=4;
        ELSE PUT "ERROR: UNCLASSIFIED OCCUPATION STATUS " Y1= X14= X4100= X4106= X7401=;
      %END;
%* JXB - new in 2013 for new CPS occ codes;
      %ELSE %IF (&YEAR GE 2013) %THEN %DO;
        IF (10 <= X7401 <=200) THEN OCCAT2=1;
        ELSE IF (220 <= X7401 <=1530) THEN OCCAT2=1;
        ELSE IF (1600 <= X7401 <=1860) THEN OCCAT2=1;
        ELSE IF (2000 <= X7401 <=3655) THEN OCCAT2=1;
        ELSE IF (1540 <= X7401 <=1560) THEN OCCAT2=2;
        ELSE IF (4700 <= X7401 <=5940) THEN OCCAT2=2;
        ELSE IF (1900 <= X7401 <=1965) THEN OCCAT2=2;
        ELSE IF (7900 <= X7401 <=7900) THEN OCCAT2=2;
        ELSE IF (3700 <= X7401 <=4320) THEN OCCAT2=2;
        ELSE IF (4400 <= X7401 <=4400) THEN OCCAT2=2;
        ELSE IF (4420 <= X7401 <=4650) THEN OCCAT2=2;
        ELSE IF (9830 <= X7401 <=9840) THEN OCCAT2=2;
        ELSE IF (6200 <= X7401 <=7855) THEN OCCAT2=3;
        ELSE IF (8330 <= X7401 <=8350) THEN OCCAT2=3;
        ELSE IF (8360 <= X7401 <=8430) THEN OCCAT2=3;
        ELSE IF (8440 <= X7401 <=8730) THEN OCCAT2=3;
        ELSE IF (8740 <= X7401 <=8760) THEN OCCAT2=3;
        ELSE IF (8800 <= X7401 <=8800) THEN OCCAT2=3;
        ELSE IF (8810 <= X7401 <=8810) THEN OCCAT2=3;
        ELSE IF (8830 <= X7401 <=9750) THEN OCCAT2=3;
        ELSE IF (4410 <= X7401 <=4410) THEN OCCAT2=3;
        ELSE IF (7920 <= X7401 <=8320) THEN OCCAT2=3;
        ELSE IF (205 <= X7401 <=205) THEN OCCAT2=3;
        ELSE IF (4340 <= X7401 <=4350) THEN OCCAT2=3;
        ELSE IF (6000 <= X7401 <=6130) THEN OCCAT2=3;
        ELSE IF (X7401=0) THEN OCCAT2=4;
        ELSE PUT "ERROR: UNCLASSIFIED OCCUPATION STATUS " Y1= X14= X4100= X4106= X7401=;
      %END;
    %END;
    %IF (&PUBLIC EQ YES) %THEN %DO;
       IF X7401=1 THEN OCCAT2=1;
       ELSE IF X7401 IN(2 3) THEN OCCAT2=2;
       ELSE IF X7401 IN(4 5 6) THEN OCCAT2=3;
       ELSE IF X7401=0 THEN OCCAT2=4;
    %END;

*   industry classifications for head: 1=mining + construction +
    manufacturing, 2=transportation + communications + utilities and
    sanitary services + wholesale trade + finance, insurance and
    real estate, 3=agriculture + retail trade + services + public
    administration;
*   NOTE: for the public version of the dataset, the categories 2
    and 3 are combined;
    %IF (&PUBLIC EQ NO) %THEN %DO;
      %IF (&YEAR LT 2004) %THEN %DO;
        IF (OCCAT1>=3) THEN INDCAT=4;
        ELSE IF  (X7402>=40 & X7402<=392) THEN INDCAT=1;
        ELSE IF ((X7402>=400 & X7402<=571)|(X7402>=700 & X7402<=712))
          THEN INDCAT=2;
        ELSE IF ((X7402>0 & X7402<40)|(X7402>=580 & X7402=<691)|
          (X7402>=721)) THEN INDCAT=3;
      %END;
      %ELSE %DO;
        IF (OCCAT1>=3) THEN INDCAT=4;
        ELSE IF  ((370 <= X7402 <= 490) | (X7402=770) | 
          (1070 <= X7402 <= 3990)) THEN INDCAT=1; 
        ELSE IF ((570 <= X7402 <= 690) | (4070 <= X7402 <= 4590) | 
          (6070 <= X7402 <= 6390) | (6470 <= X7402 <= 6780) | 
          (6870 <= X7402 <= 7190)) THEN INDCAT=2; 
        ELSE IF ((170 <= X7402 <= 290) | (4670 <= X7402 <= 5790) |
          (7270 <= X7402 <= 9890)) THEN INDCAT=3; 
      %END;
    %END;
    %ELSE %DO;
      IF (OCCAT1>=3) THEN INDCAT=4;
      ELSE IF (X7402 IN (2 3)) THEN INDCAT=1;
      ELSE INDCAT=2;
    %END;
    
*   Census regions: 1=northeast, 2=north central, 3=south, 4=west;
    %IF (&YEAR NE 1989 AND &PUBLIC NE YES) %THEN %DO;
      REGION=X30022;
    %END;
    %ELSE %IF (&PUBLIC NE YES) %THEN %DO;
      REGION=X40083;
    %END;

*   Urbanicity: 1=MSA, 2=non-MSA;
    %IF (&PUBLIC NE YES) %THEN %DO;
      IF X8460<100000 THEN URBAN=1;
      ELSE URBAN=2;
    %END;

*   Annualized food spending - not asked prior to 2004;
    %IF &YEAR LE 2001 %THEN %DO;
       FOODHOME=0;
       FOODAWAY=0;
       FOODDELV=0;
    %END;
    %ELSE %IF &YEAR GE 2004 %THEN %DO;
       FOODHOME=((MAX(0,X3024))*%MCONV(X3025))*12;
       FOODAWAY=((MAX(0,X3029))*%MCONV(X3030))*12;
       FOODDELV=((MAX(0,X3027))*%MCONV(X3028))*12;
    %END;
 
***************************************************************************;
*   income variables;
    
*   HH income in previous calendar year;
*   NOTE: For 2004 forward, IRA and withdrawals from tax-
    deferred pension accounts added to INCOME below;
    INCOME=MAX(0,X5729);

*   HH income components in previous calendar year;
*   NOTE: Components of income may not sum to INCOME and in the public
*   data X5704/X5714/X5724 may have a value of -9 if X5729 was
*   negative and X5704/X5714/X5724 was also negative;
    WAGEINC=X5702;
    BUSSEFARMINC=X5704+X5714;
    INTDIVINC=X5706+X5708+X5710;
    KGINC=X5712;
    SSRETINC=X5722;
    TRANSFOTHINC=X5716+X5718+X5720+X5724;
    RENTINC=X5714;
*   for 2004 and beyond, add in the amount of withdrawals from IRAs
    and tax-deferred pension accounts (already included in earlier
    years). need to convert pension withdrawals to annual frequency,
    using macro that converts amounts to monthly, then multiply by 12;;
*   in 2010, the 5th and 6th current and future pension grids are droppped;
    %IF (&YEAR GE 2010) %THEN %DO;
      PENACCTWD=X6558+X6566+X6574+MAX(0,(X6464*%MCONV(X6465))*12)
        +MAX(0,(X6469*%MCONV(X6470))*12)+MAX(0,(X6474*%MCONV(X6475))*12)
        +MAX(0,(X6479*%MCONV(X6480))*12)+MAX(0,(X6965*%MCONV(X6966))*12)
        +MAX(0,(X6971*%MCONV(X6972))*12)+MAX(0,(X6977*%MCONV(X6978))*12)
        +MAX(0,(X6983*%MCONV(X6984))*12);
      INCOME=INCOME+PENACCTWD;
      SSRETINC=SSRETINC+PENACCTWD;
    %END;
    %ELSE %IF (&YEAR GE 2004) %THEN %DO;
      PENACCTWD=X6558+X6566+X6574+MAX(0,(X6464*%MCONV(X6465))*12)
        +MAX(0,(X6469*%MCONV(X6470))*12)+MAX(0,(X6474*%MCONV(X6475))*12)
        +MAX(0,(X6479*%MCONV(X6480))*12)+MAX(0,(X6484*%MCONV(X6485))*12)
        +MAX(0,(X6489*%MCONV(X6490))*12)+MAX(0,(X6965*%MCONV(X6966))*12)
        +MAX(0,(X6971*%MCONV(X6972))*12)+MAX(0,(X6977*%MCONV(X6978))*12)
        +MAX(0,(X6983*%MCONV(X6984))*12)+MAX(0,(X6989*%MCONV(X6990))*12)
        +MAX(0,(X6995*%MCONV(X6996))*12);
      INCOME=INCOME+PENACCTWD;
      SSRETINC=SSRETINC+PENACCTWD;
    %END;
    %ELSE %DO;
      PENACCTWD=0;
    %END;

*   normal income;
    NORMINC=.B;
    %IF (&YEAR GE 2004) %THEN %DO;
      IF (X7650^=3) THEN NORMINC=MAX(0,X7362)+PENACCTWD;
      ELSE NORMINC=INCOME;
    %END;
    %ELSE %IF (&YEAR GE 1995) %THEN %DO;
      IF (X7650^=3) THEN NORMINC=MAX(0,X7362);
      ELSE NORMINC=INCOME;
    %END;
    
*   if ADJINC=YES, adjust actual/normal income to level of survey year;
    %IF (&ADJINC EQ YES) %THEN %DO;
      ARRAY PYRINC {*} INCOME WAGEINC BUSSEFARMINC INTDIVINC KGINC 
        SSRETINC TRANSFOTHINC RENTINC;
      DO I=1 TO DIM(PYRINC);
        PYRINC{I}=PYRINC{I}*&CPILAG;
      END; 
      %IF (&YEAR GE 1995) %THEN %DO;
        NORMINC=NORMINC*&CPILAG;
      %END;
    %END;

***************************************************************************;
*   attitudinal and related variables;
    
*   adjusting for durables purchases/investments, spent
    more/same/less than income in past year;
*   WSAVED: 1=spending exceeded income, 2=spending equaled income,
    3=spending less than income;
*   SAVED: 1=spent less than income, 0=all others;
    %IF (&YEAR GE 1992) %THEN %DO;
      IF (X7508>0) THEN WSAVED=X7508;
      ELSE IF (X7510=2 & X7509=1) THEN WSAVED=3;
      ELSE WSAVED=X7510;
      SAVED=(WSAVED=3);
    %END;
    
*   reasons for saving: 1=cant save, 2=education, 3=family, 4=home,
    5=purchases, 6=retirement, 7=liquidity/the future, 8=investment,
    9=no particular reason;
*   NOTE: multiple saving reasons may be reported: here choosing only
    first (most important) reason;
    ARRAY SAVRES {*} SAVRES1-SAVRES9;
    ARRAY RES {*} X3006 /*X3007*/ ;
    DO I=1 TO DIM(SAVRES);
      SAVRES{I}=0;
    END;
    DO I=1 TO DIM(RES);
      IF (RES{I} IN (-2 -1)) THEN SAVRES{1}=1;
      ELSE IF (RES{I} IN (1 2)) THEN SAVRES{2}=1;
      ELSE IF (RES{I} IN (3 5 6)) THEN SAVRES{3}=1;
      ELSE IF (RES{I} EQ 11) THEN SAVRES{4}=1;
      ELSE IF (RES{I} IN (12 13 14 15 16 27 29 30 9 18 20 41)) THEN
        SAVRES{5}=1;
      ELSE IF (RES{I} IN (17 22)) THEN SAVRES{6}=1;
      ELSE IF (RES{I} IN (23 24 25 32 92 93)) THEN SAVRES{7}=1;
      ELSE IF (RES{I} IN (21 26 28)) THEN SAVRES{8}=1;
      ELSE IF (RES{I} IN (31 33 40 90 91 -7)) THEN SAVRES{9}=1;
      ELSE PUT "UNCLASSIFIED SAVING REASON! " &ID= RES{I}=;
    END;
/*
    1.  Children's education; education of grandchildren
    2.  Own education; spouse's education; education -- NA for whom
    3.  "For the children/family"  -- NFS; "to help the kids
        out"; estate
    5.  Wedding, Bar Mitzvah, and other ceremonies (except 17)
    6.  To have children/a family
    9.  To move (except 11)
   11.  Buying own house (code summer cottage in 12)
   12.  Purchase of cottage or second home for own use
   13.  Buy a car, boat or other vehicle
   14.  Home improvements/repairs
   15.  To travel; take vacations; take other time off
   16.  Buy durable household goods, appliances, home
        furnishings; hobby and recreational items; for other 
        purchases not codable above or not further specified; 
        "buy things when we need/want them"; moving/special occasions 
   17.  Burial/funeral expenses
   18.  Charitable or religious contributions
   20.  "To enjoy life"
   21.  Buying (investing in) own business/farm; equipment for
        business/farm
   22.  Retirement/old age
   23.  Reserves in case of unemployment
   24.  In case of illness; medical/dental expenses
   25.  Emergencies; "rainy days"; other unexpected needs; for
       "security" and independence
   26.  Investments reasons (to get interest, to be
        diversified, to buy other forms of assets)
   27.  To meet contractual commitments (debt repayment,
        insurance, taxes, etc.), to pay off house
   28.  "To get ahead;" to advance standard of living
   29.  Ordinary living expenses/bills
   31.  No reason (except 90, 91, 92)
   32.  "For the future"
   90.  Had extra income; saved becaused had the money left
        over -- no other purpose specified
   91.  Wise/prudent thing to do; good discipline to save; habit
   92.  Liquidity; to have cash available/on hand
   -1.  Don't/can't save; "have no money"
*/


*   R would spend more if assets appreciated in value:
    1=agree strongly, 2=agree somewhat, 3=neither agree nor
    disagree, 4=disagree somewhat, 5=disagree strongly;
    %IF (&YEAR GE 1998) %THEN %DO;
      SPENDMOR=X6789;
    %END;

*   R would spend less if assets depreciated in value:
    1=agree strongly, 2=agree somewhat, 3=neither agree nor
    disagree, 4=disagree somewhat, 5=disagree strongly;
    %IF (&YEAR GE 2010) %THEN %DO;
      SPENDLESS=X7492;
    %END;

*   Households overall expenses over last 12 months:
		1=unusually high, 2=unusually low, 3=normal;
    %IF (&YEAR GE 2010) %THEN %DO;
      EXPENSHILO=X7491;
    %END;
            
*   household had any late payments in last year;
    LATE=(X3004=5);
*   household had any payments more than 60 days past due in last year;
    LATE60=(X3005=1);
*   turndown for credit, feared denial, or either;
    TURNDOWN=(X407=1);
    FEARDENIAL=(X409=1);
    TURNFEAR=(TURNDOWN=1 | FEARDENIAL=1);
*   have a payday loan;
    %IF &YEAR GE 2007 %THEN %DO;
       HPAYDAY=(X7063=1);
    %END;
    %ELSE %DO;
       HPAYDAY=0;
    %END;
*  bankruptcy in the last five years;
   %IF &YEAR GE 1998 %THEN %DO;
      BNKRUPLAST5=(X6774 >= (&YEAR-5));
   %END;
   %ELSE %DO;
      BNKRUPLAST5=0;  
   %END;
   

***************************************************************************;
*  shopping for financial services;

    %IF (&YEAR GE 1995) %THEN %DO;   
      BSHOPNONE=(X7100=1);
      BSHOPGRDL=(X7100=5);
      BSHOPMODR=(X7100 IN (2 3 4));

      ISHOPNONE=(X7111=1);
      ISHOPGRDL=(X7111=5);
      ISHOPMODR=(X7111 IN (2 3 4));

    * information sources used in borrowing and investment;
      ARRAY SRCS {*} BCALL BMAGZNEWS BMAILADTV BINTERNET BFRIENDWORK BFINPRO BSELF 
         BDONT BOTHER ICALL IMAGZNEWS IMAILADTV IINTERNET IFRIENDWORK IFINPRO ISELF 
         IDONT IOTHER BFINPLAN IFINPLAN;
      DO J=1 TO DIM(SRCS);
        SRCS{J}=0;
      END;
      %IF (&YEAR GE 1998) %THEN %DO;
        %IF (&YEAR EQ 1998) %THEN %DO;
          ARRAY BORRSRC {*} X7101-X7110 X6849;	   
          ARRAY INVSRC {*} X7112-X7121;			   
        %END;
        %ELSE %DO;
          ARRAY BORRSRC {*} X7101-X7110 X6849 X6861-X6864;	   
          ARRAY INVSRC {*} X7112-X7121 X6865-X6869;			   
        %END;
        DO I=1 TO DIM(BORRSRC);						   
          IF BORRSRC{I}=1 THEN BCALL=1;
          ELSE IF BORRSRC{I}=2 THEN BMAGZNEWS=1;			   
          ELSE IF BORRSRC{I} IN (3 6 4 32) THEN BMAILADTV=1;
          ELSE IF BORRSRC{I}=5 THEN BINTERNET=1;	   
          ELSE IF BORRSRC{I} IN (7 18) THEN BFRIENDWORK=1;
          ELSE IF BORRSRC{I} IN (10 11 20 21 23 24) THEN BFINPRO=1;
          ELSE IF (BORRSRC{I} IN (8 9 12)) THEN BFINPLAN=1;
          ELSE IF BORRSRC{I} IN (13 17 19 22) THEN BSELF=1;
          ELSE IF BORRSRC{I} IN (14 16) THEN BDONT=1;	
          ELSE IF BORRSRC{I} IN (-7) THEN BOTHER=1;			   
          ELSE IF BORRSRC{I}^=0 THEN PUT "WARNING: UNCLASSIFIED BORROWING INFO SOURCE: " BORRSRC{I}=;
        END;
        DO I=1 TO DIM(INVSRC);
          IF INVSRC{I}=1 THEN ICALL=1;   
          ELSE IF INVSRC{I}=2 THEN IMAGZNEWS=1;			   
          ELSE IF INVSRC{I} IN (3 6 4 32) THEN IMAILADTV=1;
          ELSE IF INVSRC{I}=5 THEN IINTERNET=1;
          ELSE IF INVSRC{I} IN (7 18 19) THEN IFRIENDWORK=1;
          ELSE IF INVSRC{I} IN (10 11 20 23 24 25) THEN IFINPRO=1;
          ELSE IF (INVSRC{I} IN (8 9 12)) THEN IFINPLAN=1;
          ELSE IF INVSRC{I} IN (13 17 21 22) THEN ISELF=1;
          ELSE IF INVSRC{I} IN (14 16) THEN IDONT=1;	
          ELSE IF INVSRC{I} IN (-7) THEN IOTHER=1;			   
          ELSE IF INVSRC{I}^=0 THEN PUT "WARNING: UNCLASSIFIED INVESTMENT INFO SOURCE: " INVSRC{I}=;    
        END;	
      %END;

      %IF (&YEAR GE 2004) %THEN %DO;
        ARRAY INSTUSE {*} X6600 X6601 X6602 X6603 X6604 X6605 X6606 X6607			 
          X6870 X6871 X6872 X6873													 
          X6608 X6609 X6610 X6611 X6612 X6613 X6614 X6615							 
          X6874 X6875 X6876 X6877													 
          X6616 X6617 X6618 X6619 X6620 X6621 X6622 X6623							 
          X6878 X6879 X6880 X6881													 
          X6624 X6625 X6626 X6627 X6628 X6629 X6630 X6631							 
          X6882 X6883 X6884 X6885													 
          X6632 X6633 X6634 X6635 X6636 X6637 X6638 X6639							 
          X6886 X6887 X6888 X6889													 
          X6640 X6641 X6642 X6643 X6644 X6645 X6646 X6647							 
          X6890 X6891 X6892 X6893													 
          X6656 X6657 X6658 X6659 X6660 X6661 X6662 X6663							 
          X6894 X6895 X6896 X6897;                                                 
      %END;
      %ELSE %IF (&YEAR GE 1995) %THEN %DO;
        ARRAY INSTUSE {*} X6600 X6601 X6602 X6603 X6604 X6605 X6606 X6607			 
          X6870 X6871 X6872 X6873													 
          X6608 X6609 X6610 X6611 X6612 X6613 X6614 X6615							 
          X6874 X6875 X6876 X6877													 
          X6616 X6617 X6618 X6619 X6620 X6621 X6622 X6623							 
          X6878 X6879 X6880 X6881													 
          X6624 X6625 X6626 X6627 X6628 X6629 X6630 X6631							 
          X6882 X6883 X6884 X6885													 
          X6632 X6633 X6634 X6635 X6636 X6637 X6638 X6639							 
          X6886 X6887 X6888 X6889													 
          X6640 X6641 X6642 X6643 X6644 X6645 X6646 X6647							 
          X6890 X6891 X6892 X6893;													 
      %END;
      INTERNET=0;
      DO I=1 TO DIM(INSTUSE);
        IF (INSTUSE{I}=12) THEN INTERNET=1;
      END;
      IF (BINTERNET=1 | IINTERNET=1) THEN INTERNET=1;
    %END;


***************************************************************************;
*   assets, debts, networth, and related varaibles;

***************************************************************************;
*   financial assets and related variables;
    
*   checking accounts other than money market;
    CHECKING=MAX(0,X3506)*(X3507=5)+MAX(0,X3510)*(X3511=5)
      +MAX(0,X3514)*(X3515=5)+MAX(0,X3518)*(X3519=5)
      +MAX(0,X3522)*(X3523=5)+MAX(0,X3526)*(X3527=5)
      +MAX(0,X3529)*(X3527=5);
*   have any checking account: 1=yes, 0=no;
    HCHECK=(((X3507=5)+(X3511=5)+(X3515=5)+(X3519=5)
      +(X3523=5)+(X3527=5)+(X3527=5))>0);
*   have no checking account: 1=no checking, 0=have checking;
*   NOTE: NOCHK=0 may include instances where R has a money market
    account that is used for checking;
    NOCHK=(X3501=5);
*   people w/o checking accounts: ever had an account?: 1=yes, 5=no;
    EHCHKG=X3502;
*   people w/o checking accounts: why have no account?:
    1=don*t write enough checks to make it worthwhile,
    2=minimum balance is too high, 3=do not like dealing with banks,
    4=service charges are too high, 5=no bank has convenient hours
    or location, 12=checkbook has been/could be lost/stolen,
    13=haven*t gotten around to it, 14=R has alternative source of
    checking services (MMA, MIA, etc) (does not include individuals
    who write checks for R), 15=R not allowed to have account (e.g.,
    asset test for welfare), 20=R does not need/want a checking
    account (NEC), 21=credit problems, bankruptcy, R does not meet
    depository*s qualifications for having an account, 95=don*t have
    (enough) money, -1=can*t manage/balance a checking account,
    -7=other, 0=inapplicable. (R has a checking account: X3501=1);
*   NOTE: codeframe varies over the survey years, so beware of
    constructing overly specific comparisons of the distribution of
    households over these categories over time;
    WHYNOCKG=X3503;

    %IF (&YEAR<=1992) %THEN %DO;
      DONTWRIT=(WHYNOCKG=1);
      MINBAL=(WHYNOCKG=2);
      DONTLIKE=(WHYNOCKG=3);
      SVCCHG=(WHYNOCKG=4);
      NOMONEY=(WHYNOCKG=95); 
      CANTMANG=(WHYNOCKG=96); 
      CREDIT=.; 
      DONTWANT=.; 
      OTHER=(WHYNOCKG NOT IN(0 1 2 3 4 95 96 21 20));
    %END;
    %ELSE %DO;
      DONTWRIT=(WHYNOCKG=1);
      MINBAL=(WHYNOCKG=2);
      DONTLIKE=(WHYNOCKG=3);
      SVCCHG=(WHYNOCKG=4);
      CANTMANG=(WHYNOCKG=-1); 
      NOMONEY=(WHYNOCKG=95); 
      CREDIT=(WHYNOCKG=21); 
      DONTWANT=(WHYNOCKG=20); 
      OTHER=(WHYNOCKG NOT IN (0 1 2 3 4 -1 95 21 20));
    %END;
    
*   reason choose main checking institution - codeframe is similar
    across years, just more codes in later years;
    * code 50 added in 2010;
    CKLOCATION=(X3530=3);	
  	CKLOWFEEBAL=(X3530=7);
  	CKMANYSVCS=(X3530=6);	
   	CKRECOMFRND=(X3530=1);
   	CKPERSONAL=(X3530=11);
   	CKCONNECTN=(X3530=35);
   	CKLONGTIME=(X3530=14);
   	CKSAFETY=(X3530=8);	
   	CKCONVPAYRL=(X3530=9);
   	CKOTHCHOOSE=(X3530 NOT IN(0 3 7 6 1 11 35 14 8 9));

*   savings accounts;
    %IF (&YEAR LE 2001) %THEN %DO;
      SAVING=MAX(0,X3804)+MAX(0,X3807)+MAX(0,X3810)+MAX(0,X3813)
        +MAX(0,X3816)+MAX(0,X3818);
    %END;
    %ELSE %DO;
      SAVING=MAX(0,X3730*(X3732 NOT IN (4 30)))
        +MAX(0,X3736*(X3738 NOT IN (4 30)))
        +MAX(0,X3742*(X3744 NOT IN (4 30)))+MAX(0,X3748*(X3750 NOT IN (4 30)))
        +MAX(0,X3754*(X3756 NOT IN (4 30)))+MAX(0,X3760*(X3762 NOT IN (4 30)))
        +MAX(0,X3765);
    %END;
*   have savings account: 1=yes, 0=no;
    HSAVING=(SAVING>0);
  
*   money market deposit accounts;
*   NOTE: includes money market accounts used for checking and other
    money market account held at commercial banks, savings and
    loans, savings banks, and credit unions;
    %IF (&YEAR LE 2001) %THEN %DO;
      MMDA=MAX(0,X3506)*((X3507=1)*(11<=X9113<=13))
        +MAX(0,X3510)*((X3511=1)*(11<=X9114<=13))
        +MAX(0,X3514)*((X3515=1)*(11<=X9115<=13))
        +MAX(0,X3518)*((X3519=1)*(11<=X9116<=13))
        +MAX(0,X3522)*((X3523=1)*(11<=X9117<=13))
        +MAX(0,X3526)*((X3527=1)*(11<=X9118<=13))
        +MAX(0,X3529)*((X3527=1)*(11<=X9118<=13))
        +MAX(0,X3706)*(11<=X9131<=13)+MAX(0,X3711)*(11<=X9132<=13)
        +MAX(0,X3716)*(11<=X9133<=13)+MAX(0,X3718)*(11<=X9133<=13);
    %END;
    %ELSE %DO;
      MMDA=MAX(0,X3506)*((X3507=1)*(11<=X9113<=13))
        +MAX(0,X3510)*((X3511=1)*(11<=X9114<=13))
        +MAX(0,X3514)*((X3515=1)*(11<=X9115<=13))
        +MAX(0,X3518)*((X3519=1)*(11<=X9116<=13))
        +MAX(0,X3522)*((X3523=1)*(11<=X9117<=13))
        +MAX(0,X3526)*((X3527=1)*(11<=X9118<=13))
        +MAX(0,X3529)*((X3527=1)*(11<=X9118<=13))
        +MAX(0,X3730*(X3732 IN (4 30))*(X9259>=11 & X9259<=13))
        +MAX(0,X3736*(X3738 IN (4 30))*(X9260>=11 & X9260<=13))
        +MAX(0,X3742*(X3744 IN (4 30))*(X9261>=11 & X9261<=13))
        +MAX(0,X3748*(X3750 IN (4 30))*(X9262>=11 & X9262<=13))
        +MAX(0,X3754*(X3756 IN (4 30))*(X9263>=11 & X9263<=13))
        +MAX(0,X3760*(X3762 IN (4 30))*(X9264>=11 & X9264<=13))
        +MAX(0,X3765*(X3762 IN (4 30))*(X9264>=11 & X9264<=13));
    %END;
*   money market mutual funds;
*   NOTE: includes money market accounts used for checking and other
    money market account held at institutions other than commercial
    banks, savings and loans, savings banks, and credit unions;
    %IF (&YEAR LE 2001) %THEN %DO;
      MMMF=MAX(0,X3506)*(X3507=1)*(X9113<11|X9113>13)
        +MAX(0,X3510)*(X3511=1)*(X9114<11|X9114>13)
        +MAX(0,X3514)*(X3515=1)*(X9115<11|X9115>13)
        +MAX(0,X3518)*(X3519=1)*(X9116<11|X9116>13)
        +MAX(0,X3522)*(X3523=1)*(X9117<11|X9117>13)
        +MAX(0,X3526)*(X3527=1)*(X9118<11|X9118>13)
        +MAX(0,X3529)*(X3527=1)*(X9118<11|X9118>13)
        +MAX(0,X3706)*(X9131<11|X9131>13)+MAX(0,X3711)*(X9132<11|X9132>13)
        +MAX(0,X3716)*(X9133<11|X9133>13)+MAX(0,X3718)*(X9133<11|X9133>13);
    %END;
    %ELSE %DO;
      MMMF=MAX(0,X3506)*(X3507=1)*(X9113<11|X9113>13)
        +MAX(0,X3510)*(X3511=1)*(X9114<11|X9114>13)
        +MAX(0,X3514)*(X3515=1)*(X9115<11|X9115>13)
        +MAX(0,X3518)*(X3519=1)*(X9116<11|X9116>13)
        +MAX(0,X3522)*(X3523=1)*(X9117<11|X9117>13)
        +MAX(0,X3526)*(X3527=1)*(X9118<11|X9118>13)
        +MAX(0,X3529)*(X3527=1)*(X9118<11|X9118>13)
        +MAX(0,X3730*(X3732 IN (4 30))*(X9259<11|X9259>13))
        +MAX(0,X3736*(X3738 IN (4 30))*(X9260<11|X9260>13))
        +MAX(0,X3742*(X3744 IN (4 30))*(X9261<11|X9261>13))
        +MAX(0,X3748*(X3750 IN (4 30))*(X9262<11|X9262>13))
        +MAX(0,X3754*(X3756 IN (4 30))*(X9263<11|X9263>13))
        +MAX(0,X3760*(X3762 IN (4 30))*(X9264<11|X9264>13))
        +MAX(0,X3765*(X3762 IN (4 30))*(X9264<11|X9264>13));
    %END;
*   all types of money market accounts;
    MMA=MMDA+MMMF;
*   have any type of money market account: 1=yes, 0=no;
    HMMA=(MMA>0);
    
*   call accounts at brokerages;
    CALL=MAX(0,X3930);
*   have call account: 1=yes, 0=no;
    HCALL=(CALL>0);
    
*   all types of transactions accounts (liquid assets);
    LIQ=CHECKING+SAVING+MMA+CALL;
*   have any types of transactions accounts: 1=yes, 0=no;
*   here include even accounts with zero reported balances;
    %IF (&YEAR GE 2004) %THEN %DO;
      HLIQ=(LIQ>0 | X3501=1 | X3727=1 | X3929=1);
    %END;
    %ELSE %DO;
      HLIQ=(LIQ>0 | X3501=1 | X3701=1 | X3801=1 | X3929=1);
    %END;  
*   include accounts with zero balances (for tabling program);  
    LIQ=MAX(HLIQ,LIQ);

*   certificates of deposit;
    CDS=MAX(0,X3721);
*   have CDs: 1=yes, 0=no;
    HCDS=(CDS>0);
    
*   mutual funds;
*   stock mutual funds;
    STMUTF=(X3821=1)*MAX(0,X3822);
*   tax-free bond mutual funds;
    TFBMUTF=(X3823=1)*MAX(0,X3824);
*   government bond mutual funds;
    GBMUTF=(X3825=1)*MAX(0,X3826);
*   other bond mutual funds;
    OBMUTF=(X3827=1)*MAX(0,X3828);
*   combination and other mutual funds;
    COMUTF=(X3829=1)*MAX(0,X3830);
    %IF (&YEAR GE 2004) %THEN %DO;
*     other mutual funds;
      OMUTF=(X7785=1)*MAX(0,X7787);
    %END;
    %ELSE %DO;
      OMUTF=0;
    %END;
*   total directly-held mutual funds, excluding MMMFs;
    NMMF=STMUTF+TFBMUTF+GBMUTF+OBMUTF+COMUTF;
    %IF (&YEAR GE 2004) %THEN %DO;
      NMMF=NMMF+OMUTF;
    %END;
*   have any mutual funds excluding MMMFs: 1=yes, 0=no;
    HNMMF=(NMMF>0);
    
*   stocks;
    STOCKS=MAX(0,X3915);
*   have stocks: 1=yes, 0=no;
    HSTOCKS=(STOCKS>0);
*   number different companies in which hold stock;
    NSTOCKS=X3914;
*   Wilshire index of stock prices;
    %IF (&YEAR GE 1998) %THEN %DO;
      WILSH=X33001;
    %END;
    %ELSE %DO;
      WILSH=.B;
    %END;

*   bonds, not including bond funds or savings bonds;
*   tax-exempt bonds (state and local bonds);
    NOTXBND=X3910;
*   mortgage-backed bonds;
    MORTBND=X3906;
*   US government and government agency bonds and bills;
    GOVTBND=X3908;
*   corporate and foreign bonds;
    %IF (&YEAR GE 1992) %THEN %DO;
      OBND=X7634+X7633;
    %END;
    %ELSE %DO;
      OBND=X3912;
    %END;
*   total bonds, not including bond funds or savings bonds;
    BOND=NOTXBND+MORTBND+GOVTBND+OBND;
*   have bonds: 1=yes, 0=no;
    HBOND=(BOND>0);
    
*   quasi-liquid retirement accounts (IRAs and thrift-type accounts);
*   individual retirement accounts/Keoghs;
    %IF (&YEAR GE 2004) %THEN %DO;
      IRAKH=SUM(0,X6551,X6559,X6567,X6552,X6560,X6568,X6553,X6561,
        X6569,X6554,X6562,X6570);
    %END;
    %ELSE %DO;
      IRAKH=MAX(0,X3610)+MAX(0,X3620)+MAX(0,X3630);
    %END;
    %IF (&YEAR LT 2004) %THEN %DO;
*     account-type pension plans (included if type is 401k, 403b,
      thrift, savings, SRA, or if participant has option to borrow or
      withdraw);
*     PENEQ counts thrift amounts invested in stock;
      ARRAY PTYPE{*} X4216 X4316 X4416 X4816 X4916 X5016;
      ARRAY PAMT{*} X4226 X4326 X4426 X4826 X4926 X5026;
      ARRAY PBOR{*} X4227 X4327 X4427 X4827 X4927 X5027;
      ARRAY PWIT{*} X4231 X4331 X4431 X4831 X4931 X5031;
      ARRAY PALL{*} X4234 X4334 X4434 X4834 X4934 X5034;
      THRIFT = 0; PENEQ =0 ;
      RTHRIFT=0; STHRIFT=0;
      REQ=0; SEQ=0;
      DO I=1 TO DIM(PTYPE);
        HOLD=MAX(0,PAMT{I})*(PTYPE{I} IN (1 2 7 11 12 18) | 
          PBOR{I}=1|PWIT{I}=1);
        IF (I<=3) THEN RTHRIFT=RTHRIFT+HOLD;
        ELSE STHRIFT=STHRIFT+HOLD;
        THRIFT=THRIFT+HOLD;
        PENEQ =PENEQ+HOLD*((PALL{I}=1)+.5*(PALL{I}=3));
        IF (PTYPE{I}<=.Z) THEN PUT &ID= PTYPE{I}=;
        IF (PAMT{I}<=.Z) THEN PUT &ID= PAMT{I}=;
        IF (PBOR{I}<=.Z) THEN PUT &ID= PBOR{I}=;
        IF (PWIT{I}<=.Z) THEN PUT &ID= PWIT{I}=;
        IF (PALL{I}<=.Z) THEN PUT &ID= PALL{I}=;
        IF (I<=3) THEN REQ=PENEQ;
        ELSE SEQ=PENEQ-REQ;
      END;
*     allocate the pension mopups;
*     where possible, use information for first three pensions to infer
      characteristics of this amount;
*     where not possible to infer whether R can borrow/make withdrawals,
      assume this is possible;
*     where not possible to determine investment direction, assume half
      in stocks;
      PMOP=.;
      IF (X4436>0) THEN DO;
        IF (PTYPE{1} IN (1 2 7 11 12 18)|PTYPE{2} IN (1 2 7 11 12 18)|
           PTYPE{3} IN (1 2 7 11 12 18)| PWIT{1}=1|PWIT{2}=1|PWIT{3}=1|
           PBOR{1}=1|PBOR{2}=1|PBOR{3}=1) 
          THEN PMOP=X4436;
        ELSE IF (PTYPE{1}^=0 & PTYPE{2}^=0 & PTYPE{3}^=0 & PWIT{1}^=0 &
          PWIT{2}^=0 & PWIT{3}^=0) THEN PMOP=0;
        ELSE PMOP=X4436;
        THRIFT=THRIFT+PMOP;
        IF (REQ>0) THEN PENEQ=PENEQ+PMOP*(REQ/RTHRIFT);
        ELSE PENEQ=PENEQ+PMOP/2;
      END;
      PMOP=.;
      IF (X5036>0) THEN DO;
        IF (PTYPE{4} IN (1 2 7 11 12 18)|PTYPE{5} IN (1 2 7 11 12 18)|
           PTYPE{6} IN (1 2 7 11 12 18)| PWIT{4}=1|PWIT{5}=1|PWIT{6}=1| 
           PBOR{4}=1|PBOR{5}=1|PBOR{6}=1)
          THEN PMOP=X5036;
        ELSE IF (PTYPE{4}^=0 & PTYPE{5}^=0 & PTYPE{6}^=0 & PWIT{4}^=0 &
          PWIT{5}^=0 & PWIT{6}^=0) THEN PMOP=0;
        ELSE PMOP=X5036;
        THRIFT=THRIFT+PMOP;
        IF (SEQ>0) THEN PENEQ=PENEQ+PMOP*(SEQ/STHRIFT);
        ELSE PENEQ=PENEQ+PMOP/2;
      END;
      DROP HOLD PMOP;
    %END;
    %ELSE %IF (&YEAR GE 2004 & &YEAR LT 2010) %THEN %DO;
*     account-type pension plans (included if type is 401k, 403b,
      thrift, savings, SRA, or if participant has option to borrow or
      withdraw);
*     PENEQ counts thrift amounts invested in stock;
      ARRAY PTYPE1{*} X11000 X11100 X11200 X11300 X11400 X11500;
      ARRAY PTYPE2{*} X11001 X11101 X11201 X11301 X11401 X11501;
      ARRAY PAMT{*} X11032 X11132 X11232 X11332 X11432 X11532;
      ARRAY PBOR{*} X11025 X11125 X11225 X11325 X11425 X11525;
      ARRAY PWIT{*} X11031 X11131 X11231 X11331 X11431 X11531;
      ARRAY PALL{*} X11036 X11136 X11236 X11336 X11436 X11536;
      ARRAY PPCT{*} X11037 X11137 X11237 X11337 X11437 X11537;
      THRIFT = 0; PENEQ =0 ;
      RTHRIFT=0; STHRIFT=0;
      REQ=0; SEQ=0;
      DO I=1 TO DIM(PTYPE1);
        HOLD=MAX(0,PAMT{I})*(PTYPE1{I} IN (5 6 10 21)|
          PTYPE2{I} IN (2 3 4 6 20 21 22 26)|PBOR{I}=1|PWIT{I}=1);
        IF (I<=3) THEN RTHRIFT=RTHRIFT+HOLD;
        ELSE STHRIFT=STHRIFT+HOLD;
        THRIFT=THRIFT+HOLD;
        PENEQ=PENEQ+HOLD*((PALL{I}=1)+(PALL{I}=3)*(MAX(0,PPCT{I})/10000));
        IF (I<=3) THEN REQ=PENEQ;
        ELSE SEQ=PENEQ-REQ;
      END;
*     allocate the pension mopups;
*     where possible, use information for first three pensions to infer
      characteristics of this amount;
*     where not possible to infer whether R can borrow/make withdrawals,
      assume this is possible;
*     where not possible to determine investment direction, assume half
      in stocks;
      HOLD=.; PMOP=.;
      IF (X11259>0) THEN DO;
        IF (PTYPE1{1} IN (5 6 10 21)|PTYPE1{2} IN (5 6 10 21)|
          PTYPE1{3} IN (5 6 10 21)| PTYPE2{1} IN (2 3 4 6 20 21 22 26) |
          PTYPE2{2} IN (2 3 4 6 20 21 22 26) | 
          PTYPE2{3} IN (2 3 4 6 20 21 22 26) | 
          PWIT{1}=1|PWIT{2}=1|PWIT{3}=1|PBOR{1}=1|PBOR{2}=1|PBOR{3}=1)
          THEN PMOP=X11259;
        ELSE IF (PTYPE1{1}^=0 & PTYPE1{2}^=0 & PTYPE1{3}^=0 & PWIT{1}^=0 &
          PWIT{2}^=0 & PWIT{3}^=0) THEN PMOP=0;
        ELSE PMOP=X11259;
        THRIFT=THRIFT+PMOP;
        IF (REQ>0) THEN PENEQ=PENEQ+PMOP*(REQ/RTHRIFT);
        ELSE PENEQ=PENEQ+PMOP/2;
      END;
      IF (X11559>0) THEN DO;
        IF (PTYPE1{4} IN (5 6 10 21)|PTYPE1{5} IN (5 6 10 21) | 
          PTYPE1{6} IN (5 6 10 21)| PTYPE2{4} IN (2 3 4 6 20 21 22 26) | 
          PTYPE2{5} IN (2 3 4 6 20 21 22 26) | 
          PTYPE2{6} IN (2 3 4 6 20 21 22 26) | 
          PWIT{4}=1|PWIT{5}=1|PWIT{6}=1|PBOR{4}=1|PBOR{5}=1|PBOR{6}=1)
          THEN PMOP=X11559;
        ELSE IF (PTYPE1{4}^=0 & PTYPE1{5}^=0 & PTYPE1{6}^=0 & PWIT{4}^=0 &
          PWIT{5}^=0 & PWIT{6}^=0) THEN PMOP=0;
        ELSE PMOP=X11559;
        THRIFT=THRIFT+PMOP;
        IF (SEQ>0) THEN PENEQ=PENEQ+PMOP*(SEQ/STHRIFT);
        ELSE PENEQ=PENEQ+PMOP/2;
      END;
      DROP HOLD PMOP RTHRIFT STHRIFT REQ SEQ;
    %END;
    %ELSE %IF (&YEAR GE 2010) %THEN %DO;
*     account-type pension plans (included if type is 401k, 403b,
      thrift, savings, SRA, or if participant has option to borrow or
      withdraw);
*     PENEQ counts thrift amounts invested in stock;
      ARRAY PTYPE1{*} X11000 X11100 X11300 X11400;
      ARRAY PTYPE2{*} X11001 X11101 X11301 X11401;
      ARRAY PAMT{*} X11032 X11132 X11332 X11432;
      ARRAY PBOR{*} X11025 X11125 X11325 X11425;
      ARRAY PWIT{*} X11031 X11131 X11331 X11431;
      ARRAY PALL{*} X11036 X11136 X11336 X11436;
      ARRAY PPCT{*} X11037 X11137 X11337 X11437;
      THRIFT = 0; PENEQ =0 ;
      RTHRIFT=0; STHRIFT=0;
      REQ=0; SEQ=0;
      DO I=1 TO DIM(PTYPE1);
        HOLD=MAX(0,PAMT{I})*(PTYPE1{I} IN (1)|
          PTYPE2{I} IN (2 3 4 6 20 21 22 26)|PBOR{I}=1|PWIT{I}=1);
        IF (I<=2) THEN RTHRIFT=RTHRIFT+HOLD;
        ELSE STHRIFT=STHRIFT+HOLD;
        THRIFT=THRIFT+HOLD;
        PENEQ=PENEQ+HOLD*((PALL{I}=1)+(PALL{I} IN(3 30))*(MAX(0,PPCT{I})/10000));
        IF (I<=2) THEN REQ=PENEQ;
        ELSE SEQ=PENEQ-REQ;
      END;
*     allocate the pension mopups;
*     where possible, use information for first three pensions to infer
      characteristics of this amount;
*     where not possible to infer whether R can borrow/make withdrawals,
      assume this is possible;
*     where not possible to determine investment direction, assume half
      in stocks;
      HOLD=.; PMOP=.;
      IF (X11259>0) THEN DO;
        IF (PTYPE1{1} IN (1)|PTYPE1{2} IN (1)|
          PTYPE2{1} IN (2 3 4 6 20 21 22 26) |
          PTYPE2{2} IN (2 3 4 6 20 21 22 26) |  
          PWIT{1}=1|PWIT{2}=1|PBOR{1}=1|PBOR{2}=1)
          THEN PMOP=X11259;
        ELSE IF (PTYPE1{1}^=0 & PTYPE1{2}^=0 & PWIT{1}^=0 &
          PWIT{2}^=0 ) THEN PMOP=0;
        ELSE PMOP=X11259;
        THRIFT=THRIFT+PMOP;
        IF (REQ>0) THEN PENEQ=PENEQ+PMOP*(REQ/RTHRIFT);
        ELSE PENEQ=PENEQ+PMOP/2;
      END;
      IF (X11559>0) THEN DO;
        IF (PTYPE1{3} IN (1)|PTYPE1{4} IN (1) | 
          PTYPE2{3} IN (2 3 4 6 20 21 22 26) | 
          PTYPE2{4} IN (2 3 4 6 20 21 22 26) |  
          PWIT{3}=1|PWIT{4}=1|PBOR{3}=1|PBOR{4}=1)
          THEN PMOP=X11559;
        ELSE IF (PTYPE1{3}^=0 & PTYPE1{4}^=0 & PWIT{3}^=0 &
          PWIT{4}^=0 ) THEN PMOP=0;
        ELSE PMOP=X11559;
        THRIFT=THRIFT+PMOP;
        IF (SEQ>0) THEN PENEQ=PENEQ+PMOP*(SEQ/STHRIFT);
        ELSE PENEQ=PENEQ+PMOP/2;
      END;
      DROP HOLD PMOP RTHRIFT STHRIFT REQ SEQ;
    %END;
*   future pensions (accumulated in an account for the R/S); 
    %IF (&YEAR GE 2010) %THEN %DO;
      FUTPEN=MAX(0,X5604)+MAX(0,X5612)+MAX(0,X5620)+MAX(0,X5628); 
    %END;
    %ELSE %DO;
    FUTPEN=MAX(0,X5604)+MAX(0,X5612)+MAX(0,X5620)+MAX(0,X5628)+
      MAX(0,X5636)+MAX(0,X5644); 
    %END;
*   NOTE: there is very little evidence that pensions with currently
    received benefits recorded in the SCFs before 2001 were any type
    of 401k or related account from which the R was making
    withdrawals:  the questions added in 2001 allow one to distinguish
    such account, and there are 55 of them in that year:
    create a second version of RETQLIQ to include this information;
    %IF (&YEAR GE 2004 AND &YEAR LE 2007) %THEN %DO;
      CURRPEN=X6462+X6467+X6472+X6477+X6482+X6487+X6957;
*     total quasi-liquid: sum of IRAs, thrift accounts, and future pensions; 
*     this version includes currently received benefits;
      RETQLIQ=IRAKH+THRIFT+FUTPEN+CURRPEN;
*     have quasi-liquid assets: 1=yes, 0=no;
      HRETQLIQ=(RETQLIQ>0);
    %END;
    %ELSE %IF (&YEAR GE 2010) %THEN %DO;
      CURRPEN=X6462+X6467+X6472+X6477+X6957;
*     total quasi-liquid: sum of IRAs, thrift accounts, and future pensions; 
*     this version includes currently received benefits;
      RETQLIQ=IRAKH+THRIFT+FUTPEN+CURRPEN;
*     have quasi-liquid assets: 1=yes, 0=no;
      HRETQLIQ=(RETQLIQ>0);
    %END;
    %ELSE %IF (&YEAR EQ 2001) %THEN %DO;
      CURRPEN=X6462+X6467+X6472+X6477+X6482+X6487;
*     total quasi-liquid: sum of IRAs, thrift accounts, and future pensions; 
*     this version includes currently received benefits;
      RETQLIQ=IRAKH+THRIFT+FUTPEN+CURRPEN;
*     have quasi-liquid assets: 1=yes, 0=no;
      HRETQLIQ=(RETQLIQ>0);
    %END;
    %ELSE %DO;
      CURRPEN=0;
*     total quasi-liquid: sum of IRAs, thrift accounts, and future pensions; 
*     assuming that account type plans with currently received
      benefits are negligible;
      RETQLIQ=IRAKH+THRIFT+FUTPEN;
*     have quasi-liquid assets: 1=yes, 0=no;
      HRETQLIQ=(RETQLIQ>0);
    %END;
*   other pension characteristics:
    ANYPEN: 1=either head or spouse/partner has any type of
    pension, 0=otherwise,
    DBPLANCJ: 1=either head or spouse/partner has a defined benefit
    pension on a current job, 0=otherwise,
    DBPLANT: 1=either head or spouse/partner has DB plan on current
    job or some type of pension from a past job to be received in
    the future, 0=otherwise,
    DCPLANCJ: 1=either head or spouse/partner has any type of
    account-based plan on a current job, 0=otherwise
    BPLANCJ: 1=either head or spouse/partner has both types of pension
    plan on a current job
    NOTE: For DCPLANCJ and DBPLANCJ checking for plans from the
    current job from which R or S/P is currently receiving benefits -
    prior to 2001 do not know if these are account or income plans, so
    assume all current job currently receiving pension plans are defined benefit;
    ANYPEN=(X4135=1|X4735=1|X5313=1|X5601=1); 
    %IF (&YEAR GE 2010) %THEN %DO;
      DBPLANCJ=((X11000=5 & X11001^=30)|X11001=1 |
        (X11100=5 & X11101^=30)|X11101=1|
        (X11300=5 & X11301^=30)|X11301=1|(X11400=5 & X11401^=30)|X11401=1|
        (X5316=1 & X6461=5) 
        | (X5324=1 & X6466=5) | (X5332=1 & X6471=5) | (X5416=1 & X6476=5)) ;
      DCPLANCJ=(X11032>0|X11132>0|X11332>0|X11432>0|
        X11032=-1|X11132=-1|X11332=-1|X11432=-1|
        (X5316=1 & X6461=1) | (X5324=1 & X6466=1)| 
        (X5332=1 & X6471=1) | (X5416=1 & X6476=1));
    %END;
    %ELSE %IF (&YEAR GE 2004) %THEN %DO;
      DBPLANCJ=((X11000=4 & X11001^=30)|X11001=1 |
        (X11100=4 & X11101^=30)|X11101=1|(X11200=4 & X11201^=30)|X11201=1|
        (X11300=4 & X11301^=30)|X11301=1|(X11400=4 & X11401^=30)|X11401=1|
        (X11500=4 & X11501^=30)|X11501=1 | (X5316=1 & X6461=5) 
        | (X5324=1 & X6466=5) | (X5332=1 & X6471=5) | (X5416=1 & X6476=5)
        | (X5424=1 & X6481=5) | (X5432=1 & X6486=5)) ;
      DCPLANCJ=(X11032>0|X11132>0|X11232>0|X11332>0|X11432>0|
        X11532>0|X11032=-1|X11132=-1|X11232=-1|X11332=-1|X11432=-1|
        X11532=-1 | (X5316=1 & X6461=1) | (X5324=1 & X6466=1) 
        | (X5332=1 & X6471=1) | (X5416=1 & X6476=1) 
        | (X5424=1 & X6481=1) | (X5432=1 & X6486=1));
    %END;
    %ELSE %IF (&YEAR EQ 2001) %THEN %DO;
      DBPLANCJ=(X4203=1|X4303=1|X4403=1|X4803=1|X4903=1|X5003=1 
        | (X5316=1 & X6461=5) | (X5324=1 & X6466=5) 
        | (X5332=1 & X6471=5) | (X5416=1 & X6476=5)
        | (X5424=1 & X6481=5) | (X5432=1 & X6486=5));
      DCPLANCJ=((X4203 IN (2 3))|(X4303 IN (2 3))|(X4403 IN (2 3))
        |(X4803 IN (2 3))|(X4903 IN (2 3))|(X5003 IN (2 3))
        | (X5316=1 & X6461=1) | (X5324=1 & X6466=1) 
        | (X5332=1 & X6471=1) | (X5416=1 & X6476=1) 
        | (X5424=1 & X6481=1) | (X5432=1 & X6486=1));
    %END;
    %ELSE %DO;
      DBPLANCJ=(X4203=1|X4303=1|X4403=1|X4803=1|X4903=1|X5003=1
        |X5316=1|X5324=1|X5332=1|X5416=1|X5424=1|X5432=1);
      DCPLANCJ=((X4203 IN (2 3))|(X4303 IN (2 3))|(X4403 IN (2 3))
        |(X4803 IN (2 3))|(X4903 IN (2 3))|(X5003 IN (2 3)));
    %END;
*     NOTE! some plans captured by the X5314>0 conditions may
      not be DB plans;
    %IF (&YEAR GE 2010) %THEN %DO;
      DBPLANT=(DBPLANCJ=1|(X6461=5|X6466=5|X6471=5|X6476=5)
        |(X5603 IN (1 3)|X5611 IN (1 3)|X5619 IN (1 3)|
        X5627 IN (1 3)));
    %END;
    %ELSE %IF (&YEAR GE 2001) %THEN %DO;
      DBPLANT=(DBPLANCJ=1|(X6461=5|X6466=5|X6471=5|X6476=5|X6481=5|
        X6486=5)|(X5603 IN (1 3)|X5611 IN (1 3)|X5619 IN (1 3)|
        X5627 IN (1 3)|X5635 IN (1 3)|X5643 IN (1 3)));
    %END;
    %ELSE %DO;
      DBPLANT =(DBPLANCJ=1|(X5314>0)|(X5603 IN (1 3)|X5611 IN (1 3)|
        X5619 IN (1 3)|X5627 IN (1 3)|X5635 IN (1 3)|X5643 IN (1 3)));
    %END;
    BPLANCJ =(DBPLANCJ=1 & DCPLANCJ=1); 
    
*   savings bonds;
    SAVBND=X3902;
*   have savings bonds: 1=yes, 0=no;
    HSAVBND=(SAVBND>0);
    
*   cash value of whole life insurance;
    CASHLI=MAX(0,X4006);
*   have cash value LI: 1=yes, 0=no;
    HCASHLI=(CASHLI>0);
    
*   other managed assets (trusts, annuities and managed investment
    accounts in which HH has equity interest);
    %IF (&YEAR GE 2004) %THEN %DO;
      ANNUIT=MAX(0,X6577);
      TRUSTS=MAX(0,X6587);
      OTHMA=ANNUIT+TRUSTS;
    %END;
    %ELSE %IF (&YEAR GE 1998 AND &YEAR LE 2001) %THEN %DO;
      ANNUIT=MAX(0,X6820);
      TRUSTS=MAX(0,X6835);
      OTHMA=ANNUIT+TRUSTS;
    %END;
    %ELSE %DO;
      OTHMA=MAX(0,X3942);
      ANNUIT=((X3935=1)/MAX(1,((X3934=1)+(X3935=1)+(X3936=1)+(X3937=1))))*MAX(0,X3942);
      TRUSTS=OTHMA-ANNUIT;
    %END;
*   have other managed assets: 1=yes, 0=no;
    HOTHMA=(OTHMA>0);
    
*   other financial assets: includes loans from the household to
    someone else, future proceeds, royalties, futures, non-public
    stock, deferred compensation, oil/gas/mineral invest., cash
    n.e.c.;
*   NOTE: because of the collapsing of categories in the public
    version of the dataset, both codes 71 (oil/gas/mineral leases or
    investments) and 72 (futures contracts, stock options) are
    combined in code 71: thus, the sum will be treated as a
    nonfinancial asset.
    Additionally, codes 77 (future lottery/prize receipts) and 79
    (other obligations to R, tax credits) are combined in code 77;
    %IF (&PUBLIC EQ NO) %THEN %DO;
      OTHFIN=X4018+X4022*(X4020 IN (61,62,63,64,65,66,67,72,73,74,77,79
        80,81,82,83,84))+
        X4026*(X4024 IN (61,62,63,64,65,66,67,72,73,74,77,79,80,81,82,83,84))+
        X4030*(X4028 IN (61,62,63,64,65,66,67,72,73,74,77,79,80,81,82,83,84));
    %END;
    %ELSE %DO;
      OTHFIN=X4018+X4022*(X4020 IN (61,62,63,64,65,66,71,72,73,74,77,
        80,81,-7))+
        X4026*(X4024 IN (61,62,63,64,65,66,71,72,73,74,77,80,81,-7))+
        X4030*(X4028 IN (61,62,63,64,65,66,71,72,73,74,77 80,81,-7));
    %END;
*   have other financial assets: 1=yes, 0=no;
    HOTHFIN=(OTHFIN>0);
    
*   financial assets invested in stock:
    1. directly-held stock 
    2. stock mutual funds: full value if described as stock mutual fund,
       1/2 value of combination mutual funds
    3. IRAs/Keoghs invested in stock: 
       full value if mostly invested in stock, 
       1/2 value if split between stocks/bonds or stocks/money market,
       1/3 value if split between stocks/bonds/money market
    4. other managed assets w/equity interest (annuities, trusts, MIAs):
       full value if mostly invested in stock,
       1/2 value if split between stocks/MFs & bonds/CDs, or
    	 "mixed/diversified,"
       1/3 value if "other"
    5. thrift-type retirement accounts invested in stock
       full value if mostly invested in stock
       1/2 value if split between stocks and interest earning assets
    6. savings accounts classified as 529 or other accounts that may
       be invested in stocks;
    %IF (&YEAR GE 2010) %THEN %DO;
      EQUITY=STOCKS+STMUTF+.5*COMUTF+OMUTF+
        (X6551+X6552+X6553+X6554)*((X6555=1)+(X6555 IN(3 30))*(MAX(0,X6556)/10000))+
        (X6559+X6560+X6561+X6562)*((X6563=1)+(X6563 IN(3 30))*(MAX(0,X6564)/10000))+
        (X6567+X6568+X6569+X6570)*((X6571=1)+(X6571 IN(3 30))*(MAX(0,X6572)/10000))+
        ANNUIT*((X6581=1)+(X6581 IN(3 30))*(MAX(0,X6582)/10000))+
        TRUSTS*((X6591=1)+(X6591 IN(3 30))*(MAX(0,X6592)/10000))+PENEQ+
        (X6461=1)*X6462*((X6933=1)+(X6933 IN(3 30))*(MAX(0,X6934)/10000))+
        (X6466=1)*X6467*((X6937=1)+(X6937 IN(3 30))*(MAX(0,X6938)/10000))+
        (X6471=1)*X6472*((X6941=1)+(X6941 IN(3 30))*(MAX(0,X6942)/10000))+
        (X6476=1)*X6477*((X6945=1)+(X6945 IN(3 30))*(MAX(0,X6946)/10000))+
        X5604*((X6962=1)+(X6962 IN(3 30))*(MAX(0,X6963)/10000))+
        X5612*((X6968=1)+(X6968 IN(3 30))*(MAX(0,X6969)/10000))+
        X5620*((X6974=1)+(X6974 IN(3 30))*(MAX(0,X6975)/10000))+
        X5628*((X6980=1)+(X6980 IN(3 30))*(MAX(0,X6981)/10000))+
        X3730*((X7074=1)+(X7074 IN(3 30))*(MAX(0,X7075)/10000))+
        X3736*((X7077=1)+(X7077 IN(3 30))*(MAX(0,X7078)/10000))+
        X3742*((X7080=1)+(X7080 IN(3 30))*(MAX(0,X7081)/10000))+
        X3748*((X7083=1)+(X7083 IN(3 30))*(MAX(0,X7084)/10000))+
        X3754*((X7086=1)+(X7086 IN(3 30))*(MAX(0,X7087)/10000))+
        X3760*((X7089=1)+(X7089 IN(3 30))*(MAX(0,X7090)/10000));

    %END;
    %ELSE %IF (&YEAR GE 2007) %THEN %DO;
      EQUITY=STOCKS+STMUTF+.5*COMUTF+OMUTF+
        (X6551+X6552+X6553+X6554)*((X6555=1)+(X6555=3)*X6556/10000)+
        (X6559+X6560+X6561+X6562)*((X6563=1)+(X6563=3)*X6564/10000)+
        (X6567+X6568+X6569+X6570)*((X6571=1)+(X6571=3)*X6572/10000)+
        ANNUIT*((X6581=1)+(X6581=3)*X6582/10000)+
        TRUSTS*((X6591=1)+(X6591=3)*X6592/10000)+PENEQ+
        (X6461=1)*X6462*((X6933=1)+(X6933=3)*X6934/10000)+
        (X6466=1)*X6467*((X6937=1)+(X6937=3)*X6938/10000)+
        (X6471=1)*X6472*((X6941=1)+(X6941=3)*X6942/10000)+
        (X6476=1)*X6477*((X6945=1)+(X6945=3)*X6946/10000)+
        (X6481=1)*X6482*((X6949=1)+(X6949=3)*X6950/10000)+
        (X6486=1)*X6487*((X6953=1)+(X6953=3)*X6954/10000)+
        X5604*((X6962=1)+(X6962=3)*X6963/10000)+
        X5612*((X6968=1)+(X6968=3)*X6969/10000)+
        X5620*((X6974=1)+(X6974=3)*X6975/10000)+
        X5628*((X6980=1)+(X6980=3)*X6981/10000)+
        X5636*((X6986=1)+(X6986=3)*X6987/10000)+
        X5644*((X6992=1)+(X6992=3)*X6993/10000)+
        X3730*((X7074=1)+(X7074 IN(3))*(MAX(0,X7075)/10000))+
        X3736*((X7077=1)+(X7077 IN(3))*(MAX(0,X7078)/10000))+
        X3742*((X7080=1)+(X7080 IN(3))*(MAX(0,X7081)/10000))+
        X3748*((X7083=1)+(X7083 IN(3))*(MAX(0,X7084)/10000))+
        X3754*((X7086=1)+(X7086 IN(3))*(MAX(0,X7087)/10000))+
        X3760*((X7089=1)+(X7089 IN(3))*(MAX(0,X7090)/10000));
    %END;
    %ELSE %IF (&YEAR EQ 2004) %THEN %DO;
      EQUITY=STOCKS+STMUTF+.5*COMUTF+OMUTF+
        (X6551+X6552+X6553+X6554)*((X6555=1)+(X6555=3)*X6556/10000)+
        (X6559+X6560+X6561+X6562)*((X6563=1)+(X6563=3)*X6564/10000)+
        (X6567+X6568+X6569+X6570)*((X6571=1)+(X6571=3)*X6572/10000)+
        ANNUIT*((X6581=1)+(X6581=3)*X6582/10000)+
        TRUSTS*((X6591=1)+(X6591=3)*X6592/10000)+PENEQ+
        (X6461=1)*X6462*((X6933=1)+(X6933=3)*X6934/10000)+
        (X6466=1)*X6467*((X6937=1)+(X6937=3)*X6938/10000)+
        (X6471=1)*X6472*((X6941=1)+(X6941=3)*X6942/10000)+
        (X6476=1)*X6477*((X6945=1)+(X6945=3)*X6946/10000)+
        (X6481=1)*X6482*((X6949=1)+(X6949=3)*X6950/10000)+
        (X6486=1)*X6487*((X6953=1)+(X6953=3)*X6954/10000)+
        X5604*((X6962=1)+(X6962=3)*X6963/10000)+
        X5612*((X6968=1)+(X6968=3)*X6969/10000)+
        X5620*((X6974=1)+(X6974=3)*X6975/10000)+
        X5628*((X6980=1)+(X6980=3)*X6981/10000)+
        X5636*((X6986=1)+(X6986=3)*X6987/10000)+
        X5644*((X6992=1)+(X6992=3)*X6993/10000);
    %END;
    %ELSE %IF (&YEAR EQ 2001) %THEN %DO;
      EQUITY=STOCKS+STMUTF+.5*COMUTF+
        IRAKH*((X3631=2) +.5*(X3631=5|X3631=6)+.3*(X3631=4))+
        ANNUIT*((X6826=1)+.5*(X6826=5|X6826=6)+.3*(X6826=-7))+
        TRUSTS*((X6841=1)+.5*(X6841=5|X6841=6)+.3*(X6841=-7))+
        PENEQ+X6462*((X6463=1)+.5*(X6463=3))+X6467*((X6468=1)+.5*(X6468=3))
        +X6472*((X6473=1)+.5*(X6473=3))+X6477*((X6478=1)+.5*(X6478=3))
        +X6482*((X6483=1)+.5*(X6483=3))+X6487*((X6488=1)+.5*(X6488=3))
        +X5604*((X6491=1)+.5*(X6491=3))+X5612*((X6492=1)+.5*(X6492=3))
        +X5620*((X6493=1)+.5*(X6493=3))+X5628*((X6494=1)+.5*(X6494=3))
        +X5636*((X6495=1)+.5*(X6495=3))+X5644*((X6496=1)+.5*(X6496=3)) ;
    %END;
    %ELSE %IF (&YEAR EQ 1998) %THEN %DO;
      EQUITY=STOCKS+STMUTF+.5*COMUTF+
        IRAKH*((X3631=2) +.5*(X3631=5|X3631=6)+.3*(X3631=4))+
        ANNUIT*((X6826=1)+.5*(X6826=5|X6826=6)+.3*(X6826=-7))+
        TRUSTS*((X6841=1)+.5*(X6841=5|X6841=6)+.3*(X6841=-7))+
        PENEQ;
    %END;
    %ELSE %DO;
      EQUITY=STOCKS+STMUTF+.5*COMUTF+
        IRAKH*((X3631=2)+.5*(X3631=5|X3631=6)+.3*(X3631=4))+
        OTHMA*((X3947=1)+.5*(X3947=5|X3947=6)+.3*(X3947=4|X3947=-7))+
        PENEQ;
    %END;
*   have stock equity: 1=yes, 0=no;
    HEQUITY=(EQUITY>0);
*   equity in directly held stocks, stock mutual funds,
    combination mutual funds, and other mutuals funds ;
    %IF &YEAR LE 2001 %THEN %DO;
       DEQ=STOCKS+STMUTF+.5*COMUTF;
    %END;
    %ELSE %IF &YEAR GE 2004 %THEN %DO;
       DEQ=STOCKS+STMUTF+.5*COMUTF+OMUTF;
    %END;
*   equity held in savings accounts such as 529s, Coverdells or other
    types with investment choice;
    %IF (&YEAR GE 2007) %THEN %DO;
      SAVEQ=X3730*((X7074=1)+(X7074 IN(3))*(MAX(0,X7075)/10000))+
          X3736*((X7077=1)+(X7077 IN(3))*(MAX(0,X7078)/10000))+
          X3742*((X7080=1)+(X7080 IN(3))*(MAX(0,X7081)/10000))+
          X3748*((X7083=1)+(X7083 IN(3))*(MAX(0,X7084)/10000))+
          X3754*((X7086=1)+(X7086 IN(3))*(MAX(0,X7087)/10000))+
          X3760*((X7089=1)+(X7089 IN(3))*(MAX(0,X7090)/10000));
    %END;
    %ELSE %DO;
      SAVEQ=0;
    %END;
*   equity in quasi-liquid retirement assets;
    %IF (&YEAR LE 1998) %THEN %DO;
      RETEQ=IRAKH*((X3631=2)+.5*(X3631 IN (5 6))+.3*(X3631=4))+PENEQ; 
    %END;
    %ELSE %IF (&YEAR EQ 2001) %THEN %DO;
      RETEQ=IRAKH*((X3631=2)+.5*(X3631 IN (5 6))+.3*(X3631=4))+PENEQ
      +X6462*((X6463=1)+.5*(X6463=3))+X6467*((X6468=1)+.5*(X6468=3))
      +X6472*((X6473=1)+.5*(X6473=3))+X6477*((X6478=1)+.5*(X6478=3))
      +X6482*((X6483=1)+.5*(X6483=3))+X6487*((X6488=1)+.5*(X6488=3))
      +X5604*((X6491=1)+.5*(X6491=3))+X5612*((X6492=1)+.5*(X6492=3))
      +X5620*((X6493=1)+.5*(X6493=3))+X5628*((X6494=1)+.5*(X6494=3))
      +X5636*((X6495=1)+.5*(X6495=3))+X5644*((X6496=1)+.5*(X6496=3)); 
    %END;
    %ELSE %IF (&YEAR GE 2010) %THEN %DO;
      RETEQ=(SUM(0,X6551,X6552,X6553,X6554)*((X6555=1)
        +(X6555 IN(3 30))*(MAX(0,X6556)/10000))+
        SUM(0,X6559,X6560,X6561,X6562)*((X6563=1)
        +(X6563 IN(3 30))*(MAX(0,X6564)/10000))+
        SUM(0,X6567,X6568,X6569,X6570)*((X6571=1)
        +(X6571 IN(3 30))*(MAX(0,X6572)/10000)))+
        PENEQ+(X6461=1)*X6462*((X6933=1)+(X6933 IN(3 30))*(MAX(0,X6934)/10000))+
        (X6466=1)*X6467*((X6937=1)+(X6937 IN(3 30))*(MAX(0,X6938)/10000))+
        (X6471=1)*X6472*((X6941=1)+(X6941 IN(3 30))*(MAX(0,X6942)/10000))+
        (X6476=1)*X6477*((X6945=1)+(X6945 IN(3 30))*(MAX(0,X6946)/10000))+
        X5604*((X6962=1)+(X6962 IN(3 30))*(MAX(0,X6963)/10000))+
        X5612*((X6968=1)+(X6968 IN(3 30))*(MAX(0,X6969)/10000))+
        X5620*((X6974=1)+(X6974 IN(3 30))*(MAX(0,X6975)/10000))+
        X5628*((X6980=1)+(X6980 IN(3 30))*(MAX(0,X6981)/10000));
    %END;
    %ELSE %DO;
      RETEQ=(SUM(0,X6551,X6552,X6553,X6554)*((X6555=1)+(X6555=3)*X6556/10000)+
        SUM(0,X6559,X6560,X6561,X6562)*((X6563=1)+(X6563=3)*X6564/10000)+
        SUM(0,X6567,X6568,X6569,X6570)*((X6571=1)+(X6571=3)*X6572/10000))+
        PENEQ+(X6461=1)*X6462*((X6933=1)+(X6933=3)*X6934/10000)+
        (X6466=1)*X6467*((X6937=1)+(X6937=3)*X6938/10000)+
        (X6471=1)*X6472*((X6941=1)+(X6941=3)*X6942/10000)+
        (X6476=1)*X6477*((X6945=1)+(X6945=3)*X6946/10000)+
        (X6481=1)*X6482*((X6949=1)+(X6949=3)*X6950/10000)+
        (X6486=1)*X6487*((X6953=1)+(X6953=3)*X6954/10000)+
        X5604*((X6962=1)+(X6962=3)*X6963/10000)+
        X5612*((X6968=1)+(X6968=3)*X6969/10000)+
        X5620*((X6974=1)+(X6974=3)*X6975/10000)+
        X5628*((X6980=1)+(X6980=3)*X6981/10000)+
        X5636*((X6986=1)+(X6986=3)*X6987/10000)+
        X5644*((X6992=1)+(X6992=3)*X6993/10000);
    %END;
*   ratio of equity to normal income;
    EQUITINC=.;
    %IF (&YEAR GE 1995) %THEN %DO;
      EQUITINC=EQUITY/MAX(100,NORMINC);
    %END;
    
*   brokerage account info;
*   have a brokerage account;
    HBROK=(X3923=1);
*   traded in the past year;
    HTRAD=(X3928>0);
*   number of trades per year; 
    %IF (&YEAR GE 1995) %THEN %DO;
      IF (X7193 NOT IN (0 -1 1 2 3 4 5 6 8 11 12 18 25)) THEN PUT
        "WARNING: UNRECOGNIZED FREQUENCY FOR X7193! " &ID= X7193=;
      PTRAD=((X7193=1)*250+(X7193=2)*52+(X7193=3)*26+(X7193=4)*12
       +(X7193=5)*4+(X7193=6)+(X7193=8)+(X7193=11)*2 +
        (X7193=12)*6+(X7193=18)*8*250+(X7193=25)/2);
      NTRAD=PTRAD*MAX(0,X3928); 
    %END;
    %ELSE %DO;
      NTRAD=MAX(0,X3928); 
    %END;
    
*   total financial assets;
    FIN=LIQ+CDS+NMMF+STOCKS+BOND+RETQLIQ+SAVBND+CASHLI+OTHMA+OTHFIN;
*   have any financial assets: 1=yes, 0=no;
    HFIN=(FIN>0);

***************************************************************************;
*   nonfinancial assets and related variables;
    
*   value of all vehicles (includes autos, motor homes, RVs, airplanes,
    boats);
    %IF (&YEAR GE 1995) %THEN %DO;
      VEHIC=MAX(0,X8166)+MAX(0,X8167)+MAX(0,X8168)+MAX(0,X8188)+
        MAX(0,X2422)+MAX(0,X2506)+MAX(0,X2606)+MAX(0,X2623);
    %END;
    %ELSE %DO;
      VEHIC=MAX(0,X8166)+MAX(0,X8167)+MAX(0,X8168)+
        MAX(0,X2422)+MAX(0,X2506)+MAX(0,X2606)+MAX(0,X2623);
    %END;
*   have any vehicles: 1=yes, 0=no;
    HVEHIC=(VEHIC>0);
*   vehicle supplied by a business;
*   have a business vehicle: 1=yes, 0=no;
    BUSVEH=(X2501=1);
*   number of business vehicles;
    NBUSVEH=X2502;
*   owned vehicles (excludes motorcycles, RVs, motor homes,
    tractors, snow blowers etc);
*   have an owned vehicle: 1=yes, 0=no;
    OWN=(X2201=1);
*   number of owned vehicles;
    NOWN=X2202;
*   value of owned vehicles;
    %IF (&YEAR GE 1995) %THEN %DO;
      VOWN=X8166+X8167+X8168+X8188;
    %END;
    %ELSE %DO;
      VOWN=X8166+X8167+X8168;
    %END;
*   leased vehicles;
*   have leased vehicle: 1=yes, 0=no;
    LEASE=(X2101=1);
*   number of leased vehicles;
    NLEASE=X2102;
*   value of leased vehicles;
    VLEASE=X8163+X8164;
*   total number of vehicles (owned and leased);
    NVEHIC=NOWN+NLEASE; 
*   new model-year car (owned or leased); 
*   NEWCAR1: 1=number of car/truck/SUV with model year no older than
    two years before the survey,
*   NEWCAR2: 1=number of car/truck/SUV with model year no older than
    one year before the survey;
    %IF (&YEAR GE 1995) %THEN %DO;
      NEWCAR1=(X2104>=%EVAL(&YEAR-2))+(X2111>=%EVAL(&YEAR-2))+
        (X2205>=%EVAL(&YEAR-2))+(X2305>=%EVAL(&YEAR-2))+
        (X2405>=%EVAL(&YEAR-2))+(X7152>=%EVAL(&YEAR-2));
      NEWCAR2=(X2104>=%EVAL(&YEAR-1))+(X2111>=%EVAL(&YEAR-1))+
        (X2205>=%EVAL(&YEAR-1))+(X2305>=%EVAL(&YEAR-1))+
        (X2405>=%EVAL(&YEAR-1))+(X7152>=%EVAL(&YEAR-1));
    %END;
    %ELSE %DO;
      NEWCAR1=(X2104>=%EVAL(&YEAR-2))+(X2111>=%EVAL(&YEAR-2))+
        (X2205>=%EVAL(&YEAR-2))+(X2305>=%EVAL(&YEAR-2))+
        (X2405>=%EVAL(&YEAR-2));
      NEWCAR2=(X2104>=%EVAL(&YEAR-1))+(X2111>=%EVAL(&YEAR-1))+
        (X2205>=%EVAL(&YEAR-1))+(X2305>=%EVAL(&YEAR-1))+
        (X2405>=%EVAL(&YEAR-1));
    %END;
    
*   primary residence;
*   for farmers, assume X507 (percent of farm used for
    farming/ranching) is maxed at 90%;
    IF (X507>9000) THEN X507=9000;
*   compute value of business part of farm net of outstanding mortgages;
    FARMBUS=0;
    IF (X507>0) THEN DO;
      FARMBUS=(X507/10000)*(X513+X526-X805-X905-X1005);
      X805=X805*((10000-X507)/10000);
      X808=X808*((10000-X507)/10000);
      X813=X813*((10000-X507)/10000);
      X905=X905*((10000-X507)/10000);
      X908=X908*((10000-X507)/10000);
      X913=X913*((10000-X507)/10000);
      X1005=X1005*((10000-X507)/10000);
      X1008=X1008*((10000-X507)/10000);
      X1013=X1013*((10000-X507)/10000);
      IF (X1103=1) THEN DO;
        FARMBUS=FARMBUS-X1108*(X507/10000);
        X1108=X1108*((10000-X507)/10000);
        X1109=X1109*((10000-X507)/10000);
      END;
      IF (X1114=1) THEN DO;
        FARMBUS=FARMBUS-X1119*(X507/10000);
        X1119=X1119*((10000-X507)/10000);
        X1120=X1120*((10000-X507)/10000);
      END;
      IF (X1125=1) THEN DO;
        FARMBUS=FARMBUS-X1130*(X507/10000);
        X1130=X1130*((10000-X507)/10000);
        X1131=X1131*((10000-X507)/10000);
      END;
      IF (X1136>0 & (X1108+X1119+X1130 > 0)) THEN DO;
        FARMBUS=FARMBUS-X1136*(X507/10000)*((X1108*(X1103=1)
        +X1119*(X1114=1)+X1130*(X1125=1))/(X1108+X1119+X1130));
        X1136=X1136*((10000-X507)/10000)*((X1108*(X1103=1)
        +X1119*(X1114=1)+X1130*(X1125=1))/(X1108+X1119+X1130));
      END;
    END;

*NOTE: needed to add END to IF (X1136>0) statement;

*   value of primary residence;
    HOUSES=SUM(0,X604,X614,X623,X716)+((10000-MAX(0,X507))/10000)*(X513+X526);
*   NOTE: if R only owns a part of the property, the values reported
    should be only Rs share;
*   have owned principal residence: 1=yes, 0=no;
    HHOUSES=(HOUSES^=0);

*   homeownership class: 1=owns ranch/farm/mobile home/house/condo/
    coop/etc., 2=otherwise;
    %IF (&YEAR GE 1995) %THEN %DO;
      IF (X508 IN (1 2)|X601 IN (1 2 3)|X701 IN (1 3 4 5 6 8)|
        (X701=-7 & X7133=1)) THEN HOUSECL=1;
    %END;
    %ELSE %DO;
      IF (X508 IN (1 2)|X601 IN (1 2 3)|X701 IN (1 3 4 5 6)) THEN
        HOUSECL=1;
    %END;
    ELSE HOUSECL=2; 
    
*   other residential real estate: includes land contracts/notes
    household has made, properties other than the principal
    residence that are coded as 1-4 family residences, time shares,
    and vacations homes;
*   JXB - no more x1405,x1409.c1505.x1509,x1619 in 2013. Replaced 
    with x1306,x1310,x1325,x1329,x1339;
    %IF (&YEAR GE 2013) %THEN %DO;
    ORESRE=MAX(X1306,X1310)+MAX(X1325,X1329)+MAX(0,X1339) 
      +(X1703 IN (12,14,21,22,25,40,41,42,43,44,49,50,52,999))*
      MAX(0,X1706)*(X1705/10000)
      +(X1803 IN (12,14,21,22,25,40,41,42,43,44,49,50,52,999))*
      MAX(0,X1806)*(X1805/10000)
      +MAX(0,X2002);
    %END;
    %ELSE %IF (&YEAR EQ 2010) %THEN %DO;
    ORESRE=MAX(X1405,X1409)+MAX(X1505,X1509)+MAX(0,X1619) 
      +(X1703 IN (12,14,21,22,25,40,41,42,43,44,49,50,52,999))*
      MAX(0,X1706)*(X1705/10000)
      +(X1803 IN (12,14,21,22,25,40,41,42,43,44,49,50,52,999))*
      MAX(0,X1806)*(X1805/10000)
      +MAX(0,X2002);
    %END;
    %ELSE %DO;
    ORESRE=MAX(X1405,X1409)+MAX(X1505,X1509)+MAX(X1605,X1609)+MAX(0,X1619) 
      +(X1703 IN (12,14,21,22,25,40,41,42,43,44,49,50,52,999))*
      MAX(0,X1706)*(X1705/10000)
      +(X1803 IN (12,14,21,22,25,40,41,42,43,44,49,50,52,999))*
      MAX(0,X1806)*(X1805/10000)
      +(X1903 IN (12,14,21,22,25,40,41,42,43,44,49,50,52,999))*
      MAX(0,X1906)*(X1905/10000)
      +MAX(0,X2002);
    %END;
*   have other residential real estate: 1=yes, 0=no;
    HORESRE=(ORESRE>0);
    
*   net equity in nonresidential real estate: real estate other than
    the principal residence, properties coded as 1-4 family
    residences, time shares, and vacation homes net of mortgages and
    other loans taken out for investment real estate;
    %IF (&YEAR GE 2010) %THEN %DO;
    NNRESRE =(X1703 IN (1 2 3 4 5 6 7 10 11 13 15 24 45 46 47 48 51 53 -7))*
      MAX(0,X1706)*(X1705/10000)
      +(X1803 IN (1 2 3 4 5 6 7 10 11 13 15 24 45 46 47 48 51 53 -7))*
      MAX(0,X1806)*(X1805/10000)
      +MAX(0,X2012)
      -(X1703 IN (1 2 3 4 5 6 7 10 11 13 15 24 45 46 47 48 51 53 -7))*
      X1715*(X1705/10000) 
      -(X1803 IN (1 2 3 4 5 6 7 10 11 13 15 24 45 46 47 48 51 53 -7))*
      X1815*(X1805/10000)
      -X2016;
    %END;  
    %ELSE %DO;
    NNRESRE =(X1703 IN (1 2 3 4 5 6 7 10 11 13 15 24 45 46 47 48 51 53 -7))*
      MAX(0,X1706)*(X1705/10000)
      +(X1803 IN (1 2 3 4 5 6 7 10 11 13 15 24 45 46 47 48 51 53 -7))*
      MAX(0,X1806)*(X1805/10000)
      +(X1903 IN (1 2 3 4 5 6 7 10 11 13 15 24 45 46 47 48 51 53 -7))*
      MAX(0,X1906)*(X1905/10000)+MAX(0,X2012)
      -(X1703 IN (1 2 3 4 5 6 7 10 11 13 15 24 45 46 47 48 51 53 -7))*
      X1715*(X1705/10000) 
      -(X1803 IN (1 2 3 4 5 6 7 10 11 13 15 24 45 46 47 48 51 53 -7))*
      X1815*(X1805/10000)
      -(X1903 IN (1 2 3 4 5 6 7 10 11 13 15 24 45 46 47 48 51 53 -7))*
      X1915*(X1905/10000)-X2016;
    %END;
*   remove installment loans for PURPOSE=78 from NNRESRE only
    where such property exists--otherwise, if ORESRE exists, include
    loan as RESDBT---otherwise, treat as installment loan;
    IF (NNRESRE^=0) THEN DO;
      FLAG781=1;
      NNRESRE=NNRESRE-X2723*(X2710=78)-X2740*(X2727=78)-X2823*(X2810=78)
        -X2840*(X2827=78)-X2923*(X2910=78)-X2940*(X2927=78);
    END;
    ELSE FLAG781=0;
*   have nonresidential real estate assets: 1=yes, 0=no;
    HNNRESRE=(NNRESRE ^=0); 
    
*   business interests;
*   for businesses where the HH has an active interest, value is net
    equity if business were sold today, plus loans from HH to
    business, minus loans from business to HH not previously
    reported, plus value of personal assets used as collateral for
    business loans that were reported earlier;
*   for businesses where the HH does not have an active interest,
    market value of the interest;
    %IF (&YEAR GE 2010) %THEN %DO;
      BUS=MAX(0,X3129)+MAX(0,X3124)-MAX(0,X3126)*(X3127=5)+
        MAX(0,X3121)*(X3122 IN (1 6))+
        MAX(0,X3229)+MAX(0,X3224)-MAX(0,X3226)*(X3227=5)+
        MAX(0,X3221)*(X3222 IN (1 6))+
        MAX(0,X3335)+FARMBUS+
        MAX(0,X3408)+MAX(0,X3412)+MAX(0,X3416)+MAX(0,X3420)+
        MAX(0,X3452)+MAX(0,X3428);
      ACTBUS=MAX(0,X3129)+MAX(0,X3124)-MAX(0,X3126)*(X3127=5)+
        MAX(0,X3121)*(X3122 IN (1 6))+
        MAX(0,X3229)+MAX(0,X3224)-MAX(0,X3226)*(X3227=5)+
        MAX(0,X3221)*(X3222 IN (1 6))+
        MAX(0,X3335)+FARMBUS;
      NONACTBUS=MAX(0,X3408)+MAX(0,X3412)+MAX(0,X3416)+MAX(0,X3420)+
        MAX(0,X3452)+MAX(0,X3428);
    %END;
    %ELSE %DO;
      BUS=MAX(0,X3129)+MAX(0,X3124)-MAX(0,X3126)*(X3127=5)+
        MAX(0,X3121)*(X3122 IN (1 6))+
        MAX(0,X3229)+MAX(0,X3224)-MAX(0,X3226)*(X3227=5)+
        MAX(0,X3221)*(X3222 IN (1 6))+
        MAX(0,X3329)+MAX(0,X3324)-MAX(0,X3326)*(X3327=5)+
        MAX(0,X3321)*(X3322 IN (1 6))+
        MAX(0,X3335)+FARMBUS+
        MAX(0,X3408)+MAX(0,X3412)+MAX(0,X3416)+MAX(0,X3420)+
        MAX(0,X3424)+MAX(0,X3428);
      ACTBUS=MAX(0,X3129)+MAX(0,X3124)-MAX(0,X3126)*(X3127=5)+
        MAX(0,X3121)*(X3122 IN (1 6))+
        MAX(0,X3229)+MAX(0,X3224)-MAX(0,X3226)*(X3227=5)+
        MAX(0,X3221)*(X3222 IN (1 6))+
        MAX(0,X3329)+MAX(0,X3324)-MAX(0,X3326)*(X3327=5)+
        MAX(0,X3321)*(X3322 IN (1 6))+
        MAX(0,X3335)+FARMBUS;
      NONACTBUS=MAX(0,X3408)+MAX(0,X3412)+MAX(0,X3416)+MAX(0,X3420)+
        MAX(0,X3424)+MAX(0,X3428);
    %END;
*   have business assets: 1=yes, 0=no;
*   definition keys off of ownership question rather than BUS variable
    to deal with businesses with a zero equity value; 
    HBUS=(X3103=1 | X3401=1);
    
*   other nonfinancial assets: defined as total value of
    miscellaneous assets minus other financial assets:
    includes gold, silver (incl. silverware), other metals or metals
    NA type, jewelry, gem stones (incl. antique), cars (antique or
    classic), antiques, furniture, art objects, paintings,
    sculpture, textile art, ceramic art, photographs, (rare) books,
    coin collections, stamp collections, guns, misc. real estate
    (exc. cemetery), cemetery plots, china, figurines,
    crystal/glassware, musical instruments, livestock, horses,
    crops, oriental rugs, furs, other collections, incl. baseball
    cards, records, wine, oil/gas/mineral leases or investments,
    computer, equipment/tools, association or exchange membership,
    and other miscellaneous assets;
    OTHNFIN=X4022+X4026+X4030-OTHFIN+X4018;
*   have other nonfinancial assets: 1=yes, 0=no;
    HOTHNFIN=(OTHNFIN>0);
    
*   total nonfinancial assets;
    NFIN=VEHIC+HOUSES+ORESRE+NNRESRE+BUS+OTHNFIN;
*   have any nonfinancial assets: 1=yes, 0=no;
    HNFIN=(NFIN^=0);
*   total nonfinancial assets excluding principal residences;
    NHNFIN=NFIN-HOUSES;
  
***************************************************************************;
*   total assets;
    ASSET=FIN+NFIN;
*   have any assets: 1=yes, 0=no;
    HASSET=(ASSET^=0);

***************************************************************************;
*   debts and related variables;
    
*   housing debt (mortgage, home equity loans and HELOCs --
    mopup LOCs divided between HE and other);
    IF (X1108+X1119+X1130)>=1 THEN DO;
      HELOC=X1108*(X1103=1)+X1119*(X1114=1)+X1130*(X1125=1)
        +MAX(0,X1136)*(X1108*(X1103=1)
        +X1119*(X1114=1)+X1130*(X1125=1))/(X1108+X1119+X1130);
      MRTHEL=X805+X905+X1005+
        X1108*(X1103=1)+X1119*(X1114=1)+X1130*(X1125=1)
        +MAX(0,X1136)*(X1108*(X1103=1)
        +X1119*(X1114=1)+X1130*(X1125=1))/(X1108+X1119+X1130);
      NH_MORT=MRTHEL-HELOC;
    END;
    ELSE DO;
      HELOC=0;
      MRTHEL=X805+X905+X1005+.5*(MAX(0,X1136))*(HOUSES>0);
      NH_MORT=MRTHEL-HELOC;
    END;

*   Home equity equals home value less all home secured debt;
    HOMEEQ=HOUSES-MRTHEL;

*   have principal residence debt: 1=yes, 0=no;
    HMRTHEL=(MRTHEL>0);
    HHELOC=(HELOC>0); 
    HNH_MORT=(NH_MORT>0);

*   have principal residence debt by type: 1=yes, 0=no;
*   have first-lien mortgage;
    HPRIM_MORT=(X805>0); 
*   have purchase loan - first mortgage; 
*   - note: X7137 only available from 1995 survey forward;
    %IF (&YEAR GE 1995) %THEN %DO;
      PURCH1=((X802>0 & X7137=0) | X7137=8); 
    %END;
*   refinanced;
    %IF (&YEAR GE 1995) %THEN %DO;
      REFIN_EVER=(X7137>0 & X7137^=8); 
    %END;
*   extracted equity from refinance;
    %IF (&YEAR GE 1995) %THEN %DO;
      HEXTRACT_EVER=(X7137 IN(2 3 4));
    %END;
*   have second/third mortgage;
    HSEC_MORT=((X905+X1005)>0);
*   have purchase loan - second/third mortgage; 
    PURCH2=(X918=1 | X1018=1); 
*   have loan used for other purpose;
    HMORT2=((X918^=0 & X918^=1) | (X1018^=0 & X1018^=1));
*   have a  HELOC;
    HELOC_YN=(X1103=1 | X1114=1 | X1125=1);

*   other lines of credit;
    IF (X1108+X1119+X1130)>=1 THEN DO;
      OTHLOC=X1108*(X1103^=1)+X1119*(X1114^=1)+X1130*(X1125^=1)+
        MAX(0,X1136)*(X1108*(X1103^=1)+X1119*(X1114^=1)+
        X1130*(X1125^=1))/(X1108+X1119+X1130);
    END;
    ELSE DO;
      OTHLOC=((HOUSES<=0)+.5*(HOUSES>0))*(MAX(0,X1136));
    END;
*   have balances on lines of credit other than HELOCs: 1=yes, 0=no;
    HOTHLOC=(OTHLOC>0);
    
*   debt for other residential property: includes land contracts,
    residential property other than the principal residence, misc
    vacation, and installment debt reported for cottage/vacation home
    code 67); 
*   NOTE: debt for nonresidential real estate is netted out of the
    corresponding assets; 
    MORT1=(X1703 IN (12,14,21,22,25,40,41,42,43,44,49,50,52,53,999))*
     X1715*(X1705/10000);
    MORT2=(X1803 IN (12,14,21,22,25,40,41,42,43,44,49,50,52,53,999))*
     X1815*(X1805/10000);
    MORT3=0;

*   JXB - in 2013 RESDBT, use 1318 for 1417, 1337 for 1517, 1342 for 1621;
    %IF (&YEAR LE 2007) %THEN %DO;
      MORT3=(X1903 IN (12,14,21,22,25,40,41,42,43,44,49,50,52,53,999))*
        X1915*(X1905/10000);
      RESDBT=X1417+X1517+X1617+X1621+MORT1+MORT2+MORT3+X2006;
    %END;
    %ELSE %IF (&YEAR EQ 2010) %THEN %DO;
      RESDBT=X1417+X1517+X1621+MORT1+MORT2+X2006;
    %END;
    %ELSE %IF (&YEAR EQ 2013) %THEN %DO;
      RESDBT=X1318+X1337+X1342+MORT1+MORT2+X2006;
    %END;
*   see note above at definition of NNRESRE;
    IF (FLAG781^=1 & ORESRE>0) THEN DO;
      FLAG782=1;
      RESDBT=RESDBT+X2723*(X2710=78)+X2740*(X2727=78) 
       +X2823*(X2810=78)+X2840*(X2827=78)
       +X2923*(X2910=78)+X2940*(X2927=78);
    END;
    ELSE FLAG782=0;
*   for parallel treatment, only inlcude PURPOSE=67 where
    ORESRE>0--otherwise, treat as installment loan;
    IF (ORESRE>0) THEN DO;
      FLAG67=1;
      RESDBT=RESDBT+X2723*(X2710=67)+X2740*(X2727=67) 
       +X2823*(X2810=67)+X2840*(X2827=67)
       +X2923*(X2910=67)+X2940*(X2927=67);
    END;
    ELSE FLAG67=0;

*   have other residential real estate debt: 1=yes, 0=no;
    HRESDBT=(RESDBT>0);
    
*   credit card debt;
*   NOTE: from 1992 forward, specific question addresses revolving
    debt at stores, and this amount is treated as credit card debt here;
*   convenience use of credit cards - NOCCBAL, excludes charge
    accounts at stores; 
   %IF (&YEAR GE 2010) %THEN %DO;
      CCBAL = MAX(0,X427)+MAX(0,X413)+MAX(0,X421)+MAX(0,X430)+
        +MAX(0,X7575);
      NOCCBAL=((MAX(0,X427)+MAX(0,X413)+MAX(0,X421)+MAX(0,X430))=0);
    %END;
    %ELSE %IF (&YEAR GE 1992) %THEN %DO;
      CCBAL = MAX(0,X427)+MAX(0,X413)+MAX(0,X421)+MAX(0,X430)+
        MAX(0,X424)+MAX(0,X7575);
      NOCCBAL=((MAX(0,X427)+MAX(0,X413)+MAX(0,X421)+MAX(0,X430)+MAX(0,X424))=0);
    %END;
    %ELSE %DO;
      CCBAL = MAX(0,X427)+MAX(0,X413)+MAX(0,X421)+MAX(0,X430)+
        MAX(0,X424);
      NOCCBAL=((MAX(0,X427)+MAX(0,X413)+MAX(0,X421)+MAX(0,X430)+MAX(0,X424))=0);
    %END;
*   have credit card balances: 1=yes, 0=no;
    HCCBAL=(CCBAL>0);
    
*   installment loans not classified elsewhere;
*   subdivide into vehicle loans, education loans, and other
    installment loans;
    %IF (&YEAR GE 1995) %THEN %DO;
      VEH_INST=X2218+X2318+X2418+X7169+X2424+X2519+X2619+X2625;
      EDN_INST=X7824+X7847+X7870+X7924+X7947+X7970+X7179+
        X2723*(X2710=83)+X2740*(X2727=83)+X2823*(X2810=83)+
        X2840*(X2827=83)+X2923*(X2910=83)+X2940*(X2927=83);
      INSTALL=X2218+X2318+X2418+X7169+X2424+X2519+X2619+X2625+X7183
        +X7824+X7847+X7870+X7924
        +X7947+X7970+X7179
        +X1044+X1215+X1219;
    %END;
    %ELSE %IF (&YEAR EQ 1992) %THEN %DO;
      VEH_INST=X2218+X2318+X2418+X2424+X2519+X2619+X2625;
      EDN_INST=X7824+X7847+X7870+X7924+X7947+X7970+
        X2723*(X2710=83)+X2740*(X2727=83)+X2823*(X2810=83)+
        X2840*(X2827=83)+X2923*(X2910=83)+X2940*(X2927=83);
      INSTALL=X2218+X2318+X2418+X2424+X2519+X2619+X2625
        +X7824+X7847+X7870+X7924+X7947+X7970
        +X1044+X1215+X1219;
    %END;
    %ELSE %DO;
      VEH_INST=X2218+X2318+X2418+X2424+X2519+X2619+X2625;
      EDN_INST=X2723*(X2710=83)+X2740*(X2727=83)+
        X2823*(X2810=83)+X2840*(X2827=83)+X2923*(X2910=83)+
        X2940*(X2927=83);
      INSTALL=X2218+X2318+X2418+X2424+X2519+X2619+X2625
        +X1044+X1215+X1219;
    %END;
*   see notes above at definitions of NNRESRE and RESDBT;
    IF (FLAG781=0 & FLAG782=0) THEN DO;
      INSTALL=INSTALL+X2723*(X2710=78)+X2740*(X2727=78)+
        X2823*(X2810=78)+X2840*(X2827=78)+X2923*(X2910=78)+
        X2940*(X2927=78);
    END;
    IF (FLAG67=0) THEN DO;
      INSTALL=INSTALL+X2723*(X2710=67)+X2740*(X2727=67)+
        X2823*(X2810=67)+X2840*(X2827=67)+X2923*(X2910=67)+
        X2940*(X2927=67);
    END;
    INSTALL=INSTALL+X2723*(X2710 NOT IN (67 78))+X2740*(X2727 NOT IN
      (67 78))+X2823*(X2810 NOT IN (67 78))+X2840*(X2827 NOT IN
      (67 78))+X2923*(X2910 NOT IN (67 78))+X2940*(X2927 NOT IN (67 78));
    OTH_INST=INSTALL-VEH_INST-EDN_INST;
    HVEH_INST=(VEH_INST>0);
    HEDN_INST=(EDN_INST>0);
    HOTH_INST=(OTH_INST>0);

*   have any installment debt: 1=yes, 0=no;
    HINSTALL=(INSTALL>0);
    
*   margin loans; 
*   except in 1995, the SCF does not ask whether the margin loan
    was reported earlier: the instruction explicitly excludes loans
    reported earlier;
    %IF (&YEAR EQ 1995) %THEN %DO;
      OUTMARG=MAX(0,X3932)*(X7194=5);
    %END;
    %ELSE %DO;
      OUTMARG=MAX(0,X3932);
    %END;
    
*   pension loans not reported earlier;
    %IF (&YEAR GE 2010) %THEN %DO;
      OUTPEN1=MAX(0,X11027)*(X11070=5);
      OUTPEN2=MAX(0,X11127)*(X11170=5);
      OUTPEN4=MAX(0,X11327)*(X11370=5);
      OUTPEN5=MAX(0,X11427)*(X11470=5);
      OUTPEN3=0;
      OUTPEN6=0;
    %END;
    %ELSE %IF (&YEAR GE 2004) %THEN %DO;
      OUTPEN1=MAX(0,X11027)*(X11070=5);
      OUTPEN2=MAX(0,X11127)*(X11170=5);
      OUTPEN3=MAX(0,X11227)*(X11270=5);
      OUTPEN4=MAX(0,X11327)*(X11370=5);
      OUTPEN5=MAX(0,X11427)*(X11470=5);
      OUTPEN6=MAX(0,X11527)*(X11570=5);
    %END;
    %ELSE %DO;
      OUTPEN1=MAX(0,X4229)*(X4230=5);
      OUTPEN2=MAX(0,X4329)*(X4330=5);
      OUTPEN3=MAX(0,X4429)*(X4430=5);
      OUTPEN4=MAX(0,X4829)*(X4830=5);
      OUTPEN5=MAX(0,X4929)*(X4930=5);
      OUTPEN6=MAX(0,X5029)*(X5030=5);
    %END;
    OUTPEN=MAX(0,SUM(OUTPEN1,OUTPEN2,OUTPEN3,OUTPEN4,OUTPEN5,OUTPEN6));
    
*   other debts (loans against pensions, loans against life insurance,
    margin loans, miscellaneous);
    %IF (&YEAR GE 2010) %THEN %DO;
    ODEBT=OUTPEN1+OUTPEN2+OUTPEN4+OUTPEN5
      +MAX(0,X4010)+MAX(0,X4032)+OUTMARG;
    %END;
    %ELSE %DO;
    ODEBT=OUTPEN1+OUTPEN2+OUTPEN3+OUTPEN4+OUTPEN5+OUTPEN6
      +MAX(0,X4010)+MAX(0,X4032)+OUTMARG;
    %END;
*   have any other debts: 1=yes, 0=no;
    HODEBT=(ODEBT>0);
    
*   total debt;
    DEBT=MRTHEL+RESDBT+OTHLOC+CCBAL+INSTALL+ODEBT;
*   have any debts: 1=yes, 0=no;
    HDEBT=(DEBT>0);

***************************************************************************;
*   total net worth;
    NETWORTH=ASSET-DEBT;
    IF (NETWORTH<=.Z) THEN PUT Y1= &PID= FIN= NFIN= DEBT= LIQ= CDS= NMMF=
      STOCKS=  BOND=  RETQLIQ=  SAVBND=  CASHLI=  OTHMA=  OTHFIN= 
      VEHIC= HOUSES= ORESRE= NNRESRE= BUS= OTHNFIN=
      MRTHEL= RESDBT= OTHLOC= CCBAL= INSTALL= ODEBT=;

* levarage ratio;
  IF (DEBT >0 & ASSET > 0) THEN LEVRATIO=(DEBT/ASSET);
  ELSE IF (DEBT > 0 & ASSET=0) THEN LEVRATIO=1;
  ELSE LEVRATIO=0;

* debt to income ratio, if no income, assign arbitrary value of 10 to
  match how computed in Bulletin article tables;
  IF (DEBT > 0 & INCOME > 0) THEN DEBT2INC=(DEBT/INCOME);
  ELSE IF (DEBT > 0 & INCOME=0) THEN DEBT2INC=10;
  ELSE DEBT2INC=0; 

***************************************************************************;
*   Capital Gains;

*   define capital gains over all assets in SCF where measure is
    possible;
  
*   principal residences: current value less original purchase price
    and less improvements;
    KGHOUSE=MAX(X513,X526,X604,X614,X623,X716)-
      MAX(X607,X617,X627+X631,X635,X717)-X1202;
    
*   other real estate: current value less purchase price adjusted
    for share owned;
    %IF (&YEAR GE 2010) %THEN %DO;
    KGORE=(X1705/10000)*(X1706-X1709)+(X1805/10000)*(X1806-X1809)+
      (X2002-X2003)+(X2012-X2013);
    %END;  
    %ELSE %DO;
    KGORE=(X1705/10000)*(X1706-X1709)+(X1805/10000)*(X1806-X1809)+
      (X1905/10000)*(X1906-X1909)+(X2002-X2003)+(X2012-X2013);
    %END;  
    
*   businesses: current value less tax basis--active and non-active
    businesses;
    %IF (&YEAR GE 2010) %THEN %DO;
    KGBUS=(X3129-X3130)+(X3229-X3230)+(X3408-X3409)+
      (X3412-X3413)+(X3416-X3417)+(X3420-X3421)+(X3452-X3453)+
      (X3428-X3429);
    %END;  
    %ELSE %DO;
    KGBUS=(X3129-X3130)+(X3229-X3230)+(X3329-X3330)+(X3408-X3409)+
      (X3412-X3413)+(X3416-X3417)+(X3420-X3421)+(X3424-X3425)+
      (X3428-X3429);
    %END;

*   adjust for capital gains on farm businesses;
    IF (X507>0) THEN DO;
      FARMBUS_KG=(X507/10000)*KGHOUSE;
      KGHOUSE=KGHOUSE-FARMBUS_KG;
      KGBUS=KGBUS+FARMBUS_KG;
    END;
    ELSE FARMBUS_KG=0;
    
*   stocks and mutual funds: current value less gains/losses;
    KGSTMF=(X3918-X3920)+(X3833-X3835);
    
*   total gains/losses;
    KGTOTAL=KGHOUSE+KGORE+KGBUS+KGSTMF;
  

*   whether have capital gains ;
    ARRAY HAVEO{*} HNETWRTH HKGTOTAL HKGHOUSE HKGORE HKGBUS HKGSTMF;
    ARRAY AMTO{*}  NETWORTH KGTOTAL KGHOUSE KGORE KGBUS KGSTMF;

    DO I=1 TO DIM(AMTO);
      HAVEO{I}=(AMTO{I}^=0);
    END;

***************************************************************************;
***************************************************************************;
*   the following section deals with characteristics of loans:
    compute loan payments and related measures, array loan balances by
    use of the borrowed funds and by type of institution from which
    the loans was obtained;

***************************************************************************;
*   compute payments on all loans on a monthly basis;

*   credit card payments;
*   use HREF assumption of 2.5% per month for credit card payments;
    CCPAY=.025*CCBAL;

*   mortgage payments;
    IF (X808>0) THEN PAYMORT1=X808*(%MCONV(F=X809));
    ELSE IF (X813>0) THEN PAYMORT1=X813*(%MCONV(F=X814));
    ELSE PAYMORT1=0;
    IF (X908>0) THEN PAYMORT2=X908*(%MCONV(F=X909));
    ELSE IF (X913>0) THEN PAYMORT2=X913*(%MCONV(F=X914));
    ELSE PAYMORT2=0;
    IF (X1008>0) THEN PAYMORT3=X1008*(%MCONV(F=X1009));
    ELSE IF (X1013>0) THEN PAYMORT3=X1013*(%MCONV(F=X1014));
    ELSE PAYMORT3=0;
    %IF (&YEAR GE 1992) %THEN %DO;
      IF (X1039>0) THEN PAYMORTO=X1039*(%MCONV(F=X7567));
        ELSE IF (X1040>0) THEN PAYMORTO=X1040*(%MCONV(F=X1041));
        ELSE PAYMORTO=0;
    %END;
    %ELSE %DO;
      IF (X1039>0) THEN PAYMORTO=X1039;
      ELSE IF (X1040>0) THEN PAYMORTO=X1040*(%MCONV(F=X1041));
      ELSE PAYMORTO=0;
    %END;
*   lines of credit;
    IF (X1109>0) THEN PAYLOC1=X1109*(%MCONV(F=X1110));
    ELSE PAYLOC1=0;
    IF (X1120>0) THEN PAYLOC2=X1120*(%MCONV(F=X1121));
    ELSE PAYLOC2=0;
    IF (X1131>0) THEN PAYLOC3=X1131*(%MCONV(F=X1132));
    ELSE PAYLOC3=0;
*   LOC mopup (payment is estimated using the ratio of payment to debt
    on the last LOC in the grid, and is assumed to have the same
    purpose and per);
    IF (X1136>0) THEN DO;
      IF (X1130^=0) THEN DO;
        IF (X1125=1) THEN HMOP=1;
        PTMOP=(MAX(0,X1131)/X1130)*X1136;
        PPMOP=X1132;
        PAYLOCO=PTMOP*(%MCONV(F=PPMOP));
      END;
      ELSE IF (X1119^=0) THEN DO;
        IF (X1114=1) THEN HMOP=1;
        PTMOP=(MAX(0,X1120)/X1119)*X1136;
        PPMOP=X1121;
        PAYLOCO=PTMOP*(%MCONV(F=PPMOP));
      END;
      ELSE IF (X1108^=0) THEN DO;
        IF (X1103=1) THEN HMOP=1;
        PTMOP=(MAX(0,X1109)/X1108)*X1136;
        PPMOP=X1110;
        PAYLOCO=PTMOP*(%MCONV(F=PPMOP));
      END;
*     if 3 LOCs in grid but no amounts outstanding, payment
      rate and per for the mopup are set to median values for all
      LOCs (median payment rate for 1995 is .0316 7/19/96); 
      ELSE DO;
        IF (X1125=1) THEN HMOP=1;
        PTMOP=(.0316)*X1136;
        PPMOP=4;
        PAYLOCO=PTMOP*(%MCONV(F=PPMOP));
      END;  
    END;
    ELSE DO;
      HMOP=0;
      PTMOP=0;
      PPMOP=0;
      PAYLOCO=0;
    END;
    IF PAYLOCO<=.Z THEN PUT 'PAYLOCO PROBLEM ' &ID= X1136= X1125= X1130=
      X1131= X1132= X1114= X1119= X1120= X1121= X1103= X1108= X1109=
      X1110=;
*   home improvement loans;
    %IF (&YEAR GE 1992) %THEN %DO;
      IF (X1210>0) THEN PAYHI1=X1210*(%MCONV(F=X7565));    
      ELSE IF (X1211>0) THEN PAYHI1=X1211*(%MCONV(F=X1212));
      ELSE PAYHI1=0;
    %END;
    %ELSE %DO;
      IF (X1210>0) THEN PAYHI1=X1210;
      ELSE  IF (X1211>0) THEN PAYHI1=X1211*(%MCONV(F=X1212));
      ELSE PAYHI1=0;
    %END;
    IF (X1220>0) THEN PAYHI2=X1220*(%MCONV(F=X1221));    
      ELSE PAYHI2=0;
*   land contracts, assumed to be interest payments on amount
    outstanding, at 8 percent APR;
*   JXB - in 2013 RESDBT, use 1318 for 1417, 1337 for 1517, 1342 for 1621;
    %IF (&YEAR LE 2010) %THEN %DO;
      PAYLC1=(1.08**(1/12)-1)*X1417;
      PAYLC2=(1.08**(1/12)-1)*X1517;
      PAYLCO=(1.08**(1/12)-1)*X1621;
    %END;
    %ELSE %IF (&YEAR GE 2013) %THEN %DO;
      PAYLC1=(1.08**(1/12)-1)*X1318;
      PAYLC2=(1.08**(1/12)-1)*X1337;
      PAYLCO=(1.08**(1/12)-1)*X1342;
    %END;
*   other residential property;
    IF (X1718>0) THEN PAYORE1=
      (X1703 IN (12,14,21,22,25,40,41,42,43,44,49,50,52,999))*
      X1718*(%MCONV(F=X1719))*(X1705/10000);    
    ELSE IF (X1723>0) THEN PAYORE1=
      (X1703 IN (12,14,21,22,25,40,41,42,43,44,49,50,52,999))*
      X1723*(%MCONV(F=X1724))*(X1705/10000);
    ELSE PAYORE1=0;
    IF (X1818>0) THEN PAYORE2=
      (X1803 IN (12,14,21,22,25,40,41,42,43,44,49,50,52,999))*
      X1818*(%MCONV(F=X1819))*(X1805/10000);    
    ELSE IF (X1823>0) THEN PAYORE2=
      (X1803 IN (12,14,21,22,25,40,41,42,43,44,49,50,52,999))*
      X1823*(%MCONV(F=X1824))*(X1805/10000);
    ELSE PAYORE2=0;

    IF (X2007>0) THEN PAYOREV=X2007*(%MCONV(F=X2008));    
    ELSE PAYOREV=0;

    %IF (&YEAR LT 2010) %THEN %DO;
    
      PAYLC3=(1.08**(1/12)-1)*X1617;

      IF (X1918>0) THEN PAYORE3=
        (X1903 IN (12,14,21,22,25,40,41,42,43,44,49,50,52,999))*
        X1918*(%MCONV(F=X1919))*(X1905/10000);    
      ELSE IF (X1923>0) THEN PAYORE3=
        (X1903 IN (12,14,21,22,25,40,41,42,43,44,49,50,52,999))*
        X1923*(%MCONV(F=X1924))*(X1905/10000);
      ELSE PAYORE3=0;
    %END;
    %ELSE %DO;
      PAYORE3=0;
    %END;

*   vehicle loans;
    %IF (&YEAR>=1992) %THEN %DO;
      IF (X2213>0) THEN PAYVEH1=X2213*(%MCONV(F=X7537));    
        ELSE IF (X2214>0) THEN PAYVEH1=X2214*(%MCONV(F=X2215));
        ELSE PAYVEH1=0;
      IF (X2313>0) THEN PAYVEH2=X2313*(%MCONV(F=X7536));    
        ELSE IF (X2314>0) THEN PAYVEH2=X2314*(%MCONV(F=X2315));
       ELSE PAYVEH2=0;
      IF (X2413>0) THEN PAYVEH3=X2413*(%MCONV(F=X7535));    
        ELSE IF (X2414>0) THEN PAYVEH3=X2414*(%MCONV(F=X2415));
        ELSE PAYVEH3=0;
    %END;
    %ELSE %DO;
      IF (X2213>0) THEN PAYVEH1=X2213;    
        ELSE IF (X2214>0) THEN PAYVEH1=X2214*(%MCONV(F=X2215));
        ELSE PAYVEH1=0;
      IF (X2313>0) THEN PAYVEH2=X2313;    
        ELSE IF (X2314>0) THEN PAYVEH2=X2314*(%MCONV(F=X2315));
        ELSE PAYVEH2=0;
      IF (X2413>0) THEN PAYVEH3=X2413;    
        ELSE IF (X2414>0) THEN PAYVEH3=X2414*(%MCONV(F=X2415));
        ELSE PAYVEH3=0;
    %END;
    %IF (&YEAR GE 1995) %THEN %DO;
      IF (X7162>0) THEN PAYVEH4=X7162*(%MCONV(F=X7163));    
        ELSE IF (X7164>0) THEN PAYVEH4=X7164*(%MCONV(F=X7165));
        ELSE PAYVEH4=0;
    %END;
    %ELSE %DO;
      PAYVEH4=0;
    %END;
    IF (X2425>0) THEN PAYVEHM=X2425*(%MCONV(F=X2426));    
      ELSE PAYVEHM=0;
    %IF (&YEAR GE 1992) %THEN %DO;
      IF (X2514>0) THEN PAYVEO1=X2514*(%MCONV(F=X7531));    
        ELSE IF (X2515>0) THEN PAYVEO1=X2515*(%MCONV(F=X2516));
        ELSE PAYVEO1=0;
      IF (X2614>0) THEN PAYVEO2=X2614*(%MCONV(F=X7530));    
        ELSE IF (X2615>0) THEN PAYVEO2=X2615*(%MCONV(F=X2616));
        ELSE PAYVEO2=0;
    %END;
    %ELSE %DO;
      IF (X2514>0) THEN PAYVEO1=X2514;    
        ELSE IF (X2515>0) THEN PAYVEO1=X2515*(%MCONV(F=X2516));
        ELSE PAYVEO1=0;
      IF (X2614>0) THEN PAYVEO2=X2614;    
        ELSE IF (X2615>0) THEN PAYVEO2=X2615*(%MCONV(F=X2616));
        ELSE PAYVEO2=0;
    %END;
    IF (X2626>0) THEN PAYVEOM=X2626*(%MCONV(F=X2627));    
      ELSE PAYVEOM=0;
*   student loans;
    %IF (&YEAR GE 1992) %THEN %DO;
      IF (X7815>0) THEN PAYEDU1=X7815*(%MCONV(F=X7816));    
      ELSE IF (X7817>0) THEN PAYEDU1=X7817*(%MCONV(F=X7818));
      ELSE PAYEDU1=0;
      IF (X7838>0) THEN PAYEDU2=X7838*(%MCONV(F=X7839));    
      ELSE IF (X7840>0) THEN PAYEDU2=X7840*(%MCONV(F=X7841));
      ELSE PAYEDU2=0;
      IF (X7861>0) THEN PAYEDU3=X7861*(%MCONV(F=X7862));    
      ELSE IF (X7863>0) THEN PAYEDU3=X7863*(%MCONV(F=X7864));
      ELSE PAYEDU3=0;
      IF (X7915>0) THEN PAYEDU4=X7915*(%MCONV(F=X7916));    
      ELSE IF (X7917>0) THEN PAYEDU4=X7917*(%MCONV(F=X7918));
      ELSE PAYEDU4=0;
      IF (X7938>0) THEN PAYEDU5=X7938*(%MCONV(F=X7939));    
      ELSE IF (X7940>0) THEN PAYEDU5=X7940*(%MCONV(F=X7941));
      ELSE PAYEDU5=0;
      IF (X7961>0) THEN PAYEDU6=X7961*(%MCONV(F=X7962));    
      ELSE IF (X7963>0) THEN PAYEDU6=X7963*(%MCONV(F=X7964));
      ELSE PAYEDU6=0;
      %IF (&YEAR GE 1995) %THEN %DO;
        IF (X7180>0) THEN PAYEDU7=X7180*(%MCONV(F=X7181));    
          ELSE PAYEDU7=0;
      %END;
      %ELSE %DO;
        PAYEDU7=0;
      %END;
    %END;
    %ELSE %DO;
      PAYEDU1=0;
      PAYEDU2=0;
      PAYEDU3=0;
      PAYEDU4=0;
      PAYEDU5=0;
      PAYEDU6=0;
      PAYEDU7=0;
    %END;
*   installment loans (loans for nonresidential property not included);
    %IF (&YEAR GE 1992) %THEN %DO;
      IF (X2718>0) THEN PAYILN1=(X2710^=78|FLAG781=0)*X2718*
        (%MCONV(F=X7527));
      ELSE IF (X2719>0) THEN PAYILN1=(X2710^=78|FLAG781=0)*X2719*
        (%MCONV(F=X2720));
      ELSE PAYILN1=0;
      IF (X2735>0) THEN PAYILN2=(X2727^=78|FLAG781=0)*X2735*
        (%MCONV(F=X7526));
      ELSE IF (X2736>0) THEN PAYILN2=(X2727^=78|FLAG781=0)*X2736*
        (%MCONV(F=X2737));
      ELSE PAYILN2=0;
      IF (X2818>0) THEN PAYILN3=(X2810^=78|FLAG781=0)*X2818*
        (%MCONV(F=X7525));    
      ELSE IF (X2819>0) THEN PAYILN3=(X2810^=78|FLAG781=0)*X2819*
        (%MCONV(F=X2820));
      ELSE PAYILN3=0;
      IF (X2835>0) THEN PAYILN4=(X2827^=78|FLAG781=0)*X2835*
        (%MCONV(F=X7524));
      ELSE IF (X2836>0) THEN PAYILN4=(X2827^=78|FLAG781=0)*X2836*
        (%MCONV(F=X2837));
      ELSE PAYILN4=0;
      IF (X2918>0) THEN PAYILN5=(X2910^=78|FLAG781=0)*X2918*
        (%MCONV(F=X7523));
      ELSE IF (X2919>0) THEN PAYILN5=(X2910^=78|FLAG781=0)*X2919*
        (%MCONV(F=X2920));
      ELSE PAYILN5=0;
      IF (X2935>0) THEN PAYILN6=(X2927^=78|FLAG781=0)*X2935*
        (%MCONV(F=X7522));
      ELSE IF (X2936>0) THEN PAYILN6=(X2927^=78|FLAG781=0)*X2936*
        (%MCONV(F=X2937));
      ELSE PAYILN6=0;
      %IF (&YEAR GE 1995) %THEN %DO;
        IF (X7184>0) THEN PAYILN7=X7184*(%MCONV(F=X7185));    
          ELSE PAYILN7=0;
      %END;
      %ELSE %DO;
        PAYILN7=0;
      %END;
    %END;
    %ELSE %DO;
      IF (X2718>0) THEN PAYILN1=(X2710^=78|FLAG781=0)*X2718;    
        ELSE IF (X2719>0) THEN PAYILN1=(X2710^=78|FLAG781=0)*X2719*
          (%MCONV(F=X2720));
        ELSE PAYILN1=0;
      IF (X2735>0) THEN PAYILN2=(X2727^=78|FLAG781=0)*X2735;    
        ELSE IF (X2736>0) THEN PAYILN2=(X2727^=78|FLAG781=0)*X2736*
          (%MCONV(F=X2737));
        ELSE PAYILN2=0;
      IF (X2818>0) THEN PAYILN3=(X2810^=78|FLAG781=0)*X2818;    
        ELSE IF (X2819>0) THEN PAYILN3=(X2810^=78|FLAG781=0)*X2819*
          (%MCONV(F=X2820));
        ELSE PAYILN3=0;
      IF (X2835>0) THEN PAYILN4=(X2827^=78|FLAG781=0)*X2835;    
        ELSE IF (X2836>0) THEN PAYILN4=(X2827^=78|FLAG781=0)*X2836*
          (%MCONV(F=X2837));
        ELSE PAYILN4=0;
      IF (X2918>0) THEN PAYILN5=(X2910^=78|FLAG781=0)*X2918;    
        ELSE IF (X2919>0) THEN PAYILN5=(X2910^=78|FLAG781=0)*X2919*
          (%MCONV(F=X2920));
        ELSE PAYILN5=0;
      IF (X2935>0) THEN PAYILN6=(X2927^=78|FLAG781=0)*X2935;    
        ELSE IF (X2936>0) THEN PAYILN6=(X2927^=78)*X2936*(%MCONV(F=X2937));
        ELSE PAYILN6=0;
        PAYILN7=0;
    %END;
*   payments on margin loans, assumed to be interest payments on
    balance outstanding, paid at 8 percent APR;
    PAYMARG=(1.08**(1/12)-1)*OUTMARG;
*   loans against insurance policies;
    PAYINS=MAX(0,X4011)*(%MCONV(F=X4012));
*   payments on other loans, set to zero;
    PAYOTH=0;
*   payments on loans against pension plans not previously reported
    (details available starting in 1998, but still only include loans
    not previously reported); 
    %IF (&YEAR GE 2010) %THEN %DO;
      IF (X11070=5) THEN PAYPEN1=X11028*(%MCONV(F=X11029));
      ELSE PAYPEN1=0;
      IF (X11170=5) THEN PAYPEN2=X11128*(%MCONV(F=X11129));
      ELSE PAYPEN2=0;
      PAYPEN3=0;
      IF (X11370=5) THEN PAYPEN4=X11328*(%MCONV(F=X11329));
      ELSE PAYPEN4=0;
      IF (X11470=5) THEN PAYPEN5=X11428*(%MCONV(F=X11429));
      ELSE PAYPEN5=0;
      PAYPEN6=0;	
    %END;
    %ELSE %IF (&YEAR GE 2004) %THEN %DO;
      IF (X11070=5) THEN PAYPEN1=X11028*(%MCONV(F=X11029));
      ELSE PAYPEN1=0;
      IF (X11170=5) THEN PAYPEN2=X11128*(%MCONV(F=X11129));
      ELSE PAYPEN2=0;
      IF (X11270=5) THEN PAYPEN3=X11228*(%MCONV(F=X11229));
      ELSE PAYPEN3=0;
      IF (X11370=5) THEN PAYPEN4=X11328*(%MCONV(F=X11329));
      ELSE PAYPEN4=0;
      IF (X11470=5) THEN PAYPEN5=X11428*(%MCONV(F=X11429));
      ELSE PAYPEN5=0;
      IF (X11570=5) THEN PAYPEN6=X11528*(%MCONV(F=X11529));
      ELSE PAYPEN6=0;
    %END;
    %ELSE %IF (&YEAR GE 1998 AND &YEAR LE 2001) %THEN %DO;
      IF (X4230=5) THEN PAYPEN1=X7211*(%MCONV(F=X7212));    
      ELSE PAYPEN1=0;
      IF (X4330=5) THEN PAYPEN2=X7220*(%MCONV(F=X7221));    
      ELSE PAYPEN2=0;
      IF (X4430=5) THEN PAYPEN3=X7229*(%MCONV(F=X7230));    
      ELSE PAYPEN3=0;
      IF (X4830=5) THEN PAYPEN4=X7278*(%MCONV(F=X7279));    
      ELSE PAYPEN4=0;
      IF (X4930=5) THEN PAYPEN5=X7287*(%MCONV(F=X7288));    
      ELSE PAYPEN5=0;
      IF (X5030=5) THEN PAYPEN6=X7296*(%MCONV(F=X7297));    
      ELSE PAYPEN6=0;
    %END;
    %ELSE %DO;
      PAYPEN1=(1.065**(1/12)-1)*MAX(0,X4229)*(X4230=5);
      PAYPEN2=(1.065**(1/12)-1)*MAX(0,X4329)*(X4330=5);
      PAYPEN3=(1.065**(1/12)-1)*MAX(0,X4429)*(X4430=5);
      PAYPEN4=(1.065**(1/12)-1)*MAX(0,X4829)*(X4830=5);
      PAYPEN5=(1.065**(1/12)-1)*MAX(0,X4929)*(X4930=5);
      PAYPEN6=(1.065**(1/12)-1)*MAX(0,X5029)*(X5030=5);
    %END;

    %IF (&YEAR GE 2010) %THEN %DO;
      ARRAY PAYMENTS {*} CCPAY PAYMORT1 PAYMORT2 PAYMORT3 PAYMORTO
        PAYLOC1 PAYLOC2 PAYLOC3 PAYLOCO PAYHI1 PAYHI2
        PAYLC1 PAYLC2 PAYLCO PAYORE1 PAYORE2 PAYOREV
        PAYVEH1 PAYVEH2 PAYVEH3 PAYVEH4 PAYVEHM PAYVEO1 PAYVEO2 PAYVEOM   
        PAYEDU1 PAYEDU2 PAYEDU3 PAYEDU4 PAYEDU5 PAYEDU6 PAYEDU7
        PAYILN1 PAYILN2 PAYILN3 PAYILN4 PAYILN5 PAYILN6 PAYILN7
        PAYMARG PAYINS PAYOTH
        PAYPEN1 PAYPEN2 PAYPEN4 PAYPEN5;
    %END;
    %ELSE %DO;
      ARRAY PAYMENTS {*} CCPAY PAYMORT1 PAYMORT2 PAYMORT3 PAYMORTO
        PAYLOC1 PAYLOC2 PAYLOC3 PAYLOCO PAYHI1 PAYHI2
        PAYLC1 PAYLC2 PAYLC3 PAYLCO PAYORE1 PAYORE2 PAYORE3 PAYOREV
        PAYVEH1 PAYVEH2 PAYVEH3 PAYVEH4 PAYVEHM PAYVEO1 PAYVEO2 PAYVEOM   
        PAYEDU1 PAYEDU2 PAYEDU3 PAYEDU4 PAYEDU5 PAYEDU6 PAYEDU7
        PAYILN1 PAYILN2 PAYILN3 PAYILN4 PAYILN5 PAYILN6 PAYILN7
        PAYMARG PAYINS PAYOTH
        PAYPEN1 PAYPEN2 PAYPEN3 PAYPEN4 PAYPEN5 PAYPEN6;
    %END;  

*   compute total monthly payments;
    TPAY=0;
    DO I=1 TO DIM(PAYMENTS);
      TPAY=TPAY+MAX(0,PAYMENTS{I});
      IF PAYMENTS{I}<=.Z THEN PUT 'ERROR? MISSING PAYMENTS AMOUNT ' &ID=
        PAYMENTS{I}=;
    END;

    %IF (&YEAR GE 2010) %THEN %DO;
      MORTPAY=PAYMORT1+PAYMORT2+PAYMORT3+PAYORE1+PAYORE2+
        PAYLC1+PAYLC2+PAYLCO+PAYOREV+PAYLOC1*(X1103=1)+PAYLOC2*(X1114=1)
        +PAYLOC3*(X1125=1);
*   total non-mortgage non-revolving consumer debt;
      CONSPAY=PAYMORTO+PAYVEH1+PAYVEH2+PAYVEH3+PAYVEH4+PAYVEHM+PAYVEO1+
        PAYVEO2+PAYVEOM+PAYEDU1+PAYEDU2+PAYEDU3+PAYEDU4+PAYEDU5+PAYEDU6+
        PAYEDU7+PAYILN1+PAYILN2+PAYILN3+PAYILN4+PAYILN5+PAYILN6+PAYILN7+
        PAYMARG+PAYINS+PAYOTH+PAYPEN1+PAYPEN2+PAYPEN4+PAYPEN5+
        +PAYHI1+PAYHI2;
    %END;	
    %ELSE %DO;
      MORTPAY=PAYMORT1+PAYMORT2+PAYMORT3+PAYORE1+PAYORE2+PAYORE3+
        PAYLC1+PAYLC2+PAYLC3+PAYLCO+PAYOREV+PAYLOC1*(X1103=1)+PAYLOC2*(X1114=1)
        +PAYLOC3*(X1125=1);
*   total non-mortgage non-revolving consumer debt;
      CONSPAY=PAYMORTO+PAYVEH1+PAYVEH2+PAYVEH3+PAYVEH4+PAYVEHM+PAYVEO1+
        PAYVEO2+PAYVEOM+PAYEDU1+PAYEDU2+PAYEDU3+PAYEDU4+PAYEDU5+PAYEDU6+
        PAYEDU7+PAYILN1+PAYILN2+PAYILN3+PAYILN4+PAYILN5+PAYILN6+PAYILN7+
        PAYMARG+PAYINS+PAYOTH+PAYPEN1+PAYPEN2+PAYPEN3+PAYPEN4+PAYPEN5+
        PAYPEN6+PAYHI1+PAYHI2;
    %END;

*   total revolving debt (excluding HELOCs);
    REVPAY=CCPAY+PAYLOC1*(X1103 IN(0 5))+PAYLOC2*(X1114 IN(0
    5))+PAYLOC3*(X1125 IN(0 5))+PAYLOCO;

*   compute ratio of monthly payments to monthly income;
*   NOTE: for ADJINC=YES, income already inflated to same dollars as
    payments;
    %LET CPIADJ89=(&CPIADJ*(1886/&CPIBASE));
    PIRTOTAL=(TPAY/MAX((INCOME/12),(100/&CPIADJ89)));
    PIRMORT=(MORTPAY/MAX((INCOME/12),(100/&CPIADJ89)));
    PIRCONS=(CONSPAY/MAX((INCOME/12),(100/&CPIADJ89)));
    PIRREV=(REVPAY/MAX((INCOME/12),(100/&CPIADJ89)));
    PIR40=(PIRTOTAL > .4);
***************************************************************************;
*   loan purpose;
*   for convenience, initialize nonexistent debt questions at zero;
    %IF (&YEAR LT 1995) %THEN %DO;
      X7169=0;
      X7179=0;
      X7183=0;
    %END;
    %IF (&YEAR EQ 1989) %THEN %DO;
      X7824=0;
      X7847=0;
      X7870=0;
      X7924=0;
      X7947=0;
      X7970=0;
    %END;
*   assign a purpose code for borrowing where the purposes of the debt
    was not asked directly as a part of the individual loan sequence;

*   credit card debt assigned to goods and services (11);
    PURPCC=(CCBAL>0)*11;
*   first mortgages assigned to home purchase (1);
*   NOTE: in 1995+ surveys collected this info if date of mortgage was
    not the same as the date of home purchase, but code 1 assumed here
    for comparability with earlier years;
    PURPMORT=(X805>0)*1;
*   other home purchase loans assigned to home purchase (1);
    PURPOP=(X1044>0)*1;
*   mop-up line of credit assigned to unclassifiable (-7);
    PURPLOC=(X1136>0)*-7;
*   home improvement loans assigned to home improvement category 3,
    but code as 300 for programming convenience;
    PURPHI1=(X1215>0)*300;
    PURPHI2=(X1219>0)*300;
*   money still owed on properties for which R holds a land 
    contract assigned to debt for residential property (67);
* JXB - changed in 2013;
    %IF (&YEAR GE 2013) %THEN %DO;
      PURPLC1=(X1318>0)*67;
      PURPLC2=(X1337>0)*67;
      PURPLC4=(X1342>0)*67;
    %END;
    %ELSE %IF (&YEAR EQ 2010) %THEN %DO;
      PURPLC1=(X1417>0)*67;
      PURPLC2=(X1517>0)*67;
      PURPLC4=(X1621>0)*67;
    %END;
    %ELSE %IF (&YEAR LE 2007) %THEN %DO;
      PURPLC1=(X1417>0)*67;
      PURPLC2=(X1517>0)*67;
      PURPLC3=(X1617>0)*67;
      PURPLC4=(X1621>0)*67;
    %END;
*   residential property other than principal residence assigned to
    debt for residential property (67);
    PURPORE1=(X1703 IN (12 14 21 22 25 40 41 42 43 44 49 50 52 999))*
      (X1715>0)*67;
    PURPORE2=(X1803 IN (12 14 21 22 25 40 41 42 43 44 49 50 52 999))*
      (X1815>0)*67;
    %IF (&YEAR LT 2010) %THEN %DO;	 
      PURPORE3=(X1903 IN (12 14 21 22 25 40 41 42 43 44 49 50 52 999))*
        (X1915>0)*67;
    %END;
    PURPORE4=(X2006>0)*67;
*   vehicles assigned to vehicle category 10 but code as 100 for
    programming convenience;;
    PURPVEH1=(X2218>0)*100;
    PURPVEH2=(X2318>0)*100;
    PURPVEH3=(X2418>0)*100;
    PURPVEH4=.;
    %IF (&YEAR GE 1995) %THEN %DO;
      PURPVEH4=(X7169>0)*10;
    %END;
    PURPVEHM=(X2424>0)*100;
    PURPVEO1=(X2519>0)*100;
    PURPVEO2=(X2619>0)*100;
    PURPVEOM=(X2625>0)*100;
*   loans collected directly as student loans assigned to education
    category (83);
    %IF (&YEAR GE 1992) %THEN %DO;
      PURPED1=(X7824>0)*83;
      PURPED2=(X7847>0)*83;
      PURPED3=(X7870>0)*83;
      PURPED4=(X7924>0)*83;
      PURPED5=(X7947>0)*83;
      PURPED6=(X7970>0)*83;
      %IF (&YEAR GE 1995) %THEN %DO;
        PURPED7=(X7179>0)*83;
      %END;
      %ELSE %DO;
        PURPED7=0;
      %END;
    %END;
    %ELSE %DO;
      PURPED1=0;
      PURPED2=0;
      PURPED3=0;
      PURPED4=0;
      PURPED5=0;
      PURPED6=0;
      PURPED7=0;
    %END;
*   purpose of installment loans (excl loans for nonresidential property,
    which are netted out of such assets); 
    IF (X2710^=78|FLAG781=0) THEN PURPILN1=(X2710^=78|FLAG781=0)*X2710;
    ELSE PURPILN1=-9;
    IF (X2727^=78|FLAG781=0) THEN PURPILN2=(X2727^=78|FLAG781=0)*X2727;
    ELSE PURPILN2=-9;
    IF (X2810^=78|FLAG781=0) THEN PURPILN3=(X2810^=78|FLAG781=0)*X2810;
    ELSE PURPILN3=-9;
    IF (X2827^=78|FLAG781=0) THEN PURPILN4=(X2827^=78|FLAG781=0)*X2827;
    ELSE PURPILN4=-9;
    IF (X2910^=78|FLAG781=0) THEN PURPILN5=(X2910^=78|FLAG781=0)*X2910;
    ELSE PURPILN5=-9;
    IF (X2927^=78|FLAG781=0) THEN PURPILN6=(X2927^=78|FLAG781=0)*X2927;
    ELSE PURPILN6=-9;
*   other installment loans assigned to unclassifiable (-7);
    %IF (&YEAR GE 1995) %THEN %DO;
      PURPILN7=(X7183>0)*-7;
    %END;
    %ELSE %DO;
      PURPILN7=-9;
    %END;
  * margin loans assigned to investment (71);
    PURPMARG=(X3932>0)*71;
  * insurance loans assigned to unclassifiable (-7);
    PURPINS=(X4010>0)*-7;
  * misc sect N loans assigned to unclassifiable (-7);
    PURPOTH=(X4031>0)*-7;
  * pension loans assigned to unclassifiable (-99) for years before 1998;
    %IF (&YEAR GE 2010) %THEN %DO;
      IF (X11070 IN (5 0)) THEN PURPPEN1=X11030;
      ELSE PURPPEN1=-9;
      IF (X11170 IN (5 0)) THEN PURPPEN2=X11130;
      ELSE PURPPEN2=-9;
      PURPPEN3=-9;
      IF (X11370 IN (5 0)) THEN PURPPEN4=X11330;
      ELSE PURPPEN4=-9;
      IF (X11470 IN (5 0)) THEN PURPPEN5=X11430;
      ELSE PURPPEN5=-9;
      PURPPEN6=-9;
    %END;   
    %ELSE %IF (&YEAR GE 2004) %THEN %DO;
      IF (X11070 IN (5 0)) THEN PURPPEN1=X11030;
      ELSE PURPPEN1=-9;
      IF (X11170 IN (5 0)) THEN PURPPEN2=X11130;
      ELSE PURPPEN2=-9;
      IF (X11270 IN (5 0)) THEN PURPPEN3=X11230;
      ELSE PURPPEN3=-9;
      IF (X11370 IN (5 0)) THEN PURPPEN4=X11330;
      ELSE PURPPEN4=-9;
      IF (X11470 IN (5 0)) THEN PURPPEN5=X11430;
      ELSE PURPPEN5=-9;
      IF (X11570 IN (5 0)) THEN PURPPEN6=X11530;
      ELSE PURPPEN6=-9;
    %END;
    %ELSE %IF (&YEAR GE 1998) %THEN %DO;
      IF (X4230 IN (5 0)) THEN PURPPEN1=X6791;
      ELSE PURPPEN1=-9;
      IF (X4330 IN (5 0)) THEN PURPPEN2=X6792;
      ELSE PURPPEN2=-9;
      IF (X4430 IN (5 0)) THEN PURPPEN3=X6793;
      ELSE PURPPEN3=-9;
      IF (X4830 IN (5 0)) THEN PURPPEN4=X6794;
      ELSE PURPPEN4=-9;
      IF (X4930 IN (5 0)) THEN PURPPEN5=X6795;
      ELSE PURPPEN5=-9;
      IF (X5030 IN (5 0)) THEN PURPPEN6=X6796;
      ELSE PURPPEN6=-9;
    %END;
    %ELSE %DO;
      IF (X4230 IN (5 0)) THEN PURPPEN1=(X4229>0)*(-99)*(X4230=5);
      ELSE PURPPEN1=-9;
      IF (X4330 IN (5 0)) THEN PURPPEN2=(X4329>0)*(-99)*(X4330=5);
      ELSE PURPPEN2=-9;
      IF (X4430 IN (5 0)) THEN PURPPEN3=(X4429>0)*(-99)*(X4430=5);
      ELSE PURPPEN3=-9;
      IF (X4830 IN (5 0)) THEN PURPPEN4=(X4829>0)*(-99)*(X4830=5);
      ELSE PURPPEN4=-9;
      IF (X4930 IN (5 0)) THEN PURPPEN5=(X4929>0)*(-99)*(X4930=5);
      ELSE PURPPEN5=-9;
      IF (X5030 IN (5 0)) THEN PURPPEN6=(X5029>0)*(-99)*(X5030=5);
      ELSE PURPPEN6=-9;
    %END;

***************************************************************************;
*   assign codes for type of lender where type is not directly known;

*   credit card debt assigned to store and credit cards;
    TYPECC=98;
*   mopup LOC assigned to unclassified ;
    TYPELOC=97;
*   home improvement loans assigned to unclassified;
    TYPEHI1=97;
    TYPEHI2=97;
*   outstanding loans on land contracts held by HH assigned to unclassified;
    TYPELC1=97;
    TYPELC2=97;
    TYPELC3=97;
    TYPELC4=97;
*   loans for vacation homes and recreational properties assigned to
    unclassified;
    TYPEORE4=97;
    TYPEORE5=97;
*   loans for remaining passenger vehicles assigned to unclassified;
    TYPEVEHM=97;
*   loans for remaining other vehicles assigned to unclassified;
    TYPEVEOM=97;
*   other education loans assigned to unclassified; 
    TYPEEDU7=97;
*   adjust institution for assignment of loans to NNRESRE;
    INSTALL_I1=X9107*(X2710^=78|FLAG781=0);
    INSTALL_I2=X9108*(X2727^=78|FLAG781=0);
    INSTALL_I3=X9109*(X2810^=78|FLAG781=0);
    INSTALL_I4=X9110*(X2827^=78|FLAG781=0);
    INSTALL_I5=X9111*(X2910^=78|FLAG781=0);
    INSTALL_I6=X9112*(X2927^=78|FLAG781=0);
*   other installment loans assigned to unclassified; 
    TYPEILN7=97;
*   margin loans assigned to brokerages;
    TYPEMARG=16;
*   loans against insurance policies assigned to insurance companies;
    TYPEINS=17;
*   other loans assigned to individual lender NEC;
    TYPEOTH=27;
*   loans against pensions assigned to borrowing against pension plans;
    TYPEPEN1=96;
    TYPEPEN2=96;
    TYPEPEN3=96;
    TYPEPEN4=96;
    TYPEPEN5=96;
    TYPEPEN6=96;

***************************************************************************;
* summarize loan balances by purpose and lender;

*   NOTE: arrays of purposes, lenders, and amounts outstanding must be
*   in same order;
*   array purposes for loans;

  %IF (&YEAR GE 2013) %THEN %DO;	 
    ARRAY PURPOSES {*} PURPCC PURPMORT X918  X1018 PURPOP
      X1106 X1117 X1128 PURPLOC PURPHI1 PURPHI2  
      PURPLC1 PURPLC2 PURPLC4
      PURPORE1 PURPORE2 PURPORE4
      PURPVEH1 PURPVEH2 PURPVEH3 PURPVEH4 PURPVEHM
      PURPVEO1 PURPVEO2 PURPVEOM
      PURPED1 PURPED2 PURPED3 PURPED4 PURPED5 PURPED6 PURPED7
      PURPILN1 PURPILN2 PURPILN3 PURPILN4 PURPILN5 PURPILN6 PURPILN7
      PURPMARG PURPINS PURPOTH       
      PURPPEN1 PURPPEN2 PURPPEN4  PURPPEN5;
*   array amounts outstanding;
    ARRAY DEBTOUT {*}  CCBAL X805 X905  X1005 X1044
      X1108 X1119 X1130 X1136 X1215 X1219 
      X1318 X1337 X1342
      MORT1   MORT2   X2006
      X2218    X2318    X2418    X7169    X2424
      X2519    X2619    X2625
      X7824   X7847   X7870   X7924   X7947   X7970   X7179
      X2723 X2740 X2823 X2840 X2923 X2940 X7183
      OUTMARG  X4010   X4032
      OUTPEN1  OUTPEN2  OUTPEN4   OUTPEN5; 
*   array types of lenders;
    ARRAY LENTYPE {*}  TYPECC X9083    X9084 X9085 X9086 
      X9087 X9088 X9089 TYPELOC X9090   TYPEHI2  
      TYPELC1 TYPELC2 TYPELC4
      X9099    X9100    TYPEORE4
      X9102    X9103    X9104    X9215    TYPEVEHM
      X9105    X9106    TYPEVEOM
      X9203   X9204   X9205   X9206   X9207   X9208   TYPEEDU7
      INSTALL_I1 INSTALL_I2 INSTALL_I3 INSTALL_I4 INSTALL_I5
      INSTALL_I6 TYPEILN7 TYPEMARG TYPEINS TYPEOTH        
      TYPEPEN1 TYPEPEN2 TYPEPEN4  TYPEPEN5;
  %END;
  %ELSE %IF (&YEAR EQ 2010) %THEN %DO;	 
    ARRAY PURPOSES {*} PURPCC PURPMORT X918  X1018 PURPOP
      X1106 X1117 X1128 PURPLOC PURPHI1 PURPHI2  
      PURPLC1 PURPLC2 PURPLC4
      PURPORE1 PURPORE2 PURPORE4
      PURPVEH1 PURPVEH2 PURPVEH3 PURPVEH4 PURPVEHM
      PURPVEO1 PURPVEO2 PURPVEOM
      PURPED1 PURPED2 PURPED3 PURPED4 PURPED5 PURPED6 PURPED7
      PURPILN1 PURPILN2 PURPILN3 PURPILN4 PURPILN5 PURPILN6 PURPILN7
      PURPMARG PURPINS PURPOTH       
      PURPPEN1 PURPPEN2 PURPPEN4  PURPPEN5;
*   array amounts outstanding;
    ARRAY DEBTOUT {*}  CCBAL X805 X905  X1005 X1044
      X1108 X1119 X1130 X1136 X1215 X1219 
      X1417 X1517 X1621
      MORT1   MORT2   X2006
      X2218    X2318    X2418    X7169    X2424
      X2519    X2619    X2625
      X7824   X7847   X7870   X7924   X7947   X7970   X7179
      X2723 X2740 X2823 X2840 X2923 X2940 X7183
      OUTMARG  X4010   X4032
      OUTPEN1  OUTPEN2  OUTPEN4   OUTPEN5; 
*   array types of lenders;
    ARRAY LENTYPE {*}  TYPECC X9083    X9084 X9085 X9086 
      X9087 X9088 X9089 TYPELOC X9090   TYPEHI2  
      TYPELC1 TYPELC2 TYPELC4
      X9099    X9100    TYPEORE4
      X9102    X9103    X9104    X9215    TYPEVEHM
      X9105    X9106    TYPEVEOM
      X9203   X9204   X9205   X9206   X9207   X9208   TYPEEDU7
      INSTALL_I1 INSTALL_I2 INSTALL_I3 INSTALL_I4 INSTALL_I5
      INSTALL_I6 TYPEILN7 TYPEMARG TYPEINS TYPEOTH        
      TYPEPEN1 TYPEPEN2 TYPEPEN4  TYPEPEN5;
  %END;
  %ELSE %IF (&YEAR LT 2010) %THEN %DO;
    ARRAY PURPOSES {*} PURPCC PURPMORT X918  X1018 PURPOP
      X1106 X1117 X1128 PURPLOC PURPHI1 PURPHI2  
      PURPLC1 PURPLC2 PURPLC3 PURPLC4
      PURPORE1 PURPORE2 PURPORE3 PURPORE4
      PURPVEH1 PURPVEH2 PURPVEH3 PURPVEH4 PURPVEHM
      PURPVEO1 PURPVEO2 PURPVEOM
      PURPED1 PURPED2 PURPED3 PURPED4 PURPED5 PURPED6 PURPED7
      PURPILN1 PURPILN2 PURPILN3 PURPILN4 PURPILN5 PURPILN6 PURPILN7
      PURPMARG PURPINS PURPOTH       
      PURPPEN1 PURPPEN2 PURPPEN3 PURPPEN4  PURPPEN5 PURPPEN6;
*   array amounts outstanding;
    ARRAY DEBTOUT {*}  CCBAL X805 X905  X1005 X1044
      X1108 X1119 X1130 X1136 X1215 X1219 
      X1417 X1517 X1617 X1621
      MORT1   MORT2   MORT3   X2006
      X2218    X2318    X2418    X7169    X2424
      X2519    X2619    X2625
      X7824   X7847   X7870   X7924   X7947   X7970   X7179
      X2723 X2740 X2823 X2840 X2923 X2940 X7183
      OUTMARG  X4010   X4032
      OUTPEN1  OUTPEN2  OUTPEN3   OUTPEN4   OUTPEN5  OUTPEN6; 
*   array types of lenders;
    ARRAY LENTYPE {*}  TYPECC X9083    X9084 X9085 X9086 
      X9087 X9088 X9089 TYPELOC X9090   TYPEHI2  
      TYPELC1 TYPELC2 TYPELC3 TYPELC4
      X9099    X9100    X9101    TYPEORE4
      X9102    X9103    X9104    X9215    TYPEVEHM
      X9105    X9106    TYPEVEOM
      X9203   X9204   X9205   X9206   X9207   X9208   TYPEEDU7
      INSTALL_I1 INSTALL_I2 INSTALL_I3 INSTALL_I4 INSTALL_I5
      INSTALL_I6 TYPEILN7 TYPEMARG TYPEINS TYPEOTH        
      TYPEPEN1 TYPEPEN2 TYPEPEN3  TYPEPEN4  TYPEPEN5 TYPEPEN6;
  %END;	 

*   aggregate balances by loan purposes:
*   NOTE: beginning with 1998, loan purpose codes were collapsed in
    the public dataset in a way that makes it impossible to recreate
    exactly the same loan purpose categories as those used in the
    1/2000 Bulletin article: for consistency across years, the 1998
    collapsed categories are used in creating loan purpose categories
    here for the public datasets;
*   NOTE: Starting in 2007, "Other unclassifiable loans" includes
    unclassifiable borrowing against pension accounts, which was 
    shown separately in prior years;
*   categories for the internal data:
    1=home purchase loans
    2=home improvement loans
    3=vehicle loans
    4=loan for purchase of goods and services+prof services
    5=loans for investments
    6=loans for education
    7=mortgage loans for other real estate
    8=other unclassifiable loans;
*   categories for the public data:
    1=loans for home purchase, cottage, vacation property
    2=home improvement loans
    3=vehicle loans
    4=loan for purchase of goods and services
    5=loans for investments and mortgage loans for other real estate
    6=loans for education and loans for professional expenses
    7=other unclassifiable loans;
*   NOTE: purpose variable for intallment loans set to zero where loan
    has been allocated to NNRESRE;

    ARRAY PLOAN{*} PLOAN1-PLOAN8;
*   initialize;
    DO I=1 TO DIM(PLOAN);
      PLOAN{I}=0;
    END;
    DO I=1 TO DIM(PURPOSES);
      %IF (&PUBLIC NE YES) %THEN %DO;
        CALL VNAME(DEBTOUT{I},VARNM);
        IF (PURPOSES{I} IN (1)) THEN PLOAN1=PLOAN1+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I} IN (3 300 4)) THEN
          PLOAN2=PLOAN2+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I} IN (10 24 61 63 65 100)) THEN PLOAN3=
          PLOAN3+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I} IN ( -1 6 11 12 13 14 15 16 17 18 20 23 25 
          26 27 29 31 32 34 35 36 40 49 50 69 80 81 82 84 85
          86 88 89 90 91 92 93 94 95 96 97)) THEN
	  PLOAN4=PLOAN4+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I} IN (71 72 73 74 75 76 79)) THEN PLOAN5=
          PLOAN5+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I} IN (83)) THEN PLOAN6=PLOAN6+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I}=67) THEN PLOAN7=PLOAN7+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I}=78) THEN PLOAN7=PLOAN7+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I} IN (-99 -7)) THEN PLOAN8=PLOAN8+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I} NOT IN (0 . -9) | (DEBTOUT{I}>0 &
          PURPOSES{I}^=-9)) THEN  PUT "ERROR? UNCLASSIFIED LOAN TYPE "
	  &ID= DEBTOUT{I}= PURPOSES{I} =;
      %END;
      %ELSE %IF (&YEAR GE 1998) %THEN %DO;
        IF (PURPOSES{I} IN (1 67)) THEN PLOAN1=PLOAN1+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I} IN (2 300)) THEN PLOAN2=PLOAN2+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I} IN (3 6 100)) THEN
          PLOAN3=PLOAN3+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I} IN (4 5 8 10 11)) THEN
          PLOAN4=PLOAN4+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I} IN (7 71)) THEN PLOAN5=PLOAN5+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I} IN (9 83)) THEN PLOAN6=PLOAN6+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I} IN (-99 -7)) THEN PLOAN7=PLOAN7+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I} NOT IN (0 . -9) | (DEBTOUT{I}>0 &
          PURPOSES{I}^=-9)) THEN  PUT "ERROR? UNCLASSIFIED LOAN TYPE "
	  &ID= DEBTOUT{I}= PURPOSES{I} =;
      %END;
      %ELSE %DO;
        IF (PURPOSES{I} IN (1 67)) THEN PLOAN1=PLOAN1+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I} IN (3 300 4)) THEN
          PLOAN2=PLOAN2+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I} IN (10 100 24 61 63 65)) THEN
          PLOAN3=PLOAN3+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I} IN (6 11 12 13 14 15 16 17 18 20 23 25 
          26 29 31 34 35 36 49 50 69 80 81 84 85
          88 89 90 91 92 93 94 95 96)) THEN
	  PLOAN4=PLOAN4+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I} IN (71 72 73 74 75 76 79)) THEN PLOAN5=
          PLOAN5+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I}=78) THEN PLOAN5=PLOAN5+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I} IN (82 83)) THEN PLOAN6=PLOAN6+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I} IN (-99 -7)) THEN PLOAN7=PLOAN7+MAX(0,DEBTOUT{I});
        ELSE IF (PURPOSES{I} NOT IN (0 . -9) | (DEBTOUT{I}>0 &
          PURPOSES{I}^=-9)) THEN  PUT "ERROR? UNCLASSIFIED LOAN TYPE "
	  &ID= DEBTOUT{I}= PURPOSES{I} =;
      %END;
    END;
   
*   calculate alternative aggregate balances by loan purpose that, for
    2004 forward, splits out purpose for first-lien equity extraction;
    ARRAY PLOANB{*} PLOANB1-PLOANB8;
    DO I=1 TO DIM(PLOAN);
      PLOANB{I}=PLOAN{I};
    END;
    %IF (&YEAR GE 2004 & &PUBLIC NE YES) %THEN %DO;
      IF X7137 IN (2 3 4) THEN DO;
        IF X7137 IN (2 3) THEN TOTEXTRACT=X7138*(X805/MAX(1,X804));
        ELSE IF X7137=4 THEN TOTEXTRACT=X805;
        PLOANB1=PLOAN1-TOTEXTRACT;

        IF (X6723 IN (1)) THEN PLOANB1=PLOANB1+MAX(0,TOTEXTRACT);
        ELSE IF (X6723 IN (3 300 4)) THEN
          PLOANB2=PLOANB2+MAX(0,TOTEXTRACT);
        ELSE IF (X6723 IN (10 24 61 63 65 100)) THEN PLOANB3=
          PLOANB3+MAX(0,TOTEXTRACT);
        ELSE IF (X6723 IN ( -1 6 11 12 13 14 15 16 17 18 20 23 25 
          26 27 29 31 32 34 35 36 40 49 50 69 80 81 82 84 85
          86 88 89 90 91 92 93 94 95 96 97)) THEN
        PLOANB4=PLOANB4+MAX(0,TOTEXTRACT);
        ELSE IF (X6723 IN (71 72 73 74 75 76 79)) THEN PLOANB5=
          PLOANB5+MAX(0,TOTEXTRACT);
        ELSE IF (X6723 IN (83)) THEN PLOANB6=PLOANB6+MAX(0,TOTEXTRACT);
        ELSE IF (X6723=67) THEN PLOANB7=PLOANB7+MAX(0,TOTEXTRACT);
        ELSE IF (X6723=78) THEN PLOANB7=PLOANB7+MAX(0,TOTEXTRACT);
        ELSE IF (X6723 IN (-99 -7)) THEN PLOANB8=PLOANB8+MAX(0,TOTEXTRACT);
        ELSE IF (X6723 NOT IN (0 . -9) | (TOTEXTRACT>0 &
          X6723^=-9)) THEN  PUT "ERROR? UNCLASSIFIED LOAN TYPE "
          &ID= TOTEXTRACT= X6723=;
      END;   
    %END;

*   as check compute total by loan purpose (TPLOAN, TPLOANB) and 
    compare with total debt;
    TPLOAN=0;
    TPLOANB=0;
    DO I=1 TO DIM(PLOAN);
      TPLOAN=TPLOAN+MAX(0,PLOAN{I});
      TPLOANB=TPLOANB+MAX(0,PLOANB{I});
    END;
    IF (ROUND(DEBT,1)^=ROUND(TPLOAN,1)) THEN PUT
      "ERROR IN TOTAL DEBT/DEBT BY PURPOSE " &ID= DEBT= TPLOAN=;
    IF (ROUND(DEBT,1)^=ROUND(TPLOANB,1)) THEN PUT
      "ERROR IN TOTAL DEBT/DEBT BY PURPOSE " &ID= DEBT= TPLOANB=;

*   aggregate total balances by institution type:
    1=commercial bank 
    2=savings and loan 
    3=credit union 
    4=finance, loan or leasing company, inc debt consolidator  
    5=brokerage, broad financial services company and life insurance
    6=real estate lender
    7=individual lender 
    8=other nonfinancial
    9=government
    10=store and credit cards
    11=Loans against pensions
    12=Other unclassifiable (inc foreign) ;
    ARRAY LLOAN{*} LLOAN1-LLOAN12;
*   initialize;
    DO I=1 TO DIM(LLOAN);
      LLOAN{I}=0;
    END;
    DO I=1 TO DIM(LENTYPE);
      CALL VNAME(DEBTOUT{I},VARNM);
        IF (LENTYPE{I} IN (11)) THEN LLOAN1=LLOAN1+MAX(0,DEBTOUT{I});
        ELSE IF (LENTYPE{I} IN (12)) THEN LLOAN2=LLOAN2+MAX(0,DEBTOUT{I});
        ELSE IF (LENTYPE{I} IN (13)) THEN LLOAN3=LLOAN3+MAX(0,DEBTOUT{I});
        ELSE IF (LENTYPE{I} IN (14 21 43 56)) THEN LLOAN4=LLOAN4+
          MAX(0,DEBTOUT{I});
        ELSE IF (LENTYPE{I} IN (16 17 92 94 29)) THEN LLOAN5=LLOAN5+
          MAX(0,DEBTOUT{I});
        ELSE IF (LENTYPE{I} IN (18 19 31)) THEN LLOAN6=LLOAN6+
          MAX(0,DEBTOUT{I});
        ELSE IF (LENTYPE{I} IN (20 26 27 85)) THEN LLOAN7=LLOAN7+
          MAX(0,DEBTOUT{I});
        ELSE IF (LENTYPE{I} IN (15 22 23 24 25 30 32 39 40 41 44 45 46 47
          57 61 80 81 95 99)) THEN LLOAN8=LLOAN8+MAX(0,DEBTOUT{I});
        ELSE IF (LENTYPE{I} IN (33 34 35 42 93)) THEN LLOAN9=
          LLOAN9+MAX(0,DEBTOUT{I});
        ELSE IF (LENTYPE{I} IN (38 50 51 52 90 98)) THEN LLOAN10=
          LLOAN10+MAX(0,DEBTOUT{I});
        ELSE IF (LENTYPE{I} IN (96)) THEN LLOAN11=LLOAN11+MAX(0,DEBTOUT{I});
        ELSE IF (LENTYPE{I} IN (-7 -1 28 37 75 97)) THEN LLOAN12=
          LLOAN12+MAX(0,DEBTOUT{I});
        ELSE IF (DEBTOUT{I}>0 & LENTYPE{I}>.Z & LENTYPE{I}^=0) THEN PUT
          'ERROR? UNCLASSIFIED LENDER TYPE ' &ID= LENTYPE{I}= DEBTOUT{I}= ;
        IF DEBTOUT{I}<=.Z THEN PUT 'BALANCE MISSING ' &ID= DEBTOUT{I}= ;
    END;
*   as check compute total by loan purpose (TLLOAN) and compare with
    total debt;
    TLLOAN=0;
    DO I=1 TO DIM(LLOAN);
      TLLOAN=TLLOAN+MAX(0,LLOAN{I});
    END;
    IF (ROUND(DEBT,1)^=ROUND(TLLOAN,1)) THEN PUT
      "ERROR IN TOTAL DEBT/DEBT BY INSTITUTION TYPE " &ID= DEBT= TLLOAN=;

***************************************************************************;
*   miscellaneous data fixes;

    %IF (&PUBLIC NE YES) %THEN %DO;
      %INCLUDE FIX&YEAR;
    %END;


***************************************************************************;
*   where REAL=YES, put key dollar terms into dollars corresponding
    the BASE year;
%*  adjust all dollar values;
    %IF (&REAL EQ YES) %THEN %DO;
      ARRAY VALUE{*} INCOME ASSET NETWORTH FIN LIQ CDS NMMF STOCKS
        BOND RETQLIQ SAVBND CASHLI OTHMA OTHFIN NFIN VEHIC HOUSES ORESRE
        NNRESRE BUS OTHNFIN DEBT MRTHEL RESDBT OTHLOC CCBAL INSTALL
        ODEBT KGTOTAL KGHOUSE KGORE KGBUS KGSTMF TPAY MORTPAY CONSPAY
        REVPAY TPLOAN PLOAN1-PLOAN8 TLLOAN LLOAN1-LLOAN12 EQUITY DEQ
        VLEASE VOWN RETEQ NORMINC CHECKING SAVING MMA CALL HOMEEQ
        IRAKH PENEQ VEH_INST EDN_INST OTH_INST HELOC NH_MORT WAGEINC 
        BUSSEFARMINC INTDIVINC KGINC SSRETINC TRANSFOTHINC RENTINC
        NHNFIN THRIFT CURRPEN FUTPEN STMUTF TFBMUTF
        GBMUTF OBMUTF COMUTF OMUTF ANNUIT TRUSTS ACTBUS NONACTBUS NOTXBND
        MORTBND GOVTBND OBND OUTPEN OUTMARG
        PAYMORT1 PAYMORT2 PAYMORT3 PAYMORTO PAYLOC1
        PAYLOC2 PAYLOC3 PAYLOCO
        PAYHI1 PAYHI2 PAYLC1 PAYLC2 PAYLCO PAYORE1 PAYORE2 PAYORE3
        PAYOREV PAYVEH1 PAYVEH2 PAYVEH3 PAYVEH4 PAYVEHM PAYVEO1 PAYVEO2
        PAYVEOM PAYEDU1 PAYEDU2 PAYEDU3 PAYEDU4 PAYEDU5 PAYEDU6 PAYEDU7
        PAYILN1 PAYILN2 PAYILN3 PAYILN4 PAYILN5 PAYILN6 PAYILN7 PAYMARG
        PAYINS PAYPEN1 PAYPEN2 PAYPEN3 PAYPEN4 PAYPEN5 PAYPEN6 FARMBUS_KG MORT1
        MORT2 MORT3 FARMBUS NHNFIN PENACCTWD MMDA MMMF FOODHOME FOODAWAY FOODDELV; 
      DO I=1 TO DIM(VALUE);
        IF VALUE{I}<=.Z & VALUE{I}^=.B THEN PUT "MISSING VALUE " &ID=
          VALUE{I}=;
        IF VALUE{I} NOT IN(-1 -2) THEN VALUE{I}=VALUE{I}*&CPIADJ;
      END;
    %END;

***************************************************************************;
*   create categorical variables based on dollar variables;
*   NOTE: need to put here so categories are fixed in real terms when
    REAL=YES;
    IF INCOME<10000 THEN INCCL2=1;
    ELSE IF 10000<=INCOME<25000 THEN INCCL2=2;
    ELSE IF 25000<=INCOME<50000 THEN INCCL2=3;
    ELSE IF 50000<=INCOME<100000 THEN INCCL2=4;
    ELSE IF INCOME>=100000 THEN INCCL2=5;

***************************************************************************;
  RUN;
***************************************************************************;
***************************************************************************;
* compute net worth and income quantiles and merge with main dataset;

  %MACRO PCTL(VAR=,PPOINTS=,TAG=);
    PROC MEANS NOPRINT DATA=SCFB&YEAR;
      VAR WGT;
      OUTPUT OUT=TOTPOP SUM=TOTPOP;
    RUN;
    DATA TOTPOP;
      SET TOTPOP;
      MERGEID=1;
    RUN;
    DATA %UNQUOTE(&TAG.CAT)(KEEP=WGT &VAR &PID &ID TOTPOP);
      MERGE SCFB&YEAR(KEEP=WGT &VAR &PID MERGEID &ID) TOTPOP;
      BY MERGEID;
    RUN;
    PROC SORT DATA=&TAG.CAT;
      BY &VAR &PID;
    RUN;
    DATA %UNQUOTE(&TAG.CAT)(KEEP=&ID &TAG.CAT);
      SET %UNQUOTE(&TAG.CAT);
      RETAIN CUMSUM 0;

*     extract number of items in percentile point list;
      %LET I=1;
      %LET NEXT=%SCAN(&PPOINTS,&I);
      %DO %UNTIL (&NEXT EQ );
        %LET I=%EVAL(&I + 1);
        %LET NEXT=%SCAN(&PPOINTS,&I);
      %END;
      %LET NCAT=%EVAL (&I-1);
        
      CUMSUM=CUMSUM+WGT;

      %LET NEXT=%SCAN(&PPOINTS,&NCAT);
      IF ((CUMSUM/TOTPOP)>=%UNQUOTE(%STR(.)&NEXT)) THEN &TAG.CAT=&NCAT;
      %LET I=%EVAL(&NCAT-1);
      %LET NEXT=%SCAN(&PPOINTS,&I));
      %DO I=1 %TO %EVAL(&NCAT-1);
        %LET WNCAT=%EVAL(&NCAT-&I);
        %LET NEXT=%SCAN(&PPOINTS,&WNCAT);
        ELSE IF ((CUMSUM/TOTPOP)>=%UNQUOTE(%STR(.)&NEXT)) THEN
          &TAG.CAT=&WNCAT;
      %END;
*     print out the empirical breakpoints;
      CATVAL=&TAG.CAT;
      IF (LAG1(CATVAL)^=CATVAL) THEN PUT &TAG.CAT= &VAR=;
    RUN;
    PROC SORT DATA=%UNQUOTE(&TAG.CAT);
      BY &ID;
    RUN;
  %MEND PCTL;
  %PCTL(VAR=NETWORTH,PPOINTS=0 25 50 75 90,TAG=NW);
  %PCTL(VAR=INCOME,PPOINTS=0 20 40 60 80 90,TAG=INC);
  %PCTL(VAR=ASSET,PPOINTS=0 20 40 60 80 90,TAG=ASSET);
  %PCTL(VAR=NORMINC,PPOINTS=0 20 40 60 80 90,TAG=NINC);
  %PCTL(VAR=NORMINC,PPOINTS=0 50 90,TAG=NINC2);
  %PCTL(VAR=NETWORTH,PPOINTS=0 10 20 30 40 50 60 70 80 90,TAG=NW10);
  %PCTL(VAR=INCOME,PPOINTS=0 10 20 30 40 50 60 70 80 90,TAG=INC10);
  %PCTL(VAR=NORMINC,PPOINTS=0 10 20 30 40 50 60 70 80 90,TAG=NINC10);
  DATA SCFB&YEAR;
    MERGE SCFB&YEAR NWCAT INCCAT ASSETCAT NINCCAT NINC2CAT NW10CAT
      INC10CAT NINC10CAT;
    BY &ID;
  RUN;

***************************************************************************;
***************************************************************************;
* check for overly influential observations;

* for convenience, compute mean of weight/key dollar values over implicates;
  PROC MEANS NOPRINT;
    WEIGHT WGT;
    VAR WGT FIN LIQ CDS NMMF STOCKS BOND RETQLIQ SAVBND
      CASHLI OTHMA OTHFIN NFIN VEHIC HOUSES ORESRE NNRESRE
      BUS OTHNFIN DEBT MRTHEL RESDBT OTHLOC CCBAL INSTALL
      ODEBT ASSET NETWORTH TPAY;
    CLASS &ID; 
    OUTPUT OUT=OUT0 MEAN=
      MWGT MFIN MLIQ MCDS MNMMF MSTOCKS MBOND MRETQLIQ MSAVBND
      MCASHLI MOTHMA MOTHFIN MNFIN MVEHIC MHOUSES MORESRE MNNRESRE
      MBUS MOTHNFIN MDEBT MMRTHEL MRESDBT MOTHLOC MCCBAL MINSTALL
      MODEBT MASSET MNETWORT MTPAY;
  RUN;
  DATA OUT0;
    SET OUT0;
    MERGEID=1; 
  RUN;
* compute weighted sum of mean over implicates for each variable;
  PROC MEANS DATA=OUT0 NOPRINT SUM;
    WEIGHT MWGT;
    VAR MFIN MLIQ MCDS MNMMF MSTOCKS MBOND MRETQLIQ MSAVBND
      MCASHLI MOTHMA MOTHFIN MNFIN MVEHIC MHOUSES MORESRE MNNRESRE
      MBUS MOTHNFIN MDEBT MMRTHEL MRESDBT MOTHLOC MCCBAL MINSTALL
      MODEBT MASSET MNETWORT MTPAY;
    OUTPUT OUT=OUT1 SUM=
      TFIN TLIQ TCDS TNMMF TSTOCKS TBOND TRETQLIQ TSAVBND
      TCASHLI TOTHMA TOTHFIN TNFIN TVEHIC THOUSES TORESRE TNNRESRE
      TBUS TOTHNFIN TDEBT TMRTHEL TRESDBT TOTHLOC TCCBAL TINSTALL
      TODEBT TASSET TNETWORT TTPAY;
  RUN;
  DATA OUT1;
    SET OUT1;
    MERGEID=1;
  RUN;
* compare share of average of implicates to totals and print out cases
  accounting for 1.5 percent or more of the totals;
  DATA OUT2;
     MERGE OUT0 OUT1;
     BY MERGEID ;

    ARRAY SVAR {*}
      SFIN SLIQ SCDS SNMMF SSTOCKS SBOND SRETQLIQ SSAVBND
      SCASHLI SOTHMA SOTHFIN SNFIN SVEHIC SHOUSES SORESRE SNNRESRE
      SBUS SOTHNFIN SDEBT SMRTHEL SRESDBT SOTHLOC SCCBAL SINSTALL
      SODEBT SASSET SNETWORT STPAY;
    ARRAY TVAR {*}
      TFIN TLIQ TCDS TNMMF TSTOCKS TBOND TRETQLIQ TSAVBND
      TCASHLI TOTHMA TOTHFIN TNFIN TVEHIC THOUSES TORESRE TNNRESRE
      TBUS TOTHNFIN TDEBT TMRTHEL TRESDBT TOTHLOC TCCBAL TINSTALL
      TODEBT TASSET TNETWORT TTPAY;
    ARRAY MVAR {*}
      MFIN MLIQ MCDS MNMMF MSTOCKS MBOND MRETQLIQ MSAVBND
      MCASHLI MOTHMA MOTHFIN MNFIN MVEHIC MHOUSES MORESRE MNNRESRE
      MBUS MOTHNFIN MDEBT MMRTHEL MRESDBT MOTHLOC MCCBAL MINSTALL
      MODEBT MASSET MNETWORT MTPAY;
    DO I=1 TO DIM(SVAR);
      SVAR{I}=MWGT*MVAR{I}/TVAR{I};
      SVAR{I}=SVAR{I}*100;
      SVAR{I}=ROUND(SVAR{I},.1);
      IF (SVAR{I}>1.5) THEN PUT &ID= MWGT= MVAR{I}= TVAR{I}= SVAR{I}=;
      IF (SVAR{I}<=.Z) THEN PUT "MISSING VALUE " &ID= MWGT= MVAR{I}=
        TVAR{I}= SVAR{I}=;
    END; 

  RUN;

***************************************************************************;
***************************************************************************;
* save final version of the dataset;
  %IF (&REAL EQ YES) %THEN %LET WORKING=WORKINGR;
  %ELSE %LET WORKING=WORKINGN;
  %IF (&PUBLIC EQ YES) %THEN %LET DSNAME=&WORKING%STR(.SCFP);
  %ELSE %LET DSNAME=&WORKING%STR(.SCFB);
  DATA &DSNAME&YEAR;
    SET SCFB&YEAR;
    DROP MERGEID;
  RUN;
  PROC CONTENTS;
  RUN;
* export data set for use in spreadsheet;
* export file for NORC;
  %IF (&REAL EQ YES) %THEN %LET NR=real;
  %ELSE %LET NR=nominal;
  %IF (&PUBLIC EQ YES) %THEN %LET PUB=SCFP;
  %ELSE %LET PUB=SCFB;
  PROC EXPORT DATA=%UNQUOTE(&DSNAME&YEAR)
    OUTFILE=
      "/mecs/scf7/analysis/bulletin/2013/data/transf/bull_&NR/&PUB&YEAR%STR(.)csv"
    DBMS=CSV REPLACE;
  RUN;


***************************************************************************;
***************************************************************************;

%MEND BULLIT;
***************************************************************************;
***************************************************************************;

* CPI-U-RS here should be September values;
* omit decimal in CPI;
* for the public datasets, only adjust past year income to the survey year; 

%BULLIT(YEAR=2013,REAL=NO,ADJINC=NO,CPIBASE=3438,PUBLIC=YES);
%BULLIT(YEAR=2010,REAL=NO,ADJINC=NO,CPIBASE=3208,PUBLIC=YES);
%BULLIT(YEAR=2007,REAL=NO,ADJINC=NO,CPIBASE=3062,PUBLIC=YES);
%BULLIT(YEAR=2004,REAL=NO,ADJINC=NO,CPIBASE=2788,PUBLIC=YES);
%BULLIT(YEAR=2001,REAL=NO,ADJINC=NO,CPIBASE=2618,PUBLIC=YES);
%BULLIT(YEAR=1998,REAL=NO,ADJINC=NO,CPIBASE=2405,PUBLIC=YES);
%BULLIT(YEAR=1995,REAL=NO,ADJINC=NO,CPIBASE=2265,PUBLIC=YES);
%BULLIT(YEAR=1992,REAL=NO,ADJINC=NO,CPIBASE=2116,PUBLIC=YES);
%BULLIT(YEAR=1989,REAL=NO,ADJINC=NO,CPIBASE=1902,PUBLIC=YES);

/*
* make real versions of public data;
%BULLIT(YEAR=2013,REAL=YES,ADJINC=YES,CPIBASE=3438,PUBLIC=YES);
%BULLIT(YEAR=2010,REAL=YES,ADJINC=YES,CPIBASE=3438,PUBLIC=YES);
%BULLIT(YEAR=2007,REAL=YES,ADJINC=YES,CPIBASE=3438,PUBLIC=YES);
%BULLIT(YEAR=2004,REAL=YES,ADJINC=YES,CPIBASE=3438,PUBLIC=YES);
%BULLIT(YEAR=2001,REAL=YES,ADJINC=YES,CPIBASE=3438,PUBLIC=YES);
%BULLIT(YEAR=1998,REAL=YES,ADJINC=YES,CPIBASE=3438,PUBLIC=YES);
%BULLIT(YEAR=1995,REAL=YES,ADJINC=YES,CPIBASE=3438,PUBLIC=YES);
%BULLIT(YEAR=1992,REAL=YES,ADJINC=YES,CPIBASE=3438,PUBLIC=YES);
%BULLIT(YEAR=1989,REAL=YES,ADJINC=YES,CPIBASE=3438,PUBLIC=YES);
*/

***************************************************************************;
***************************************************************************;
*ENDSAS;