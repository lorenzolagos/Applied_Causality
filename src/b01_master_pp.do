/*
This do file leaves creates master files for the probationary period analysis

The data is stored under MasterPP_SS.dta files.

Run after clean_RAIS.do
Run before before xxx

NEED A HIRE DATE
*/

cd "C:\Users\lfl2121\Documents\RAIS_new\ProbPer"
*cd "D:\RAIS_new\ProbPer"

cap log close

clear all

set more off
pause on

set seed 12345

local qgrp "RR"
*local qgrp "RR AC AP TO PI SE AL RO PB AM RN MA MS MT ES PA DF CE GO PE BA SC RS PR MG RJ SP"
local m : word count `qgrp'
forvalues k = 1/`m' {
local q : word `k' of `qgrp'

cd "C:\Users\lfl2121\Documents\RAIS_new\ProbPer/`q'"
*cd "D:\RAIS_new\ProbPer/`q'"

log using MasterPP_`q'_log, replace
timer clear 1
timer on 1


local agrp "2003 2004 2005 2006 2007 2008 2009 2010"
local n : word count `agrp'
forvalues i = 1/`n' {
	local a : word `i' of `agrp'

	cd "C:\Users\lfl2121\Documents\RAIS_new\Stata_files/`q'"
	*cd "D:\RAIS_new\Stata_files/`q'"
	
	use `q'`a'_RAIS, clear
	
	cd "C:\Users\lfl2121\Documents\RAIS_new\ProbPer/`q'"
	*cd "D:\RAIS_new\ProbPer/`q'"

***SAMPLE OF IDENTIFIABLE WORKERS***
	
	display "Unidentifiable IDs (drop x==0)"
	gen x=0 if pis<1000 | pis==. | estabid==0 | remunavgnom==0 | remunavgminw==0
	replace x=1 if x==.
	tab x
	drop if x==0
	drop x
	

***GENERATE YEAR***

	gen year = `a' 
	
	keep pis estabid year municipality wagetype hiredate emptype hiretype sepcause tenure concla cnae10class cbo02
	order pis estabid year municipality wagetype hiredate emptype hiretype sepcause tenure concla cnae10class cbo02
	
***REMOVE DUPLICATES***
	
	display "Duplicates (drop x==0)"
	egen x=tag(pis estabid year hiredate)
	tab x
	drop if x==0
	drop x

***SAVE INTO A MASTERPP_SS FILE***

	if `i'==1{
		save MasterPP_`q', replace
	}
	else{
		append using MasterPP_`q'
		save MasterPP_`q', replace
	}
}


***MODE VARIABLES***

sort pis year, stable

egen _mode=mode(concla), by(estabid) minmode
replace concla=_mode
drop _mode

egen _mode=mode(municipality), by(estabid) minmode
replace municipality=_mode
drop _mode

egen _mode=mode(cnae10class), by(estabid) minmode
replace cnae10class=_mode
drop _mode

mdesc


***REHIRED (has been observed at the same firm before hiring)***

sort pis year hiredate, stable

*people being hired within a year who have been observed at firm previously
egen x=tag(pis estabid)
gen rehired=(x==0)&(hiretype!=0)
*people being hired at end of year but fired shortly into next year who have been observed at firm previously
bysort pis: replace x=1 if x[_n]==0 & x[_n-1]==1 & estabid[_n]==estabid[_n-1] & year[_n]==year[_n-1]+1
replace rehired=1 if x==0 & hiretype==0 & tenure<=4
drop x

***REOCCUPATION (has been observed at the same occupation before hiring)***

sort pis year hiredate, stable

*people being hired within a year who have been observed at occupation previously
egen x=tag(pis cbo02)
gen reocc=(x==0)&(hiretype!=0)
*people being hired at end of year but fired shortly into next year who have been observed at occupation previously
bysort pis: replace x=1 if x[_n]==0 & x[_n-1]==1 & estabid[_n]==estabid[_n-1] & year[_n]==year[_n-1]+1
replace reocc=1 if x==0 & hiretype==0 & tenure<=4
drop x


***MAX MODE & MIN***

display "Drop Rehired (drop rehired==1)"
tab rehired
keep if rehired==0

display "Drop Reoccupation (drop reocc==1)"
tab reocc
keep if reocc==0

display "Urban + Open Ended (drop x==0)"
gen x=(emptype==10)
tab x
drop if x==0
drop x

display "End of Contract + Layoff (drop x==0)"
gen x=(sepcause==12) | (sepcause==11)
tab x
drop if x==0
drop x

display "Normal Hire (drop x==0)"
gen x=(hiretype==0) | (hiretype==1) | (hiretype==2)
tab x
drop if x==0
drop x

display "Monthly Wage (drop x==0)"
gen x=(wagetype==1)
tab x
drop if x==0
drop x

display "Drop Anomalous End of Contract (drop x==1)"
gen x=(sepcause==12) & (tenure>3)
tab x
drop if x==1
drop x

display "Drop Unnecessary Layoffs (drop x==1)"
gen x=(sepcause==11) & (tenure>4)
tab x
drop if x==1
drop x

keep pis estabid year sepcause tenure cbo02 concla municipality cnae10class

replace tenure=2.9 if tenure==3 & sepcause==12
replace tenure=1.9 if tenure==2 & sepcause==12
replace tenure=1.4 if tenure==1.5 & sepcause==12
replace tenure=0.9 if tenure==1 & sepcause==12

*maximum tenure for end of contract
egen ppmax=max(tenure) if sepcause==12, by (estabid cbo02 year)
replace tenure=. if round(tenure,0.1)==round(ppmax,0.1) & sepcause==12
*mode tenure for end of contract (ignoring the maximum)
egen ppmode=mode(tenure) if sepcause==12, by (estabid cbo02 year) 
gen x=1 if tenure!=. & sepcause==12
egen countmode=sum(x), by(estabid cbo02 year)
replace countmode=. if countmode==0
drop x
replace tenure=round(ppmax,0.1) if tenure==. & sepcause==12
*minimum tenure for layoff
egen ppmin=min(tenure) if sepcause==11, by (estabid cbo02 year)

bysort estabid cbo02 year sepcause: gen countmax=_N if sepcause==12
bysort estabid cbo02 year sepcause: gen countmin=_N if sepcause==11

***FINALIZE DATASET*** 

*collapse to the (firm, occupation, year) level
collapse ppmax ppmode ppmin countmax countmode countmin, by(estabid cbo02 year concla municipality cnae10class)

replace countmax=0 if countmax==.
replace countmode=0 if countmode==.
replace countmin=0 if countmin==.
replace ppmax=round(ppmax,0.1)
replace ppmode=round(ppmode,0.1)

*descriptive statistics
mdesc
tab ppmax
tab ppmax if countmax>=5
tab ppmode
tab ppmode if countmax>=5 & countmode>=5

*end of contract (top)... must have 5 observations
gen pp1obs=ppmax if countmax>=5
gen pp1=pp1obs
replace pp1=3 if pp1obs>2 & pp1obs<=3
replace pp1=2 if pp1obs>1.5 & pp1obs<=2
replace pp1=1.5 if pp1obs>1 & pp1obs<=1.5
replace pp1=1 if pp1obs<=1

*end of contract (middle)... must have 5 observations
gen pp2obs=ppmode if countmax>=5 & countmode>=5
gen pp2=pp2obs
replace pp2=3 if pp2obs>2 & pp2obs<=3
replace pp2=2 if pp2obs>1.5 & pp2obs<=2
replace pp2=1.5 if pp2obs>1 & pp2obs<=1.5
replace pp2=1 if pp2obs<=1
replace pp2=0 if pp2obs==. & countmax>=5

sort cbo02 estabid year
order cbo02 estabid year municipality concla cnae10class pp1 pp2 pp1obs pp2obs

keep if pp1!=.

label data "RAIS `q' Master PP"
save MasterPP_`q', replace

timer off 1
timer list 1
log close

}
