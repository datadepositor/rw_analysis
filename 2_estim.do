********************************************************************************
*				program: estimation			   	   							   *
*				author: Ziwei Ye 											   *
********************************************************************************
global root "D:"
global data "$root\1_data"

global ds "$data\source"
global dd "$data\data"
global dp "$data\proc"
global dt "$data\temp"
global do "$data\out"



*====================used================================*
use $dd\estim.dta, clear

gl vce scrd 

gen pattern=lag1_cornshare_corn/100 
gen hist_scrd=lag2_rw_rate_scrd/100
gen hist_countyfips=lag2_rw_rate_countyfips/100
for var prec*: gen Xcm=X/100

keep if year>=2005 & year<=2016
drop if soilinsecticide==1 

egen miss=rowmiss(rw bt hist_scrd prec4cm-prec8cm tMin1 tMin2 tMin3 seedtrtrate locid)
drop if miss

save $dd\estim_use.dta, replace

*-------------------regression analysis----------------------------------------*
//fractional response model
//crd level Bt rate
qui fracreg probit rw  i.bt##c.hist_scrd prec4cm-prec8cm i.seedtrtrate  i.locid , vce(cluster $vce)
est store f010
margins, dydx(*) post
est store f01

qui fracreg probit rw  i.bt##c.hist_scrd prec4cm-prec8cm i.seedtrtrate pattern i.locid, vce(cluster $vce)
est store f020
margins, dydx(*) post
est store f02

qui fracreg probit rw  i.bt##c.hist_scrd prec4cm-prec8cm i.seedtrtrate tMin1 tMin2 tMin3 i.locid, vce(cluster $vce)
est store f030
margins, dydx(*) post
est store f03


esttab  f010 f020 f030 using $do\results_v6.rtf, ///
	star(* 0.10 ** 0.05 *** 0.01) ///
	keep(1.bt hist_scrd 1.bt#c.hist_scrd prec4cm prec5cm prec6cm prec7cm prec8cm 1.seedtrtrate pattern tMin1 tMin2 tMin3) ///
	title ("Table S2. Fractional response model results (coefficient estimates) for root injury, using Crop Reporting District level Bt planting rate") ///
	se b(%9.2f) label ///
	replace
	
//county level Bt rate
qui fracreg probit rw  i.bt##c.hist_countyfips prec4cm-prec8cm i.seedtrtrate  i.locid, vce(cluster $vce)
est store f110
margins, dydx(*) post
est store f11

qui fracreg probit rw  i.bt##c.hist_countyfips prec4cm-prec8cm i.seedtrtrate  pattern i.locid, vce(cluster $vce)
est store f120
margins, dydx(*) post
est store f12

qui fracreg probit rw i.bt##c.hist_countyfips prec4cm-prec8cm i.seedtrtrate tMin1 tMin2 tMin3 i.locid, vce(cluster $vce)
est store f130
margins, dydx(*) post
est store f13

esttab  f110  f120  f130  using $do\results_v6.rtf, ///
	star(* 0.10 ** 0.05 *** 0.01) ///
	keep(1.bt hist_countyfips 1.bt#c.hist_countyfips prec4cm prec5cm prec6cm prec7cm prec8cm 1.seedtrtrate pattern tMin1 tMin2 tMin3) ///
	title ("Table S3. Fractional response model results (coefficient estimates) for root injury, using county-level Bt planting rate") ///
	se b(%9.2f) label ///
	append

//linear comparison
//crd level Bt rate
qui reg rw  i.bt##c.hist_scrd prec4cm-prec8cm i.seedtrtrate  i.locid, vce(cluster $vce)
est store n01

qui reg rw  i.bt##c.hist_scrd prec4cm-prec8cm i.seedtrtrate pattern  i.locid, vce(cluster $vce)
est store n02

qui reg rw  i.bt##c.hist_scrd prec4cm-prec8cm i.seedtrtrate tMin1 tMin2 tMin3 i.locid, vce(cluster $vce)
est store n03

esttab  n01 n02  n03 using $do\results_v6.rtf, ///
	star(* 0.10 ** 0.05 *** 0.01) ///
	keep(1.bt hist_scrd 1.bt#c.hist_scrd prec4cm prec5cm prec6cm prec7cm prec8cm 1.seedtrtrate pattern tMin1 tMin2 tMin3) ///
	title ("Table S4. Linear regression results for root injury, using Crop Reporting District level Bt planting rate") ///
	se b(%9.2f) label ///
	append

//county-level Bt rate
qui reg rw  i.bt##c.hist_countyfips prec4cm-prec8cm i.seedtrtrate  i.locid, vce(cluster $vce)
est store n11

qui reg rw  i.bt##c.hist_countyfips prec4cm-prec8cm i.seedtrtrate pattern i.locid, vce(cluster $vce)
est store n12

qui reg rw  i.bt##c.hist_countyfips prec4cm-prec8cm i.seedtrtrate tMin1 tMin2 tMin3 i.locid, vce(cluster $vce)
est store n13

esttab  n11 n12  n13 using $do\results_v6.rtf, ///
	star(* 0.10 ** 0.05 *** 0.01) ///
	keep(1.bt hist_countyfips 1.bt#c.hist_countyfips prec4cm prec5cm prec6cm prec7cm prec8cm 1.seedtrtrate pattern tMin1 tMin2 tMin3) ///
	title ("Table S5. Linear regression results for root injury, using county-level Bt planting rate") ///
	se b(%9.2f) label ///
	append	

	
*------------predicted root injury table in the main manuscript (Table 1)------*
*using estimates from the fractional response model with crd level Bt planting rate
est restore f010
margins bt, at(hist_scrd=(0 1)) post
lincom _b[1._at#1.bt]-_b[1._at#0.bt]
lincom _b[2._at#1.bt]-_b[2._at#0.bt]
lincom _b[2._at#1.bt]-_b[1._at#1.bt]
lincom _b[2._at#0.bt]-_b[1._at#0.bt]
lincom _b[2._at#1.bt]-_b[2._at#0.bt]-(_b[1._at#1.bt]-_b[1._at#0.bt])






