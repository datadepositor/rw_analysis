********************************************************************************
*				program: policy analysis						       		   *
*				author: Ziwei Ye											   *
********************************************************************************
global root "D:"
global data "$root\1_data"

global ds "$data\source"
global dd "$data\data"
global dp "$data\proc"
global dt "$data\temp"
global do "$data\out"

*====================initialize fixed parameters===============================*
 
*****parameters, fixed
global a=0.47 //bt treatment effect 
gl b=0.33 //selection effect
gl r=0.30 //suppression effect


*=====================compile parameters from dataset===========================*
*****maize cash price, state level
use $dd\estim.dta, clear
keep state* east
duplicates drop state, force
expand 18
bys state: gen year=2003+_n-1 //because trial field rw in 2004 represents the area rw in 2003.
save $dt\sample_statefips.dta, replace

use D:\Dropbox\1_Research\8_climate_biotech\1_data\temp\cash.dta, clear
//evelator-level maize cash price data; see "Data and materials availability" for data availability

duplicates drop year statefips, force
keep year statefips cash_statefips cash_na
merge 1:1 year statefips using $dt\sample_statefips.dta, keep(match using) nogen
save $dt\par_pcash.dta, replace
//aggregated state-level maize cash price data; see "Data and materials availability" for data availability

*****obtain yield potential, seed premium, Bt rate, and m.
use $dd\estim.dta, clear
keep if soilinsecticide==0 &seedtrtrate==0 

gen premium_acre_statefips=ps_acre_statefips1-ps_acre_statefips0
keep year east state* countyfips  bt rw premium_acre* *rw_rate_scrd* def16
collapse (mean) premium_acre* *rw_rate* rw def16, by(year east state* countyfips bt)
reshape wide rw, i(year east stateabbr statefips countyfips *rw_rate_scrd premium_acre_statefips def16) j(bt)

xtset countyfips year
for var rw0 rw1: rename X X_tminus1
gen rw0=f.rw0_tminus1
gen rw1=f.rw1_tminus1
drop rw*_tminus1

keep if year>=2014 & year<=2016
drop if rw0==. & rw1==.

replace rw1=rw0-${a}+${b}*lag1_rw_rate_scrd/100 if rw1==. 
replace rw0=rw1+${a}-${b}*lag1_rw_rate_scrd/100 if rw0==. 
for var rw0 rw1: replace X=0 if X<0

save $dt\threeyear.dta, replace



//match with county level yield.
import delimited $ds\yield_NASS\yield_2010to16_county.csv, clear
//county-level maize yield data from NASS, USDA; see "Data and materials availability" for data availability

keep if period=="YEAR"
keep year stateansi countyansi value
rename value yield
gen countyfips=stateansi*1000+countyansi
drop if countyfips==.
drop *ansi

merge 1:1 countyfips year using $dt\threeyear.dta, keep(match using) nogen

//back out yield potential
gen rw_rate=rw_rate_scrd/100
gen rw_loss=(rw_rate*rw1+(100-rw_rate)*rw0)/100 
gen yield_potential=yield/(1-0.45*rw_loss)

preserve
collapse (mean) mean1= yield mean2=yield_potential ///
	(sd) sd1= yield sd2=yield_potential ///
	(count) n1=yield n2=yield_potential ///
	, by(state stateabbr statefips)
save $do\yield_state.dta, replace
restore

preserve
collapse (mean) mean1= yield mean2=yield_potential ///
	(sd) sd1= yield sd2=yield_potential ///
	(count) n1=yield n2=yield_potential ///
	, by(east)
tostring east, force replace
replace east="West" if east=="0"
replace east="East" if east=="1"
rename east state
save $do\yield_east.dta, replace
restore


//calculate c/m
merge m:1 year statefips using $dt\par_pcash.dta, keep(match master) nogen

replace cash_statefips=cash_na if cash_statefips==.
gen m=3*0.15*cash_statefips*yield_potential/def16
gen c=premium_acre_statefips/def16
gen cbym=c/m
collapse (mean) c m cash_statefips yield_potential cbym rw_rate, by(east)
sort east


//store c/m for east and west.
gl cbym0=cbym[1] //+${b}-${r} - the two terms are added in the old version for benefit calculation
gl cbym1=cbym[2] //+${b}-${r}
di "$cbym0, $cbym1"

save $dt\par.dta, replace


*====================initialize changing parameters and find equilibriums===============================*
loc taolist "0.70 0.78 0.95"
loc ulist "0.05 0.15 0.25 0.50"
loc intlist "private group"

foreach int of local intlist{
	foreach u of local ulist{
	foreach tao of local taolist{
	
	
macro drop u tao int
gl u=`u' //choose spillover parameter
gl tao=`tao' //choose discount factor


gl btilda=${b}*${u}/(1-${u})
gl rtilda=${r}*${u}/(1-${u})


*****calculate k(P) and g(P)

forv i=0/1{
gl i=`i' //region 




//adjust C to obtain solutions to private interest versus group interest.
if "`int'"=="private"{
gl C`i'=${cbym${i}}+${tao}*(${b}-${r}) 
}

if "`int'"=="group"{
gl C`i'=${cbym${i}}+${tao}*(${b}-${r}) +${tao}*(${btilda}-${rtilda}) //add the second part for group interest solution. 
}

//set up 
clear
set obs 101
gen P=(_n-1)/100


//k(P): match with ecdf to find out the approximated F value for each argument value Fx

gen Fx_0=. //the x inside F() for s(t-1)=0
replace Fx_0=${C`i'}+${rtilda}*P
gen Fx_1=. //the x inside F() for s(t-1)=1
replace Fx_1=${C`i'}+${rtilda}*P+${r}


gen id=_n
expand 100
sort id
bys id: gen percentile=_n 
merge m:1 percentile using $dt/ecdf.dta, nogen
sort  id percentile

gen dif_0=abs(rw${i}-Fx_0)
gen dif_1=abs(rw${i}-Fx_1)
bysort id (dif_0): gen is_minimum_0 = dif_0 == dif_0[1] 
//for each P, find the percentile whose argument is the closest to the CDF argument Fx by comparing the absolute difference between Fx and rw
bysort id (dif_1): gen is_minimum_1 = dif_1 == dif_1[1] 

preserve
keep if is_minimum_0==1
keep id percentile cdf${i}
collapse (max) percentile cdf${i}, by(id)
for var percentile cdf${i}: rename X X_0
save $dt\temp_0.dta, replace
restore
preserve
keep if is_minimum_1==1
keep id percentile cdf${i}
collapse (max) percentile cdf${i}, by(id)
for var percentile cdf${i}: rename X X_1
save $dt\temp_1.dta, replace
restore

keep id P Fx_*
duplicates drop id, force
merge 1:1 id using $dt\temp_0.dta, nogen
merge 1:1 id using $dt\temp_1.dta, nogen

gen fn_0=0 //the indicator function for s(t-1)=0
replace fn_0=1 if ${a}-${btilda}*P- ${C`i'}>=0
gen fn_1=0 //the indicator function for s(t-1)=1
replace fn_1=1 if ${a}-${btilda}*P-${b}-${C`i'}>=0

gen k${i}=(1-P)*fn_0*(1-cdf${i}_0)+P*fn_1*(1-cdf${i}_1)

//g(P)
gen g=P

//find P s.t. g=k
gen dif=abs(k${i}-g)
sort dif
gen is_minimum${i}=1 in 1
scalar Pstar=P in 1
di Pstar

keep P k${i} g is_minimum${i}
save $dt\temp${i}.dta, replace
}


use $dt\temp0.dta, clear
merge 1:1 P using $dt\temp1.dta, nogen
save $dt\function.dta, replace

//figure for (tao, u, int) combination


loc TAO=100*${tao}
loc U=100*${u}

	use  $dt\function.dta, clear
	sort is_minimum0
	scalar Pstar0=P in 1
	sort is_minimum1
	scalar Pstar1=P in 1

	use  $dt\function.dta, clear
	sort is_minimum0
	scalar Pstar0=P in 1
	sort is_minimum1
	scalar Pstar1=P in 1
	graph twoway (line g P , sort lc(black) lw(0.5)) ///
		(line k1 P , sort lc(red*1.5)  lp(shortdash)  lw(0.5)) ///
		(line k0 P , sort lc(green*1.5) lp(longdash)  lw(0.5)) ///
		(scatteri `=Pstar0' `=Pstar0' , legend(off)  mcolor(red) ) ///
		(scatteri `=Pstar1' `=Pstar1' , legend(off)  mcolor(red) ), ///
		xlabel(0 (0.2) 1.0, format(%03.1f))  ///
		ylabel(0 (0.2) 1.0, format(%03.1f))  ///
		title("{bf:u=`U'%}", size(medium)) /// 
		leg(region(lwidth(none) color(none)) on order(3 2 1) ring(0) position(12) row(1)  size(small) label(1 "g(P)") label(2 "k(P) East") label(3 "k(P) West" )) ///
		xtitle("Bt planting rate (P)", size(medium) margin(t=3)) /// 
		ytitle("Functions")   ///
		saving($dt\functions_tao`TAO'_u`U'_`int'.gph, replace) nodraw

	
	
		
*====================cost and benefit analysis==============================================*
scalar drop _all

*****status-quo(Pbar) and equilibrium(Pstar) coverage from previous section 
use  $dt\function.dta, clear
sort is_minimum0
scalar Pstar0=P in 1
sort is_minimum1
scalar Pstar1=P in 1

use $dt\par.dta, clear
sort east
scalar Pbar0=rw_rate in 1
scalar Pbar1=rw_rate in 2

di Pstar0
di Pstar1
di Pbar0
di Pbar1



*****calculate cost and benefit for one period 
//OPTIMUM

forv rg=0/1{
	local scenario "star"
	foreach sc of local scenario{ 
	
	//obtain F(C+rtilda*P) and F(C+rtilda*P+r) 
	scalar temp01_`sc'`rg'=${C`rg'}+${rtilda}*P`sc'`rg'
	scalar temp11_`sc'`rg'=${C`rg'}+${rtilda}*P`sc'`rg'+${r}
	
	use $dt\ecdf.dta, clear
	gen dif=abs(rw`rg'-temp01_`sc'`rg')
	sort dif
	gen is_minimum`rg'=1 in 1
	scalar F01_`sc'`rg'=cdf`rg' in 1 //F(C+rtilda*P)
	
	use $dt\ecdf.dta, clear
	gen dif=abs(rw`rg'-temp11_`sc'`rg')
	sort dif
	gen is_minimum`rg'=1 in 1
	scalar F11_`sc'`rg'=cdf`rg' in 1 //F(C+rtilda*P+r)
	
	//obtain P(st=1,s(t-1)=0) and P(st=1,s(t-1)=1)
	scalar P01_`sc'`rg'=cond(${a}-${btilda}*P`sc'`rg'>=${C`rg'},1,0)*(1-F01_`sc'`rg')
	scalar P11_`sc'`rg'=cond(${a}-${btilda}*P`sc'`rg'-${b}>=${C`rg'},1,0)*(1-F11_`sc'`rg')
	
	//obtain the probability weight for each combinations of s(t-1) and st
	scalar Pw00_`sc'`rg'=(1-P`sc'`rg')*(1-P01_`sc'`rg') //s(t-1)=0,st=0
	scalar Pw01_`sc'`rg'=(1-P`sc'`rg')*P01_`sc'`rg' //s(t-1)=0,st=1
	scalar Pw10_`sc'`rg'=P`sc'`rg'*(1-P11_`sc'`rg') //s(t-1)=1,st=0
	scalar Pw11_`sc'`rg'=P`sc'`rg'*P11_`sc'`rg' //s(t-1)=1,st=1
	
	//obtain unconstrained benefit for each combinations of s(t-1) and st
	scalar phiP00_`sc'`rg'=${rtilda}*P`sc'`rg'
	scalar phiP01_`sc'`rg'=${a}+${rtilda}*P`sc'`rg'-${btilda}*P`sc'`rg'
	scalar phiP10_`sc'`rg'=${r}+${rtilda}*P`sc'`rg'
	scalar phiP11_`sc'`rg'=${a}+${r}-${b}+${rtilda}*P`sc'`rg'-${btilda}*P`sc'`rg'
	
	
	loc sslist "00 01 10 11"
	foreach ss of local sslist{
		use $dt\ecdf.dta, clear
		
	//obtain F(phi(P)) for each combinations of s(t-1) and st
		gen dif=abs(rw`rg'-phiP`ss'_`sc'`rg')
		sort dif
		gen is_minimum`rg'=1 in 1
		scalar FphiP`ss'_`sc'`rg'=cdf`rg' in 1 //F(phi(P)) for region=i
		scalar pc`ss'_`sc'`rg'=percentile in 1 //the integer percentile corresponding to phiP, to be used for graphical approximation in the next step
	
	//obtain F(v) integral over 0 to phi(P)
		drop dif is_minimum
		sort percentile

		tsset percentile
		gen l_rw`rg'=l.rw`rg'
		gen l_cdf`rg'=l.cdf`rg'
		replace l_rw`rg'=0 in 1
		replace l_cdf`rg'=0 in 1

		scalar integral`ss'_`sc'`rg'=0
		forv i=1(1)`=pc`ss'_`sc'`rg''{
			sort percentile
			loc a1=rw`rg' in `i'
			loc a2=l_rw`rg' in `i'
			loc b1=cdf`rg' in `i'
			loc b2=l_cdf`rg' in `i'

			scalar integral`ss'_`sc'`rg'=integral`ss'_`sc'`rg'+(`b1'+`b2')*(`a1'-`a2')/(2)

		}
		
	//obtain Emin
		scalar Emin`ss'_`sc'`rg'=phiP`ss'_`sc'`rg'-integral`ss'_`sc'`rg'
		 
	}
	
	//obtain expected payoff (in terms of root injury) over v
	scalar E_`sc'`rg'=Emin00_`sc'`rg'*Pw00_`sc'`rg'+Emin01_`sc'`rg'*Pw01_`sc'`rg'+Emin10_`sc'`rg'*Pw10_`sc'`rg'+Emin11_`sc'`rg'*Pw11_`sc'`rg'-${cbym`rg'}*P`sc'`rg'

	//obtain expected payoff (in terms of money, that is, Benefit) for one period
	use $dt\par.dta, clear
	sort east
	loc temp=`rg'+1
	scalar Benefit_`sc'`rg'=m[`temp']*E_`sc'`rg'
    }

//STATUS QUO
	local scenario "bar"
	foreach sc of local scenario{ 
	
	
	//obtain the probability weight for each combinations of s(t-1) and st
	scalar Pw00_`sc'`rg'=(1-P`sc'`rg')*(1-P`sc'`rg') //s(t-1)=0,st=0
	scalar Pw01_`sc'`rg'=(1-P`sc'`rg')*P`sc'`rg' //s(t-1)=0,st=1
	scalar Pw10_`sc'`rg'=P`sc'`rg'*(1-P`sc'`rg') //s(t-1)=1,st=0
	scalar Pw11_`sc'`rg'=P`sc'`rg'*P`sc'`rg' //s(t-1)=1,st=1
	
	//obtain unconstrained benefit for each combinations of s(t-1) and st
	scalar phiP00_`sc'`rg'=${rtilda}*P`sc'`rg'
	scalar phiP01_`sc'`rg'=${a}+${rtilda}*P`sc'`rg'-${btilda}*P`sc'`rg'
	scalar phiP10_`sc'`rg'=${r}+${rtilda}*P`sc'`rg'
	scalar phiP11_`sc'`rg'=${a}+${r}-${b}+${rtilda}*P`sc'`rg'-${btilda}*P`sc'`rg'
	
	
	loc sslist "00 01 10 11"
	foreach ss of local sslist{
		use $dt\ecdf.dta, clear
		
	//obtain F(phi(P)) for each combinations of s(t-1) and st
		gen dif=abs(rw`rg'-phiP`ss'_`sc'`rg')
		sort dif
		gen is_minimum`rg'=1 in 1
		scalar FphiP`ss'_`sc'`rg'=cdf`rg' in 1 //F(phi(P)) for region=i
		scalar pc`ss'_`sc'`rg'=percentile in 1 //the integer percentile corresponding to phiP, to be used for graphical approximation in the next step
	
	//obtain F(v) integral over 0 to phi(P)
		drop dif is_minimum
		sort percentile

		tsset percentile
		gen l_rw`rg'=l.rw`rg'
		gen l_cdf`rg'=l.cdf`rg'
		replace l_rw`rg'=0 in 1
		replace l_cdf`rg'=0 in 1

		scalar integral`ss'_`sc'`rg'=0
		forv i=1(1)`=pc`ss'_`sc'`rg''{
			sort percentile
			loc a1=rw`rg' in `i'
			loc a2=l_rw`rg' in `i'
			loc b1=cdf`rg' in `i'
			loc b2=l_cdf`rg' in `i'

			scalar integral`ss'_`sc'`rg'=integral`ss'_`sc'`rg'+(`b1'+`b2')*(`a1'-`a2')/(2)

		}
		
	//obtain Emin
		scalar Emin`ss'_`sc'`rg'=phiP`ss'_`sc'`rg'-integral`ss'_`sc'`rg'
		 
	}
	
	//obtain expected payoff (in terms of root injury) over v
	scalar E_`sc'`rg'=Emin00_`sc'`rg'*Pw00_`sc'`rg'+Emin01_`sc'`rg'*Pw01_`sc'`rg'+Emin10_`sc'`rg'*Pw10_`sc'`rg'+Emin11_`sc'`rg'*Pw11_`sc'`rg'-${cbym`rg'}*P`sc'`rg'

	//obtain expected payoff (in terms of money, that is, Benefit) for one period
	use $dt\par.dta, clear
	sort east
	loc temp=`rg'+1
	scalar Benefit_`sc'`rg'=m[`temp']*E_`sc'`rg'
    }
}


	
*****obtain matrix of Benefit by T, tao, u and region. Note equilibirum coverage is dependent on tao, u and regin. 
//obtain the discount multiplier for tao and T period
clear
loc length "10"
set obs `length'
egen T=seq(),f(1) t(`length')
gen multiplier=.

forv T=1/`length'{
	replace multiplier=(1-${tao}^`T')/(1-${tao}) in `T'
}

//obtain cost and benefit for the range of T
forv rg=0/1{

	local scenario "star bar"
	foreach sc of local scenario{ 

	gen Benefit_`sc'`rg'=`=Benefit_`sc'`rg''*multiplier
	gen P`sc'`rg'=`=P`sc'`rg''
	}
	
	gen DIFF_`rg'=Benefit_star`rg'-Benefit_bar`rg'
    }


loc TAO=100*${tao}
loc U=100*${u}
gen TAO=`TAO'
gen U=`U'
gen interest="`int'"
save $dt\net_tao`TAO'_u`U'_`int'.dta, replace

}
}
}

*****obtain the complete dataset for cost-benefit analysis
clear
loc taolist "0.70 0.78 0.95"
loc ulist "0.05 0.15 0.25 0.50"
loc intlist "private group"

foreach int of local intlist{
	foreach u of local ulist{
	foreach tao of local taolist{
		loc TAO=100*`tao'
		loc U=100*`u'
		append using $dt\net_tao`TAO'_u`U'_`int'.dta

	}
	}
}
save $dt\net_figure.dta, replace

*==============================outputs: figures==================================*
graph set window fontface "Helvetica"
global psize "10pt" //panel text size
global tsize "7pt" //figure test size

*------Fig.3 regional heterogeneity figure----------------------------*
//CDF figure
use $dt\ecdf.dta, clear
graph twoway (line cdf0 rw0 , sort(cdf0) lc("162 92 115") lp(longdash)  lw(0.6)) ///
	(line cdf1 rw1 , sort(cdf1) lc("119 143 165") lp(shortdash)  lw(0.6)), ///
	xlabel(0 (0.2) 1.0, format(%03.1f) labsize("$tsize"))  ylabel(0 (0.2) 1.0, format(%03.1f)  labsize("$tsize")) ///
	leg(region(lwidth(none)) on order(1 2) ring(0) position(5) col(1)  size("$tsize") label(1 "West") label(2 "East" )) ///
	xtitle("{bf:Root injury}", size("$tsize") margin(t=3)) /// 
	ytitle("{bf:Cumulative distribution function}", margin(r=3) size("$tsize"))  ///
	title("{bf:A}",size("$psize") margin(b=3) ring(1) position(11)) ///
	saving($do\ecdf.gph, replace) 


//yield potential figure
use $do\yield_state.dta, clear
append using $do\yield_east.dta
save $dt\temp.dta, replace

keep if state=="West"|state=="East"

forv i=1/2{
	gen hi`i'=mean`i'+invttail(n`i'-1,0.025)*(sd`i' / sqrt(n`i'))
	gen lo`i'=mean`i'-invttail(n`i'-1,0.025)*(sd`i' / sqrt(n`i'))
}

reshape long mean sd n hi lo, i(state) j(yieldvar) //yieldvar=1 for yield, =2 for yield_potential
sort statefips

save $dt\temp.dta, replace


loc bwidth=3
loc xlab1=1+`bwidth'/2
loc xlab2=`xlab1'+`bwidth'
sort state yieldvar
loc m1=round(mean[1],0.1) //yield, east
if round(`m1', 0.1) == int(`m1') {
    local m1 `m1'.0
}
di "`m1'"
loc m2=round(mean[2],0.1) //yield potential, east
loc m3=round(mean[3],0.1) //yield, west
loc m4=round(mean[4],0.1) //yield potential, west

loc d1=round(mean[2]-mean[1],0.1)
di "`d1'"
loc d0=round(mean[4]-mean[3],0.1)
di "`d0'"

use $dt\temp.dta, clear
gen statevar=yieldvar+1 if state=="West"
replace statevar=yieldvar+1+`bwidth' if state=="East"
graph twoway ///
	(bar mean statevar if yieldvar==1 & state=="West",  barw(1) fcolor("195 149 164") lcolor("162 92 115") lwidth(0.6)) ///
	(bar mean statevar if yieldvar==1 & state=="East",  barw(1) fcolor("169 184 198") lcolor("119 143 165") lwidth(0.6) ) ///
	(bar mean statevar if yieldvar==2 & state=="West", barw(1) fcolor("250 230 240") lcolor("195 149 164") lwidth(0.6)) ///
	(bar mean statevar if yieldvar==2 & state=="East", barw(1) fcolor("231 239 250")  lcolor("169 184 198") lwidth(0.6) ) ///
	(rcap hi lo statevar if yieldvar==1 & state=="West", color("162 92 115") lwidth(0.6)) ///
	(rcap hi lo statevar if yieldvar==2 & state=="West", color("195 149 164") lwidth(0.6)) ///
	(rcap hi lo statevar if yieldvar==1 & state=="East", color("119 143 165") lwidth(0.6)) ///
	(rcap hi lo statevar if yieldvar==2 & state=="East", color("169 184 198") lwidth(0.6)) ///
	, ///
	text(285 2.5 "{bf:Loss=`d0'}",color("162 92 115") size("$tsize") box bcolor("250 230 240") margin(medsmall)) ///
	text(285 5.5 "{bf:Loss=`d1'}",color("119 143 165") size("$tsize") box bcolor("231 239 250") margin(medsmall)) ///
	text(200 1.8 "Y=`m3'", color("162 92 115") size("$tsize")) ///
	text(257 3.0 "P=`m4'", color("195 149 164") size("$tsize")) ///
	text(202 4.8 "Y=`m1'", color("119 143 165") size("$tsize")) ///	
	text(220 6.0 "P=`m2'", color("169 184 198") size("$tsize")) ///	
	xscale(range(1 7)) yscale(range(100 300)) ///
	legend(region(lstyle(none)) off) ///
	xlabel( `xlab1' "West"  `xlab2' "East",labsize("$tsize") ) ///
	ylabel(100 150 200 250,labsize("$tsize") ) ///
	xtitle("{bf:Region}", size("$tsize") margin(t=3)) /// 
	ytitle("{bf:Yield (bushel/acre)}", size("$tsize") margin(r=2)) ///
	title("{bf:B}",size("$psize") margin(b=3) ring(1) position(11)) ///
	saving($do/yield_compare.gph, replace)

cd $do
graph combine ecdf.gph yield_compare.gph, row(1)  iscale(1) ///
	xsize(4.75) ysize(2.5)
	
graph export $do\hetero.tif, replace width(1425) //horizontal pixels=4.75inch*300dots per inch
graph export $do\hetero.png, replace width(1425)
	



*----------Fig.4 equilibrium optimmum and cost-benefit in one figure-----------------------*

*****obtain dataset for making one figure
use $dt\net_figure.dta, clear
sort U TAO T
order U, before(T)
order TAO, after(U)
keep if T==1|T==5|T==10
keep P* TAO U T interest DIFF_*

reshape wide P* DIFF*, i(TAO U T) j(interest ) string
sort   T U TAO
gen u=U/100
gen tao=TAO/100
drop  Pbar0group Pbar1group
order U TAO T Pbar0private Pstar0private Pstar0group Pbar1private Pstar1private Pstar1group ///
	DIFF_0private DIFF_0group DIFF_1private DIFF_1group
save $dt\temp.dta, replace

use $dt\temp.dta, clear
keep U TAO T DIFF_0private DIFF_0group DIFF_1private DIFF_1group
tostring T, force replace
reshape wide DIFF*, i(TAO U ) j(T) string
save $dt\temp1.dta, replace

use $dt\temp.dta, clear
keep  U TAO u tao Pbar0private Pstar0private Pstar0group Pbar1private Pstar1private Pstar1group
duplicates drop  U TAO, force
merge 1:1 U TAO using $dt\temp1.dta, nogen

gen tao_private=tao-0.015
gen tao_group=tao+0.015
	
save $dt\one_figure.dta, replace

*****figure for Pvalues
set scheme s1color
loc lw "0.3"

//west
use $dt\one_figure.dta, clear
sort U TAO
local bar=Pbar0private[1]
local up1=Pstar0private[1]
local up2=Pstar0private[2]
local up3=Pstar0private[3]
local lp1=Pstar0private[10]
local lp2=Pstar0private[11]
local lp3=Pstar0private[12]

local ug1=Pstar0group[1]
local ug2=Pstar0group[2]
local ug3=Pstar0group[3]
local lg1=Pstar0group[10]
local lg2=Pstar0group[11]
local lg3=Pstar0group[12]

twoway	(scatteri `bar' 0.65 `bar' 1.0, recast(line) lpattern(shortdash) lcolor(gray)) ///
		(scatteri `up1' 0.685 `lp1' 0.685,recast(line)  lc(red*1.5) lw(`lw')) ///
		(scatteri `ug1' 0.715 `lg1' 0.715,recast(line)  lc(black) lw(`lw')) ///
		(scatteri `up2' 0.765 `lp2' 0.765,recast(line) lc(red*1.5) lw(`lw')) ///
		(scatteri `ug2' 0.795 `lg2' 0.795,recast(line) lc(black) lw(`lw')) ///
		(scatteri `up3' 0.935 `lp3' 0.935,recast(line) lc(red*1.5) lw(`lw')) ///
		(scatteri `ug3' 0.965 `lg3' 0.965,recast(line) lc(black) lw(`lw')) ///
		(scatter Pstar0private tao_private if U==5, msymbol(circle) mfcolor(white) lcolor(black) mcolor(red*1.5) msize(medium) mlw(`lw') ) ///
		(scatter Pstar0private tao_private if U==15, msymbol(smtriangle) mfcolor(white) lcolor(black) mcolor(red*1.5) mlw(`lw') ) ///
		(scatter Pstar0private tao_private if U==25, msymbol(X) mfcolor(white) lcolor(black) mcolor(red*1.5) mlw(`lw')) ///
		(scatter Pstar0private tao_private if U==50, msymbol(smdiamond) mfcolor(white) lcolor(black) mcolor(red*1.5) mlw(`lw'))  ///
		(scatter Pstar0group tao_group if U==5, msymbol(circle) mfcolor(white)  lcolor(black) mcolor(black) mlw(`lw')) ///
		(scatter Pstar0group tao_group if U==15, msymbol(smtriangle) mfcolor(white) lcolor(black) mcolor(black) mlw(`lw')) ///
		(scatter Pstar0group tao_group if U==25, msymbol(X) mfcolor(white) lcolor(black) mcolor(black) mlw(`lw')) ///
		(scatter Pstar0group tao_group if U==50, msymbol(smdiamond) mfcolor(white) lcolor(black) mcolor(black) mlw(`lw')) ///
		(rbar  DIFF_0private10 DIFF_0private1 tao_private if U==25, barw(0.025) fcolor("255 204 204") lw(0) yaxis(2)) ///
		(scatter DIFF_0private1 tao_private if U==25, msymbol(P) mcolor("255 153 153") msize(tiny) yaxis(2) ) ///
		(scatter DIFF_0private5 tao_private if U==25, msymbol(P) mcolor("255 0 0") msize(tiny) yaxis(2)) ///
		(scatter DIFF_0private10 tao_private if U==25, msymbol(P) mcolor("153 0 0") msize(tiny) yaxis(2) ) ///
		(rbar  DIFF_0group10 DIFF_0group1 tao_group if U==25, barw(0.025) fcolor("240 240 240") lw(0)  yaxis(2)) ///
		(scatter DIFF_0group1 tao_group if U==25, msymbol(P) mcolor("224 224 224") msize(tiny) yaxis(2)) ///
		(scatter DIFF_0group5 tao_group if U==25, msymbol(P) mcolor("128 128 128") msize(tiny) yaxis(2)) ///		
		(scatter DIFF_0group10 tao_group if U==25, msymbol(P) mcolor("32 32 32") msize(tiny) yaxis(2)) ///
		, ///
		yline(0, lcolor(black)) ///
		xlabel(0.70 0.78 0.95, format(%03.2f) labsize("$tsize")) ///
		ylabel(0 (0.2) 0.8, format(%03.1f) labsize("$tsize")) ///
		ylabel(none, axis(2)) ///
		yscale(range(-0.4 0.8) axis(1)) yscale(range(0 200) axis(2)) ///
		leg(region(lwidth(none) color(none)) on ring(0) position(6) col(1)  size("$tsize") ///
			order(1 "{bf:Status quo}" ///
			2 "{bf:Individual optimum under}" "{bf: (group optimum in black)}" ///
			8 "{it:u}=5%" 9 "{it:u}=15%" 10 "{it:u}=25%" 11 "{it:u}=50%" ///
			16 "{bf:Benefits of shifting to}" "{bf:individual optimum over}" "{bf:(group optimum in gray)}" ///
			17 "1 year" 18 "5 years"  19 "10 years")) ///
		title("{bf:A.West}",size("$psize") margin(b=3) position(11)) ///
		ytitle("{bf:                           Bt planting rate}",size("$tsize") margin(r=3)) ///
		xtitle("{bf:Discount factor}",size("$tsize") margin(t=3)) ///
		saving($dt\Pwest.gph, replace) ///
		fxsize(80)
		
//east
use $dt\one_figure.dta, clear
sort U TAO
local bar=Pbar1private[1]
local up1=Pstar1private[1]
local up2=Pstar1private[2]
local up3=Pstar1private[3]
local lp1=Pstar1private[10]
local lp2=Pstar1private[11]
local lp3=Pstar1private[12]

local ug1=Pstar1group[1]
local ug2=Pstar1group[2]
local ug3=Pstar1group[3]
local lg1=Pstar1group[10]
local lg2=Pstar1group[11]
local lg3=Pstar1group[12]

twoway	(scatteri `bar' 0.65 `bar' 1.0, recast(line) lpattern(shortdash) lcolor(gray)) ///
		(scatteri `up1' 0.685 `lp1' 0.685,recast(line)  lc(red*1.5) lw(`lw')) ///
		(scatteri `ug1' 0.715 `lg1' 0.715,recast(line)  lc(black) lw(`lw')) ///
		(scatteri `up2' 0.765 `lp2' 0.765,recast(line) lc(red*1.5) lw(`lw')) ///
		(scatteri `ug2' 0.795 `lg2' 0.795,recast(line) lc(black) lw(`lw')) ///
		(scatteri `up3' 0.935 `lp3' 0.935,recast(line) lc(red*1.5) lw(`lw')) ///
		(scatteri `ug3' 0.965 `lg3' 0.965,recast(line) lc(black) lw(`lw')) ///
		(scatter Pstar1private tao_private if U==5, msymbol(circle) lw(vvthin) mfcolor(white) lcolor(black) mcolor(red*1.5)  mlw(`lw')) ///
		(scatter Pstar1private tao_private if U==15, msymbol(X) mfcolor(white) lcolor(black) mcolor(red*1.5)  mlw(`lw')) ///
		(scatter Pstar1private tao_private if U==25, msymbol(smtriangle) mfcolor(white) lcolor(black) mcolor(red*1.5) mlw(`lw') ) ///
		(scatter Pstar1private tao_private if U==50, msymbol(smdiamond) mfcolor(white) lcolor(black) mcolor(red*1.5)  mlw(`lw'))  ///
		(scatter Pstar1group tao_group if U==5, msymbol(circle) mfcolor(white)  lcolor(black) mcolor(black)  mlw(`lw')) ///
		(scatter Pstar1group tao_group if U==15, msymbol(X) mfcolor(white) lcolor(black) mcolor(black)  mlw(`lw')) ///
		(scatter Pstar1group tao_group if U==25, msymbol(smtriangle) mfcolor(white) lcolor(black) mcolor(black)  mlw(`lw')) ///
		(scatter Pstar1group tao_group if U==50, msymbol(smdiamond) mfcolor(white) lcolor(black) mcolor(black)  mlw(`lw')) ///
		(rbar  DIFF_1private10 DIFF_1private1 tao_private if U==25, barw(0.025) fcolor("255 204 204") lw(0) yaxis(2)) ///
		(scatter DIFF_1private1 tao_private if U==25, msymbol(P) mcolor("255 153 153") msize(tiny) yaxis(2) ) ///
		(scatter DIFF_1private5 tao_private if U==25, msymbol(P) mcolor("255 0 0") msize(tiny) yaxis(2)) ///
		(scatter DIFF_1private10 tao_private if U==25, msymbol(P) mcolor("153 0 0") msize(tiny) yaxis(2) ) ///
		(rbar  DIFF_1group10 DIFF_1group1 tao_group if U==25, barw(0.025) fcolor("240 240 240") lw(0)  yaxis(2)) ///
		(scatter DIFF_1group1 tao_group if U==25, msymbol(P) mcolor("224 224 224") msize(tiny) yaxis(2)) ///
		(scatter DIFF_1group5 tao_group if U==25, msymbol(P) mcolor("128 128 128") msize(tiny) yaxis(2)) ///		
		(scatter DIFF_1group10 tao_group if U==25, msymbol(P) mcolor("32 32 32") msize(tiny) yaxis(2)) ///
		, ///
		yline(0, lcolor(black)) ///
		xlabel(0.70 0.78 0.95, format(%03.2f) labsize("$tsize")) ///
		ylabel(none, axis(1)) ///
		ylabel(0 (25) 50, format(%2.0f) axis(2) labsize("$tsize")) ///
		yscale(range(-0.4 0.8) axis(1)) yscale(range(0 200) axis(2)) ///
		leg(off) ///
		title("{bf:B.East}",size("$psize") margin(b=3) position(11)) ///
		ytitle("{bf:Benefits ($/acre)                                      }",size("$tsize") margin(l=3) axis(2)) ///
		xtitle("{bf:Discount factor}",size("$tsize") margin(t=3)) ///
		saving($dt\Peast.gph, replace) ///
		fxsize(80)		
		
grc1leg2  $dt\Pwest.gph  $dt\Peast.gph  ,  rows(1) legendfrom( $dt\Pwest.gph) ///
	position(3) ring(1) saving($do\onefigure.gph, replace)  imargin(tiny) ///
	noauto iscale(1) 
graph dis,	xsize(4.75) ysize(2.5)

graph export $do\onefigure.tif, replace width(1425) //horizontal pixels=4.75inch*300dots per inch
graph export $do\onefigure.png, replace width(1425)


*------fig.S5-S7 equilibrium visualization figures------------------------------------*
cd $dt
loc TAOlist "70 78 95"
foreach TAO of local TAOlist{
grc1leg ///
	functions_tao`TAO'_u5_private.gph functions_tao`TAO'_u15_private.gph ///
	functions_tao`TAO'_u25_private.gph functions_tao`TAO'_u50_private.gph ///
	functions_tao`TAO'_u5_group.gph functions_tao`TAO'_u15_group.gph ///
	functions_tao`TAO'_u25_group.gph functions_tao`TAO'_u50_group.gph, ///
	row(2) legendfrom(functions_tao`TAO'_u5_private.gph) 
graph display, xsize(6) ysize(3.8)
graph export $do\equi_tao`TAO'.png, replace width(2000) 
}	

*------fig.S8-S9 cost-benefit full figures-------------------------------------*
*****cost-benefit figures
//east
use $dt\net_figure.dta, clear
keep if interest=="private"
local Ulist "5 15 25 50"
foreach U of local Ulist{
graph twoway (line Benefit_bar1 T if TAO==95&U==`U', lcolor(black%50) lpattern(shortdash)) (line Benefit_star1 T if TAO==95&U==`U', lcolor(red*1.5) lpattern(shortdash)) ///
	(line Benefit_bar1 T if TAO==70&U==`U', lcolor(black%50) lpattern(longdash)) (line Benefit_star1 T if TAO==70&U==`U', lcolor(red*1.5) lpattern(longdash)) ///
	(line Benefit_bar1 T if TAO==78&U==`U', lcolor(black%50) lpattern(solid)) (line Benefit_star1 T if TAO==78&U==`U', lcolor(red*1.5) lpattern(solid)), ///
	xlabel(1 5 10 )  ///
	ylabel(-10(20)50) ///
	title("{bf:u=`U'%}",size(medium) margin(b=3)) ///
		xtitle("") /// 
	ytitle("Benefit ($/acre)",margin(r=3))   ///
		xtitle("Number of years (T)", size(medium) margin(t=3)) /// 
		leg(region(lwidth(none) color(none)) on ring(0) position(11) col(2)  size(small) ///
		order(1 "Status quo: {&tau}=0.95" 2 "Ind. opt.: {&tau}=0.95"  5 "Status quo: {&tau}=0.78" ///
		6 "Ind. opt.: {&tau}=0.78" 3 "Status quo: {&tau}=0.70" 4 "Ind. opt.: {&tau}=0.70")) ///
	saving($dt\net_u`U'.gph, replace) ///
	fxsize(90) nodraw
}

cd $dt
grc1leg net_u5.gph net_u15.gph net_u25.gph net_u50.gph ,  rows(2) legendfrom(net_u5.gph) ring(1) position(6)
graph display, xsize(5) ysize(6) 
graph export $do\net_east.png, replace width(2000) 

//west
use $dt\net_figure.dta, clear
keep if interest=="private"
local Ulist "5 15 25 50"
foreach U of local Ulist{
graph twoway (line Benefit_bar0 T if TAO==95&U==`U', lcolor(black%50) lpattern(shortdash)) (line  Benefit_star0 T if TAO==95&U==`U', lcolor(red*1.5) lpattern(shortdash)) ///
	(line  Benefit_bar0 T if TAO==70&U==`U', lcolor(black%50) lpattern(longdash)) (line  Benefit_star0 T if TAO==70&U==`U', lcolor(red*1.5) lpattern(longdash)) ///
	(line  Benefit_bar0 T if TAO==78&U==`U', lcolor(black%50) lpattern(solid)) (line  Benefit_star0 T if TAO==78&U==`U', lcolor(red*1.5) lpattern(solid)), ///
	xlabel(1 5 10)  ///
	ylabel(0(250)750) ///
	title("{bf:u=`U'%}",size(medium) margin(b=3)) ///
		xtitle("") /// 
	ytitle("Benefit ($/acre)",margin(r=3))   ///
		xtitle("Number of years (T)", size(medium) margin(t=3)) /// 
		leg(region(lwidth(none) color(none)) on ring(0) position(11) col(2)  size(small) ///
		order(1 "Status quo: {&tau}=0.95" 2 "Ind. opt.: {&tau}=0.95"  5 "Status-quo: {&tau}=0.78" ///
		6 "Ind. opt.: {&tau}=0.78" 3 "Status quo: {&tau}=0.70" 4 "Ind. opt.: {&tau}=0.70")) ///
	saving($dt\net_u`U'.gph, replace) ///
	fxsize(90) nodraw
}

cd $dt
grc1leg net_u5.gph net_u15.gph net_u25.gph net_u50.gph ,  rows(2) legendfrom(net_u5.gph) ring(1) position(6)
graph display, xsize(5) ysize(6) 
graph export $do\net_west.png, replace width(2000) 




*------------------table.S6-S7 obtain table that corresponds to Fig.4---------*
//Pvalues
use $dt\net_figure.dta, clear
gen u=U/100
gen tao=TAO/100
keep if T==1
keep P* tao u interest
reshape wide P*, i(tao u) j(interest) string
order u tao Pbar0private Pstar0private Pstar0group Pbar1private Pstar1private Pstar1group
drop Pbar*group
cd $dt
asdoc list, save(Pvalues.doc) replace dec(2) tzok

//Benefits
use $dt\net_figure.dta, clear
keep TAO U T Benefit* DIFF* interest
keep if U==25
keep if T==1|T==5|T==10
gen tao=TAO/100
order  tao, before(T)
order interest, before(tao)
drop TAO U
replace interest="individual" if interest=="private"
cd $dt
asdoc list, save(Benefits.doc) replace dec(2) tzok


/*
*====================obtain ECDF==============================================*

*****ECDF
use $dd\estim.dta, clear
keep if bt==0
keep if soilinsecticide==0
keep if seedtrtrate==0
rename year year_trial
gen year=year_trial-1
keep if year>=2014 & year<=2016
keep year east rw

forv i=0/1{
preserve
//hist rw if east==`i'&bt==0

cumul rw if east==`i', gen(cdf`i') 

keep rw cdf`i' 
keep if cdf`i'!=.
sort cdf`i'
gen cdf`i'_per=cdf`i'*100
tostring cdf`i'_per, force replace
split cdf`i'_per, gen(part) parse(".")
sort cdf`i'_per
bysort part1 (cdf`i'): gen is_minimum = cdf`i' == cdf`i'[1] //quantile level at interval=1
keep if is_minimum==1
keep rw part1 cdf`i'
destring part1, force replace
drop if part1==.
sort part1
rename part1 percentile 
rename rw rw`i'
save $dt\temp`i'.dta, replace
restore
}


clear
set obs 100
gen percentile=_n
merge 1:1 percentile using $dt\temp0.dta, nogen
merge 1:1 percentile using $dt\temp1.dta, nogen
for var cdf*: replace X=percentile/100 //use the closest percentile so don't need to insert missing values for cdf0/cdf1.
save $dt\ecdf.dta, replace

*/