********************************************************************************
*				program: descriptives									   	   *
*				author: Ziwei Ye 											   *
********************************************************************************
global root "D:"
global data "$root\1_data"

global ds "$data\source"
global dd "$data\data"
global dp "$data\proc"
global dt "$data\temp"
global do "$data\out"

graph set window fontface "Helvetica"
global psize "10pt" //panel text size
global tsize "7pt" //figure test size

*------------------Fig.2 trend analysis----------------------------------------* 
//panel A
cd $do
use $dd\estim_use.dta, clear
keep if soilinsecticide==0 &seedtrtrate==0  
drop if stateabbr=="OH" | stateabbr=="MN"

preserve
collapse (mean) hist_scrd, by(year)
save $dt\bt_year.dta, replace
restore

collapse (mean) mean= rw  ///
	(sd) sd= rw ///
	(count) n=rw ///
	, by(bt year)
gen hi=mean+invttail(n-1,0.025)*(sd / sqrt(n))
gen lo=mean-invttail(n-1,0.025)*(sd / sqrt(n))

reshape wide mean sd n hi lo, i(year) j(bt)
gen mean_d=mean0-mean1
merge 1:1 year using $dt\bt_year.dta


graph twoway ///
(rarea hi1 lo1 year , fcolor("250 230 240") lw(0)) ///
 (connected mean1 year ,msymbol(i) lw(1) color("195 149 164")) ///
 (rarea hi0 lo0 year , fcolor("231 239 250") lw(0)) ///
 (connected mean0 year ,msymbol(i) lw(1) color("169 184 198")) ///
 (connected mean_d year,msymbol(O) lpattern(solid) lcolor(black) mcolor(black) mfcolor(white) lw(1) yaxis(1)) ///
(bar hist_scrd year,bcolor("green%15") lw(0) yaxis(1)) , ///
	xlabel(2004 (4) 2016,labsize("$tsize"))  ///
	ylabel(0 (0.2) 1, labsize("$tsize") axis(1) format(%03.1f))  ///
	legend(ring(0) position(11) region(lwidth(none)) size("$tsize") col(2) ///
	order( 2 "Bt injury" 4 "Non-Bt injury" 5 "Bt efficacy" 6 "Bt planting" "history")) ///
	xtitle("{bf:Trial year}",size("$tsize") margin(t=3)) ///
	ytitle("{bf:Bt planting history / Root injury}", axis(1) margin(r=3) size("$tsize")) ///	
	title("{bf:A}",size("$psize") ring(1) position(11)  margin(b=3)) ///
	saving($dt\bteffect_year.gph, replace) 
graph display, xsize(2.25) ysize(2.25)
graph save $dt\bteffect_year.gph, replace


//panel B
use $dd\estim_use.dta, clear
keep if soilinsecticide==0 &seedtrtrate==0  
drop if stateabbr=="OH" | stateabbr=="MN"

local interv=10/100
local begin=`interv'/2
egen rwcut = cut(hist_scrd), at(0 (`interv') 1)
collapse (mean) rw, by(bt rwcut)
reshape wide rw, i(rwcut) j(bt)
gen rw_d=rw0-rw1
gen xlab=rwcut+`interv'/2
keep if rwcut<0.75

graph twoway ///
	(lfitci rw0 xlab,  lcolor("169 184 198") clw(1) lpattern(shortdash) acolor("231 239 250") alwidth(0)  ) ///
	(lfitci rw1 xlab,  lcolor("195 149 164") clw(1) lpattern(shortdash) acolor("250 230 240") alwidth(0)) ///
	(lfitci rw_d xlab,  lcolor(black) clw(1) lpattern(shortdash)  acolor(grey%10) alwidth(0)  ) ///
	(scatter rw_d xlab,msymbol(O)  mcolor(black) lpattern(solid)  mfcolor(white) lw(1) msize(medium)) ///
	(scatter rw0 xlab, msymbol(X) mcolor("169 184 198") lpattern(solid) lw(2) msize(large) ) ///
	(scatter rw1 xlab,msymbol(X) mcolor("195 149 164") lpattern(solid)  lw(2) msize(large) ), ///
	xlabel(0 (0.2) 0.8,labsize("$tsize") format(%03.1f)) ///
	ylabel(0 (0.2) 1,labsize("$tsize")format(%03.1f) ) ///
	title("{bf:B}",size("$size") ring(1) position(11) margin(b=3)) /// 
	leg(region(lwidth(none)) on  ring(0) position(11) col(1)  size("$tsize") ///
	order(4 "Bt injury" 2 "non-Bt injury" 6 "Bt efficacy")) ///
	xtitle("{bf:Bt planting history}", size("$tsize") margin(t=3)) ///
	ytitle("{bf:Root injury}", margin(r=3) size("$tsize") )  ///
	saving($dt\bteffect_coverage_CC.gph, replace) 
graph display, xsize(2.25) ysize(2.25)
graph save $dt\bteffect_coverage_CC.gph, replace

cd $dt
graph combine bteffect_year.gph bteffect_coverage_CC.gph, iscale(1) ///
	col(2) rows(1) ///
	xsize(4.75) ysize(2.5)
graph export $do\bteffect.tif, replace width(1425)
graph export $do\bteffect.png, replace width(1425)



*-----------------Fig.1B RW trial sites map-------------------------------------------------------*
use $dd\estim_use.dta, clear
keep year location stateabbr _X _Y
bys location: gen nobs_avg=_N/12 

collapse (mean) nobs_avg, by(stateabbr location  _X _Y)

egen size=cut(nobs_avg), at(0,2,4,6,8,10,12,14,16)

geo2xy _Y _X,proj(albers) replace
save $dt\RWmap.dta, replace

*keep Midwest area only
cd "D:\Dropbox\0_Graduate Program\计量\stata\graph\examples" 
use usa_state_shp, clear //usa state shapefile, publicly accessible
merge m:1 _ID using usa_state  // to get the identifiers
destring STATEFP, replace 
keep if inlist(STATEFP,19,17,18,39,38,46,31,20,27,29,55,26)

*get the mainland map
drop _CX- _merge 
sort _ID shape_order
geo2xy _Y _X, proj(albers) replace
save MidWest_shp_clean.dta, replace //the basemap dataset to be used.

use usa_state_shp, clear
merge m:1 _ID using usa_state  // to get the identifiers
destring STATEFP, replace 
keep if inlist(STATEFP,19,17,18,39,38,46,31,20,27,29,55,26)

//twoway (scatter _Y _X) //get the map of the mainland states
use usa_state, clear //the master data must contain the basemap identifier
spmap using MidWest_shp_clean, id(_ID) fcolor(white) ///basemap defined as mainland shapefile in stata format
	point(data($dt\RWmap.dta) x(_X) y(_Y) fcolor(Greens) by(size) ///
	shape(triangle triangle triangle triangle triangle triangle triangle) ///
	ocolor(green*1.7 green*1.7 green*1.7 green*1.7 green*1.7 green*1.7 green*1.7) ///
	legenda(off))

graph dis, xsize(3) ysize(2.5)
graph export $do/RWmap_v6.png, replace width(900)
graph export $do/RWmap_v6.tif, replace width(900)
graph export $do/RWmap_v6.pdf, replace 

*--------------------fig.S1 state maps-----------------------------------------*
*****state level
//prepare background map data for cornbelt region
use "D:\Dropbox\1_Research\9_halo\1_data\source\maptile_geographies\state_database_clean.dta",clear 
//the dta. file can be obtained here: https://michaelstepner.com/maptile/geographies/
merge m:1 statefips using $dt/cornbeltstates.dta
keep if cornbelt==1
save $dt/state_database_clean_cornbelt.dta, replace

use "D:\Dropbox\1_Research\9_halo\1_data\source\maptile_geographies\state_coords_clean.dta",clear
//the dta. file can be obtained here: https://michaelstepner.com/maptile/geographies/
gen oid=_n
rename _ID _polygonid
merge m:1 _polygonid using $ds\maptile_geographies\state_database_clean.dta, keepusing(statefips) nogen
merge m:1 statefips using $dt\cornbeltstates.dta, nogen keepusing(cornbelt)
sort oid

keep if cornbelt==1
rename _polygonid _ID
drop statefips cornbelt
save $dt\state_coords_clean_cornbelt.dta, replace

use $dd\estim.dta, clear
keep if bt==0
keep if soilinsecticide==0 &seedtrtrate==0 
keep if year<=2016 & year>=2014

gen rw_scrd=rw_rate_scrd/100
gen rw_statefips=rw_rate_statefips/100
collapse (mean) rw rw_scrd rw_statefips, by(statefips)

capture drop _merge
merge 1:1 statefips using $dt\cornbeltstates.dta, keepusing(cornbelt) nogen

keep if cornbelt==1

merge 1:1 statefips using $dt\state_database_clean_cornbelt.dta, nogen
drop _merge


save $dt\rwxy.dta, replace

spmap rw using $ds\maptile_geographies\state_coords_clean.dta, id(_polygonid)  ///
	clmethod(custom) clbreaks(0 0.15 0.3 0.45 1) fcolor(Blues) ///
	legend(lab(1 "No data") lab(2 "0.00 - 0.15") lab(3 "0.15 - 0.30") lab(4 "0.30 - 0.45") lab(5 "0.45 - 1.00"))  ///
	legend(ring(1) position(3) size(medium))  saving($do\threeyear_rw, replace) ///
	title("{bf:A. Root injury}",size(medium) pos(11))

spmap rw_scrd using $ds\maptile_geographies\state_coords_clean.dta, id(_polygonid)  ///
	clmethod(custom) clbreaks(0 0.20 0.4 0.6 1) fcolor(Reds) ///
	legend(lab(1 "No data") lab(2 "0.0 - 0.2") lab(3 "0.2 - 0.4") lab(4 "0.4 - 0.6") lab(5 "0.6 - 1.00"))  ///
	legend(ring(1) position(3) size(medium))  saving($do\threeyear_rw_scrd, replace) ///
	title("{bf:B. Bt planting rate}",size(medium) pos(11))


graph combine $do\threeyear_rw.gph $do\threeyear_rw_scrd.gph, col(1)
graph display, xsize(4) ysize(5.5) scale(1.4)
graph export "$do/threeyear.png", wid(2000) replace  
graph export "$do/threeyear.tif", wid(2000) replace  


*----------------fig.S2 rotation-----------------------------------------*
use $dd\estim_use, clear
collapse (mean) cornshare_corn  share_corncorn share_corn, by(year stateabbr)
for var cornshare_corn  share_corncorn share_corn: replace X=X/100
graph twoway (bar cornshare_corn year,bcolor(gray%30) lw(0) yaxis(1)) ///
	(scatter share_corncorn year, msymbol(O) mcolor(black) mfcolor(white) msize(1)) ///
	(scatter share_corn year,msymbol(O) mcolor(black) msize(1)) ///
	,by(stateabbr, col(5)) ///
	xlabel(2004 (6) 2016, labsize(medium)) ylabel(0 (0.25) 0.50, labsize(medium) format(%03.1f) ) ///
	legend(region(lwidth(none)) size(small) col(1) ///
	lab(1 "Ratio of maize-to-maize acres to maize acres") ///
	lab(2 "Ratio of maize-to-maize acres to total acres") ///
	lab(3 "Ratio of maize acres to total acres")) ///
	xtitle("Year", size(small) margin(t=3)) ytitle("Ratio", size(small) margin(r=3)) 
//manually edit box color and so on
graph export "$do\cornshare.tif", wid(2000) replace  
graph export "$do\cornshare.png", wid(2000) replace  


*----------------fig.S3 root injury and Bt planting history by state---------------*

use $dd\estim_use.dta, clear
keep if bt==0
keep if soilinsecticide==0 &seedtrtrate==0  &year>=2005 
collapse (mean) rw hist_scrd, by(stateabbr year)

graph twoway (bar hist_scrd year, bcolor(green%30))(scatter rw year, yaxis(2)) ///
	(lfit rw year, yaxis(2)  ) ///
	, by(stateabbr, col(5)) ///
	legend(region(lwidth(none)) size(small) col(3) label(1 "Bt planting rate") label(2 "Root injury") label(3 "Fitted value") ) ///
	ylabel(0 (0.4) 0.8, labsize(medium) format(%03.1f) ) ylabel(0 (0.3) 0.6, labsize(medium) format(%03.1f) axis(2) )  ///
	xlabel(2004 (6) 2016, labsize(medium))  ///
	xtitle("Trial year", size(small)) ytitle("Bt planting rate",axis(1) size(small)) ytitle("Root injury",axis(2) size(small)) 
//manually edit box color and so on
graph export "$do\bystate.tif", wid(2000) replace  
graph export "$do\bystate.png", wid(2000) replace  
		
*----------------------------------fig.S4 time trends by trait------------------*
use $dd\estim_use.dta, clear
duplicates drop year scrd, force
collapse (sum) *acres_scrd, by(year)

foreach v of varlist only_* two_*{
	renvars `v',postdrop(11) display
}

foreach v of varlist only_* two_*{
	gen `v'_rate=`v'/acres_scrd
}

gen h1=only_Cry3Bb1_rate
loc t=1
foreach v of varlist only_Cry3435A_rate only_mCry3A_rate two_Cry3Bb1_Cry3435A_rate  two_mCry3A_Cry3435A_rate two_mCry3A_eCry31Ab_rate{
	loc s=`t'
	loc t=`t'+1
	gen h`t'=h`s'+`v'
	loc s=`s'+1
}

sort year

//red to blue pallette
loc r1 "219 049 036"
loc r2 "252 140 090"
loc r3 "255 223 146"
loc b3 "230 241 243"
loc b2 "144 190 224"
loc b1 "075 116 178"


loc barw=0.8
loc lw=0.2
graph twoway (bar h6 h5 h4 h3 h2 h1 year,barw(`barw' `barw' `barw' `barw' `barw' `barw') ///
	lw(`lw' `lw' `lw' `lw' `lw' `lw') lcolor(black black black black black black) ///
	color("`b1'" "`b2'" "`b3'" "`r3'"  "`r2'"  "`r1'")), ///
	yscale(range(0 0.8))  ///
	ylabel(0 (0.2) 0.8 ,format(%03.1f)  angle(horizontal)) ///
	xlabel(2004 (4) 2016) ///
	xtitle("Year",margin(t=3)) ytitle("Bt planting rate",margin(r=3)) ///
	saving($dt\rwstacking.gph, replace)   ///
	leg(region(lwidth(none) color(none)) ring(0) position(11) col(1)  size(medsmall) ///
	lab(1 "mCry3A+eCry3.1Ab") lab(2 "Cry34/35A+mCry3A") lab(3 "Cry3Bb1+Cry34/35A") ///
	lab(4 "mCry3A") lab(5 "Cry34/35A") lab(6 "Cry3Bb1")) 
	
graph dis,	xsize(6) ysize(5.5) scale(0.9)
graph export "$do\rwstacking.tif", wid(2000) replace  
graph export "$do\rwstacking.png", wid(2000) replace  

*------------------------fig.S10 below-ground trait stacking---------------------------*
use $dd\estim.dta, clear
keep year *_rate_nation
duplicates drop year, force
keep if year>=2000 & year<=2016

gen h1=rwonly_rate_nation/100
gen h2=h1+rwcbonly_rate_nation/100
gen h3=h2+rwhtonly_rate_nation/100
gen h4=h3+rwcbht_rate_nation/100


sort year

//red to blue pallette
loc r1 "219 049 036"
loc r2 "252 140 090"
loc r3 "255 223 146"
loc b3 "230 241 243"
loc b2 "144 190 224"
loc b1 "075 116 178"

loc barw=0.8
loc lw=0.2
graph twoway (bar   h4 h3 h2 h1 year,barw(`barw' `barw' `barw' `barw') lw(`lw' `lw' `lw' `lw') lc(black black black black) color("`b3'" "`b2'" "`b1'" "`r1'")), ///
	yscale(range(0 0.6))  ///
	ylabel(0 (0.2) 0.6 ,format(%03.1f)  angle(horizontal)) ///
	xlabel(2000 (4) 2016) ///
	xtitle("Year",margin(t=3)) ytitle("Planting rate",margin(r=3)) ///
	saving($dt\gestacking.gph, replace)   ///
	leg(region(lwidth(none) color(none)) ring(0) position(11) col(1)  size(medsmall) ///
	lab(4 "RW") lab(3 "RW+CB") lab(2 "RW+HT") lab(1 "RW+CB+HT"))

graph dis,	xsize(6) ysize(5.5) scale(0.9)
graph export "$do\gestacking.tif", wid(2000) replace  
graph export "$do\gestacking.png", wid(2000) replace  	



*-----------------table. S1 cross validation between TraitTrak and USDA NASS---------------------*
//obtain TraitTrak rates
import delimited "$ds\ISUTraitTrakData_NEW.csv", clear 
//The TraitTrak data from GfK Kynetec company; see "Data and materials availability" for data availability
keep if crop=="Corn"

gen bt=(strpos(seedtraitcollapsed, "RW")>0|strpos(seedtraitcollapsed, "CB")>0)
gen ht=(strpos(seedtraitcollapsed, "Gly Tol")>0|strpos(seedtraitcollapsed, "LL")>0|strpos(seedtraitcollapsed, "HT other")>0)
gen btonly=(bt==1&ht==0)
gen htonly=(bt==0&ht==1)
gen stacked=(bt==1&ht==1)

gl vlist "bt ht btonly htonly stacked "
for var $vlist :gen X_acres=acres*X

collapse (sum) *acres, by(year statefips state)

foreach v of global vlist{
	gen `v'_rate_TT=`v'_acres/acres //Bt rate in Kynetec TraitTrak data.
}

save $dt/bt_TraitTrak.dta, replace

//obtain NASS rates
import delimited "$ds/BiotechCropsAllTables2023.csv", clear 
//State-level biotech crop adoption data from NASS, USDA; see "Data and materials availability" for data availability
gen variety="stacked" if attribute=="Stacked gene varieties (percent of all corn planted)"
replace variety="btonly" if attribute=="Insect-resistant (Bt) only (percent of all corn planted) "
replace variety="htonly" if attribute=="Herbicide-tolerant (HT) only (percent of all corn planted) "
keep if variety!=""
destring value, force replace
gen rate=value/100
keep year state rate variety
reshape wide rate, i(year state) j(variety) string
gen ratebt=ratebtonly+ratestacked
gen rateht=ratehtonly+ratestacked

//merge
merge 1:1 year state using $dt/bt_TraitTrak.dta, keep(matched)

//correlation
corr ratebt bt_rate_TT
corr rateht ht_rate_TT



