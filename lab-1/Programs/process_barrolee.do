use Data/BL2013_MF_v1.3.dta, clear
replace WBcode="MDA" if WBcode=="ROM"
replace WBcode="SRB" if WBcode=="SER"
ren WBcode countrycode

keep countrycode year agefrom ageto yr_sch pop
order countrycode year agefrom ageto yr_sch pop
gen tyr=yr_sch*pop
drop if agefrom==15 & ageto==999
drop if agefrom==25 & ageto==999
collapse (sum) pop tyr, by(countrycode year)
gen yr_sch_lf=tyr/pop
drop tyr pop
save temp, replace

*
use Data/BL2013_MF_v1.3.dta, clear
replace WBcode="MDA" if WBcode=="ROM"
replace WBcode="SRB" if WBcode=="SER"
ren WBcode countrycode

keep if ageto==999
keep if agefrom==25 | agefrom==15
gen a15=(agefrom==15)
keep countrycode year yr_sch a15
reshape wide yr_sch, i(countrycode year) j(a15)
ren (yr_sch0 yr_sch1) (yr_sch25up yr_sch15up)
merge 1:1 countrycode year using temp
erase temp.dta

keep countrycode year yr_sch15up
ren yr_sch15up yr_sch
label var yr_sch "Average years of schooling of the 15+ population"
save Data/bl_data, replace
