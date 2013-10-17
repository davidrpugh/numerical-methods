version 12
* Compute price levels for CIG in every benchmark 
* Note that basic heading data for ICP 2005 cannot be provided. 
* However if available through World Bank channels, the program can easily be ammended

local bms="1970 1975 1980 1985 1996"
foreach x of local bms {
	clear mata
	import excel using Data/icp`x'.xls, sheet("StataReady") firstrow clear
	run level_analysis
	replace year=`x'
	save ppp_icp_`x', replace
	}

* Append earlier years
append using ppp_icp_1985
append using ppp_icp_1980
append using ppp_icp_1975
append using ppp_icp_1970

foreach x of local bms {
	erase ppp_icp_`x'.dta
	}

* Order variables and save results
order countrycode year pl_*_geks pl_*_gk
save ppp_icp_ex2005, replace
