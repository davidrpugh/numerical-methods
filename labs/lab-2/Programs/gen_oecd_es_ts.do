version 12

use Data/ppp_oecd, clear
drop if countrycode=="DEU"

append using Data/ppp_es

merge 1:1 countrycode year using Data/na_data
keep if year>=1995
bysort countrycode: egen check=total(pl_c_geks)
drop if check==0
drop check _merge

local var="c i g gdp"
foreach x of local var {
	gen p_`x'=v_`x'/q_`x'
	}
ren (pl_c_$index1 pl_i_$index1 pl_g_$index1) (pl_c pl_i pl_g)
keep countrycode year pl_c pl_i pl_g xr p_c p_i p_g p_gdp

* Change XR to Germany base
gen temp=xr if countrycode=="DEU"
bysort year: egen xr_base=total(temp)
replace xr=xr/xr_base
drop temp xr_base

* Generate base country price indexes and normalize PPPs to p_.._base
local exp="c i g"
foreach x of local exp {
	gen temp=p_`x' if countrycode=="DEU"
	bysort year: egen p_`x'_base=total(temp)
	drop temp
	gen ppp_`x'=pl_`x'*xr // *p_`x'_base
	gen adj`x'=ppp_`x'/p_`x'
	}
gen i_bm_cig=adjc~=.

* Interpolate and extrapolate
foreach x of local exp {
	bysort countrycode: ipolate adj`x' year, generate(adj`x'2)
	}
gen i_ip_cig=(adjc2~=.)-i_bm_cig
foreach x of local exp {
	local loop=0
	while `loop' < (2005-1995) {
		bysort countrycode (year): replace adj`x'2=adj`x'2[_n+1] if adj`x'2==. // Move backwards		
		bysort countrycode (year): replace adj`x'2=adj`x'2[_n-1] if adj`x'2==. // Move forwards
		local loop=`loop'+1
		}
	}
gen i_ep_cig=(adjc2~=.)-i_ip_cig-i_bm_cig

* Calculate PPPs for the full time series
foreach x of local exp {
	replace ppp_`x'=p_`x'*adj`x'2
	}

* Change to base USA
drop p_c_base p_i_base p_g_base
local p="xr ppp_c ppp_i ppp_g p_c p_i p_g"
foreach x of local p {
	gen temp=`x' if countrycode=="USA"
	bysort year: egen `x'_base=total(temp)
	replace `x'=`x'/`x'_base
	drop temp
	}
*
gen pl_c_eso=ppp_c/xr
gen pl_i_eso=ppp_i/xr
gen pl_g_eso=ppp_g/xr
ren (i_bm_cig i_ip_cig) (i_bm_cig_eso i_ip_cig_eso)
keep countrycode year pl_c_eso pl_i_eso pl_g_eso i_bm_cig_eso i_ip_cig_eso

save Data/ppp_oecd_es_ts, replace
