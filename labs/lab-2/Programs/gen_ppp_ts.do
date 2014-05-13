clear mata
local cig="c i g"
local cigxm="c i g x m"
local pl_set="c i g x1 m1 x2 m2 x3 m3 x4 m4 x5 m5 x6 m6"
local p="c i g x1 m1 x2 m2 x3 m3 x4 m4 x5 m5 x6 m6 x m r gdp"

* Load data and generate price series
use Data/na_data, clear
fillin countrycode year
drop _fillin

if "$chn"=="icp" {
	drop if countrycode=="CH2"
	}
	else {
	drop if countrycode=="CHN"
	replace countrycode="CHN" if countrycode=="CH2" // Replace standard NA data by Maddison/Wu NA data
	}

foreach x of local cigxm {
	gen p_`x'=v_`x'/q_`x'
	}
gen p_gfcf=v_gfcf/q_gfcf
replace p_i=p_gfcf if p_i<0 // Three times when v_i is negative but q_i positive due to inventory changes
replace v_i=v_gfcf if countrycode=="SAU" & year==1973

keep countrycode year p_c p_i p_g p_x p_m v_c v_gdp v_i v_g v_x v_m xr xr2 pop

* Merge in the US trade prices by BEC (applied to all countries, later combined with country-specific price trends)
merge m:1 year using Data/us_trade_prices
drop if _merge==2
drop _merge

* Merge in the benchmark relative CIG price data
if "$chn"=="icp" {
	local filename="ppp_icp_CHN_ICPbh.dta"
	}
else {
	local filename="ppp_icp"
	}
merge m:1 countrycode year using "Data/`filename'", keepusing(pl_c_$index1 pl_i_$index1 pl_g_$index1)
foreach x of local cig {
	ren pl_`x'_$index1 pl_`x'
	}
drop _merge

* Merge in the benchmark relative trade price data
merge 1:1 countrycode year using Data/xm_bec_data_q
drop if _merge==2
drop _merge
	
* Merge in the total merchandise trade data
merge 1:1 countrycode year using Data/tot_xm
drop _merge

* Integrate OECD/Eurostat block for 1995-2010
merge 1:1 countrycode year using Data/ppp_oecd_es_ts
drop if _merge==2
drop _merge
foreach x of local cig {
	replace pl_`x'=pl_`x'_eso if pl_`x'_eso~=.
	}
drop pl_*_eso

* Keep only 2005 benchmark results if global bm is set to 2005
foreach x of local pl_set {
	replace pl_`x'=. if year~=$bm & $bm~=.
	}

* Initialize base country & year
gen base_country=countrycode=="USA"
gen base_year=(year==2005)

* Generate shares of GDP at current national prices and implicit quantities at constant national prices
gen v_da=v_c+v_i+v_g
gen v_t=v_c+v_i+v_g+v_x-v_m
foreach x of local cig {
	gen sh_`x'_da=v_`x'/v_da
	}
foreach x of local cigxm {
	gen sh_`x'=v_`x'/v_t
	gen q_`x'=v_`x'/p_`x'
	}

* Generate Fisher GDP at constant national prices (NA)
bysort countrycode (year): gen gp_gdp=log((((p_c[_n]*q_c[_n-1]+p_i[_n]*q_i[_n-1]+p_g[_n]*q_g[_n-1]+ ///
										       p_x[_n]*q_x[_n-1]-p_m[_n]*q_m[_n-1])/v_t[_n-1])* ///
								   (v_t[_n]/(p_c[_n-1]*q_c[_n]+p_i[_n-1]*q_i[_n]+p_g[_n-1]*q_g[_n]+ ///
										       p_x[_n-1]*q_x[_n]-p_m[_n-1]*q_m[_n])))^(0.5))
bysort countrycode (year): gen gq_gdp=log(v_t[_n]/v_t[_n-1])-gp_gdp
run index.do
gen price_base=1
index price_base gp_gdp base_year countrycode p_gdp

* Extrapolate total merchandise trade using total NA trade
replace tot_x=tot_x*xr
replace tot_m=tot_m*xr
encode countrycode, gen(c)
xtset c year
local final=r(tmaxs)
forval t=1950(1)`final' {
	replace tot_x=L.tot_x*v_x/L.v_x if tot_x==. & year==`t'
	replace tot_m=L.tot_m*v_m/L.v_m if tot_m==. & year==`t'
}	
forval t=`final'(-1)1950 {
	replace tot_x=F.tot_x*v_x/F.v_x if tot_x==. & year==`t'
	replace tot_m=F.tot_m*v_m/F.v_m if tot_m==. & year==`t'
	}

* Compare Comtrade totals with NA totals. If id2=1, Comtrade totals are used
gen id=0
replace id=1 if tot_x<v_x & tot_m<v_m
replace id=. if tot_x==.
by c: egen tot_id=total(id) if year>=1962
by c: egen nobs_na=count(v_c) if year>=1962
gen temp=tot_id/nobs_na
by c: egen av_id=mean(temp)
drop temp
gen id2=0
replace id2=1 if av_id>.75

* Check whether C+I+G+Comtrade X-M is positive; otherwise replace id2 with zero for those years
gen pos_check=v_c+v_i+v_g+tot_x-tot_m
replace id2=0 if pos_check<0

* Separate exports & imports into merchandise trade (v_x and v_m) and other trade balance (v_r) if id2==1
gen v_r=0
gen q_r=0
replace v_r=. if v_c==.
replace q_r=. if v_c==.
replace v_r=(v_x-tot_x)-(v_m-tot_m) if id2==1
replace q_r=((v_x/p_x)-(tot_x/p_x))-((v_m/p_m)-(tot_m/p_m)) if id2==1
gen p_r=v_r/q_r
replace p_r=1 if v_r==0

replace v_x=tot_x if id2==1
replace v_m=tot_m if id2==1
replace q_x=v_x/p_x
replace q_m=v_m/p_m

* Extrapolate BEC-1 trade using v_x and v_m
forval i=1(1)6 {
	replace x`i'=1e-6*x`i'*xr
	replace m`i'=1e-6*m`i'*xr
	by c: ipolate x`i' year, gen(x`i'2)
	replace x`i'=x`i'2
	by c: ipolate m`i' year, gen(m`i'2)
	replace m`i'=m`i'2
	drop x`i'2 m`i'2
	local loop=0
	* Move backwards
	forval t=2005(-1)1950 {
		replace x`i'=F.x`i'*v_x/F.v_x if x`i'==. & year==`t'
		replace m`i'=F.m`i'*v_m/F.v_m if m`i'==. & year==`t'
		}
	* Move forwards
	forval t=2005(1)`final' {
		replace x`i'=L.x`i'*v_x/L.v_x if x`i'==. & year==`t'
		replace m`i'=L.m`i'*v_m/L.v_m if m`i'==. & year==`t'
		}
}
* Constrain x1..x6 and m1...m6 to sum to v_x and v_m
gen sum_x=x1+x2+x3+x4+x5+x6
gen sum_m=m1+m2+m3+m4+m5+m6
forval i=1(1)6 {
	gen v_x`i'=x`i'*v_x/sum_x
	gen v_m`i'=m`i'*v_m/sum_m
	drop x`i' m`i'
	}

* Generate base country price indexes
foreach x of local p {
	gen temp=p_`x' if countrycode=="USA"
	bysort year: egen p_`x'_base=total(temp)
	drop temp
	}

* Adjust BEC price indexes for country-level inflation and compute implict quantities
forval i=1(1)6 {
	replace p_x`i'=p_x`i'*p_x/p_x_base
	replace p_m`i'=p_m`i'*p_m/p_m_base
	gen q_x`i'=v_x`i'/p_x`i'
	gen q_m`i'=v_m`i'/p_m`i'
	}

* Normalize PPPs to p_GDP(us) (or to 1 depending on $norm) and interpolate and extrapolate
xtset c year
foreach x of local pl_set {
	if $norm==0 {
		gen ppp_`x'=pl_`x'*xr2*p_`x'_base
		}
	else if $norm==1 {
		gen ppp_`x'=pl_`x'*xr2
	}
	gen adj`x'=ppp_`x'/p_`x'
	}
gen i_bm_cig=adjc~=.
gen i_bm_xm=adjx1~=.
foreach x of local pl_set {
	by c: ipolate adj`x' year, generate(adj`x'2)
	}
gen i_ip_cig=(adjc2~=.)-i_bm_cig
gen i_ip_xm=(adjx12~=.)-i_bm_xm
foreach x of local pl_set {
	forval t=2005(-1)1950 {
		replace adj`x'2=F.adj`x'2 if adj`x'2==. & year==`t'
		}
	forval t=1970(1)`final' {
		replace adj`x'2=L.adj`x'2 if adj`x'2==. & year==`t'
	}
}
gen i_ep_cig=(adjc2~=.)-i_ip_cig-i_bm_cig
gen i_ep_xm=(adjx12~=.)-i_ip_xm-i_bm_xm
drop if i_bm_cig==0&i_ip_cig==0&i_ep_cig==0
replace i_bm_cig=i_bm_cig_eso if i_bm_cig_eso~=.
replace i_ip_cig=i_ip_cig_eso if i_ip_cig_eso~=.

* Calculate PPPs for the full time series
foreach x of local pl_set {
	replace ppp_`x'=p_`x'*adj`x'2
	}

* Fill in relative prices and trade if missing for the entire period
forval i=1(1)6 {
	replace ppp_x`i'=xr2 if ppp_x`i'==. & ppp_c~=.
	replace ppp_m`i'=xr2 if ppp_m`i'==. & ppp_c~=.
	replace v_x`i'=0 if v_x`i'==. & ppp_c~=.
	replace v_m`i'=0 if v_m`i'==. & ppp_c~=.
	}

* Generate exchange rate indicator variable and replace market by market+estimated set
gen i_xr=0
replace i_xr=1 if xr~=xr2
drop xr
ren xr2 xr

* Divide values and PPPs by the exchange rate
foreach x of local pl_set {
	replace v_`x'=v_`x'/xr
	replace ppp_`x'=ppp_`x'/xr
	}

* Run a GK/GEKS for every year
su year
local tmax=r(max)
mata: p_gdp=J(0,1,.)
mata: p_cig=J(0,1,.)
mata: p_x=J(0,1,.)
mata: p_m=J(0,1,.)
mata: c_list=J(0,1,"")
mata: y_list=J(0,1,.)
forval t=1950(1)`tmax' {
	gen present=(ppp_c~=. & year==`t')
	mata: st_sview(c=0,.,"countrycode","present")
	mata: st_view(present=0,.,"present","present")
	mata: base=sum(strmatch(c,"USA"):*(1..colsum(present))')
	mata: st_view(p=0,.,("ppp_c","ppp_i","ppp_g","ppp_x1","ppp_x2","ppp_x3","ppp_x4","ppp_x5","ppp_x6","ppp_m1","ppp_m2","ppp_m3","ppp_m4","ppp_m5","ppp_m6"),"present")
	mata: st_view(v=0,.,("v_c","v_i","v_g","v_x1","v_x2","v_x3","v_x4","v_x5","v_x6","v_m1","v_m2","v_m3","v_m4","v_m5","v_m6"),"present")
	mata: v=v:+1e-6
	mata: v=(v[,1..9],-v[,10..15])'
	mata: p=p'
	mata: p_y_gdp=$index2(v,p,base)
	if "$index2"=="geks" {
		mata: if (missing(p_y_gdp)>0) p_y_gdp=geks2(v,p,base) ;;
		mata: p_y_cig=$index2(v[1..3,.],p[1..3,.],base)
		mata: p_y_x=$index2(v[4..9,.],p[4..9,.],base)
		mata: p_y_m=$index2(-v[10..15,.],p[10..15,.],base)
	}
	mata: p_gdp=(p_gdp\p_y_gdp')
	if "$index2"=="geks" {
		mata: p_cig=(p_cig\p_y_cig')
		mata: p_x=(p_x\p_y_x')
		mata: p_m=(p_m\p_y_m')
	}
	mata: c_list=(c_list\c)
	mata: y_list=(y_list\J(colsum(present),1,`t'))
	drop present
}
preserve
	clear
	mata
		st_addvar("str10","countrycode")
		st_addvar("double",("year", "ppp_gdp_o"))
		st_addobs(rows(c_list))
		st_sstore(.,"countrycode",c_list)
		st_store(.,("year","ppp_gdp_o"),(y_list,p_gdp))
	end
	if "$index2"=="geks" {
		mata: st_addvar("double",("ppp_cig","ppp_x","ppp_m"))
		mata: st_store(.,("ppp_cig","ppp_x","ppp_m"),(p_cig,p_x,p_m))
	}
	save temp_p, replace
restore
merge 1:1 countrycode year using temp_p
drop _merge
erase temp_p.dta

* Compute price levels (USA GDP 2005=1), values in (XR-converted) dollars, implicit quantities, and reference prices
if $norm==0 {
	gen pl_gdpo=ppp_gdp_o*p_gdp_base
	}
else if $norm==1 {
	gen pl_gdpo=ppp_gdp_o
	}

foreach x of local pl_set {
	gen q_`x'_pi=v_`x'/ppp_`x'
	bysort year: egen calc_pi_`x'1=total(v_`x'/pl_gdpo)
	bysort year: egen calc_pi_`x'2=total(q_`x'_pi)
	gen pi_`x'= calc_pi_`x'1/calc_pi_`x'2
	replace pl_`x'=v_`x'/(pi_`x'*q_`x'_pi)
	drop calc_pi_`x'1 calc_pi_`x'2
	}
gen v_cig=v_c+v_i+v_g
replace v_gdp=v_gdp/xr
replace v_x=v_x/xr
replace v_m=v_m/xr
replace v_r=v_gdp-(v_c+v_i+v_g+v_x-v_m) // v_r is now the sum of the residual trade balance and any statistical discrepancy

if "$index2"=="geks" {
	gen pl_cig=ppp_cig
	gen pl_x=ppp_x
	gen pl_m=ppp_m
}
else {
	gen pl_cig=v_cig/(pi_c*q_c_pi+pi_i*q_i_pi+pi_g*q_g_pi)
	gen pl_x=v_x/(pi_x1*q_x1_pi+pi_x2*q_x2_pi+pi_x3*q_x3_pi+pi_x4*q_x4_pi+pi_x5*q_x5_pi+pi_x6*q_x6_pi)
	gen pl_m=v_m/(pi_m1*q_m1_pi+pi_m2*q_m2_pi+pi_m3*q_m3_pi+pi_m4*q_m4_pi+pi_m5*q_m5_pi+pi_m6*q_m6_pi)
}
sort countrycode year
