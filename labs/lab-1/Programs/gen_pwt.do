clear all
* Set parameters 
global index1="geks" 	// Determines aggregation method for aggregating basic headings, set to geks or gk
global index2="gk" 		// Determines aggregation method for aggregating main expenditure categories, set to geks or gk
global bm=. 			// Set to . to use all benchmark data; set to 2005 to keep only 2005 benchmark results (mimicking PWT7); alternative values: 1970, 1975, 1980, 1985 and 1996
global chn="pwt" 		// Set to "pwt" for adjusted basic headings and NA time series; set to "icp" for original basic headings and NA time series
global norm=0 			// Set to 0 to normalise to USA=GDP deflator; set to 1 to normalise to USA=1 in every year
* For original PWT8.0, set parameters to: "geks" "gk" . "pwt" 0

* Run preparatory programs
** Run PPP Mata programs (needs to run only once)
 run "Mata functions/laspeyres_func_create.do"
 run "Mata functions/paasche_func_create.do"
 run "Mata functions/fischer_func_create.do"
 run "Mata functions/geks_func_create.do"
 run "Mata functions/gk_func_create.do"

** Compute CIG prices for all benchmarks [not feasible unless ICP 2005 basic heading data is available; only a recalculation of older benchmarks is possible]
* run gen_bm_level
** Load indexing program into memory
 run index
** Generate capital, labor, schooling and labor share data
 run process_barrolee
** Generate PPP time series based on benchmarks and National Accounts
 run gen_oecd_es_ts
 run gen_ppp_ts

* Define expenditure categories
local exp="c i g x1 m1 x2 m2 x3 m3 x4 m4 x5 m5 x6 m6"
replace pop=1e-6*pop

* Compute implicit trade balance quantities for GDPe and GDPo computation
ren pl_cig pl_gdpe
gen q_x_e_pi=v_x/pl_gdpe
gen q_m_e_pi=v_m/pl_gdpe
gen q_r_e_pi=v_r/pl_gdpe
gen q_r_pi=v_r/pl_gdpo

* Calculate CGDPe, CGDPo and expenditure shares
gen cgdpe=(pi_c*q_c_pi)+(pi_i*q_i_pi)+(pi_g*q_g_pi)+q_x_e_pi-q_m_e_pi+q_r_e_pi
gen cgdpo=(pi_c*q_c_pi)+(pi_i*q_i_pi)+(pi_g*q_g_pi)+ ///
		  (pi_x1*q_x1_pi)+(pi_x2*q_x2_pi)+(pi_x3*q_x3_pi)+(pi_x4*q_x4_pi)+ ///
		  (pi_x5*q_x5_pi)+(pi_x6*q_x6_pi)-(pi_m1*q_m1_pi)-(pi_m2*q_m2_pi)- ///
		  (pi_m3*q_m3_pi)-(pi_m4*q_m4_pi)-(pi_m5*q_m5_pi)-(pi_m6*q_m6_pi)+q_r_pi
foreach x of local exp {
	gen csh_`x'=pi_`x'*q_`x'_pi/cgdpo
	}
forval i=1(1)6 {
	replace csh_m`i'=-csh_m`i'
	}
gen csh_x=csh_x1+csh_x2+csh_x3+csh_x4+csh_x5+csh_x6
gen csh_m=csh_m1+csh_m2+csh_m3+csh_m4+csh_m5+csh_m6
gen csh_r=q_r_pi/cgdpo
if "$index2"=="geks" {
	replace cgdpe=v_gdp/pl_gdpe
	replace cgdpo=v_gdp/pl_gdpo
	}

* Calculate RGDPe
bysort countrycode (year): gen cgdpe_ly=(pi_c[_n-1]*q_c_pi)+(pi_i[_n-1]*q_i_pi) ///
										+(pi_g[_n-1]*q_g_pi)+(q_x_e_pi-q_m_e_pi+q_r_e_pi)
bysort countrycode (year): gen cgdpe_ny=(pi_c[_n+1]*q_c_pi)+(pi_i[_n+1]*q_i_pi) ///
										+(pi_g[_n+1]*q_g_pi)+(q_x_e_pi-q_m_e_pi+q_r_e_pi)
bysort countrycode (year): gen grgdpe=log(((cgdpe_ly[_n]/cgdpe[_n-1])*(cgdpe[_n]/cgdpe_ny[_n-1]))^(0.5))
index cgdpe grgdpe base_year countrycode rgdpe

* Calculate RGDPo
bysort countrycode (year): gen cgdpo_ly=(pi_c[_n-1]*q_c_pi)+(pi_i[_n-1]*q_i_pi) ///
										+(pi_g[_n-1]*q_g_pi)+q_r_pi ///
										+(pi_x1[_n-1]*q_x1_pi)+(pi_x2[_n-1]*q_x2_pi) ///
										+(pi_x3[_n-1]*q_x3_pi)+(pi_x4[_n-1]*q_x4_pi) ///
										+(pi_x5[_n-1]*q_x5_pi)+(pi_x6[_n-1]*q_x6_pi) ///
										-(pi_m1[_n-1]*q_m1_pi)-(pi_m2[_n-1]*q_m2_pi) ///
										-(pi_m3[_n-1]*q_m3_pi)-(pi_m4[_n-1]*q_m4_pi) ///
										-(pi_m5[_n-1]*q_m5_pi)-(pi_m6[_n-1]*q_m6_pi)
bysort countrycode (year): gen cgdpo_ny=(pi_c[_n+1]*q_c_pi)+(pi_i[_n+1]*q_i_pi) ///
										  +(pi_g[_n+1]*q_g_pi)+q_r_pi ///
										  +(pi_x1[_n+1]*q_x1_pi)+(pi_x2[_n+1]*q_x2_pi) ///
										  +(pi_x3[_n+1]*q_x3_pi)+(pi_x4[_n+1]*q_x4_pi) ///
										  +(pi_x5[_n+1]*q_x5_pi)+(pi_x6[_n+1]*q_x6_pi) ///
										  -(pi_m1[_n+1]*q_m1_pi)-(pi_m2[_n+1]*q_m2_pi) ///
										  -(pi_m3[_n+1]*q_m3_pi)-(pi_m4[_n+1]*q_m4_pi) ///
										  -(pi_m5[_n+1]*q_m5_pi)-(pi_m6[_n+1]*q_m6_pi)
bysort countrycode (year): gen grgdpo=log(((cgdpo_ly[_n]/cgdpo[_n-1])*(cgdpo[_n]/cgdpo_ny[_n-1]))^(0.5))
index cgdpo grgdpo base_year countrycode rgdpo

* Get RGDPna
merge 1:1 countrycode year using Data/na_data, keepus(q_gdp)
keep if _merge==3
drop _merge
ren q_gdp rgdpna
gen temp1=rgdpna if year==2005
gen temp2=cgdpo if year==2005
bysort countrycode: egen rgdpna_05=total(temp1)
bysort countrycode: egen cgdpo_05=total(temp2)
replace rgdpna=rgdpna*cgdpo_05/rgdpna_05

* Get the employment data
merge 1:1 countrycode year using Data/pwt_emp_data
keep if _merge~=2
drop _merge
merge m:1 year using Data/basic_ch2_data, keepus(emp_wu)
replace emp=emp_wu if countrycode=="CHN"
drop emp_wu _merge
replace emp=1e-3*emp

* Compute CTFP
merge 1:1 countrycode year using Data/capital_results
drop _merge
if $norm==1 {
	gen temp=pl_k if countrycode=="USA"
	bysort year: egen pl_k_usa=total(temp)
	replace pl_k=pl_k/pl_k_usa
	replace ck=ck*pl_k_usa
	drop temp pl_k_usa
	}

merge 1:1 countrycode year using Data/lab_share_data, keepus(labsh)
drop _merge

merge 1:1 countrycode year using Data/bl_data
drop if _merge==2
drop _merge

drop c
encode countrycode, gen(c)
xtset c year
by c: ipolate yr_sch year, gen(s)

* Piece-wise linear return to education from Psacharopoulos
gen hc=exp(.134*min(s,4)+.101*max(0,min(4,(s-4)))+.068*max(0,s-8))
replace hc=. if s==.
replace hc=L.hc if year==2011
replace emp=1e-3*emp

gen yo=cgdpo/pop
gen l=emp/pop
gen k=ck/pop

local var="ye yo l k hc labsh"
foreach x of local var {
	gen temp=`x' if countrycode=="USA"
	bysort year: egen `x'_us=total(temp)
	drop temp
	}
gen ryo=yo/yo_us
gen avlsh=.5*(labsh+labsh_us)
gen rq=exp(avlsh*log(l/l_us)+avlsh*log(hc/hc_us)+(1-avlsh)*log(k/k_us))
gen ctfp=ryo/rq

* Compute RTFPna
merge 1:1 countrycode year using Data/rkna_data
drop _merge

xtset c year
gen avlsh2=.5*(labsh+L.labsh)
gen grtfpna=log(rgdpna/L.rgdpna)-avlsh2*log(emp/L.emp)-avlsh2*log(hc/L.hc)-(1-avlsh2)*log(rkna/L.rkna)
index price_base grtfpna base_year countrycode rtfpna 
replace rtfpna=. if hc==. | emp==. | rkna==. | labsh==.
replace hc=. if xr==.

* Construct identifyer variable for CIG
gen i_cig=0
replace i_cig=1 if i_bm_cig==1
replace i_cig=2 if i_ip_cig==1
replace i_cig=. if xr==.

lab def i 0 "Extrapolated" 1 "Benchmark" 2 "Interpolated"
lab val i_cig i

* Construct identifyer variable for XM
gen i_xm=0
replace i_xm=1 if i_bm_xm==1
replace i_xm=2 if i_ip_xm==1
replace i_xm=. if xr==.
lab val i_xm i

* Label and clean up identifyer variable for XR
replace i_xr=. if xr==.
lab def x 0 "Market" 1 "Estimated"
lab val i_xr x

* Mark outliers, discussed in "Outliers in PWT80"
gen i_outlier=0
replace i_outlier=1 if countrycode=="BRN" & year>=1970 & year<=1985
replace i_outlier=1 if countrycode=="BDI" & year>=1960 & year<=1999
replace i_outlier=1 if countrycode=="BMU" & year>=1998 & year<=2004
replace i_outlier=1 if countrycode=="GNB" & year>=1960 & year<=1995
replace i_outlier=1 if countrycode=="MOZ" & year>=1960 & year<=1991
replace i_outlier=1 if countrycode=="SLV" & year>=1986
replace i_outlier=. if xr==.
lab def y 0 "Regular" 1 "Outlier"
lab val i_outlier y

* Get expenditure share correlation and statistical capacity indicator
merge 1:1 countrycode year using Data/esh_cor_data
drop _merge
merge 1:1 countrycode year using Data/statcap_data
drop if _merge==2
drop _merge

* Drop GDPe and GDPo series if there are no national accounts for those years
replace rgdpe=. if cgdpo==.
replace rgdpo=. if cgdpo==.
replace rgdpna=. if cgdpo==.
replace rtfpna=. if cgdpo==.

* Change format to get full precision for export
format %14.3f cgdpo
format %14.3f rgdpo
format %14.3f cgdpe
format %14.3f rgdpe
format %14.3f rgdpna

* Get the country name and currency unit
merge m:1 countrycode using Data/names_units
keep if _merge==3
drop _merge

* Generate final format, including labeling and variable selection
run gen_final_pwt.do
