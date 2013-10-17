* Identifier variables
label var country "Country name"
label var currency_unit "Currency unit"

* Real GDP, employment and population levels
label var rgdpe "Expenditure-side real GDP at chained PPPs (in mil. 2005US$)"
label var rgdpo "Output-side real GDP at chained PPPs (in mil. 2005US$)"
label var pop 	"Population (in millions)"
label var emp	"Number of persons engaged (in millions)"
label var avh	"Average annual hours worked by persons engaged"
label var hc	"Human capital index, based on years of schooling (Barro/Lee, 2010) and returns to education (Psacharopoulos, 1994)"

* Current price GDP, capital and TFP
label var cgdpe "Expenditure-side real GDP at current PPPs (in mil. 2005US$)"
label var cgdpo "Output-side real GDP at current PPPs (in mil. 2005US$)"
label var ck "Capital stock at current PPPs (in mil. 2005US$)"
label var ctfp "TFP level at current PPPs (USA=1)"

* National accounts-based variables
label var rgdpna "Real GDP at constant 2005 national prices (in mil. 2005US$)"
label var rkna "Capital stock at constant 2005 national prices (in mil. 2005US$)"
label var rtfpna "TFP at constant national prices (2005=1)"
label var labsh "Share of labour compensation in GDP at current national prices"

* Exchange rates and GDP price levels
label var xr "Exchange rate, national currency/USD (market+estimated)"label var pl_gdpe "Price level of CGDPe (PPP/XR), price level of USA GDPo in 2005=1"label var pl_gdpo "Price level of CGDPo (PPP/XR),  price level of USA GDPo in 2005=1"

* Data information variables
label var i_cig "0/1/2: relative price data for consumption, investment and government is extrapolated (0), benchmark (1), or interpolated (2)"
label var i_xm "0/1/2: relative price data for exports and imports is extrapolated (0), benchmark (1), or interpolated (2)"
label var i_xr "0/1: the exchange rate is market-based (0) or estimated (1)"
label var i_outlier "0/1: the observation on pl_gdpe or pl_gdpo is not an outlier (0) or an outlier (1)"
label var cor_exp "Correlation between expenditure shares of the country and the US (benchmark observations only)"
label var statcap "Statistical capacity indicator (source: World Bank, developing countries only)"

* Shares in CGDPo
label var csh_c "Share of household consumption at current PPPs"
label var csh_i "Share of gross capital formation at current PPPs"
label var csh_g "Share of government consumption at current PPPs"
label var csh_x "Share of merchandise exports at current PPPs"
label var csh_m "Share of merchandise imports at current PPPs"
label var csh_r "Share of residual trade and GDP statistical discrepancy at current PPPs"

* Price levels, expenditure categories and capital
label var pl_c "Price level of household consumption,  price level of USA GDPo in 2005=1"
label var pl_i "Price level of capital formation, price level of USA GDPo in 2005=1"
label var pl_g "Price level of government consumption, price level of USA GDPo in 2005=1"
label var pl_x "Price level of exports, price level of USA GDPo in 2005=1"
label var pl_m "Price level of imports, price level of USA GDPo in 2005=1"
label var pl_k "Price level of the capital stock, price level of USA 2005=1"

* Detailed trade series
label var csh_x1 "Share of food and beverages exports at current PPPs"
label var csh_x2 "Share of industrial supplies exports at current PPPs"
label var csh_x3 "Share of fuels and lubricants exports at current PPPs"
label var csh_x4 "Share of capital goods exports at current PPPs"
label var csh_x5 "Share of transport equipment exports at current PPPs"
label var csh_x6 "Share of consumer goods exports at current PPPs"
label var csh_m1 "Share of food and beverages imports at current PPPs"
label var csh_m2 "Share of industrial supplies imports at current PPPs"
label var csh_m3 "Share of fuels and lubricants imports at current PPPs"
label var csh_m4 "Share of capital goods imports at current PPPs"
label var csh_m5 "Share of transport equipment imports at current PPPs"
label var csh_m6 "Share of consumer goods imports at current PPPs"

label var pl_x1 "Price level of food and beverages exports, USA GDPo in 2005=1"
label var pl_x2 "Price level of industrial supplies exports, USA GDPo in 2005=1"
label var pl_x3 "Price level of fuels and lubricants exports, USA GDPo in 2005=1"
label var pl_x4 "Price level of capltal goods exports, USA GDPo in 2005=1"
label var pl_x5 "Price level of transport equipment exports, USA GDPo in 2005=1"
label var pl_x6 "Price level of consumer goods exports, USA GDPo in 2005=1"
label var pl_m1 "Price level of food and beverages imports, USA GDPo in 2005=1"
label var pl_m2 "Price level of industrial supplies imports, USA GDPo in 2005=1"
label var pl_m3 "Price level of fuels and lubricants imports, USA GDPo in 2005=1"
label var pl_m4 "Price level of capltal goods imports, USA GDPo in 2005=1"
label var pl_m5 "Price level of transport equipment imports, USA GDPo in 2005=1"
label var pl_m6 "Price level of consumer goods imports, USA GDPo in 2005=1"

keep countrycode country currency_unit year rgdpe rgdpo pop emp avh hc cgdpe cgdpo ck ctfp ///
	rgdpna rkna rtfpna labsh xr pl_gdpe pl_gdpo i_cig i_xm i_xr i_outlier cor_exp ///
	statcap csh_c csh_i csh_g csh_x csh_m csh_r pl_c pl_i pl_g pl_x pl_m pl_k
	
order countrycode country currency_unit year rgdpe rgdpo pop emp avh hc cgdpe cgdpo ck ctfp ///
	rgdpna rkna rtfpna labsh xr pl_gdpe pl_gdpo i_cig i_xm i_xr i_outlier cor_exp ///
	statcap csh_c csh_i csh_g csh_x csh_m csh_r pl_c pl_i pl_g pl_x pl_m pl_k
  
sort countrycode year
compress

* Drop variables that make no sense in the GEKS data
if "$index2"=="geks" {
	drop csh* r*
	}

save pwt_output, replace
