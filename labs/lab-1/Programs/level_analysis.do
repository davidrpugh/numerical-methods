version 12
mata
	k=4								// Number of variables at the top that are not v or p
	c=st_nvar()-k					// Number of countries
end	

// Generate indicator variable
gen value=(var=="v")				// where the value data can be found
gen price=(var=="p")				// where the price data can be found
gen xrate=(var=="xr")				// where the exchange rate data can be found
gen pop=(var=="pop")				// where the population data can be found
gen con=(cig=="c")					// consumption categories
gen inv=(cig=="i")					// investment categories
gen gov=(cig=="g")					// government categories
gen tb=(cig=="x")					// trade balance
gen save_var=0

// Run the analysis
mata
// Preparation
	range_c=k+1..k+c				// Give the range where the country data can be found (variable numbers x..y) 
	base=st_varindex("USA")-k		// Define the base country

	st_view(p=0,.,range_c,"price") 	// View on the price data
	st_view(v=0,.,range_c,"value") 	// View on the value data
	st_view(xr=0,.,range_c,"xrate")	// View on exchange rate data
	st_view(pop=0,.,range_c,"pop")	// View on population data

	range_bh=(1..rows(p))'						// Range of basic headings
	st_view(con=0,range_bh,st_varindex("con"))	// con identifies consumption headings
	st_view(inv=0,range_bh,st_varindex("inv"))	// inv identifies investment headings
	st_view(gov=0,range_bh,st_varindex("gov"))	// gov identifies government headings
	st_view(tb=0,range_bh,st_varindex("tb"))	// tb identifies inventory & trade balance headings

	p=p:/xr
	v=v:/xr
	
	p_c=select(p,con)
	p_i=select(p,inv)
	p_g=select(p,gov)

	v_c=select(v,con)
	v_i=select(v,inv)
	v_g=select(v,gov)

	cig=J(rows(p),1,1)-tb
	v_ex=select(v,cig)
	p_ex=select(p,cig)

	// GEKS
	pl_c_geks=geks(v_c,p_c,base)					// Estimate GEKS for C
		st_addobs(1) 								// Add a new observation for results
		st_store(st_nobs(),range_c,pl_c_geks)		// Write results
		stata(`"replace var="pl_c_geks" in -1"')	// Name aggregate variable
		stata(`"replace save_var=1 in -1"')			// Include this to save the resulting variable
	pl_i_geks=geks(v_i,p_i,base)					// Estimate GEKS for I
		st_addobs(1)
		st_store(st_nobs(),range_c,pl_i_geks)		
		stata(`"replace var="pl_i_geks" in -1"')	
		stata(`"replace save_var=1 in -1"')			
	pl_g_geks=geks(v_g,p_g,base)					// Estimate GEKS for G
		st_addobs(1) 								
		st_store(st_nobs(),range_c,pl_g_geks)		
		stata(`"replace var="pl_g_geks" in -1"')	
		stata(`"replace save_var=1 in -1"')			
	pl_gdp_geks=geks(v_ex,p_ex,base)				// Estimate GEKS for CIG
		st_addobs(1) 								
		st_store(st_nobs(),range_c,pl_gdp_geks)	
		stata(`"replace var="pl_geks" in -1"')		

	// GK
	pl_gdp_gk=gk(v_ex,p_ex,base)					// Estimate GK for CIG
		st_addobs(1) 								
		st_store(st_nobs(),range_c,pl_gdp_gk)		
		stata(`"replace var="pl_gk" in -1"')		
	ref_p=rowsum(v_ex:/pl_gdp_gk):/rowsum(v_ex:/p_ex)	// Calculate the reference prices corresponding to the PPPs
	q=v:/p											// Notional quantities
	pl_c_gk=colsum(v_c):/colsum(select(ref_p,select(con,cig)):*select(q,con))
		st_addobs(1) 								
		st_store(st_nobs(),range_c,pl_c_gk)			
		stata(`"replace var="pl_c_gk" in -1"')		
		stata(`"replace save_var=1 in -1"')			
	pl_i_gk=colsum(v_i):/colsum(select(ref_p,select(inv,cig)):*select(q,inv))
		st_addobs(1) 								
		st_store(st_nobs(),range_c,pl_i_gk)			
		stata(`"replace var="pl_i_gk" in -1"')		
		stata(`"replace save_var=1 in -1"')			
	pl_g_gk=colsum(v_g):/colsum(select(ref_p,select(gov,cig)):*select(q,gov))
		st_addobs(1) 								
		st_store(st_nobs(),range_c,pl_g_gk)			
		stata(`"replace var="pl_g_gk" in -1"')		
		stata(`"replace save_var=1 in -1"')			
end

keep if save_var==1									// Get rid of variables that are not of interest
drop cig bh_code bh_name value price xrate pop con inv gov tb save_var // Drop all indicator variables

* Transpose the dataset
mata
	data=st_data(.,2..c+1)'							// Grab all the data
	country_names=st_varname(2..st_nvar())'			// Grab the country names
	var_names=st_sdata(.,1)'						// Grab the variable names
	stata("clear")									// Clear the data from Stata
	st_addvar("str10",("countrycode"))				// Make the new variable countrycode
	st_addvar("double",("year",var_names))			// Make the year and other variables
	st_addobs(c)									// Add enough observations to include all countries
	st_sstore(.,1,country_names)					// Fill the country names
	st_store(.,3..cols(var_names)+2,data)			// Fill the variables
end
