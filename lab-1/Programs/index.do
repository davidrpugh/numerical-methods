capture program drop index
program index
	// Argument list
	// 1: value
	// 2: growth rate
	// 3: base year
	// 4: by group variable
	// 5: new variable
	tempvar cum_growth val_base cum_growth_base
	bysort `4': gen `cum_growth'=exp(sum(`2'))
	bysort `4': gen `val_base'=sum(`1'*`3')
	bysort `4': replace `val_base'=`val_base'[_N]
	bysort `4': gen `cum_growth_base'=sum(`cum_growth'*`3')
	bysort `4': replace `cum_growth_base' = `cum_growth_base'[_N]
	bysort `4': gen `5'=`val_base'*`cum_growth'/`cum_growth_base'
end
