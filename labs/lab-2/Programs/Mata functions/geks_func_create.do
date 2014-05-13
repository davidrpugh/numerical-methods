clear
mata
mata clear

function geks(v,p,base)
{
	P=J(cols(p),cols(p),.)
	for(b=1;b<=cols(p);b++) P[b,.]=fischer(v,p,b)
	ml=mean(log(P)')'
	ML=ml*J(1,cols(p),1)
	GEKS=exp(ML-ML')
	p_agg=GEKS[base,]
	return(p_agg)
}

mata mosave geks(), replace

end
