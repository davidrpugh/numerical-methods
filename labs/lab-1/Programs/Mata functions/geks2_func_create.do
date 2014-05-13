clear
mata
mata clear

function geks2(v,p,base)
{
	P=J(cols(p),cols(p),.)
	for(b=1;b<=cols(p);b++) P[b,.]=tornquist(v,p,b)
	ml=mean(log(P)')'
	ML=ml*J(1,cols(p),1)
	GEKS=exp(ML-ML')
	p_agg=GEKS[base,]
	return(p_agg)
}

mata mosave geks2(), replace

end
