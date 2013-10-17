clear *
mata

function laspeyres(v,p,base)
{
	p=p:/p[,base]
	p_lasp=colsum(v[,base]:*p):/colsum(v[,base])
	return(p_lasp)
}

mata mosave laspeyres(), replace

end

