clear *
mata

function tornquist(v,p,base)
{
	sh=v:/colsum(v)
	avsh=.5*(sh:+sh[.,base])
	p_t=exp(colsum(avsh:*log(p:/p[.,base])))
	return(p_t)
}

mata mosave tornquist(), replace

end

