clear *
mata

function fischer(v,p,base)
{
	p_lasp=laspeyres(v,p,base)
	p_paasche=paasche(v,p,base)
	p_fischer=sqrt(p_lasp:*p_paasche)
	return(p_fischer)
}

mata mosave fischer(), replace

end

