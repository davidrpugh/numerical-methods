clear *
mata

function paasche(v,p,base)
{
	p=p:/p[,base]
	p_paasche=colsum(v):/colsum(v:/p)
	return(p_paasche)
}

mata mosave paasche(), replace

end

