clear
mata
mata clear

function gk(v,p,base)
{
	p_gk=J(100,cols(p),.)
	p_gk[1,]=J(1,cols(p),1)
	i=1
	dif=1000
	q=v:/p
	while(dif>1e-8)
		{
		ref_p=rowsum(v:/p_gk[i,]):/rowsum(q)
		i++
		p_gk[i,]=colsum(v):/colsum(ref_p:*q)
		p_gk[i,]=p_gk[i,]:/p_gk[i,base]
		dif=sum(abs(log(p_gk[i,]:/p_gk[i-1,])))
		}
	p_gk=p_gk[i,]
	return(p_gk)
}

mata mosave gk(), replace

end

