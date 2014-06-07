%
%program to work out model odds from LL values
%
%
%number of models
%
n=3
%
%LL vector (GH H Z)
%

LL=[-127.338556 	-93.612833 	-251.341045 ]
sum=1.0;
for i=2:n
      rp(i)=exp(LL(i)-LL(1));
      sum=sum+rp(i);
end
p(1)=1/sum;
sump=p(1);
for i=2:n
      p(i)=rp(i)*p(1);
      sump=sump+p(i);
end
p
sump

