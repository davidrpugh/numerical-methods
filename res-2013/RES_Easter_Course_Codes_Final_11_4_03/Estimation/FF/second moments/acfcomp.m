function acf = acfcomp(z)
% ACFCOMP(z) computes sample ACF
% z should be a column vector
n=length(z);
w=2/sqrt(n);
acvf=autocov(z);
na=length(acvf)-1;
acf=acvf/acvf(1);

