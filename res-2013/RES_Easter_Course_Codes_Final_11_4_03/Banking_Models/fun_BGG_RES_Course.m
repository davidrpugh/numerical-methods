% ------------------------------------------------------
% This program defines the conditions needed to find the
% steady state value of hours worked using fsolve
% ------------------------------------------------------

function y=fun_BGG_RES_Course(x,PIE)
global M_
 
%% DO NOT CHANGE THIS PART.
%%
%% Here we load the values of the deep parameters in a loop.
%%
NumberOfParameters = M_.param_nbr;                            % Number of deep parameters.
for i = 1:NumberOfParameters                                  % Loop...
  paramname = deblank(M_.param_names(i,:));                   %    Get the name of parameter i. 
  eval([ paramname ' = M_.params(' int2str(i) ');']);         %    Get the value of parameter i.
end                                                           % End of the loop.  
check = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Derived variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varrho=exp(x(1))/(1+exp(x(1)));
K=x(2);
sigmaE=exp(x(3))/(1+exp(x(3)));
Apsi=x(4);
hF=exp(x(5))/(1+exp(x(5)));
psi=x(6);
mu=exp(x(7))/(1+exp(x(7))); 
xiE=exp(x(8))/(1+exp(x(8)));
betta=exp(x(9))/(1+exp(x(9)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Define only the steady state relationships needed to find the value of hours worked that supports
%the resource constraint relationship Y=C+I+G
%
PIETILDE=PIE^(1-gamp);
Delta=((1-xi)*(((1-xi*PIETILDE^(zzeta-1))/(1-xi))^(1/(1-zzeta)))^(-zzeta))/(1-xi*PIETILDE^zzeta);
h=hss;
A=Ass;
Q=1;
Rn=Rnss;
Rex=(1+g)^((1-varrho)*(sigma_c-1)+1)/betta;
DD=1/Rex;

MC=(1-1/zzeta)*(1-xi*betta*PIETILDE^zzeta*(1+g)^((1-varrho)*(1-sigma_c)))...
/(1-xi*betta*PIETILDE^(zzeta-1)*(1+g)^((1-varrho)*(1-sigma_c)))*(((1-xi*PIETILDE^(zzeta-1))...
/(1-xi))^(1/(1-zzeta)));
PWP=MC;

YW=A*h^(alp)*(K/(1+g))^(1-alp)/Delta;
I=(delta+g)*K/(1+g);
Y=(1-c)*YW;
G=gy*Y;
WP=alp*PWP*YW/h;
C=WP*((1-varrho)*(1-h))/((1-hab/(1+g))*varrho);
tax=G/(WP*h);
%
%Banking Sector
%
Z=(1-alp)*PWP*YW/(K/(1+g));
Rk=Z+1-delta;
p=1/(2*Apsi)*(psi-1+epsilonA+Apsi);
fnGam=1/(4*Apsi)*(psi^2-(1-epsilonA-epsilonA-Apsi)^2)+psi*(1-p);
fnG=1/(4*Apsi)*(psi^2-(1-epsilonA-Apsi)^2);
DGam=1/(2*Apsi)*(1-epsilonA-psi)+1/2;
DG=psi/(2*Apsi);
rho=DGam/((fnGam-mu*fnG)*DGam+(1-fnGam)*(DGam-mu*DG));
NW=(sigmaE+xiE)*(1-fnGam)*Rk*Q*K/(1+g);
phi=Q*K/NW;
spread=Rk-Rex;
CE=(1-xiE)*(1-sigmaE)*(1-fnGam)*Rk*Q*K/(1+g);

%
%Flexi-Price 
%
RexF=Rex;
MCF=(1-1/zzeta);
PWPF=MCF;

KYF=(1-alp)*PWPF/(RexF-1+delta)*(1+g);

YWF=A*hF*(KYF/(1+g))^((1-alp)/alp);

KF=KYF*YWF;
IF=(delta+g)*KF/(1+g);
YF=(1-c)*YWF;
GF=gy*YF;
WPF=alp*PWPF*YWF/hF;
CF=WPF*((1-varrho)*(1-hF))/((1-hab/(1+g))*varrho);
taxF=GF/(WPF*hF);

%Resource constraint 
%this is the equation for which fsolve will find the solution
y=[Y-C-CE-I-G-mu*fnG*Rk*Q*K/(1+g); YF-CF-IF-GF; 
   rho-rhocalib;  p-pcalib; phi-phicalib; CE/Y-CEcalib;
   Rk-rho*Rex; phi*Rk*(fnGam-mu*fnG)-Rex*(phi-1)
   Rex-Rn/PIE;];



