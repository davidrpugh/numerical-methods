% ------------------------------------------------------
% This program defines the conditions needed to find the
% steady state value of hours worked using fsolve
% ------------------------------------------------------

function y=fun_BGG_Course_3(x,PIE,params)

                                                         % End of the loop.  
check = 0;
%%

gy  = params(1);
alp = params(2);
c = params(3);
zzeta = params(4);
delta = params(5);
sigma_c = params(6); 
rhoA = params(7);
rhoG = params(8);
rhoMS = params(9);
rhoMPS = params(10);
Ass = params(11);
phiX = params(12);
xi = params(13);
alpha_r = params(14);
alpha_pie = params(15);
alpha_y = params(16);
pcalib = params(17);
rhocalib = params(18);
phicalib = params(19);
CEcalib = params(20);
hab = params(21);
g = params(22);
gamp = params(23);
PIEss = params(24);
hss = params(25);
rhocapqual = params(26);
alpha_q  = params(27);
epsilonA = params(28);
Rnss= params(29);
trend = params(30);
conspie = params(31);
consr = params(32);
consrkn= params(33);

%%Define hours as unknown
varrho=exp(x(1))/(1+exp(x(1)));
K=x(2);
sigmaE=exp(x(3))/(1+exp(x(3)));
Apsi=x(4);
hF=exp(x(5))/(1+exp(x(5)));
psi=x(6);
mu=exp(x(7))/(1+exp(x(7))); 
xiE=exp(x(8))/(1+exp(x(8)));
betta=exp(x(9))/(1+exp(x(9)));

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
%
%PWP=1-1/zzeta;
%
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

%%%%%%%%%%%YW=A*h*KY^((1-alp)/alp);
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



