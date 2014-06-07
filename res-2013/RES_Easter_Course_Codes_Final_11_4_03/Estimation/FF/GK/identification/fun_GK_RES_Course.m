% ------------------------------------------------------
% This program defines the conditions needed to find the
% steady state value of hours worked using fsolve
% ------------------------------------------------------

function y=fun_GK_Course_F_4(x,PIE,params)

check = 0;

gy = params(1);
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
sigmaB = params(17);
lev = params(18);
creditspread = params(19);
hab = params(20);
g = params(21);
gamp = params(22);
PIEss = params(23);
hss = params(24);
rhocapqual = params(25);
Rnss = params(26);
trend = params(27);
conspie = params(28);
consr = params(29);
consrkn = params(30);

%%Define hours as unknown
varrho=exp(x(1))/(1+exp(x(1)));
K=x(2);
ThetaB=x(3);
xiB=x(4);
hF=exp(x(5))/(1+exp(x(5)));
betta=exp(x(6))/(1+exp(x(6)));


%Define only the steady state relationships needed to find the value of hours worked that supports
%the resource constraint relationship Y=C+I+G
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

S=K;
Z=(1-alp)*PWP*YW/(K/(1+g));
phiB=lev;
NW=(Q*S)/phiB;
Dep=(Q*S)-NW;
omega=1-sigmaB+sigmaB*ThetaB*phiB;
nuB=omega;
muB=(phiB*ThetaB-nuB)/phiB;
Rk=((muB)/(DD*omega)+(Rex));
spread=Rk-Rex;
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
y=[Y-C-I-G; YF-CF-IF-GF; 
   Z+(1-delta)-Rk; 
   creditspread-spread;
 NW-((((sigmaB+xiB)*(Z+(1-delta))*Q*S/(1+g))-(sigmaB*(Rex)*Q*S/(1+g)))/(1-sigmaB*(Rex)/(1+g)));
%NW-(sigmaB+xiB)*(Z+(1-delta))*Q*S/(1+g)+sigmaB*Rex*Dep/(1+g);
Rex-Rn/PIE;];




