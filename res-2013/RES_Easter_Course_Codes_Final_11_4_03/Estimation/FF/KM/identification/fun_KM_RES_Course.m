% ------------------------------------------------------
% This program defines the conditions needed to find the
% steady state value of hours worked using fsolve
% ------------------------------------------------------

function y=fun_NK_Course_Revised_4_F(x,PIE,params)

check = 0;

gy=params(1); 
alp=params(2);
c=params(3);
zzeta=params(4);
delta=params(5);
sigma_c=params(6);
sigma_entre=params(7);
rhoA=params(8);
rhoG=params(9);
rhoMS=params(10);
rhoMPS=params(11);
Ass=params(12);
phiX=params(13);
xi=params(14);
alpha_y=params(15);
hab=params(16);
hss=params(17);
g=params(18);
gamp=params(19);
PIEss=params(20);
Rnss=params(21);
alpha_r=params(22);
alpha_pie=params(23);
rhocapqual=params(24);
habE=params(25);
zzetaL=params(26);
CEcalib=params(27);
spreadcalib=params(28);
m=params(29);
trend=params(30);
conspie=params(31);
consr=params(32);
consrkn=params(33);



%%Define hours as unknown
varrho=exp(x(1))/(1+exp(x(1)));
hF=exp(x(2))/(1+exp(x(2)));
K=exp(x(3));
bettaE=exp(x(4))/(1+exp(x(4)));
TE=x(5);
betta=exp(x(6))/(1+exp(x(6)));
%Define only the steady state relationships needed to find the value of hours worked that supports
%the resource constraint relationship Y=C+I+G

%%%%%%%%%%%define QQ=J/H
%QQ=((1-xi*PIE^(zzeta-1))/(1-xi))^(1/(1-zzeta));
%%%%%%%%%%%dispersion
PIETILDE=PIE^(1-gamp);
Delta=((1-xi)*(((1-xi*PIETILDE^(zzeta-1))/(1-xi))^(1/(1-zzeta)))^(-zzeta))/(1-xi*PIETILDE^zzeta);
h=hss;
A=Ass;
Rex=(1+g)^((1-varrho)*(sigma_c-1)+1)/betta;
Rn=Rex*PIE;
DD=1/Rex;
%%%%%%%%%%%PWP=1-1/zzeta;
MC=(1-1/zzeta)*(1-xi*betta*PIETILDE^zzeta*(1+g)^((1-varrho)*(1-sigma_c)))...
/(1-xi*betta*PIETILDE^(zzeta-1)*(1+g)^((1-varrho)*(1-sigma_c)))*(((1-xi*PIETILDE^(zzeta-1))...
/(1-xi))^(1/(1-zzeta)));
PWP=MC;


%%%%%%%%%%%YW=A*h*KY^((1-alp)/alp);

YW=A*h^(alp)*(K/(1+g))^(1-alp)/Delta;

I=(delta+g)*K/(1+g);
Y=(1-c)*YW;
G=gy*Y;
WP=alp*PWP*YW/h;
C=WP*((1-varrho)*(1-h))/((1-hab/(1+g))*varrho);
tax=G/(WP*h);
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
Z=(1-alp)*PWP*YW/(K/(1+g));
Q=1;
DDE=bettaE/(1+g)^sigma_entre;
RL=1/(1-1/zzetaL)*Rn;
ThetaE=1/RL-DDE/PIE;
Rk=(Z+(1-delta)*Q)/Q;
RLex=RL/PIE;
L=m*PIE*Q*(1-delta)*K/RL;
CE=(Rk/(1+g)-1)*Q*K+L*(1-RLex/(1+g))-TE;
spread=Rk-Rex;


y=[Y-C-CE-I-G; YF-CF-IF-GF;...
Rk-(1-ThetaE*m*PIE*Q*(1-delta))/DDE;...
spread-spreadcalib; CE/Y-CEcalib; Rn-Rnss];



