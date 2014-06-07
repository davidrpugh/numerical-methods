%compute the steady state of NK Model 

function [ys,check]=KM_Course_6_MCMC_steadystate(ys,exe);
global M_

% add this for ACES only - DO NOT CHANGE THIS PART.
%global lq_instruments 
%%%%%%%%%%%% add these for LQ-ACEs
%iPI= strmatch(lq_instruments.seedvar,M_.endo_names, 'exact');
%eval([ lq_instruments.seedvar ' = ys(iPI);']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DO NOT CHANGE THIS PART.
%%
%% Here we load the values of the deep parameters in a loop.
% %%
% NumberOfParameters = M_.param_nbr;                            % Number of deep parameters.
% for i = 1:NumberOfParameters                                  % Loop...
%   paramname = deblank(M_.param_names(i,:));                   %    Get the name of parameter i. 
%   eval([ paramname ' = M_.params(' int2str(i) ');']);         %    Get the value of parameter i.
% end                                                           % End of the loop.  
check = 0;
%%


%% THIS BLOCK IS MODEL SPECIFIC.
%%
%% Here the user has to define the steady state.
%%


params = M_.params;

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
%


Stochg=g;
PIE=PIEss;
PIETILDE=PIE^(1-gamp);


%This part calls the function fun_RBC.m
 %Use fsolve to find the steady state values of hours worked
 %initial values

rho=0.8; 
x0=[log(rho/(1-rho)),log(hss/(1-hss)), 1, 0.5, 0, 0.001];

 %options = optimset('TolFun',1e-16);
 %options = optimset('MaxIter',10);
[x,fval] =fsolve(@fun_KM_RES_Course,x0,optimset('Display','off'),PIE,params);
%[x,fval] =fsolve(@fun_RBC_Habit3_Ex_Growth_calib,x0,options,PIE);
%pause

%Derived variables 
varrho=exp(x(1))/(1+exp(x(1)));
hF=exp(x(2))/(1+exp(x(2)));
K=exp(x(3));
bettaE=exp(x(4))/(1+exp(x(4)));
TE=x(5);
betta=exp(x(6))/(1+exp(x(6)));
%
%
h=hss;
A=Ass;
Rex=(1+g)^((1-varrho)*(sigma_c-1)+1)/betta;
Rn=Rex*PIE;
DD=1/Rex;
DDL=DD;
Delta=((1-xi)*(((1-xi*PIETILDE^(zzeta-1))/(1-xi))^(1/(1-zzeta)))^(-zzeta))/(1-xi*PIETILDE^zzeta);
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

%Post recursive Steady state relationship
LAMBDA=1/(1-sigma_c)*((C*(1-hab/(1+g)))^((1-varrho)*(1-sigma_c))*(1-h)^(varrho*(1-sigma_c))-1);
LAMBDAC=(1-varrho)*((C*(1-hab/(1+g)))^((1-varrho)*(1-sigma_c)-1)*(1-h)^(varrho*(1-sigma_c)));
XX=(Rex)*LAMBDAC;
Q=1;
X=1+g;
Z1=2.0*phiX*(X-1-g)*X^2*Q/(Rex);
Z2=(1-alp)*PWP*YW/(K/(1+g))+(1-delta)*Q;

H=Y*LAMBDAC/(1-betta*xi*PIETILDE^(zzeta-1)*(1+g)^((1-varrho)*(1-sigma_c)));


Htilde=(PIETILDE^(zzeta-1))*H;

%J=H*(((1-xi*PIETILDE^(zzeta-1))/(1-xi))^(1/(1-zzeta)));
J=(1/(1-1/zzeta))*Y*LAMBDAC*MC/(1-betta*xi*PIETILDE^(zzeta)*(1+g)^((1-varrho)*(1-sigma_c)));

Jtilde=PIETILDE^zzeta*J;





Rn=PIE*(1+g)^((1-varrho)*(sigma_c-1)+1)/betta;


INVPIE=1/PIE;

Z=(1-alp)*PWP*YW/(K/(1+g));
DDE=bettaE/(1+g)^sigma_entre;
QK=Q*K;
RL=1/(1-1/zzetaL)*Rn;
ThetaE=1/RL-DDE/PIE;
Rk=(Z+(1-delta)*Q)/Q;
L=m*PIE*Q*(1-delta)*K/RL;
RLex=RL/PIE;
CE=(Rk/(1+g)-1)*Q*K+L*(1-RLex/(1+g))-TE;
LAMBDAE=((CE*(1-habE/(1+Stochg)))^(1-sigma_entre)-1)/(1-sigma_entre);
LAMBDACE=(CE*(1-habE/(1+Stochg)))^(-sigma_entre);
DDERL=DDE*RLex;
DDERk=DDE*Rk;
phi=Q*K/(Q*K-L);
phiphi=1;
YY=1;
CC=1;
CECE=1;
hh=1;
WPWP=1; 
II=1; 
KK=1;
RR=1;
RnRn=1;
QQ=1;
PIEPIE=1;
spread=Rk-Rex;




MS=1;
MPS=1;

%
%Flexi-Price ss
%

RF=Rex;
MCF=(1-1/zzeta);
PWPF=MCF;
KYF=(1-alp)*PWPF/(RF-1+delta)*(1+g);


YWF=A*hF*(KYF/(1+g))^((1-alp)/alp);

KF=KYF*YWF;
IF=(delta+g)*KF/(1+g);
YF=(1-c)*YWF;
WPF=alp*PWPF*YWF/hF;
CF=WPF*((1-varrho)*(1-hF))/((1-hab/(1+g))*varrho);
GF=YF*G/Y;
taxF=GF/(WPF*hF);

LAMBDAF=1/(1-sigma_c)*((CF*(1-hab/(1+g)))^((1-varrho)*(1-sigma_c))*(1-hF)^(varrho*(1-sigma_c))-1);
LAMBDACF=(1-varrho)*((CF*(1-hab/(1+g)))^((1-varrho)*(1-sigma_c)-1)*(1-hF)^(varrho*(1-sigma_c)));
XF=X;
QF=Q; 
XXF=1; 
DDF=DD;
DDFL=DDF;
Z1F=2.0*phiX*(XF-1-Stochg)*XF^2*QF*DDF;
Z2F=(1-alp)*PWPF*YWF/(KF/(1+g))+(1-delta)*QF;
RkF=Z2F/QF; 
spreadF=0; 
RRF=RR;
YYF=YY;
CCF=CC;
hhF=hh;
WPWPF=WPWP;
IIF=II;
KKF=KK;
QQF=1;
OUTGAP=YF/Y;
OUTGAPOUTGAP=1;
capqual=1;
OMEGA=LAMBDA;
OMEGAOMEGA=-1;
RkRk=1;
spreadspread=0;
DDEL=DDE;

dy=trend;
pinfobs=log(PIEPIE)+conspie;
robs=log(RnRn)+consr;
rkn_obs = log(RkRk*PIEPIE)+consrkn;

ys=[LAMBDA OMEGA LAMBDAC Rex h WP YW Y  PWP K I  tax  C A  G X Q XX ...
Z1 Z2 MC H Htilde J Jtilde PIE Rn INVPIE Rk spread ...
RR RnRn YY CC hh WPWP II KK Delta  MS MPS ...
 varrho PIETILDE  PIEPIE QQ DD DDL ...
LAMBDAF LAMBDACF RF hF WPF YWF YF  PWPF KF IF  taxF  CF  XF QF ...
Z1F Z2F MCF  RkF spreadF  OUTGAP ...
RRF YYF CCF hhF WPWPF IIF KKF  QQF OUTGAPOUTGAP OMEGAOMEGA ...
 GF Stochg capqual LAMBDAE LAMBDACE CE DDE DDEL RL RLex DDERL DDERk QK L TE ThetaE bettaE RkRk spreadspread ...
phi phiphi CECE betta dy pinfobs robs rkn_obs]';
