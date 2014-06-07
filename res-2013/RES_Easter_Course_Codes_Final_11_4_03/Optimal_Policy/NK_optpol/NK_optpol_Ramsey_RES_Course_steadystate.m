%compute the steady state of NK Model 

function [ys,check]=NK_optpol_RES_Course(ys,exe);
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
%%
%NumberOfParameters = M_.param_nbr;                            % Number of deep parameters.
%for i = 1:NumberOfParameters                                  % Loop...
%  paramname = deblank(M_.param_names(i,:));                   %    Get the name of parameter i. 
%  eval([ paramname ' = M_.params(' int2str(i) ');']);         %    Get the value of parameter i.
%end                                                           % End of the loop.  
check = 0;
%%
params = M_.params;

gy=params(1);
alp=params(2);
c=params(3);
zzeta=params(4);
delta=params(5);
sigma_c=params(6);
rhoA=params(7);
rhoG=params(8);
rhoMS=params(9);
rhoMPS=params(10);
Ass=params(11);
phiX=params(12);
xi=params(13);
hab=params(14);
hss=params(15);
g=params(16);
gamp=params(17);
PIEss=params(18);
Rnss=params(19);
rhocapqual=params(20);
ppsi=params(21);
trend=params(22);
conspie=params(23);
consr=params(24);
varrho =params(25);
betta=params(26);
wr =params(27);
alpha_r =params(28);
alpha_y=params(29);
alpha_pie=params(30);
alpha_A =params(31);
alpha_G=params(32); 
alpha_MS =params(33);
alpha_capqual=params(34);
alpha_trend=params(35);


%% THIS BLOCK IS MODEL SPECIFIC.
%%
%% Here the user has to define the steady state.
%%

%
Stochg=g;
%


%This part calls the function fun_RBC.m
 %Use fsolve to find the steady state values of hours worked
 %initial values


x0=[log(hss/(1-hss)),log(hss/(1-hss))];
PIE=PIEss;
PIETILDE=PIE^(1-gamp);
 %options = optimset('TolFun',1e-16);
 %options = optimset('MaxIter',10);
[x,fval] =fsolve(@fun_NK_optpol_Ramsey_RES_Course,x0,optimset('Display','off'),PIE);
%[x,fval] =fsolve(@fun_RBC_Habit3_Ex_Growth_calib,x0,options,PIE);
%pause

%Derived variables (cut and paste from fun_RBC.m)
h=exp(x(1))/(1+exp(x(1)));
hF=exp(x(2))/(1+exp(x(2)));
bettag=betta*(1+g)^((1-varrho)*(1-sigma_c));
%h=hss;
A=Ass;
Rex=(1+g)^((1-varrho)*(sigma_c-1)+1)/betta;
Rn=PIE*Rex;
DD=1/Rex;
DDL=DD;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calvo Contracts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Delta=((1-xi)*(((1-xi*PIETILDE^(zzeta-1))/(1-xi))^(1/(1-zzeta)))^(-zzeta))/(1-xi*PIETILDE^zzeta);
MC=(1-1/zzeta)*(1-xi*betta*PIETILDE^zzeta*(1+g)^((1-varrho)*(1-sigma_c)))...
/(1-xi*betta*PIETILDE^(zzeta-1)*(1+g)^((1-varrho)*(1-sigma_c)))*(((1-xi*PIETILDE^(zzeta-1))...
/(1-xi))^(1/(1-zzeta)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rotemberg Contracts
%
%Replace above with with
%Delta=1;
%MC=1-(1-ppsi*(PIE-1)*PIE*(1-DD*(1+g)))/zzeta;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PWP=MC;

KY=(1-alp)*PWP/(Rex-1+delta)*(1+g);


YW=A*h*(KY/(1+g))^((1-alp)/alp)/Delta^(1/alp);

K=KY*YW;
I=(delta+g)*K/(1+g);
Y=(1-c)*YW;
G=gy*Y;
WP=alp*PWP*YW/h;
C=WP*((1-varrho)*(1-h))/((1-hab/(1+g))*varrho);
tax=G/(WP*h);

%Post recursive Steady state relationship
LAMBDA=((((1+Stochg)*C-hab*C)^(1-varrho)*(1-h)^varrho)^(1-sigma_c))/(1-sigma_c);
LAMBDAC=(1-varrho)*((C*(1-hab/(1+g)))^((1-varrho)*(1-sigma_c)-1)*(1-h)^(varrho*(1-sigma_c)));
XX=(Rex)*LAMBDAC;
Q=1;
X=1+g;
Z1=2.0*phiX*(X-1-g)*X^2*Q/(Rex);
Z2=(1-alp)*PWP*YW/(K/(1+g))+(1-delta)*Q;

H=Y*LAMBDAC/(1-betta*xi*PIETILDE^(zzeta-1)*(1+g)^((1-varrho)*(1-sigma_c)));
Htilde=(PIETILDE^(zzeta-1))*H;
J=(1/(1-1/zzeta))*Y*LAMBDAC*MC/(1-betta*xi*PIETILDE^(zzeta)*(1+g)^((1-varrho)*(1-sigma_c)));
Jtilde=PIETILDE^zzeta*J;

INVPIE=1/PIE;
YY=1;
CC=1;
hh=1;
WPWP=1; 
II=1; 
KK=1;
RR=1;
RnRn=1;
QQ=1;
PIEPIE=1;
Rk=Rex;
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
Z1F=Z1;
Z2F=Z2;
RkF=Rk; 
spreadF=spread; 
RRF=RR;
YYF=YY;
CCF=CC;
hhF=hh;
WPWPF=WPWP;
IIF=II;
KKF=KK;
QQF=1;
DDF=DD;
DDFL=DDF;
OUTGAP=YF/Y;
OUTGAPOUTGAP=1;
capqual=1;
OMEGA=LAMBDA;
OMEGAscaled=100*LAMBDA;
LAMBDAtrue=LAMBDA;
OMEGAtrue=LAMBDAtrue;
OMEGAtruescaled=100*LAMBDAtrue;
OMEGAOMEGA=-1;
CE=((((1+Stochg-hab)*C*1.01)^(1-varrho)*(1-h)^varrho)^(1-sigma_c))/(1-sigma_c)-LAMBDAtrue;

RkRk=1;
spreadspread=0;
phiphi=1;
NWNW=1;
dy=trend;
pinfobs=log(PIEPIE)+conspie;
robs=log(RnRn)+consr;

%% END OF THE MODEL SPECIFIC BLOCK.


%% DO NOT CHANGE THIS PART.
%%
%% Here we define the steady state vZNues of the endogenous variables of
%% the model.
%%
NumberOfEndogenousVariables = M_.endo_nbr;                    % Number of endogenous variables.
ys = zeros(NumberOfEndogenousVariables,1);                    % Initialization of ys (steady state).
for i = 1:NumberOfEndogenousVariables                         % Loop...
  varname = deblank(M_.endo_names(i,:));                      %    Get the name of endogenous variable i.                     
  eval(['ys(' int2str(i) ') = ' varname ';']);                %    Get the steady state vZNue of this variable.
end                                                           % End of the loop.
%%


