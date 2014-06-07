%compute the steady state of NK Model 

function [ys,check]=NK_Rotemberg_RES_Course_2(ys,exe);
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
NumberOfParameters = M_.param_nbr;                            % Number of deep parameters.
for i = 1:NumberOfParameters                                  % Loop...
  paramname = deblank(M_.param_names(i,:));                   %    Get the name of parameter i. 
  eval([ paramname ' = M_.params(' int2str(i) ');']);         %    Get the value of parameter i.
end                                                           % End of the loop.  
check = 0;
%%


%% THIS BLOCK IS MODEL SPECIFIC.
%%
%% Here the user has to define the steady state.
%%

%
Stochg=g;

%This part calls the function fun_RBC.m
 %Use fsolve to find the steady state values of hours worked
 %initial values
rho=0.8;
x0=[log(rho/(1-rho)),log(hss/(1-hss)), log(rho/(1-rho))];
PIE=PIEss;
 %options = optimset('TolFun',1e-16);
 %options = optimset('MaxIter',10);
[x,fval] =fsolve(@fun_NK_Rotemberg_RES_Course,x0,optimset('Display','off'),PIE);
%[x,fval] =fsolve(@fun_RBC_Habit3_Ex_Growth_calib,x0,options,PIE);
%pause

%Derived variables (cut and paste from fun_RBC.m)
varrho=exp(x(1))/(1+exp(x(1)));
hF=exp(x(2))/(1+exp(x(2)));
betta=exp(x(3))/(1+exp(x(3)));
h=hss;
A=Ass;
Rex=(1+g)^((1-varrho)*(sigma_c-1)+1)/betta;
DD=1/Rex;
MC=1-(1-ppsi*(PIE-1)*PIE*(1-DD*(1+g)))/zzeta;
PWP=MC;
KY=(1-alp)*PWP/(Rex-1+delta)*(1+g);

%%%%%%%%%%%YW=A*h*KY^((1-alp)/alp);
YW=A*h*(KY/(1+g))^((1-alp)/alp);

K=KY*YW;
I=(delta+g)*K/(1+g);
Y=(1-c)*YW;
G=gy*Y;
WP=alp*PWP*YW/h;
C=WP*((1-varrho)*(1-h))/((1-hab/(1+g))*varrho);
tax=G/(WP*h);

%Post recursive Steady state relationship
LAMBDA=1/(1-sigma_c)*((C*(1-hab/(1+g)))^((1-varrho)*(1-sigma_c))*(1-h)^(varrho*(1-sigma_c))-1);
LAMBDAC=(1-varrho)*((C*(1-hab/(1+g)))^((1-varrho)*(1-sigma_c)-1)*(1-h)^(varrho*(1-sigma_c)));
OMEGA=LAMBDA;
XX=(Rex)*LAMBDAC;
Q=1;
X=1+g;
Z1=2.0*phiX*(X-1-g)*X^2*Q/(Rex);
Z2=(1-alp)*PWP*YW/(K/(1+g))+(1-delta)*Q;



%%%%%%%%%%%Rn=1/betta; - Rn should be exogenous 
Rn=PIE*(1+g)^((1-varrho)*(sigma_c-1)+1)/betta;

%%%%%%%%%%%INVPIE=1;
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


%%%%%%%%%%%Delta=1;

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
OUTGAP=YF/Y;
OUTGAPOUTGAP=1;
capqual=1;

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


