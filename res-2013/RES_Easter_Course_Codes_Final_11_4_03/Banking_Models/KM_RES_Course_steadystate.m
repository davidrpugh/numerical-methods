%compute the steady state of NK Model 

function [ys,check]=KM_Course_1(ys,exe);
global M_

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
%


PIE=PIEss;
PIETILDE=PIE^(1-gamp);
 %options = optimset('TolFun',1e-16);
 %options = optimset('MaxIter',10);
 rho=0.8; 
x0=[log(rho/(1-rho)),log(hss/(1-hss)), 1, 0.5, 0, 0.001];
%x0=[2.2729,-0.7838, 1.8302, 4.6272, -0.0550, 6.3426];
[x,fval] =fsolve(@fun_KM_RES_Course,x0,optimset('Display','off'),PIE);

%
%to reset initial values uncomment:
%x
%pause
%
%and replace x0 above

%Derived variables
%
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


