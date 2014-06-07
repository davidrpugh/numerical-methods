/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Kiyotaki-Moore-Collatoralized Constraint Model
///////////////////////////////////////////////////////////////////////////////////////////////// 
//(c) CIMS Univeristy of Surrey
//The Science and  Art of DSGE Modelling: Construction, Calibration, Estimation and Policy
//////////////////////////////////////////////////////////////////////////////////////////////////
//This is a free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.  See <http://www.gnu.org/licenses/> for more information.
////////////////////////////////////////////////////////////////////////////////////////////////// 
//
///
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF ENDOGENOUS VARIABLES%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var LAMBDA OMEGA LAMBDAC Rex h WP YW Y  PWP K I  tax  C A  G X Q XX 
Z1 Z2 MC H Htilde J Jtilde PIE Rn INVPIE Rk spread 
RR RnRn YY CC hh WPWP II KK Delta  MS MPS
 varrho PIETILDE  PIEPIE QQ DD DDL
LAMBDAF LAMBDACF RF hF WPF YWF YF  PWPF KF IF  taxF  CF  XF QF 
Z1F Z2F MCF  RkF spreadF  OUTGAP 
RRF YYF CCF hhF WPWPF IIF KKF  QQF OUTGAPOUTGAP OMEGAOMEGA
 GF Stochg capqual LAMBDAE LAMBDACE CE DDE DDEL RL RLex DDERL DDERk QK L TE ThetaE bettaE RkRk spreadspread
phi phiphi CECE betta;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF OBSERVABLE VARIABLES%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var dy pinfobs robs rkn_obs; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF EXOGENOUS VARIABLES%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varexo epsA epsG epsMS epsMPS epsAtrend epscapqual;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF PARAMETERS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters gy  alp c zzeta  delta sigma_c sigma_entre rhoA rhoG rhoMS rhoMPS
Ass  phiX xi alpha_y hab  hss g gamp PIEss Rnss
alpha_r $\alpha_r$
alpha_pie $\alpha_{\pi}$ rhocapqual
habE   zzetaL CEcalib spreadcalib m ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%ADD OBSERVATION TRENDS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters trend conspie consr consrkn;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PARAMETERS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PIEss=1.0063;
Rnss=1.013142;
g=0.0046;
gamp=0.1;
gy=0.2;
hss=0.35;
alp=0.70;
zzeta=7.0;
c=1/zzeta;
delta=0.0250;
sigma_c=2.0;
phiX=1.24;
xi=0.75;
///////////////////
//Flexi-Price Case
////////////////////
//xi=0.00001;
/////////////////////////////
//MP rule
//////////////////////////////
alpha_r=0.7;
alpha_pie=2.0;
alpha_y=0.5;

hab=0.7;
///////////////////////
//Choice of Units
//////////////////////
Ass=1;
//////////////////////
//
////////////////////
///shock persistence
//////////////////////
rhoA=0.7;
rhoG=0.7;
rhoMS=0.7;
rhoMPS=0;
rhocapqual=0.7;
//////////////////////
//
// FF Parameters
/////////////////////
habE=hab;
sigma_entre=1.5;
m=0.25;
zzetaL=500;
spreadcalib=0.01/4;
CEcalib=0.1;//as in BGG

trend=g*100;
consr=(Rnss-1)*100; % ss nominal interest rate
conspie=100*(PIEss-1); % quarterly ss inflation rate
consrkn=2.1955;% ss nominal BAA corporate bonds

//
// ----------------------------
// *** DSGE-Model-equations ***
// ----------------------------
model;

%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%ccalibration of varrho and betta - eqns 1, 2
%%%%%%%%%%%%%%%%%%%%%%%%%
varrho=STEADY_STATE(varrho);
betta=STEADY_STATE(betta);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Stochastic Trend
%
log(1+Stochg)=log(1+g)+epsAtrend;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Euler Eqn, Single period utility and Marginal utilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

XX=(Rn(-1)/PIE)*LAMBDAC;
LAMBDAC=betta*(1+Stochg(+1))^((1-varrho)*(1-sigma_c)-1)*XX(+1);
LAMBDA=((((C-hab*C(-1)/(1+Stochg))^(1-varrho))*((1-h)^varrho))^(1-sigma_c)-1)/(1-sigma_c);
LAMBDAC=(1-varrho)*((C-hab*C(-1)/(1+Stochg))^((1-varrho)*(1-sigma_c)-1))*((1-h)^(varrho*(1-sigma_c)));
%
%Stochastic Discount Factor =correct date on growth
%
DD=betta*(LAMBDAC(+1)/LAMBDAC)*(1+Stochg(+1))^((1-varrho)*(1-sigma_c)-1);
DDL=betta*(LAMBDAC/LAMBDAC(-1))*(1+Stochg)^((1-varrho)*(1-sigma_c)-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Intertemporal welfare En 8
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OMEGA=(1-betta)*LAMBDA+betta*OMEGA(+1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Household FOC for labour Supply
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varrho*((C-hab*C(-1)/(1+Stochg))^((1-varrho)*(1-sigma_c)))*((1-h)^(varrho*(1-sigma_c)-1))/LAMBDAC=WP;
%%%%%%%%%%%%%%%%
%Retail Output
%%%%%%%%%%%%%%%%
Y=(1-c)*YW;
%%%%%%%%%%%%%%%%%%%%%%%%%
%Wholesale Firm FOC
%%%%%%%%%%%%%%%%%%%%%%%%% 
YW=((A*h)^alp)*(K(-1)/(1+Stochg))^(1-alp)/Delta;
PWP*(alp*YW)/h=WP;
%%%%%%%%%%%%%%%%%%%%
%Capital Producers Eqn 13
%%%%%%%%%%%%%%%%%%%
X=I/I(-1)*(1+Stochg);
K=capqual(+1)*((I*(1-phiX*(X-1-Stochg)^2.0)+(1-delta)*K(-1)/(1+Stochg)));
Z1=2.0*phiX*(X-1-Stochg)*X^2*Q*DDL;
Q*(1- phiX*(X-1-Stochg)^2.0-2.0*X*phiX*(X-1-Stochg))+Z1(+1)=1;
INVPIE=1/PIE;
Rex=Rn(-1)*INVPIE;
Z2=(1-alp)*PWP*YW/(K(-1)/(1+Stochg))+(1-delta)*Q;
Rk=capqual*Z2/Q(-1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%REPLACE THE FOLLOWING WITH CC RELATIONSHIPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rn*INVPIE(+1)=Rk(+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Resource Constraint (Output Equilibrium)- 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Y=C+G+I;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CC MODEL Eq 21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bettaE=STEADY_STATE(bettaE);
LAMBDAE=((CE-habE*CE(-1)/(1+Stochg))^(1-sigma_entre)-1)/(1-sigma_entre);
LAMBDACE=(CE-habE*CE(-1)/(1+Stochg))^(-sigma_entre);
%
%Stochastic Discount Factor =correct date on growth
%
DDE=bettaE*(LAMBDACE(+1)/LAMBDACE)*(1+Stochg(+1))^(-sigma_entre);
DDEL=bettaE*(LAMBDACE/LAMBDACE(-1))*(1+Stochg)^(-sigma_entre);
RL=1/(1-1/zzetaL)*Rn;
RLex=RL(-1)/PIE;
DDERL=DDEL*RLex;
DDERL(+1)+ThetaE*RL=1;
DDERk=DDEL*Rk;
//
DDERk(+1)+ThetaE*m*(1-delta)*PIE(+1)*Q(+1)/Q=1;
TE=STEADY_STATE(TE);
%
%plus binding constraint 32
%
QK=Q*K(-1);
RL*L=m*(1-delta)*QK(+1)*PIE(+1);
%
%plus entrepreneur's BC -eqn 34
%
CE+Q*K+RLex*L(-1)/(1+Stochg)+TE=Rk*K(-1)*Q(-1)/(1+Stochg)+L;
phi=Q*K/(Q*K-L);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Resource Constraint (Output Equilibrium) -eqn 35
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
Y=C+CE+G+I;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%New variables: LAMBDAE LAMBDACE CE DDE RL DDERL DDERk m L TE RLex ThetaE
%New parameters: habE  BettaE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spread=Rk-Rex;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Government Balanced Budget
%%%%%%%%%%%%%%%%%%%%%%%%%%% 
G=h*WP*tax;
%%%%%%%%%%%%%%%%%%%
%Price Setting
%%%%%%%%%%%%%%%%%%%
H-xi*betta*Htilde(+1)*(1+Stochg(+1))^((1-varrho)*(1-sigma_c))=Y*LAMBDAC;
J-xi*betta*Jtilde(+1)*(1+Stochg(+1))^((1-varrho)*(1-sigma_c))=(1/(1-(1/zzeta)))*Y*LAMBDAC*MC*MS;
PIETILDE=PIE/PIE(-1)^gamp;
Htilde=(PIETILDE^(zzeta-1))*H;
Jtilde=PIETILDE^zzeta*J;
%
%Eqn 45-47
%
1=xi*(PIETILDE^(zzeta-1))+(1-xi)*((J/H)^(1-zzeta));
Delta=xi*(PIETILDE^zzeta)*Delta(-1)+(1-xi)*(J/H)^(-zzeta);
MC=PWP;

%%%%%%%%%%%%%%%
%%Taylor rule%% eqn 48
%%%%%%%%%%%%%%%
%
%Implementable Rule
%
//log(Rn/(STEADY_STATE(Rn)))=alpha_r*log((Rn(-1))/(STEADY_STATE(Rn)))
//+alpha_pie*(1-alpha_r)*log((PIE)/(STEADY_STATE(PIE)))+alpha_y*(1-alpha_r)*log(Y/STEADY_STATE(Y))+epsMPS;
%
%Conventional Taylor Rule
%
log((Rn)/(STEADY_STATE(Rn)))=alpha_r*log((Rn(-1))/(STEADY_STATE(Rn)))
+alpha_pie*(1-alpha_r)*log((PIE)/(STEADY_STATE(PIE)))+
alpha_y*(1-alpha_r)*log((Y/STEADY_STATE(Y))/(YF/STEADY_STATE(YF)))+log(MPS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FLEXI-PRICE, NO FF MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Euler Eqn, Single period utility and Marginal utilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LAMBDAF=((((CF-hab*CF(-1)/(1+Stochg))^(1-varrho))*((1-hF)^varrho))^(1-sigma_c)-1)/(1-sigma_c);
LAMBDACF=(1-varrho)*((CF-hab*CF(-1)/(1+Stochg))^((1-varrho)*(1-sigma_c)-1))*((1-hF)^(varrho*(1-sigma_c)));
%
%Stochastic Discount Factor Eqn 50
%
#DDF=betta*(LAMBDACF(+1)/LAMBDACF)*(1+Stochg(+1))^((1-varrho)*(1-sigma_c)-1);
%
% Euler Eqn 51
%
DDF*RF=1;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Household FOC for labour Supply
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varrho*((CF-hab*CF(-1)/(1+Stochg))^((1-varrho)*(1-sigma_c)))*((1-hF)^(varrho*(1-sigma_c)-1))/LAMBDACF=WPF;
%%%%%%%%%%%%%%%%
%Retail Output eq 54
%%%%%%%%%%%%%%%%
YF=(1-c)*YWF;
%%%%%%%%%%%%%%%%%%%%%%%%%
%Wholesale Firm FOC eqn 55
%%%%%%%%%%%%%%%%%%%%%%%%% 
YWF=((A*hF)^alp)*(KF(-1)/(1+Stochg))^(1-alp);
PWPF*(alp*YWF)/hF=WPF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Resource Constraint (Output Equilibrium)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GF=YF*G/Y;
YF=CF+GF+IF;
%%%%%%%%%%%%%%%%%%%%
%Capital Producers eq 59
%%%%%%%%%%%%%%%%%%%
XF=IF/IF(-1)*(1+Stochg);
KF=capqual(+1)*((IF*(1-phiX*(XF-1-Stochg)^2.0)+(1-delta)*KF(-1)/(1+Stochg)));
#DDFL=betta*(LAMBDACF/LAMBDACF(-1))*(1+Stochg)^((1-varrho)*(1-sigma_c)-1);
Z1F=2.0*phiX*(XF-1-Stochg)*XF^2*QF*DDFL;
QF*(1- phiX*(XF-1-Stochg)^2.0-2.0*XF*phiX*(XF-1-Stochg))+Z1F(+1)=1;
Z2F=(1-alp)*PWPF*YWF/(KF(-1)/(1+Stochg))+(1-delta)*QF;
RkF=capqual*Z2F/QF(-1);
DDF*RF=DDF*RkF(+1);
spreadF=RkF-RF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Government Balanced Budget 67
%%%%%%%%%%%%%%%%%%%%%%%%%%% 
GF=hF*WPF*taxF;
%%%%%%%%%%%%%%%%%%%
%Price Setting
%%%%%%%%%%%%%%%%%%%
PWPF=MCF;
MCF=(1-1/zzeta)/MS;
OUTGAP=YF/Y;

%%%%%%%%%%%%%%%%%%%
%%Shock processes%%
%%%%%%%%%%%%%%%%%%%
log(A)-log(STEADY_STATE(A))=rhoA*(log(A(-1))-log(STEADY_STATE(A)))-epsA;
log(G)-log(STEADY_STATE(G))=rhoG*(log(G(-1))-log(STEADY_STATE(G)))-epsG;
log(MS)=rhoMS*log(MS(-1))+epsMS;
log(MPS)=rhoMPS*log(MPS(-1))+epsMPS;
log(capqual)=rhocapqual*log(capqual(-1))-epscapqual;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Variables in deviation form%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
YY=Y/STEADY_STATE(Y);
KK=K/STEADY_STATE(K);
II=I/STEADY_STATE(I);
//CC=C/STEADY_STATE(C);
CC=(C+CE)/(STEADY_STATE(C)+STEADY_STATE(CE));
CECE=CE/STEADY_STATE(CE);
WPWP=WP/STEADY_STATE(WP);
hh=h/STEADY_STATE(h);
RR=(Rex)/(STEADY_STATE(Rex));
RnRn=(Rn)/(STEADY_STATE(Rn));
RkRk=(Rk)/(STEADY_STATE(Rk));
PIEPIE=PIE/(STEADY_STATE(PIE));
QQ=Q/(STEADY_STATE(QQ));
//
//
YYF=YF/STEADY_STATE(YF);
KKF=KF/STEADY_STATE(KF);
IIF=IF/STEADY_STATE(IF);
CCF=CF/STEADY_STATE(CF);
WPWPF=WPF/STEADY_STATE(WPF);
hhF=hF/STEADY_STATE(hF);
RRF=(RF)/(STEADY_STATE(RF));
QQF=QF/(STEADY_STATE(QQF));
OUTGAPOUTGAP=OUTGAP/(STEADY_STATE(OUTGAP));
OMEGAOMEGA=-OMEGA/(STEADY_STATE(OMEGA));
spreadspread=RkRk(+1)-RR(+1);
phiphi=phi/(STEADY_STATE(phi));

%%%%%%%%%%%%%%%%%%%%%%%%
%%Measurment equations%%
%%%%%%%%%%%%%%%%%%%%%%%%
// In the latest Dynare, one could define PIE(conspie)
dy=log(YY)-log(YY(-1))+trend+epsAtrend;
pinfobs = log(PIEPIE)+conspie; //PIEss-1
robs = log(RnRn)+consr;
rkn_obs = log(RkRk*PIEPIE)+consrkn;
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%SPECIFICATION OF SHOCKS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
shocks;
%
%for irfs
%
var epsA; stderr 1;
var epsG; stderr 1;
var epsMS; stderr 1;
var epsMPS; stderr 1;
var epsAtrend; stderr 0.1;
var epscapqual; stderr 1.0;
end;
/*
steady;
check;

stoch_simul(periods=0,order=1,irf=40) QQ PIEPIE RnRn YY CC CECE II hh WPWP RR OUTGAPOUTGAP 
OMEGAOMEGA spreadspread phiphi;

*/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%BAYESIAN ESTIMATION STARTS HERE%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estimated_params;
// PARAM NAME, INITVAL, LB, UB, PRIOR_SHAPE, PRIOR_P1, PRIOR_P2, PRIOR_P3, PRIOR_P4, JSCALE
// PRIOR_SHAPE: BETA_PDF, GAMMA_PDF, NORMAL_PDF, INV_GAMMA_PDF
// Priors used as in SW07
//
stderr epsA, INV_GAMMA_PDF,0.1,2;      //technology
stderr epsG, INV_GAMMA_PDF,0.5,2;      //government spending - rescaled stderr/gy
stderr epsMPS,INV_GAMMA_PDF,0.1,2;      //interest rate rule
stderr epsMS, INV_GAMMA_PDF,0.1,2;      //mark-up
stderr epsAtrend,INV_GAMMA_PDF,0.1,2;      //trend shock
stderr epscapqual,INV_GAMMA_PDF,0.1,2;      // capital quality shock



rhoA,     BETA_PDF,     0.5,0.20;   //AR1 technology
rhoG,     BETA_PDF,     0.5,0.20;   //AR1 government spending
rhoMS,    BETA_PDF,     0.5,0.20;   //AR1 mark-up
rhocapqual, BETA_PDF,     0.5,0.20;   //AR1 capital quality


phiX,    NORMAL_PDF,   2.00, 1.50; //investment adj cost
sigma_c,  NORMAL_PDF,   1.50, 0.375;//consumption utility
hab,     BETA_PDF,     0.70, 0.10; //habit
xi,        BETA_PDF,     0.50, 0.10; //calvo prices
gamp,      BETA_PDF,     0.50, 0.15; //indexation prices

alpha_pie,  NORMAL_PDF,   2.00, 0.25; //feedback inflation
alpha_r,   BETA_PDF,     0.75, 0.10; //lagged interest rate
alpha_y,     NORMAL_PDF,   0.125,0.05; //feedback output gap 


end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF OBSERVABLE VARIABLES%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
varobs dy pinfobs robs rkn_obs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%sample periods 84:1-08:2
%%4 quarters for initialisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
options_.plot_priors=1;
//estimation 1st stage - posterior mode computation
estimation(datafile=us_data,mode_compute=4,first_obs=5,presample=4, 
              prefilter=0,mh_replic=0,mh_nblocks=0,mh_jscale=0.40,mh_drop=0.2);

