/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Dynare set-up:  GK Banking Model
//
//Add banking sector to NK_RES_Course with Calvo Contracts, capital and costs of investment and  non-distortionary taxes
//With external habit and price indexing
//Exogenous non-zero growth and inflation
//varrho calibrated to target steady state hours, hss
//betta calibrated to target steady state nominal interest rate, Rnss 
//
//Five exogenous stochastic processes: temporary technology (A), gov spending (G), 
//price mark-up (MS), monetary policy shock (MPS), stochastic trend and capital quality shock

//
//
//Banking sector calibrated to target spread and leverage
//
//
//Note ss of BGG model with is not the same as the flexi-price
//
//(c) CIMS Univeristy of Surrey
//The Science and  Art of DSGE Modelling: Construction, Calibration, Estimation and Policy
//////////////////////////////////////////////////////////////////////////////////////////////////
//This is a free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.  See <http://www.gnu.org/licenses/> for more information.
////////////////////////////////////////////////////////////////////////////////////////////////// 
//
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF ENDOGENOUS VARIABLES%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var LAMBDA LAMBDAC Rex  h WP YW Y PWP K I  tax  C A  G X Q XX 
Z1 MC H Htilde J  Jtilde PIE Rn INVPIE RR RnRn YY CC hh WPWP II KK QQ
S NW phiB nuB muB omega Rk Z Dep DD DDL spread  Delta
NWNW spreadspread  MS MPS PIETILDE PIEPIE 
LAMBDAF LAMBDACF RF hF WPF YWF YF  PWPF KF IF  taxF  CF  XF QF 
Z1F Z2F MCF  RkF spreadF  OUTGAP 
RRF YYF CCF hhF WPWPF IIF KKF  QQF OUTGAPOUTGAP
GF Stochg capqual RkRk phiphi varrho betta ThetaB xiB;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF OBSERVABLE VARIABLES%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var dy pinfobs robs rkn_obs; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF EXOGENOUS VARIABLES%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varexo epsA epsG epsMS epsMPS epsAtrend epscapqual mes_y mes_pie mes_r mes_rkn;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF PARAMETERS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters gy  alp c zzeta delta sigma_c  
rhoA rhoG rhoMS rhoMPS Ass phiX xi alpha_r theta_pi theta_y 
sigmaB lev creditspread hab g gamp PIEss hss rhocapqual Rnss;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%ADD OBSERVATION TRENDS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters trend conspie consr consrkn;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PARAMETERS%% set to agree with Gertler-Karadi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//
//Common Parameter Values for NK, BGG, GK and KM Simulations Only
//Anticipates Bayesian estimation results
////////////////////////////////////////////////////
//Steady State Inflation, growth and nominal interest rate
//
PIEss=1.0063;
//PIEss=1.0;//zero inflation
g=0.0046;
//g=0;//zero growth
Rnss=1.013142;
//Rnss=1.013142/1.0063;// to calibrate betta for zero inflation
////////////////////////////////////////////////////
//
gy=0.2;
hss=0.35;
alp=0.70;
zzeta=7.0;
c=1/zzeta;
delta=0.0250;
phiX=2.0;
sigma_c=2.0;
hab=0.7;
gamp=0.2;
//gamp=0.0;//no indexing
xi=0.75;
///////////////////
//Flexi-Price Case
////////////////////
//xi=0.00001;
/////////////////////////////
//
//MP rule
//
alpha_r=0.7;
theta_pi=2.0;
theta_y=0.2;
//theta_y=0.0;//inflation targeting alone
//
////////////////////
///shock persistence
//////////////////////
//
rhoA=0.75;
rhoG=0.75;
rhoMS=0.75;
rhoMPS=0.75;
rhocapqual=0.75;
///////////////////////
//Choice of Units
//////////////////////
Ass=1;
//////////////////////
//
//
//////////////////////////
//Indeterminacy Analysis
//////////////////////////
//
//alpha_r=parameter2(i);
//theta_pi=parameter1(j);
//
//
//PIEss=parameter2(i);
//theta_pi=parameter1(j);
//
/////////////////////////////

////////////////////////////
//Banking Sector Parameters 
///////////////////////////
//sigmaB=0.9688;
//sigmaB=0.950;// SP not satisfied
//sigmaB=0.94;
//sigmaB=0.93;
//////////////
//GK Handbook
/////////////
sigmaB=0.975;
lev=4;
creditspread=0.01/4;

//
//


trend=g*100;
consr=(Rnss-1)*100; % ss nominal interest rate
conspie=100*(PIEss-1); % quarterly ss inflation rate
consrkn=2.1955;% ss nominal BAA corporate bonds

////////////////////////////////
//
// ----------------------------
// *** DSGE-Model-equations ***
// ----------------------------
//
model;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%calibration of varrho 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
varrho=STEADY_STATE(varrho);
betta=STEADY_STATE(betta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Stochastic Trend
%
log(1+Stochg)=log(1+g)+epsAtrend;
%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%Single period utility%%
%%%%%%%%%%%%%%%%%%%%%%%%%
LAMBDA=((((C-hab*C(-1)/(1+Stochg))^(1-varrho))*((1-h)^varrho))^(1-sigma_c)-1)/(1-sigma_c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Marginal utility of consumption%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LAMBDAC=(1-varrho)*((C-hab*C(-1)/(1+Stochg))^((1-varrho)*(1-sigma_c)-1))*((1-h)^(varrho*(1-sigma_c)));

%%%%%%%%%%%%%%%%%%
%%Euler equation%%
%%%%%%%%%%%%%%%%%%
LAMBDAC=betta*(1+Stochg)^((1-varrho)*(1-sigma_c)-1)*XX(+1);
XX=Rn(-1)/(PIE)*LAMBDAC;



%%%%%%%%%%%%%%%%%%%%%
%%Labour supply foc%%
%%%%%%%%%%%%%%%%%%%%%
varrho*((C-hab*C(-1)/(1+Stochg))^((1-varrho)*(1-sigma_c)))*((1-h)^(varrho*(1-sigma_c)-1))/LAMBDAC=WP;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Wholesale and retail sector relation%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y=(1-c)*YW;

%%%%%%%%%%%%%%%%%%%%%%
%%Price Dispersion%%%% 
%%%%%%%%%%%%%%%%%%%%%%
Delta=xi*(PIETILDE^zzeta)*Delta(-1)+(1-xi)*(J/H)^(-zzeta);

%%%%%%%%%%%%%%%%%%%%%%%
%%Production Function%% 
%%%%%%%%%%%%%%%%%%%%%%%
YW=((A*h)^alp)*(K(-1)/(1+Stochg))^(1-alp)/Delta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Wholesale firms FOC for labour%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PWP*(alp*YW)/h=WP;

%%%%%%%%%%%%%%%%%%%%%%%
%%Resource constraint%% 
%%%%%%%%%%%%%%%%%%%%%%%
Y=C+G+I;
%
%%%%%%%%%%%%%%%%%%%%
%Capital Producers 
%%%%%%%%%%%%%%%%%%%
X=I/I(-1)*(1+Stochg);
K=capqual(+1)*((I*(1-phiX*(X-1-Stochg)^2.0)+(1-delta)*K(-1)/(1+Stochg)));
Z1=2.0*phiX*(X-1-Stochg)*X^2*Q*DDL;
Q*(1- phiX*(X-1-Stochg)^2.0-2.0*X*phiX*(X-1-Stochg))+Z1(+1)=1;

INVPIE=1/PIE;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Real ex-post interest rate%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rex=(Rn(-1))*INVPIE;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Balance budget constraint%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G=h*WP*tax; 

%%%%%%%%%%%%%%%%%%%%%%
%%Inflation Dynamics%%
%%%%%%%%%%%%%%%%%%%%%%
H-xi*betta*Htilde(+1)*(1+Stochg(+1))^((1-varrho)*(1-sigma_c))=Y*LAMBDAC;
J-xi*betta*Jtilde(+1)*(1+Stochg(+1))^((1-varrho)*(1-sigma_c))=(1/(1-(1/zzeta)))*Y*LAMBDAC*MC*MS;
PIETILDE=PIE/PIE(-1)^gamp;
Htilde=(PIETILDE^(zzeta-1))*H;
Jtilde=PIETILDE^zzeta*J;
1=xi*(PIETILDE^(zzeta-1))+(1-xi)*((J/H)^(1-zzeta));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Mark-up Monopolistic pricing%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
MC=PWP;

%%%%%%%%%%%%%%%%%%
%%GK Banking Sector%%
%%%%%%%%%%%%%%%%%%
%
%Calibration of banking parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ThetaB=STEADY_STATE(ThetaB);
xiB=STEADY_STATE(xiB);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=((I*(1-phiX*(X-1-Stochg)^2.0)+(1-delta)*capqual*S(-1)/(1+Stochg)));
Q*S=phiB*NW;
phiB=nuB/(ThetaB-muB);
NW=(sigmaB+xiB)*(Z+(1-delta)*Q)*S(-1)*capqual/(1+Stochg)-(Rex)*sigmaB*Dep(-1)/(1+Stochg);
Dep=Q*S-NW;
DD=betta*(LAMBDAC(+1)/LAMBDAC)*(1+Stochg(+1))^((1-varrho)*(1-sigma_c)-1);
DDL=betta*(LAMBDAC/LAMBDAC(-1))*(1+Stochg)^((1-varrho)*(1-sigma_c)-1);
%%%%%%%%%%%%%%%%%%%%%%%
%
%nu_{d,t} in notes
%
nuB=DD*omega(+1)*(Rex(+1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%
spread=Rk-Rex;
%
%mu_{s,t} in notes
%
muB=DD*omega(+1)*spread(+1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega=1-sigmaB+sigmaB*ThetaB*phiB;

Z=(1-alp)*PWP*YW/(K(-1)/(1+Stochg));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%Needs rescaling for zero credit-spread  case
% 
//
//0.1*Rk=0.1*capqual*((Z+(1-delta)*Q)/Q(-1));
//
Rk=capqual*((Z+(1-delta)*Q)/Q(-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%
%%Taylor rule%%
%%%%%%%%%%%%%%%
%
%Implementable Rule
%
//log((Rn)/(STEADY_STATE(Rn)))=alpha_r*log(Rn(-1)/STEADY_STATE(Rn))
//+theta_pi*(1-alpha_r)*log(PIE/STEADY_STATE(PIE))
//+theta_y*(1-alpha_r)*log(Y/STEADY_STATE(Y))+epsMPS;
%
%Conventional Taylor Rule
%

log((Rn)/(STEADY_STATE(Rn)))=alpha_r*log((Rn(-1))/(STEADY_STATE(Rn)))
+theta_pi*(1-alpha_r)*log((PIE)/(STEADY_STATE(PIE)))
+theta_y*(1-alpha_r)*log((Y/STEADY_STATE(Y))/(YF/STEADY_STATE(YF)))+log(MPS);
%
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
%Stochastic Discount Factor 
%
#DDF=betta*(LAMBDACF(+1)/LAMBDACF)*(1+Stochg(+1))^((1-varrho)*(1-sigma_c)-1);
%
% Euler Eqn 
%
DDF*RF=1;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Household FOC for labour Supply
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varrho*((CF-hab*CF(-1)/(1+Stochg))^((1-varrho)*(1-sigma_c)))*((1-hF)^(varrho*(1-sigma_c)-1))/LAMBDACF=WPF;
%%%%%%%%%%%%%%%%
%Retail Output 
%%%%%%%%%%%%%%%%
YF=(1-c)*YWF;
%%%%%%%%%%%%%%%%%%%%%%%%%
%Wholesale Firm FOC 
%%%%%%%%%%%%%%%%%%%%%%%%% 
YWF=((A*hF)^alp)*(KF(-1)/(1+Stochg))^(1-alp);
PWPF*(alp*YWF)/hF=WPF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Resource Constraint (Output Equilibrium)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GF=YF*G/Y;
YF=CF+GF+IF;
%%%%%%%%%%%%%%%%%%%%
%Capital Producers 
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
%Government Balanced Budget 
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
CC=C/STEADY_STATE(C);
WPWP=WP/STEADY_STATE(WP);
hh=h/STEADY_STATE(h);
RR=Rex/(STEADY_STATE(Rex));
RkRk=Rk/(STEADY_STATE(Rk));
RnRn=(Rn)/(STEADY_STATE(Rn));
QQ=Q/(STEADY_STATE(QQ));
PIEPIE=PIE/(STEADY_STATE(PIE));
spreadspread=RkRk(+1)-RR(+1);
NWNW=NW/STEADY_STATE(NW);
phiphi=phiB/STEADY_STATE(phiB);
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

%%%%%%%%%%%%%%%%%%%%%%%%
%%Measurment equations%%
%%%%%%%%%%%%%%%%%%%%%%%%
// In the latest Dynare, one could define PIE(conspie)
dy=log(YY)-log(YY(-1))+trend+epsAtrend+mes_y;
pinfobs = log(PIEPIE)+conspie+mes_pie; //PIEss-1
robs = log(RnRn)+consr+mes_r;
rkn_obs = log(RkRk*PIEPIE)+consrkn+mes_rkn;
end;


/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%INITIAL GUESSES FOR STEADY-STATE COMPUTATION%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initval;
end;*/

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%SPECIFICATION OF SHOCKS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
shocks;
var epsA; stderr 1.0;
var epsG; stderr 1.0;
var epsMS; stderr 1.0;
var epsMPS; stderr 1.0;
var epsAtrend; stderr 1.0;
var epscapqual; stderr 1.0;
%
%for second moment mean calculations; common for NK, BGG, GK, KM
% 

//var epsA; stderr 0.01;
//var epsG; stderr 0.01;
//var epsMS; stderr 0.01;
//var epsMPS; stderr 0.01;
//var epsAtrend; stderr 0.01;
//var epscapqual; stderr 0.01;
end;

steady;

check;


stoch_simul(order=1,irf=40) QQ PIEPIE RnRn YY CC  II hh WPWP RR RkRk OUTGAPOUTGAP
spreadspread NWNW phiphi;
//
//No irfs - theoretical moments
//
//stoch_simul(order=2, irf=0) Q PIE Rn Y C I h WP Rex OUTGAP spreadspread;
//
