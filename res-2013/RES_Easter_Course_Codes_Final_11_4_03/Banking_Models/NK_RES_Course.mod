/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Dynare set-up:  NK Core Model with Calvo Contracts, capital and costs of investment and  non-distortionary taxes
//
//With external habit and price indexing
//Exogenous non-zero growth and inflation
//
//varrho calibrated to target steady state hours, hss
//betta calibrated to target steady state nominal interest rate, Rnss 

//
//with Rotemberg contract versus Calvo options
//Flexi-Price Sector and Taylor Rule with Output gap
//
//Five exogenous stochastic processes: temporary technology (A), gov spending (G), 
//price mark-up (MS), monetary policy shock (MPS), stochastic trend and capital quality shock
//Note ss of NK model with non-zero inflation is not the same as the flexi-price
//
//Intertemporal welfare OMEGA in Bellman equation form for zero growth case
// 
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
Z1 Z2 MC  Rn PIE INVPIE Rk spread RR RnRn YY CC hh WPWP II KK   MS MPS
  PIEPIE QQ  DD DDL DDF DDFL LAMBDAF LAMBDACF RF hF WPF YWF YF  PWPF KF IF  taxF  CF  XF QF 
Z1F Z2F MCF  RkF spreadF  OUTGAP RRF YYF CCF hhF WPWPF IIF KKF  QQF OUTGAPOUTGAP OMEGAOMEGA
 GF Stochg capqual  RkRk spreadspread phiphi NWNW  Delta H Htilde J Jtilde PIETILDE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF CALIBRATED PARAMETERS AS ENDOGENOUS VARIABLES%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var   betta varrho;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF OBSERVABLE VARIABLES%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var dy pinfobs robs; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF EXOGENOUS VARIABLES%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varexo epsA epsG epsMS epsMPS epsAtrend epscapqual mes_y mes_pie mes_r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF PARAMETERS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters gy  alp c zzeta  delta sigma_c rhoA rhoG rhoMS rhoMPS
Ass  phiX xi theta_y hab  hss g gamp PIEss Rnss
alpha_r $\alpha_r$
theta_pi $\theta_{\pi}$ rhocapqual ppsi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%ADD OBSERVATION TRENDS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters trend conspie consr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PARAMETERS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Set varrho, betta as exogenous parameters
%
%varrho=0.924456;
%betta=0.998162;
//
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
//////////////////////////
//Indeterminacy Analysis - remove // and comment out stoch_simul at end
//////////////////////////
//
//alpha_r=parameter2(i);
//theta_pi=parameter1(j);
//
// or
//
//PIEss=parameter2(i);
//theta_pi=parameter1(j);
//
//
//Link price adjustment costs in Rotemberg with Calvo; for zero inflation and growth
//
ppsi=(zzeta-1)*xi/((1-xi)*(1-betta*xi));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define parameters for observation equqtions
%(only needed with esimation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trend=g*100;
consr=(Rnss-1)*100; % ss nominal interest rate
conspie=100*(PIEss-1); % quarterly ss inflation rate
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// ----------------------------
// *** DSGE-Model-equations ***
// ----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model;

%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%calibration of varrho and betta - now out
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
LAMBDA=((((C-hab*C(-1)/(1+Stochg))^(1-varrho))*((1-h)^varrho))^(1-sigma_c)-1)/(1-sigma_c);
LAMBDAC=(1-varrho)*((C-hab*C(-1)/(1+Stochg))^((1-varrho)*(1-sigma_c)-1))*((1-h)^(varrho*(1-sigma_c)));
%
%
XX=(Rn(-1)/PIE)*LAMBDAC;
LAMBDAC=betta*(1+Stochg(+1))^((1-varrho)*(1-sigma_c)-1)*XX(+1);

%Stochastic Discount Factor
%
DD=betta*(LAMBDAC(+1)/LAMBDAC)*(1+Stochg(+1))^((1-varrho)*(1-sigma_c)-1);
DDL=betta*(LAMBDAC/LAMBDAC(-1))*(1+Stochg)^((1-varrho)*(1-sigma_c)-1);
%
%Alternative Euler eqn formulation
%
//1=DD*Rn/PIE(+1);


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Intertemporal welfare (zero ss growth case only)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OMEGA=(1-betta)*LAMBDA+betta*OMEGA(+1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Resource Constraint (Output Equilibrium)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y=C+G+I;
%%%%%%%%%%%%%%%%%%%%
%Capital Producers 
%%%%%%%%%%%%%%%%%%%
X=I/I(-1)*(1+Stochg);
K=capqual(+1)*((I*(1-phiX*(X-1-Stochg)^2.0)+(1-delta)*K(-1)/(1+Stochg)));
Z1=2.0*phiX*(X-1-Stochg)*X^2*Q*DDL;
Q*(1- phiX*(X-1-Stochg)^2.0-2.0*X*phiX*(X-1-Stochg))+Z1(+1)=1;
INVPIE=1/PIE;
Rex=Rn(-1)*INVPIE;
Z2=(1-alp)*PWP*YW/(K(-1)/(1+Stochg))+(1-delta)*Q;
Rk=capqual*Z2/Q(-1);
DD*Rn*INVPIE(+1)=DD*Rk(+1);
spread=Rk(+1)-Rex(+1);
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
%
%
1=xi*(PIETILDE^(zzeta-1))+(1-xi)*((J/H)^(1-zzeta));
Delta=xi*(PIETILDE^zzeta)*Delta(-1)+(1-xi)*(J/H)^(-zzeta);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Price Setting -Rotemberg Contracts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Replace Calvo price-setting above
%
//H-xi*betta*Htilde(+1)=Y*LAMBDAC;
//J-xi*betta*Jtilde(+1)=(1/(1-(1/zzeta)))*Y*LAMBDAC*MC*MS;
//PIETILDE=PIE/PIE(-1)^gamp;
//Htilde=(PIETILDE^(zzeta-1))*H;
//Jtilde=PIETILDE^zzeta*J;
//1=xi*(PIETILDE^(zzeta-1))+(1-xi)*((J/H)^(1-zzeta));
//Delta=xi*(PIETILDE^zzeta)*Delta(-1)+(1-xi)*(J/H)^(-zzeta);
%
%with 
%
//1-ppsi*(PIE-1)*PIE+ppsi*DD*((PIE(+1)-1)*PIE(+1)*(Y(+1)/Y)*(1+Stochg(+1)))=(1-MC*MS)*zzeta;
//Delta=1;
//
% and remove H, J, PIETILDE, Htilde, Jtilde from variable list
% 
%In steady state and fun files replace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Delta=((1-xi)*(((1-xi*PIETILDE^(zzeta-1))/(1-xi))^(1/(1-zzeta)))^(-zzeta))/(1-xi*PIETILDE^zzeta);
%MC=(1-1/zzeta)*(1-xi*betta*PIETILDE^zzeta)/(1-xi*betta*PIETILDE^(zzeta-1))*(((1-xi*PIETILDE^(zzeta-1))/(1-xi))^(1/(1-zzeta)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%with
%Delta=1;
%MC=1-(1-ppsi*(PIE-1)*PIE*(1-DD*(1+g)))/zzeta;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Continue with NK Model (Calvo or Rotemberg)

MC=PWP;

%%%%%%%%%%%%%%%
%%Taylor rule%% 
%%%%%%%%%%%%%%%
%
%Implementable Rule
%
//log((Rn)/(STEADY_STATE(Rn)))=alpha_r*log((Rn(-1))/(STEADY_STATE(Rn)))
//+theta_pi*(1-alpha_r)*log((PIE)/(STEADY_STATE(PIE)))+theta_y*(1-alpha_r)*log(Y/STEADY_STATE(Y))+epsMPS;
%
%Conventional Taylor Rule
%
log((Rn)/(STEADY_STATE(Rn)))=alpha_r*log((Rn(-1))/(STEADY_STATE(Rn)))
+theta_pi*(1-alpha_r)*log((PIE)/(STEADY_STATE(PIE)))
+theta_y*(1-alpha_r)*log((Y/STEADY_STATE(Y))/(YF/STEADY_STATE(YF)))+log(MPS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FLEXI-PRICE MODEL
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
DDF=betta*(LAMBDACF(+1)/LAMBDACF)*(1+Stochg(+1))^((1-varrho)*(1-sigma_c)-1);
%
% F Euler Eqn 
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
DDFL=betta*(LAMBDACF/LAMBDACF(-1))*(1+Stochg)^((1-varrho)*(1-sigma_c)-1);
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
RR=(Rex)/(STEADY_STATE(Rex));
RnRn=(Rn)/(STEADY_STATE(Rn));
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
RkRk=Rk/(STEADY_STATE(Rk));
spreadspread=RkRk(+1)-RR(+1);
//
//leverage and net worth - to compare irfs with banking models
//
phiphi=1;
NWNW=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Measurment equations - for estimation%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy=log(YY)-log(YY(-1))+trend+epsAtrend+mes_y;
pinfobs = log(PIEPIE)+conspie+mes_pie; //PIEss-1
robs = log(RnRn)+consr+mes_r;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%SPECIFICATION OF SHOCKS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
shocks;
%
%for irfs and first order solution only
%
var epsA; stderr 1;
var epsG; stderr 1;
var epsMS; stderr 1;
var epsMPS; stderr 1;
var epsAtrend; stderr 1;
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
//
//comment out simul for stability analysis
//
//irfs 
//
stoch_simul(order=1,irf=40) QQ PIEPIE RnRn YY CC II hh WPWP RR 
OUTGAPOUTGAP OMEGAOMEGA spreadspread phiphi NWNW;
//
//No irfs - theoretical moments
//
//stoch_simul(order=2, irf=0) Q PIE Rn Y C I h WP Rex OUTGAP OMEGA;
//


