/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Dynare set-up:  GK Banking Model
//
//with unconventional monetary policy
//
//
//Corrected, stationarized H and J price dynamics
//
//
//Price Dynamics Corrected
//
//Built from NK Course Model Revised  with Calvo Contracts, capital and costs of investment and  non-distortionary taxes
// With external habit and price indexing
//Exogenous non-zero growth and inflation
//Lagged variables defined for use with ACES
//Calibrated to target hours worked
//
//Banking Model has an expected and steady state spread
//
//Calibrated to hit targets for spread and leverage
//
//Add Flexi-Price Sector and Taylor Rule with Output gap
//varrho calibrated to target hss in NK model
//
//Note ss of GK model with non-zero inflation is not the same as the flexi-price
//Add trend shock epsAtrend so g becomes Stochg (version 2)
//
//Add negative capital quality shock (version 3) - enters at the end of period t 
// 
//
//
//Remove lags used in Aces-Dynare
// 
//Correction to habit lagged term with growth and growth in  SDF should be dated t+1
// 
//
//Correction to discounting by price setters; betta goes to SDF
// 
//Correction to SDF with growth timed at +1
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
S NW phiB nuB muB omega Rk Z Dep DD DDL spread ThetaB xiB  Delta
NWNW spreadspread  MS MPS MRS PIETILDE PIEPIE varrho
LAMBDAF LAMBDACF RF hF WPF YWF YF  PWPF KF IF  taxF  CF  XF QF 
Z1F Z2F MCF  RkF spreadF  OUTGAP 
RRF YYF CCF hhF WPWPF IIF KKF  QQF OUTGAPOUTGAP OMEGAOMEGA
GF Stochg capqual RkRk phiphi betta bettag OMEGA OMEGAtrue LAMBDAtrue CE GAMMA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF OBSERVABLE VARIABLES%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var dy pinfobs robs rkn_obs; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF EXOGENOUS VARIABLES%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varexo epsA epsG epsMS epsMPS epsMRS epsAtrend epscapqual mes_y mes_pie mes_r mes_rkn;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF PARAMETERS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters gy  alp c zzeta delta sigma_c  
rhoA rhoG rhoMS rhoMPS rhoMRS Ass phiX xi 
sigmaB lev creditspread hab g gamp PIEss hss rhocapqual Rnss;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%ADD OBSERVATION TRENDS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters trend conspie consr consrkn;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%optimized parameters%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters wr tau alpha_r alpha_pie alpha_y alpha_gam alpha_gam_spread;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PARAMETERS%% set to agree with Gertler-Karadi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//PIEss=1.0063;
PIEss=1.0;
g=0.0046;
//g=0;
Rnss=1.013142;
gamp=0.2;
gy=0.2;
hss=0.35;
//hss=0.40;
alp=0.70;
zzeta=7.0;
c=1/zzeta;
delta=0.0250;
sigma_c=2.0;
hab=0.7;
phiX=2.0;
xi=0.75;
///////////////////
//Flexi-Price Case
////////////////////
//xi=0.00001;
/////////////////////////////
//////////////////////////////////
//MP rule
alpha_r=0.7;
alpha_pie=2.0;
alpha_y=0.5;
//
wr=0.0;
//
//Rule in form suitable for optimzation
//
alpha_r=0.7;
alpha_r=0.0;
alpha_pie=2.0*(1-alpha_r);
alpha_y=0.5*(1-alpha_r);
alpha_gam=0.0;
alpha_gam_spread=0.0;
%
%Efficiency Cost Parameter
%
tau=0.1;
tau=0.01;
//tau=0.001;
//
//////////////////////////
//Indeterminacy Analysis
//////////////////////////
//
//alpha_r=parameter2(i);
//alpha_pie=parameter1(j);
//
//
//PIEss=parameter2(i);
//alpha_pie=parameter1(j);
//
/////////////////////////////
%
%
%

////////////////////////////
//Choice of Units
////////////////////////////
Ass=1;
//////////////////////////
//shock persistence
//////////////////////////
rhoA=0.75;
rhoG=0.75;
rhoMS=0.75;
rhoMPS=0.75;
rhoMRS=0.75;
rhocapqual=0.75;
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
//creditspread=0.0088125/PIEss;
creditspread=0.01/4;
//creditspread=0.0;
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
//
//--------------------------------------------------------------
//define auxiliary variables for instrument and lagged variables
//Only Necessary for ACES-DYNARE
//
//--------------------------------------------------------------
//

%%%%%%%%%%%%%%%%%%%%%%%%%%
%calibration of varrho (Eqn 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
varrho=STEADY_STATE(varrho);
betta=STEADY_STATE(betta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Stochastic Trend
%
log(1+Stochg)=log(1+g)+epsAtrend;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Single-Period and intertemporal welfare
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bettag=betta*(1+Stochg)^((1-varrho)*(1-sigma_c));
LAMBDA=((((1+Stochg)*C-hab*C(-1))^(1-varrho)*(1-h)^varrho)^(1-sigma_c))/(1-sigma_c) -wr*(Rn-STEADY_STATE(Rn))^2;
LAMBDAtrue=((((1+Stochg)*C-hab*C(-1))^(1-varrho)*(1-h)^varrho)^(1-sigma_c))/(1-sigma_c);
%
%
%CE Calculation
%
CE=((((1+Stochg)*C*1.01-hab*C(-1)*1.01)^(1-varrho)*(1-h)^varrho)^(1-sigma_c))/(1-sigma_c)-LAMBDAtrue;
%
OMEGA=(1-bettag)*LAMBDA+bettag*OMEGA(+1);
OMEGAtrue=(1-bettag)*LAMBDAtrue+bettag*OMEGAtrue(+1);
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Marginal utility of consumption%% Eqn 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LAMBDAC=(1-varrho)*((C-hab*C(-1)/(1+Stochg))^((1-varrho)*(1-sigma_c)-1))*((1-h)^(varrho*(1-sigma_c)));

%%%%%%%%%%%%%%%%%%
%%Euler equation%%
%%%%%%%%%%%%%%%%%%
LAMBDAC=betta*(1+Stochg(+1))^((1-varrho)*(1-sigma_c)-1)*XX(+1);
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
%%Price Dispersion%%%% Eqn 9
%%%%%%%%%%%%%%%%%%%%%%
Delta=xi*(PIETILDE^zzeta)*Delta(-1)+(1-xi)*(J/H)^(-zzeta);

%%%%%%%%%%%%%%%%%%%%%%%
%%Production Function%% Eqn 10
%%%%%%%%%%%%%%%%%%%%%%%
YW=((A*h)^alp)*(K(-1)/(1+Stochg))^(1-alp)/Delta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Wholesale firms FOC for labour%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PWP*(alp*YW)/h=WP;

%%%%%%%%%%%%%%%%%%%%%%%
%%Resource constraint%% (Eqn 12)
%%%%%%%%%%%%%%%%%%%%%%%
Y=C+G+I+tau*GAMMA*Q*K;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Capital law of motion%%%%%
%%Costs of investment case%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=I/I(-1)*(1+Stochg);
//
// add capital quality shock
//
K=capqual(+1)*((I*(1-phiX*(X-1-Stochg)^2.0)+(1-delta)*K(-1)/(1+Stochg)));
//K=capqual*((I*(1-phiX*(X-1-Stochg)^2.0)+(1-delta)*K(-1)/(1+Stochg)));


%%%%%%%%%%%%%%
%%Investment%%
%%%%%%%%%%%%%%
Z1=2.0*phiX*(X-1-Stochg)*X^2*Q*DDL;
Q*(1- phiX*(X-1-Stochg)^2.0-2.0*X*phiX*(X-1-Stochg))+Z1(+1)=1;

%%%%%%%%%%%%%%%
%%Fischer Eqn%%
%%%%%%%%%%%%%%%
INVPIE=1/PIE;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Real ex-post interest rate%% Eqn 34
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rex=(Rn(-1))*INVPIE;

%%%%%%%%%%%%%
%%Tobin's Q%%
%%%%%%%%%%%%%
//Z2=(1-alp)*PWP*YW/K(-1)+(1-delta)*Q;
//(Rn)*INVPIE(+1)=Z2(+1)/Q;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Balance budget constraint%%% Eqn 35
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G=h*WP*tax; 

%%%%%%%%%%%%%%%%%%%%%%
%%Inflation Dynamics%%
%%%%%%%%%%%%%%%%%%%%%%
//DDHtilde=DD(-1)*Htilde;
//DDJtilde=DD(-1)*Jtilde;
H-xi*betta*Htilde(+1)*(1+Stochg(+1))^((1-varrho)*(1-sigma_c))=Y*LAMBDAC;
J-xi*betta*Jtilde(+1)*(1+Stochg(+1))^((1-varrho)*(1-sigma_c))=(1/(1-(1/zzeta)))*Y*LAMBDAC*MC*MS;
PIETILDE=PIE/PIE(-1)^gamp;
Htilde=(PIETILDE^(zzeta-1))*H;
Jtilde=PIETILDE^zzeta*J;
1=xi*(PIETILDE^(zzeta-1))+(1-xi)*((J/H)^(1-zzeta));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Mark-up Monopolistic pricing%% Eqn 28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
MC=PWP;

%%%%%%%%%%%%%%%%%%
%%Banking Sector%%
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
//S=K/capqual(+1);   
Q*S=phiB*NW+GAMMA*Q*S;
phiB=nuB/(ThetaB-muB);
NW=(sigmaB+xiB)*(Z+(1-delta)*Q)*S(-1)*capqual/(1+Stochg)-(Rex)*sigmaB*Dep(-1)/(1+Stochg);
//NW=(sigmaB+xiB)*(Rk)*QL*S(-1)-(Rex)*sigmaB*DepL;
Dep=Q*S-NW;
DD=betta*(LAMBDAC(+1)/LAMBDAC)*(1+Stochg(+1))^((1-varrho)*(1-sigma_c)-1);
DDL=betta*(LAMBDAC/LAMBDAC(-1))*(1+Stochg)^((1-varrho)*(1-sigma_c)-1);
nuB=DD*omega(+1)*Rex(+1);
spread=Rk-Rex;
//spread=Rk(+1)-Rex(+1);
muB=DD*omega(+1)*spread(+1);
omega=1-sigmaB+sigmaB*ThetaB*phiB;

Z=(1-alp)*PWP*YW/(K(-1)/(1+Stochg));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Eq 42
%Needs rescaling for zero credit-spread  case
% 
//0.1*Rk=0.1*(Z+(1-delta)*Q)/Q(-1);
%
//
// add capital quality shock
//
//0.1*Rk=0.1*capqual*((Z+(1-delta)*Q)/Q(-1));
Rk=capqual*((Z+(1-delta)*Q)/Q(-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%
%%Taylor rule%%
%%%%%%%%%%%%%%%
%
%Implementable Rule. Not suitable for optimized rule
%log((Rn)/(STEADY_STATE(Rn)))=alpha_r*log(Rn(-1)/STEADY_STATE(Rn))
%+alpha_pie*(1-alpha_r)*log(PIE/STEADY_STATE(PIE))
%+alpha_y*(1-alpha_r)*log(Y/STEADY_STATE(Y))+epsMPS;
%
%Conventional Taylor Rule. Not suitable for optimized rule
%

%log((Rn)/(STEADY_STATE(Rn)))=alpha_r*log((Rn(-1))/(STEADY_STATE(Rn)))
%+alpha_pie*(1-alpha_r)*log((PIE)/(STEADY_STATE(PIE)))
%+alpha_y*(1-alpha_r)*log((Y/STEADY_STATE(Y))/(YF/STEADY_STATE(YF)))+log(MPS);
%
%Implementable Rule. Suitable for optimized rule
%
log((Rn)/(STEADY_STATE(Rn)))=alpha_r*log(Rn(-1)/STEADY_STATE(Rn))
+alpha_pie*log(PIE/STEADY_STATE(PIE))
+alpha_y*log(Y/STEADY_STATE(Y))+epsMPS;
%
%
%Conventional Taylor Rule in form suitable for optimized rule
%
%log((Rn)/(STEADY_STATE(Rn)))=alpha_r*log((Rn(-1))/(STEADY_STATE(Rn)))
%+alpha_pie*log((PIE)/(STEADY_STATE(PIE)))
%+alpha_y*log((Y/STEADY_STATE(Y))/(YF/STEADY_STATE(YF)))+log(MPS);
%
%Unconventional Policy
%
(GAMMA-STEADY_STATE(GAMMA))=alpha_gam* (GAMMA(-1)-STEADY_STATE(GAMMA))+
50*alpha_gam_spread*(spreadspread-STEADY_STATE(spreadspread));
%
%Only Conventional Policy Case
%
%GAMMA=STEADY_STATE(GAMMA);
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
//Z1F=2.0*phiX*(XF-1-Stochg)*XF^2*QF/RF(-1);
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
log(MRS)=rhoMRS*log(MRS(-1))+epsMRS;
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
//spreadspread=spread/STEADY_STATE(spread);
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
OMEGAOMEGA=-OMEGA/(STEADY_STATE(OMEGA));

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
//var epsA; stderr 1.0;
//var epsG; stderr 1.0;
//var epsMS; stderr 1.0;
//var epsMPS; stderr 1.0;
//var epsAtrend; stderr 0.1;
//var epscapqual; stderr 1.0;
%
%for second moment mean calculations
% 
//var epsA; stderr 1;//to check SO solutions (negative mean C)
var epsA; stderr 0.01;
var epsG; stderr 0.01;
var epsMS; stderr 0.01;
var epsMPS; stderr 0.01;
var epsMRS; stderr 0.01;
var epsAtrend; stderr 0.01;
var epscapqual; stderr 0.01;
end;

steady;

check;



//
//comment out check 
//
//irfs 
//
//Compute welfare-optimized rules in dynare
//'stoch_simul' not required
//'steady not' required (without steady it produces zero deterministic welfare)
//Needs optwelrule.m, stoch_simul_wel.m, disp_th_moments_wel.m in the working dir.
//NOTE: when declaring parameters put the feedbacks to the last: 
//e.g. parameters a, b, c, ..., alpha_r alpha_pie alpha_y;
//For optimising it takes their assigned values as the starting values
//May have problem with indeter/instab: the usual error is that the iteration won't
//start if the starting values cause indeter - I'll need to think about this
//This command optwelrule trigers the complutation and returns WELFARE FUNCTION and Opt rule
//Arguments explained:
//   strvcat('OMEGA') - tells it the variable name of inter-temp welfare
//   3 - the number of feedbacks chosen to be optimised (takes e.g. 3 2 or 1 in this model)
//       e.g. if 1 is chosen it only optimises alpha_y (the last parameter)
//   [] - for lower and upper bounds respectively (being empty takes -inf and inf)
//----------------------------------------------------------------------------
//
optwelrule(strvcat('OMEGA'),5,[0 0 0 0 0],[1 5 5 1 1]);
//stoch_simul(order=1,irf=40) QQ PIEPIE RnRn YY CC II hh WPWP RR 
//OUTGAPOUTGAP OMEGAOMEGA spreadspread phiphi NWNW;
//
//No irfs - theoretical moments
//
stoch_simul(order=2, irf=0) Q PIE Rn Y C I h WP Rex OUTGAP OMEGA OMEGAtrue;
//



