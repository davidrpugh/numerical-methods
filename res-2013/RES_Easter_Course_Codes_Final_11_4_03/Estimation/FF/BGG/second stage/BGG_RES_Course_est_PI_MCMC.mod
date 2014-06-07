/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Dynare set-up:  BGG Banking Model
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
 NW Rk Z  DD DDL spread sigmaE Apsi  mu
psi rho p fnGam fnG DGam DG rhoRex CE phi  rhorho  psipsi Delta
NWNW spreadspread  MS MPS PIETILDE PIEPIE varrho
LAMBDAF LAMBDACF RF hF WPF YWF YF  PWPF KF IF  taxF  CF  XF QF 
Z1F Z2F MCF  RkF spreadF  OUTGAP 
RRF YYF CCF hhF WPWPF IIF KKF  QQF OUTGAPOUTGAP
 GF Stochg capqual RkRk  S phiphi xiE betta;
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
parameters gy  alp c zzeta delta sigma_c  
rhoA rhoG rhoMS rhoMPS Ass phiX xi alpha_r alpha_pie alpha_y 
 pcalib rhocalib phicalib CEcalib hab g gamp PIEss hss rhocapqual alpha_q  epsilonA Rnss;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%ADD OBSERVATION TRENDS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters trend conspie consr consrkn;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PARAMETERS%% set to agree with Gertler-Karadi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PIEss=1.0063;
g=0.0046;
Rnss=1.013142;
gamp=0.3;
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
//////////////////////////////////
//MP rule
alpha_r=0.7;
alpha_pie=2.0;
alpha_y=0.5;
alpha_q=0.0;

hab=0.7;
////////////////////////////
//Choice of Units
////////////////////////////
Ass=1;
//////////////////////////
//shock persistence
//////////////////////////
rhoA=0.7;
rhoG=0.7;
rhoMS=0.7;
rhoMPS=0;
rhocapqual=0.7;
////////////////////////////
//Banking Sector Parameters
///////////////////////////
epsilonA=0.0;//Remove this feature
pcalib=0.02;
rhocalib=1.005;
phicalib=2.0;
CEcalib=0.1;
//
//
////////////////////////////////
//

trend=g*100;
consr=(Rnss-1)*100; % ss nominal interest rate
conspie=100*(PIEss-1); % quarterly ss inflation rate
consrkn=2.1955;% ss nominal BAA corporate bonds

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
%%%%%%%%%%%%%%%%%%%%%%%%%
%%Single period utility%%
%%%%%%%%%%%%%%%%%%%%%%%%%
LAMBDA=((((C-hab*C(-1)/(1+Stochg))^(1-varrho))*((1-h)^varrho))^(1-sigma_c)-1)/(1-sigma_c);

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
YW=(A*(h)^alp)*(K(-1)/(1+Stochg))^(1-alp)/Delta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Wholesale firms FOC for labour%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PWP*(alp*YW)/h=WP;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Capital law of motion%%%%%
%%Costs of investment case%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=I/I(-1)*(1+Stochg);

K=capqual(+1)*((I*(1-phiX*(X-1-Stochg)^2.0)+(1-delta)*K(-1)/(1+Stochg)));


%%%%%%%%%%%%%%
%%Investment%%
%%%%%%%%%%%%%%
Z1=2.0*phiX*(X-1-Stochg)*X^2*Q*DD(-1);
Q*(1- phiX*(X-1-Stochg)^2.0-2.0*X*phiX*(X-1-Stochg))+Z1(+1)=1;

%%%%%%%%%%%%%%%
%%Fischer Eqn%%
%%%%%%%%%%%%%%%
INVPIE=1/PIE;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Real ex-post interest rate%% Eqn 17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rex=(Rn(-1))*INVPIE;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Balance budget constraint%%% Eqn 18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G=h*WP*tax; 

%%%%%%%%%%%%%%%%%%%%%%
%%Inflation Dynamics%%
%%%%%%%%%%%%%%%%%%%%%%
DD=betta*(LAMBDAC(+1)/LAMBDAC)*(1+Stochg(+1))^((1-varrho)*(1-sigma_c)-1);
DDL=betta*(LAMBDAC/LAMBDAC(-1))*(1+Stochg)^((1-varrho)*(1-sigma_c)-1);
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
%%BGG Banking Sector
%%%%%%%%%%%%%%%%%%
%
%Calibration of banking parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigmaE=STEADY_STATE(sigmaE);
Apsi=STEADY_STATE(Apsi);
mu=STEADY_STATE(mu);
xiE=STEADY_STATE(xiE);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Eqn 31-44
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=((I*(1-phiX*(X-1-Stochg)^2.0)+(1-delta)*capqual*S(-1)/(1+Stochg)));
rho=DGam/((fnGam-mu*fnG)*DGam+(1-fnGam)*(DGam-mu*DG));
p=1/(2*Apsi)*(psi-1+epsilonA+Apsi);
fnGam=1/(4*Apsi)*(psi^2-(1-epsilonA-Apsi)^2)+psi*(1-p);
fnG=1/(4*Apsi)*(psi^2-(1-epsilonA-Apsi)^2);
DGam=1/(2*Apsi)*(1-epsilonA-psi)+1/2;
DG=psi/(2*Apsi);
rhoRex=rho*Rex;
Z=(1-alp)*PWP*YW/(K(-1)/(1+Stochg));
Rk=capqual*((Z+(1-delta)*Q)/Q(-1));
Rk(+1)=rhoRex(+1);
NW=(sigmaE+xiE)*(1-fnGam)*Rk*Q(-1)*K(-1)/(1+Stochg);
phi=Q*K/NW;
phi(-1)*Rk*(fnGam-mu*fnG)=Rex*(phi(-1)-1);
spread=Rk-Rex;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%Resource constraint%% Eqn 45-46
%%%%%%%%%%%%%%%%%%%%%%%
Y=C+CE+G+I+mu*fnG*Rk*Q(-1)*K(-1)*capqual/(1+Stochg);
CE=(1-xiE)*(1-sigmaE)*(1-fnGam)*Rk*Q(-1)*K(-1)*capqual/(1+Stochg);
%%%%%%%%%%%%%%%
%%Taylor rule%%
%%%%%%%%%%%%%%%
%
%Implementable Rule
%
//log((Rn)/(STEADY_STATE(Rn)))=alpha_r*log(Rn(-1)/STEADY_STATE(Rn))
//+alpha_pie*(1-alpha_r)*log(PIE/STEADY_STATE(PIE))
//+alpha_y*(1-alpha_r)*log(Y/STEADY_STATE(Y))+epsMPS;
//
//+0.1*log(phi/STEADY_STATE(phi));
%
%Conventional Taylor Rule
%

log((Rn)/(STEADY_STATE(Rn)))=alpha_r*log((Rn(-1))/(STEADY_STATE(Rn)))
+alpha_pie*(1-alpha_r)*log((PIE)/(STEADY_STATE(PIE)))
+alpha_y*(1-alpha_r)*log((Y/STEADY_STATE(Y))/(YF/STEADY_STATE(YF)))+log(MPS);
//
//+alpha_q*(1-alpha_r)*log((Q/STEADY_STATE(Q))/(QF/STEADY_STATE(QF)))
//+alpha_q*(1-alpha_r)*log(Q/STEADY_STATE(Q))

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
rhorho=rho/STEADY_STATE(rho);
psipsi=psi/STEADY_STATE(psi);
phiphi=phi/STEADY_STATE(phi);
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
dy=log(YY)-log(YY(-1))+trend+epsAtrend;
pinfobs = log(PIEPIE)+conspie; //PIEss-1
robs = log(RnRn)+consr;
rkn_obs = log(RkRk*PIEPIE)+consrkn;
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
var epsAtrend; stderr 0.1;
var epscapqual; stderr 1.0;
end;
/*
steady;

check;

stoch_simul(periods=0, order=1,irf=40)QQ PIEPIE RnRn YY CC  II hh WPWP RR RkRk OUTGAPOUTGAP
spreadspread NWNW rhorho p psipsi phiphi;

*/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%BAYESIAN ESTIMATION STARTS HERE%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estimated_params;
// PARAM NAME, INITVAL, LB, UB, PRIOR_SHAPE, PRIOR_P1, PRIOR_P2, PRIOR_P3, PRIOR_P4, JSCALE
// PRIOR_SHAPE: BETA_PDF, GAMMA_PDF, NORMAL_PDF, INV_GAMMA_PDF
// Priors used as in SW07
//
stderr epsA,  INV_GAMMA_PDF,0.1,2;      //technology
stderr epsG,  INV_GAMMA_PDF,0.5,2;      //government spending - rescaled stderr/gy
stderr epsMPS, INV_GAMMA_PDF,0.1,2;      //interest rate rule
stderr epsMS, INV_GAMMA_PDF,0.1,2;      //mark-up
stderr epsAtrend,  INV_GAMMA_PDF,0.1,2;      //trend shock
stderr epscapqual,INV_GAMMA_PDF,0.1,2;      // capital quality shock




rhoA,            BETA_PDF,     0.5,0.20;   //AR1 technology
rhoG,          BETA_PDF,     0.5,0.20;   //AR1 government spending
rhoMS,          BETA_PDF,     0.5,0.20;   //AR1 mark-up
rhocapqual,     BETA_PDF,     0.5,0.20;   //AR1 capital quality

phiX,           NORMAL_PDF,   2.00, 1.50; //investment adj cost
sigma_c,        NORMAL_PDF,   1.50, 0.375;//consumption utility
hab,           BETA_PDF,     0.70, 0.10; //habit
xi,             BETA_PDF,     0.50, 0.10; //calvo prices
gamp,          BETA_PDF,     0.50, 0.15; //indexation prices


alpha_pie,      NORMAL_PDF,   2.00, 0.25; //feedback inflation
alpha_r,        BETA_PDF,     0.75, 0.10; //lagged interest rate
alpha_y,         NORMAL_PDF,   0.125,0.05; //feedback output gap 


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
//estimation 2nd stage - MCMC
estimation(datafile=us_data,mode_compute=0,first_obs=5,presample=4,mode_file=BGG_RES_Course_est_PI_mode,
              prefilter=0,mh_replic=1000,mh_nblocks=2,mh_jscale=0.30,mh_drop=0.2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%POST-ESTIMATION SIMULATIONS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
// generate IRfunctions and Moments of obs based on the post. distribution
// simulate the model by sampling shocks from their distribution
// (with the given variabilities) to see how the system reacts for
// every period and generates the simulated endogenous var, and runs an IRF,
// similarly by sampling shocks from the distribution (with 1 s.d.),
// to see how the system reacts for the given 20 periods.
stoch_simul(irf=20,ar=10,conditional_variance_decomposition=[1 4 10 100]) dy pinfobs robs rkn_obs;

