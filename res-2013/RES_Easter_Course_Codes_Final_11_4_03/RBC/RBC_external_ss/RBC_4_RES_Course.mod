//Dynare set-up:  One-Sector RBC Model without costs of investment and with non-distortionary taxes
//RBC model with varrho calibrated to hit hours, delta to hit investment share, beta to hit real interest rate   
//(c) CIMS Univeristy of Surrey
//The Science and  Art of DSGE Modelling: Construction, Calibration, Estimation and Policy
//////////////////////////////////////////////////////////////////////////////////////////////////
//This is a free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.  See <http://www.gnu.org/licenses/> for more information.
//////////////////////////////////////////////////////////////////////////////////////////////////
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF ENDOGENOUS VARIABLES%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var   LAMBDA LAMBDAC R C  WP h Y
YW PWP K  I  tax   A  G
RR   YY CC hh WPWP II KK
varrho delta betta; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF EXOGENOUS VARIABLES%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varexo epsA epsG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF PARAMETERS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters gy alp c zzeta  sigma_c rhoA rhoG Ass hobs iobs Robs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PARAMETER SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gy=0.2;
alp=0.70;
zzeta=7.0;
c=1/zzeta;
sigma_c=2.0;
hobs=0.35;
iobs=0.2;
Robs=1.01;
%
%Choice of Units
Ass=1;
%shock persistence
rhoA=0.75;
rhoG=0.75;



% ----------------------------
% *** DSGE-Model-equations ***
% ----------------------------
model;
%------------------------------
%Calibration of varrho, delta and betta
%------------------------------
%

varrho=STEADY_STATE(varrho);
delta=STEADY_STATE(delta);
betta=STEADY_STATE(betta);

%%%%%%%%%%%%%%%%%%%%%%%%%
%%Single period utility%%
%%%%%%%%%%%%%%%%%%%%%%%%%
LAMBDA=(((C^(1-varrho))*((1-h)^varrho))^(1-sigma_c)-1)/(1-sigma_c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Marginal utility of consumption%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LAMBDAC=(1-varrho)*(C^((1-varrho)*(1-sigma_c)-1))*((1-h)^(varrho*(1-sigma_c)));

%%%%%%%%%%%%%%%%%%
%%Euler equation%%
%%%%%%%%%%%%%%%%%%
%
%R set in period t to yield interest over period t+1
%
LAMBDAC=betta*(R)*LAMBDAC(+1);

%%%%%%%%%%%%%%%%%%%%%
%%Labour supply foc%%
%%%%%%%%%%%%%%%%%%%%%
varrho*(C^((1-varrho)*(1-sigma_c)))*((1-h)^(varrho*(1-sigma_c)-1))/LAMBDAC=WP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Wholesale and retail sector relation%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y=(1-c)*YW;

%%%%%%%%%%%%%%%%%%%%%%%
%%Production Function%%
%%%%%%%%%%%%%%%%%%%%%%%
%
%Note: Capital stock is end-of-period
%
YW=((A*h)^alp)*(K(-1))^(1-alp);

%%%%%%%%%%%%%%%%%%%%%%%%
%%Wholesale firms FOCs%%
%%%%%%%%%%%%%%%%%%%%%%%%
(1-alp)*PWP(+1)*YW(+1)/K=R-1+delta;
PWP*(alp*YW)/h=WP;

%%%%%%%%%%%%%%%%%%%%%%%
%%Resource constraint%%
%%%%%%%%%%%%%%%%%%%%%%%
Y=C+G+I;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%Capital law of motion%% 
%%%%%%%%%%%%%%%%%%%%%%%%%
K=I+(1-delta)*K(-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Balance budget constraint%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G=h*WP*tax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Mark-up Monopolistic pricing%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
PWP=(1-1/zzeta);

%%%%%%%%%%%%%%%%%%%
%%Shock processes%%
%%%%%%%%%%%%%%%%%%%
log(A)-log(STEADY_STATE(A))=rhoA*(log(A(-1))-log(STEADY_STATE(A)))+epsA;
log(G)-log(STEADY_STATE(G))=rhoG*(log(G(-1))-log(STEADY_STATE(G)))+epsG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Variables in deviation form for IRFs%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
YY=Y/STEADY_STATE(Y);
KK=K/STEADY_STATE(K);
II=I/STEADY_STATE(I);
CC=C/STEADY_STATE(C);
WPWP=WP/STEADY_STATE(WP);
hh=h/STEADY_STATE(h);
RR=(R)/(STEADY_STATE(R));
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
var epsA; stderr 1;
var epsG; stderr 1;
end;

steady;

check;


stoch_simul(order=1,irf=40) YY CC II hh WPWP RR;
