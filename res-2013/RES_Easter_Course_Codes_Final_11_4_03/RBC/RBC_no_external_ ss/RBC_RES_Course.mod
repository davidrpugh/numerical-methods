//Dynare set-up:  One-Sector RBC Model without costs of investment and with non-distortionary taxes
//The basic RBC model with all parameters exogenous
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
var   LAMBDA LAMBDAC R C  WP h Y YW PWP K  I  tax   A  G
RR   YY CC hh WPWP II KK KYW IY CY;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF EXOGENOUS VARIABLES%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varexo epsA epsG ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF PARAMETERS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters gy varrho alp c zzeta betta delta sigma_c rhoA rhoG Ass;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%CALIBRATION OF PARAMETERS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gy=0.2;
alp=0.70;
zzeta=7.0;
c=1/zzeta;
betta=0.99;
delta=0.0250;
sigma_c=2.0;
varrho=0.684;
%Choice of Units
Ass=1;
%shock persistence
rhoA=0.75;
rhoG=0.75;





% ----------------------------
% *** DSGE-Model-equations ***
% ----------------------------
model;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Variables used in steady state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KYW=K/YW;
IY=I/Y;
CY=C/Y;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%INITIAL GUESSES FOR STEADY-STATE COMPUTATION%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initval;
A=Ass;
R=1.0/betta;
PWP=1-1/zzeta;
KYW=PWP*(1-alp)/(R-1+delta);
IY=(1-alp)*delta/(R-1+delta);
CY=1-IY-gy;
h=(1-varrho)*alp/(varrho*CY+(1-varrho)*alp);
YW=A*h*KYW^((1-alp)/alp);
Y=(1-c)*YW;
K=KYW*YW;
I=IY*Y;
G=gy*Y;
C=Y-I-G;
WP=alp*PWP*YW/h;
tax=G/(WP*h);
LAMBDA=1/(1-sigma_c)*(C^((1-varrho)*(1-sigma_c))*(1-h)^(varrho*(1-sigma_c))-1);
LAMBDAC=(1-varrho)*C^((1-varrho)*(1-sigma_c)-1)*(1-h)^(varrho*(1-sigma_c));

YY=1;
CC=1;
hh=1;
WPWP=1; 
II=1; 
KK=1;
RR=1;
end;

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
