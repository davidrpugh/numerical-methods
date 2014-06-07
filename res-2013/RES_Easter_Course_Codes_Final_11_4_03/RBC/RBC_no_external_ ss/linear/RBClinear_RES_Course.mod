//Dynare set-up:  One-Sector RBC Model without costs of investment and with non-distortionary taxes
//Linear version
//(c) CIMS Univeristy of Surrey
//The Science and  Art of DSGE Modelling: Construction, Calibration, Estimation and Policy
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF ENDOGENOUS VARIABLES%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var LAMBDAC LAMBDAL R  h WP YW  
Y K  I  tax  C  A  G ; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF EXOGENOUS VARIABLES%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varexo epsA epsG ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF PARAMETERS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters gy varrho alp c zzeta betta delta sigma_c rhoA rhoG cy iy hss Rss Ass;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%CALIBRATION OF PARAMETERS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gy=0.2;
varrho=0.684;
alp=0.70;
zzeta=7.0;
c=1/zzeta;
betta=0.99;
delta=0.0250;
sigma_c=2.0;
%Choice of Units
Ass=1;
%shock persistence
rhoA=0.75;
rhoG=0.75;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%STEADY STATE RELATIONSHIPS FROM NON-LINEAR SET_UP%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rss=1.0/betta;
iy=((1-alp)*delta)/((Rss-1+delta));
cy=1-iy-gy;
hss=alp*(1-varrho)/(cy*varrho+alp*(1-varrho));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// ----------------------------
// *** DSGE-Model-equations ***
// ----------------------------
model(linear);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Marginal utility of consumption%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LAMBDAC=((1-varrho)*(1-sigma_c)-1)*C+ varrho*(sigma_c-1)*(hss/(1-hss))*h;
%
%deliberate mistake in linearization
%
%LAMBDAC=((1-varrho)*(1-sigma_c)-1)*C+ (sigma_c-1)*(hss/(1-hss))*h; 
%

%%%%%%%%%%%%%%%%%%
%%Euler equation%%
%%%%%%%%%%%%%%%%%%
LAMBDAC(+1)=LAMBDAC-R; 

%%%%%%%%%%%%%%%%%%%%%
%%Labour supply foc%%
%%%%%%%%%%%%%%%%%%%%%
LAMBDAL=LAMBDAC+C+(hss/(1-hss))*h;
WP=LAMBDAL-LAMBDAC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Wholesale and retail sector relation%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y=YW;

%%%%%%%%%%%%%%%%%%%%%%%
%%Production Function%%
%%%%%%%%%%%%%%%%%%%%%%%
YW=alp*A+alp*h+(1-alp)*K(-1);

%%%%%%%%%%%%%%%%%%%%%%%%
%%Wholesale firms FOCs%%
%%%%%%%%%%%%%%%%%%%%%%%%
Y(+1)-K=(Rss)/(Rss-1+delta)*R;
WP=YW-h;

%%%%%%%%%%%%%%%%%%%%%%%
%%Resource constraint%%
%%%%%%%%%%%%%%%%%%%%%%%
Y=cy*C+iy*I+gy*G;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%Capital law of motion%% 
%%%%%%%%%%%%%%%%%%%%%%%%%
K=delta*I +(1-delta)*K(-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Non-distortionart taxes%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tax=G;


%%%%%%%%%%%%%%%%%%%
%%Shock processes%%
%%%%%%%%%%%%%%%%%%%
A=rhoA*A(-1)+epsA;
G=rhoG*G(-1)+epsG;
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


stoch_simul(order=1,irf=40) Y C I h WP R;
