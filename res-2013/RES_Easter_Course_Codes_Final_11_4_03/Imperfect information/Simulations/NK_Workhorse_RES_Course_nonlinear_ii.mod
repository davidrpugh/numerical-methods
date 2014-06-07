//------------------------------------------------
//filename; lppy2011_optpol.mod
//Levine at el (2011) model, CDMA WP1002
//perfect-imperfect info. 
//optimal policy excercies
//outputs saved in the 'filename' folder
//only applicable to Dynare 4.1.2 for now
//this version - 23 Aug 2011
//------------------------------------------------

var A Z MUC PI R MC N Y C MUN H J D MStot MSper MStem G ZR YL PIL RL;
var YY CC RR PIPI NN;
varexo eps_a eps_msper eps_mstem eps_g; 

parameters beta xi hc sigma rho_a rho_g rho_ms rho_r thetap eta zeta varrho gy Ass; 
beta  = 0.99;
//estimated parameters------
xi    = 0.6319;
hc    = 0.8817;
sigma = 2.0066; 
rho_a = 0.9407;
rho_g = 0.9587;
rho_ms= 0.4997;
rho_r = 0.7316;
thetap= 0.6800;
//--------------------------
eta   = 6.00;
zeta  = 7.67;
gy    = 0.30;
Ass   = 1.00;
//------------------------
//calibration of varrho - with wd (=0.4) and cy
//varrho =(1-wd)/(1+wd*(cy*(1-hc)/((1-1/zeta)*(1-1/eta))-1))
//------------------------
varrho =(1-0.4)/(1+0.4*((1-gy)*(1-hc)/((1-1/zeta)*(1-1/eta))-1));


model; 
//--------------------------------------------------------------
//define auxiliary variables for instrument and lagged variables
//--------------------------------------------------------------
ZR=R;
RL=ZR(-1);
YL=Y(-1);
PIL=PI(-1);
Z=C(-1);

//-----------------
//structural model
//-----------------
1 = beta*(R)*(MUC(+1)/(MUC*PI(+1)));
MUC = (1-varrho)*(C-hc*Z)^((1-varrho)*(1-sigma)-1)*(1-N)^(varrho*(1-sigma));
MUN = -(C-hc*Z)^((1-varrho)*(1-sigma))*varrho*(1-N)^(varrho*(1-sigma)-1);
MC  = -MUN/MUC/A;
N   = Y*D/A;
PI^(1-zeta)*H - xi*beta*H(+1)= Y*MUC;
PI^(-zeta)*J - xi*beta*J(+1)= (zeta/(zeta-1))*MC*MStot*Y*MUC*(eta/(eta-1));
1 = xi*PI^(zeta-1)+(1-xi)*(J/H/PI)^(1-zeta);
D=xi*(PI^zeta)*D(-1)+(1-xi)*((J/H/PI)^(-zeta));
Y = C+G;
MStot=MSper*MStem;
log(A) = rho_a*log(A(-1))+eps_a;
log(G/STEADY_STATE(G))  = rho_g*log(G(-1)/STEADY_STATE(G))+eps_g;
log(MSper) = rho_ms*log(MSper(-1))+eps_msper;
log(MStem) = eps_mstem;

//-----------------------------------------------
//policy rule
//add a tag to identify the rule in the Jacobian
//feedback coeff. must not be expressions (func.)
//-----------------------------------------------
//[tag1 = 'aceslq_sim_rule' , tag2 = 'inflation']
log((R)/(STEADY_STATE(R)))=rho_r*log(RL/STEADY_STATE(RL))+thetap*log(PI/STEADY_STATE(PI));

YY=Y/STEADY_STATE(Y);
CC=C/STEADY_STATE(C);
RR=(R)/(STEADY_STATE(R));
PIPI=PI / STEADY_STATE(PI);
NN=(N)/(STEADY_STATE(N));
end;

//-----------------------------------------------------------------------------------
//initial values for the Ramsey steady state
//only need to specify a value for seedvar if exists external steady program, e.g. PI
//-----------------------------------------------------------------------------------
initval;
A     		 =1;
Z     		 =0.28;
MUC   		 =4.39305;
PI    		 =1;
R     		 =1.0101;
MC    		 =0.724685;
N     		 =0.4;
Y     		 =0.4;
C     		 =0.28;
MUN   		 =-3.18358;
H     		 =4.6932;
J     		 =4.6932;
D     		 =1;
MStot 		 =1;
MSper 		 =1;
MStem 		 =1;
G     		 =0.12;
ZR    		 =1.0101;
YL    		 =0.4;
PIL   		 =1;
RL    		 =1.0101;
YY           =1;
CC           =1;
RR           =1;
PIPI         =1;
NN           =1;
end;

shocks;
//var eps_a; stderr 0.7535;
//var eps_g; stderr 2.3105;
//var eps_msper; stderr 0.05/((1-beta*xi)*(1-xi)/xi);
//var eps_mstem; stderr 0.05/((1-beta*xi)*(1-xi)/xi);
var eps_a; stderr 1;
var eps_g; stderr 1;
var eps_msper; stderr 1;
var eps_mstem; stderr 1;
end;

steady;
check;
varobs YY MSper PIPI CC;
stoch_simul(partial_information,irf=30) PIPI NN CC YY RR; 
