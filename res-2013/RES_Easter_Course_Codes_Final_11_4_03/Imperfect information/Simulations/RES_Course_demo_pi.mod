//Simple demo model (RES School 16/4/2013)
//In the notes: Imperfect Information for Simulation and Estimation in Dynare
//Standard Dynare assumption: agants have perfect information
//agents observe [pi x w]'


var pi x;
varexo eps w; 

parameters beta rho; 
beta = 0.99;
rho  = 0.90;

model(linear);
pi(+1) = 1/beta*pi+x+w;
x      = rho*x(-1)+eps;
end;

shocks;
var eps; stderr 1;
var w;   stderr 1; 
end;


stoch_simul(irf=10) pi;

//***********************************************
// lines 19 - 25 equiv to the following:
//shocks;
//var eps; stderr 1;
//var w;   stderr 1;
//end;
//varobs pi x;
//stoch_simul(partial_information, irf=10) pi;
//**********************************************
//
//OR -- the reduced-form solution of this particular example
//shows that perfect info.(AI) is equiv to setting var(w)=0 under II
//
//**********************************************
// lines 19 - 25 equiv to the following:
//shocks;
//var eps; stderr 1;
//var w;   stderr 0;
//end;
//varobs pi;
//stoch_simul(partial_information, irf=10) pi;
//*********************************************