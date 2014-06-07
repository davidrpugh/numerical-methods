//Simple demo model (RES School 16/4/2013)
//In the notes: Imperfect Information for Simulation and Estimation in Dynare
//Relaxed assumption: agants have imperfect information (II) - possible symmetric info.
//agents observe [pi]' so have one observable and two shocks (eps and w)
//rational learning about the unobservable shock via Kalman updating

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
var w;   stderr 1.41;
end;

varobs pi;
stoch_simul(partial_information, irf=10) pi;