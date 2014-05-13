// following are the 'endogenous' variables:
var lc, lk, ly, a;
// following are the shocks
varexo e;

parameters beta, alpha, gamma, rho, delta;

alpha = 0.36;
gamma = 2;
rho   = 0.95;
beta  = 0.99;
delta = .02;

[css,kss,yss,zss] = ssrbc(alpha, beta, delta);

model;
    1/exp(gamma*(lc)) = beta/exp(gamma*(lc(+1)))*(alpha*exp(ly(+1)-lk) + 1 - delta);
    exp(ly) = exp(alpha*lk(-1))*exp(a); 
    a = rho*a(-1) + e;
    exp(ly) = exp(lc) + exp(lk)-(1-delta)*exp(lk(-1));
end;

initval;
lc=log(css);
lk=log(kss);
ly=log(css+delta*kss);
a=0;
end;

%verifies that the steady state computed by rbc_steadystate.m is correct
steady;

shocks;
var e; 
stderr 0.01;
end;

stoch_simul(order=2,irf=20,nograph) lc lk ly a;