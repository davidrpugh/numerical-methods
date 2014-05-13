function [ggk, gga, gge, ggsig, ggkk, ggaa, ggee, ggsigsig, ggka, ggke, ggksig, ggae, ggasig, ggesig, kss] = ...
    second(alpha,beta,delta,gamma,rho,Ve,lg)

%This finds the second order expansion of the policy rule which solves the 
%model, max(sum, t=0,...,infinity) beta^t * u(c), subject to
%c+kp-(1-delta)*k=k^alpha*exp((1-alpha)*a(t)), where a(t) = rho*a(t-1) +
%epsilon(t), where Eepsilon(t)^2 = Ve. Also, u(c)=c^(1-gamma)/(1-gamma).
%
%The policy rule is:
%
% kp = g(k,a,epsilon;sig)
%
% where
%
% g = kss + gk*(k - kss) + ga*a + ge*epsilon + gsig*sig ...
%     + (1/2)*(gkk*(k - kss)^2 + gaa*a^2 + gee*epsilon^2 + gsigsig*sig^2) ...
%     + gka*(k - kss)*a + gke*(k - kss)*epsilon + gksig*(k - kss)*sig ...
%     + gae*a*epsilon + gasig*a*sig + gesig*epsilon*sig ;
%
%In the above expression, 
%   kss is steady state capital
%   k is beginning of time t capital
%   kp is the beginning of time t+1 capital, chosen in time t
%   a is a(t-1)
%   epsilon is epsilon(t)
%   sig is a constant, which should be set to unity
%
%If lg = 1, then kss, k and kp are the log of the capital stock. Otherwise,
%they are the level.
%
%The strategy used here follows the one described in Swanson, Anderson, Levin,
%'Higher-Order Perturbation Solutions to Dynamic, Discrete-Time Rational
%Expectations Models', manuscript, Federal Reserve Board, July 8, 2005.
%This code is meant for class purposes only, and is very very slow (the code below
%would be substantially faster if subs were replaced by eval). For superior code, 
%which is much, much faster and much, much
%more general, see Swanson-Anderson-Levin's PerturbationAIM software.
    
%steady state:
kss=(alpha*beta/(1-(1-delta)*beta))^(1/(1-alpha));
if lg == 1;kss=log(kss);end

ge=4;%this is undone by the subsequent line, and is present to deal with an incomprehensible quirk in Matlab

syms c k a epsilon x n kp g gk ga ge gsig
syms gkk gaa gee gsigsig gka gke gksig gae gasig gesig sig epsilonp

%this is the Taylor series expansion of the policy rule about nonstochastic steady state, whose coefficients
%we will compute using the implicit function theorem
g = kss + gk*(k - kss) + ga*a + ge*epsilon + gsig*sig ...
    + (1/2)*(gkk*(k - kss)^2 + gaa*a^2 + gee*epsilon^2 + gsigsig*sig^2) ...
    + gka*(k - kss)*a + gke*(k - kss)*epsilon + gksig*(k - kss)*sig ...
    + gae*a*epsilon + gasig*a*sig + gesig*epsilon*sig ;

%g=g(k,a,epsilon,sig), k=k at time t; a = a at time t-1; epsilon = eps at time t, sig = constant

%here, f is the time t production function and fk is the time t marginal product of capital
if lg == 1
    f= exp(alpha*k) * exp(rho*a + epsilon) + (1-delta) * exp(k) ;
    fk = alpha * exp(k*(alpha-1)) * exp(rho*a + epsilon) + (1-delta);
else
    f= k^alpha * exp(rho*a + epsilon) + (1-delta) * k ;
    fk = alpha * k^(alpha-1) * exp(rho*a + epsilon) + (1-delta);
end

%c is time t consumption, obtained from the resource constraint, with t+1 capital left as a variable, kp
if lg == 1
    c = f - exp(kp);
else
    c = f - kp;
end

%u is the time t level of the utility function, and uc is the time t marginal utility
u = (c^(1-gamma))/(1-gamma);
up= c^(-gamma);%this is a function of a, epsilon, k, kp               
uc = subs(up,kp,g);%replace k at time t+1 (i.e., kp) with the policy rule, g

%
%now go for the time t+1 marginal utility of consumption, ucp
ucp = subs(up,[epsilon],[sig*epsilonp]);%replace epsilon in up by sig*epsilonp (epsilonp is t+1 epsilon)
ucp = subs(ucp,[a],[rho*a+epsilon]);%replace a by rho*a+epsilon
ucp = subs(ucp,k,g);%replace period t capital by t+1 capital, g

%
%go for the period t+1 capital decision
gp = subs(g,[epsilon],[sig*epsilonp]);
gp = subs(gp,[a],[rho*a+epsilon]);
gp = subs(gp,k,g);

%
%the final step in obtaining t+1 marginal utility of consumption: replace period t+1 
%capital with period t+1 capital decision
ucp = subs(ucp,kp,gp);

%
%go after the period t+1 marginal product of capital
fkp = subs(fk,[epsilon],[sig*epsilonp]);
fkp = subs(fkp,[a],[rho*a+epsilon]);
fkp = subs(fkp,k,g);

%
%construct the period t intertemporal error function. Applying the implicit
%function theorem to this will allow us to solve for the coefficients in g
R = uc - beta*ucp*fkp;

Rk=diff(R,k);%differentiate R with respect to k
RRk=subs(Rk,[k a epsilon sig epsilonp ],[kss 0 0 0 0]);%evaluate derivative in steady state and replace 
                                                       %epsilonp by 0. This is because the only way 
                                                       %epsilonp enters in Rk is linearly, and so this 
                                                       %goes to zero upon application of the expectation operator.
gg=roots(sym2poly(RRk));%RRk is now a polynomial in gk, and sym2poly recovers the coefficients of this polynomial
                        %roots then finds the roots of this polynomial.
                        %Save the root that is less than unity in absolute
                        %value
ggk=gg(1);%the notation used here is that gij denotes a variable name and ggij is the corresponding value
if abs(ggk) > 1
    ggk=gg(2);
end
if abs(ggk) > 1
    error('fatal(second) too many explosive eigenvalues')
end

%the following operations are the analog of the ones used to obtain ggk
Ra = diff(R,a);
RRa = subs(Ra,[k a epsilon sig epsilonp gk],[kss 0 0 0 0 ggk]);
gga = roots(sym2poly(RRa));

Re = diff(R,epsilon);
RRe = subs(Re,[k a epsilon sig epsilonp gk ga],[kss 0 0 0 0 ggk gga]);
gge = roots(sym2poly(RRe));

Rsig=diff(R,sig);
RRsig = subs(Rsig,[k a epsilon sig epsilonp gk],[kss 0 0 0 0 ggk]);
ggsig = roots(sym2poly(RRsig));

sprintf('Following are the terms in the first-order expansion:\n')
sprintf('gk = %8.7f, ga = %8.7f, ge = %8.7f, gsig = %8.7f \n',ggk,gga,gge,ggsig)

Rkk=diff(Rk,k);
Rkk = subs(Rkk,[k a epsilon sig epsilonp gk ga ge gsig],[kss 0 0 0 0 ggk gga gge ggsig]);
ggkk=roots(sym2poly(Rkk));

Rka=diff(Rk,a);
Rka = subs(Rka,[k a epsilon sig epsilonp gk ga ge gsig gkk],[kss 0 0 0 0 ggk gga gge ggsig ggkk]);
ggka=roots(sym2poly(Rka));

Rke=diff(Rk,epsilon);
Rke = subs(Rke,[k a epsilon sig epsilonp gk ga ge gsig gkk gka],[kss 0 0 0 0 ggk gga gge ggsig ggkk ggka]);
ggke = roots(sym2poly(Rke));

Rksig=diff(Rk,sig);
Rksig = subs(Rksig,[k a epsilon sig epsilonp gk ga ge gsig],[kss 0 0 0 0 ggk gga gge ggsig]);
ggksig = roots(sym2poly(Rksig));

Raa=diff(Ra,a);
Raa = subs(Raa,[k a epsilon sig epsilonp gk ga ge gsig gka gkk],[kss 0 0 0 0 ggk gga gge ggsig ggka ggkk]);
ggaa = roots(sym2poly(Raa));

Rae=diff(Ra,epsilon);
Rae = subs(Rae,[k a epsilon sig epsilonp gk ga ge gsig gaa gka gkk],[kss 0 0 0 0 ggk gga gge ggsig ggaa ggka ggkk]);
ggae = roots(sym2poly(Rae));

Rasig=diff(Ra,sig);
Rasig = subs(Rasig,[k a epsilon sig epsilonp gk ga ge gsig gksig],[kss 0 0 0 0 ggk gga gge ggsig ggksig]);
ggasig = roots(sym2poly(Rasig));


Ree=diff(Re,epsilon);
Ree = subs(Ree,[k a epsilon sig epsilonp gk ga ge gsig gaa gka gkk],[kss 0 0 0 0 ggk gga gge ggsig ggaa ggka ggkk]);
ggee = roots(sym2poly(Ree));

Resig=diff(Re,sig);
Resig = subs(Resig,[k a epsilon sig epsilonp gk ga ge gsig gasig gksig],[kss 0 0 0 0 ggk gga gge ggsig ggasig ggksig]);
ggesig = roots(sym2poly(Resig));

Rsigsig=diff(Rsig,sig);
Rsigsig2=diff(Rsigsig,epsilonp,2);%take the second derivative of Rsigsig w.r.t epsilonp, 
                                  %to isolate coefficient on epsilonp^2
Rsigsig2= subs(Rsigsig2,[k a epsilon sig epsilonp gk ga ge gsig gee],[kss 0 0 0 0 ggk gga gge ggsig ggee]);%evaluate the second derivative in ss
Rsigsig = subs(Rsigsig,[k a epsilon sig epsilonp gk ga ge gsig gee],[kss 0 0 0 0 ggk gga gge ggsig ggee]);
Rsigsignew = Rsigsig +(1/2)*Rsigsig2*Ve;%make sure that E[epsilonp^2] is replaced by Ve 
ggsigsig = roots(sym2poly(Rsigsignew));

%
sprintf('Following are the second-order terms in the expansion \n')
sprintf('gkk = %8.7f, gka = %8.7f, gke = %8.7f, gksig = %8.7f, gaa = %8.7f \n',ggkk,ggka,ggke,ggksig,ggaa)
sprintf('gae = %8.7f, gasig = %8.7f, gee = %8.7f, gesig = %8.7f, gsigsig = %8.7f \n',ggae, ggasig, ggee, ggesig,ggsigsig)
