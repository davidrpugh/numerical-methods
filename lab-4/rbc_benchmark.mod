////////// Declare variables //////////

/* 

///// Endogenous variables /////
   
   K: capital
   C: consumption
   Z: productivity
   Y: output
   I: investment
   W: real wage
   r: net interest rate
   L: labor supply
   zero_profit: zero profit condition

*/
var Y, K, L, C, I, r, W, Z, zero_profit;
///// Exogenous variables /////

// eps_z: productivity shock 
varexo eps_z;
////////// Declare parameters //////////
parameters b, alpha, beta, delta, sigma, theta, omega, rho_z, KSS, LSS, YSS;
// discount factor
# insert your code here!

// coefficient of relative risk aversion
theta = 2.5;

// inverse elasticity of substitution for leisure
# insert your code here!

// weight for leisure in utility
b = 2.0;

// elasticity of substitution between capital and labor
# insert your code here!

// relative weight of capital in production
alpha = 0.40;

// depreciation rate of capital
# insert your code here!

// persistence of productivity process
rho_z = 0.95;

// steady state labor supply (computed for theta=omega=delta=sigma=1)
LSS = (1 - alpha) / (1 - alpha + b * (1 - alpha * beta));

// steady state capital (computed for theta=omega=delta=sigma=1)
KSS =  (alpha * beta)^(1 / (1 - alpha)) * LSS;

// steady state output (computed for theta=omega=delta=sigma=1)
YSS = KSS^alpha * LSS^(1 - alpha);

////////// Model equations //////////

model;

// production (i.e., equation 3.1)
# insert your code here!

// productivity process (i.e., equation 3.2)
Z = rho_z * Z(-1) + eps_z;

// real wage (i.e., equation 3.3)
# insert your code here!

// rental rate (i.e., equation 3.4)
r = ((alpha *K(-1)^(-1 / sigma)) / (alpha * K(-1)^((sigma - 1) / sigma) + (1 - alpha) * L^((sigma - 1) / sigma))) * Y;

// check that zero profit condition holds (i.e., equation 3.5)
zero_profit = Y - W * L - r * K(-1);

// household budget constraint (i.e., equation 3.6a)
# insert your code here!

// household evolution of capital (i.e., equation 3.6b)
K = (1 - delta) * K(-1) + I;

// consumption Euler equation (i.e., equation 3.7)
# insert your code here!

// intra-temporal consumption/labor trade-off (i.e., equation 3.8)
b * (1 - L)^(-omega) = W * C^(-theta);

end;
////////// Initial values for computing steady state //////////

initval;
K = KSS;
L = LSS;
Y = YSS;
C = (1 - alpha * beta) * YSS;
Z = 0.0;
I = alpha * beta * YSS;
W = (1 - alpha) * (YSS / LSS);
r = (1 / beta);
zero_profit = 0.0;
end;
////////// Variance covariance matrix //////////
vcov = [0.007];