//------------------------------------------------
//filename: lppy2011_optpol.mod - linear
//Levine at el (2011) model, CDMA WP1002
//perfect-imperfect info. 
//optimal policy excercies
//outputs saved in the 'filename' folder
//only applicable to Dynare 4.1.2 for now
//this version - 26 July 2012
//------------------------------------------------

options_.usePartInfo=1;

var pi mc mun muc c y n r g a ms mstot;
varexo eps_a eps_m eps_ms eps_g eps_e; 

parameters beta xi hc wd sigma gamma rho_g rho_a rho_r rho_ms thetap cy alpha; 
parameters eta zeta; // these are the ones in the nonlinear model
// calibrated parameters
beta = 0.99;
wd = 0.40;
cy = 0.70; // cy = Cbar/Ybar
//varrho =(1-wd)/(1+wd*(cy*(1-hc)/((1-1/zeta)*(1-1/eta))-1));
////////////
xi    = 0.76;
hc    = 0.0;
sigma = 2.30;
rho_a = 0.95;
rho_g = 0.92;
rho_ms= 0.42;
rho_r = 0.63;
thetap= 1.90*(1-rho_r);
eta   = 6.00;
zeta  = 7.67;
gamma = 0.00;
alpha = 1.00;

model(linear);
//--------------------------------------
//calibration of varrho - with wd and cy
//--------------------------------------
#varrho =(1-wd)/(1+wd*(cy*(1-hc)/((1-1/zeta)*(1-1/eta))-1));

//-----------------
//structual model
//-----------------
pi = (beta/(1+beta*gamma))*pi(+1)+(gamma/(1+beta*gamma))*pi(-1)+(((1-beta*xi)*(1-xi))/((1+beta*gamma)*xi))*(mc+mstot);
mc = mun-muc-a+(1-alpha)*n;
mun = (c-hc*c(-1))/(1-hc)+wd*n/(1-wd)+muc;
muc = ((1-varrho)*(1-sigma)-1)*(c-hc*c(-1))/(1-hc)-wd*varrho*(1-sigma)*n/(1-wd);
muc(+1) = muc-(r-pi(+1));
y = cy*c+(1-cy)*g;
y=alpha*n+a;

//exogenous processes
//mark-up shock process specified as in SW 2007:
r = rho_r*r(-1)+thetap*pi+eps_e;
g = rho_g*g(-1)+eps_g;
a = rho_a*a(-1)+eps_a;
ms = rho_ms*ms(-1)+eps_ms; // persistent mark-up shock  //SW07-choice
mstot = ms+eps_m; // persistent+transient components where eps_m iid-Normal price mark-up shock
end;

shocks;
var eps_g;  stderr 1;
var eps_a;  stderr 1;
var eps_ms; stderr 1;
var eps_m;  stderr 1;
var eps_e;  stderr 1;
end;

check;
//stoch_simul(linear,irf=30)pi n c y r;

estimated_params;
// priors obtained from the happiness paper (i.e. Levin et al. (2005)) and SW 2007
sigma,  normal_pdf, 1.50,0.375;   // SW07-chioce
xi,     beta_pdf,   0.50, 0.10;   // SW07-chioce
//hc,     beta_pdf,   0.50, 0.10;
thetap, normal_pdf, 1.50, 0.25;   // SW07-chioce 

// AR1 parameters
rho_r,   beta_pdf, 0.80, 0.10;
rho_a,   beta_pdf, 0.85, 0.10;
rho_g,   beta_pdf, 0.85, 0.10;
rho_ms,  beta_pdf, 0.50 ,0.20; // persistent price-markup

// stds of AR(1) innovations
stderr eps_a,    inv_gamma_pdf, 0.125, 2.00;
stderr eps_g,    inv_gamma_pdf, 0.50, 2.00;
stderr eps_ms,   inv_gamma_pdf, 0.10, 2.00; // shock price markup AR1 - SW07

// stds of iid shocks:
stderr eps_m,    inv_gamma_pdf, 0.10, 2.00; // price markup iid-normal - SW07
stderr eps_e,    inv_gamma_pdf, 0.10, 2.00;  
end;

// Data 1966Q1 - 2006Q4
// Data 1980Q1 - 2006Q4
varobs pi y r;
options_.plot_priors=0;
estimation(datafile=data_lin_81Q1_06Q4, first_obs=1, mode_compute=4,
           mh_replic=0, mh_drop=0.3, mh_jscale=0.3);
//estimation(datafile=data_lin_81Q1_06Q4, first_obs=1, mode_compute=0, mode_file=mod_z_ai_linear_subdat2_est_mode,
//           mh_replic=100000, mh_drop=0.3, mh_jscale=0.25);
//estimation(datafile=data_lin_81Q1_06Q4, first_obs=1, mode_compute=4,
//           mh_replic=0, mh_drop=0.3, mh_jscale=0.3);

//stoch_simul(irf=60) pi y r;
