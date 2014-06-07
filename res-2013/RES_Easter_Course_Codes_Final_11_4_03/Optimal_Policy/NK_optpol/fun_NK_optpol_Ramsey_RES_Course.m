% ------------------------------------------------------
% This program defines the conditions needed to find the
% steady state value of hours worked using fsolve
% ------------------------------------------------------

function y=fun_NK_optpol_RES_Course(x,PIE)

global M_
 
%% DO NOT CHANGE THIS PART.
%%
%% Here we load the values of the deep parameters in a loop.
%%
%NumberOfParameters = M_.param_nbr;                            % Number of deep parameters.
%for i = 1:NumberOfParameters                                  % Loop...
%  paramname = deblank(M_.param_names(i,:));                   %    Get the name of parameter i. 
%  eval([ paramname ' = M_.params(' int2str(i) ');']);         %    Get the value of parameter i.
%end                                                           % End of the loop.  
check = 0;
%%
params = M_.params;

gy=params(1);
alp=params(2);
c=params(3);
zzeta=params(4);
delta=params(5);
sigma_c=params(6);
rhoA=params(7);
rhoG=params(8);
rhoMS=params(9);
rhoMPS=params(10);
Ass=params(11);
phiX=params(12);
xi=params(13);
hab=params(14);
hss=params(15);
g=params(16);
gamp=params(17);
PIEss=params(18);
Rnss=params(19);
rhocapqual=params(20);
ppsi=params(21);
trend=params(22);
conspie=params(23);
consr=params(24);
varrho =params(25);
betta=params(26);
wr =params(27);
alpha_r =params(28);
alpha_y=params(29);
alpha_pie=params(30);
alpha_A =params(31);
alpha_G=params(32); 
alpha_MS =params(33);
alpha_capqual=params(34);
alpha_trend=params(35);

%%Define hours as unknown for both sticky and flexi-price parts
h=exp(x(1))/(1+exp(x(1)));
hF=exp(x(2))/(1+exp(x(2)));

PIETILDE=PIE^(1-gamp);
%h=hss;
A=Ass;
Rex=(1+g)^((1-varrho)*(sigma_c-1)+1)/betta;
Rn=PIE*Rex;
DD=1/Rex;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calvo Contracts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Delta=((1-xi)*(((1-xi*PIETILDE^(zzeta-1))/(1-xi))^(1/(1-zzeta)))^(-zzeta))/(1-xi*PIETILDE^zzeta);
MC=(1-1/zzeta)*(1-xi*betta*PIETILDE^zzeta*(1+g)^((1-varrho)*(1-sigma_c)))...
/(1-xi*betta*PIETILDE^(zzeta-1)*(1+g)^((1-varrho)*(1-sigma_c)))*(((1-xi*PIETILDE^(zzeta-1))...
/(1-xi))^(1/(1-zzeta)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Rotemberg Contracts
%
%Replace with
%Delta=1;
%MC=1-(1-ppsi*(PIE-1)*PIE*(1-DD*(1+g)))/zzeta;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


PWP=MC;
KY=(1-alp)*PWP/(Rex-1+delta)*(1+g);
YW=A*h*(KY/(1+g))^((1-alp)/alp)/Delta^(1/alp);
K=KY*YW;
I=(delta+g)*K/(1+g);
Y=(1-c)*YW;
G=gy*Y;
WP=alp*PWP*YW/h;
C=WP*((1-varrho)*(1-h))/((1-hab/(1+g))*varrho);
tax=G/(WP*h);
%
%Flexi-Price 
%
RexF=Rex;
MCF=(1-1/zzeta);
PWPF=MCF;

KYF=(1-alp)*PWPF/(RexF-1+delta)*(1+g);
YWF=A*hF*(KYF/(1+g))^((1-alp)/alp);

KF=KYF*YWF;
IF=(delta+g)*KF/(1+g);
YF=(1-c)*YWF;
GF=gy*YF;
WPF=alp*PWPF*YWF/hF;
CF=WPF*((1-varrho)*(1-hF))/((1-hab/(1+g))*varrho);
taxF=GF/(WPF*hF);
%Resource constraints
%this is the equation for which fsolve will find the solutions for h and hF
y=[Y-C-I-G; YF-CF-IF-GF];



