% ------------------------------------------------------
% This program defines the conditions needed to find the
% steady state value of hours worked using fsolve
% ------------------------------------------------------

function y=fun_GK_RES_Course(x,PIE)

global M_
 
%% DO NOT CHANGE THIS PART.
%%
%% Here we load the values of the deep parameters in a loop.
%%
NumberOfParameters = M_.param_nbr;                            % Number of deep parameters.
for i = 1:NumberOfParameters                                  % Loop...
  paramname = deblank(M_.param_names(i,:));                   %    Get the name of parameter i. 
  eval([ paramname ' = M_.params(' int2str(i) ');']);         %    Get the value of parameter i.
end                                                           % End of the loop.  
check = 0;
%%

%%Define hours as unknown
varrho=exp(x(1))/(1+exp(x(1)));
K=x(2);
ThetaB=exp(x(3))/(1+exp(x(3)));
xiB=x(4);
hF=exp(x(5))/(1+exp(x(5)));
betta=exp(x(6))/(1+exp(x(6)));


%Define only the steady state relationships needed to find the value of hours worked that supports
%the resource constraint relationship Y=C+I+G
PIETILDE=PIE^(1-gamp);
Delta=((1-xi)*(((1-xi*PIETILDE^(zzeta-1))/(1-xi))^(1/(1-zzeta)))^(-zzeta))/(1-xi*PIETILDE^zzeta);
h=hss;
A=Ass;
Q=1;
Rn=Rnss;
Rex=(1+g)^((1-varrho)*(sigma_c-1)+1)/betta;
DD=1/Rex;
MC=(1-1/zzeta)*(1-xi*betta*PIETILDE^zzeta*(1+g)^((1-varrho)*(1-sigma_c)))...
/(1-xi*betta*PIETILDE^(zzeta-1)*(1+g)^((1-varrho)*(1-sigma_c)))*(((1-xi*PIETILDE^(zzeta-1))...
/(1-xi))^(1/(1-zzeta)));
PWP=MC;

YW=A*h^(alp)*(K/(1+g))^(1-alp)/Delta;
I=(delta+g)*K/(1+g);
Y=(1-c)*YW;
G=gy*Y;
WP=alp*PWP*YW/h;
C=WP*((1-varrho)*(1-h))/((1-hab/(1+g))*varrho);
tax=G/(WP*h);

S=K;
Z=(1-alp)*PWP*YW/(K/(1+g));
phiB=lev;
NW=(Q*S)/phiB;
Dep=(Q*S)-NW;
omega=1-sigmaB+sigmaB*ThetaB*phiB;
nuB=omega;
muB=(phiB*ThetaB-nuB)/phiB;
Rk=((muB)/(DD*omega)+(Rex));
spread=Rk-Rex;
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

%Resource constraint 
%this is the equation for which fsolve will find the solution
y=[Y-C-I-G; YF-CF-IF-GF; 
   Z+(1-delta)-Rk; 
   creditspread-spread;
 NW-((((sigmaB+xiB)*(Z+(1-delta))*Q*S/(1+g))-(sigmaB*(Rex)*Q*S/(1+g)))/(1-sigmaB*(Rex)/(1+g)));
Rex-Rn/PIE;];



