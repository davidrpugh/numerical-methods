% ------------------------------------------------------
% This program defines the conditions needed to find the
% steady state value of hours worked using fsolve
% ------------------------------------------------------

function y=fun_KM_RES_Course(x,PIE)

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

%%Derived variables
varrho=exp(x(1))/(1+exp(x(1)));
hF=exp(x(2))/(1+exp(x(2)));
K=exp(x(3));
bettaE=exp(x(4))/(1+exp(x(4)));
TE=x(5);
betta=exp(x(6))/(1+exp(x(6)));
PIETILDE=PIE^(1-gamp);
Delta=((1-xi)*(((1-xi*PIETILDE^(zzeta-1))/(1-xi))^(1/(1-zzeta)))^(-zzeta))/(1-xi*PIETILDE^zzeta);
h=hss;
A=Ass;
Rex=(1+g)^((1-varrho)*(sigma_c-1)+1)/betta;
Rn=Rex*PIE;
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


Z=(1-alp)*PWP*YW/(K/(1+g));
Q=1;
DDE=bettaE/(1+g)^sigma_entre;
RL=1/(1-1/zzetaL)*Rn;
ThetaE=1/RL-DDE/PIE;
Rk=(Z+(1-delta)*Q)/Q;
RLex=RL/PIE;
L=m*PIE*Q*(1-delta)*K/RL;
CE=(Rk/(1+g)-1)*Q*K+L*(1-RLex/(1+g))-TE;
spread=Rk-Rex;
%
%this is the equation for which fsolve will find the solution
%
y=[Y-C-CE-I-G; YF-CF-IF-GF;...
Rk-(1-ThetaE*m*PIE*Q*(1-delta))/DDE;...
spread-spreadcalib; CE/Y-CEcalib; Rn-Rnss];



