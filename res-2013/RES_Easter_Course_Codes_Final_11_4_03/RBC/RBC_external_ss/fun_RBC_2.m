% ------------------------------------------------------
% This program defines the conditions needed to find the
% steady state value of hours worked using fsolve
% ------------------------------------------------------

function y=fun_RBC(x)

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

%%Define varrho as unknown endogenous variable
varrho=exp(x(1))/(1+exp(x(1)));
%Define only the steady state relationships needed to find the value of varrho to satisfy C=WP*((1-varrho)*(1-h))/varrho
A=Ass;
h=hobs;
R=1.0/betta;
PWP=1-1/zzeta;
KY=(1-alp)*PWP/(R-1+delta);
YW=A*h*KY^((1-alp)/alp);
K=KY*YW;
I=(delta)*K;
Y=(1-c)*YW;
G=gy*Y;
WP=alp*PWP*YW/h;
C=Y-I-G;
tax=G/(WP*h);
%Resource constraint 
%this is the equation for which fsolve will find the solution
y=[C-WP*((1-varrho)*(1-h))/varrho];



