%compute the steady state of RBC Model 
function[ys,check]=RBC_3_Course_steadystate(ys,exe)
global M_
 
%% DO NOT CHANGE THIS PART.
%%
%% Here we load the values of the deep parameters in a loop.
%%
NumberOfParameters = M_.param_nbr;                            % Number of deep parameters.
for i = 1:NumberOfParameters                                  % Loop...
  paramname = deblank(M_.param_names(i,:)); %    Get the name of parameter i.
 eval([ paramname ' = M_.params(' int2str(i) ');']);   %    Get the value of parameter i.
end                                                           % End of the loop.  
check = 0;
%%


%% THIS BLOCK IS MODEL SPECIFIC.
%%
%% Here the user has to define the steady state.
%%



 %This part calls the function fun_RBC.m
 %Use fsolve to find the steady state values of hours worked
 %initial values
 x0=[0.5 0.0 2.0];
 [x,fval] =fsolve(@fun_RBC_4,x0,optimset('Display','off'));
 
%Derived variables (cut and paste from fun_RBC.m)
varrho=exp(x(1))/(1+exp(x(1)));
delta=exp(x(2))/(1+exp(x(2)));
betta=exp(x(3))/(1+exp(x(3)));
A=Ass;
h=hobs;
R=1.0/betta;
PWP=1-1/zzeta;
KY=(1-alp)*PWP/(R-1+delta);
YW=A*h*KY^((1-alp)/alp);
K=KY*YW;
Y=(1-c)*YW;
I=iobs*Y;
G=gy*Y;
WP=alp*PWP*YW/h;
C=WP*((1-varrho)*(1-h))/varrho;
tax=G/(WP*h);
%Post recursive Steady state relationship
LAMBDA=1/(1-sigma_c)*(C^((1-varrho)*(1-sigma_c))*(1-h)^(varrho*(1-sigma_c))-1);
LAMBDAC=(1-varrho)*C^((1-varrho)*(1-sigma_c)-1)*(1-h)^(varrho*(1-sigma_c));

YY=1;
CC=1;
hh=1;
WPWP=1; 
II=1; 
KK=1;
RR=1;
%% END OF THE MODEL SPECIFIC BLOCK.


%% DO NOT CHANGE THIS PART.
%%
%% Here we define the steady state vZNues of the endogenous variables of
%% the model.
%%
NumberOfEndogenousVariables = M_.endo_nbr;                    % Number of endogenous variables.
ys = zeros(NumberOfEndogenousVariables,1);                    % Initialization of ys (steady state).
for i = 1:NumberOfEndogenousVariables                         % Loop...
  varname = deblank(M_.endo_names(i,:));                      %    Get the name of endogenous variable i.                     
  eval(['ys(' int2str(i) ') = ' varname ';']);                %    Get the steady state vZNue of this variable.
end                                                           % End of the loop.
%%