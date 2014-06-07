%compute the steady state of RBC Model 
function[ys,check]=RBC_Course_steadystate(ys,exe)
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


%% THIS BLOCK IS MODEL SPECIFIC.
%%
%% Here the user has to define the steady state.
%%



 %Uses Analytical Solution without calling fun_RBC.m
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


A=Ass;
R=1.0/betta;
PWP=1-1/zzeta;
KYW=PWP*(1-alp)/(R-1+delta);
IY=(1-alp)*delta/(R-1+delta);
CY=1-IY-gy;
h=(1-varrho)*alp/(varrho*CY+(1-varrho)*alp);
YW=A*h*KYW^((1-alp)/alp);
Y=(1-c)*YW;
K=KYW*YW;
I=IY*Y;
G=gy*Y;
C=Y-I-G;
WP=alp*PWP*YW/h;
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