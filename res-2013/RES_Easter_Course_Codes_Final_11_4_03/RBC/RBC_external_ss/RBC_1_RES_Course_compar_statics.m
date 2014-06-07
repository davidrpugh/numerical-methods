%
% Comparative statics program
%
%compute the steady state of RBC Model 
%function[ys,check]=RBC_Course_steadystate(ys,exe)
%
% global statement needed if fun and fsolve used; repeat statement in fun
%
global gy alp zzeta c betta delta sigma_c varrho Ass
 
%% DO NOT CHANGE THIS PART.
%%
%% Here we load the values of the deep parameters in a loop.
%%
%NumberOfParameters = M_.param_nbr;                            % Number of deep parameters.
%for i = 1:NumberOfParameters                                  % Loop...
%  paramname = deblank(M_.param_names(i,:));                   %    Get the name of parameter i. 
%  eval([ paramname ' = M_.params(' int2str(i) ');']);         %    Get the value of parameter i.
%end                                                           % End of the loop.  
%check = 0;
%%


%% THIS BLOCK IS MODEL SPECIFIC.
%%
%% Here the user has to define the steady state.
%%
%
%For use of steady state program as a stand-alone replace all above with:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PARAMETER SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gy=0.2;
alp=0.70;
zzeta=7.0;
%zzeta=1000000000;
c=1/zzeta;
betta=0.9871;
delta=0.0250;
sigma_c=2.0;
varrho=0.6073;

%
%
%Choice of Units
%
Ass=1;
%
%shock persistence
%
rhoA=0.7;
rhoG=0.7;

%
%loop for comparative stattics
%
N=200;
for i=1:N
    varrho=0.5+(i-1)*0.001;
    varrhovec(i)=varrho;
 %   
 %Uses Analytical Solution without calling fun_RBC.m
 %
 
A=Ass;
R=1.0/betta-1;
PWP=1-1/zzeta;
KYW=PWP*(1-alp)/(R+delta);
IYW=delta*KYW;
IY=IYW/(1-c);
CY=1-IY-gy;
X=(1-varrho)*alp/(varrho*CY);
h=X/(1+X);
%
hvec(i)=h;
% 
YW=A*h*KYW^((1-alp)/alp);
WP=alp*PWP*YW/h;
K=KYW*YW;
I=(delta)*K;
Y=(1-c)*YW;
G=gy*Y;
C=Y-I-G;
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

end
figure
plot(hvec, varrhovec);
xlabel('hours');
ylabel('Preference Parameter \rho');
%% END OF THE MODEL SPECIFIC BLOCK.
%
%rest not needed for comparative ststics
%
%
%% DO NOT CHANGE THIS PART.
%%
%% Here we define the steady state vZNues of the endogenous variables of
%% the model.
%%
%NumberOfEndogenousVariables = M_.endo_nbr;                    % Number of endogenous variables.
%ys = zeros(NumberOfEndogenousVariables,1);                    % Initialization of ys (steady state).
%for i = 1:NumberOfEndogenousVariables                         % Loop...
%  varname = deblank(M_.endo_names(i,:));                      %    Get the name of endogenous variable i.                     
%  eval(['ys(' int2str(i) ') = ' varname ';']);                %    Get the steady state vZNue of this variable.

%end                                                           % End of the loop.
%%