function [xfeedback,objval] = optwelrule(welvar,nf,lb,ub)
% uses the fmincon function to find an optimized simple rule (and the expected welfare) 
% in dynare wrt the feedback parameters
% developed by B Yang, Jan 2013
% 08/01/13 - Version 1
% Copyright (C) - CIMS


global M_ oo_

fh = str2func('stoch_simul_wel');
optim_options = optimset('display','iter','LargeScale','off', ...
                                 'MaxFunEvals',100000,'TolFun',1e-8,'TolX',1e-6);
x0=[];
param_nbr=M_.param_nbr;
param_names=M_.param_names;
fb_i=[param_nbr-nf+1:param_nbr];
fbnames=param_names(fb_i,:);

if isempty(lb)
    lb=-inf(1,nf);
end    
if isempty(ub)
    ub=inf(1,nf);
end 

%find the index of omega and it's deterministic steady state
welvar_i = strmatch(welvar, M_.endo_names, 'exact'); 
determwel = oo_.steady_state(welvar_i);

%x0=[0.7;5;0.25];

%use the feedback parameters imposed to initilise
temp=M_.params(end:-1:1);
for i=1:nf
    x0(nf+1-i)=temp(i);
end    
x0=x0';

%funciton value at initial parameters 
fval=stoch_simul_wel(x0,welvar);
disp(' ')
disp('Compute welfare-optimized rules...')
disp(' ')
title='';
    labels = char('Starting values');
    headers = char('',fbnames);
    lh = size(labels,2)+2;
    dyntable(title,headers,labels,x0',lh,8,4);
disp(sprintf('Intitial ojective function: %f',-fval));

%optimal feedbacks and MINUS objective
[xfeedback,objval]=fmincon(fh,x0,[],[],[],[],lb,ub,[],optim_options,welvar);
disp(' ')
disp(sprintf('Objective function: %f',-objval));

w=[-objval determwel];
 title='WELFARE FUNCTION';
    headers=char('','Expected','Deterministic');
    labels = deblank(welvar);
    lh = size(labels,2)+2;
    dyntable(title,headers,labels,w,lh,8,4);
if determwel==0
    disp(' ')
    warning('Zero deterministic function: need to find the deterministic steady state of the model'); 
end    

 title='OPTIMIZED RULE';
    labels = char('Coeff.');
    headers = char('',fbnames);
    lh = size(labels,2)+2;
    dyntable(title,headers,labels,xfeedback',lh,8,4);