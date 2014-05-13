function [ys,check] = rbc_steadystate(ys,exe)
%input vectors, ys and exe not used. Still, they must be included because
%the call to rbc_steadystate has these vectors as arguments.

%notes on the steady state
%suppose your mod file is called fname.mod. 
%if you supply a file called fname_steadystate.mod, dynare uses this file to
%compute the steady state of the model.
%if you include the steady; command, dynare will get the steady state
%supplied by fname_steadystate and verify that the equations in
%fname_static.m are satisfied. If they are not satisfied then Dynare
%crashes, with a report of which equations are not satisfied. For this
%purpose, the equations are in the order that they appear in the model
%statements.
%if steady; is not invoked, then the steady state in fname_steadystate.mod
%is used without any check.

%if you only know how to approximately compute the steady state and you
%want dynare to refine your calculation, then don't include
%fname_steadystate.m in the directory and instead enter your guess of the
%steady state in the initval command. If the guess is good enough, dynare
%will be able to refine it. Otherwise, Dynare will give you a message
%saying that it cannot find the steady state.

global M_

check = 0;

% because some of the parameter names in this example can be confused with function values,
% it is a good idea to preset the names of the parameters to ensure they
% are not confused for functions:
alpha=0;beta=0;delta=0;

%now, read in the values of the parameters:
Np = M_.param_nbr;                                            % Number of objects in the parameter statement
for i = 1:Np                                                  %
    paramname = deblank(M_.param_names(i,:));                   % Get the name of parameter i.
    eval([ paramname ' = M_.params(' int2str(i) ');']);         % assign the value to the parameter
end                                                           %

%compute the model steady state:
k    = (1/alpha*(1/beta - (1-delta)))^(1/(alpha-1));
c     = k^alpha - delta*k;
y     =  k^alpha;
z     = 1;

err(1) =  1/c - beta/c*(alpha*y(+1)/k + 1 - delta);
err(2) =  y - (c + k - (1-delta)*k);
err(3) =  y - (k^alpha)*(z^(1-alpha));


if max(abs(err))>1e-9
    error('fatal (ssrbc1) failed to compute steady state')
end

%assign values to all the variables in the var command in the mod file
lc=log(c);
lk=log(k);
ly=log(y);
a=log(z);
a=10;
%put the steady states in the vector ys, making sure to get everything in
%the right order (i.e., the order in which they appear in the var command)
Ne=M_.orig_endo_nbr;
ys = zeros(Ne,1);
for i = 1:Ne                                
    varname = deblank(M_.endo_names(i,:));  
        eval(['ys(' int2str(i) ') = ' varname ';']);                % Get the steady state value of this variable.
end                                                        
