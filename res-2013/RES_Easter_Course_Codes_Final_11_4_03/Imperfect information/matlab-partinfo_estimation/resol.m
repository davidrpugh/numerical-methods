function [dr,info]=resol(steady_state_0,check_flag)
% function [dr,info]=resol(steady_state_0,check_flag)
% Computes first and second order approximations
%
% INPUTS
%    steady_state_0: vector of variables in steady state
%    check_flag=0:   all the approximation is computed
%    check_flag=1:   computes only the eigenvalues
%
% OUTPUTS
%    dr:             structure of decision rules for stochastic simulations
%    info=1:         the model doesn't determine the current variables '...' uniquely
%    info=2:         MJDGGES returns the following error code'
%    info=3:         Blanchard Kahn conditions are not satisfied: no stable '...' equilibrium
%    info=4:         Blanchard Kahn conditions are not satisfied:'...' indeterminacy
%    info=5:         Blanchard Kahn conditions are not satisfied:'...' indeterminacy due to rank failure
%    info=6:         The jacobian evaluated at the steady state is complex.
%    info=19:        The steadystate file did not compute the steady state (inconsistent deep parameters).
%    info=20:        can't find steady state info(2) contains sum of sqare residuals
%    info=21:        steady state is complex 
%                               info(2) contains sum of sqare of
%                               imaginary part of steady state
%    info=30:        Variance can't be computed
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2001-2011 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

global M_ options_ oo_
global it_

jacobian_flag = 0; 

if isfield(oo_,'dr');
    dr = oo_.dr;
end

options_ = set_default_option(options_,'jacobian_flag',1);
info = 0;

it_ = M_.maximum_lag + 1 ;

if M_.exo_nbr == 0
    oo_.exo_steady_state = [] ;
end

% check if steady_state_0 is steady state
tempex = oo_.exo_simul;
oo_.exo_simul = repmat(oo_.exo_steady_state',M_.maximum_lag+M_.maximum_lead+1,1);
if M_.exo_det_nbr > 0 
    tempexdet = oo_.exo_det_simul;
    oo_.exo_det_simul = repmat(oo_.exo_det_steady_state',M_.maximum_lag+M_.maximum_lead+1,1);
end
steady_state = steady_state_0;
check1 = 0;
% testing for steadystate file
if (~options_.bytecode)
    fh = str2func([M_.fname '_static']);
end;

if options_.steadystate_flag
    [steady_state,check1] = feval([M_.fname '_steadystate'],steady_state,...
                           [oo_.exo_steady_state; ...
                        oo_.exo_det_steady_state]);
    if size(steady_state,1) < M_.endo_nbr 
        if length(M_.aux_vars) > 0
            steady_state = add_auxiliary_variables_to_steadystate(steady_state,M_.aux_vars,...
                                                           M_.fname,...
                                                           oo_.exo_steady_state,...
                                                           oo_.exo_det_steady_state,...
                                                           M_.params,...
                                                           options_.bytecode);
        else
            error([M_.fname '_steadystate.m doesn''t match the model']);
        end
    end

else
    % testing if steady_state_0 isn't a steady state or if we aren't computing Ramsey policy
    if  options_.ramsey_policy == 0
        if options_.linear == 0
            % nonlinear models
            if (options_.block == 0 && options_.bytecode == 0)
                if max(abs(feval(fh,steady_state,[oo_.exo_steady_state; ...
                                        oo_.exo_det_steady_state], M_.params))) > options_.dynatol
                    [steady_state,check1] = dynare_solve(fh,steady_state,options_.jacobian_flag,...
                                                  [oo_.exo_steady_state; ...
                                        oo_.exo_det_steady_state], M_.params);
                end
            else
                [steady_state,check1] = dynare_solve_block_or_bytecode(steady_state,...
                                                                [oo_.exo_steady_state; ...
                                    oo_.exo_det_steady_state], M_.params);
            end;
        else
            % linear models
            [fvec,jacob] = feval(fh,steady_state,[oo_.exo_steady_state;...
                                oo_.exo_det_steady_state], M_.params);
            if max(abs(fvec)) > 1e-12
                steady_state = steady_state-jacob\fvec;
            end
        end
    end
end
% testing for problem
dr.ys = steady_state;

if check1
    if options_.steadystate_flag
        info(1)= 19;
        resid = check1 ;
    else
        info(1)= 20;
        resid = feval(fh,steady_state_0,oo_.exo_steady_state, M_.params);
    end
    info(2) = resid'*resid ;
    return
end

if ~isreal(steady_state)
    info(1) = 21;
    info(2) = sum(imag(steady_state).^2);
    steady_state = real(steady_state);
    dr.ys = steady_state;
    return
end

if options_.block
    [dr,info,M_,options_,oo_] = dr_block(dr,check_flag,M_,options_,oo_);
elseif((options_.usePartInfo==1)  | (options_.partial_information ==1) | (options_.useACES==1))%&& (check_flag == 0)
    [dr,info,M_,options_,oo_] = dr1_PI(dr,check_flag,M_,options_,oo_);
else
    [dr,info,M_,options_,oo_] = dr1(dr,check_flag,M_,options_,oo_);
end
if info(1)
    return
end

if M_.exo_det_nbr > 0
    oo_.exo_det_simul = tempexdet;
end
oo_.exo_simul = tempex;
tempex = [];

% 01/01/2003 MJ added dr_algo == 1
% 08/24/2001 MJ uses Schmitt-Grohe and Uribe (2001) constant correction
%               in dr.ghs2 
% 05/26/2003 MJ added temporary values for oo_.exo_simul
