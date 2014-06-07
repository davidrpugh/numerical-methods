function [fval,th_mean, info] = stoch_simul_wel(x0,var_list)

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

global M_ options_ oo_ it_

%impose the feedback parameters
temp=M_.params(end:-1:1);
for i=1:size(x0,1)
    temp(i)=x0(size(x0,1)+1-i);
end    
M_.params=temp(end:-1:1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_for_deep_parameters_calibration(M_);

options_old = options_;
if options_.linear
    options_.order = 1;
end
if options_.order == 1
    options_.replic = 1;
elseif options_.order == 3
    options_.k_order_solver = 1;
end

if isempty(options_.qz_criterium)
    options_.qz_criterium = 1+1e-6;
end

%if options_.partial_information == 1 || options_.ACES_solver == 1
%    PI_PCL_solver = 1;
%    if options_.order ~= 1
%        warning('STOCH_SIMUL: forcing order=1 since you are using partial_information or ACES solver')
%        options_.order = 1;
%    end
%else
    PI_PCL_solver = 0;
%end

TeX = options_.TeX;

if size(var_list,1) == 0
    var_list = M_.endo_names(1:M_.orig_endo_nbr, :);
end

[i_var,nvar] = varlist_indices(var_list,M_.endo_names);

iter_ = max(options_.periods,1);
if M_.exo_nbr > 0
    oo_.exo_simul= ones(iter_ + M_.maximum_lag + M_.maximum_lead,1) * oo_.exo_steady_state';
end

check_model;

if PI_PCL_solver
    [oo_.dr, info] = PCL_resol(oo_.steady_state,0);
else
    [oo_.dr, info] = resol(oo_.steady_state,0);
end

if info(1)
    options_ = options_old;
   % print_info(info, options_.noprint);
   % return
end

if options_.periods > 0 && ~PI_PCL_solver
    if options_.periods <= options_.drop
        disp(['STOCH_SIMUL error: The horizon of simulation is shorter' ...
              ' than the number of observations to be DROPed'])
        options_ =options_old;
        return
    end
    oo_.endo_simul = simult(oo_.dr.ys,oo_.dr);
    dyn2vec;
end

if options_.nomoments == 0
    if PI_PCL_solver
        PCL_Part_info_moments (0, PCL_varobs, oo_.dr, i_var);
    elseif options_.periods == 0
        [th_mean, th_sd] = disp_th_moments_wel(oo_.dr,var_list); 
        if size(var_list,1)==1
           fval=-th_mean;
        else
           fval=th_sd(2,:);
        end
    else
        disp_moments(oo_.endo_simul,var_list);
    end
end

