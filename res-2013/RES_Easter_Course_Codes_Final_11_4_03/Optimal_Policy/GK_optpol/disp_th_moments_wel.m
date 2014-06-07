function [m sd] = disp_th_moments_wel(dr,var_list)
% Display theoretical moments of variables

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

global M_ oo_ options_

if size(var_list,1) == 0
    var_list = M_.endo_names(1:M_.orig_endo_nbr, :);
end
nvar = size(var_list,1);
ivar=zeros(nvar,1);
for i=1:nvar
    i_tmp = strmatch(var_list(i,:),M_.endo_names,'exact');
    if isempty(i_tmp)
        error (['One of the variable specified does not exist']) ;
    else
        ivar(i) = i_tmp;
    end
end

[oo_.gamma_y,stationary_vars] = th_autocovariances(dr,ivar,M_,options_);
m = dr.ys(ivar);
non_stationary_vars = setdiff(1:length(ivar),stationary_vars);
m(non_stationary_vars) = NaN;

i1 = find(abs(diag(oo_.gamma_y{1})) > 1e-12);
s2 = diag(oo_.gamma_y{1});
sd = sqrt(s2);
% in the case for welfare 2nd order is always needed
%if options_.order == 2
    m = m+oo_.gamma_y{options_.ar+3};
%end