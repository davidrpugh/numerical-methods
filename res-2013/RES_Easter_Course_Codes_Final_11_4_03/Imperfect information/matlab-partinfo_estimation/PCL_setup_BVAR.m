function [GAM, COV_P, VV, DELTA, ivarobs]=PCL_setup_BVAR( H, varobs)
% sets up parameters and calls part-info kalman filter
% developed  from notes by Prof. Joe Pearlman to 
% suit partial information RE solution in accordance with, and based on, the 
% Pearlman, Currie and Levine 1986 solution.
% 22/10/06 - Version 2 for new Riccati with 4 params instead 5 

% Copyright (C) 2006-2010 Dynare Team
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
    warning_old_state = warning;
    warning off
    dr=oo_.dr;
    TeX = options_.TeX;
    
    [PDIDPDRD,PP,CCCC,AA]=PCL_setup_core ( H, varobs);

    MSIG=disclyap_fast(CCCC, CCCC*PDIDPDRD*PP*CCCC');

    COV_P=[ PP, PP; PP, PP+MSIG]; % P0
    FL_RANK=dr.PI_FL_RANK;
    ss=size(dr.PI_ghx,1);
    dr.PI_GG=[CCCC (AA-CCCC)*(eye(ss-FL_RANK)-PDIDPDRD); zeros(ss-FL_RANK) AA*(eye(ss-FL_RANK)-PDIDPDRD)];

    GAM= [ AA*(eye(ss-FL_RANK)-PDIDPDRD) zeros(ss-FL_RANK); (AA-CCCC)*(eye(ss-FL_RANK)-PDIDPDRD),  CCCC];

    VV = [  dr.PI_TT1 dr.PI_TT2];
    nn=size(VV,1);

    %DELTA
    COV_OMEGA= COV_P( end-nn+1:end, end-nn+1:end);
    COV_YR0= VV*COV_OMEGA*VV';
    diagCovYR0=diag(COV_YR0);
    diagSqrtCovYR0=sqrt(diagCovYR0);
    DELTA=inv(diag(diagSqrtCovYR0));

    if ~isfield(options_,'varobs')
      error(' No observables for DSGE-VAR');
    else
      if isempty(options_.varobs)
        error(' No observables for DSGE-VAR');
      else
        varobs=options_.varobs;
        nobs = size(varobs,1);
      end
    end
    ivarobs=zeros(nobs,1);
    if TeX
      varobs_listTeX = [];
    end
    for i=1:nobs
      i_tmp = strmatch(varobs(i,:),M_.endo_names,'exact');
      if isempty(i_tmp)
        error (['One of the specified observed variables does not exist']) ;
      else
        ivarobs(i) = i_tmp;
        if TeX
          varobs_listTeX = strvcat(varobs_listTeX,deblank(M_.endo_names_tex(i_tmp,:)));
        end
      end
    end

