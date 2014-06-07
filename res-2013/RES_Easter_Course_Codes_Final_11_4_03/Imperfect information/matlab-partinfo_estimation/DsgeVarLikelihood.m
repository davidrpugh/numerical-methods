function [fval,cost_flag,info,PHI,SIGMAu,iXX,prior] = DsgeVarLikelihood(xparam1,gend)
% function [fval,cost_flag,info,PHI,SIGMAu,iXX] = DsgeVarLikelihood(xparam1,gend)
% Evaluates the posterior kernel of the bvar-dsge model. 
% 
% INPUTS 
%   o xparam1       [double]     Vector of model's parameters.
%   o gend          [integer]    Number of observations (without conditionning observations for the lags).
%  
% OUTPUTS 
%   o fval          [double]     Value of the posterior kernel at xparam1.
%   o cost_flag     [integer]    Zero if the function returns a penalty, one otherwise.
%   o info          [integer]    Vector of informations about the penalty.
%   o PHI           [double]     Stacked BVAR-DSGE autoregressive matrices (at the mode associated to xparam1).
%   o SIGMAu        [double]     Covariance matrix of the BVAR-DSGE (at the mode associated to xparam1).
%   o iXX           [double]     inv(X'X).
%   o prior         [double]     a matlab structure describing the dsge-var prior.  
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 2006-2008 Dynare Team
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

global bayestopt_ estim_params_ M_ options_

nvx = estim_params_.nvx;
nvn = estim_params_.nvn;
ncx = estim_params_.ncx;
ncn = estim_params_.ncn;
np  = estim_params_.np;
nx = nvx+nvn+ncx+ncn+np;
ns = nvx+nvn+ncx+ncn;

NumberOfObservedVariables = size(options_.varobs,1);
NumberOfLags = options_.dsge_varlag;
NumberOfParameters = NumberOfObservedVariables*NumberOfLags ;
if ~options_.noconstant
    NumberOfParameters = NumberOfParameters + 1;
end

mYY = evalin('base', 'mYY');
mYX = evalin('base', 'mYX');
mXY = evalin('base', 'mXY');
mXX = evalin('base', 'mXX');

fval = [];
cost_flag = 1;

if options_.mode_compute ~= 1 & any(xparam1 < bayestopt_.lb)
    k = find(xparam1 < bayestopt_.lb);
    fval = bayestopt_.penalty+sum((bayestopt_.lb(k)-xparam1(k)).^2);
    cost_flag = 0;
    info = 41;
    return;
end

if options_.mode_compute ~= 1 & any(xparam1 > bayestopt_.ub)
    k = find(xparam1 > bayestopt_.ub);
    fval = bayestopt_.penalty+sum((xparam1(k)-bayestopt_.ub(k)).^2);
    cost_flag = 0;
    info = 42;
    return;
end

Q = M_.Sigma_e;
for i=1:estim_params_.nvx
    k = estim_params_.var_exo(i,1);
    Q(k,k) = xparam1(i)*xparam1(i);
end
offset = estim_params_.nvx;
if estim_params_.nvn
    disp('DsgeVarLikelihood :: Measurement errors are implemented!')
    return
end 
if estim_params_.ncx
    disp('DsgeVarLikelihood :: Correlated structural innovations are not implemented!')
    return
end

M_.params(estim_params_.param_vals(:,1)) = xparam1(offset+1:end);
M_.Sigma_e = Q;

%% Weight of the dsge prior:
dsge_prior_weight = M_.params(strmatch('dsge_prior_weight',M_.param_names));
% Is the DSGE prior proper?
if dsge_prior_weight<(NumberOfParameters+NumberOfObservedVariables)/gend;
    fval = bayestopt_.penalty+abs(gend*dsge_prior_weight-(NumberOfParameters+NumberOfObservedVariables));
    cost_flag = 0;
    info = 51;
    return;
end


%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------
[T,R,SteadyState,info] = dynare_resolve('restrict');

if info(1) == 1 || info(1) == 2 || info(1) == 5
    fval = bayestopt_.penalty+1;
    cost_flag = 0;
    return
elseif info(1) == 3 || info(1) == 4 || info(1) == 19 || info(1) == 20 || info(1) == 21
    fval = bayestopt_.penalty+info(2);
    cost_flag = 0;
    return
end

if ~options_.noconstant
    if options_.loglinear 
        constant = transpose(log(SteadyState(bayestopt_.mfys)));
    else
        constant = transpose(SteadyState(bayestopt_.mfys));
    end 
else
    constant = zeros(1,NumberOfObservedVariables);
end
if bayestopt_.with_trend == 1
    disp('DsgeVarLikelihood :: Linear trend is not yet implemented!')
    return
end

%------------------------------------------------------------------------------
% 3. theoretical moments (second order)
%------------------------------------------------------------------------------
if ~(options_.usePartInfo==1  | (options_.partial_information ==1))

tmp0 = lyapunov_symm(T,R*Q*R',options_.qz_criterium,options_.lyapunov_complex_threshold);% I compute the variance-covariance matrix
mf  = bayestopt_.mf1;          % of the restricted state vector.

% Get the non centered second order moments
TheoreticalAutoCovarianceOfTheObservedVariables = ...
    zeros(NumberOfObservedVariables,NumberOfObservedVariables,NumberOfLags+1);
TheoreticalAutoCovarianceOfTheObservedVariables(:,:,1) = tmp0(mf,mf)+constant'*constant;
for lag = 1:NumberOfLags
    tmp0 = T*tmp0;
    TheoreticalAutoCovarianceOfTheObservedVariables(:,:,lag+1) = tmp0(mf,mf) ...
        + constant'*constant;
end

else % if(options_.usePartInfo==1  | (options_.partial_information ==1))
  [GAM, COV_P, VV, DELTA, ivarobs]=PCL_setup_BVAR ( M_.H, options_.varobs);
  nn=size(VV,1);
  options_ = set_default_option(options_,'ar',5);
  COV_YRk= zeros(nn,nn);
  COV_OMEGA= COV_P( end-nn+1:end, end-nn+1:end);
  COV_YRk = VV*COV_OMEGA*VV';
  TheoreticalAutoCovarianceOfTheObservedVariables = ...
    zeros(NumberOfObservedVariables,NumberOfObservedVariables,NumberOfLags+1);
  TheoreticalAutoCovarianceOfTheObservedVariables(:,:,1) = COV_YRk(ivarobs,ivarobs)+constant'*constant;
  for lag = 1:NumberOfLags
    COV_P=GAM*COV_P;
    COV_OMEGA= COV_P( end-nn+1:end, end-nn+1:end);
    COV_YRk = VV*COV_OMEGA*VV';
    %AutoCOR_YRkMAT=DELTA*COV_YRk*DELTA;
    %TheoreticalAutoCovarianceOfTheObservedVariables(:,:,lag+1) = AutoCOR_YRkMAT(ivarobs,ivarobs)+constant'*constant;
    TheoreticalAutoCovarianceOfTheObservedVariables(:,:,lag+1) = COV_YRk(ivarobs,ivarobs)+constant'*constant;
  end
end
% Build the theoretical "covariance" between Y and X
GYX = zeros(NumberOfObservedVariables,NumberOfParameters);
for i=1:NumberOfLags
    GYX(:,(i-1)*NumberOfObservedVariables+1:i*NumberOfObservedVariables) = ...
        TheoreticalAutoCovarianceOfTheObservedVariables(:,:,i+1);
end
if ~options_.noconstant
    GYX(:,end) = constant';
end
% Build the theoretical "covariance" between X and X
GXX = kron(eye(NumberOfLags), ...
           TheoreticalAutoCovarianceOfTheObservedVariables(:,:,1));
for i = 1:NumberOfLags-1
    tmp1 = diag(ones(NumberOfLags-i,1),i); 
    tmp2 = diag(ones(NumberOfLags-i,1),-i);
    GXX = GXX + kron(tmp1,TheoreticalAutoCovarianceOfTheObservedVariables(:,:,i+1));
    GXX = GXX + kron(tmp2,TheoreticalAutoCovarianceOfTheObservedVariables(:,:,i+1)');
end

if ~options_.noconstant
    % Add one row and one column to GXX
    GXX = [GXX , kron(ones(NumberOfLags,1),constant') ; [  kron(ones(1,NumberOfLags),constant) , 1] ];
end

GYY = TheoreticalAutoCovarianceOfTheObservedVariables(:,:,1);

assignin('base','GYY',GYY);
assignin('base','GXX',GXX);
assignin('base','GYX',GYX);

if ~isinf(dsge_prior_weight)
    tmp0 = dsge_prior_weight*gend*TheoreticalAutoCovarianceOfTheObservedVariables(:,:,1) + mYY ;
    tmp1 = dsge_prior_weight*gend*GYX + mYX;
    tmp2 = inv(dsge_prior_weight*gend*GXX+mXX);
    SIGMAu = tmp0 - tmp1*tmp2*tmp1'; clear('tmp0');
    if ~ispd(SIGMAu)
        v = diag(SIGMAu);
        k = find(v<0);
        fval = bayestopt_.penalty + sum(v(k).^2);
        info = 52;
        cost_flag = 0;
        return;
    end
    SIGMAu = SIGMAu / (gend*(1+dsge_prior_weight));
    PHI = tmp2*tmp1'; clear('tmp1');
    prodlng1 = sum(gammaln(.5*((1+dsge_prior_weight)*gend- ...
                               NumberOfObservedVariables*NumberOfLags ...
                               +1-(1:NumberOfObservedVariables)')));
    prodlng2 = sum(gammaln(.5*(dsge_prior_weight*gend- ...
                               NumberOfObservedVariables*NumberOfLags ...
                               +1-(1:NumberOfObservedVariables)')));  
    lik = .5*NumberOfObservedVariables*log(det(dsge_prior_weight*gend*GXX+mXX)) ...
          + .5*((dsge_prior_weight+1)*gend-NumberOfParameters)*log(det((dsge_prior_weight+1)*gend*SIGMAu)) ...
          - .5*NumberOfObservedVariables*log(det(dsge_prior_weight*gend*GXX)) ...
          - .5*(dsge_prior_weight*gend-NumberOfParameters)*log(det(dsge_prior_weight*gend*(GYY-GYX*inv(GXX)*GYX'))) ...
          + .5*NumberOfObservedVariables*gend*log(2*pi)  ...
          - .5*log(2)*NumberOfObservedVariables*((dsge_prior_weight+1)*gend-NumberOfParameters) ...
          + .5*log(2)*NumberOfObservedVariables*(dsge_prior_weight*gend-NumberOfParameters) ...
          - prodlng1 + prodlng2;
else
    iGXX = inv(GXX);
    SIGMAu = GYY - GYX*iGXX*transpose(GYX);
    PHI = iGXX*transpose(GYX);
    lik = gend * ( log(det(SIGMAu)) + NumberOfObservedVariables*log(2*pi) +  ...
                   trace(inv(SIGMAu)*(mYY - transpose(mYX*PHI) - mYX*PHI + transpose(PHI)*mXX*PHI)/gend));
    lik = .5*lik;% Minus likelihood
end      

lnprior = priordens(xparam1,bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4);
fval = (lik-lnprior);

if (nargout == 6)
    if isinf(dsge_prior_weight)
        iXX = iGXX;
    else
        iXX = tmp2;
    end
end

if (nargout==7)
    if isinf(dsge_prior_weight)
        iXX = iGXX;
    else
        iXX = tmp2;
    end
    iGXX = inv(GXX); 
    prior.SIGMAstar = GYY - GYX*iGXX*GYX';
    prior.PHIstar = iGXX*transpose(GYX);
    prior.ArtificialSampleSize = fix(dsge_prior_weight*gend);
    prior.DF = prior.ArtificialSampleSize - NumberOfParameters - NumberOfObservedVariables;
    prior.iGXX = iGXX;
end
