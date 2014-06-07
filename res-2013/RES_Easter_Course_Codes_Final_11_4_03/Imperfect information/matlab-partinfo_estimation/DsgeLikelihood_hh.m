function [fval,llik,cost_flag,ys,trend_coeff,info] = DsgeLikelihood_hh(xparam1,gend,data,data_index,number_of_observations,no_more_missing_observations)
% function [fval,cost_flag,ys,trend_coeff,info] = DsgeLikelihood(xparam1,gend,data,data_index,number_of_observations,no_more_missing_observations)
% Evaluates the posterior kernel of a dsge model. 
% 
% INPUTS 
%   xparam1                        [double]   vector of model parameters.
%   gend                           [integer]  scalar specifying the number of observations.
%   data                           [double]   matrix of data
%   data_index                     [cell]     cell of column vectors
%   number_of_observations         [integer]
%   no_more_missing_observations   [integer] 
% OUTPUTS 
%   fval        :     value of the posterior kernel at xparam1.
%   cost_flag   :     zero if the function returns a penalty, one otherwise.
%   ys          :     steady state of original endogenous variables
%   trend_coeff :
%   info        :     vector of informations about the penalty:
%                     41: one (many) parameter(s) do(es) not satisfied the lower bound
%                     42: one (many) parameter(s) do(es) not satisfied the upper bound
%               
% SPECIAL REQUIREMENTS
%

% Copyright (C) 2004-2009 Dynare Team
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

  global bayestopt_ estim_params_ options_ trend_coeff_ M_ oo_
  fval		= [];
  ys		= [];
  trend_coeff	= [];
  llik = NaN;
  cost_flag  	= 1;
  nobs 		= size(options_.varobs,1);
  llik=NaN; 
  %------------------------------------------------------------------------------
  % 1. Get the structural parameters & define penalties
  %------------------------------------------------------------------------------
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
  H = M_.H;
  for i=1:estim_params_.nvx
    k =estim_params_.var_exo(i,1);
    Q(k,k) = xparam1(i)*xparam1(i);
  end
  offset = estim_params_.nvx;
  if estim_params_.nvn
    for i=1:estim_params_.nvn
      k = estim_params_.var_endo(i,1);
      H(k,k) = xparam1(i+offset)*xparam1(i+offset);
    end
    offset = offset+estim_params_.nvn;
  end	
  if estim_params_.ncx
    for i=1:estim_params_.ncx
      k1 =estim_params_.corrx(i,1);
      k2 =estim_params_.corrx(i,2);
      Q(k1,k2) = xparam1(i+offset)*sqrt(Q(k1,k1)*Q(k2,k2));
      Q(k2,k1) = Q(k1,k2);
    end
    [CholQ,testQ] = chol(Q);
    if testQ 	%% The variance-covariance matrix of the structural innovations is not definite positive.
		%% We have to compute the eigenvalues of this matrix in order to build the penalty.
		a = diag(eig(Q));
		k = find(a < 0);
		if k > 0
		  fval = bayestopt_.penalty+sum(-a(k));
		  cost_flag = 0;
		  info = 43;
		  return
		end
    end
    offset = offset+estim_params_.ncx;
  end
  if estim_params_.ncn 
    for i=1:estim_params_.ncn
      k1 = options_.lgyidx2varobs(estim_params_.corrn(i,1));
      k2 = options_.lgyidx2varobs(estim_params_.corrn(i,2));
      H(k1,k2) = xparam1(i+offset)*sqrt(H(k1,k1)*H(k2,k2));
      H(k2,k1) = H(k1,k2);
    end
    [CholH,testH] = chol(H);
    if testH
      a = diag(eig(H));
      k = find(a < 0);
      if k > 0
	fval = bayestopt_.penalty+sum(-a(k));
	cost_flag = 0;
	info = 44;
	return
      end
    end
    offset = offset+estim_params_.ncn;
  end
  if estim_params_.np > 0
      M_.params(estim_params_.param_vals(:,1)) = xparam1(offset+1:end);
  end
  M_.Sigma_e = Q;
  M_.H = H;
  %------------------------------------------------------------------------------
  % 2. call model setup & reduction program
  %------------------------------------------------------------------------------
  [T,R,SteadyState,info] = dynare_resolve('restrict');
  if info(1) == 1 || info(1) == 2 || info(1) == 5
    fval = bayestopt_.penalty+1;
    cost_flag = 0;
    return
  elseif info(1) == 3 || info(1) == 4 || info(1)==6 ||info(1) == 19 || info(1) == 20 || info(1) == 21
    fval = bayestopt_.penalty+info(2);%^2; % penalty power raised in DR1.m and resol already. GP July'08
    cost_flag = 0;
    return
  end
  bayestopt_.mf = bayestopt_.mf1;
  if ~options_.noconstant
    if options_.loglinear == 1
      constant = log(SteadyState(bayestopt_.mfys));
    else
      constant = SteadyState(bayestopt_.mfys);
    end
  else
    constant = zeros(nobs,1);
  end
  if bayestopt_.with_trend == 1
    trend_coeff = zeros(nobs,1);
    t = options_.trend_coeffs;
    for i=1:length(t)
      if ~isempty(t{i})
	trend_coeff(i) = evalin('base',t{i});
      end
    end
    trend = repmat(constant,1,gend)+trend_coeff*[1:gend];
  else
    trend = repmat(constant,1,gend);
  end
  start = options_.presample+1;
  np    = size(T,1);
  mf    = bayestopt_.mf;
  no_missing_data_flag = (number_of_observations==gend*nobs);
  %------------------------------------------------------------------------------
  % 3. Initial condition of the Kalman filter
  %------------------------------------------------------------------------------
  PartInfoKF=11; %'PIKF';%-1;
  if(options_.usePartInfo==1  | (options_.partial_information ==1))
        kalman_algo = PartInfoKF;
        options_.lik_init=PartInfoKF; % if lik_init not changed, kalman_algo may be changed

        NY=M_.endo_nbr;  % number of endogenous vars.
        NX=M_.exo_nbr; % no of exogenous varexo shock variables.
        FL_RANK=oo_.dr.PI_FL_RANK;
        %        NETA=oo_.dr.nfwrd+oo_.dr.nboth; % no of exp. errors  set to no of forward looking equations
        [NOBS,sample]=size(data);  % number of observations =nn and of num of samples =T  


        % Mapping the observations to the variables e.g.
        % LL(inCC, CC) =1;    % consumption
        % LL(inEE, EE) =1;	% Employment
        % ...
        % LL(inYFD, PI) =rhoYFDPI;  auxiliary, noisy relationship

        LL=zeros(NOBS,NY);
        % observation vector indices
        % mapping to endogenous variables.
        for i = 1: NOBS % size(bayestopt_.mfys,2)
            LL(i,bayestopt_.mfys(:,i))=1;
        end
        %LL=[ eye(NOBS, NX) LL zeros(NOBS,NETA)]
        % Transform measurement (observation) vector so that it is equal to
        % L1*s(t)+L2*x(t)

        L1=LL*oo_.dr.PI_TT1;
        L2=LL*oo_.dr.PI_TT2;

        SDX=Q^0.5; % =SD,not V-COV, of Exog shocks or M_.Sigma_e^0.5 num_exog x num_exog matrix
        VV=H; % V-COV of observation errors.

        G1=oo_.dr.PI_ghx;  % main transformation matrix
        impact=oo_.dr.PI_ghu;  % control impact matrix
        nmat=oo_.dr.PI_nmat;  % PI N matrix
        CC=oo_.dr.PI_CC;    % any constant(s)
        MM=impact*SDX; % R*(Q^0.5) in standard KF notation
        ss=size(G1,1); % size of state space
        shat=zeros(ss,1); % initial Sigma (P) "hat"

else

  kalman_algo = options_.kalman_algo;
  if options_.lik_init == 1		% Kalman filter
      if kalman_algo ~= 2
          kalman_algo = 1;
      end
      Pstar = lyapunov_symm(T,R*Q*R',options_.qz_criterium,options_.lyapunov_complex_threshold);
      Pinf	= [];
  elseif options_.lik_init == 2	% Old Diffuse Kalman filter
      if kalman_algo ~= 2
          kalman_algo = 1;
      end
      Pstar = options_.Harvey_scale_factor*eye(np);
      Pinf = [];
  elseif options_.lik_init == 3	% Diffuse Kalman filter
      if kalman_algo ~= 4
          kalman_algo = 3;
      end
    [Z,ST,R1,QT,Pstar,Pinf] = schur_statespace_transformation(mf,T,R,Q,options_.qz_criterium);
  end
end % part info  
  kalman_tol = options_.kalman_tol;
  riccati_tol = options_.riccati_tol;
  mf = bayestopt_.mf1;
  Y   = data-trend;
  %------------------------------------------------------------------------------
  % 4. Likelihood evaluation
  %------------------------------------------------------------------------------
  if (kalman_algo==1)% Multivariate Kalman Filter
      if no_missing_data_flag
          [LIK, lik] = kalman_filter(T,R,Q,H,Pstar,Y,start,mf,kalman_tol,riccati_tol); 
      else
          [LIK, lik] = ...
              missing_observations_kalman_filter(T,R,Q,H,Pstar,Y,start,mf,kalman_tol,riccati_tol, ...
                                                 data_index,number_of_observations,no_more_missing_observations);
      end
      if isinf(LIK)
          kalman_algo = 2;
      end
  end
  if (kalman_algo==2)% Univariate Kalman Filter
    no_correlation_flag = 1;
    if isequal(H,0)
        H = zeros(nobs,1);
    else
        if all(all(abs(H-diag(diag(H)))<1e-14))% ie, the covariance matrix is diagonal...
            H = diag(H);
        else
            no_correlation_flag = 0;
        end
    end
    if no_correlation_flag
      [LIK, lik] = univariate_kalman_filter(T,R,Q,H,Pstar,Y,start,mf,kalman_tol,riccati_tol,data_index,number_of_observations,no_more_missing_observations);
    else
      [LIK, lik] = univariate_kalman_filter_corr(T,R,Q,H,Pstar,Y,start,mf,kalman_tol,riccati_tol,data_index,number_of_observations,no_more_missing_observations);
    end
  end
  if (kalman_algo==3)% Multivariate Diffuse Kalman Filter
    if no_missing_data_flag
        [LIK, lik] = diffuse_kalman_filter(ST,R1,Q,H,Pinf,Pstar,Y,start,Z,kalman_tol, ...
                                           riccati_tol);
    else
        [LIK, lik] = missing_observations_diffuse_kalman_filter(ST,R1,Q,H,Pinf, ...
                                                          Pstar,Y,start,Z,kalman_tol,riccati_tol,...
                                                          data_index,number_of_observations,...
                                                          no_more_missing_observations);
    end
    if isinf(LIK)
        kalman_algo = 4;
    end
  end
  if (kalman_algo==4)% Univariate Diffuse Kalman Filter
    no_correlation_flag = 1;
    if isequal(H,0)
        H = zeros(nobs,1);
    else
        if all(all(abs(H-diag(diag(H)))<1e-14))% ie, the covariance matrix is diagonal...
            H = diag(H);
        else
            no_correlation_flag = 0;
        end
    end
    if no_correlation_flag
        [LIK, lik] = univariate_diffuse_kalman_filter(ST,R1,Q,H,Pinf,Pstar,Y, ...
                                                      start,Z,kalman_tol,riccati_tol,data_index,...
                                                      number_of_observations,no_more_missing_observations);
    else
        [LIK, lik] = univariate_diffuse_kalman_filter_corr(ST,R1,Q,H,Pinf,Pstar, ...
                                                          Y,start,Z,kalman_tol,riccati_tol,...
                                                          data_index,number_of_observations,...
                                                          no_more_missing_observations);
    end
  end
  if  (kalman_algo == PartInfoKF)
      VV=H; % V-COV of observation errors.
      [xhatnew,shatnew,signew,LIK,lik]=pt_info_kf_shell ...
            ( data',CC,SDX,G1,VV,MM,nmat,[],L1,L2,shat,sample,NY,NX,FL_RANK, start);
      if LIK == Inf
         LIK = bayestopt_.penalty;
      end
  end
  if imag(LIK) ~= 0
      likelihood = bayestopt_.penalty;
  else
      likelihood = LIK;
  end
  % ------------------------------------------------------------------------------
  % Adds prior if necessary
  % ------------------------------------------------------------------------------
  lnprior = priordens(xparam1,bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4);
  fval    = (likelihood-lnprior);
  options_.kalman_algo = kalman_algo;
  lik=lik(start:end,:);
  llik=[-lnprior; lik(:)];
  %llik=[-lnprior; lik(start:end)];
  
