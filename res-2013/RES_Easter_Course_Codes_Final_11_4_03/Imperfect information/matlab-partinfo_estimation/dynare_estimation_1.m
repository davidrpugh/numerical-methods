function dynare_estimation_1(var_list_,dname)
% function dynare_estimation_1(var_list_,dname)
% runs the estimation of the model
%  
% INPUTS
%   var_list_:  selected endogenous variables vector
%   dname:      alternative directory name
%  
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2003-2011 Dynare Team
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

global M_ options_ oo_ estim_params_ bayestopt_

if nargin > 1
  [data,rawdata,xparam1] = dynare_estimation_init(var_list_,dname);
else
  [data,rawdata,xparam1] = dynare_estimation_init(var_list_);
end

%% Set various options.
options_ = set_default_option(options_,'mh_nblck',2); 
options_ = set_default_option(options_,'nodiagnostic',0);

%% Set number of observations
gend = options_.nobs;
%% Set the number of observed variables.
n_varobs = size(options_.varobs,1);
%% Get the number of parameters to be estimated. 
nvx = estim_params_.nvx;  % Variance of the structural innovations (number of parameters).
nvn = estim_params_.nvn;  % Variance of the measurement innovations (number of parameters).
ncx = estim_params_.ncx;  % Covariance of the structural innovations (number of parameters).
ncn = estim_params_.ncn;  % Covariance of the measurement innovations (number of parameters).
np  = estim_params_.np ;  % Number of deep parameters.
nx  = nvx+nvn+ncx+ncn+np; % Total number of parameters to be estimated.
%% Set the names of the priors.
pnames = ['     ';'beta ';'gamm ';'norm ';'invg ';'unif ';'invg2'];
%% Set parameters bounds
lb = bayestopt_.lb;
ub = bayestopt_.ub;

dr = oo_.dr;

%% load mode file is necessary
if ~isempty(options_.mode_file) && ~options_.mh_posterior_mode_estimation
    load(options_.mode_file);
end

%% Compute the steady state: 
if ~isempty(estim_params_)
    set_parameters(xparam1);
end
if options_.steadystate_flag% if the *_steadystate.m file is provided.
    [ys,tchek] = feval([M_.fname '_steadystate'],...
                       [zeros(M_.exo_nbr,1);...
                        oo_.exo_det_steady_state]);
    if size(ys,1) < M_.endo_nbr 
        if length(M_.aux_vars) > 0
            ys = add_auxiliary_variables_to_steadystate(ys,M_.aux_vars,...
                                                        M_.fname,...
                                                        zeros(M_.exo_nbr,1),...
                                                        oo_.exo_det_steady_state,...
                                                        M_.params,...
                                                        options_.bytecode);
        else
            error([M_.fname '_steadystate.m doesn''t match the model']);
        end
    end
    oo_.steady_state = ys;
else% if the steady state file is not provided.
    [dd,info] = resol(oo_.steady_state,0);
    oo_.steady_state = dd.ys; clear('dd');
end
if all(abs(oo_.steady_state(bayestopt_.mfys))<1e-9)
    options_.noconstant = 1;
else
    options_.noconstant = 0;
end


%% compute sample moments if needed (bvar-dsge)
if options_.dsge_var
    if options_.missing_data
        error('I cannot estimate a DSGE-VAR model with missing observations!')
    end
    if options_.noconstant
        evalin('base',...
               ['[mYY,mXY,mYX,mXX,Ydata,Xdata] = ' ...
                'var_sample_moments(options_.first_obs,' ... 
                'options_.first_obs+options_.nobs-1,options_.dsge_varlag,-1,' ...
                'options_.datafile, options_.varobs,options_.xls_sheet,options_.xls_range);'])
    else% The steady state is non zero ==> a constant in the VAR is needed!
        evalin('base',['[mYY,mXY,mYX,mXX,Ydata,Xdata] = ' ...
                       'var_sample_moments(options_.first_obs,' ...
                       'options_.first_obs+options_.nobs-1,options_.dsge_varlag,0,' ...
                       'options_.datafile, options_.varobs,options_.xls_sheet,options_.xls_range);'])
    end
end

[data_index,number_of_observations,no_more_missing_observations] = describe_missing_data(data,gend,n_varobs);
missing_value = ~(number_of_observations == gend*n_varobs);

initial_estimation_checks(xparam1,gend,data,data_index,number_of_observations,no_more_missing_observations);

if isequal(options_.mode_compute,0) && isempty(options_.mode_file) && options_.mh_posterior_mode_estimation==0
    if options_.smoother == 1
        [atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,T,R,P,PK,decomp] = DsgeSmoother(xparam1,gend,data,data_index,missing_value);
        oo_.Smoother.SteadyState = ys;
        oo_.Smoother.TrendCoeffs = trend_coeff;
        if options_.filter_covariance
            oo_.Smoother.variance = P;
        end
        i_endo = bayestopt_.smoother_saved_var_list;
        if options_.nk ~= 0
            oo_.FilteredVariablesKStepAhead = ...
                aK(options_.filter_step_ahead,i_endo,:);
            if ~isempty(PK)
                oo_.FilteredVariablesKStepAheadVariances = ...
                    PK(options_.filter_step_ahead,i_endo,i_endo,:);
            end
            if ~isempty(decomp)
                oo_.FilteredVariablesShockDecomposition = ...
                    decomp(options_.filter_step_ahead,i_endo,:,:);
            end
        end
        for i=bayestopt_.smoother_saved_var_list'
            i1 = dr.order_var(bayestopt_.smoother_var_list(i));
            eval(['oo_.SmoothedVariables.' deblank(M_.endo_names(i1,:)) ' = atT(i,:)'';']);
            eval(['oo_.FilteredVariables.' deblank(M_.endo_names(i1,:)) ' = squeeze(aK(1,i,:));']);
            eval(['oo_.UpdatedVariables.' deblank(M_.endo_names(i1,:)) ' = updated_variables(i,:)'';']);
        end
        for i=1:M_.exo_nbr
            eval(['oo_.SmoothedShocks.' deblank(M_.exo_names(i,:)) ' = innov(i,:)'';']);
        end
    end

    return;
end

if isequal(options_.mode_compute,6)
    % Erase previously computed optimal mh scale parameter.
    delete([M_.fname '_optimal_mh_scale_parameter.mat'])
end


%% Estimation of the posterior mode or likelihood mode
if ~isequal(options_.mode_compute,0) && ~options_.mh_posterior_mode_estimation
    if ~options_.dsge_var
        fh=str2func('DsgeLikelihood');
    else
        fh=str2func('DsgeVarLikelihood');
    end
    switch options_.mode_compute
      case 1
        optim_options = optimset('display','iter','LargeScale','off', ...
                                 'MaxFunEvals',100000,'TolFun',1e-8,'TolX',1e-6);
        if isfield(options_,'optim_opt')
            eval(['optim_options = optimset(optim_options,' options_.optim_opt ');']);
        end
        if ~options_.dsge_var
            [xparam1,fval,exitflag,output,lamdba,grad,hessian_fmincon] = ...
                fmincon(fh,xparam1,[],[],[],[],lb,ub,[],optim_options,gend,data,data_index,number_of_observations,no_more_missing_observations);
        else
            [xparam1,fval,exitflag,output,lamdba,grad,hessian_fmincon] = ...
                fmincon(fh,xparam1,[],[],[],[],lb,ub,[],optim_options,gend);
        end
      case 2
        error('ESTIMATION: mode_compute=2 option (Lester Ingber''s Adaptive Simulated Annealing) is no longer available')
      case 3
        optim_options = optimset('display','iter','MaxFunEvals',100000,'TolFun',1e-8,'TolX',1e-6);
        if isfield(options_,'optim_opt')
            eval(['optim_options = optimset(optim_options,' options_.optim_opt ');']);
        end
        if ~options_.dsge_var
            [xparam1,fval,exitflag] = fminunc(fh,xparam1,optim_options,gend,data,data_index,number_of_observations,no_more_missing_observations);
        else
            [xparam1,fval,exitflag] = fminunc(fh,xparam1,optim_options,gend);
        end
      case 4
        H0 = 1e-4*eye(nx);
        crit = 1e-7;
        nit = 1000;
        verbose = 2;
        if ~options_.dsge_var
            [fval,xparam1,grad,hessian_csminwel,itct,fcount,retcodehat] = ...
                csminwel1('DsgeLikelihood',xparam1,H0,[],crit,nit,options_.gradient_method,options_.gradient_epsilon,gend,data,data_index,number_of_observations,no_more_missing_observations);
            disp(sprintf('Objective function at mode: %f',fval))
            disp(sprintf('Objective function at mode: %f',DsgeLikelihood(xparam1,gend,data,data_index,number_of_observations,no_more_missing_observations)))
        else
            [fval,xparam1,grad,hessian_csminwel,itct,fcount,retcodehat] = ...
                csminwel1('DsgeVarLikelihood',xparam1,H0,[],crit,nit,options_.gradient_method,options_.gradient_epsilon,gend);
            disp(sprintf('Objective function at mode: %f',fval))
            disp(sprintf('Objective function at mode: %f',DsgeVarLikelihood(xparam1,gend)))
        end
      case 5
        if isfield(options_,'hess')
            flag = options_.hess;
        else
            flag = 1;
        end
        if ~exist('igg','var'),  % by M. Ratto
            hh=[];
            gg=[];
            igg=[];
        end   % by M. Ratto
        if isfield(options_,'ftol')
            crit = options_.ftol;
        else
            crit = 1.e-5;
        end
        if isfield(options_,'nit')
            nit = options_.nit;
        else
            nit=1000;
        end
        if ~options_.dsge_var
            [xparam1,hh,gg,fval,invhess] = newrat('DsgeLikelihood',xparam1,hh,gg,igg,crit,nit,flag,gend,data,data_index,number_of_observations,no_more_missing_observations);
        else
            [xparam1,hh,gg,fval,invhess] = newrat('DsgeVarLikelihood',xparam1,hh,gg,igg,crit,nit,flag,gend);
        end
        parameter_names = bayestopt_.name;
        save([M_.fname '_mode.mat'],'xparam1','hh','gg','fval','invhess','parameter_names');
      case 6
        if ~options_.dsge_var
            fval = DsgeLikelihood(xparam1,gend,data,data_index,number_of_observations,no_more_missing_observations);
        else
            fval = DsgeVarLikelihood(xparam1,gend);
        end
        OldMode = fval;
        if ~exist('MeanPar','var')
            MeanPar = xparam1;
        end
        if exist('hh','var')
            CovJump = inv(hh);
        else% The covariance matrix is initialized with the prior
            % covariance (a diagonal matrix) %%Except for infinite variances ;-)
            varinit = 'prior';
            if strcmpi(varinit,'prior')  
                stdev = bayestopt_.p2;
                indx = find(isinf(stdev));
                stdev(indx) = ones(length(indx),1)*sqrt(10);
                vars = stdev.^2;
                CovJump = diag(vars);
            elseif strcmpi(varinit,'eye')
                vars = ones(length(bayestopt_.p2),1)*0.1;  
                CovJump = diag(vars);          
            else
                disp('gmhmaxlik :: Error!')
                return
            end
        end
        OldPostVar = CovJump;
        Scale = options_.mh_jscale;
        for i=1:options_.Opt6Iter
            if i == 1
                if options_.Opt6Iter > 1
                    flag = '';
                else
                    flag = 'LastCall';
                end
                if ~options_.dsge_var
                    [xparam1,PostVar,Scale,PostMean] = ...
                        gmhmaxlik('DsgeLikelihood',xparam1,[lb ub],options_.Opt6Numb,Scale,flag,MeanPar,CovJump,gend,data,...
                                  data_index,number_of_observations,no_more_missing_observations);
                    fval = DsgeLikelihood(xparam1,gend,data,data_index,number_of_observations,no_more_missing_observations);
                else
                    [xparam1,PostVar,Scale,PostMean] = ...
                        gmhmaxlik('DsgeVarLikelihood',xparam1,[lb ub],options_.Opt6Numb,Scale,flag,MeanPar,CovJump,gend);
                    fval = DsgeVarLikelihood(xparam1,gend);
                end
                options_.mh_jscale = Scale;
                mouvement = max(max(abs(PostVar-OldPostVar)));
                disp(' ')
                disp('========================================================== ')
                disp(['   Change in the covariance matrix = ' num2str(mouvement) '.'])
                disp(['   Mode improvement = ' num2str(abs(OldMode-fval))])
                disp(['   New value of jscale = ' num2str(Scale)])
                disp('========================================================== ')
                OldMode = fval;
            else
                OldPostVar = PostVar;
                if i<options_.Opt6Iter
                    flag = '';
                else
                    flag = 'LastCall';
                end
                if ~options_.dsge_var
                    [xparam1,PostVar,Scale,PostMean] = ...
                        gmhmaxlik('DsgeLikelihood',xparam1,[lb ub],...
                                  options_.Opt6Numb,Scale,flag,PostMean,PostVar,gend,data,data_index,number_of_observations,no_more_missing_observations);
                    fval = DsgeLikelihood(xparam1,gend,data,data_index,number_of_observations,no_more_missing_observations);
                else
                    [xparam1,PostVar,Scale,PostMean] = ...
                        gmhmaxlik('DsgeVarLikelihood',xparam1,[lb ub],...
                                  options_.Opt6Numb,Scale,flag,PostMean,PostVar,gend);
                    fval = DsgeVarLikelihood(xparam1,gend);          
                end
                options_.mh_jscale = Scale;
                mouvement = max(max(abs(PostVar-OldPostVar)));
                fval = DsgeLikelihood(xparam1,gend,data,data_index,number_of_observations,no_more_missing_observations);
                disp(['Change in the covariance matrix = ' num2str(mouvement) '.'])
                disp(['Mode improvement = ' num2str(abs(OldMode-fval))])
                OldMode = fval;
            end
            hh = inv(PostVar);
            save([M_.fname '_mode.mat'],'xparam1','hh');
            save([M_.fname '_optimal_mh_scale_parameter.mat'],'Scale');
            bayestopt_.jscale = ones(length(xparam1),1)*Scale;%??!
        end    
      case 7
        optim_options = optimset('display','iter','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6);
        if isfield(options_,'optim_opt')
            eval(['optim_options = optimset(optim_options,' options_.optim_opt ');']);
        end
        if ~options_.dsge_var
            [xparam1,fval,exitflag] = fminsearch(fh,xparam1,optim_options,gend,data,data_index,number_of_observations,no_more_missing_observations);
        else
            [xparam1,fval,exitflag] = fminsearch(fh,xparam1,optim_options,gend);
        end
        
      case 101
        myoptions=soptions;
        myoptions(2)=1e-6; %accuracy of argument
        myoptions(3)=1e-6; %accuracy of function (see Solvopt p.29)
        myoptions(5)= 1.0;
        
        [xparam1,fval]=solvopt(xparam1,fh,[],myoptions,gend,data);
      case 102
        %simulating annealing
        %        LB=zeros(size(xparam1))-20;
        % UB=zeros(size(xparam1))+20;
        LB = xparam1 - 1;
        UB = xparam1 + 1;
        neps=10;
        %  Set input parameters. 
        maxy=0;
        epsilon=1.0e-9;
        rt_=.10;
        t=15.0;
        ns=10;
        nt=10;
        maxevl=100000000;
        idisp =1;
        npar=length(xparam1);
        
        disp(['size of param',num2str(length(xparam1))])    
        c=.1*ones(npar,1);
        %*  Set input values of the input/output parameters.*
        
        vm=1*ones(npar,1);
        disp(['number of parameters= ' num2str(npar) 'max= '  num2str(maxy) 't=  ' num2str(t)]);
        disp(['rt_=  '  num2str(rt_) 'eps=  '  num2str(epsilon) 'ns=  '  num2str(ns)]);
        disp(['nt=  '  num2str(nt) 'neps= '   num2str(neps) 'maxevl=  '  num2str(maxevl)]);
        %      disp(['iprint=   '   num2str(iprint) 'seed=   '   num2str(seed)]);
        disp '  ';
        disp '  ';
        disp(['starting values(x)  ' num2str(xparam1')]);
        disp(['initial step length(vm)  '  num2str(vm')]);
        disp(['lower bound(lb)', 'initial conditions', 'upper bound(ub)' ]);
        disp([LB xparam1 UB]);
        disp(['c vector   ' num2str(c')]);
        
        %  keyboard 
        if ~options_.dsge_var
            [xparam1, fval, nacc, nfcnev, nobds, ier, t, vm] = sa(fh,xparam1,maxy,rt_,epsilon,ns,nt ...
                                                              ,neps,maxevl,LB,UB,c,idisp ,t,vm,gend,data,data_index,number_of_observations,no_more_missing_observations);
        else
            [xparam1, fval, nacc, nfcnev, nobds, ier, t, vm] = sa(fh,xparam1,maxy,rt_,epsilon,ns,nt ...
                                                              ,neps,maxevl,LB,UB,c,idisp ,t,vm,gend);
        end
      otherwise
        if ischar(options_.mode_compute)
            if options_.dsge_var
                [xparam1, fval, retcode ] = feval(options_.mode_compute,fh,xparam1,gend,data);
            else
                [xparam1, fval, retcode ] = feval(options_.mode_compute, ...
                                                  fh,xparam1,gend,data,data_index,number_of_observations,no_more_missing_observations);
            end
        else
            error(['ESTIMATION: mode_compute=' int2str(options_.mode_compute) ...
                   ' option is unknown!'])
        end
    end
    if ~isequal(options_.mode_compute,6)
        if options_.cova_compute == 1
            if ~options_.dsge_var
                hh = reshape(hessian('DsgeLikelihood',xparam1, ...
                                     options_.gstep,gend,data,data_index,number_of_observations,...
                                     no_more_missing_observations),nx,nx);
            else
                hh = reshape(hessian('DsgeVarLikelihood',xparam1,options_.gstep,gend),nx,nx);
            end
        end
    end
    parameter_names = bayestopt_.name;
    if options_.cova_compute
        save([M_.fname '_mode.mat'],'xparam1','hh','parameter_names');
    else
        save([M_.fname '_mode.mat'],'xparam1','parameter_names');
    end
end

if options_.cova_compute == 0
    hh = NaN(length(xparam1),length(xparam1));
end

if ~options_.mh_posterior_mode_estimation && options_.cova_compute
    try
        chol(hh);
    catch
        disp(' ')
        disp('POSTERIOR KERNEL OPTIMIZATION PROBLEM!')
        disp(' (minus) the hessian matrix at the "mode" is not positive definite!')
        disp('=> posterior variance of the estimated parameters are not positive.')
        disp('You should  try  to change the initial values of the parameters using')
        disp('the estimated_params_init block, or use another optimization routine.')
        warning('The results below are most likely wrong!');
    end
end

if options_.mode_check == 1 && ~options_.mh_posterior_mode_estimation && options_.cova_compute
    if options_.cova_compute
        mode_check(xparam1,0,hh,gend,data,lb,ub,data_index,number_of_observations,no_more_missing_observations);
    else
        mode_check(xparam1,0,[],gend,data,lb,ub,data_index,number_of_observations,no_more_missing_observations);
    end
end

if ~options_.mh_posterior_mode_estimation
    if options_.cova_compute
        invhess = inv(hh);
        stdh = sqrt(diag(invhess));
    end
else
    variances = bayestopt_.p2.*bayestopt_.p2;
    idInf = isinf(variances);
    variances(idInf) = 1;
    invhess = options_.mh_posterior_mode_estimation*diag(variances);
    xparam1 = bayestopt_.p5;
    idNaN = isnan(xparam1);
    xparam1(idNaN) = bayestopt_.p1(idNaN);
    xparam1 = transpose(xparam1);
end

if ~options_.cova_compute
    stdh = NaN(length(xparam1),1);
end

if any(bayestopt_.pshape > 0) && ~options_.mh_posterior_mode_estimation
    disp(' ')
    disp('RESULTS FROM POSTERIOR MAXIMIZATION')
    tstath = zeros(nx,1);
    for i = 1:nx
        tstath(i) = abs(xparam1(i))/stdh(i);
    end
    
    header_width = row_header_width(M_,estim_params_,bayestopt_);
    
    tit1 = sprintf('%-*s %7s %8s %7s %6s %4s %6s\n',header_width-2,' ','prior mean', ...
                   'mode','s.d.','t-stat','prior','pstdev');
    if np
        ip = nvx+nvn+ncx+ncn+1;
        disp('parameters')
        disp(tit1)
        for i=1:np
            name = bayestopt_.name{ip};
            disp(sprintf('%-*s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', ...
                         header_width,name, ...
                         bayestopt_.p1(ip),xparam1(ip),stdh(ip),tstath(ip), ...
                         pnames(bayestopt_.pshape(ip)+1,:), ...
                         bayestopt_.p2(ip)));
            eval(['oo_.posterior_mode.parameters.' name ' = xparam1(ip);']);
            eval(['oo_.posterior_std.parameters.' name ' = stdh(ip);']); 
            ip = ip+1;
        end
    end
    if nvx
        ip = 1;
        disp('standard deviation of shocks')
        disp(tit1)
        for i=1:nvx
            k = estim_params_.var_exo(i,1);
            name = deblank(M_.exo_names(k,:));
            disp(sprintf('%-*s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', ...
                         header_width,name,bayestopt_.p1(ip),xparam1(ip), ...
                         stdh(ip),tstath(ip),pnames(bayestopt_.pshape(ip)+1,:), ...
                         bayestopt_.p2(ip))); 
            M_.Sigma_e(k,k) = xparam1(ip)*xparam1(ip);
            eval(['oo_.posterior_mode.shocks_std.' name ' = xparam1(ip);']);
            eval(['oo_.posterior_std.shocks_std.' name ' = stdh(ip);']); 
            ip = ip+1;
        end
    end
    if nvn
        disp('standard deviation of measurement errors')
        disp(tit1)
        ip = nvx+1;
        for i=1:nvn
            name = deblank(options_.varobs(estim_params_.var_endo(i,1),:));
            disp(sprintf('%-*s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', ...
                         header_width,name,bayestopt_.p1(ip), ...
                         xparam1(ip),stdh(ip),tstath(ip), ...
                         pnames(bayestopt_.pshape(ip)+1,:), ...
                         bayestopt_.p2(ip)));
            eval(['oo_.posterior_mode.measurement_errors_std.' name ' = xparam1(ip);']);
            eval(['oo_.posterior_std.measurement_errors_std.' name ' = stdh(ip);']); 
            ip = ip+1;
        end
    end
    if ncx
        disp('correlation of shocks')
        disp(tit1)
        ip = nvx+nvn+1;
        for i=1:ncx
            k1 = estim_params_.corrx(i,1);
            k2 = estim_params_.corrx(i,2);
            name = [deblank(M_.exo_names(k1,:)) ',' deblank(M_.exo_names(k2,:))];
            NAME = [deblank(M_.exo_names(k1,:)) '_' deblank(M_.exo_names(k2,:))];
            disp(sprintf('%-*s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', name, ...
                         header_width,bayestopt_.p1(ip),xparam1(ip),stdh(ip),tstath(ip),  ...
                         pnames(bayestopt_.pshape(ip)+1,:), bayestopt_.p2(ip)));
            M_.Sigma_e(k1,k2) = xparam1(ip)*sqrt(M_.Sigma_e(k1,k1)*M_.Sigma_e(k2,k2));
            M_.Sigma_e(k2,k1) = M_.Sigma_e(k1,k2);
            eval(['oo_.posterior_mode.shocks_corr.' NAME ' = xparam1(ip);']);
            eval(['oo_.posterior_std.shocks_corr.' NAME ' = stdh(ip);']); 
            ip = ip+1;
        end
    end
    if ncn
        disp('correlation of measurement errors')
        disp(tit1)
        ip = nvx+nvn+ncx+1;
        for i=1:ncn
            k1 = estim_params_.corrn(i,1);
            k2 = estim_params_.corrn(i,2);
            name = [deblank(M_.endo_names(k1,:)) ',' deblank(M_.endo_names(k2,:))];
            NAME = [deblank(M_.endo_names(k1,:)) '_' deblank(M_.endo_names(k2,:))];
            disp(sprintf('%-*s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', name, ...
                         header_width,bayestopt_.p1(ip),xparam1(ip),stdh(ip),tstath(ip), ...
                         pnames(bayestopt_.pshape(ip)+1,:), bayestopt_.p2(ip)));
            eval(['oo_.posterior_mode.measurement_errors_corr.' NAME ' = xparam1(ip);']);
            eval(['oo_.posterior_std.measurement_errors_corr.' NAME ' = stdh(ip);']); 
            ip = ip+1;
        end
    end
    %% Laplace approximation to the marginal log density:
    if options_.cova_compute
        estim_params_nbr = size(xparam1,1);
        scale_factor = -sum(log10(diag(invhess)));
        log_det_invhess = -estim_params_nbr*log(scale_factor)+log(det(scale_factor*invhess));
        if ~options_.dsge_var
            md_Laplace = .5*estim_params_nbr*log(2*pi) + .5*log_det_invhess ...
                - DsgeLikelihood(xparam1,gend,data,data_index,number_of_observations,no_more_missing_observations);
        else
            md_Laplace = .5*estim_params_nbr*log(2*pi) + .5*log_det_invhess ...
                - DsgeVarLikelihood(xparam1,gend);
        end
        oo_.MarginalDensity.LaplaceApproximation = md_Laplace;    
        disp(' ')
        disp(sprintf('Log data density [Laplace approximation] is %f.',md_Laplace))
        disp(' ')
    end
elseif ~any(bayestopt_.pshape > 0) && options_.mh_posterior_mode_estimation
    disp(' ')
    disp('RESULTS FROM MAXIMUM LIKELIHOOD')
    tstath = zeros(nx,1);
    for i = 1:nx
        tstath(i) = abs(xparam1(i))/stdh(i);
    end
    header_width = row_header_width(M_,estim_params_,bayestopt_);
    tit1 = sprintf('%-*s %10s %7s %6s\n',header_width-2,' ','Estimate','s.d.','t-stat');
    if np
        ip = nvx+nvn+ncx+ncn+1;
        disp('parameters')
        disp(tit1)
        for i=1:np
            name = bayestopt_.name{ip};
            disp(sprintf('%-*s %8.4f %7.4f %7.4f', ...
                         header_width,name,xparam1(ip),stdh(ip),tstath(ip)));
            eval(['oo_.mle_mode.parameters.' name ' = xparam1(ip);']);
            eval(['oo_.mle_std.parameters.' name ' = stdh(ip);']); 
            ip = ip+1;
        end
    end 
    if nvx
        ip = 1;
        disp('standard deviation of shocks')
        disp(tit1)
        for i=1:nvx
            k = estim_params_.var_exo(i,1);
            name = deblank(M_.exo_names(k,:));
            disp(sprintf('%-*s %8.4f %7.4f %7.4f',header_width,name,xparam1(ip),stdh(ip),tstath(ip)));
            M_.Sigma_e(k,k) = xparam1(ip)*xparam1(ip);
            eval(['oo_.mle_mode.shocks_std.' name ' = xparam1(ip);']);
            eval(['oo_.mle_std.shocks_std.' name ' = stdh(ip);']); 
            ip = ip+1;
        end
    end
    if nvn
        disp('standard deviation of measurement errors')
        disp(tit1)
        ip = nvx+1;
        for i=1:nvn
            name = deblank(options_.varobs(estim_params_.var_endo(i,1),:));
            disp(sprintf('%-*s %8.4f %7.4f %7.4f',header_width,name,xparam1(ip),stdh(ip),tstath(ip)))
            eval(['oo_.mle_mode.measurement_errors_std.' name ' = xparam1(ip);']);
            eval(['oo_.mle_std.measurement_errors_std.' name ' = stdh(ip);']);      
            ip = ip+1;
        end
    end
    if ncx
        disp('correlation of shocks')
        disp(tit1)
        ip = nvx+nvn+1;
        for i=1:ncx
            k1 = estim_params_.corrx(i,1);
            k2 = estim_params_.corrx(i,2);
            name = [deblank(M_.exo_names(k1,:)) ',' deblank(M_.exo_names(k2,:))];
            NAME = [deblank(M_.exo_names(k1,:)) '_' deblank(M_.exo_names(k2,:))];
            disp(sprintf('%-*s %8.4f %7.4f %7.4f', header_width,name,xparam1(ip),stdh(ip),tstath(ip)));
            M_.Sigma_e(k1,k2) = xparam1(ip)*sqrt(M_.Sigma_e(k1,k1)*M_.Sigma_e(k2,k2));
            M_.Sigma_e(k2,k1) = M_.Sigma_e(k1,k2);
            eval(['oo_.mle_mode.shocks_corr.' NAME ' = xparam1(ip);']);
            eval(['oo_.mle_std.shocks_corr.' NAME ' = stdh(ip);']);      
            ip = ip+1;
        end
    end
    if ncn
        disp('correlation of measurement errors')
        disp(tit1)
        ip = nvx+nvn+ncx+1;
        for i=1:ncn
            k1 = estim_params_.corrn(i,1);
            k2 = estim_params_.corrn(i,2);
            name = [deblank(M_.endo_names(k1,:)) ',' deblank(M_.endo_names(k2,:))];
            NAME = [deblank(M_.endo_names(k1,:)) '_' deblank(M_.endo_names(k2,:))];
            disp(sprintf('%-*s %8.4f %7.4f %7.4f',header_width,name,xparam1(ip),stdh(ip),tstath(ip)));
            eval(['oo_.mle_mode.measurement_error_corr.' NAME ' = xparam1(ip);']);
            eval(['oo_.mle_std.measurement_error_corr.' NAME ' = stdh(ip);']);
            ip = ip+1;
        end
    end
end


OutputDirectoryName = CheckPath('Output');

if any(bayestopt_.pshape > 0) && options_.TeX %% Bayesian estimation (posterior mode) Latex output
    if np
        filename = [OutputDirectoryName '/' M_.fname '_Posterior_Mode_1.TeX'];
        fidTeX = fopen(filename,'w');
        fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
        fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (parameters)\n');
        fprintf(fidTeX,['%% ' datestr(now,0)]);
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,'{\\tiny \n');
        fprintf(fidTeX,'\\begin{table}\n');
        fprintf(fidTeX,'\\centering\n');
        fprintf(fidTeX,'\\begin{tabular}{l|lcccc} \n');
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Posterior mode & s.d. \\\\ \n');
        fprintf(fidTeX,'\\hline \\\\ \n');
        ip = nvx+nvn+ncx+ncn+1;
        for i=1:np
            fprintf(fidTeX,'$%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n',...
                    M_.param_names_tex(estim_params_.param_vals(i,1),:),...
                    deblank(pnames(bayestopt_.pshape(ip)+1,:)),...
                    bayestopt_.p1(ip),...
                    bayestopt_.p2(ip),...
                    xparam1(ip),...
                    stdh(ip));
            ip = ip + 1;    
        end
        fprintf(fidTeX,'\\hline\\hline \n');
        fprintf(fidTeX,'\\end{tabular}\n ');    
        fprintf(fidTeX,'\\caption{Results from posterior parameters (parameters)}\n ');
        fprintf(fidTeX,'\\label{Table:Posterior:1}\n');
        fprintf(fidTeX,'\\end{table}\n');
        fprintf(fidTeX,'} \n');
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
    end
    if nvx
        TeXfile = [OutputDirectoryName '/' M_.fname '_Posterior_Mode_2.TeX'];
        fidTeX = fopen(TeXfile,'w');
        fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
        fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (standard deviation of structural shocks)\n');
        fprintf(fidTeX,['%% ' datestr(now,0)]);
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,'{\\tiny \n');
        fprintf(fidTeX,'\\begin{table}\n');
        fprintf(fidTeX,'\\centering\n');
        fprintf(fidTeX,'\\begin{tabular}{l|lcccc} \n');
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Posterior mode & s.d. \\\\ \n');
        fprintf(fidTeX,'\\hline \\\\ \n');
        ip = 1;
        for i=1:nvx
            k = estim_params_.var_exo(i,1);
            fprintf(fidTeX,[ '$%s$ & %4s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n'],...
                    deblank(M_.exo_names_tex(k,:)),...
                    deblank(pnames(bayestopt_.pshape(ip)+1,:)),...
                    bayestopt_.p1(ip),...
                    bayestopt_.p2(ip),...
                    xparam1(ip), ...
                    stdh(ip)); 
            ip = ip+1;
        end
        fprintf(fidTeX,'\\hline\\hline \n');
        fprintf(fidTeX,'\\end{tabular}\n ');    
        fprintf(fidTeX,'\\caption{Results from posterior parameters (standard deviation of structural shocks)}\n ');
        fprintf(fidTeX,'\\label{Table:Posterior:2}\n');
        fprintf(fidTeX,'\\end{table}\n');
        fprintf(fidTeX,'} \n');
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
    end
    if nvn
        TeXfile = [OutputDirectoryName '/' M_.fname '_Posterior_Mode_3.TeX'];
        fidTeX  = fopen(TeXfile,'w');
        fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
        fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (standard deviation of measurement errors)\n');
        fprintf(fidTeX,['%% ' datestr(now,0)]);
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,'\\begin{table}\n');
        fprintf(fidTeX,'\\centering\n');
        fprintf(fidTeX,'\\begin{tabular}{l|lcccc} \n');
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. &  Posterior mode & s.d. \\\\ \n')
        fprintf(fidTeX,'\\hline \\\\ \n');
        ip = nvx+1;
        for i=1:nvn
            idx = strmatch(options_.varobs(estim_params_.var_endo(i,1),:),M_.endo_names);
            fprintf(fidTeX,'$%s$ & %4s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n',...
                    deblank(M_.endo_names_tex(idx,:)), ...
                    deblank(pnames(bayestopt_.pshape(ip)+1,:)), ...        
                    bayestopt_.p1(ip), ...
                    bayestopt_.p2(ip),...        
                    xparam1(ip),...
                    stdh(ip)); 
            ip = ip+1;
        end
        fprintf(fidTeX,'\\hline\\hline \n');
        fprintf(fidTeX,'\\end{tabular}\n ');    
        fprintf(fidTeX,'\\caption{Results from posterior parameters (standard deviation of measurement errors)}\n ');
        fprintf(fidTeX,'\\label{Table:Posterior:3}\n');
        fprintf(fidTeX,'\\end{table}\n');
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
    end
    if ncx
        TeXfile = [OutputDirectoryName '/' M_.fname '_Posterior_Mode_4.TeX'];
        fidTeX = fopen(TeXfile,'w');
        fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
        fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (correlation of structural shocks)\n');
        fprintf(fidTeX,['%% ' datestr(now,0)]);
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,'\\begin{table}\n');
        fprintf(fidTeX,'\\centering\n');
        fprintf(fidTeX,'\\begin{tabular}{l|lcccc} \n');
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. &  Posterior mode & s.d. \\\\ \n')
        fprintf(fidTeX,'\\hline \\\\ \n');
        ip = nvx+nvn+1;
        for i=1:ncx
            k1 = estim_params_.corrx(i,1);
            k2 = estim_params_.corrx(i,2);
            fprintf(fidTeX,[ '$%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n'],...
                    [deblank(M_.exo_names_tex(k1,:)) ',' deblank(M_.exo_names_tex(k2,:))], ...
                    deblank(pnames(bayestopt_.pshape(ip)+1,:)), ...
                    bayestopt_.p1(ip), ...
                    bayestopt_.p2(ip), ...
                    xparam1(ip), ...
                    stdh(ip));
            ip = ip+1;
        end
        fprintf(fidTeX,'\\hline\\hline \n');
        fprintf(fidTeX,'\\end{tabular}\n ');    
        fprintf(fidTeX,'\\caption{Results from posterior parameters (correlation of structural shocks)}\n ');
        fprintf(fidTeX,'\\label{Table:Posterior:4}\n');
        fprintf(fidTeX,'\\end{table}\n');
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
    end
    if ncn
        TeXfile = [OutputDirectoryName '/' M_.fname '_Posterior_Mode_5.TeX'];
        fidTeX = fopen(TeXfile,'w');
        fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
        fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (correlation of measurement errors)\n');
        fprintf(fidTeX,['%% ' datestr(now,0)]);
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,'\\begin{table}\n');
        fprintf(fidTeX,'\\centering\n');
        fprintf(fidTeX,'\\begin{tabular}{l|lcccc} \n');
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. &  Posterior mode & s.d. \\\\ \n')
        fprintf(fidTeX,'\\hline \\\\ \n');
        ip = nvx+nvn+ncx+1;
        for i=1:ncn
            k1 = estim_params_.corrn(i,1);
            k2 = estim_params_.corrn(i,2);
            fprintf(fidTeX,'$%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n',...
                    [deblank(M_.endo_names_tex(k1,:)) ',' deblank(M_.endo_names_tex(k2,:))], ...
                    pnames(bayestopt_.pshape(ip)+1,:), ...
                    bayestopt_.p1(ip), ...
                    bayestopt_.p2(ip), ...
                    xparam1(ip), ...
                    stdh(ip));
            ip = ip+1;
        end
        fprintf(fidTeX,'\\hline\\hline \n');
        fprintf(fidTeX,'\\end{tabular}\n ');    
        fprintf(fidTeX,'\\caption{Results from posterior parameters (correlation of measurement errors)}\n ');
        fprintf(fidTeX,'\\label{Table:Posterior:5}\n');
        fprintf(fidTeX,'\\end{table}\n');
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
    end
end

if np > 0
    pindx = estim_params_.param_vals(:,1);
    save([M_.fname '_params.mat'],'pindx');
end

if (any(bayestopt_.pshape  >0 ) && options_.mh_replic) || ...
        (any(bayestopt_.pshape >0 ) && options_.load_mh_file)  %% not ML estimation
    bounds = prior_bounds(bayestopt_);
    bounds(:,1)=max(bounds(:,1),lb);
    bounds(:,2)=min(bounds(:,2),ub);
    bayestopt_.lb = bounds(:,1);
    bayestopt_.ub = bounds(:,2);
    if any(xparam1 < bounds(:,1)) || any(xparam1 > bounds(:,2))
        find(xparam1 < bounds(:,1))
        find(xparam1 > bounds(:,2))
        error('Mode values are outside prior bounds. Reduce prior_trunc.')
    end
    % runs MCMC
    if options_.mh_replic
        if options_.load_mh_file && options_.use_mh_covariance_matrix
            invhess = compute_mh_covariance_matrix;
        end
        if options_.cova_compute
            if options_.dsge_var
                feval(options_.posterior_sampling_method,'DsgeVarLikelihood',options_.proposal_distribution,xparam1,invhess,bounds,gend);
            else
                feval(options_.posterior_sampling_method,'DsgeLikelihood',options_.proposal_distribution,xparam1,invhess,bounds,gend,data,...
                      data_index,number_of_observations,no_more_missing_observations);
            end
        else
            error('I Cannot start the MCMC because the hessian of the posterior kernel at the mode was not computed.')
        end
    end
    if options_.mh_posterior_mode_estimation
        CutSample(M_, options_, estim_params_);
        return
    else
        if ~options_.nodiagnostic && options_.mh_replic > 1000 && options_.mh_nblck > 1
            McMCDiagnostics(options_, estim_params_, M_);
        end
        %% Here i discard first half of the draws:
        CutSample(M_, options_, estim_params_);
        %% Estimation of the marginal density from the Mh draws:
        if options_.mh_replic
            [marginal,oo_] = marginal_density(M_, options_, estim_params_, oo_);
            oo_ = GetPosteriorParametersStatistics(estim_params_, M_, options_, bayestopt_, oo_);
            oo_ = PlotPosteriorDistributions(estim_params_, M_, options_, bayestopt_, oo_);
        else
            load([M_.fname '_results'],'oo_');
        end
        metropolis_draw(1);
        if options_.bayesian_irf
            PosteriorIRF('posterior');
        end
        if options_.moments_varendo
            oo_ = compute_moments_varendo('posterior',options_,M_,oo_,var_list_);
        end
        if options_.smoother || ~isempty(options_.filter_step_ahead) || options_.forecast
            prior_posterior_statistics('posterior',data,gend,data_index,missing_value);
        end
        xparam = get_posterior_parameters('mean');
        set_all_parameters(xparam);
    end
end

if (~((any(bayestopt_.pshape > 0) && options_.mh_replic) || (any(bayestopt_.pshape ...
                                                      > 0) && options_.load_mh_file)) ...
    || ~options_.smoother ) && M_.endo_nbr^2*gend < 1e7 ...
         && options_.usePartInfo==0 && options_.partial_information == 0 % to be fixed   
    %% ML estimation, or posterior mode without metropolis-hastings or metropolis without bayesian smooth variable
    [atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,T,R,P,PK,decomp] = DsgeSmoother(xparam1,gend,data,data_index,missing_value);
    oo_.Smoother.SteadyState = ys;
    oo_.Smoother.TrendCoeffs = trend_coeff;
    oo_.Smoother.variance = P;
    i_endo = bayestopt_.smoother_saved_var_list;
    if options_.nk ~= 0
        oo_.FilteredVariablesKStepAhead = aK(options_.filter_step_ahead, ...
                                             i_endo,:);
        if isfield(options_,'kalman_algo')
            if ~isempty(PK)
                oo_.FilteredVariablesKStepAheadVariances = ...
                    PK(options_.filter_step_ahead,i_endo,i_endo,:);
            end
            if ~isempty(decomp)
                oo_.FilteredVariablesShockDecomposition = ...
                    decomp(options_.filter_step_ahead,i_endo,:,:);
            end
        end
    end
    for i=bayestopt_.smoother_saved_var_list'
        i1 = dr.order_var(bayestopt_.smoother_var_list(i));
        eval(['oo_.SmoothedVariables.' deblank(M_.endo_names(i1,:)) ' = atT(i,:)'';']);
        eval(['oo_.FilteredVariables.' deblank(M_.endo_names(i1,:)) ' = squeeze(aK(1,i,:));']);
        eval(['oo_.UpdatedVariables.' deblank(M_.endo_names(i1,:)) ...
              ' = updated_variables(i,:)'';']);
    end
    [nbplt,nr,nc,lr,lc,nstar] = pltorg(M_.exo_nbr);
    if options_.TeX
        fidTeX = fopen([M_.fname '_SmoothedShocks.TeX'],'w');
        fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation.m (Dynare).\n');
        fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
        fprintf(fidTeX,' \n');
    end    
    if nbplt == 1
        hh = figure('Name','Smoothed shocks');
        NAMES = [];
        if options_.TeX, TeXNAMES = []; end
        for i=1:M_.exo_nbr
            subplot(nr,nc,i);
            plot(1:gend,innov(i,:),'-k','linewidth',1)
            hold on
            plot([1 gend],[0 0],'-r','linewidth',.5)
            hold off
            xlim([1 gend])
            name    = deblank(M_.exo_names(i,:));
            if isempty(NAMES)
                NAMES = name;
            else
                NAMES = char(NAMES,name);
            end
            if ~isempty(options_.XTick)
                set(gca,'XTick',options_.XTick)
                set(gca,'XTickLabel',options_.XTickLabel)
            end
            if options_.TeX
                texname = M_.exo_names_tex(i,1);
                if isempty(TeXNAMES)
                    TeXNAMES = ['$ ' deblank(texname) ' $'];
                else
                    TeXNAMES = char(TeXNAMES,['$ ' deblank(texname) ' $']);
                end
            end
            title(name,'Interpreter','none')
            eval(['oo_.SmoothedShocks.' deblank(M_.exo_names(i,:)) ' = innov(i,:)'';']);
        end
        eval(['print -depsc2 ' M_.fname '_SmoothedShocks' int2str(1) '.eps']);
        if ~exist('OCTAVE_VERSION')
            eval(['print -dpdf ' M_.fname '_SmoothedShocks' int2str(1)]);
            saveas(hh,[M_.fname '_SmoothedShocks' int2str(1) '.fig']);
        end
        if options_.nograph, close(hh), end
        if options_.TeX
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            for jj = 1:M_.exo_nbr
                fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
            end
            fprintf(fidTeX,'\\centering \n');
            fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedShocks%s}\n',M_.fname,int2str(1));
            fprintf(fidTeX,'\\caption{Smoothed shocks.}');
            fprintf(fidTeX,'\\label{Fig:SmoothedShocks:%s}\n',int2str(1));
            fprintf(fidTeX,'\\end{figure}\n');
            fprintf(fidTeX,'\n');
            fprintf(fidTeX,'%% End of TeX file.\n');
            fclose(fidTeX);
        end
    else
        for plt = 1:nbplt-1
            hh = figure('Name','Smoothed shocks');
            set(0,'CurrentFigure',hh)
            NAMES = [];
            if options_.TeX, TeXNAMES = []; end
            for i=1:nstar
                k = (plt-1)*nstar+i;
                subplot(nr,nc,i);
                plot([1 gend],[0 0],'-r','linewidth',.5)
                hold on
                plot(1:gend,innov(k,:),'-k','linewidth',1)
                hold off
                name = deblank(M_.exo_names(k,:));
                if isempty(NAMES)
                    NAMES = name;
                else
                    NAMES = char(NAMES,name);
                end
                if ~isempty(options_.XTick)
                    set(gca,'XTick',options_.XTick)
                    set(gca,'XTickLabel',options_.XTickLabel)
                end
                xlim([1 gend])
                if options_.TeX
                    texname = M_.exo_names_tex(k,:);
                    if isempty(TeXNAMES)
                        TeXNAMES = ['$ ' deblank(texname) ' $'];
                    else
                        TeXNAMES = char(TeXNAMES,['$ ' deblank(texname) ' $']);
                    end
                end    
                title(name,'Interpreter','none')
                eval(['oo_.SmoothedShocks.' deblank(name) ' = innov(k,:)'';']);
            end
            eval(['print -depsc2 ' M_.fname '_SmoothedShocks' int2str(plt) '.eps']);
            if ~exist('OCTAVE_VERSION')
                eval(['print -dpdf ' M_.fname '_SmoothedShocks' int2str(plt)]);
                saveas(hh,[M_.fname '_SmoothedShocks' int2str(plt) '.fig']);
            end
            if options_.nograph, close(hh), end
            if options_.TeX
                fprintf(fidTeX,'\\begin{figure}[H]\n');
                for jj = 1:nstar
                    fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
                end    
                fprintf(fidTeX,'\\centering \n');
                fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedShocks%s}\n',M_.fname,int2str(plt));
                fprintf(fidTeX,'\\caption{Smoothed shocks.}');
                fprintf(fidTeX,'\\label{Fig:SmoothedShocks:%s}\n',int2str(plt));
                fprintf(fidTeX,'\\end{figure}\n');
                fprintf(fidTeX,'\n');
            end    
        end
        hh = figure('Name','Smoothed shocks');
        set(0,'CurrentFigure',hh)
        NAMES = [];
        if options_.TeX, TeXNAMES = []; end
        for i=1:M_.exo_nbr-(nbplt-1)*nstar
            k = (nbplt-1)*nstar+i;
            if lr ~= 0
                subplot(lr,lc,i);
            else
                subplot(nr,nc,i);
            end    
            plot([1 gend],[0 0],'-r','linewidth',0.5)
            hold on
            plot(1:gend,innov(k,:),'-k','linewidth',1)
            hold off
            name     = deblank(M_.exo_names(k,:));
            if isempty(NAMES)
                NAMES = name;
            else
                NAMES = char(NAMES,name);
            end
            if ~isempty(options_.XTick)
                set(gca,'XTick',options_.XTick)
                set(gca,'XTickLabel',options_.XTickLabel)
            end
            xlim([1 gend])
            if options_.TeX
                texname  = M_.exo_names_tex(k,:);
                if isempty(TeXNAMES)
                    TeXNAMES = ['$ ' deblank(texname) ' $'];
                else
                    TeXNAMES = char(TeXNAMES,['$ ' deblank(texname) ' $']);
                end
            end
            title(name,'Interpreter','none')
            eval(['oo_.SmoothedShocks.' deblank(name) ' = innov(k,:)'';']);
        end
        eval(['print -depsc2 ' M_.fname '_SmoothedShocks' int2str(nbplt) '.eps']);
        if ~exist('OCTAVE_VERSION')
            eval(['print -dpdf ' M_.fname '_SmoothedShocks' int2str(nbplt)]);
            saveas(hh,[M_.fname '_SmoothedShocks' int2str(nbplt) '.fig']);
        end
        if options_.nograph, close(hh), end
        if options_.TeX
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            for jj = 1:size(NAMES,1);
                fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
            end    
            fprintf(fidTeX,'\\centering \n');
            fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedShocks%s}\n',M_.fname,int2str(nbplt));
            fprintf(fidTeX,'\\caption{Smoothed shocks.}');
            fprintf(fidTeX,'\\label{Fig:SmoothedShocks:%s}\n',int2str(nbplt));
            fprintf(fidTeX,'\\end{figure}\n');
            fprintf(fidTeX,'\n');
            fprintf(fidTeX,'%% End of TeX file.\n');
            fclose(fidTeX);
        end    
    end
    %%
    %%  Smooth observational errors...
    %%
    if options_.noconstant
        yf = zeros(n_varobs,gend);
    else
        if options_.prefilter == 1
            yf = atT(bayestopt_.mf,:)+repmat(bayestopt_.mean_varobs,1,gend);
        elseif options_.loglinear == 1
            yf = atT(bayestopt_.mf,:)+repmat(log(ys(bayestopt_.mfys)),1,gend)+...
                 trend_coeff*[1:gend];
        else
            yf = atT(bayestopt_.mf,:)+repmat(ys(bayestopt_.mfys),1,gend)+...
                 trend_coeff*[1:gend];
        end
    end
    if nvn
        number_of_plots_to_draw = 0;
        index = [];
        for i=1:n_varobs
            if max(abs(measurement_error(10:end))) > 0.000000001
                number_of_plots_to_draw = number_of_plots_to_draw + 1;
                index = cat(1,index,i);
            end
            eval(['oo_.SmoothedMeasurementErrors.' deblank(options_.varobs(i,:)) ...
                  ' = measurement_error(i,:)'';']);
        end
        [nbplt,nr,nc,lr,lc,nstar] = pltorg(number_of_plots_to_draw);
        if options_.TeX
            fidTeX = fopen([M_.fname '_SmoothedObservationErrors.TeX'],'w');
            fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation.m (Dynare).\n');
            fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
            fprintf(fidTeX,' \n');
        end
        if nbplt == 1
            hh = figure('Name','Smoothed observation errors');
            set(0,'CurrentFigure',hh)
            NAMES = [];
            if options_.TeX, TeXNAMES = []; end
            for i=1:number_of_plots_to_draw
                subplot(nr,nc,i);
                plot(1:gend,measurement_error(index(i),:),'-k','linewidth',1)
                hold on
                plot([1 gend],[0 0],'-r','linewidth',.5)
                hold off
                name    = deblank(options_.varobs(index(i),:));
                if isempty(NAMES)
                    NAMES = name;
                else
                    NAMES = char(NAMES,name);
                end
                if ~isempty(options_.XTick)
                    set(gca,'XTick',options_.XTick)
                    set(gca,'XTickLabel',options_.XTickLabel)
                end
                if options_.TeX
                    idx = strmatch(options_.varobs(indx(i),:),M_.endo_names,'exact');
                    texname = M_.endo_names_tex(idx,:);
                    if isempty(TeXNAMES)
                        TeXNAMES = ['$ ' deblank(texname) ' $'];
                    else
                        TeXNAMES = char(TeXNAMES,['$ ' deblank(texname) ' $']);
                    end
                end
                title(name,'Interpreter','none')
            end
            eval(['print -depsc2 ' M_.fname '_SmoothedObservationErrors' int2str(1) '.eps']);
            if ~exist('OCTAVE_VERSION')
                eval(['print -dpdf ' M_.fname '_SmoothedObservationErrors' int2str(1)]);
                saveas(hh,[M_.fname '_SmoothedObservationErrors' int2str(1) '.fig']);
            end
            if options_.nograph, close(hh), end
            if options_.TeX
                fprintf(fidTeX,'\\begin{figure}[H]\n');
                for jj = 1:number_of_plots_to_draw
                    fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
                end    
                fprintf(fidTeX,'\\centering \n');
                fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedObservationErrors%s}\n',M_.fname,int2str(1));
                fprintf(fidTeX,'\\caption{Smoothed observation errors.}');
                fprintf(fidTeX,'\\label{Fig:SmoothedObservationErrors:%s',int2str(1));
                fprintf(fidTeX,'\\end{figure}\n');
                fprintf(fidTeX,'\n');
                fprintf(fidTeX,'%% End of TeX file.\n');
                fclose(fidTeX);
            end
        else
            for plt = 1:nbplt-1
                hh = figure('Name','Smoothed observation errors');
                set(0,'CurrentFigure',hh)
                NAMES = [];
                if options_.TeX, TeXNAMES = []; end
                for i=1:nstar
                    k = (plt-1)*nstar+i;
                    subplot(nr,nc,i);
                    plot([1 gend],[0 0],'-r','linewidth',.5)
                    hold on
                    plot(1:gend,measurement_error(index(k),:),'-k','linewidth',1)
                    hold off
                    name = deblank(options_.varobs(index(k),:));
                    if isempty(NAMES)
                        NAMES = name;
                    else
                        NAMES = char(NAMES,name);
                    end
                    if ~isempty(options_.XTick)
                        set(gca,'XTick',options_.XTick)
                        set(gca,'XTickLabel',options_.XTickLabel)
                    end
                    if options_.TeX
                        idx = strmatch(options_.varobs(k),M_.endo_names,'exact');
                        texname = M_.endo_names_tex(idx,:);
                        if isempty(TeXNAMES)
                            TeXNAMES = ['$ ' deblank(texname) ' $'];
                        else
                            TeXNAMES = char(TeXNAMES,['$ ' deblank(texname) ' $']);
                        end
                    end
                    title(name,'Interpreter','none')
                end
                eval(['print -depsc2 ' M_.fname '_SmoothedObservationErrors' int2str(plt) '.eps']);
                if ~exist('OCTAVE_VERSION')
                    eval(['print -dpdf ' M_.fname '_SmoothedObservationErrors' int2str(plt)]);
                    saveas(hh,[M_.fname '_SmoothedObservationErrors' int2str(plt) '.fig']);
                end
                if options_.nograph, close(hh), end
                if options_.TeX
                    fprintf(fidTeX,'\\begin{figure}[H]\n');
                    for jj = 1:nstar
                        fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
                    end    
                    fprintf(fidTeX,'\\centering \n');
                    fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedObservationErrors%s}\n',M_.fname,int2str(plt));
                    fprintf(fidTeX,'\\caption{Smoothed observation errors.}');
                    fprintf(fidTeX,'\\label{Fig:SmoothedObservationErrors:%s}\n',int2str(plt));
                    fprintf(fidTeX,'\\end{figure}\n');
                    fprintf(fidTeX,'\n');
                end    
            end
            hh = figure('Name','Smoothed observation errors');
            set(0,'CurrentFigure',hh)
            NAMES = [];
            if options_.TeX, TeXNAMES = []; end
            for i=1:number_of_plots_to_draw-(nbplt-1)*nstar
                k = (nbplt-1)*nstar+i;
                if lr ~= 0
                    subplot(lr,lc,i);
                else
                    subplot(nr,nc,i);
                end    
                plot([1 gend],[0 0],'-r','linewidth',0.5)
                hold on
                plot(1:gend,measurement_error(index(k),:),'-k','linewidth',1)
                hold off
                name     = deblank(options_.varobs(index(k),:));
                if isempty(NAMES)
                    NAMES = name;
                else
                    NAMES = char(NAMES,name);
                end
                if ~isempty(options_.XTick)
                    set(gca,'XTick',options_.XTick)
                    set(gca,'XTickLabel',options_.XTickLabel)
                end
                if options_.TeX
                    idx = strmatch(options_.varobs(index(k)),M_.endo_names,'exact');
                    texname = M_.endo_names_tex(idx,:);
                    if isempty(TeXNAMES)
                        TeXNAMES = ['$ ' deblank(texname) ' $'];
                    else
                        TeXNAMES = char(TeXNAMES,['$ ' deblank(texname) ' $']);
                    end
                end
                title(name,'Interpreter','none');
            end
            eval(['print -depsc2 ' M_.fname '_SmoothedObservationErrors' int2str(nbplt) '.eps']);
            if ~exist('OCTAVE_VERSION')
                eval(['print -dpdf ' M_.fname '_SmoothedObservationErrors' int2str(nbplt)]);
                saveas(hh,[M_.fname '_SmoothedObservationErrors' int2str(nbplt) '.fig']);
            end
            if options_.nograph, close(hh), end
            if options_.TeX
                fprintf(fidTeX,'\\begin{figure}[H]\n');
                for jj = 1:size(NAMES,1);
                    fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
                end    
                fprintf(fidTeX,'\\centering \n');
                fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedObservedErrors%s}\n',M_.fname,int2str(nbplt));
                fprintf(fidTeX,'\\caption{Smoothed observed errors.}');
                fprintf(fidTeX,'\\label{Fig:SmoothedObservedErrors:%s}\n',int2str(nbplt));
                fprintf(fidTeX,'\\end{figure}\n');
                fprintf(fidTeX,'\n');
                fprintf(fidTeX,'%% End of TeX file.\n');
                fclose(fidTeX);
            end    
        end
    end 
    %%
    %%  Historical and smoothed variabes
    %%
    [nbplt,nr,nc,lr,lc,nstar] = pltorg(n_varobs);
    if options_.TeX
        fidTeX = fopen([M_.fname '_HistoricalAndSmoothedVariables.TeX'],'w');
        fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation.m (Dynare).\n');
        fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
        fprintf(fidTeX,' \n');
    end    
    if nbplt == 1
        hh = figure('Name','Historical and smoothed variables');
        NAMES = [];
        if options_.TeX, TeXNAMES = []; end
        for i=1:n_varobs
            subplot(nr,nc,i);
            plot(1:gend,yf(i,:),'-r','linewidth',1)
            hold on
            plot(1:gend,rawdata(:,i),'-k','linewidth',1)
            hold off
            name    = deblank(options_.varobs(i,:));
            if isempty(NAMES)
                NAMES = name;
            else
                NAMES = char(NAMES,name);
            end
            if ~isempty(options_.XTick)
                set(gca,'XTick',options_.XTick)
                set(gca,'XTickLabel',options_.XTickLabel)
            end
            xlim([1 gend])
            if options_.TeX
                idx = strmatch(options_.varobs(i,:),M_.endo_names,'exact');
                texname = M_.endo_names_tex(idx,:);
                if isempty(TeXNAMES)
                    TeXNAMES = ['$ ' deblank(texname) ' $'];
                else
                    TeXNAMES = char(TeXNAMES,['$ ' deblank(texname) ' $']);
                end
            end
            title(name,'Interpreter','none')
        end
        eval(['print -depsc2 ' M_.fname '_HistoricalAndSmoothedVariables' int2str(1) '.eps']);
        if ~exist('OCTAVE_VERSION')
            eval(['print -dpdf ' M_.fname '_HistoricalAndSmoothedVariables' int2str(1)]);
            saveas(hh,[M_.fname '_HistoricalAndSmoothedVariables' int2str(1) '.fig']);
        end
        if options_.nograph, close(hh), end
        if options_.TeX
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            for jj = 1:n_varobs
                fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
            end    
            fprintf(fidTeX,'\\centering \n');
            fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_HistoricalAndSmoothedVariables%s}\n',M_.fname,int2str(1));
            fprintf(fidTeX,'\\caption{Historical and smoothed variables.}');
            fprintf(fidTeX,'\\label{Fig:HistoricalAndSmoothedVariables:%s}\n',int2str(1));
            fprintf(fidTeX,'\\end{figure}\n');
            fprintf(fidTeX,'\n');
            fprintf(fidTeX,'%% End of TeX file.\n');
            fclose(fidTeX);
        end    
    else
        for plt = 1:nbplt-1
            hh = figure('Name','Historical and smoothed variables');
            set(0,'CurrentFigure',hh)
            NAMES = [];
            if options_.TeX, TeXNAMES = []; end
            for i=1:nstar
                k = (plt-1)*nstar+i;
                subplot(nr,nc,i);
                plot(1:gend,yf(k,:),'-r','linewidth',1)
                hold on
                plot(1:gend,rawdata(:,k),'-k','linewidth',1)
                hold off
                name = deblank(options_.varobs(k,:));
                if isempty(NAMES)
                    NAMES = name;
                else
                    NAMES = char(NAMES,name);
                end
                if ~isempty(options_.XTick)
                    set(gca,'XTick',options_.XTick)
                    set(gca,'XTickLabel',options_.XTickLabel)
                end
                xlim([1 gend])
                if options_.TeX
                    idx = strmatch(options_.varobs(k,:),M_.endo_names,'exact');
                    texname = M_.endo_names_tex(idx,:);
                    if isempty(TeXNAMES)
                        TeXNAMES = ['$ ' deblank(texname) ' $'];
                    else
                        TeXNAMES = char(TeXNAMES,['$ ' deblank(texname) ' $']);
                    end
                end    
                title(name,'Interpreter','none')
            end
            eval(['print -depsc2 ' M_.fname '_HistoricalAndSmoothedVariables' int2str(plt) '.eps']);
            if ~exist('OCTAVE_VERSION')
                eval(['print -dpdf ' M_.fname '_HistoricalAndSmoothedVariables' int2str(plt)]);
                saveas(hh,[M_.fname '_HistoricalAndSmoothedVariables' int2str(plt) '.fig']);
            end
            if options_.nograph, close(hh), end
            if options_.TeX
                fprintf(fidTeX,'\\begin{figure}[H]\n');
                for jj = 1:nstar
                    fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
                end    
                fprintf(fidTeX,'\\centering \n');
                fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_HistoricalAndSmoothedVariables%s}\n',M_.fname,int2str(plt));
                fprintf(fidTeX,'\\caption{Historical and smoothed variables.}');
                fprintf(fidTeX,'\\label{Fig:HistoricalAndSmoothedVariables:%s}\n',int2str(plt));
                fprintf(fidTeX,'\\end{figure}\n');
                fprintf(fidTeX,'\n');
            end    
        end
        hh = figure('Name','Historical and smoothed variables');
        set(0,'CurrentFigure',hh)
        NAMES = [];
        if options_.TeX, TeXNAMES = []; end
        for i=1:n_varobs-(nbplt-1)*nstar
            k = (nbplt-1)*nstar+i;
            if lr ~= 0
                subplot(lr,lc,i);
            else
                subplot(nr,nc,i);
            end    
            plot(1:gend,yf(k,:),'-r','linewidth',1)
            hold on
            plot(1:gend,rawdata(:,k),'-k','linewidth',1)
            hold off
            name = deblank(options_.varobs(k,:));
            if isempty(NAMES)
                NAMES = name;
            else
                NAMES = char(NAMES,name);
            end
            if ~isempty(options_.XTick)
                set(gca,'XTick',options_.XTick)
                set(gca,'XTickLabel',options_.XTickLabel)
            end
            xlim([1 gend])
            if options_.TeX
                idx = strmatch(options_.varobs(i,:),M_.endo_names,'exact');
                texname = M_.endo_names_tex(idx,:);
                if isempty(TeXNAMES)
                    TeXNAMES = ['$ ' deblank(texname) ' $'];
                else
                    TeXNAMES = char(TeXNAMES,['$ ' deblank(texname) ' $']);
                end
            end
            title(name,'Interpreter','none');
        end
        eval(['print -depsc2 ' M_.fname '_HistoricalAndSmoothedVariables' int2str(nbplt) '.eps']);
        if ~exist('OCTAVE_VERSION')
            eval(['print -dpdf ' M_.fname '_HistoricalAndSmoothedVariables' int2str(nbplt)]);
            saveas(hh,[M_.fname '_HistoricalAndSmoothedVariables' int2str(nbplt) '.fig']);
        end
        if options_.nograph, close(hh), end
        if options_.TeX
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            for jj = 1:size(NAMES,1);
                fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
            end    
            fprintf(fidTeX,'\\centering \n');
            fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_HistoricalAndSmoothedVariables%s}\n',M_.fname,int2str(nbplt));
            fprintf(fidTeX,'\\caption{Historical and smoothed variables.}');
            fprintf(fidTeX,'\\label{Fig:HistoricalAndSmoothedVariables:%s}\n',int2str(nbplt));
            fprintf(fidTeX,'\\end{figure}\n');
            fprintf(fidTeX,'\n');
            fprintf(fidTeX,'%% End of TeX file.\n');
            fclose(fidTeX);
        end
    end
end

if options_.forecast > 0 && options_.mh_replic == 0 && ~options_.load_mh_file 
    forecast(var_list_,'smoother');
end

if np > 0
    pindx = estim_params_.param_vals(:,1);
    save([M_.fname '_pindx.mat'] ,'pindx');
end

