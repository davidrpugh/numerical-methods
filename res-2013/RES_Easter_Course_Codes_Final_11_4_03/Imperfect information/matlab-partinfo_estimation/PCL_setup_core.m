function  [PDIDPDRD,PP,CCCC, AA]=PCL_setup_core ( H, varobs)
% sets up parameters and calls part-info kalman filter
% developed by G Perendia, July 2006 for implementation from notes by Prof. Joe Pearlman to 
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
    [junk,OBS] = ismember(varobs,M_.endo_names,'rows');
    NOBS = length(OBS);
    
    G1=dr.PI_ghx;
    impact=dr.PI_ghu;
    nmat=dr.PI_nmat;
    CC=dr.PI_CC;
    NX=M_.exo_nbr; % no of exogenous varexo shock variables.
    FL_RANK=dr.PI_FL_RANK;
    NY=M_.endo_nbr;
    LL = sparse(1:NOBS,OBS,ones(NOBS,1),NY,NY);

    if exist( 'irfpers')==1
        if ~isempty(irfpers)
            if irfpers<=0, irfpers=20, end;
        else
            irfpers=20;
        end
     else
        irfpers=20;
    end      

    ss=size(G1,1);

    pd=ss-size(nmat,1);
    SDX=M_.Sigma_e^0.5; % =SD,not V-COV, of Exog shocks or M_.Sigma_e^0.5 num_exog x num_exog matrix
    if isempty(H)
        H=M_.H;
    end
    VV=H; % V-COV of observation errors.
    MM=impact*SDX; % R*(Q^0.5) in standard KF notation
    % observation vector indices
    % mapping to endogenous variables.

    L1=LL*dr.PI_TT1;
    L2=LL*dr.PI_TT2;

    MM1=MM(1:ss-FL_RANK,:);
    U11=MM1*MM1';
    % SDX
    U22=0;
    % determine K1 and K2 observation mapping matrices
    % This uses the fact that measurements are given by L1*s(t)+L2*x(t)
    % and s(t) is expressed in the dynamics as
    % H1*eps(t)+G11*s(t-1)+G12*x(t-1)+G13*x(t).  
    % Thus the observations o(t) can be written in the form
    % o(t)=K1*[eps(t)' s(t-1)' x(t-1)']' + K2*x(t) where
    % K1=[L1*H1 L1*G11 L1*G12] K2=L1*G13+L2

    G12=G1(NX+1:ss-2*FL_RANK,:);
    KK1=L1*G12;
    K1=KK1(:,1:ss-FL_RANK);
    K2=KK1(:,ss-FL_RANK+1:ss)+L2;

    %pre calculate time-invariant factors
    A11=G1(1:pd,1:pd);
    A22=G1(pd+1:end, pd+1:end);
    A12=G1(1:pd, pd+1:end);
    A21=G1(pd+1:end,1:pd);
    Lambda= nmat*A12+A22;
    I_L=inv(Lambda);
    BB=A12*inv(A22);
    FF=K2*inv(A22);       
    QQ=BB*U22*BB' + U11;        
    UFT=U22*FF';
    % kf_param structure:
    AA=A11-BB*A21;
    CCCC=A11-A12*nmat; % F in new notation
    DD=K1-FF*A21; % H in new notation
    EE=K1-K2*nmat;
    RR=FF*UFT+VV;
    if ~any(RR) 
        % if zero add some dummy measurement err. variance-covariances 
        % with diagonals 0.000001. This would not be needed if we used
        % the slow solver, or the generalised eigenvalue approach, 
        % but these are both slower.
        RR=eye(size(RR,1))*1.0e-6;
    end
    SS=BB*UFT;
    VKLUFT=VV+K2*I_L*UFT;
    ALUFT=A12*I_L*UFT;
    FULKV=FF*U22*I_L'*K2'+VV;
    FUBT=FF*U22*BB';
    nmat=nmat;
    % initialise pshat
    AQDS=AA*QQ*DD'+SS;
    DQDR=DD*QQ*DD'+RR;
    I_DQDR=inv(DQDR);
    AQDQ=AQDS*I_DQDR;
    ff=AA-AQDQ*DD;
    hh=AA*QQ*AA'-AQDQ*AQDS';%*(DD*QQ*AA'+SS');
    rr=DD*QQ*DD'+RR;
    ZSIG0=disc_riccati_fast(ff,DD,rr,hh);
    PP=ZSIG0 +QQ;

    exo_names=M_.exo_names(M_.exo_names_orig_ord,:);

    DPDR=DD*PP*DD'+RR;
    I_DPDR=inv(DPDR);
    PDIDPDRD=PP*DD'*I_DPDR*DD;
