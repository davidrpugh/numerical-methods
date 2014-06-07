function [xhatnew,zhatnew,signew,lh,yhat,e]=pt_info_kf(y,zhat,sig)
% Modified by Joe April 2011
% Kalman Filter for Part Information estimation 
% function [xhatnew,zhatnew,signew,lh,yhat,e]=pt_info_kf(y,zhat,sig)
% it has also: nmat,AA,BB,CCC,DD,EE,FF,I_L,UFT,VKLUFT,ALUFT,QQ,RR,QU,BUFT,FUBT,FULKV)
% passed in global structure kf_param.
% the plant equation is z(t)=CCC*z(t-1)+Me, where e is N(0,I).  
% The observation equation is y(t)=AA*s(t).  
% The prior distribution for z is N(zhat,sig).  
% The posterior on the new state is N(zhatnew,signew) and lh is a two-dimensional vector containing
% the increments to the two component terms of the log likelihood function.  
% 
% June 06: written by G Perendia from notes by Prof. Joe Pearlman and original KF.m by A Justiniano
% to suit partial information, RE solution in accordance with, and based on, 
% the Pearlman, Currie and Levine 1986 solution.
% 
%	 	nmat (i.e. N matrix) = Z22'\Z12'

%            %pre calculated time-invariant factors - by the calling function
%            A11=G1(1:pd,1:pd);
%            A22=G1(pd+1:end, pd+1:end);
%            A12=G1(1:pd, pd+1:end);
%            A21=G1(pd+1:end,1:pd);
%            BB=A12*inv(A22);
%            K2A22i=K2*inv(A22);
%            QQ=U11;

%		precalculated members of globally defined structure kf_param:
%            AA=A11-BB*A21;
%            FF=A11-A12*nmat;
%            HH=K1-K2A22i*A21;
%            EE=K1-K2*nmat;
%            RR=K2A22i*UFT+VV;
%            SS=BB*UFT;
%            VKLUFT=VV+K2*I_L*UFT;
%            ALUFT=A12*I_L*UFT;
%            FULKV=VV;
%            nmat


    global kf_param;

    HH=kf_param.HH;
    PP=kf_param.PP;
    EE=kf_param.EE;
    yhat=kf_param.EE*zhat;
    e=y-yhat;

    % Def. Esig= (EPE'+V) i.e.
    Esig=EE*sig*EE'+kf_param.VV; 
    if any(Esig)
        try
            I_ESIG=inv(Esig); 
        catch
            lh=Inf;
            xhatnew=[];signew=[];zhatnew=[];
            return;
        end
    else
        I_ESIG=Inf;
    end
    FF=kf_param.FF;
    CDFT=FF*sig*HH';
    if ~any(CDFT)| ~any(e)
        zhatnew= FF*zhat;
    elseif  isinf(I_ESIG)% == Inf
        %error ('Inverting Zero or Singular matrix');
        warning ('Inverting Zero or Singular matrix');
        lh=Inf;
        xhatnew=[];signew=[];zhatnew=[];
        lastwarn('')
        return
    else
        % Old Setup zhatnew= FF*zhat+ CDFT*I_ESIG*e;
        %New setup
       % zhatnew= FF*zhat+ FF*sig*HH'*I_ESIG*e; 
        zhatnew= FF*zhat+ FF*sig*EE'*I_ESIG*e; 
    end

    % Def. DsigR=(DPD'+FU22F'+V), i.e.
    DsigR=HH*PP*HH'+kf_param.RR; 
    if any(DsigR) % i.e. not null matrix
        I_DSIG=inv(DsigR); 
    else
        I_DSIG=Inf;
    end

    % New: Zt+1 = FZtF' + PH' (HPH' + V )^(-1)HP - FZtE' (EZtE' + V )^(-1)EZtF'     
       signew= FF*sig*FF' + PP*HH' * I_DSIG* HH*PP - ...
          FF*sig*EE'*inv(EE*sig*EE')*EE*sig*FF' ;% +  kf_param.QQ;
    

    xhatnew=-kf_param.nmat*zhatnew;
    
  %  DEKV=(HH*sig*kf_param.EE'+kf_param.FULKV);
  %  if ~any(DEKV) | ~any(Esig) %==0
  %      CovE= zeros(size(DEKV*Esig));
  %  elseif  isinf(I_DSIG)% == Inf
  %      warning ('Inverting Zero or Singular matrix');
  %      lh=Inf;
  %      xhatnew=[];signew=[];zhatnew=[];
  %      lastwarn('')
  %      return
  %  else
         CovE= Esig; % *I_DSIG*DEKV; % Mistake CovE= DEKV*I_DSIG*Esig;
  %  end
    invCovE=inv(CovE);
    [LastWarningTxt LastWarningID]=lastwarn;
    if strcmp('MATLAB:nearlySingularMatrix',LastWarningID) | ...
        strcmp('MATLAB:illConditionedMatrix',LastWarningID) | ...
        strcmp('MATLAB:singularMatrix',LastWarningID)
     %   disp(['PI_KF:  ' LastWarningTxt]);
        lh=Inf;
        xhatnew=[];signew=[];zhatnew=[];
        lastwarn('')
        return
    end 
    if ~any(e) % i.e. e==0
        lh= -log(det(CovE));
    elseif  ~any(CovE) %==Inf
    %     danger of Inverting Zero  matrix');
    %    lh= - 1e10; % very small
            lh=Inf;
    elseif  isinf(invCovE) %==Inf
    %    error ('Inverting Zero or Singular matrix');
        warning ('Inverting Zero or Singular matrix CovE:');
    %    lh= - det(CovE);
    %    lh= - 1e10
            lh=Inf;
    else
        lh= -log(det(CovE)) - (e'*invCovE*e) ;
    end
    lh=lh/2;
    if isnan(lh) | isinf(lh)
          warning (' likelihood not a number: %d ', lh );
    end 
    
