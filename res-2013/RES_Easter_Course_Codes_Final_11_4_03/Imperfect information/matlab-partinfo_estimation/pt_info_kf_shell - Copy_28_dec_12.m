function [xhatnew,shatnew,signew,LIK,lht]=pt_info_kf_shell( Y,CC,SDX,G1,VV,MM,nmat,zmat,L1,L2,shat,T,NY,NX,FL_RANK,start)
% Modified by Joe April 2011
% sets up parameters and calls part-info kalman filter
% developed by G Perendia, July 2006 for implementation from notes by Prof. Joe Pearlman to 
% suit partial information RE solution in accordance with, and based on, the 
% Pearlman, Currie and Levine 1986 solution.
% 22/10/06 - Version 2 for new Riccati with 4 params instead 5 

        % Recall that the state space is given by the 
        % predetermined variables s(t-1), x(t-1) 
        % and the jump variables x(t).
        % The jump variables have dimension NETA
        
            global kf_param;
            global options_;
            global M_;
		ss=size(G1,1);
        
        pd=ss-size(nmat,1);
        
        %MM1=MM(1:ss-FL_RANK,:);
        %U11=MM1*MM1';
        % SDX
	
        %U11=SDX*SDX; % SDX - Cov(error terms)
        %U11=[U11 zeros(size(U11,1),pd-size(U11,1));zeros(pd-size(U11,1),pd)];
        U22=0;
        % determine K1 and K2 observation mapping matrices
        % This uses the fact that measurements are given by L1*s(t)+L2*x(t)
        % and s(t) is expressed in the dynamics as
        % H1*eps(t)+G11*s(t-1)+G12*x(t-1)+G13*x(t).  
        % Thus the observations o(t) can be written in the form
        % o(t)=K1*[eps(t)' s(t-1)' x(t-1)']' + K2*x(t) where
        % K1=[L1*H1 L1*G11 L1*G12] K2=L1*G13+L2
        
 % New definition       
		G12=G1(NX+1:ss-2*FL_RANK,:);
 % End of new definition       
		KK1=L1*G12;
        K1=KK1(:,1:ss-FL_RANK);
        K2=KK1(:,ss-FL_RANK+1:ss)+L2;
                                                %KK=[KK1;L2];
                                                %K1=KK(:,1:ss-NETA);
                                                %K2=KK(:,ss-NETA+1:end);

                                                
        % red_form added BY to fix Ricc
        E2=[K1,K2];
        E5=zeros(size(L1,1),0);
        DD=zeros(size(G1,2),0);
        W=zeros(size(G1,2));
        num_inst=0;
        ACES_M=size(G1,2)-FL_RANK;
        disp_flag=0;
        [G1,HB,E2,junk,DD,MM,FL_RANK,ACES_M]=red_form_PIest(G1,MM,DD,E2,E5,W,num_inst,FL_RANK,ACES_M,M_,options_,disp_flag);
      
        ss=size(G1,1); 
        K1=E2(:,1:ss-FL_RANK);
        K2=E2(:,ss-FL_RANK+1:ss);
       
        MM1=MM(1:ss-FL_RANK,:);
        U11=MM1*MM1';
       
        [m,lam] = eig(G1.');                                         
        m=conj(m);
        mT=m';       
        
        mtemp=mT(find(diag(lam)>1),:);
        
        m1=mtemp(:,1:ss-FL_RANK);
        m2=mtemp(:,ss-FL_RANK+1:ss);
        nmat=inv(m2)*m1;
        nmat=real(nmat);
        pd=ss-size(nmat,1);
        
        %pre calculate time-invariant factors
        A11=G1(1:pd,1:pd);
        A22=G1(pd+1:end, pd+1:end);
        A12=G1(1:pd, pd+1:end);
        A21=G1(pd+1:end,1:pd);
        %A11_A12Nmat= A11-A12*nmat % test 
        BB=A12*inv(A22);
        K2A22i=K2*inv(A22);       
        QQ=U11;      % =U11 
        % kf_param structure:
        kf_param.AA=A11-BB*A21;
        kf_param.FF=A11-A12*nmat;
        kf_param.HH=K1-K2A22i*A21;
        HH=kf_param.HH;
        kf_param.EE=K1-K2*nmat;
        kf_param.RR=VV;  %=VV
        %if ~any(kf_param.RR) 
            % if zero add some dummy measurement err. variance-covariances 
            % with diagonals 0.000001. This would not be needed if we used
            % the slow solver, or the generalised eigenvalue approach, 
            % but these are both slower.
        %    kf_param.RR=eye(size(kf_param.RR,1))*0;% 1.0e-6;
        %end
        kf_param.VV=VV;  %=VV
        kf_param.FULKV=VV;
        kf_param.QQ=QQ;
        kf_param.nmat=nmat;

        % initialise pshat
        AQDS=kf_param.AA*QQ*HH';
        DQDR=HH*QQ*HH'+kf_param.RR;
        
        %RR=kf_param.RR
        %HH=kf_param.HH
        %U11
        %QQ
        %RRR = chol(DQDR)
        %d = eig(DQDR)
        
        smallnumber=1.0e-8;
    
        lastwarn('');
        %if any(DQDR) I_DQDR=inv(DQDR); else I_DQDR=Inf; end;
        if any(DQDR) 
          warning off all;
          I_DQDR=inv(DQDR);
          warning on all;
          if any(find(I_DQDR==Inf)) && any(~any(kf_param.RR))
            lastwarn('');
            kf_param.RR=eye(size(kf_param.RR,1))*1.0e-8;%1.0e-6;%1.0e-8;%0;%
            DQDR=HH*QQ*HH'+kf_param.RR;
            I_DQDR=inv(DQDR);
          end;
        elseif (~any(kf_param.RR))
          %kf_param.RR=eye(size(kf_param.RR,1))*1.0e-8;%1.0e-6;%1.0e-8;%0;%  
          %DQDR=HH*QQ*HH'+kf_param.RR;
           DQDR=HH*QQ*HH';
          for i=1:size(DQDR,1)
             if DQDR(i,i)<=smallnumber
                DQDR(i,i)=smallnumber;
             end
          end 
          I_DQDR=inv(DQDR);
        else
          I_DQDR=Inf
        end;
      
        %I_DQDR
        
        [LastWarningTxt LastWarningID]=lastwarn;
        if strcmp('MATLAB:nearlySingularMatrix',LastWarningID) | ...
           strcmp('MATLAB:illConditionedMatrix',LastWarningID) | ...
            strcmp('MATLAB:singularMatrix',LastWarningID)
        %    disp(['PI_KF:  ' LastWarningTxt]);
        
        %   SDD=SDX
        
            lht=Inf;
            xhatnew=[];signew=[];shatnew=[]; LIK=lht;
            lastwarn('');
            return
        end 
        %end
        if ~any(AQDS) %==0 bypass division by possibly non-invertble matrix
            AQDQ=0;
            ff=kf_param.AA;
            hh=kf_param.AA*QQ*kf_param.AA';%zeros(size(A11));
        elseif  isinf(I_DQDR)% e.g == Inf  i.e. could not bypass it!!
           % error ('Inverting Zero or Singular matrix');
            warning ('Inverting Zero or Singular matrix');
            lht=Inf;
            xhatnew=[];signew=[];shatnew=[]; LIK=lht;
            return
        else  
%           AQDQ=(kf_param.AA*QQ*kf_param.HH'+kf_param.SS)*inv(kf_param.HH*QQ*kf_param.HH'+kf_param.RR);
            AQDQ=AQDS*I_DQDR;
            ff=kf_param.AA-AQDQ*HH;
            hh=kf_param.AA*QQ*kf_param.AA'-AQDQ*AQDS';
        end
        try
            %New Riccati routine having 4 inputs F,D,R,H 
            %solving P=FPF' - FPD'inv(DPD'+rr)DPF' + H
            %where rr=H*Q*H'+R; - moved out from old Riccati, i.e.:
            %rr=HH*QQ*HH'+kf_param.RR;
            rr=HH*QQ*HH';
            for i=1:size(rr,1)
               if rr(i,i)<=smallnumber
                  rr(i,i)=smallnumber;
               end
            end 

            PP=disc_riccati_fast(ff,HH,rr,hh);
            %PP=kalman_steady_state(ff,kf_param.HH,rr,hh);%(T,QQ,Z,H)
            %PP=kalman_steady_state(ff,hh,kf_param.HH',rr)';
            % Test
            %ZSIG0_DIF=ZSIG0_OLD-ZSIG0;
            %if any(ZSIG0_DIF)
            %    ZSIG0_DIF
            %end
        catch
            try
                lerror=lasterror;
                disp(['PT_Info_KF_Shell: ' lerror.message]);
                if ~strfind(lerror.message,'Riccati not converged fast enough!')
                    disp '** Unexpected Error PT_Info_KF_Shell:disc_riccati_fast: ** :';
                    button=questdlg('Continue Y/N?',...
                        'Unexpected Error in PT_Info_KF_Shell','No','Yes','Yes'); 
                    switch button 
                    case 'No' 
                        error ('Terminated')
                    %case 'Yes'

                    end
                end
                lht=Inf;
                xhatnew=[];signew=[];shatnew=[]; LIK=lht;
                return
            catch
                disp '** Unexpected Error in PT_Info_KF_Shell:disc_riccati_fast ** :';
                disp lerror.message;
                button=questdlg('Continue Y/N?', ...
                    'Unexpected Error in PT_Info_KF_Shell','No','Yes','Yes'); 
                switch button 
                case 'No' 
                    error ('Terminated') 
                case 'Yes' 
                    lht=Inf;
                    xhatnew=[];signew=[];shatnew=[]; LIK=lht;
                    return
                end
            end
        end
        
        %P0 = Z+Q
        SIG0=PP +QQ; 
        kf_param.PP=SIG0;
        
        FDT=SIG0*HH';
        DSDR=HH*SIG0*HH'+kf_param.RR;
        lastwarn('')
        if ~any(DSDR) I_DSDR = Inf; else I_DSDR =inv(DSDR ); end
        [LastWarningTxt LastWarningID]=lastwarn;
        if strcmp('MATLAB:nearlySingularMatrix',LastWarningID) | ...
            strcmp('MATLAB:illConditionedMatrix',LastWarningID) | ...
            strcmp('MATLAB:singularMatrix',LastWarningID)
       %     disp(['PI_KF:  ' LastWarningTxt]);
            lht=Inf;
            xhatnew=[];signew=[];shatnew=[]; LIK=lht;
            lastwarn('');
            return
        end 
        if ~any(FDT)  %==0
            MSIG=disclyap_fast(kf_param.FF,zeros(size(kf_param.FF)));
        elseif isinf(I_DSDR)% == Inf
           % error ('Inverting Zero or Singular matrix');
            warning ('Inverting Zero or Singular matrix');
            lht=Inf;
            xhatnew=[];signew=[];shatnew=[]; LIK=lht;
            return
        else
            MSIG=disclyap_fast(kf_param.FF, FDT*I_DSDR*FDT');
        end
        if isnan(MSIG)
          xhatnew=[];signew=[];shatnew=[]; lht=Inf;LIK=Inf;
          display(' disclyap_fast returned NaN!');
          return;
        end
        pshat=MSIG;
        lht=zeros(T,1);
		xhatnew=zeros(FL_RANK,T);
        shat =shat(1:pd);
        ii=1;   
        for ii=1:T;
            WW=(Y(ii,:)'-CC);
            [xhatnewt,shatnew,signew,lik]=pt_info_kf( WW ,shat,pshat);
            lht(ii)=lik; 
            if lik==Inf 
              xhatnew=[];signew=[];shatnew=[]; LIK=Inf;
              return; 
            end
            shat=shatnew; 
            pshat=signew;
    %        if ii==50
    %            disp('pshat at 50')
    %            pshat
    %        end
            xhatnew(:,ii)=xhatnewt; 
        end
        %save('PI_LIKS','PIldf','PIlefe');
        %xhatnewt;
        %shatnew;
        %signew%(1:4,:)
        smpl= T; % =size(Y,1);
        pp   = size(Y,2);
        % likelihoods are negative
%        lht(smpl+1)= -0.5*smpl*pp*log(2*pi);
%        LIK=-(sum(lht(start:end))-(start-1)*lht(smpl+1)/smpl);
        % likelihoods are positive
%        lht(smpl+1)= 0.5*smpl*pp*log(2*pi);
%        LIK=(sum(lht(start:end))-(start-1)*lht(smpl+1)/smpl);
        penalty= 0.5*smpl*pp*log(2*pi);
        LIK=-sum(lht(start:end))+penalty -(start-1)*penalty/smpl;

        clear kf_param;
