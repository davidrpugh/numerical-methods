function [HAG,HB,CG,W_red,D_red,C_red,nm,mm] = red_form_PIest(aa,cc,dd,e2,e5,W,num_inst,nm,mm,M_,options_,disp_flag)

global lq_instruments;
% system defined as x(t+1)=aa*x(t)+BB*w(t); 
%                   y(t)=Chat*x(t)


BB=[dd, cc]; %define BB by augmenting D and C 
Chat=e2;

% reduce Chat
if isfield(lq_instruments,'optvar')
    num_optv=size(lq_instruments.optvar,1);
    for i=1:num_optv
        %find indices of variables of interest 
        i_optv(i) = strmatch(deblank(lq_instruments.optvar(i,:)),...
                    M_.endo_names(lq_instruments.m_var,:),'exact');
    end
    Chat=Chat(i_optv,:); %extract elements in E2 based on i_optv
    E5=e5(i_optv,:); %%extract elements in E5 based on i_optv
end
%if isempty(num_optv)
%        error('Need to specifiy the names of variables for computing the optimal policy, or the Ricatti interation may run into problems')
%end

% partition W - [QQ, UU; UU',RR]
i_varW=[];
    for m=1:num_inst
        i_varW(m)=size(e2,2)+m; %find indices of instrument(s)
    end
QQ=W([1:size(e2,2)],[1:size(e2,2)]); %(cross) terms of xt
UU=W([1:size(e2,2)],[i_varW:end]); %cross terms of xt/wt
% this is a test on setting cross terms to zeros
%UU=zeros(size(BB,1),num_inst);
RR=W(i_varW,i_varW); %terms of wt


%% Controllable form
max_iter=size(aa,2);
mTnew=[];
% testing using different tol
% one changes it in the mod file
%if isfield(lq_instruments,'red_tol')
   %tol=lq_instruments.red_tol;
%else
   tol=1e-15;
%end

penalty=0;

for n=1:max_iter
    % find the left eigenvectors of aa
    [m,lam] = eig(aa.');
    m=conj(m);
    mT=m';
    %additional condition to exit
    %penaltyold=zeros(n,1);
    %if penalty>0 & n>1 
    %    penaltyold(n)=penalty; % store pen
    %    if n-2>0 & (penaltyold(n)-penaltyold(n-1))==(penaltyold(n-1)-penaltyold(n-2)) %if change(pen)=constant
    %      if disp_flag    
    %        disp('system is now controllable');
    %      end
    %    break
    %    end
    %end
    %if penalty==0 & n>2
       temp=mT*BB; temp2=[];
       for i=1:size(aa,2)
           temp2(i)=any(abs(temp(i,:))>tol); %find if there are any zero rows in mT*BB
%        if isempty(temp2)
%           emp=i; 
%        end
       end 
%      if size(find(temp2),2)==size(aa,2) % only works for isreal(mT)
       if size(find(temp2(1:size(aa,1))),2)==size(aa,1) %if these is no more zero elements 
          if disp_flag
            disp('system is now controllable');
          end
       break
       end
    %end   
    %else    
        for i=1:size(mT,1)
            temp=mT(i,:)*BB; %examine each row of mT
            tempD=mT(i,:)*dd; %examine each row of mT re instrument
            if abs(tempD)<tol % check if mTD=0
               if isreal(temp)
                  mTnew=mT(i,:); %define mTnew=mT(i,:) for which mT*BB=0
                  n_new=size(mTnew,2); %current size of the system at each reduction
                  critval=n_new-nm; %theshold index distinguishing BL and FL vars
                  tem_i=find(abs(mTnew)>tol); %find the indices of nonzero elements in mT
                  %if all elements of tem_i do not fall into the BL region
                  %(<=critval) or FL region (>critval)
                  if isempty(find(all(tem_i<=critval)))&isempty(find(all(tem_i>critval)))
                     if abs(lam(i,i))<1 
                        penalty=penalty+1;
                     end      
                  elseif abs(temp)<tol %check if mTB=0 when all elements of tem_i do fall into either BL region
                                       %(<=critval) or FL region (>critval)   
                      if all(tem_i<=critval) %if tem_i fall into the BL region(<=critval)
                         mm=mm-1; %reduce BL
                         %now find mtrices F T inv(T) and D
                         F=eye(size(aa));
                         %pick the first nonzero in mTnew
                         [junk,j,val]=find(mTnew);
                         for i=1:size(j,2)
                            if abs(val(i))>tol & j(i)>M_.exo_nbr
                               j=j(i);
                            break
                            end
                         end   
                         F(j,:)=[]; T=[F;mTnew]; D=inv(T);
                         D(:,size(F,1)+size(mTnew,1))=[]; 
                         %find the controllable Q, U, A and B
                         QQ=D'*QQ*D; UU=(UU'*D)'; aa=F*aa*D; BB=F*BB;
                         %need to reduce dd and cc as well for ismember !!!!!!
                         dd=F*dd; cc=F*cc;
                         %matrix CD
                         Chat=Chat*D;
                         break 
                      else
                         nm=nm-1; %reduce FL
                         %now find mtrices F T inv(T) and D
                         F=eye(size(aa));
                         %pick the first nonzero in mTnew
                         [junk,j,val]=find(mTnew);
                         for i=1:size(j,2)
                            if abs(val(i))>tol & j(i)>M_.exo_nbr
                               j=j(i);
                            break
                            end
                         end    
                         F(j,:)=[]; T=[F;mTnew]; D=inv(T);
                         D(:,size(F,1)+size(mTnew,1))=[]; 
                         %find the controllable Q, U, A and B
                         QQ=D'*QQ*D; UU=(UU'*D)'; aa=F*aa*D; BB=F*BB;
                         %need to reduce dd and cc as well for ismember !!!!!!
                         dd=F*dd; cc=F*cc;
                         %matrix CD
                         Chat=Chat*D;
                         break 
                      end      
                  end
               else % if this temp is complex
                  mTnew=mT(i,:); 
                  m1T=real(mTnew);
                  m2T=imag(mTnew);
                  mTnew=[m1T;m2T];
                  F=eye(size(aa));
                  %j = find_jcomplex(mTnew,tol);% !!!!!!!!!!needs to be changed, see below
                  [j, mTnew]= find_jcomplex(mTnew, m1T, tol);
                  % reducing n, nm and m when complex
                  n_new=size(mTnew,2);
                  critval=n_new-nm;
                  tem_i=find(abs(mTnew)>tol);
                  [junk, tem_i]=find(abs(mTnew)>tol);
                  tem_i=unique(tem_i);
                  if isempty(find(all(tem_i<=critval)))|isempty(find(all(tem_i>critval)))&abs(lam(i,i))<1 
                     penalty=penalty+1;
                  end   
                  if abs(temp)<tol %mTB=0   
                      if all(tem_i<=critval) %if tem_i fall into the BL region(<=critval)
                         mm=mm-size(j,2);
                      elseif all(tem_i>critval)|(isempty(find(all(tem_i<=critval)))&isempty(find(all(tem_i>critval)))&abs(lam(i,i))>1)%if tem_i fall into the FL region(>critval)
                             nm=nm-size(j,2);
                      end 
                  elseif (isempty(find(all(tem_i<=critval)))&isempty(find(all(tem_i>critval)))&abs(lam(i,i))>1)          
                         nm=nm-size(j,2);    
                  %now find mtrices F T inv(T) and D
                 F(j,:)=[];
                 T=[F;mTnew];
                 D=inv(T);
                 % this is an error (for the complex numbers)
                 %D(:,size(F,1)+size(mTnew,1))=[];
                 D=D(:,1:size(F,1));
                 % find the controllable Q, U, A and B
                 QQ=D'*QQ*D;
                 UU=(UU'*D)';
                 aa=F*aa*D;
                 BB=F*BB;
                 % need to reduce dd and cc as well for ismember !!!!!!!
                 dd=F*dd;
                 cc=F*cc;
                 % matrix CD
                 Chat=Chat*D;
               break 
                  end
               end
           end    
        end
    %end  
end
if n==max_iter
    if disp_flag
       disp('Controllable reduction: max_iter is reached');
       penalty
    end
end    
% find red-form matrices FAD and FB (control)
FAD=aa;
FB=BB;
ChatD=Chat;
QQ_red=QQ;
UU_redT=UU';

%% Observable form - new version - 27/06/12
% define Cbar – stack C (optvar by n), Q (n by n) and U (number of instru by n)
Cbar=[ChatD;QQ_red;UU_redT];
% max_iter is the size of A after controllability (=n)
max_iter=size(FAD,2);
snew=[];
emp=[];
temp2=[];
r=0;
p=0;
ss=0;
nstate=[];
SB=[];
SF=[];

for n=1:max_iter
    % find the eigenvectors of A – s is a full matrix whose columns are the e-vectors (mu is the 
    %  e-values and is never used)
    [s,mu] = eig(FAD);
    % temp is Cbar*s - (optvar+n+number_of_instru) by n
    %temp=Cbar*s;
  
    %if n>1 && (isempty(SB) && isempty(SF))
    %     if disp_flag    
    %     disp('system is now observable');
    %     end
    %     break
    %end
    
    tempstate=size(nstate,2);
    if tempstate-2>0 && n<max_iter && (nstate(tempstate)-nstate(tempstate-1))==0 %if change(nstate)=constant 
          if disp_flag    
            disp('system is now observable');
          end
        break
    end 
    
    %------------------------------
    %find the BL and FL regions
    n_new=size(FAD,2);% n_new is the current size of A
    % critval is n_new-(number of FL var) so from 1 to critval (including crival) is the 
    % BL region and remaining is FL area
    critval=n_new-nm;           
    %------------------------------           
    for ii=1:size(s,2)
        temp3=Cbar*s(:,ii); % Cbar*each column of s (starting from the ii’th)
        snew=s(:,ii); % record the ii’th column of s (i.e. ii’th e-vector)
        tem2=find(abs(snew)>tol);% tem2 has the indices of nonzero elements in current snew
        % proceed if temp3=0 AND if all nonzero indices in tem2 do not fall into
        % both the BL and FL regions (otherwise go back to “for ii=1:size(s,2)”)
        if all(abs(temp3)<tol) && (isempty(find(all(tem2<=critval))) + isempty(find(all(tem2>critval)))==1)  
            if isempty(find(all(tem2<=critval))) && isempty(find(all(tem2>critval)))
            disp('contains both forward-looking and backward-looking, procedure fails');
            return
            elseif all(tem2<=critval)
                %mm=mm-1;
                r=r+1;
                SB(:,r)=snew;
            end     
        end
    end
    if ~isreal(SB)
       SB=[real(SB) imag(SB)];
    end   
    SB(mm+1:end,:)=[];
    SBO=orth(SB);
    SBN=null(SB');
    if ~isempty(SB) || size(SBO,2)<r
       SBNbar=[SBN' zeros(size(SBN,2),nm); zeros(nm, size(SBN,1)) eye(nm)];
       FAD=SBNbar*FAD*SBNbar';
       FB=SBNbar*FB;
       %Cbar=Cbar*SN;
       QQ_red=SBNbar*QQ_red*SBNbar';
       UU=SBNbar*UU;
       UU_redT=UU';
       ChatD=ChatD*SBNbar';
       Cbar=[ChatD;QQ_red;UU_redT];
       mm=size(SBN,2);
       r=0;
       SB=[];
    else 
        for nn=1:max_iter
            [s,mu] = eig(FAD);
            n_new=size(FAD,2);% n_new is the current size of A
            critval=n_new-nm; 
            for ii=1:size(s,2)
                temp3=Cbar*s(:,ii); % Cbar*each column of s (starting from the ii’th)
                snew=s(:,ii); % record the ii’th column of s (i.e. ii’th e-vector)
                tem2=find(abs(snew)>tol);% tem2 has the indices of nonzero elements in current snew
                if all(abs(temp3)<tol) && (isempty(find(all(tem2<=critval))) + isempty(find(all(tem2>critval)))==1)  
                  if isempty(find(all(tem2<=critval))) && isempty(find(all(tem2>critval)))
                  disp('contains both forward-looking and backward-looking, procedure fails');
                  return
                elseif all(tem2>critval)
                      %nm=nm-1;
                      p=p+1;
                     SF(:,p)=snew;
                  end     
                end
            end
            if ~isreal(SF)
               SF=[real(SF) imag(SF)];
            end
            SF=SF(critval+1:end,:);
            SFO=orth(SF);
            SFN=null(SF');
            if ~isempty(SF) || size(SFO,2)<p
               SFNbar=[eye(critval) zeros(critval,size(SFN,1)); zeros(size(SFN,2), critval)  SFN'];
               FAD=SFNbar*FAD*SFNbar';
               FB=SFNbar*FB;
               %Cbar=Cbar*SN;
               QQ_red=SFNbar*QQ_red*SFNbar';
               UU=SFNbar*UU;
               UU_redT=UU';
               ChatD=ChatD*SFNbar';
               Cbar=[ChatD;QQ_red;UU_redT];
               nm=nm-p;
               p=0;
               SF=[];
            else
               % record the size of the state
               ss=ss+1;
               nstate(ss)=size(FAD,1);
            break    
            end   
        end
    end
end                
                                
              
% find red-form matrices HAG,HB,CG (obsevable)
HAG=FAD;
HB=FB;
CG=Cbar(1:size(ChatD,1),:);
%UU_redT=UU';
mm=size(HAG,1)-nm;

% decompose HB to C and D: this may be a problem - check again
i_d=find(all(ismember(BB,dd)));
i_c=find(all(ismember(BB,cc)));
i_d=[1:size(dd,2)];
i_c=[size(dd,2)+1:size(HB,2)];
D_red=HB(:,i_d); 
C_red=HB(:,i_c);

% rescale w_red if terms of instruments are zeros - applicable for opt 
scale=1e-4; 
%if isempty(find(RR)) %%%%% fix this !!!!!!!!!!!!!!!!!!
for i=1:num_inst
if abs(RR(i,i))<tol
RR(i,i)=RR(i,i)+scale;      
end
end
% resemble W - [QQ, UU; UU',RR]
W_red=[QQ_red, UU_redT'; UU_redT, RR;];
%W_red=zeros(size(QQ_red,1)+size(UU_redT,1));
%W_red(i_varW,i_varW)=RR;
%W_red([1:size(QQ_red,2)],[1:size(QQ_red,2)])=QQ_red;
%W_red([1:size(e2,2)],[i_varW:end]);
%W_red=[QQ_red]

