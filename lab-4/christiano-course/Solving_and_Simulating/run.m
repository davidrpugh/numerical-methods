%this computes and studies the 2nd order approximation to the policy rule
%for capital in the neoclassical model.

close all
clear all

figure
k=0:.001:1.5;
KK=1.5:.001:2.5;
kk=[k KK];
%y=a+b*x+c*x^2
%KK(1)^.36=a+b*KK(1)+c*KK(1)^2
%derivative: b+2*c*x
%.36*KK(1)^(.36-1)
%.36*KK(1)^(.36-1)=b+2*c*KK(1)

c=2;
%fix c and then 
b=0.36*KK(1)^(0.36-1)-2*c*KK(1);
a = KK(1)^.36-b*KK(1)-c*KK(1)^2;
tt=1:length(k);
plot(k,k.^.36,kk,kk,KK,a+b*KK+c*KK.^2)
axis tight
xlabel('\it K_{t}','Fontsize',28)
ylabel('\it K_{t+1}','Fontsize',28)


alpha=0.36;
beta=0.99;
delta=.02;
rho=0.95;
Ve=.01^2;
lg=0;

gamma=2;
[gklin, galin, gelin, gsiglin, gkklin, gaalin, geelin, gsigsiglin, gkalin, gkelin, gksiglin, gaelin, ...
        gasiglin, gesiglin, ksslin] = second(alpha,beta,delta,gamma,rho,Ve,lg);
    
lg=1;
gamma=2;
[gk, ga, ge, gsig, gkk, gaa, gee, gsigsig, gka, gke, gksig, gae, ...
        gasig, gesig, kss] = second(alpha,beta,delta,gamma,rho,Ve,lg);
    
yy=[];
yp=[];
yplin=[];
yy=-500:.1:5;
for ii = 1:length(yy)
    yplin(ii) = gklin*yy(ii) + gsiglin + (1/2)*gkklin*(yy(ii))^2 + gsigsiglin ...
        + gksiglin*yy(ii) ;
end
figure    
plot(yy,yplin-yy,yy,yy-yy)
axis tight

figure
    
yy=[];
yp=[];
yplin=[];
yy=-7:.001:.1;
yy=-5:.001:.1;
for ii = 1:length(yy)
    yp(ii) = gk*yy(ii) + gsig + (1/2)*(gkk*yy(ii)^2 + gsigsig) + gksig*yy(ii) ;
    yplin(ii) = gklin*(exp(yy(ii)+kss)-ksslin) + gsiglin + (1/2)*gkklin*(exp(yy(ii)+kss)-ksslin)^2 + gsigsiglin ...
        + gksiglin*(exp(yy(ii)+kss)-ksslin) ;
end


figure
plot(exp(yy+kss),exp(yp+kss),exp(yy+kss),yplin+ksslin,exp(yy+kss),exp(yy+kss),'--')
legend('log linear','linear','45 degree')
title('Comparing Two Second Order Approximations','Fontsize',28)
ylabel('actual capital stock, K'', next period','Fontsize',28)
xlabel('actual K','Fontsize',28)
axis tight

yy=[];
yp=[];
yplin=[];
yy=2:.001:3.01;
for ii = 1:length(yy)
    yp(ii) = gk*yy(ii) + gsig + (1/2)*(gkk*yy(ii)^2 + gsigsig) + gksig*yy(ii) ;
    yplin(ii) = gklin*(exp(yy(ii)+kss)-ksslin) + gsiglin + (1/2)*gkklin*(exp(yy(ii)+kss)-ksslin)^2 + gsigsiglin ...
        + gksiglin*(exp(yy(ii)+kss)-ksslin) ;
end

figure
plot(exp(yy+kss),exp(yp+kss),exp(yy+kss),yplin+ksslin,exp(yy+kss),exp(yy+kss),'--')
legend('log linear','linear','45 degree')
title('Comparing Two Second Order Approximations','Fontsize',28)
ylabel('actual capital stock, K'', next period','Fontsize',28)
xlabel('actual K','Fontsize',28)
axis tight

yy=[];
yp=[];
yplin=[];
yy=-.1:.001:3;
a=2*sqrt(Ve/(1-rho^2));
aa=4*sqrt(Ve/(1-rho^2));
for ii = 1:length(yy)
    ypa(ii) = gk*yy(ii) + ge*a + (1/2)*(gkk*yy(ii)^2 + gee*a^2 + gsigsig) + gksig*yy(ii) ;
    ypaa(ii) = gk*yy(ii) + ge*aa + (1/2)*(gkk*yy(ii)^2 + gee*aa^2 + gsigsig) + gksig*yy(ii) ;
    yp(ii) = gk*yy(ii) + (1/2)*(gkk*yy(ii)^2 + gsigsig) + gksig*yy(ii) ;
    yplin(ii) = gklin*(exp(yy(ii)+kss)-ksslin) + gsiglin + (1/2)*(gkklin*(exp(yy(ii)+kss)-ksslin)^2 + gsigsiglin) ...
        + gksiglin*(exp(yy(ii)+kss)-ksslin) ;
end

[Ylin,Ilin]=min(abs(yplin-(exp(yy+kss)-ksslin)));
[Y,II]=min(abs((yp-yy)./yy));
ix=0;
ixa=0;
ixaa=0;
IF=[];
IFa=[];
IFaa=[];
for ii = 2:length(yp)
    if (yp(ii)-yy(ii))*(yp(ii-1)-yy(ii-1))<0
        ix=ix+1;
        IF(ix)=ii;
    end
    if (ypa(ii)-yy(ii))*(ypa(ii-1)-yy(ii-1))<0
        ixa=ixa+1;
        IFa(ixa)=ii;
    end
    if (ypaa(ii)-yy(ii))*(ypaa(ii-1)-yy(ii-1))<0
        ixaa=ixaa+1;
        IFaa(ixaa)=ii;
    end
end

if ix > 2
    error('fatal (run): something very wrong')
end

figure
tt=1:length(yp);
plot(yy,ypa-yy,yy,ypaa-yy,yy,yp-yy,yy,zeros(length(yy),1),yy(IF),yp(IF)-yy(IF),'*',yy(IFa),ypa(IFa)-yy(IFa),'*',yy(IFaa),ypaa(IFaa)-yy(IFaa),'*')
axis tight
title('Second order approximation of policy rule for log-linear capital','Fontsize',28)
ylabel('\it k_{t+1} - k_{t}','Fontsize',28)
xlabel('\it k_{t} - k*','Fontsize',28)

figure
ts=1:20;
numsim=30;

e=zeros(length(ts),1);
ys=zeros(length(ts),numsim);
yf=zeros(length(ts),numsim);
a=zeros(length(ts),numsim);
ct=6;
plt=0;
for jj = 1:numsim
    a(1,jj)=ct*sqrt(Ve);
    ys(1,jj) = ge*a(1,jj) + (1/2)*(gee*a(1,jj)^2 + gsigsig) ;
    yf(1,jj) = ge*a(1,jj)  ;
    for ii = 2:length(ts)
        a(ii,jj)=rho*a(ii-1,jj)+ct*randn*sqrt(Ve);
        ys(ii,jj) = gk*ys(ii-1,jj) + ge*a(ii,jj) + (1/2)*(gkk*ys(ii-1,jj)^2 + gee*a(ii,jj)^2 + gsigsig) + gksig*ys(ii-1,jj) ;
        yf(ii,jj) = gk*yf(ii-1,jj) + ge*a(ii,jj)  ;
    end
    if plt == 1
    plot(ts,ys(:,jj),ts,yf(:,jj),'*')
    legend('second order','first order')
    axis tight
    pause(1)
    end
end

mys=mean(ys');
myf=mean(yf');

plot(exp(yy+kss),exp(yp+kss),exp(yy+kss),yplin+ksslin,exp(yy+kss),exp(yy+kss),'--')
legend('log linear','linear','45 degree')
title('Comparing Two Second Order Approximations','Fontsize',28)
ylabel('actual capital stock, K'', next period','Fontsize',28)
xlabel('actual K','Fontsize',28)
axis tight

   
gamma=20;
[gk2, ga2, ge2, gsig2, gkk2, gaa2, gee2, gsigsig2, gka2, gke2, gksig2, gae2, ...
        gasig2, gesig2, kss2] = second(alpha,beta,delta,gamma,rho,Ve,lg);
   
fprintf('gk = %6.2f, ga = %6.2f, ge = %6.2f, gsig = %6.2f\n',gk,ga,ge,gsig)
fprintf('gkk = %6.2f, gaa = %6.2f, gee = %6.2f, gsigsig = %6.2f, gka = %6.2f\n',gkk,gaa,gee,gsigsig,gka)
fprintf('gke = %6.2f, gksig = %6.2f, gae = %6.2f, gasig = %6.2f, gesig = %6.2f, kss = %6.2f\n',gke,gksig ...
    ,gae,gasig,gesig,kss)

%The solution is of this form:
%
% g = kss + gk*(k - kss) + ga*a_1 + ge*epsilon + gsig*sig ...
%     + (1/2)*(gkk*(k - kss)^2 + gaa*a_1^2 + gee*epsilon^2 + gsigsig*sig^2) ...
%     + gka*(k - kss)*a_1 + gke*(k - kss)*epsilon + gksig*(k - kss)*sig ...
%     + gae*a_1*epsilon + gasig*a_1*sig + gesig*epsilon*sig ;
%
% Under this arrangement, the state includes epsilon and a_1. 
% However, intuition suggests that only current a needs to be in the state.
% We now show how to confirm this.
%
% Here, a=rho*a_1 + epsilon.
%
% if ga=rho*ge, then ga*a_1 + ge*epsilon = ge*(rho*a_1 + epsilon)=ge*a,
% so we check this condition across ga, ge, rho.
%
% Consider:
%
%    gee*a^2 = gee*(rho*a_1+epsilon)^2 = gee*(rho^2*a_1^2 + epsilon^2 +
%       2*rho*a_1*epsilon) 
%        = gee*rho^2*a_1^2 + gee*epsilon^2 + 2*gee*rho*a_1*epsilon.
% so,
%     (1/2)*gee*a^2 = (1/2)*gee*rho^2*a_1^2 + (1/2)*gee*epsilon^2 + gee*rho*a_1*epsilon.
%
% Thus, we expect gee*rho^2=gaa, gae=gee*rho
% Also, suppose gka=rho*gke. Then,
%  gka*(k - kss)*a_1 + gke*(k - kss)*epsilon = rho*gke*(k - kss)*a_1 + gke*(k - kss)*epsilon = gke*((k - kss)*rho*a_1 + (k - kss)*epsilon)
%  = gke*(k - kss)*a
% Finally, suppose gasig=rho*gesig. Then,
%     gasig*a_1*sig + gesig*epsilon*sig = gesig*(rho*a_1*sig + epsilon*sig)=gesig*a*sig


if abs(ga/ge-rho) > .1e-10 || abs(gee*rho^2-gaa) > .1e-10 || abs(gae-gee*rho) > .1e-10 || abs(gka-rho*gke) > .1e-10 || abs(gasig-rho*gesig) > .1e-10
    error('fatal (run) policy rule parameters do not satisfy expected conditions')
end

% thus, we can write:
%
% g = kss + gk*(k - kss) + ge*a + gsig*sig ...
%     + (1/2)*(gkk*(k - kss)^2 + gee*a^2 + gsigsig*sig^2) ...
%     + gke*(k - kss)*a + gksig*(k - kss)*sig + gesig*a*sig
%
% or, setting y=k-kss, yp = y in next period:
%
% yp = gk*y + ge*a + gsig*sig + (1/2)*(gkk*y^2 + gee*a^2 + gsigsig*sig^2) + gke*y*a + gksig*y*sig + gesig*a*sig

figure
yy=[];
yp=[];
yy=-.00007:.000001:.001;
for ii = 1:length(yy)
    yp(ii) = gk*yy(ii) + gsig + (1/2)*(gkk*yy(ii)^2 + gsigsig) + gksig*yy(ii) ;
end

[Y,I]=min(abs(yy));
[YY,II]=min(abs(yy-yp));

tt=[];
tt=I-5:I+5;
tt=1:length(yy);
plot(exp(yy(tt)+kss),exp(yp(tt)+kss),exp(yy(tt)+kss),exp(yy(tt)+kss),'--',exp(yy(I)+kss),exp(yp(I)+kss),'*',exp(yy(II)+kss),exp(yp(II)+kss),'*')
legend('policy rule','45 degree line')
title(' Policy Rule in Neighborhood of Steady State','Fontsize',28)
ylabel('{exp(k'')} ','Fontsize',28)
xlabel('{exp(k)}','Fontsize',28)
axis tight

figure
yy=[];
yp=[];
yy=-.00007:.00001:.001;
yy=-10:.01:2.8;
for ii = 1:length(yy)
    yp(ii) = gk*yy(ii) + gsig + (1/2)*(gkk*yy(ii)^2 + gsigsig) + gksig*yy(ii) ;
end

[Y,I]=min(abs(yy));
[YY,II]=min(abs(yy-yp));
tt=[];
tt=I-5:I+5;
tt=1:length(yy);
%plot(exp(yy(tt)+kss),exp(yp(tt)+kss),exp(yy(I)+kss),exp(yp(I)+kss),'*',exp(yy(tt)+kss),exp(yy(tt)+kss))
plot(exp(yy(tt)+kss),exp(yp(tt)+kss)-exp(yy(tt)+kss),exp(yy(I)+kss),exp(yp(I)+kss)-exp(yy(I)+kss),'*',exp(yy(tt)+kss),zeros(length(tt),1))
ylabel('exp(k'') - exp(k)','Fontsize',28)
xlabel('exp(k)','Fontsize',28)
title(' exp(k'') - exp(k) against exp(k)','Fontsize',28)
axis tight

figure
subplot(211)
yy=[];
yp=[];
yy=-.00007:.00001:.001;
yy=-5:.001:3;
yy=-.5:.001:3;
yy=-48:.001:1;
for ii = 1:length(yy)
    yp(ii) = gk*yy(ii) + gsig + (1/2)*(gkk*yy(ii)^2 + gsigsig) + gksig*yy(ii) ;
end
plot(yy+kss,yp+kss,yy+kss,yy+kss)
ylabel('\it k_{t+1} - \it k_{t}','Fontsize',28)
xlabel('\it k_{t}','Fontsize',28)
title('Policy rule, holding {\it a_t = 0}, for capital stock, exp({\it k_t}), and log, capital stock, {\it k_t}','Fontsize',28)
axis tight
subplot(212)
plot(exp(yy+kss),exp(yp)-exp(yy),exp(yy+kss),yy-yy)
ylabel('exp(\it k_{t+1})-exp(\it k_{t})','Fontsize',28)
xlabel('exp(\it k_{t})','Fontsize',28)
axis tight

pltt(gkk,gsigsig,gkk2,gsigsig2,gee,gee2);

Ve=.10^2;

[gk, ga, ge, gsig, gkk, gaa, gee, gsigsig, gka, gke, gksig, gae, ...
        gasig, gesig, kss] = second(alpha,beta,delta,gamma,rho,Ve,lg);

fprintf('gk = %6.2f, ga = %6.2f, ge = %6.2f, gsig = %6.2f\n',gk,ga,ge,gsig)
fprintf('gkk = %6.2f, gaa = %6.2f, gee = %6.2f, gsigsig = %6.2f, gka = %6.2f\n',gkk,gaa,gee,gsigsig,gka)
fprintf('gke = %6.2f, gksig = %6.2f, gae = %6.2f, gasig = %6.2f, gesig = %6.2f, kss = %6.2f\n',gke,gksig ...
    ,gae,gasig,gesig,kss)



