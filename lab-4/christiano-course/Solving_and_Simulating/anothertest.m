% This studies the simple nonlinear expression,
% y(t) = rho*y(t-1) + alph*y(t-1)^2 + epsilon(t)
%this is a stylized representation of the 2nd order approximation
%to the neoclassical model.
clear all
close all
randn('seed',1);
T=500;
sig=0.1;
epsilon=randn(T,1)*sig;
y(1)=epsilon(1);
y1(1)=epsilon(1);
y2(1)=epsilon(1);
rho=.8;
alph=.5;
yy=-.6:.0001:.9;
for ii = 1:length(yy)
    yp(ii)=rho*yy(ii)+alph*yy(ii)^2;
    yph(ii)=rho*yy(ii)+alph*yy(ii)^2+2*sig;
    ypl(ii)=rho*yy(ii)+alph*yy(ii)^2-2*sig;
end
figure
plot(yy,yp,yy,yy,yy,yph,yy,ypl)
axis tight

figure
for ii = 2:T
    y(ii)=rho*y(ii-1)+alph*y(ii-1)^2+epsilon(ii);
    y1(ii)=rho*y1(ii-1)+epsilon(ii);
    y2(ii)=rho*y2(ii-1)+alph*y1(ii-1)^2+epsilon(ii);
end

tt=1:T;
plot(tt,y(tt),tt,y1(tt),'*-')
legend('second order difference equation','simulations ignoring quadratic term, (\alpha = 0)')
%axis([tt(1) tt(end) -1 2]);
axis tight

figure

tt=1:T;
plot(tt,y2(tt),tt,y1(tt),'*-')
legend('pruned solution to second order difference equation','simulations ignoring quadratic term, (\alpha = 0)')
%axis([tt(1) tt(end) -1 2]);
axis tight

break

subplot(211)
plot(tt,y(tt),tt,y2(tt),'*')
legend('actual','pruned')
axis([tt(1) tt(end) -1 2]);
subplot(212)
plot(tt,y1(tt),tt,y2(tt),'*')
legend('linear','pruned')
axis([tt(1) tt(end) -1 1]);

