clear; close all;

load GK_RES_Course_est_PI_results_150000_draws;
oo1=oo_;
load BGG_RES_Course_est_PI_results_150000_draws;
oo2=oo_;
load KM_RES_Course_est_PI_results_150000_draws;
oo3=oo_;

lag = (1:1:10);

%% A
F1=figure(1);
set(F1, 'numbertitle','off')
set(F1, 'name', 'Impulse response functions (TFP)')
  h1 = area(1:10);
  set(h1,'FaceColor',[.9 .9 .9]);
subplot(2,2,1)
plot(lag,oo1.irfs.pinfobs_epsA(:,[1: 10])','b','linewidth',2);
hold on
plot(lag,oo2.irfs.pinfobs_epsA(:,[1: 10]),'--r','linewidth',2); 
plot(lag,oo3.irfs.pinfobs_epsA(:,[1: 10]), '-.g','linewidth',2);
hold off
title('Inflation Rate')
%ylabel('Output')

subplot(2,2,2)
plot(lag,oo1.irfs.dy_epsA(:,[1: 10])','b','linewidth',2);
hold on
plot(lag,oo2.irfs.dy_epsA(:,[1: 10]),'--r','linewidth',2); 
plot(lag,oo3.irfs.dy_epsA(:,[1: 10]), '-.g','linewidth',2);
hold off
title('Output')

subplot(2,2,3)
plot(lag,oo1.irfs.robs_epsA(:,[1: 10])','b','linewidth',2);
hold on
plot(lag,oo2.irfs.robs_epsA(:,[1: 10]),'--r','linewidth',2); 
plot(lag,oo3.irfs.robs_epsA(:,[1: 10]), '-.g','linewidth',2);
hold off
title('Nominal Interest Rate')

subplot(2,2,4)
plot(lag,oo1.irfs.rkn_obs_epsA(:,[1: 10])','b','linewidth',2);
hold on
plot(lag,oo2.irfs.rkn_obs_epsA(:,[1: 10]),'--r','linewidth',2); 
plot(lag,oo3.irfs.rkn_obs_epsA(:,[1: 10]), '-.g','linewidth',2);
hold off
title('BAA')

legend('Model GK','Model BGG','Model KM','Location','Best','Orientation','vertical');
