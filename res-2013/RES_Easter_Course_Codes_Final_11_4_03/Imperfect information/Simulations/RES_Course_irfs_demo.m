clear; close all;


load RES_Course_demo_pi_results;
oo1=oo_;
load RES_Course_demo_ii_results;
oo2=oo_;
load RES_Course_demo_ii2_results;
oo3=oo_;

F1=figure(1);
set(F1, 'numbertitle','off')
set(F1, 'name', 'Impulse response functions')
plot(oo1.irfs.pi_eps,'r','linewidth',1.5) 
hold on
plot(oo2.irfs.pi_eps, '--b','linewidth',1.5);
plot(oo3.irfs.pi_eps, ':g','linewidth',1.5);
hold off
xlabel('TIME (QUARTERS)')
ylabel('INFLATION  \pi_t')

legend('PI (\sigma^2_w=0)','II (\sigma^2_w=1)','II (\sigma^2_w=2)','Location','Best');