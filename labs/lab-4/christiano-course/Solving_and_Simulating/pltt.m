function pltt(gkk,gsigsig,gkk2,gsigsig2,gee,gee2)
kk=-.2:.01:.2;
z=0.5*gkk*kk.^2+.5*gsigsig;
z2=0.5*gkk2*kk.^2+.5*gsigsig2;
subplot(121)
plot(100*kk,100*z,100*kk,100*z2,'*-')
legend('\gamma = 2','\gamma = 20')
ylabel('100*( k_{t+1} (2^{nd} order) - k_t+1 (1^{st} order) )','FontSize',18)
xlabel('100*( k_t - k^* ), percent deviation of initial capital from steady state','FontSize',18)
axis([100*kk(1) 100*kk(end) 0 max(100*z2)*1.2])

z2=[];
z1=[];
aa=-.2:.01:.2;
z=0.5*gee*aa.^2+.5*gsigsig;
z2=0.5*gee2*aa.^2+.5*gsigsig2;
subplot(122)
plot(100*aa,100*z,100*aa,100*z2,'*-')
xlabel('100*a_t, percent deviation of initial shock from steady state','FontSize',18)
ylabel('100*( k_{t+1} (2^{nd} order) - k_t+1 (1^{st} order) )','FontSize',18)
axis([100*aa(1) 100*aa(end) 0 max(100*z2)])
