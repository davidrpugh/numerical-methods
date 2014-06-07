%%This code produces a set of graphs that plots the IRFs produced by the linear vs
%%the non linear version of the dynare code
%%(c) CIMS Univeristy of Surrey
%%The Science and  Art of DSGE Modelling: Construction, Calibration, Estimation and Policy
 

clear;
close all;
clc;

%Variables names to plot
names = char('Output', 'Consumption', 'Investment', 'Hours worked', 'Real wage', 'Real interest rate', 'Nominal interest rate', 'Tobin Q', 'Inflation');


%load non linear simulations
%load NK_free_parameters_RES_Course_results
%load NK_Rotemberg_free_parameters_RES_Course_results
load NK_RES_Course_2_results

%Rename the IRFs for each variable of interest in the non linear model
Y_epsA     = oo_.irfs.YY_epsA;
C_epsA     = oo_.irfs.CC_epsA;
I_epsA     = oo_.irfs.II_epsA;
h_epsA     = oo_.irfs.hh_epsA;
WP_epsA    = oo_.irfs.WPWP_epsA;
R_epsA     = oo_.irfs.RR_epsA;
Rn_epsA     = oo_.irfs.RnRn_epsA;
Q_epsA     = oo_.irfs.QQ_epsA;
PIE_epsA     = oo_.irfs.PIEPIE_epsA;

Y_epsG     = oo_.irfs.YY_epsG;
C_epsG     = oo_.irfs.CC_epsG;
I_epsG     = oo_.irfs.II_epsG;
h_epsG     = oo_.irfs.hh_epsG;
WP_epsG    = oo_.irfs.WPWP_epsG;
R_epsG     = oo_.irfs.RR_epsG;
Rn_epsG     = oo_.irfs.RnRn_epsG;
Q_epsG     = oo_.irfs.QQ_epsG;
PIE_epsG     = oo_.irfs.PIEPIE_epsG;

Y_epsMS     = oo_.irfs.YY_epsMS;
C_epsMS     = oo_.irfs.CC_epsMS;
I_epsMS     = oo_.irfs.II_epsMS;
h_epsMS     = oo_.irfs.hh_epsMS;
WP_epsMS    = oo_.irfs.WPWP_epsMS;
R_epsMS     = oo_.irfs.RR_epsMS;
Rn_epsMS     = oo_.irfs.RnRn_epsMS;
Q_epsMS     = oo_.irfs.QQ_epsMS;
PIE_epsMS     = oo_.irfs.PIEPIE_epsMS;

Y_epsMPS     = oo_.irfs.YY_epsMPS;
C_epsMPS     = oo_.irfs.CC_epsMPS;
I_epsMPS     = oo_.irfs.II_epsMPS;
h_epsMPS     = oo_.irfs.hh_epsMPS;
WP_epsMPS    = oo_.irfs.WPWP_epsMPS;
R_epsMPS     = oo_.irfs.RR_epsMPS;
Rn_epsMPS     = oo_.irfs.RnRn_epsMPS;
Q_epsMPS     = oo_.irfs.QQ_epsMPS;
PIE_epsMPS     = oo_.irfs.PIEPIE_epsMPS;

Y_epsAtrend     = oo_.irfs.YY_epsAtrend;
C_epsAtrend     = oo_.irfs.CC_epsAtrend;
I_epsAtrend     = oo_.irfs.II_epsAtrend;
h_epsAtrend     = oo_.irfs.hh_epsAtrend;
WP_epsAtrend    = oo_.irfs.WPWP_epsAtrend;
R_epsAtrend     = oo_.irfs.RR_epsAtrend;
Rn_epsAtrend     = oo_.irfs.RnRn_epsAtrend;
Q_epsAtrend     = oo_.irfs.QQ_epsAtrend;
PIE_epsAtrend     = oo_.irfs.PIEPIE_epsAtrend;

%Create a vector of IRFs per each shock
irfs_NK1_epsA =[Y_epsA;  C_epsA;  I_epsA; h_epsA;  WP_epsA;  R_epsA; Rn_epsA; Q_epsA; PIE_epsA];
irfs_NK1_epsG =[Y_epsG;  C_epsG;  I_epsG; h_epsG;  WP_epsG;  R_epsG; Rn_epsG; Q_epsG; PIE_epsG];
irfs_NK1_epsMS =[Y_epsMS;  C_epsMS;  I_epsMS; h_epsMS;  WP_epsMS;  R_epsMS; Rn_epsMS; Q_epsMS; PIE_epsMS];
irfs_NK1_epsMPS =[Y_epsMPS;  C_epsMPS;  I_epsMPS; h_epsMPS;  WP_epsMPS;  R_epsMPS; Rn_epsMPS; Q_epsMPS; PIE_epsMPS];
irfs_NK1_epsAtrend =[Y_epsAtrend;  C_epsAtrend;  I_epsAtrend; h_epsAtrend;  WP_epsAtrend;  R_epsAtrend; Rn_epsAtrend; Q_epsAtrend; PIE_epsAtrend];

%load linear simulations

%load non linear simulations
%load NK_Rotemberg_RES_Course_2_results
%load NK_RES_Course_results
load NK_Rotemberg_RES_Course_2_results
%load NK_free_parameters_RES_Course_results
%load NK_Rotemberg_free_parameters_RES_Course_results
%Rename the IRFs for each variable of interest in the non linear model
Y_epsA     = oo_.irfs.YY_epsA;
C_epsA     = oo_.irfs.CC_epsA;
I_epsA     = oo_.irfs.II_epsA;
h_epsA     = oo_.irfs.hh_epsA;
WP_epsA    = oo_.irfs.WPWP_epsA;
R_epsA     = oo_.irfs.RR_epsA;
Rn_epsA     = oo_.irfs.RnRn_epsA;
Q_epsA     = oo_.irfs.QQ_epsA;
PIE_epsA     = oo_.irfs.PIEPIE_epsA;

Y_epsG     = oo_.irfs.YY_epsG;
C_epsG     = oo_.irfs.CC_epsG;
I_epsG     = oo_.irfs.II_epsG;
h_epsG     = oo_.irfs.hh_epsG;
WP_epsG    = oo_.irfs.WPWP_epsG;
R_epsG     = oo_.irfs.RR_epsG;
Rn_epsG     = oo_.irfs.RnRn_epsG;
Q_epsG     = oo_.irfs.QQ_epsG;
PIE_epsG     = oo_.irfs.PIEPIE_epsG;

Y_epsMS     = oo_.irfs.YY_epsMS;
C_epsMS     = oo_.irfs.CC_epsMS;
I_epsMS     = oo_.irfs.II_epsMS;
h_epsMS     = oo_.irfs.hh_epsMS;
WP_epsMS    = oo_.irfs.WPWP_epsMS;
R_epsMS     = oo_.irfs.RR_epsMS;
Rn_epsMS     = oo_.irfs.RnRn_epsMS;
Q_epsMS     = oo_.irfs.QQ_epsMS;
PIE_epsMS     = oo_.irfs.PIEPIE_epsMS;

Y_epsMPS     = oo_.irfs.YY_epsMPS;
C_epsMPS     = oo_.irfs.CC_epsMPS;
I_epsMPS     = oo_.irfs.II_epsMPS;
h_epsMPS     = oo_.irfs.hh_epsMPS;
WP_epsMPS    = oo_.irfs.WPWP_epsMPS;
R_epsMPS     = oo_.irfs.RR_epsMPS;
Rn_epsMPS     = oo_.irfs.RnRn_epsMPS;
Q_epsMPS     = oo_.irfs.QQ_epsMPS;
PIE_epsMPS     = oo_.irfs.PIEPIE_epsMPS;

Y_epsAtrend     = oo_.irfs.YY_epsAtrend;
C_epsAtrend     = oo_.irfs.CC_epsAtrend;
I_epsAtrend     = oo_.irfs.II_epsAtrend;
h_epsAtrend     = oo_.irfs.hh_epsAtrend;
WP_epsAtrend    = oo_.irfs.WPWP_epsAtrend;
R_epsAtrend     = oo_.irfs.RR_epsAtrend;
Rn_epsAtrend     = oo_.irfs.RnRn_epsAtrend;
Q_epsAtrend     = oo_.irfs.QQ_epsAtrend;
PIE_epsAtrend     = oo_.irfs.PIEPIE_epsAtrend;

%Create a vector of IRFs per each shock
irfs_NK2_epsA =[Y_epsA;  C_epsA;  I_epsA; h_epsA;  WP_epsA;  R_epsA; Rn_epsA; Q_epsA; PIE_epsA];
irfs_NK2_epsG =[Y_epsG;  C_epsG;  I_epsG; h_epsG;  WP_epsG;  R_epsG; Rn_epsG; Q_epsG; PIE_epsG];
irfs_NK2_epsMS =[Y_epsMS;  C_epsMS;  I_epsMS; h_epsMS;  WP_epsMS;  R_epsMS; Rn_epsMS; Q_epsMS; PIE_epsMS];
irfs_NK2_epsMPS =[Y_epsMPS;  C_epsMPS;  I_epsMPS; h_epsMPS;  WP_epsMPS;  R_epsMPS; Rn_epsMPS; Q_epsMPS; PIE_epsMPS];
irfs_NK2_epsAtrend =[Y_epsAtrend;  C_epsAtrend;  I_epsAtrend; h_epsAtrend;  WP_epsAtrend;  R_epsAtrend; Rn_epsAtrend; Q_epsAtrend; PIE_epsAtrend];

%load linear simulations




%Options for the plot
h=figure('Position', [600, 0, 1000, 900]);
axes ('position', [0, 0, 1, 1]);

%Figure for the technology shock
%Loop over the number of endogenous variables to plot
F1=figure(1)
set(F1, 'numbertitle','off')
set(F1, 'name', 'Impulse response functions (Technology Shock)')
for j = 1:9;
    subplot(3,3,j), plot(irfs_NK1_epsA (j,:),'-k', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','k'); hold on;
                    plot(irfs_NK2_epsA(j,:),'--r', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','r'); hold on;
                    xlabel('Quarters');
                    ylabel('% dev from SS');
                    grid on
                    title(names(j,:),'FontSize',10)
axis tight;                   
end;

legend('NK (Calvo)', 'NK (Rotemberg)', 'Best', 'Orientation', 'horizontal')

%Options for the plot
h=figure('Position', [600, 0, 1000, 900]);
axes ('position', [0, 0, 1, 1]);

%Figure for the government spending shock
%Loop over the number of endogenous variables to plot
F1=figure(2)
set(F1, 'numbertitle','off')
set(F1, 'name', 'Impulse response functions (Gov Spending Shock)')
for j = 1:9;
    subplot(3,3,j), plot(irfs_NK1_epsG (j,:),'-k', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','k'); hold on;
                    plot(irfs_NK2_epsG(j,:),'--r', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','r'); hold on;
                    
                    xlabel('Quarters');
                    ylabel('% dev from SS');
                    grid on
                    title(names(j,:),'FontSize',10)
axis tight;                   
end;

legend('NK (Calvo)', 'NK (Rotemberg)' , 'Best', 'Orientation', 'horizontal')


%Options for the plot
h=figure('Position', [600, 0, 1000, 900]);
axes ('position', [0, 0, 1, 1]);

%Figure for the Monetary policy shock
%Loop over the number of endogenous variables to plot
F1=figure(3)
set(F1, 'numbertitle','off')
set(F1, 'name', 'Impulse response functions (Mark-up Shock)')
for j = 1:9;
    subplot(3,3,j), plot(irfs_NK1_epsMS (j,:),'-k', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','k'); hold on;
                    plot(irfs_NK2_epsMS(j,:),'--r', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','r'); hold on;
                    xlabel('Quarters');
                    ylabel('% dev from SS');
                    grid on
                    title(names(j,:),'FontSize',10)
axis tight;                   
end;

legend('NK (Calvo)', 'NK (Rotemberg)', 'Location', 'Best', 'Orientation', 'horizontal')

F1=figure(4)
set(F1, 'numbertitle','off')
set(F1, 'name', 'Impulse response functions (Monetary Policy Shock)')
for j = 1:9;
    subplot(3,3,j), plot(irfs_NK1_epsMPS (j,:),'-k', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','k'); hold on;
                    plot(irfs_NK2_epsMPS(j,:),'--r', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','r'); hold on;
                    
                   xlabel('Quarters');
                    ylabel('% dev from SS');
                    grid on
                    title(names(j,:),'FontSize',10)
axis tight;                   
end;

legend('NK (Calvo)', 'NK (Rotemberg)', 'Location', 'Best', 'Orientation', 'horizontal')

F1=figure(5)
set(F1, 'numbertitle','off')
set(F1, 'name', 'Impulse response functions (Trend Shock)')
for j = 1:9;
    subplot(3,3,j), plot(irfs_NK1_epsAtrend (j,:),'-k', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','k'); hold on;
                    plot(irfs_NK2_epsAtrend(j,:),'--r', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','r'); hold on;
                    
                   xlabel('Quarters');
                    ylabel('% dev from SS');
                    grid on
                    title(names(j,:),'FontSize',10)
axis tight;                   
end;

legend('NK (Calvo)', 'NK (Rotemberg)', 'Location', 'Best', 'Orientation', 'horizontal')

