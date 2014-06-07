%%This code produces a set of graphs that plots the IRFs produced by the linear vs
%%the non linear version of the dynare code
%%(c) CIMS Univeristy of Surrey
%%The Science and  Art of DSGE Modelling: Construction, Calibration, Estimation and Policy
 

clear;
close all;
clc;

%Variables names to plot
names = char('Output', 'Consumption', 'Investment', 'Hours worked', 'Real wage', 'Real interest rate');


%load model 1 (can be any model)
load RBC_RES_Course_results

%Rename the IRFs for each variable of interest in the non linear model
Y_epsA     = oo_.irfs.YY_epsA;
C_epsA     = oo_.irfs.CC_epsA;
I_epsA     = oo_.irfs.II_epsA;
h_epsA     = oo_.irfs.hh_epsA;
WP_epsA    = oo_.irfs.WPWP_epsA;
R_epsA     = oo_.irfs.RR_epsA;


Y_epsG     = oo_.irfs.YY_epsG;
C_epsG     = oo_.irfs.CC_epsG;
I_epsG     = oo_.irfs.II_epsG;
h_epsG     = oo_.irfs.hh_epsG;
WP_epsG    = oo_.irfs.WPWP_epsG;
R_epsG     = oo_.irfs.RR_epsG;

%Create a vector of IRFs per each shock
irfs_RBC_1_epsA =[Y_epsA;  C_epsA;  I_epsA; h_epsA;  WP_epsA;  R_epsA];
irfs_RBC_1_epsG =[Y_epsG;  C_epsG;  I_epsG; h_epsG;  WP_epsG;  R_epsG];

%load model 2 (can be any other model)
load RBC_Inv_Costs_RES_Course_results

%Rename the IRFs for each variable of interest in the non linear model
Y_epsA     = oo_.irfs.YY_epsA;
C_epsA     = oo_.irfs.CC_epsA;
I_epsA     = oo_.irfs.II_epsA;
h_epsA     = oo_.irfs.hh_epsA;
WP_epsA    = oo_.irfs.WPWP_epsA;
R_epsA     = oo_.irfs.RR_epsA;


Y_epsG     = oo_.irfs.YY_epsG;
C_epsG     = oo_.irfs.CC_epsG;
I_epsG     = oo_.irfs.II_epsG;
h_epsG     = oo_.irfs.hh_epsG;
WP_epsG    = oo_.irfs.WPWP_epsG;
R_epsG     = oo_.irfs.RR_epsG;

%Create a vector of IRFs per each shock
irfs_RBC_2_epsA =[Y_epsA;  C_epsA;  I_epsA; h_epsA;  WP_epsA;  R_epsA];
irfs_RBC_2_epsG =[Y_epsG;  C_epsG;  I_epsG; h_epsG;  WP_epsG;  R_epsG];




%Options for the plot
h=figure('Position', [600, 0, 1000, 900]);
axes ('position', [0, 0, 1, 1]);

%Figure for the technology shock
%Loop over the number of endogenous variables to plot
F1=figure(1)
set(F1, 'numbertitle','off')
set(F1, 'name', 'Impulse response functions (Technology Shock)')
for j = 1:6;
    subplot(3,2,j), plot(irfs_RBC_1_epsA (j,:),'-k', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','k'); hold on;
                    plot(irfs_RBC_2_epsA(j,:),'--r', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','r'); hold on;
                    xlabel('Quarters');
                    ylabel('% dev from SS');
                    grid on
                    title(names(j,:),'FontSize',10)
axis tight;                   
end;

legend('No Investment Costs', 'Investment Costs', 'Location', 'Best', 'Orientation', 'horizontal')

%Options for the plot
h=figure('Position', [600, 0, 1000, 900]);
axes ('position', [0, 0, 1, 1]);

%Figure for the government spending shock
%Loop over the number of endogenous variables to plot
F1=figure(2)
set(F1, 'numbertitle','off')
set(F1, 'name', 'Impulse response functions (Gov Spending Shock)')
for j = 1:6;
    subplot(3,2,j), plot(irfs_RBC_1_epsG (j,:),'-k', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','k'); hold on;
                    plot(irfs_RBC_2_epsG(j,:),'--r', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','r'); hold on;
                    xlabel('Quarters');
                    ylabel('% dev from SS');
                    grid on
                    title(names(j,:),'FontSize',10)
axis tight;                   
end;

legend('No Investment Costs', 'Investment Costs', 'Location', 'Best', 'Orientation', 'horizontal')




