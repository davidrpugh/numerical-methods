%%This code produces a set of graphs that plots the IRFs produced by the linear vs
%%the non linear version of the dynare code
%%(c) CIMS Univeristy of Surrey
%%The Science and  Art of DSGE Modelling: Construction, Calibration, Estimation and Policy
 

clear;
close all;
clc;

%Variables names to plot
names = char('Output', 'Consumption', 'Investment', 'Hours worked', 'Leverage', 'Spread', 'Nominal interest rate', 'Tobin Q', 'Inflation');


%load non linear simulations
load GK_RES_Course_results


%Rename the IRFs for each variable of interest in the non linear model
Y_epsA     = oo_.irfs.YY_epsA;
C_epsA     = oo_.irfs.CC_epsA;
I_epsA     = oo_.irfs.II_epsA;
h_epsA     = oo_.irfs.hh_epsA;
Leverage_epsA    = oo_.irfs.phiphi_epsA;
R_epsA     = oo_.irfs.spreadspread_epsA;
Rn_epsA     = oo_.irfs.RnRn_epsA;
Q_epsA     = oo_.irfs.QQ_epsA;
PIE_epsA     = oo_.irfs.PIEPIE_epsA;

Y_epsG     = oo_.irfs.YY_epsG;
C_epsG     = oo_.irfs.CC_epsG;
I_epsG     = oo_.irfs.II_epsG;
h_epsG     = oo_.irfs.hh_epsG;
Leverage_epsG    = oo_.irfs.phiphi_epsG;
R_epsG     = oo_.irfs.spreadspread_epsG;
Rn_epsG     = oo_.irfs.RnRn_epsG;
Q_epsG     = oo_.irfs.QQ_epsG;
PIE_epsG     = oo_.irfs.PIEPIE_epsG;

Y_epsMS     = oo_.irfs.YY_epsMS;
C_epsMS     = oo_.irfs.CC_epsMS;
I_epsMS     = oo_.irfs.II_epsMS;
h_epsMS     = oo_.irfs.hh_epsMS;
Leverage_epsMS    = oo_.irfs.phiphi_epsMS;
R_epsMS     = oo_.irfs.spreadspread_epsMS;
Rn_epsMS     = oo_.irfs.RnRn_epsMS;
Q_epsMS     = oo_.irfs.QQ_epsMS;
PIE_epsMS     = oo_.irfs.PIEPIE_epsMS;

Y_epsMPS     = oo_.irfs.YY_epsMPS;
C_epsMPS     = oo_.irfs.CC_epsMPS;
I_epsMPS     = oo_.irfs.II_epsMPS;
h_epsMPS     = oo_.irfs.hh_epsMPS;
Leverage_epsMPS    = oo_.irfs.phiphi_epsMPS;
R_epsMPS     = oo_.irfs.spreadspread_epsMPS;
Rn_epsMPS     = oo_.irfs.RnRn_epsMPS;
Q_epsMPS     = oo_.irfs.QQ_epsMPS;
PIE_epsMPS     = oo_.irfs.PIEPIE_epsMPS;

Y_epsAtrend     = oo_.irfs.YY_epsAtrend;
C_epsAtrend     = oo_.irfs.CC_epsAtrend;
I_epsAtrend     = oo_.irfs.II_epsAtrend;
h_epsAtrend     = oo_.irfs.hh_epsAtrend;
Leverage_epsAtrend    = oo_.irfs.phiphi_epsAtrend;
R_epsAtrend     = oo_.irfs.spreadspread_epsAtrend;
Rn_epsAtrend     = oo_.irfs.RnRn_epsAtrend;
Q_epsAtrend     = oo_.irfs.QQ_epsAtrend;
PIE_epsAtrend     = oo_.irfs.PIEPIE_epsAtrend;


Y_epscapqual     = oo_.irfs.YY_epscapqual;
C_epscapqual     = oo_.irfs.CC_epscapqual;
I_epscapqual     = oo_.irfs.II_epscapqual;
h_epscapqual     = oo_.irfs.hh_epscapqual;
Leverage_epscapqual    = oo_.irfs.phiphi_epscapqual;
R_epscapqual     = oo_.irfs.spreadspread_epscapqual;
Rn_epscapqual     = oo_.irfs.RnRn_epscapqual;
Q_epscapqual     = oo_.irfs.QQ_epscapqual;
PIE_epscapqual     = oo_.irfs.PIEPIE_epscapqual;

%Create a vector of IRFs per each shock
irfs_NK1_epsA =[Y_epsA;  C_epsA;  I_epsA; h_epsA;  Leverage_epsA;  R_epsA; Rn_epsA; Q_epsA; PIE_epsA];
irfs_NK1_epsG =[Y_epsG;  C_epsG;  I_epsG; h_epsG;  Leverage_epsG;  R_epsG; Rn_epsG; Q_epsG; PIE_epsG];
irfs_NK1_epsMS =[Y_epsMS;  C_epsMS;  I_epsMS; h_epsMS;  Leverage_epsMS;  R_epsMS; Rn_epsMS; Q_epsMS; PIE_epsMS];
irfs_NK1_epsMPS =[Y_epsMPS;  C_epsMPS;  I_epsMPS; h_epsMPS;  Leverage_epsMPS;  R_epsMPS; Rn_epsMPS; Q_epsMPS; PIE_epsMPS];
irfs_NK1_epsAtrend =[Y_epsAtrend;  C_epsAtrend;  I_epsAtrend; h_epsAtrend;  Leverage_epsAtrend;  R_epsAtrend; Rn_epsAtrend; Q_epsAtrend; PIE_epsAtrend];
irfs_NK1_epscapqual =[Y_epscapqual;  C_epscapqual;  I_epscapqual; h_epscapqual;  Leverage_epscapqual;  R_epscapqual; Rn_epscapqual; Q_epscapqual; PIE_epscapqual];
%load linear simulations

%load non linear simulations
load BGG_RES_Course_results

%Rename the IRFs for each variable of interest in the non linear model
%Rename the IRFs for each variable of interest in the non linear model
Y_epsA     = oo_.irfs.YY_epsA;
C_epsA     = oo_.irfs.CC_epsA;
I_epsA     = oo_.irfs.II_epsA;
h_epsA     = oo_.irfs.hh_epsA;
Leverage_epsA    = oo_.irfs.phiphi_epsA;
R_epsA     = oo_.irfs.spreadspread_epsA;
Rn_epsA     = oo_.irfs.RnRn_epsA;
Q_epsA     = oo_.irfs.QQ_epsA;
PIE_epsA     = oo_.irfs.PIEPIE_epsA;

Y_epsG     = oo_.irfs.YY_epsG;
C_epsG     = oo_.irfs.CC_epsG;
I_epsG     = oo_.irfs.II_epsG;
h_epsG     = oo_.irfs.hh_epsG;
Leverage_epsG    = oo_.irfs.phiphi_epsG;
R_epsG     = oo_.irfs.spreadspread_epsG;
Rn_epsG     = oo_.irfs.RnRn_epsG;
Q_epsG     = oo_.irfs.QQ_epsG;
PIE_epsG     = oo_.irfs.PIEPIE_epsG;

Y_epsMS     = oo_.irfs.YY_epsMS;
C_epsMS     = oo_.irfs.CC_epsMS;
I_epsMS     = oo_.irfs.II_epsMS;
h_epsMS     = oo_.irfs.hh_epsMS;
Leverage_epsMS    = oo_.irfs.phiphi_epsMS;
R_epsMS     = oo_.irfs.spreadspread_epsMS;
Rn_epsMS     = oo_.irfs.RnRn_epsMS;
Q_epsMS     = oo_.irfs.QQ_epsMS;
PIE_epsMS     = oo_.irfs.PIEPIE_epsMS;

Y_epsMPS     = oo_.irfs.YY_epsMPS;
C_epsMPS     = oo_.irfs.CC_epsMPS;
I_epsMPS     = oo_.irfs.II_epsMPS;
h_epsMPS     = oo_.irfs.hh_epsMPS;
Leverage_epsMPS    = oo_.irfs.phiphi_epsMPS;
R_epsMPS     = oo_.irfs.spreadspread_epsMPS;
Rn_epsMPS     = oo_.irfs.RnRn_epsMPS;
Q_epsMPS     = oo_.irfs.QQ_epsMPS;
PIE_epsMPS     = oo_.irfs.PIEPIE_epsMPS;

Y_epsAtrend     = oo_.irfs.YY_epsAtrend;
C_epsAtrend     = oo_.irfs.CC_epsAtrend;
I_epsAtrend     = oo_.irfs.II_epsAtrend;
h_epsAtrend     = oo_.irfs.hh_epsAtrend;
Leverage_epsAtrend    = oo_.irfs.phiphi_epsAtrend;
R_epsAtrend     = oo_.irfs.spreadspread_epsAtrend;
Rn_epsAtrend     = oo_.irfs.RnRn_epsAtrend;
Q_epsAtrend     = oo_.irfs.QQ_epsAtrend;
PIE_epsAtrend     = oo_.irfs.PIEPIE_epsAtrend;


Y_epscapqual     = oo_.irfs.YY_epscapqual;
C_epscapqual     = oo_.irfs.CC_epscapqual;
I_epscapqual     = oo_.irfs.II_epscapqual;
h_epscapqual     = oo_.irfs.hh_epscapqual;
Leverage_epscapqual    = oo_.irfs.phiphi_epscapqual;
R_epscapqual     = oo_.irfs.spreadspread_epscapqual;
Rn_epscapqual     = oo_.irfs.RnRn_epscapqual;
Q_epscapqual     = oo_.irfs.QQ_epscapqual;
PIE_epscapqual     = oo_.irfs.PIEPIE_epscapqual;

%Create a vector of IRFs per each shock
irfs_NK2_epsA =[Y_epsA;  C_epsA;  I_epsA; h_epsA;  Leverage_epsA;  R_epsA; Rn_epsA; Q_epsA; PIE_epsA];
irfs_NK2_epsG =[Y_epsG;  C_epsG;  I_epsG; h_epsG;  Leverage_epsG;  R_epsG; Rn_epsG; Q_epsG; PIE_epsG];
irfs_NK2_epsMS =[Y_epsMS;  C_epsMS;  I_epsMS; h_epsMS;  Leverage_epsMS;  R_epsMS; Rn_epsMS; Q_epsMS; PIE_epsMS];
irfs_NK2_epsMPS =[Y_epsMPS;  C_epsMPS;  I_epsMPS; h_epsMPS;  Leverage_epsMPS;  R_epsMPS; Rn_epsMPS; Q_epsMPS; PIE_epsMPS];
irfs_NK2_epsAtrend =[Y_epsAtrend;  C_epsAtrend;  I_epsAtrend; h_epsAtrend;  Leverage_epsAtrend;  R_epsAtrend; Rn_epsAtrend; Q_epsAtrend; PIE_epsAtrend];
irfs_NK2_epscapqual =[Y_epscapqual;  C_epscapqual;  I_epscapqual; h_epscapqual;  Leverage_epscapqual;  R_epscapqual; Rn_epscapqual; Q_epscapqual; PIE_epscapqual];

%load linear simulations

%load non linear simulations
load KM_RES_Course_results

%Rename the IRFs for each variable of interest in the non linear model
%Rename the IRFs for each variable of interest in the non linear model
Y_epsA     = oo_.irfs.YY_epsA;
C_epsA     = oo_.irfs.CC_epsA;
I_epsA     = oo_.irfs.II_epsA;
h_epsA     = oo_.irfs.hh_epsA;
Leverage_epsA    = oo_.irfs.phiphi_epsA;
R_epsA     = oo_.irfs.spreadspread_epsA;
Rn_epsA     = oo_.irfs.RnRn_epsA;
Q_epsA     = oo_.irfs.QQ_epsA;
PIE_epsA     = oo_.irfs.PIEPIE_epsA;

Y_epsG     = oo_.irfs.YY_epsG;
C_epsG     = oo_.irfs.CC_epsG;
I_epsG     = oo_.irfs.II_epsG;
h_epsG     = oo_.irfs.hh_epsG;
Leverage_epsG    = oo_.irfs.phiphi_epsG;
R_epsG     = oo_.irfs.spreadspread_epsG;
Rn_epsG     = oo_.irfs.RnRn_epsG;
Q_epsG     = oo_.irfs.QQ_epsG;
PIE_epsG     = oo_.irfs.PIEPIE_epsG;

Y_epsMS     = oo_.irfs.YY_epsMS;
C_epsMS     = oo_.irfs.CC_epsMS;
I_epsMS     = oo_.irfs.II_epsMS;
h_epsMS     = oo_.irfs.hh_epsMS;
Leverage_epsMS    = oo_.irfs.phiphi_epsMS;
R_epsMS     = oo_.irfs.spreadspread_epsMS;
Rn_epsMS     = oo_.irfs.RnRn_epsMS;
Q_epsMS     = oo_.irfs.QQ_epsMS;
PIE_epsMS     = oo_.irfs.PIEPIE_epsMS;

Y_epsMPS     = oo_.irfs.YY_epsMPS;
C_epsMPS     = oo_.irfs.CC_epsMPS;
I_epsMPS     = oo_.irfs.II_epsMPS;
h_epsMPS     = oo_.irfs.hh_epsMPS;
Leverage_epsMPS    = oo_.irfs.phiphi_epsMPS;
R_epsMPS     = oo_.irfs.spreadspread_epsMPS;
Rn_epsMPS     = oo_.irfs.RnRn_epsMPS;
Q_epsMPS     = oo_.irfs.QQ_epsMPS;
PIE_epsMPS     = oo_.irfs.PIEPIE_epsMPS;

Y_epsAtrend     = oo_.irfs.YY_epsAtrend;
C_epsAtrend     = oo_.irfs.CC_epsAtrend;
I_epsAtrend     = oo_.irfs.II_epsAtrend;
h_epsAtrend     = oo_.irfs.hh_epsAtrend;
Leverage_epsAtrend    = oo_.irfs.phiphi_epsAtrend;
R_epsAtrend     = oo_.irfs.spreadspread_epsAtrend;
Rn_epsAtrend     = oo_.irfs.RnRn_epsAtrend;
Q_epsAtrend     = oo_.irfs.QQ_epsAtrend;
PIE_epsAtrend     = oo_.irfs.PIEPIE_epsAtrend;


Y_epscapqual     = oo_.irfs.YY_epscapqual;
C_epscapqual     = oo_.irfs.CC_epscapqual;
I_epscapqual     = oo_.irfs.II_epscapqual;
h_epscapqual     = oo_.irfs.hh_epscapqual;
Leverage_epscapqual    = oo_.irfs.phiphi_epscapqual;
R_epscapqual     = oo_.irfs.spreadspread_epscapqual;
Rn_epscapqual     = oo_.irfs.RnRn_epscapqual;
Q_epscapqual     = oo_.irfs.QQ_epscapqual;
PIE_epscapqual     = oo_.irfs.PIEPIE_epscapqual;

%Create a vector of IRFs per each shock
irfs_NK3_epsA =[Y_epsA;  C_epsA;  I_epsA; h_epsA;  Leverage_epsA;  R_epsA; Rn_epsA; Q_epsA; PIE_epsA];
irfs_NK3_epsG =[Y_epsG;  C_epsG;  I_epsG; h_epsG;  Leverage_epsG;  R_epsG; Rn_epsG; Q_epsG; PIE_epsG];
irfs_NK3_epsMS =[Y_epsMS;  C_epsMS;  I_epsMS; h_epsMS;  Leverage_epsMS;  R_epsMS; Rn_epsMS; Q_epsMS; PIE_epsMS];
irfs_NK3_epsMPS =[Y_epsMPS;  C_epsMPS;  I_epsMPS; h_epsMPS;  Leverage_epsMPS;  R_epsMPS; Rn_epsMPS; Q_epsMPS; PIE_epsMPS];
irfs_NK3_epsAtrend =[Y_epsAtrend;  C_epsAtrend;  I_epsAtrend; h_epsAtrend;  Leverage_epsAtrend;  R_epsAtrend; Rn_epsAtrend; Q_epsAtrend; PIE_epsAtrend];
irfs_NK3_epscapqual =[Y_epscapqual;  C_epscapqual;  I_epscapqual; h_epscapqual;  Leverage_epscapqual;  R_epscapqual; Rn_epscapqual; Q_epscapqual; PIE_epscapqual];


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
                    plot(irfs_NK3_epsA(j,:),':b', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','b'); hold on;
                    xlabel('Quarters');
                    ylabel('% dev from SS');
                    grid on
                    title(names(j,:),'FontSize',10)
axis tight;                   
end;

legend('GK', 'BGG ', 'KM','Best', 'Orientation', 'horizontal')

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
                    plot(irfs_NK3_epsG(j,:),':b', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','b'); hold on;
                    
                    xlabel('Quarters');
                    ylabel('% dev from SS');
                    grid on
                    title(names(j,:),'FontSize',10)
axis tight;                   
end;

legend('GK', 'BGG', 'KM', 'Best', 'Orientation', 'horizontal')


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
                    plot(irfs_NK3_epsMS(j,:),':b', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','b'); hold on;
                    
                    xlabel('Quarters');
                    ylabel('% dev from SS');
                    grid on
                    title(names(j,:),'FontSize',10)
axis tight;                   
end;

legend('GK', 'BGG ','KM', 'Location', 'Best', 'Orientation', 'horizontal')

F1=figure(4)
set(F1, 'numbertitle','off')
set(F1, 'name', 'Impulse response functions (Monetary Policy Shock)')
for j = 1:9;
    subplot(3,3,j), plot(irfs_NK1_epsMPS (j,:),'-k', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','k'); hold on;
                    plot(irfs_NK2_epsMPS(j,:),'--r', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','r'); hold on;
                    plot(irfs_NK3_epsMPS(j,:),':b', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','b'); hold on;
                    
                   xlabel('Quarters');
                    ylabel('% dev from SS');
                    grid on
                    title(names(j,:),'FontSize',10)
axis tight;                   
end;

legend('GK', 'BGG','KM', 'Location', 'Best', 'Orientation', 'horizontal')

F1=figure(5)
set(F1, 'numbertitle','off')
set(F1, 'name', 'Impulse response functions (Trend Shock)')
for j = 1:9;
    subplot(3,3,j), plot(irfs_NK1_epsAtrend (j,:),'-k', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','k'); hold on;
                    plot(irfs_NK2_epsAtrend(j,:),'--r', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','r'); hold on;
                    plot(irfs_NK3_epsAtrend(j,:),':b', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','b'); hold on;
                    
                   xlabel('Quarters');
                    ylabel('% dev from SS');
                    grid on
                    title(names(j,:),'FontSize',10)
axis tight;                   
end;

legend('GK', 'BGG', 'KM','Location', 'Best', 'Orientation', 'horizontal')

F1=figure(6)
set(F1, 'numbertitle','off')
set(F1, 'name', 'Impulse response functions (Capital Quality Shock)')
for j = 1:9;
    subplot(3,3,j), plot(irfs_NK1_epscapqual (j,:),'-k', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','k'); hold on;
                    plot(irfs_NK2_epscapqual(j,:),'--r', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','r'); hold on;
                    plot(irfs_NK3_epscapqual(j,:),':b', 'LineWidth',2, 'MarkerSize', 5,'MarkerEdgeColor','k', 'MarkerFaceColor','b'); hold on;
                    
                   xlabel('Quarters');
                    ylabel('% dev from SS');
                    grid on
                    title(names(j,:),'FontSize',10)
axis tight;                   
end;

legend('GK', 'BGG','KM', 'Location', 'Best', 'Orientation', 'horizontal')