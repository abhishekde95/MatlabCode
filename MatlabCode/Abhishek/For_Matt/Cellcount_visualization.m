% Plotting the cell count from the CSV files

% Closing all figures 
close all; 

% Clearing all variables 
clearvars;

% Loading the data from Green Channel - PV signal 
T_PV = readtable('C:\Users\setup\Google Drive\UW\DLX_20X_monkey\V1_V2_border_Virusy3\PVab_counts_part2_10X.csv');

% Loading the data from the Red Channel - mDLX signal
T_ChR2 = readtable('C:\Users\setup\Google Drive\UW\DLX_20X_monkey\V1_V2_border_Virusy3\mDLX_counts_part2_10X.csv');

% Plotting the figure 
figure(1);
plot(T_PV.X,3666-T_PV.Y,'o','MarkerFaceColor',[0 1 0],'MarkerSize',5); hold on;
plot(T_ChR2.X,3666-T_ChR2.Y,'+','MarkerEdgeColor',[0 0 0],'MarkerSize',5);
plot(T_ChR2.X,3666-T_ChR2.Y,'o','MarkerFaceColor',[1 0 0],'MarkerSize',3); axis equal; 
set(gca,'Xlim',[0 4924],'Ylim',[0 3666],'Tickdir','out');
xlabel('X'), ylabel('Y'); legend('PV','','ChR2'); hold off;

% Calculating the selectivity and the associated error in measurement
p2 = 468/543;
err2 = 1.96*sqrt(p*(1-p)/543);
