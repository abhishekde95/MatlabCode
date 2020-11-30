% Script for comparing DLX vs AAV-PV1
% IHC analyses and cell counting of a 10X tissue, Virusy 3, 10X, Part 2
% Author - Abhishek De, 5/20

if ~exist('plot_counter')
    plot_counter = 1;
end
T_PV = readtable('C:\Users\setup\Google Drive\UW\DLX_20X_monkey_IHC\V1_V2_border_Virusy3_comparing_DLX_AAVPV1\PVab_counts_part2_10X.csv');
T_ChR2 = readtable('C:\Users\setup\Google Drive\UW\DLX_20X_monkey_IHC\V1_V2_border_Virusy3_comparing_DLX_AAVPV1\mDLX_counts_part2_10X.csv');
T_AAVPV1 = readtable('C:\Users\setup\Google Drive\UW\DLX_20X_monkey_IHC\V1_V2_border_Virusy3_comparing_DLX_AAVPV1\AAV_PV1_counts_part2_10X.csv');

figure(plot_counter); set(gcf,'Name','DLX vs AAV PV1')
subplot(221); plot(T_PV.X,3666-T_PV.Y,'o','MarkerFaceColor',[0 1 0],'MarkerSize',5,'MarkerEdgeColor',[1 1 1]); hold on;
plot(T_ChR2.X,3666-T_ChR2.Y,'+','MarkerEdgeColor',[0 0 0],'MarkerSize',5);
plot(T_ChR2.X,3666-T_ChR2.Y,'o','MarkerFaceColor',[1 0 0],'MarkerSize',3,'MarkerEdgeColor',[1 1 1]); axis equal;
set(gca,'Xlim',[0 4924],'Ylim',[0 3666],'Tickdir','out'); xlabel('X'), ylabel('Y'); legend('PV','','ChR2'); hold off;

subplot(222); plot(T_PV.X,3666-T_PV.Y,'o','MarkerFaceColor',[0 1 0],'MarkerSize',5,'MarkerEdgeColor',[1 1 1]); hold on;
plot(T_AAVPV1.X,3666-T_AAVPV1.Y,'+','MarkerEdgeColor',[0 0 0],'MarkerSize',5);
plot(T_AAVPV1.X,3666-T_AAVPV1.Y,'o','MarkerFaceColor',[0 0.5 1.0],'MarkerSize',3,'MarkerEdgeColor',[1 1 1]); axis equal;
set(gca,'Xlim',[0 4924],'Ylim',[0 3666],'Tickdir','out'); xlabel('X'), ylabel('Y'); legend('PV','','AAVPV1'); hold off;

subplot(223); plot(T_AAVPV1.X,3666-T_AAVPV1.Y,'o','MarkerFaceColor',[0 0.5 1.0],'MarkerSize',4,'MarkerEdgeColor',[1 1 1]); hold on;
plot(T_ChR2.X,3666-T_ChR2.Y,'o','MarkerFaceColor',[1 0 0],'MarkerSize',3,'MarkerEdgeColor',[1 1 1]); axis equal;
set(gca,'Xlim',[0 4924],'Ylim',[0 3666],'Tickdir','out'); xlabel('X'), ylabel('Y'); legend('PV','AAVPV1'); hold off;

subplot(224); plot(T_PV.X,3666-T_PV.Y,'o','MarkerFaceColor',[0 1 0],'MarkerSize',5,'MarkerEdgeColor',[1 1 1]); hold on;
plot(T_AAVPV1.X,3666-T_AAVPV1.Y,'o','MarkerFaceColor',[0 0.5 1.0],'MarkerSize',4,'MarkerEdgeColor',[1 1 1]);
plot(T_ChR2.X,3666-T_ChR2.Y,'o','MarkerFaceColor',[1 0 0],'MarkerSize',3,'MarkerEdgeColor',[1 1 1]); axis equal;
set(gca,'Xlim',[0 4924],'Ylim',[0 3666],'Tickdir','out'); xlabel('X'), ylabel('Y'); legend('PV','AAVPV1','ChR2'); hold off;
plot_counter = plot_counter + 1;

% Comaprison of dlx and PV
p1 = 468/543;
err1 = 1.96*sqrt(p1*(1-p1)/543);
pval1 = 1-binocdf(468,543,0.75);

% Comparison of AAV-PV1 and PV 
p2 = 459/615;
err2 = 1.96*sqrt(p2*(1-p2)/615);
pval2 = 1-binocdf(459,615,0.75);