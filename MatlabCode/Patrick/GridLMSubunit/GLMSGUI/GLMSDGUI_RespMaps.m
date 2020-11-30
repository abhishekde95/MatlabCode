function [] = GLMSDGUI_RespMaps(~,~)
global GLMSD

figure(10); clf;
set(gcf,'Position',[100 100 700 700])

MapsFig = get(gcf,'UserData');
    
% Setting this manually, but should be put into datafile
whichsub = 1;

% Set up axes
MapsFig.axes.PsychData = axes('Parent',gcf,'Units','normalized',...
    'Position',[.1 .55 .35 .35],'CameraPosition',[-.9 -.9 5],...
    'XLim',[-max(GLMSD.GLMSdata.rho) max(GLMSD.GLMSdata.rho)],...
    'YLim',[-max(GLMSD.GLMSdata.rho) max(GLMSD.GLMSdata.rho)],...
    'ZLim',[0 1.1],'XGrid','on','YGrid','on','ZGrid','on');
title('Neurometric Data','FontSize',14,'FontWeight','bold'); hold on;
xlabel('Lcc'); ylabel('Mcc'); zlabel('Predicted Percentage Correct')

MapsFig.axes.AllData = axes('Parent',gcf,'Units','normalized',...
    'Position',[.55 .55 .35 .35],'CameraPosition',[-.9 -.9 5],...
    'XLim',[-max(GLMSD.rho) max(GLMSD.rho)],...
    'YLim',[-max(GLMSD.rho) max(GLMSD.rho)],...
    'ZLim',[0 1.1],'XGrid','on','YGrid','on','ZGrid','on');
title('All Psychometric Data','FontSize',14,'FontWeight','bold'); hold on;
xlabel('Lcc'); ylabel('Mcc'); zlabel('Percentage Correct')

MapsFig.axes.RF1 = axes('Parent',gcf,'Units','normalized',...
    'Position',[.1 .1 .35 .35],'CameraPosition',[-.9 -.9 5],...
    'XLim',[-max(GLMSD.RF1.rho) max(GLMSD.RF1.rho)],...
    'YLim',[-max(GLMSD.RF1.rho) max(GLMSD.RF1.rho)],...
    'ZLim',[0 1.1],'XGrid','on','YGrid','on','ZGrid','on');
title('RF1 Data','FontSize',14,'FontWeight','bold'); hold on;
xlabel('Lcc'); ylabel('Mcc'); zlabel('Percentage Correct')

MapsFig.axes.RF2 = axes('Parent',gcf,'Units','normalized',...
    'Position',[.55 .1 .35 .35],'CameraPosition',[-.9 -.9 5],...
    'XLim',[-max(GLMSD.RF2.rho) max(GLMSD.RF2.rho)],...
    'YLim',[-max(GLMSD.RF2.rho) max(GLMSD.RF2.rho)],...
    'ZLim',[0 1.1],'XGrid','on','YGrid','on','ZGrid','on');
title('RF2 Data','FontSize',14,'FontWeight','bold'); hold on;
xlabel('Lcc'); ylabel('Mcc'); zlabel('Percentage Correct')    
    
% Neurometric Data
% GLMSD.AUC = nan(numel(GLMSD.GLMSdata.subunit{whichsub}.uniqueLcc),1);
% for s = 1:numel(GLMSD.GLMSdata.subunit{whichsub}.uniqueLcc)
%     idx = GLMSD.GLMSdata.subunit{whichsub}.uniqueIdx{s};
%     spHist = GLMSD.GLMSdata.subunit{whichsub}.fr(idx)';
%     blHist = GLMSD.GLMSdata.subunit{whichsub}.blfr';
%     thresh = unique([spHist blHist]);
%     TPR = nan(1,length(thresh));
%     FPR = TPR;
%     for n=1:length(thresh)
%         TPR(n) = sum(spHist >= thresh(n))/length(spHist);
%         FPR(n) = sum(blHist >= thresh(n))/length(blHist);
%     end
%     GLMSD.AUC(s) = trapz(fliplr(FPR),fliplr(TPR));
% end
x = GLMSD.GLMSdata.subunit{whichsub}.uniqueLcc;
y = GLMSD.GLMSdata.subunit{whichsub}.uniqueMcc;
z =  GLMSD.AUC;
F = TriScatteredInterp(x,y,z);
[qx,qy] = meshgrid(min(x):.01:max(x),min(y):.01:max(y));
qz = F(qx,qy);

axes(MapsFig.axes.PsychData); hold on;
MapsFig.data.Psych.surf = surfc(qx,qy,qz);
MapsFig.data.Psych.pts = plot3(x,y,z,'k*');


% All Psychometric data
x = GLMSD.uniqueLcc;
y = GLMSD.uniqueMcc;
z = GLMSD.pCorrect;
F = TriScatteredInterp(x,y,z);
[qx,qy] = meshgrid(min(x):.01:max(x),min(y):.01:max(y));

axes(MapsFig.axes.AllData)
MapsFig.data.AllData.pts = plot3(x,y,z,'k*');
MapsFig.data.AllData.surf = surfc(qx,qy,F(qx,qy));


% RF1
x = GLMSD.RF1.uniqueLcc;
y = GLMSD.RF1.uniqueMcc;
z = GLMSD.RF1.pCorrect;
F = TriScatteredInterp(x,y,z);
[qx,qy] = meshgrid(min(x):.01:max(x),min(y):.01:max(y));

axes(MapsFig.axes.RF1)
MapsFig.data.RF1.pts = plot3(x,y,z,'k*');
MapsFig.data.RF1.surf = surfc(qx,qy,F(qx,qy));

% RF2
x = GLMSD.RF2.uniqueLcc;
y = GLMSD.RF2.uniqueMcc;
z = GLMSD.RF2.pCorrect;
F = TriScatteredInterp(x,y,z);
[qx,qy] = meshgrid(min(x):.01:max(x),min(y):.01:max(y));

axes(MapsFig.axes.RF2)
MapsFig.data.RF2.pts = plot3(x,y,z,'k*');
MapsFig.data.RF2.surf = surfc(qx,qy,F(qx,qy));


% Save User Data
set(gcf,'UserData',MapsFig)

end