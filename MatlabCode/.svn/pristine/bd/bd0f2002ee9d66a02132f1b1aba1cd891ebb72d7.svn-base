function GLMSDGUI_Surfaces(~,~)
global GLMSD

% This analysis compares GridLMSubunitDetection surfaces to the
% corresponding GridLMSubunit neurometric surfaces.
% Aug 2014      Created.        JPW

figure(20); clf;
set(gcf,'Position',[200 200 1000 750],...
    'NumberTitle','off','Name',['Surface Analysis (' GLMSD.datafile ')'])
SurfFig = get(gcf,'UserData');

% Set up disp panels
SurfFig.disp.psychpanel = uipanel('Parent',gcf,'Units','normalized',...
    'position',[.025 .35 .3 .625],'title','Psychometric Data');
SurfFig.disp.neuropanel = uipanel('parent',gcf,'units','normalized',...
    'position',[.35 .35 .3 .625],'title','Neurometric Data');
SurfFig.disp.diffpanel = uipanel('parent',gcf,'units','normalized',...
    'position',[.675 .35 .3 .625],'title','Difference');

% Set up stats panels
SurfFig.stats.psych = uipanel('parent',gcf,'units','normalized',...
    'position',[.025 .025 .3 .3],'title','Psychometric Stats');
SurfFig.stats.neuro = uipanel('parent',gcf,'units','normalized',...
    'position',[.35 .025 .3 .3],'title','Neurometric Stats');

% set up control panel
SurfFig.conpanel = uipanel('parent',gcf,'units','normalized',...
    'position',[.675 .025 .3 .3],'title','Control Panel');

% Set up disppanel axes
% psych axes
psychpanel.axes.oneD = axes('Parent',SurfFig.disp.psychpanel,'Units','normalized',...
    'Position',[.2 .55 .7 .35],'CameraPosition',[-.9 -.9 5],...
    'XLim',[-max(GLMSD.rho) max(GLMSD.rho)],...
    'YLim',[-max(GLMSD.rho) max(GLMSD.rho)],...
    'ZLim',[0 1.1],'XGrid','on','YGrid','on','ZGrid','on','box','on');
title('1D Fit','FontSize',14,'FontWeight','bold')
xlabel('Lcc'); ylabel('Mcc'); zlabel('Percentage Correct')

psychpanel.axes.twoD = axes('Parent',SurfFig.disp.psychpanel,'Units','normalized',...
    'Position',[.2 .1 .7 .35],'CameraPosition',[-.9 -.9 5],...
    'XLim',[-max(GLMSD.rho) max(GLMSD.rho)],...
    'YLim',[-max(GLMSD.rho) max(GLMSD.rho)],...
    'ZLim',[0 1.1],'XGrid','on','YGrid','on','ZGrid','on','box','on');
title('2D Fit','FontSize',14,'FontWeight','bold')
xlabel('Lcc'); ylabel('Mcc'); zlabel('Percentage Correct')

% neuro axes
neuropanel.axes.oneD = axes('Parent',SurfFig.disp.neuropanel,'Units','normalized',...
    'Position',[.2 .55 .7 .35],'CameraPosition',[-.9 -.9 5],...
    'XLim',[-max(GLMSD.GLMSdata.rho) max(GLMSD.GLMSdata.rho)],...
    'YLim',[-max(GLMSD.GLMSdata.rho) max(GLMSD.GLMSdata.rho)],...
    'ZLim',[0 1.1],'XGrid','on','YGrid','on','ZGrid','on','box','on');
title('1D Fit','FontSize',14,'FontWeight','bold')
xlabel('Lcc'); ylabel('Mcc'); zlabel('Predicted Percentage Correct')

neuropanel.axes.twoD = axes('Parent',SurfFig.disp.neuropanel,'Units','normalized',...
    'Position',[.2 .1 .7 .35],'CameraPosition',[-.9 -.9 5],...
    'XLim',[-max(GLMSD.GLMSdata.rho) max(GLMSD.GLMSdata.rho)],...
    'YLim',[-max(GLMSD.GLMSdata.rho) max(GLMSD.GLMSdata.rho)],...
    'ZLim',[0 1.1],'XGrid','on','YGrid','on','ZGrid','on','box','on');
title('2D Fit','FontSize',14,'FontWeight','bold')
xlabel('Lcc'); ylabel('Mcc'); zlabel('Predicted Percentage Correct')

% Diff axes
diffpanel.axes.oneD = axes('Parent',SurfFig.disp.diffpanel,'Units','normalized',...
    'Position',[.2 .55 .7 .35],'CameraPosition',[-.9 -.9 5],...
    'XLim',[-max(GLMSD.GLMSdata.rho) max(GLMSD.GLMSdata.rho)],...
    'YLim',[-max(GLMSD.GLMSdata.rho) max(GLMSD.GLMSdata.rho)],...
    'ZLim',[0 1.1],'XGrid','on','YGrid','on','ZGrid','on','box','on');
title('Error between 1D Fits','FontSize',14,'FontWeight','bold')
xlabel('Lcc'); ylabel('Mcc'); zlabel('2D - 1D')

diffpanel.axes.twoD = axes('Parent',SurfFig.disp.diffpanel,'Units','normalized',...
    'Position',[.2 .1 .7 .35],'CameraPosition',[-.9 -.9 5],...
    'XLim',[-max(GLMSD.GLMSdata.rho) max(GLMSD.GLMSdata.rho)],...
    'YLim',[-max(GLMSD.GLMSdata.rho) max(GLMSD.GLMSdata.rho)],...
    'ZLim',[0 1.1],'XGrid','on','YGrid','on','ZGrid','on','box','on');
title('Error Between 2D Fits','FontSize',14,'FontWeight','bold')
xlabel('Lcc'); ylabel('Mcc'); zlabel('2D - 1D')


%%% Set up statspanel axes
% Axes for displaying parameter names
textvals{1} = 'A =';
textvals{2} = 'Sigma =';
textvals{3} = 'Exponent =';
textvals{4} = 'Baseline =';
textvals{5} = 'Angle =';
textvals{6} = [];
textvals{7} = 'Conf Int =';
textvals = textvals';
psychstats.paramNames = uicontrol('style','text','parent',SurfFig.stats.psych,...
    'units','normalized','position',[.01 .4 .25 .5],'string',textvals,...
    'fontsize',12,'horizontalalignment','center');
neurostats.paramName = uicontrol('style','text','parent',SurfFig.stats.neuro,...
    'units','normalized','position',[.01 .4 .25 .5],'string',textvals,...
    'fontsize',12,'horizontalalignment','center');

% Set up panel for displaying 1D parameter values
psychstats.oneD.paramPanel = uipanel('parent',SurfFig.stats.psych,'units','normalized',...
    'position',[.27 .4 .3 .59],'title','1D');
neurostats.oneD.paramPanel = uipanel('parent',SurfFig.stats.neuro,'units','normalized',...
    'position',[.27 .4 .3 .59],'title','1D');

% Set up edit box for displaying 1D parameter values
psychstats.oneD.paramDisp = uicontrol('style','text','parent',psychstats.oneD.paramPanel,...
    'units','normalized','position',[.01 .01 .98 .9],...
    'fontsize',12,'horizontalalignment','center');
neurostats.oneD.paramDisp = uicontrol('style','text','parent',neurostats.oneD.paramPanel,...
    'units','normalized','position',[.01 .01 .98 .9],...
    'fontsize',12,'horizontalalignment','center');

% Set up panel for displaying 2D parameter values
psychstats.twoD.paramPanel = uipanel('parent',SurfFig.stats.psych,'units','normalized',...
    'position',[.58 .4 .41 .59],'title','2D');
neurostats.twoD.paramPanel = uipanel('parent',SurfFig.stats.neuro,'units','normalized',...
    'position',[.58 .4 .41 .59],'title','2D');

% Set up edit box for displaying 2D parameter values
psychstats.twoD.paramDisp = uicontrol('style','text','parent',psychstats.twoD.paramPanel,...
    'units','normalized','position',[.01 .01 .98 .9],...
    'fontsize',12,'horizontalalignment','center');
neurostats.twoD.paramDisp = uicontrol('style','text','parent',neurostats.twoD.paramPanel,...
    'units','normalized','position',[.01 .01 .98 .9],...
    'fontsize',12,'horizontalalignment','center');

% Set up panel for displaying comparison values
psychstats.compare.axes = uipanel('parent',SurfFig.stats.psych,'units','normalized',...
    'position',[.01 .01 .98 .38],'title','1D vs 2D');
neurostats.compare.axes = uipanel('parent',SurfFig.stats.neuro,'units','normalized',...
    'position',[.01 .01 .98 .38],'title','1D vs 2D');

% Set up editable box for displaying comparison values
psychstats.compare.params = uicontrol('style','text','parent',psychstats.compare.axes,...
    'units','normalized','position',[.01 .01 .98 .98],...
    'fontsize',12,'horizontalalignment','left');
neurostats.compare.params = uicontrol('style','text','parent',neurostats.compare.axes,...
    'units','normalized','position',[.01 .01 .98 .98],...
    'fontsize',12,'horizontalalignment','left');

% Set up conpanel axes
conpanel.uicontrols.run = uicontrol('parent',SurfFig.conpanel,'units','pixels',...
    'position',[200 10 80 40],'style','pushbutton',...
    'string','Run','callback',@RunAnalysis);

% Save user variables
set(SurfFig.disp.psychpanel,'UserData',psychpanel)
set(SurfFig.disp.neuropanel,'UserData',neuropanel)
set(SurfFig.disp.diffpanel,'UserData',diffpanel)
set(SurfFig.stats.psych,'UserData',psychstats)
set(SurfFig.stats.neuro,'UserData',neurostats)
set(SurfFig.conpanel,'UserData',conpanel)
set(gcf,'UserData',SurfFig)

end


function RunAnalysis(~,~)

disp('Starting Analysis...')

%OneDPsych();
TwoDPsych();
%Compare1D2D('psych');
OneDNeuro();
%TwoDNeuro();
%Compare1D2D('neuro');
%ComputeDiff;
%CompareBestFits();
CompareBestAxis();

disp('done')

end


function OneDPsych()
global GLMSD

disp('Fitting a 1D function to Psychometric Data...')

% Load Figure Variables
SurfFig = get(gcf,'UserData');
conpanel = get(SurfFig.conpanel,'UserData');
psychpanel = get(SurfFig.disp.psychpanel,'UserData');
psychstats = get(SurfFig.stats.psych,'UserData');

% Set up x,y,z vals
xvals = GLMSD.Lcc;
yvals = GLMSD.Mcc;
zvals = GLMSD.AnsCorrect;

% Set up some variables
angs = linspace(0,pi,9);
angs(end) = [];
%nrots = numel(angs);
%GOF = nan(nrots,1);
vlb = [1 .001 .001 .01  .5 -pi];
vub = [1     50   50   5   .5   pi];
options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
Aguess = 1;
sigguess = [.04 .3 50];
expguess = 2;
blguess = .5;
guessIdx = fullfact([numel(sigguess) numel(sigguess) numel(angs)]);
GOF = nan(size(guessIdx,1),1);
params = nan(numel(guessIdx),numel(vlb));
    
% Rotate through angles
for rot = 1:size(guessIdx,1)
    
    % Generating an initial guess
    %params0 = [upperA, sigmaguess, sigmaguess, expguess, baselineguess, angs(rot)];
    params0 = [Aguess sigguess(guessIdx(rot,1)) sigguess(guessIdx(rot,2))...
        expguess blguess angs(guessIdx(rot,3))];
    
    % Fit all rotated points
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',params0,[],[],[],[],vlb,vub,[],options,[xvals yvals],zvals,'surface7','bernoulli');
    params(rot,:) = f1;
    GOF(rot) = fval;
    
end

% Find the best fitting surface
[~,bestIdx] = min(GOF);
params1 = params(bestIdx,:);

% Determine confidence intervals
thanks2greg = @(a) FitNakaRushtonFunJPW(a,[xvals yvals],zvals,'surface7','bernoulli');
[hessval,~] = hessian(thanks2greg,params1);
parvar = 1/hessval(end,end);
confint = (2* sqrt(parvar))/pi*180;
psychstats.conf.oneD = confint;

% Display surface and pts
x = linspace(-max(xvals),max(xvals),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params1,[xx(:) yy(:)],'surface7');
surface = reshape(surface,size(xx));
axes(psychpanel.axes.oneD); cla; hold on;
%psychpanel.surf.oneD = surfc(xx,yy,surface);
psychpanel.surf.oneD = contour(xx,yy,surface);
alpha(.4)
psychpanel.pts.oneD = plot3(GLMSD.uniqueLcc,GLMSD.uniqueMcc,GLMSD.pCorrect,'k*');
drawnow;

% Pull out parameters for display
statstext{1} = rnddec(params1(1),2);
statstext{2} = [num2str(rnddec(params1(2),2)) '  ' num2str(rnddec(params1(3),2))];
%statstext{3} = [num2str(rnddec(params1(4),2)) '  'num2str(rnddec(params1(5),2))]; % For surface 6
statstext{3} = num2str(rnddec(params1(4),2));
statstext{4} = rnddec(params1(5),2);
statstext{5} = rnddec(params1(6)/pi*180,2);
statstext{6} = [];
%predpts = ComputeNakaRushtonJPW(params1,[xvals yvals],'surface7');
%chisq = sum((zvals - predpts).^2);
statstext{7} = confint;
statstext = statstext';

% Display and save stats
psychstats.oneD.params = params1;
psychstats.oneD.negLL = GOF(bestIdx);
set(psychstats.oneD.paramDisp,'string',statstext);

% Set up a 'right click enlarge' option
hcmenu = uicontextmenu;
uimenu(hcmenu,'label','enlarge','callback',@Popout);
axes(psychpanel.axes.oneD);
set(gca,'uicontextmenu',hcmenu);
set(get(gca,'children'),'hittest','off')

% Save user variables
set(SurfFig.disp.psychpanel,'UserData',psychpanel)
set(SurfFig.stats.psych,'UserData',psychstats)
set(SurfFig.conpanel,'UserData',conpanel)
set(gcf,'UserData',SurfFig)

end

function OneDNeuro()
global GLMSD

disp('Fitting a 1D function to Neurometric Data...')


% Load Figure Variables
SurfFig = get(gcf,'UserData');
conpanel = get(SurfFig.conpanel,'UserData');
neuropanel = get(SurfFig.disp.neuropanel,'UserData');
neurostats = get(SurfFig.stats.neuro,'UserData');

% Set x,y,z values
xvals = GLMSD.GLMSdata.uniqueLcc;
yvals = GLMSD.GLMSdata.uniqueMcc;
zvals = GLMSD.AUC;
    
% Set up some variables
angs = linspace(0,pi,9);
angs(end) = [];
%nrots = numel(angs);
%GOF = nan(nrots,1);
vlb = [1 .001 .001 .001  .5  -pi];
vub = [1     10    10   10  .5   pi];
options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
Aguess = 1;
sigguess = 0.2;
expguess = 2;
blguess = .5;
guessIdx = fullfact([numel(sigguess) numel(sigguess) numel(angs)]);
GOF = nan(size(guessIdx,1),1);
params = nan(numel(guessIdx),numel(vlb));

% Rotate through angles
for rot = 1:size(guessIdx,1)
    
    % Generating an initial guess
    %params0 = [upperA, sigmaguess, sigmaguess, expguess, baselineguess, angs(rot)];
    params0 = [Aguess sigguess(guessIdx(rot,1)) sigguess(guessIdx(rot,2))...
        expguess blguess angs(guessIdx(rot,3))];
    
    % Fit all rotated points
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',params0,[],[],[],[],vlb,vub,[],options,[xvals yvals],zvals,'surface7','bernoulli');
    params(rot,:) = f1;
    GOF(rot) = fval;
    
end

% Find the best fitting surface
[~,bestIdx] = min(GOF);
params1 = params(bestIdx,:);

% Determine confidence intervals
thanks2greg = @(a) FitNakaRushtonFunJPW(a,[xvals yvals],zvals,'surface7','bernoulli');
[hessval,~] = hessian(thanks2greg,params1);
parvar = 1/hessval(end,end);
confint = (2* sqrt(parvar))/pi*180;
neurostats.conf.oneD = confint;

% Display surface and pts
x = linspace(-max(xvals),max(xvals),50);
[xx yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params1,[xx(:) yy(:)],'surface7');
surface = reshape(surface,size(xx));
axes(neuropanel.axes.oneD); cla; hold on;
neuropanel.pts.oneD = plot3(xvals,yvals,zvals,'k*');
%neuropanel.surf.oneD = surfc(xx,yy,surface);
neuropanel.surf.oneD = contour(xx,yy,surface);
alpha(.4); drawnow

% Pull out parameters for display
statstext{1} = rnddec(params1(1),2);
statstext{2} = [num2str(rnddec(params1(2),2)) '  ' num2str(rnddec(params1(3),2))];
%statstext{3} = [num2str(rnddec(params1(4),2)) '  ' num2str(rnddec(params1(5),2))];
statstext{3} = num2str(rnddec(params1(4),2));
statstext{4} = rnddec(params1(5),2);
statstext{5} = rnddec(params1(6)/pi*180,2);
statstext{6} = [];
%predpts = ComputeNakaRushtonJPW(params1,[xvals yvals],'surface7');
%chisq = sum((zvals - predpts).^2);
statstext{7} = confint;
statstext = statstext';

% Display and save stats
neurostats.oneD.params = params1;
neurostats.oneD.negLL = GOF(bestIdx);
set(neurostats.oneD.paramDisp,'string',statstext);

% Set up a 'right click enlarge' option
hcmenu = uicontextmenu;
uimenu(hcmenu,'label','enlarge','callback',@Popout);
axes(neuropanel.axes.oneD);
set(gca,'uicontextmenu',hcmenu);
set(get(gca,'children'),'hittest','off')

% Save user variables
set(SurfFig.disp.neuropanel,'UserData',neuropanel)
set(SurfFig.stats.neuro,'UserData',neurostats)
set(SurfFig.conpanel,'UserData',conpanel)
set(gcf,'UserData',SurfFig)

end


function TwoDPsych()
global GLMSD

disp('Fitting a 2D function to Psychometric Data...')

% Load Figure Variables
SurfFig = get(gcf,'UserData');
conpanel = get(SurfFig.conpanel,'UserData');
psychpanel = get(SurfFig.disp.psychpanel,'UserData');
psychstats = get(SurfFig.stats.psych,'UserData');

% Set up x,y,z vals
xvals = GLMSD.Lcc;
yvals = GLMSD.Mcc;
zvals = GLMSD.AnsCorrect;
    
% Set up some variables
angs = linspace(0,pi,5);
angs(end) = [];
%nrots = numel(angs);
%GOF = nan(nrots,1);
vlb = [1 0.001 0.001 0.001 0.001 .001  .5 -pi];
vub = [1  10    10    10    10    10   .5  pi];
options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
Aguess = 1;
sigguess = [.04 0.4 10];
expguess = 2;
blguess = .5;
guessIdx = fullfact([numel(sigguess) numel(sigguess) numel(sigguess) numel(sigguess) numel(angs)]);
GOF = nan(size(guessIdx,1),1);
params = nan(size(guessIdx,1),numel(vlb));
    
% Rotate through angles
for rot = 1:size(guessIdx,1)
    
    % Generating an initial guess
    params0 = [Aguess sigguess(guessIdx(rot,1)) sigguess(guessIdx(rot,2))...
        sigguess(guessIdx(rot,3)) sigguess(guessIdx(rot,4))...
        expguess blguess angs(guessIdx(rot,5))];
    
    % Fit all rotated points
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',params0,[],[],[],[],vlb,vub,[],options,[xvals yvals],zvals,'surface8','bernoulli');
    params(rot,:) = f1;
    GOF(rot) = fval;
    
end
    
% Find the best fitting surface
[~,bestIdx] = min(GOF);
params1 = params(bestIdx,:);

% Determine confidence intervals
thanks2greg = @(a) FitNakaRushtonFunJPW(a,[xvals yvals],zvals,'surface8','bernoulli');
[hessval,~] = hessian(thanks2greg,params1);
parvar = 1/hessval(end,end);
confint = (2* sqrt(parvar))/pi*180;
psychstats.conf.twoD = confint;

% Display surface and pts
x = linspace(-max(xvals),max(xvals),50);
[xx yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params1,[xx(:) yy(:)],'surface8');
surface = reshape(surface,size(xx));
axes(psychpanel.axes.twoD); cla; hold on;
psychpanel.pts.twoD = plot3(GLMSD.uniqueLcc,GLMSD.uniqueMcc,GLMSD.pCorrect,'k*');
%psychpanel.surf.twoD = surfc(xx,yy,surface);
psychpanel.surf.twoD = contour(xx,yy,surface);
alpha(.4); drawnow;

% Pull out parameters for display
statstext{1} = rnddec(params1(1),2);
statstext{2} = [num2str(rnddec(params1(2),2)) '  '...
    num2str(rnddec(params1(3),2)) '  '...
    num2str(rnddec(params1(4),2)) '  '...
    num2str(rnddec(params1(5),2))];
statstext{3} = num2str(rnddec(params1(6),2));
statstext{4} = rnddec(params1(7),2);
statstext{5} = rnddec(params1(8)/pi*180,2);
statstext{6} = [];
%predpts = ComputeNakaRushtonJPW(params1,[xvals yvals],'surface8');
%chisq = sum((zvals - predpts).^2);
statstext{7} = confint;
statstext = statstext';

% Display and save stats
psychstats.twoD.params = params1;
psychstats.twoD.negLL = GOF(bestIdx);
set(psychstats.twoD.paramDisp,'string',statstext);

% Set up a 'right click enlarge' option
hcmenu = uicontextmenu;
uimenu(hcmenu,'label','enlarge','callback',@Popout);
axes(psychpanel.axes.twoD);
set(gca,'uicontextmenu',hcmenu);
set(get(gca,'children'),'hittest','off')

% Save user variables
set(SurfFig.disp.psychpanel,'UserData',psychpanel)
set(SurfFig.stats.psych,'UserData',psychstats)
set(SurfFig.conpanel,'UserData',conpanel)
set(gcf,'UserData',SurfFig)

end


function TwoDNeuro()
global GLMSD

disp('Fitting a 2D function to Neurometric Data...')

% Load Figure Variables
SurfFig = get(gcf,'UserData');
conpanel = get(SurfFig.conpanel,'UserData');
neuropanel = get(SurfFig.disp.neuropanel,'UserData');
neurostats = get(SurfFig.stats.neuro,'UserData');

% Set up x,y,z vals
xvals = GLMSD.GLMSdata.uniqueLcc;
yvals = GLMSD.GLMSdata.uniqueMcc;
zvals = GLMSD.AUC;

% Set up some variables
angs = linspace(0,pi,5);
angs(end) = [];
%nrots = numel(angs);
%GOF = nan(nrots,1);
vlb = [1 0.001 0.001 0.001 0.001 .001  .5 -pi];
vub = [1  10    10    10    10    10   .5  pi];
options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
Aguess = 1;
sigguess = [.04 0.4 10];
expguess = 2;
blguess = .5;
guessIdx = fullfact([numel(sigguess) numel(sigguess) numel(sigguess) numel(sigguess) numel(angs)]);
GOF = nan(size(guessIdx,1),1);
params = nan(size(guessIdx,1),numel(vlb));
    
% Rotate through angles
for rot = 1:size(guessIdx,1)
    
    % Generating an initial guess
    params0 = [Aguess sigguess(guessIdx(rot,1)) sigguess(guessIdx(rot,2))...
        sigguess(guessIdx(rot,3)) sigguess(guessIdx(rot,4))...
        expguess blguess angs(guessIdx(rot,5))];
    
    % Fit all rotated points
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',params0,[],[],[],[],vlb,vub,[],options,[xvals yvals],zvals,'surface8','bernoulli');
    params(rot,:) = f1;
    GOF(rot) = fval;
    
end
    
% Find the best fitting surface
[~,bestIdx] = min(GOF);
params1 = params(bestIdx,:);

% Determine confidence intervals
thanks2greg = @(a) FitNakaRushtonFunJPW(a,[xvals yvals],zvals,'surface8','bernoulli');
[hessval,~] = hessian(thanks2greg,params1);
parvar = 1/hessval(end,end);
confint = (2* sqrt(parvar))/pi*180;
neurostats.conf.twoD = confint;
%set(conpanel.paramvals.oneD.conf,'string',confint);

% Display surface and pts
x = linspace(-max(xvals),max(xvals),50);
[xx yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params1,[xx(:) yy(:)],'surface8');
surface = reshape(surface,size(xx));
axes(neuropanel.axes.twoD); cla; hold on;
neuropanel.pts.twoD = plot3(xvals,yvals,zvals,'k*');
%neuropanel.surf.twoD = surfc(xx,yy,surface);
neuropanel.surf.twoD = contour(xx,yy,surface);
alpha(.4); drawnow;
    
% Pull out parameters for display
statstext{1} = rnddec(params1(1),2);
statstext{2} = [num2str(rnddec(params1(2),2)) '  '...
    num2str(rnddec(params1(3),2)) '  '...
    num2str(rnddec(params1(4),2)) '  '...
    num2str(rnddec(params1(5),2))];
statstext{3} = num2str(rnddec(params1(6),2));
statstext{4} = rnddec(params1(7),2);
statstext{5} = rnddec(params1(8)/pi*180,2);
statstext{6} = [];
%predpts = ComputeNakaRushtonJPW(params1,[xvals yvals],'surface8');
%chisq = sum((zvals - predpts).^2);
statstext{7} = confint;
statstext = statstext';

% Display and save stats
neurostats.twoD.params = params1;
neurostats.twoD.negLL = GOF(bestIdx);
set(neurostats.twoD.paramDisp,'string',statstext);
    
% Set up a 'right click enlarge' option
hcmenu = uicontextmenu;
uimenu(hcmenu,'label','enlarge','callback',@Popout);
axes(neuropanel.axes.twoD);
set(gca,'uicontextmenu',hcmenu);
set(get(gca,'children'),'hittest','off')

% Save user variables
set(SurfFig.disp.neuropanel,'UserData',neuropanel)
set(SurfFig.stats.neuro,'UserData',neurostats)
set(SurfFig.conpanel,'UserData',conpanel)
set(gcf,'UserData',SurfFig)

end


function Popout(a,~)
global GLMSD

% Load Figure Variables
SurfFig = get(gcf,'UserData');
psychpanel = get(SurfFig.disp.psychpanel,'UserData');
neuropanel = get(SurfFig.disp.neuropanel,'UserData');
diffpanel = get(SurfFig.disp.diffpanel,'UserData');

% Set up figure
figure(21); clf; hold on;
set(gcf,'position',[400 200 700 700])
set(gca,'CameraPosition',[-.9 -.9 5],'box','on','fontsize',16,...
    'XLim',[-max(GLMSD.rho) max(GLMSD.rho)],...
    'YLim',[-max(GLMSD.rho) max(GLMSD.rho)],...
    'ZLim',[0 1.1],'XGrid','on','YGrid','on','ZGrid','on');
xlabel('Lcc'); ylabel('Mcc');

% figure out where Popout was called
if get(a,'parent') == get(psychpanel.axes.oneD,'uicontextmenu')
    copyobj(psychpanel.pts.oneD,gca);
    copyobj(psychpanel.surf.oneD,gca);
    title('1D Fit to Psychometric Data','FontSize',18,'FontWeight','bold')
    zlabel('Percentage Correct')
    
elseif get(a,'parent') == get(neuropanel.axes.oneD,'uicontextmenu')
    copyobj(neuropanel.pts.oneD,gca)
    copyobj(neuropanel.surf.oneD,gca)
    title('1D Fit to Neurometric Data','FontSize',14,'FontWeight','bold')
    zlabel('Predicted Percentage Correct')
    
elseif get(a,'parent') == get(diffpanel.axes.oneD,'uicontextmenu')
    copyobj(diffpanel.surf.oneD,gca)
    title('Difference Between 1D Fits','FontSize',14,'FontWeight','bold')
    zlabel('Error (psych - neuro) ')
    zvals = get(diffpanel.surf.oneD,'zdata');
    zvals = zvals{1};
    set(gca,'zlim',[min(zvals(:)) max(zvals(:))]);
    
elseif get(a,'parent') == get(psychpanel.axes.twoD,'uicontextmenu')
    copyobj(psychpanel.pts.twoD,gca);
    copyobj(psychpanel.surf.twoD,gca);
    title('2D Fit to Psychometric Data','FontSize',18,'FontWeight','bold')
    zlabel('Percentage Correct')
    
elseif get(a,'parent') == get(neuropanel.axes.twoD,'uicontextmenu')
    copyobj(neuropanel.pts.twoD,gca);
    copyobj(neuropanel.surf.twoD,gca);
    title('2D Fit to Neurometric Data','FontSize',18,'FontWeight','bold')
    zlabel('Predicted Percentage Correct')
    
elseif get(a,'parent') == get(diffpanel.axes.twoD,'uicontextmenu')
    %copyobj(diffpanel.pts.twoD,gca)
    %copyobj(diffpanel.surf.twoD,gca)
    %title('Difference Between 2D Fits','FontSize',14,'FontWeight','bold')
    %zlabel('Error (psych - neuro)')
    %zvals = get(diffpanel.surf.twoD,'zdata');
    %zvals = zvals{1};
    %set(gca,'zlim',[min(zvals(:)) max(zvals(:))]);
    set(gca,'view',[0 90],'xlimmode','auto','ylim',[.4 1.1])
    copyobj(diffpanel.line.neuro,gca)
    copyobj(diffpanel.line.psych,gca)
    legend('Neuro','Psycho','location','southeast')
    xlabel('Contrast')
    ylabel('Percentage Correct')
    
end

end

function ComputeDiff()
global GLMSD

% Load Figure Variables
SurfFig = get(gcf,'UserData');
diffpanel = get(SurfFig.disp.diffpanel,'UserData');
psychstats = get(SurfFig.stats.psych,'UserData');
neurostats = get(SurfFig.stats.neuro,'UserData');

% Display difference between 1D psych and neuro
params1 = psychstats.oneD.params;
params2 = neurostats.oneD.params;

x = linspace(-max(GLMSD.rho),max(GLMSD.rho),50);
[xx yy] = meshgrid(x,x);
surface1 = ComputeNakaRushtonJPW(params1,[xx(:) yy(:)],'surface7');
surface2 = ComputeNakaRushtonJPW(params2,[xx(:) yy(:)],'surface7');
surface1 = reshape(surface1,size(xx));
surface2 = reshape(surface2,size(xx));

surfDiff = (surface2 - surface1);
axes(diffpanel.axes.oneD); cla; hold on;
diffpanel.surf.oneD = surfc(xx,yy,surfDiff);
drawnow;
set(gca,'zlim',[min(surfDiff(:))*1.1 max(surfDiff(:))*1.1])

% Set up a 'right click enlarge' option
hcmenu = uicontextmenu;
uimenu(hcmenu,'label','enlarge','callback',@Popout);
axes(diffpanel.axes.oneD)
set(gca,'uicontextmenu',hcmenu);
set(get(gca,'children'),'hittest','off')

% Display difference between 2D psych and neuro
params1 = psychstats.twoD.params;
params2 = neurostats.twoD.params;

x = linspace(-max(GLMSD.rho),max(GLMSD.rho),50);
[xx yy] = meshgrid(x,x);
surface1 = ComputeNakaRushtonJPW(params1,[xx(:) yy(:)],'surface8');
surface2 = ComputeNakaRushtonJPW(params2,[xx(:) yy(:)],'surface8');
surface1 = reshape(surface1,size(xx));
surface2 = reshape(surface2,size(xx));

surfDiff = (surface2 - surface1);
axes(diffpanel.axes.twoD); cla; hold on;
diffpanel.surf.twoD = surfc(xx,yy,surfDiff);
drawnow;
set(gca,'zlim',[min(surfDiff(:))*1.1 max(surfDiff(:))*1.1])

% Set up a 'right click enlarge' option
hcmenu = uicontextmenu;
uimenu(hcmenu,'label','enlarge','callback',@Popout);
axes(diffpanel.axes.twoD);
set(gca,'uicontextmenu',hcmenu);
set(get(gca,'children'),'hittest','off')

% Save User Variables
set(SurfFig.disp.diffpanel,'UserData',diffpanel);
set(gcf,'UserData',SurfFig);

end

function CompareBestFits()
global GLMSD

% Load Figure Variables
SurfFig = get(gcf,'UserData');
diffpanel = get(SurfFig.disp.diffpanel,'UserData');
psychstats = get(SurfFig.stats.psych,'UserData');
neurostats = get(SurfFig.stats.neuro,'UserData');

% Actual statistical test
pval = psychstats.compare.pval;
if pval < .05
    params1 = psychstats.twoD.params;
    x = linspace(-max(GLMSD.rho),max(GLMSD.rho),50);
    [xx yy] = meshgrid(x,x);
    surface1 = ComputeNakaRushtonJPW(params1,[xx(:) yy(:)],'surface8');
else
    params1 = psychstats.oneD.params;
    x = linspace(-max(GLMSD.rho),max(GLMSD.rho),50);
    [xx yy] = meshgrid(x,x);
    surface1 = ComputeNakaRushtonJPW(params1,[xx(:) yy(:)],'surface7');
end

pval = neurostats.compare.pval;
if pval < .05
    params1 = neurostats.twoD.params;
    x = linspace(-max(GLMSD.rho),max(GLMSD.rho),50);
    [xx yy] = meshgrid(x,x);
    surface2 = ComputeNakaRushtonJPW(params1,[xx(:) yy(:)],'surface8');
else
    params1 = neurostats.oneD.params;
    x = linspace(-max(GLMSD.rho),max(GLMSD.rho),50);
    [xx yy] = meshgrid(x,x);
    surface2 = ComputeNakaRushtonJPW(params1,[xx(:) yy(:)],'surface7');
end

surface1 = reshape(surface1,size(xx));
surface2 = reshape(surface2,size(xx));

surfDiff = (surface1 - surface2);
axes(diffpanel.axes.oneD); cla; hold on;
diffpanel.surf.oneD = surfc(xx,yy,surfDiff);
drawnow;
set(gca,'zlim',[min(surfDiff(:))*1.1 max(surfDiff(:))*1.1])
zlabel('Psych - Neuro')
axes(diffpanel.axes.twoD); cla;

% Set up a 'right click enlarge' option
hcmenu = uicontextmenu;
uimenu(hcmenu,'label','enlarge','callback',@Popout);
axes(diffpanel.axes.oneD);
set(gca,'uicontextmenu',hcmenu);
set(get(gca,'children'),'hittest','off')

% Save user variables
set(SurfFig.disp.diffpanel,'UserData',diffpanel)
set(SurfFig.stats.psych,'UserData',psychstats)
set(SurfFig.stats.neuro,'UserData',neurostats)
set(gcf,'UserData',SurfFig)

end

function CompareBestAxis()
global GLMSD
disp('Comparing Best Neuro Axis to Behavior...')

% Load Figure Variables
SurfFig = get(gcf,'UserData');
diffpanel = get(SurfFig.disp.diffpanel,'UserData');
psychstats = get(SurfFig.stats.psych,'UserData');
neurostats = get(SurfFig.stats.neuro,'UserData');


paramsN = neurostats.oneD.params;
paramsP = psychstats.twoD.params;

rot = paramsN(end);
rotMat = [cos(rot) sin(rot); -sin(rot) cos(rot)];
tempRotPts = [rotMat * [GLMSD.Lcc GLMSD.Mcc]']';
maxx = max(tempRotPts(:,1));
xvals = linspace(-maxx,maxx,50);
yvals = zeros(size(xvals));
zvalsN = ComputeNakaRushtonJPW(paramsN,[xvals' yvals'],'surface7');
zvalsP = ComputeNakaRushtonJPW(paramsP,[xvals' yvals'],'surface8');

axes(diffpanel.axes.twoD); cla; hold on;
set(diffpanel.axes.twoD,'view',[0 90])
xlim([-maxx maxx]); ylim([.4 1.1]);
title('Best 1D Neuro Axis')
diffpanel.line.neuro = plot(xvals,zvalsN,'r');
diffpanel.line.psych = plot(xvals,zvalsP,'g');
legend('Neuro','Psycho','location','southeast')
drawnow;
xlabel('Contrast')
zlabel('Percentage Correct')

% Set up a 'right click enlarge' option
hcmenu = uicontextmenu;
uimenu(hcmenu,'label','enlarge','callback',@Popout);
axes(diffpanel.axes.twoD);
set(gca,'uicontextmenu',hcmenu);
set(get(gca,'children'),'hittest','off')

% Save user variables
set(SurfFig.disp.diffpanel,'UserData',diffpanel)
set(SurfFig.stats.psych,'UserData',psychstats)
set(SurfFig.stats.neuro,'UserData',neurostats)
set(gcf,'UserData',SurfFig)



end


function Compare1D2D(which)
global GLMSD

disp('Comparing 1D and 2D fits...')

% Load Figure Variables
SurfFig = get(gcf,'UserData');
psychstats = get(SurfFig.stats.psych,'UserData');
neurostats = get(SurfFig.stats.neuro,'UserData');

if strcmp(which,'psych')
    
    % Prepare parameters
    nfp1 = numel(psychstats.oneD.params);
    nfp2 = numel(psychstats.twoD.params);
    LL1 = -psychstats.oneD.negLL;
    LL2 = -psychstats.twoD.negLL;
    
    % Actual statistical test
    diffLL = 2 * (LL2 - LL1);
    pval = 1 - chi2cdf(diffLL,nfp2-nfp1);
    psychstats.compare.pval = pval;
    %set(conpanel.paramvals.pval,'string',pval);

    stats{1} = [];
    stats{2} = ['df1 = ' num2str(nfp1) '          ' 'df2 = ' num2str(nfp2)];
    stats{3} = [];
    stats{4} = ['p-value = ' num2str(rnddec(pval,3))];
    stats = stats';
    set(psychstats.compare.params,'string',stats);

    % Save values
    psychstats.compare.pval = pval;
    %psychstats.compare.Fstat = Fstats;
    psychstats.compare.df1 = nfp1;
    psychstats.compare.df2 = nfp2;

    stats{1} = [];
    stats{2} = ['nfp1 = ' num2str(nfp1) '       nfp2 = ' num2str(nfp2)];
    stats{3} = ['LL1 = ' num2str(LL1) '       LL2 = ' num2str(LL2)];
    stats{4} = ['P-value = ' num2str(psychstats.compare.pval)];
    stats = stats';
    set(psychstats.compare.params,'string',stats)

elseif strcmp(which,'neuro')
    
    % Prepare parameters
    nfp1 = numel(neurostats.oneD.params);
    nfp2 = numel(neurostats.twoD.params);
    LL1 = -neurostats.oneD.negLL;
    LL2 = -neurostats.twoD.negLL;
    
    % Actual statistical test
    diffLL = 2 * (LL2 - LL1);
    pval = 1 - chi2cdf(diffLL,nfp2-nfp1);
    neurostats.compare.pval = pval;
    %set(conpanel.paramvals.pval,'string',pval);
    
    stats{1} = [];
    stats{2} = ['nfp1 = ' num2str(nfp1) '       nfp2 = ' num2str(nfp2)];
    stats{3} = ['LL1 = ' num2str(neurostats.oneD.negLL) '       LL2 = ' num2str(neurostats.twoD.negLL)];
    stats{4} = ['P-value = ' num2str(neurostats.compare.pval)];
    stats = stats';
    set(neurostats.compare.params,'string',stats)
    
end

% Save user variables
set(SurfFig.stats.psych,'UserData',psychstats)
set(SurfFig.stats.neuro,'UserData',neurostats)
set(gcf,'UserData',SurfFig)

    
end