function [] = GLMSGUI_CRF(~,~)
global GLMP


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  CRF Analysis Functions  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sets up figure for analysis of contrast response functions along a single axis.

% Set up figure
figure(40); clf;
set(gcf,'numbertitle','off','name',['GLMP CRFs (' GLMP.datafile ')'],'pos',[150 100 600 700])
CRF.controls = uipanel('parent',gcf,'units','normalized',...
    'pos',[.025 .65 .95 .325]);
cpanel = get(CRF.controls,'UserData');

% Set up stimulus map
cpanel.axes.stimmap = axes('parent',CRF.controls,'units','pixels',...
    'pos',[20 10 200 200]);
cpanel.axes.stim = polar(GLMP.subunit{1}.uniquetheta,GLMP.subunit{1}.uniquerho,'ro',...
    'parent',cpanel.axes.stimmap);
set(cpanel.axes.stim,'ButtonDownFcn',@CRFStimulusSelectionSingle)

% Set up subunit buttons
strs = {'Subunit #1'};
if ~isempty(GLMP.subunit{2})
    strs = cat(1,strs,'Subunit #2');
end
if ~isempty(GLMP.subunit{3})
    if GLMP.subunit{3}.gridX{1} == 100
        strs = cat(1,strs,'Gabor');
    else
        strs = cat(1,strs,'Subunit #3');
    end
end
cpanel.uicontrols.subunitmenu = uicontrol('parent',CRF.controls,...
    'style','popupmenu','string',strs,'fontsize',12,...
    'units','normalized','pos',[.45 .8 .225 .1],...
    'Callback',@CRFsubsel);

% Set up buttons to compare single stim to many stim
strs = {'Single Stimulus' 'Compare Stimuli'};
cpanel.uicontrols.comparemenu = uicontrol('parent',CRF.controls,...
    'style','popupmenu','string',strs,'fontsize',12,...
    'units','normalized','pos',[.45 .6 .225 .1],...
    'Callback',@singleVmulti);

% For selecting which stimuli to compare
strs = {'Sttimulus #1' 'Stimulus #2' 'Stimulus #3' 'Stimulus #4'};
cpanel.uicontrols.multistimmenu = uicontrol('parent',CRF.controls,...
    'style','popupmenu','string',strs,'fontsize',12,...
    'units','normalized','pos',[.45 .4 .225 .1],...
    'string',strs,'fontsize',12,'enable','off','UserData',nan(numel(strs),1));

% Set up display for fit values
cpanel.parvals.A = uicontrol('parent',CRF.controls,'style','edit',...
    'units','normalized','pos',[.85 .85 .1 .1],'fontsize',12);
cpanel.parvals.sig = uicontrol('parent',CRF.controls,'style','edit',...
    'units','normalized','pos',[.85 .65 .1 .1],'fontsize',12);
cpanel.parvals.exp = uicontrol('parent',CRF.controls,'style','edit',...
    'units','normalized','pos',[.85 .45 .1 .1],'fontsize',12);
cpanel.parvals.bl = uicontrol('parent',CRF.controls,'style','edit',...
    'units','normalized','pos',[.85 .25 .1 .1],'fontsize',12);
cpanel.parvals.kappa = uicontrol('parent',CRF.controls,'style','edit',...
    'units','normalized','pos',[.85 .05 .1 .1],'fontsize',12);

% Set up labels for fit values
cpanel.parlabels.A = uicontrol('parent',CRF.controls,'style','text',...
    'units','normalized','pos',[.7 .85 .125 .1],...
    'string','A =','fontsize',12,'horizontalalignment','center');
cpanel.parlabels.sig = uicontrol('parent',CRF.controls,'style','text',...
    'units','normalized','pos',[.7 .65 .125 .1],...
    'string','Sigma =','fontsize',12,'horizontalalignment','center');
cpanel.parlabels.exp = uicontrol('parent',CRF.controls,'style','text',...
    'units','normalized','pos',[.7 .45 .125 .1],...
    'string','Exponent =','fontsize',12,'horizontalalignment','center');
cpanel.parlabels.bl = uicontrol('parent',CRF.controls,'style','text',...
    'units','normalized','pos',[.7 .25 .125 .1],...
    'string','Basline =','fontsize',12,'horizontalalignment','center');
cpanel.parlabels.kappa = uicontrol('parent',CRF.controls,'style','text',...
    'units','normalized','pos',[.7 .05 .125 .1],...
    'string','Kappa =','fontsize',12,'horizontalalignment','center');

% Set up CRF Panel
CRF.display = uipanel('parent',gcf,'units','normalized',...
    'pos',[.025 .025 .95 .6]);
dpanel = get(CRF.display,'UserData');
dpanel.CRF = axes('parent',CRF.display,'units','normalized',...
    'pos',[.05 .05 .9 .9],'Visible','on','box','on');
dpanel.single.allstim = nan;
dpanel.single.meanstim = nan;
dpanel.single.NKRfun = nan;
cpanel.single.stimsel = nan(1,2);
dpanel.multi.meanstim = nan(4,1);
dpanel.multi.NKRfun = nan(4,1);
cpanel.multi.stimsel = nan(4,2);

% Save all settings
set(CRF.controls,'UserData',cpanel)
set(CRF.display,'UserData',dpanel)
set(gcf,'UserData',CRF)

end


function singleVmulti(~,~)
% Switch between the two options.

% Load variables
CRF = get(gcf,'UserData');
cpanel = get(CRF.controls,'UserData');

% 1 is single, 2 is compare
if cpanel.uicontrols.comparemenu.Value == 1
    resetCRFsingle();
else
    resetCRFmulti();
end

% Not saving variables here bc they are manipulated inside the 'reset'
% functions.  If saved again here, they would be saved as the
% previous, pre-reset variables.

end

function resetCRFsingle()
% Reset the figure to display only 1 stimulus PSTH and its individual rasters

CRF = get(gcf,'UserData');
cpanel = get(CRF.controls,'UserData');
dpanel = get(CRF.display,'UserData');

% Disenable multi stimulus comparison buttons
set(cpanel.uicontrols.multistimmenu,'enable','off')

% Show previously chosen single CRF and stimuli
if isobject(dpanel.single.allstim)
    set(dpanel.single.allstim,'visible','on')
    set(dpanel.single.meanstim,'visible','on')
    set(dpanel.single.NKRfun,'visible','on')
    set(cpanel.single.selstim,'visible','on')
end

% Hide previously chosen multi CRF and stimuli
L = ~isnan(dpanel.multi.meanstim);
set(dpanel.multi.meanstim(L),'visible','off')
set(dpanel.multi.NKRfun(L),'visible','off')
set(cpanel.multi.stimsel(L,:),'visible','off')
set(cpanel.parvals.A,'string',[])
set(cpanel.parvals.sig,'string',[])
set(cpanel.parvals.exp,'string',[])
set(cpanel.parvals.bl,'string',[])
set(cpanel.parvals.kappa,'string',[])

% Set the correct callback function
set(cpanel.axes.stim,'ButtonDownFcn',@CRFStimulusSelectionSingle)

% Save variables
set(CRF.controls,'UserData',cpanel);
set(CRF.display,'UserData',dpanel);
set(gcf,'UserData',CRF);

end

function resetCRFmulti()
% Reset the figure to display only 1 stimulus PSTH and its individual rasters

CRF = get(gcf,'UserData');
cpanel = get(CRF.controls,'UserData');
dpanel = get(CRF.display,'UserData');

% Disenable multi stimulus comparison buttons
set(cpanel.uicontrols.multistimmenu,'enable','on')

% Hide previously chosen single CRF and stimuli
if isobject(dpanel.single.allstim)
    set(dpanel.single.allstim,'visible','off')
    set(dpanel.single.meanstim,'visible','off')
    set(dpanel.single.NKRfun,'visible','off')
    set(cpanel.single.selstim,'visible','off')
    set(cpanel.parvals.A,'string',[])
    set(cpanel.parvals.sig,'string',[])
    set(cpanel.parvals.exp,'string',[])
    set(cpanel.parvals.bl,'string',[])
    set(cpanel.parvals.kappa,'string',[])
end

% Show previously chosen mulit CRF and stimuli
L = ~isnan(dpanel.multi.meanstim);
try
    set(dpanel.multi.meanstim(L),'visible','on')
    set(dpanel.multi.NKRfun(L),'visible','on')
    set(cpanel.multi.stimsel(L,:),'visible','on')
catch
    disp('Error! Please fix.')
    keyboard
end

% Set the correct callback function
set(cpanel.axes.stim,'ButtonDownFcn',@CRFStimulusSelectionMulti)

% Save variables
set(CRF.controls,'UserData',cpanel);
set(CRF.display,'UserData',dpanel);
set(gcf,'UserData',CRF);

end


function CRFsubsel(~,~)
global GLMP
% Subunit selection

% Load Variables
CRF = get(gcf,'UserData');
cpanel = get(CRF.controls,'UserData');
dpanel = get(CRF.display,'UserData');

temp = cpanel.uicontrols.subunitmenu.Value;
subval = str2double(cpanel.uicontrols.subunitmenu.String{temp}(end));
if isnan(subval)
    subval = 3;
end

cols = {'ro' 'go' 'bo'};
axes(cpanel.axes.stimmap); cla;
cpanel.axes.stim = polar(GLMP.subunit{subval}.uniquetheta,...
    GLMP.subunit{subval}.uniquerho,cols{subval},'parent',cpanel.axes.stimmap);
cpanel.subselect = subval;

if cpanel.uicontrols.comparemenu.Value == 1
    set(cpanel.axes.stim,'ButtonDownFcn',@CRFStimulusSelectionSingle)
else
    set(cpanel.axes.stim,'ButtonDownFcn',@CRFStimulusSelectionMulti)
end

% Delete CRFs of previous subunit
axes(dpanel.CRF); cla;
dpanel.single.allstim = nan;
dpanel.single.meanstim = nan;
dpanel.single.NKRfun = nan;
cpanel.single.stimsel = nan(1,2);
dpanel.multi.meanstim = nan(4,1);
dpanel.multi.NKRfun = nan(4,1);
cpanel.multi.stimsel = nan(4,2);

% Save variables
set(CRF.controls,'UserData',cpanel)
set(CRF.display,'UserData',dpanel)
set(gcf,'UserData',CRF)

end


function CRFStimulusSelectionSingle(~,~)
global GLMP

% Grab current point before anything else
h = gca;
whichpt = get(h,'CurrentPoint');
whichpt = whichpt(1,[1 2]);

% Load variables
CRF = get(gcf,'UserData');
cpanel = get(CRF.controls,'UserData');
dpanel = get(CRF.display,'UserData');

% Grab other variables
%whichsub = cpanel.subselect;
temp = cpanel.uicontrols.subunitmenu.Value;
whichsub = str2double(cpanel.uicontrols.subunitmenu.String{temp}(end));
if isnan(whichsub)
    whichsub = 3;
end
[theta,~] = cart2pol(whichpt(1),whichpt(2));
[~,idx] = min(abs(GLMP.subunit{whichsub}.theta - theta));
idx = find(GLMP.subunit{whichsub}.theta==GLMP.subunit{whichsub}.theta(idx));
theta = GLMP.subunit{whichsub}.theta(idx(1));

% Hide compare stimuli
L = ~isnan(dpanel.multi.meanstim);
if any(L)
    set(dpanel.multi.meanstim(L),'visible','off')
    set(dpanel.multi.NKRfun(L),'visible','off')
    set(cpanel.multi.stimsel(L,:),'visible','off')
end

% plot smimulus selection
axes(cpanel.axes.stimmap)
cla; hold on;
cpanel.single.selstim(1) = polar([theta theta],[0 max(GLMP.subunit{whichsub}.rho)],'-k');
cpanel.single.selstim(2) = polar(GLMP.subunit{whichsub}.theta(idx),GLMP.subunit{whichsub}.rho(idx),'k*');
cols = [1 0 0; 0 1 0; 0 0 1; 1 0 1];
cpanel.axes.stim = polar(GLMP.subunit{whichsub}.uniquetheta,GLMP.subunit{whichsub}.uniquerho,...
    'o','parent',cpanel.axes.stimmap);
set(cpanel.axes.stim,'color',cols(whichsub,:))
set(cpanel.axes.stim,'ButtonDownFcn',@CRFStimulusSelectionSingle)
set(gca,'UserData',cpanel)

% Set up CRF plot
axes(dpanel.CRF);
cla; hold on; grid on;
jtramt = mean(diff(GLMP.subunit{whichsub}.rho(idx)));
dpanel.single.allstim = scatter(GLMP.subunit{whichsub}.rho(idx),GLMP.subunit{whichsub}.nspikes(idx),'ko',...
    'jitter','on','jitterAmount',jtramt);
scatter(zeros(size(GLMP.subunit{whichsub}.blnspikes)),GLMP.subunit{whichsub}.blnspikes,'ko',...
    'jitter','on','jitterAmount',jtramt);
idx2 = find(softEq(GLMP.subunit{whichsub}.uniquetheta,theta,1));
dpanel.single.meanstim = plot(GLMP.subunit{whichsub}.uniquerho(idx2),GLMP.subunit{whichsub}.meannspikes(idx2),'k*');
plot(0,mean(GLMP.subunit{whichsub}.blnspikes),'k*')
x = cat(1,0,GLMP.subunit{whichsub}.uniquerho(idx2));
y = cat(1,mean(GLMP.subunit{whichsub}.blnspikes),GLMP.subunit{whichsub}.meannspikes(idx2));
er = cat(1,std(GLMP.subunit{whichsub}.blnspikes),sqrt(GLMP.subunit{whichsub}.varnspikes(idx2)));
errorbar(x,y,er,'rx');
ylim([0 max(GLMP.subunit{whichsub}.nspikes)])
xlim([0 max(GLMP.subunit{whichsub}.rho(idx))])

%Fit a Naka-Rushton to the datapoints
Aguess = max(GLMP.subunit{whichsub}.nspikes(idx));
sigguess = max(GLMP.subunit{whichsub}.rho(idx))/2;
expguess = 3;
blguess = mean(GLMP.subunit{whichsub}.blnspikes);
kappaguess = 0;
paramsguess = [Aguess sigguess expguess blguess kappaguess];
blnsp = GLMP.subunit{whichsub}.blnspikes;
blcont = zeros(size(blnsp));
contrast = cat(1,blcont,GLMP.subunit{whichsub}.rho(idx));
response = cat(1,blnsp,GLMP.subunit{whichsub}.nspikes(idx));
FITSTR = 'symmetric';
vlb = [.001           .005  1         0         0];
vub = [ max(response)  50   5  max(response)    5];
options = optimset('Algorithm','interior-point','MaxFunEvals',5000,...
    'MaxIter',5000,'FinDiffType','central','Hessian','bfgs','display','off',...
    'TolFun',10.^-9);
[f,~] = fmincon('FitNakaRushtonFunJPW',paramsguess,[],[],[],[],vlb,vub,[],options,contrast,response,FITSTR,'NegativeBinomial');

% Display fit values
set(cpanel.parvals.A,'string',f(1),'foregroundcolor','k');
set(cpanel.parvals.sig,'string',f(2),'foregroundcolor','k');
set(cpanel.parvals.exp,'string',f(3),'foregroundcolor','k');
set(cpanel.parvals.bl,'string',f(4),'foregroundcolor','k');
set(cpanel.parvals.kappa,'string',f(5),'foregroundcolor','k');

% Plot NakaRushton
NRxvals = linspace(0,max(GLMP.subunit{whichsub}.rho(idx)),100);
NRyvals = ComputeNakaRushtonJPW(f,NRxvals,FITSTR);
dpanel.single.NKRfun = plot(NRxvals,NRyvals,'--m');

% Save variables
set(CRF.controls,'UserData',cpanel)
set(CRF.display,'UserData',dpanel)
set(gcf,'UserData',CRF)

end


function CRFStimulusSelectionMulti(~,~)
global GLMP

% Grab current point before anything else
h = gca;
whichpt = get(h,'CurrentPoint');
whichpt = whichpt(1,[1 2]);

% Load variables
CRF = get(gcf,'UserData');
cpanel = get(CRF.controls,'UserData');
dpanel = get(CRF.display,'UserData');

% Grab other variables
temp = cpanel.uicontrols.subunitmenu.Value;
whichsub = str2double(cpanel.uicontrols.subunitmenu.String{temp}(end));
if isnan(whichsub)
    whichsub = 3;
end
[theta,~] = cart2pol(whichpt(1),whichpt(2));
[~,idx] = min(abs(GLMP.subunit{whichsub}.theta - theta));
idx = find(GLMP.subunit{whichsub}.theta == GLMP.subunit{whichsub}.theta(idx));
theta = GLMP.subunit{whichsub}.theta(idx(1));
idx2 = find(softEq(GLMP.subunit{whichsub}.theta,theta,1));

% Determine which stimulus button is pressed
whichstim = cpanel.uicontrols.multistimmenu.Value;
cpanel.uicontrols.multistimmenu.UserData(whichstim) = theta;

% Hide previous stimuli, if any
if ~isnan(dpanel.multi.NKRfun(whichstim))
    set(cpanel.multi.stimsel(whichstim,:),'visible','off')
    set(dpanel.multi.meanstim(whichstim),'visible','off')
    set(dpanel.multi.NKRfun(whichstim),'visible','off')
end

% plot stimulus map
cols = [1 0 0; 0 1 0; 0 0 1; 1 0 1];
hold on;
cpanel.multi.stimsel(whichstim,1) = polar([theta theta],[0 max(GLMP.subunit{whichsub}.rho)],'-');
cpanel.multi.stimsel(whichstim,2) = polar(GLMP.subunit{whichsub}.theta(idx),GLMP.subunit{whichsub}.rho(idx),'*');
set(cpanel.multi.stimsel(whichstim,:),'Color',cols(whichstim,:))
set(cpanel.axes.stim,'ButtonDownFcn',@CRFStimulusSelectionMulti)
set(gca,'UserData',cpanel)


% Plot CRF
axes(dpanel.CRF);
hold on; grid on;
dpanel.multi.meanstim(whichstim) = plot(GLMP.subunit{whichsub}.rho(idx2),GLMP.subunit{whichsub}.nspikes(idx2),'*','Color',cols(whichstim,:));
xlim([0 max(GLMP.subunit{whichsub}.rho)])
ylim([0 max(GLMP.subunit{whichsub}.nspikes)])

%Fit a Naka-Rushton to the datapoints
Aguess = max(GLMP.subunit{whichsub}.nspikes(idx));
sigguess = max(GLMP.subunit{whichsub}.rho(idx))/2;
expguess = 3;
blguess = mean(GLMP.subunit{whichsub}.blnspikes);
kappaguess = 0;
paramsguess = [Aguess sigguess expguess blguess kappaguess];
blnsp = GLMP.subunit{whichsub}.blnspikes;
blcont = zeros(size(blnsp));
contrast = cat(1,blcont,GLMP.subunit{whichsub}.rho(idx));
response = cat(1,blnsp,GLMP.subunit{whichsub}.nspikes(idx));
FITSTR = 'symmetric';
vlb = [.001           .005  1      .001         0];
vub = [ max(response)  50   5  max(response)    5];
options = optimset('Algorithm','interior-point','MaxFunEvals',5000,...
    'MaxIter',5000,'FinDiffType','central','Hessian','bfgs','display','off',...
    'TolFun',10.^-9);
[f,~] = fmincon('FitNakaRushtonFunJPW',paramsguess,[],[],[],[],vlb,vub,[],options,contrast,response,FITSTR,'NegativeBinomial');


% Display fit values
set(cpanel.parvals.A,'string',f(1),'ForegroundColor',cols(whichstim,:));
set(cpanel.parvals.sig,'string',f(2),'ForegroundColor',cols(whichstim,:));
set(cpanel.parvals.exp,'string',f(3),'ForegroundColor',cols(whichstim,:));
set(cpanel.parvals.bl,'string',f(4),'ForegroundColor',cols(whichstim,:));
set(cpanel.parvals.kappa,'string',f(5),'ForegroundColor',cols(whichstim,:));

% Plot NakaRushton
NRxvals = linspace(0,max(GLMP.subunit{whichsub}.rho(idx)),100);
NRyvals = ComputeNakaRushtonJPW(f,NRxvals,FITSTR);
dpanel.multi.NKRfun(whichstim) = plot(NRxvals,NRyvals,'--','color',cols(whichstim,:));

% Save variables
set(CRF.controls,'UserData',cpanel)
set(CRF.display,'UserData',dpanel)
set(gcf,'UserData',CRF)

end

