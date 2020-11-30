function [] = GLMSGUI_CRF_highlum(~,~)
global GLMP


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  CRF Analysis Functions  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sets up figure for analysis of contrast response functions along a single axis.

figure(140); clf;
set(gcf,'numbertitle','off','name',['GLMP CRFs -high lum- (' GLMP.datafile ')'],'pos',[150 100 600 700])
CRF.controls = uipanel('parent',gcf,'units','normalized',...
    'pos',[.025 .65 .95 .325]);
cpanel = get(CRF.controls,'UserData');
cpanel.uicontrols.subbuttons = uibuttongroup('Parent',CRF.controls,...
    'units','normalized','pos',[.75 .525 .24 .45],...
    'SelectionChangeFcn',@CRFsubsel);
cpanel.uicontrols.sub1 = uicontrol('style','radiobutton',...
    'parent',cpanel.uicontrols.subbuttons,'pos',[15 65 100 25],...
    'string','Subunit #1','fontsize',12);
cpanel.uicontrols.sub2 = uicontrol('style','radiobutton',...
    'parent',cpanel.uicontrols.subbuttons,'pos',[15 37 100 25],...
    'string','Subunit #2','fontsize',12);
cpanel.uicontrols.sub3 = uicontrol('style','radiobutton',...
    'parent',cpanel.uicontrols.subbuttons,'pos',[15 10 100 25],...
    'string','Subunit #3','fontsize',12);
cpanel.subselect = 1;
if numel(GLMP.subunit) == 1
    set(cpanel.uicontrols.sub2,'enable','off')
    set(cpanel.uicontrols.sub3,'enable','off')
end

cpanel.axes.stimmap = axes('parent',CRF.controls,'units','pixels',...
    'pos',[20 10 200 200]);
cpanel.axes.stim = polar(GLMP.subunit{1}.specialcases.highlum.uniquetheta,GLMP.subunit{1}.specialcases.highlum.uniquerho,'ro',...
    'parent',cpanel.axes.stimmap);
set(cpanel.axes.stim,'ButtonDownFcn',@CRFStimulusSelectionSingle)

% Set up buttons to compare single stim to many stim
cpanel.uicontrols.comparebuttons = uibuttongroup('parent',CRF.controls,...
    'units','normalized','pos',[.45 .7 .28 .275],...
    'SelectionChangeFcn',@singleVmulti);
cpanel.uicontrols.singlestimbutton = uicontrol('style','radiobutton',...
    'parent',cpanel.uicontrols.comparebuttons,'pos',[15 30 150 25],...
    'string','Single Stimulus','fontsize',12);
cpanel.uicontrols.multistimbutton = uicontrol('style','radiobutton',...
    'parent',cpanel.uicontrols.comparebuttons,'pos',[15 7 150 25],...
    'string','Compare Stimuli','fontsize',12);


% For selecting which stimuli to compare
cpanel.uicontrols.multiStimSel = uibuttongroup('parent',CRF.controls,...
    'units','normalized','pos',[.45 .05 .28 .6]);
cpanel.uicontrols.stim1 = uicontrol('style','radiobutton',...
    'parent',cpanel.uicontrols.multiStimSel,'pos',[15 100 150 25],...
    'string','Stimulus #1','fontsize',12,'enable','off');
cpanel.uicontrols.stim2 = uicontrol('style','radiobutton',...
    'parent',cpanel.uicontrols.multiStimSel,'pos',[15 70 150 25],...
    'string','Stimulus #2','fontsize',12,'enable','off');
cpanel.uicontrols.stim3 = uicontrol('style','radiobutton',...
    'parent',cpanel.uicontrols.multiStimSel,'pos',[15 40 150 25],...
    'string','Stimulus #3','fontsize',12,'enable','off');
cpanel.uicontrols.stim4 = uicontrol('style','radiobutton',...
    'parent',cpanel.uicontrols.multiStimSel,'pos',[15 10 150 25],...
    'string','Stimulus #4','fontsize',12,'enable','off');


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


function singleVmulti(~,eventdata)
% Switch between the two options.

% Load variables
CRF = get(gcf,'UserData');
cpanel = get(CRF.controls,'UserData');
dpanel = get(CRF.display,'UserData');

% if isfield(dpanel,'axes') && ~isempty(dpanel.axes.NKRfun)
%     keyboard
%     %delete(dpanel.axes.allstim);
%     %delete(dpanel.axes.meanstim);
%     %delete(dpanel.axes.NKRfun);
%     %dpanel.axes.allstim = [];
%     %dpanel.axes.meanstim = [];
%     %dpanel.axes.NKRfun = [];
%     if isfield(cpanel.axes,'selstim')
%         delete(cpanel.axes.selstim);
%         cpanel.axes.selstim = [];
%     end
%     
%     % Save variables
%     set(CRF.display,'UserData',dpanel);
%     set(CRF.controls,'UserData',cpanel);
%     set(gcf,'UserData',CRF);
% end

% Reset displays
if eventdata.NewValue == cpanel.uicontrols.singlestimbutton
    resetCRFsingle();
elseif eventdata.NewValue == cpanel.uicontrols.multistimbutton
    resetCRFmulti();
else
    keyboard
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
set(cpanel.uicontrols.stim1,'enable','off')
set(cpanel.uicontrols.stim2,'enable','off')
set(cpanel.uicontrols.stim3,'enable','off')
set(cpanel.uicontrols.stim4,'enable','off')

% Show previously chosen single CRF and stimuli
if ~isnan(dpanel.single.meanstim)
    set(dpanel.single.allstim,'visible','on')
    set(dpanel.single.meanstim,'visible','on')
    set(dpanel.single.NKRfun,'visible','on')
    set(cpanel.single.selstim,'visible','on')
end

% Hide previously chosen multi CRF and stimuli
if any(~isnan(dpanel.multi.meanstim))
    L = ~isnan(dpanel.multi.meanstim);
    set(dpanel.multi.meanstim(L),'visible','off')
    set(dpanel.multi.NKRfun(L),'visible','off')
    set(cpanel.multi.stimsel(L,:),'visible','off')
end

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
set(cpanel.uicontrols.stim1,'enable','on')
set(cpanel.uicontrols.stim2,'enable','on')
set(cpanel.uicontrols.stim3,'enable','on')
set(cpanel.uicontrols.stim4,'enable','on')

% Hide previously chosen single CRF and stimuli
if ~isnan(dpanel.single.meanstim)
    set(dpanel.single.allstim,'visible','off')
    set(dpanel.single.meanstim,'visible','off')
    set(dpanel.single.NKRfun,'visible','off')
    set(cpanel.single.selstim,'visible','off')
end

% Show previously chosen mulit CRF and stimuli
if any(~isnan(dpanel.multi.meanstim))
    L = ~isnan(dpanel.multi.meanstim);
    set(dpanel.multi.meanstim(L),'visible','on')
    set(dpanel.multi.NKRfun(L),'visible','on')
    set(cpanel.multi.stimsel(L,:),'visible','on')
end

% Set the correct callback function
set(cpanel.axes.stim,'ButtonDownFcn',@CRFStimulusSelectionMulti)

% Save variables
set(CRF.controls,'UserData',cpanel);
set(CRF.display,'UserData',dpanel);
set(gcf,'UserData',CRF);

end


function CRFsubsel(~,eventdata)
global GLMP
% Subunit selection

% Load Variables
CRF = get(gcf,'UserData');
cpanel = get(CRF.controls,'UserData');
dpanel = get(CRF.display,'UserData');

subval = get(eventdata.NewValue,'string');
axes(cpanel.axes.stimmap); cla;
if strcmp(subval,'Subunit #1')
    cpanel.axes.stim = polar(GLMP.subunit{1}.specialcases.highlum.uniquetheta,GLMP.subunit{1}.specialcases.highlum.uniquerho,'ro',...
        'parent',cpanel.axes.stimmap);
    cpanel.subselect = 1;
elseif strcmp(subval,'Subunit #2')
    cpanel.axes.stim = polar(GLMP.subunit{2}.specialcases.highlum.uniquetheta,GLMP.subunit{2}.specialcases.highlum.uniquerho,'go',...
        'parent',cpanel.axes.stimmap);
    cpanel.subselect = 2;
elseif strcmp(subval,'Subunit #3')
    cpanel.axes.stim = polar(GLMP.subunit{3}.specialcases.highlum.uniquetheta,GLMP.subunit{3}.specialcases.highlum.uniquerho,'bo',...
        'parent',cpanel.axes.stimmap);
    cpanel.subselect = 3;
end
set(cpanel.axes.stim,'ButtonDownFcn',@CRFStimulusSelectionSingle)


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
whichsub = cpanel.subselect;
[theta,~] = cart2pol(whichpt(1),whichpt(2));
[~,idx] = min(abs(GLMP.subunit{whichsub}.specialcases.highlum.theta - theta));
idx = find(GLMP.subunit{whichsub}.specialcases.highlum.theta==GLMP.subunit{whichsub}.specialcases.highlum.theta(idx));
theta = GLMP.subunit{whichsub}.specialcases.highlum.theta(idx(1));

% plot figure
cla; hold on;
cpanel.single.selstim(1) = polar([theta theta],[0 max(GLMP.subunit{whichsub}.specialcases.highlum.rho)],'-k');
cpanel.single.selstim(2) = polar(GLMP.subunit{whichsub}.specialcases.highlum.theta(idx),GLMP.subunit{whichsub}.specialcases.highlum.rho(idx),'k*');
if whichsub == 1
    cpanel.axes.stim = polar(GLMP.subunit{1}.specialcases.highlum.uniquetheta,GLMP.subunit{1}.specialcases.highlum.uniquerho,'ro',...
        'parent',cpanel.axes.stimmap);
elseif whichsub == 2
    cpanel.axes.stim = polar(GLMP.subunit{2}.specialcases.highlum.uniquetheta,GLMP.subunit{2}.specialcases.highlum.uniquerho,'go',...
        'parent',cpanel.axes.stimmap);
elseif whichsub == 3.
    cpanel.axes.stim = polar(GLMP.subunit{3}.specialcases.highlum.uniquetheta,GLMP.subunit{3}.specialcases.highlum.uniquerho,'bo',...
        'parent',cpanel.axes.stimmap);
end
set(cpanel.axes.stim,'ButtonDownFcn',@CRFStimulusSelectionSingle)
set(gca,'UserData',cpanel)

axes(dpanel.CRF);
cla; hold on; grid on;
dpanel.single.allstim = plot(GLMP.subunit{whichsub}.specialcases.highlum.rho(idx),GLMP.subunit{whichsub}.specialcases.highlum.nspikes(idx),'ko');
idx2 = find(softEq(GLMP.subunit{whichsub}.specialcases.highlum.uniquetheta,theta,1));
dpanel.single.meanstim = plot(GLMP.subunit{whichsub}.specialcases.highlum.uniquerho(idx2),GLMP.subunit{whichsub}.specialcases.highlum.meannspikes(idx2),'k*');
ylim([0 max(GLMP.subunit{whichsub}.specialcases.highlum.nspikes)])
xlim([0 max(GLMP.subunit{whichsub}.specialcases.highlum.rho(idx))])

%Fit a Naka-Rushton to the datapoints
UB = max(GLMP.subunit{whichsub}.specialcases.highlum.nspikes(idx));
LB = min(GLMP.subunit{whichsub}.specialcases.highlum.nspikes(idx));
sig = (UB-LB)/max(GLMP.subunit{whichsub}.specialcases.highlum.rho(idx));
exp = 2;
paramsguess = [UB sig exp LB];
contrast = GLMP.subunit{whichsub}.specialcases.highlum.rho(idx);
response = GLMP.subunit{whichsub}.specialcases.highlum.nspikes(idx);
FITSTR = 'symmetric';
vlb = [0  0.0001 0  0];
vub = [1000 100 10 0];
options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
[f1,fval] = fmincon('FitNakaRushtonFunJPW',paramsguess,[],[],[],[],vlb,vub,[],options,contrast,response,FITSTR);

% Plot NakaRushton
NRxvals = linspace(0,max(GLMP.subunit{whichsub}.specialcases.highlum.rho(idx)),100);
NRyvals = ComputeNakaRushtonJPW(f1,NRxvals,FITSTR);
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
whichsub = cpanel.subselect;
[theta,~] = cart2pol(whichpt(1),whichpt(2));
[~,idx] = min(abs(GLMP.subunit{whichsub}.specialcases.highlum.theta - theta));
idx = find(GLMP.subunit{whichsub}.specialcases.highlum.theta == GLMP.subunit{whichsub}.specialcases.highlum.theta(idx));
theta = GLMP.subunit{whichsub}.specialcases.highlum.theta(idx(1));
idx2 = find(softEq(GLMP.subunit{whichsub}.specialcases.highlum.uniquetheta,theta,1));


% Determine which stimulus button is pressed
tempstim = get(cpanel.uicontrols.multiStimSel,'SelectedObject');
if tempstim == cpanel.uicontrols.stim1
    whichstim = 1;
    set(cpanel.uicontrols.stim1,'UserData',theta)
elseif tempstim == cpanel.uicontrols.stim2
    whichstim = 2;
    set(cpanel.uicontrols.stim2,'UserData',theta)
elseif tempstim == cpanel.uicontrols.stim3
    whichstim = 3;
    set(cpanel.uicontrols.stim3,'UserData',theta)
elseif tempstim == cpanel.uicontrols.stim4
    whichstim = 4;
    set(cpanel.uicontrols.stim4,'UserData',theta)
end

% Hide previous stimuli, if any
if ~isnan(dpanel.multi.NKRfun(whichstim))
    set(cpanel.multi.stimsel(whichstim,:),'visible','off')
    set(dpanel.multi.meanstim(whichstim),'visible','off')
    set(dpanel.multi.NKRfun(whichstim),'visible','off')
end

% plot stimulus map
cols = [1 0 0; 0 1 0; 0 0 1; 1 0 1];
hold on;
cpanel.multi.stimsel(whichstim,1) = polar([theta theta],[0 max(GLMP.subunit{whichsub}.specialcases.highlum.rho)],'-');
cpanel.multi.stimsel(whichstim,2) = polar(GLMP.subunit{whichsub}.specialcases.highlum.theta(idx),GLMP.subunit{whichsub}.specialcases.highlum.rho(idx),'*');
set(cpanel.multi.stimsel(whichstim,:),'Color',cols(whichstim,:))
set(cpanel.axes.stim,'ButtonDownFcn',@CRFStimulusSelectionMulti)
set(gca,'UserData',cpanel)


% Plot CRF
axes(dpanel.CRF);
hold on; grid on;

dpanel.multi.meanstim(whichstim) = plot(GLMP.subunit{whichsub}.specialcases.highlum.uniquerho(idx2),GLMP.subunit{whichsub}.specialcases.highlum.meannspikes(idx2),'*');
set(dpanel.multi.meanstim(whichstim),'Color',cols(whichstim,:))
%xlim([0 max(GLMP.subunit{whichsub}.specialcases.highlum.rho(idx))])
xlim([0 max(GLMP.subunit{whichsub}.specialcases.highlum.rho)])
ylim([0 max(GLMP.subunit{whichsub}.specialcases.highlum.nspikes)])

%Fit a Naka-Rushton to the datapoints
UB = max(GLMP.subunit{whichsub}.specialcases.highlum.nspikes(idx));
LB = min(GLMP.subunit{whichsub}.specialcases.highlum.nspikes(idx));
sig = (UB-LB)/max(GLMP.subunit{whichsub}.specialcases.highlum.rho(idx));
exp = 2;
paramsguess = [UB sig exp LB];
contrast = GLMP.subunit{whichsub}.specialcases.highlum.rho(idx);
response = GLMP.subunit{whichsub}.specialcases.highlum.nspikes(idx);
FITSTR = 'symmetric';
vlb = [0  0.0001 0  0];
vub = [1000 100 10 0];
options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
[f1,fval] = fmincon('FitNakaRushtonFunJPW',paramsguess,[],[],[],[],vlb,vub,[],options,contrast,response,FITSTR);

% Plot NakaRushton
%NRxvals = linspace(0,max(GLMP.subunit{whichsub}.specialcases.highlum.rho(idx)),100);
NRxvals = linspace(0,max(GLMP.subunit{whichsub}.specialcases.highlum.rho),100);
NRyvals = ComputeNakaRushtonJPW(f1,NRxvals,FITSTR);
dpanel.multi.NKRfun(whichstim) = plot(NRxvals,NRyvals,'--');
set(dpanel.multi.NKRfun(whichstim),'Color',cols(whichstim,:));

% Save variables
set(CRF.controls,'UserData',cpanel)
set(CRF.display,'UserData',dpanel)
set(gcf,'UserData',CRF)

end

