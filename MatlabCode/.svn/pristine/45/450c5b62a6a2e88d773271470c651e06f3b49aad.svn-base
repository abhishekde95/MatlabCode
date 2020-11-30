function GLMSGUI_Special(~,~)
global GLMP
% This function deals with special cases of the GLMP dataset. The first
% special case analysis to be constructed is the high luminance case.


figure(70); clf;
set(gcf,'units','pixels','pos',[700 300 350 300],'NumberTitle','off',...
    'Name',['Special Cases Analysis (' GLMP.datafile ')']);
specialfig = get(gcf,'UserData');
specialfig.glmppanel = uipanel('Pos',[.025 .025 .95 .95],'Parent',gcf);

% Set up figure variables
glmppanel = get(specialfig.glmppanel,'UserData');


% Set up control panel
glmppanel.uicontrols.specialcases = uicontrol('style','popupmenu',...
    'parent',specialfig.glmppanel,'units','normalized','pos',[.1 .8 .8 .125],...
    'string',fieldnames(GLMP.subunit{1}.specialcases));
glmppanel.uicontrols.PSTH = uicontrol('Style','pushbutton','string','PSTH',...
    'parent',specialfig.glmppanel,'units','pixels','pos',[20 20 75 40],...
    'Callback',@GLMSGUI_PSTH_highlum);
glmppanel.uicontrols.CRF = uicontrol('Style','pushbutton','string','CRF',...
    'parent',specialfig.glmppanel,'units','pixels','pos',[20 80 75 40],...
    'Callback',@GLMSGUI_CRF_highlum);
glmppanel.uicontrols.surf = uicontrol('Style','pushbutton','string','Surface',...
    'parent',specialfig.glmppanel,'units','pixels','pos',[20 140 75 40],...
    'Callback',@GLMSGUI_Surface_highlum);

% conpanel.labels.specialcases = uicontrol('style','text','parent',specialfig.conpanel,...
%     'units','normalized','pos',[.1 .85 .8 .1],'fontsize',12,'string','Special Cases');
% conpanel.labels.specialcases = uicontrol('style','text','parent',specialfig.conpanel,...
%     'units','normalized','pos',[.05 .5 .9 .1],'string','Analyses','fontsize',12);
% conpanel.uicontrols.analyses = uibuttongroup('parent',specialfig.conpanel,...
%     'units','normalized','pos',[.05 .1 .9 .4],...
%     'selectionChangeFcn',@SpecialAnalSel);
% conpanel.analyses.CRF = uicontrol('style','radiobutton',...
%     'parent',conpanel.uicontrols.analyses,'pos',[15 60 100 25],...
%     'string','CRF','fontsize',12);
% conpanel.analyses.PSTH = uicontrol('style','radiobutton',...
%     'parent',conpanel.uicontrols.analyses,'pos',[15 35 100 25],...
%     'string','PSTH','fontsize',12);
% conpanel.analyses.surface = uicontrol('style','radiobutton',...
%     'parent',conpanel.uicontrols.analyses,'pos',[15 10 100 25],...
%     'string','Surface','fontsize',12);

% Save all settings
set(specialfig.glmppanel,'UserData',glmppanel)
set(gcf,'UserData',specialfig)



end

function SpecialAnalSel(~,b)
% This function sets up the disppanel and analpanel for the various analyses

% Get figure variables
specialfig = get(gcf,'UserData');
conpanel = get(specialfig.conpanel,'UserData');
disppanel = get(specialfig.disppanel,'UserData');
analpanel = get(specialfig.analpanel,'UserData');

if b.NewValue == conpanel.analyses.CRF
    SetupSpecialCRF()
elseif b.NewValue == conpanel.analyses.PSTH
    %SetupSpecialPSTH()
elseif b.NewValue == conpanel.analyses.surface
    %SetupSpecaialSurface()
else
    disp('Analysis not yet built.')
end

% Not saving variables here bc they are manipulated inside the 'reset'
% functions.  If saved again here, they would be saved as the
% previous, pre-reset variables.

end

function SetupSpecialCRF()
global GLMP

% Get figure variables
specialfig = get(gcf,'UserData');
conpanel = get(specialfig.conpanel,'UserData');
disppanel = get(specialfig.disppanel,'UserData');
analpanel = get(specialfig.analpanel,'UserData');

analpanel.uicontrols.subbuttons = uibuttongroup('Parent',specialfig.analpanel,...
    'units','normalized','pos',[.5 .6 .45 .375],...
    'SelectionChangeFcn',@Specialsubsel);
analpanel.uicontrols.sub1 = uicontrol('style','radiobutton',...
    'parent',analpanel.uicontrols.subbuttons,'pos',[15 65 100 25],...
    'string','Subunit #1','fontsize',12);
analpanel.uicontrols.sub2 = uicontrol('style','radiobutton',...
    'parent',analpanel.uicontrols.subbuttons,'pos',[15 37 100 25],...
    'string','Subunit #2','fontsize',12);
analpanel.uicontrols.sub3 = uicontrol('style','radiobutton',...
    'parent',analpanel.uicontrols.subbuttons,'pos',[15 10 100 25],...
    'string','Subunit #3','fontsize',12);
analpanel.subselect = 1;
if numel(GLMP.subunit) == 1
    set(analpanel.uicontrols.sub2,'enable','off')
    set(analpanel.uicontrols.sub3,'enable','off')
end

% Set up buttons to compare single stim to many stim
analpanel.uicontrols.comparebuttons = uibuttongroup('parent',specialfig.analpanel,...
    'units','normalized','pos',[.025 .7 .45 .275],...
    'SelectionChangeFcn',@singleVmulti);
analpanel.uicontrols.singlestimbutton = uicontrol('style','radiobutton',...
    'parent',analpanel.uicontrols.comparebuttons,'pos',[15 35 150 25],...
    'string','Single Stimulus','fontsize',12);
analpanel.uicontrols.multistimbutton = uicontrol('style','radiobutton',...
    'parent',analpanel.uicontrols.comparebuttons,'pos',[15 10 150 25],...
    'string','Compare Stimuli','fontsize',12);

% For selecting which stimuli to compare
analpanel.uicontrols.multiStimSel = uibuttongroup('parent',specialfig.analpanel,...
    'units','normalized','pos',[.025 .15 .45 .525]);
analpanel.uicontrols.stim1 = uicontrol('style','radiobutton',...
    'parent',analpanel.uicontrols.multiStimSel,'pos',[15 100 150 25],...
    'string','Stimulus #1','fontsize',12,'enable','off');
analpanel.uicontrols.stim2 = uicontrol('style','radiobutton',...
    'parent',analpanel.uicontrols.multiStimSel,'pos',[15 70 150 25],...
    'string','Stimulus #2','fontsize',12,'enable','off');
analpanel.uicontrols.stim3 = uicontrol('style','radiobutton',...
    'parent',analpanel.uicontrols.multiStimSel,'pos',[15 40 150 25],...
    'string','Stimulus #3','fontsize',12,'enable','off');
analpanel.uicontrols.stim4 = uicontrol('style','radiobutton',...
    'parent',analpanel.uicontrols.multiStimSel,'pos',[15 10 150 25],...
    'string','Stimulus #4','fontsize',12,'enable','off');

% Set up disppanel
disppanel.axes.stimmap = axes('parent',specialfig.disppanel,'units','pixels',...
    'pos',[20 100 200 200]);
alltheta = cat(1,GLMP.subunit{1}.uniquetheta,GLMP.subunit{1}.specialcases.highlum.uniquetheta);
allrho = cat(1,GLMP.subunit{1}.uniquerho,GLMP.subunit{1}.specialcases.highlum.uniquerho);
disppanel.axes.stim = polar(alltheta,allrho,'ro','parent',disppanel.axes.stimmap);
set(disppanel.axes.stim,'ButtonDownFcn',@CRFStimulusSelectionSingle)


% Save figure variables
set(specialfig.conpanel,'UserData',conpanel)
set(specialfig.analpanel,'UserData',analpanel)
set(specialfig.disppanel,'UserData',disppanel)
set(gcf,'UserData',specialfig)

end

function singleVmulti(~,eventdata)
% Switch between the two options.

% Load figure variables
specialfig = get(gcf,'UserData');
conpanel = get(specialfig.conpanel,'UserData');
disppanel = get(specialfig.disppanel,'UserData');
analpanel = get(specialfig.analpanel,'UserData');


% Reset displays
if eventdata.NewValue == analpanel.uicontrols.singlestimbutton
    resetCRFsingle();
elseif eventdata.NewValue == analpanel.uicontrols.multistimbutton
    resetCRFmulti();
else
    keyboard
end

% Not saving variables here bc they are manipulated inside the 'reset'
% functions.  If saved again here, they would be saved as the
% previous, pre-reset variables.

end

% cpanel.axes.stimmap = axes('parent',CRF.controls,'units','pixels',...
%     'pos',[20 10 200 200]);
% cpanel.axes.stim = polar(GLMP.subunit{1}.uniquetheta,GLMP.subunit{1}.uniquerho,'ro',...
%     'parent',cpanel.axes.stimmap);
% set(cpanel.axes.stim,'ButtonDownFcn',@CRFStimulusSelectionSingle)
% 
% 
% 
% 
% % Set up CRF Panel
% CRF.display = uipanel('parent',gcf,'units','normalized',...
%     'pos',[.025 .025 .95 .6]);
% dpanel = get(CRF.display,'UserData');
% dpanel.CRF = axes('parent',CRF.display,'units','normalized',...
%     'pos',[.05 .05 .9 .9],'Visible','on','box','on');
% dpanel.single.allstim = nan;
% dpanel.single.meanstim = nan;
% dpanel.single.NKRfun = nan;
% cpanel.single.stimsel = nan(1,2);
% dpanel.multi.meanstim = nan(4,1);
% dpanel.multi.NKRfun = nan(4,1);
% cpanel.multi.stimsel = nan(4,2);
% 
% % Save all settings
% set(CRF.controls,'UserData',cpanel)
% set(CRF.display,'UserData',dpanel)
% set(gcf,'UserData',CRF)
% 
