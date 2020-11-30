function GLMSGUI_GLMPTimeCourse
global GLMP


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  GLMP Time Course Analysis functions  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(41); clf;
set(gcf,'numbertitle','off','name',['GLMP Time Course of Responses (' GLMP.datafile ')'],'pos',[50 100 600 700])
StimTC = get(gcf,'UserData');
StimTC.controls = uipanel('parent',gcf,'units','normalized',...
    'pos',[.025 .65 .95 .325]);
cpanel = get(StimTC.controls,'UserData');

% Subunit Selection
cpanel.uicontrols.subbuttons = uibuttongroup('Parent',StimTC.controls,...
    'units','normalized','pos',[.75 .525 .24 .45],...
    'SelectionChangeFcn',@StimTCsubsel);
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

end



function StimTCsubsel(~,eventdata)
global GLMP
% Set up subunit selection

% Load variables
StimTC = get(gcf,'UserData');
cpanel = get(StimTC.controls,'UserData');
%spanel = get(PSTH.spikes,'UserData');

subval = get(eventdata.NewValue,'string');
axes(cpanel.axes.stimmap); cla;
if strcmp(subval,'Subunit #1')
    cpanel.axes.stim = polar(GLMP.subunit{1}.uniquetheta,GLMP.subunit{1}.uniquerho,'ro',...
        'parent',cpanel.axes.stimmap);
    cpanel.subselect = 1;
elseif strcmp(subval,'Subunit #2')
    cpanel.axes.stim = polar(GLMP.subunit{2}.uniquetheta,GLMP.subunit{2}.uniquerho,'go',...
        'parent',cpanel.axes.stimmap);
    cpanel.subselect = 2;
elseif strcmp(subval,'Subunit #3')
    cpanel.axes.stim = polar(GLMP.subunit{3}.uniquetheta,GLMP.subunit{3}.uniquerho,'bo',...
        'parent',cpanel.axes.stimmap);
    cpanel.subselect = 3;
end
set(cpanel.axes.stim,'ButtonDownFcn',@PSTHStimulusSelectionSingle)

% Save variables
set(PSTH.controls,'UserData',cpanel);
set(PSTH.spikes,'UserData',spanel);
set(gcf,'UserData',PSTH);

end



