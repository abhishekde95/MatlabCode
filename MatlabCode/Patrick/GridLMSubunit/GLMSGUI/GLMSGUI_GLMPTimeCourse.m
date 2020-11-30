function GLMSGUI_GLMPTimeCourse
global GLMP


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  GLMP Time Course Analysis functions  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(41); clf;
set(gcf,'numbertitle','off','name',['GLMP Time Course of Responses (' GLMP.datafile ')'],'pos',[50 100 800 500])
StimTC = get(gcf,'UserData');
StimTC.controls = uipanel('parent',gcf,'units','normalized',...
    'pos',[.025 .025 .35 .475]);
StimTC.stimmap = uipanel('parent',gcf,'units','normalized',...
    'pos',[.025 .525 .35 .45]);
StimTC.movie = uipanel('parent',gcf,'units','normalized',...
    'pos',[.4 .025 .575 .95]);

% Load variables
cpanel = get(StimTC.controls,'UserData');
mappanel = get(StimTC.movie,'UserData');
movpanel = get(StimTC.stimmap,'UserData');

% Subunit Selection
cpanel.uicontrols.subbuttons = uibuttongroup('Parent',StimTC.controls,...
    'units','normalized','pos',[.025 .525 .5 .45],...
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

% Draw stimulus map for stim #1 (as an initial condition)
mappanel.axes.map = axes('parent',StimTC.stimmap,...
    'units','normalized','pos',[.1 .1 .8 .8]);
markSize = (((GLMP.subunit{1}.meanfr) /(max(GLMP.subunit{1}.meanfr)))+.1)*50;
for t = 1:numel(GLMP.subunit{1}.uniqueLcc)
    plot(GLMP.subunit{1}.uniqueLcc(t),GLMP.subunit{1}.uniqueMcc(t),'ro','Markersize',markSize(t));
end
axis equal tight;
xlabel('Lcc')
ylabel('Mcc')


% Save variables
set(StimTC.controls,'UserData',cpanel);
set(StimTC.stimmap,'UserData',mappanel);
set(StimTC.movie,'UserData',movpanel);
set(gcf,'UserData',StimTC);

end


function StimTCsubsel(~,eventdata)
global GLMP
% Set up subunit selection

% Load variables
StimTC = get(gcf,'UserData');
cpanel = get(StimTC.controls,'UserData');
mappanel = get(StimTC.mappanel,'UserData');

subval = get(eventdata.NewValue,'string');
axes(mappanel.axes.map); cla;
if strcmp(subval,'Subunit #1')
    cpanel.subselect = 1;
    mfc = 'ro';
elseif strcmp(subval,'Subunit #2')
    cpanel.subselect = 2;
    mfc = 'go';
elseif strcmp(subval,'Subunit #3')
    cpanel.subselect = 3;
    mfc = 'bo';
end
%set(cpanel.axes.stim,'ButtonDownFcn',@StimTCStimulusSelection)

markSize = (((GLMP.subunit{cpanel.subselect}.meanfr) /(max(GLMP.subunit{cpanel.subselect}.meanfr)))+.1)*50;
grid on; hold on;
for t = 1:numel(GLMP.subunit{cpanel.subselect}.uniqueLcc)
    h = plot(GLMP.subunit{cpanel.subselect}.uniqueLcc(t),...
        GLMP.subunit{cpanel.subselect}.uniqueMcc(t),...
        mfc,'Markersize',markSize(t));
end
axis equal tight;
xlabel('Lcc')
ylabel('Mcc')

% Save variables
set(StimTC.controls,'UserData',cpanel);
set(StimTC.map,'UserData',mappanel);
set(gcf,'UserData',StimTC);

end



