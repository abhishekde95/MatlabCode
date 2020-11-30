function [] = GLMSGUI_PSTH(~,~)
global GLMP
% This function performs a peri-stimulus timing histogram

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PSTH Analysis functions  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function sets up the PSTH analysis window.  
figure(32); clf;
set(gcf,'numbertitle','off','name',['GLMP PSTHs (' GLMP.datafile ')'],'pos',[50 100 600 700])
PSTH = get(gcf,'UserData');
PSTH.controls = uipanel('parent',gcf,'units','normalized',...
    'pos',[.025 .65 .95 .325]);
cpanel = get(PSTH.controls,'UserData');

% Subunit Selection
cpanel.uicontrols.subbuttons = uibuttongroup('Parent',PSTH.controls,...
    'units','normalized','pos',[.75 .525 .24 .45],...
    'SelectionChangeFcn',@PSTHsubsel);
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

% Stimulus Map
cpanel.axes.stimmap = axes('parent',PSTH.controls,'units','pixels',...
    'pos',[20 10 200 200]);
cpanel.axes.stim = polar(GLMP.subunit{1}.uniquetheta,GLMP.subunit{1}.uniquerho,'ro',...
    'parent',cpanel.axes.stimmap);
set(cpanel.axes.stim,'ButtonDownFcn',@PSTHStimulusSelectionSingle)

% Set up buttons to compare single stim to many stim
cpanel.uicontrols.comparebuttons = uibuttongroup('parent',PSTH.controls,...
    'units','normalized','pos',[.45 .7 .28 .275],...
    'SelectionChangeFcn',@singleVmulti);
cpanel.uicontrols.singlestimbutton = uicontrol('style','radiobutton',...
    'parent',cpanel.uicontrols.comparebuttons,'pos',[15 30 150 25],...
    'string','Single Stimulus','fontsize',12);
cpanel.uicontrols.multistimbutton = uicontrol('style','radiobutton',...
    'parent',cpanel.uicontrols.comparebuttons,'pos',[15 7 150 25],...
    'string','Compare Stimuli','fontsize',12);

% For selecting which stimuli to compare
cpanel.uicontrols.multiStimSel = uibuttongroup('parent',PSTH.controls,...
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

% For all PSTH's in order of acquisition
cpanel.uicontrols.alltrials = uicontrol('style','pushbutton',...
    'parent',PSTH.controls,'units','normalized','pos',[.75 .05 .24 .2],...
    'string','All PSTH"s','fontsize',12,'callback',@PlotAllRasters);

% Set up Spikes Panel and save User Data
PSTH.spikes = uipanel('parent',gcf,'units','normalized',...
    'pos',[.025 .025 .95 .6]);

% Save variables
set(PSTH.controls,'UserData',cpanel)
set(gcf,'UserData',PSTH)

end


function singleVmulti(~,eventdata)
% Switch between the two options.

% Load variables
PSTH = get(gcf,'UserData');
cpanel = get(PSTH.controls,'UserData');
spanel = get(PSTH.spikes,'UserData');

if isfield(spanel,'axes')
    delete(spanel.axes);
    spanel.axes = [];
    if isfield(cpanel.axes,'selstim')
        delete(cpanel.axes.selstim);
        cpanel.axes.selstim = [];
    end
    
    % Save variables
    set(PSTH.spikes,'UserData',spanel);
    set(PSTH.controls,'UserData',cpanel);
end

% Reset displays
if eventdata.NewValue == cpanel.uicontrols.singlestimbutton
    resetPSTHsingle();
elseif eventdata.NewValue == cpanel.uicontrols.multistimbutton
    resetPSTHmulti();
else
    keyboard
end

% Not saving variables here bc they are manipulated inside the 'reset'
% functions.  If saved again here, they would be saved as the
% previous, pre-reset variables.

end


function resetPSTHsingle()
global GLMP
% Reset the figure to display only 1 stimulus PSTH and its individual rasters

PSTH = get(gcf,'UserData');
cpanel = get(PSTH.controls,'UserData');
spanel = get(PSTH.spikes,'UserData');

% Disenable multi stimulus comparison buttons
set(cpanel.uicontrols.stim1,'enable','off')
set(cpanel.uicontrols.stim2,'enable','off')
set(cpanel.uicontrols.stim3,'enable','off')
set(cpanel.uicontrols.stim4,'enable','off')

% Set the correct callback function
set(cpanel.axes.stim,'ButtonDownFcn',@PSTHStimulusSelectionSingle)

% Determine which subunit is selected
if get(cpanel.uicontrols.subbuttons,'SelectedObject') == cpanel.uicontrols.sub1
    whichsub = 1;
elseif get(cpanel.uicontrols.subbuttons,'SelectedObject') == cpanel.uicontrols.sub2
    whichsub = 2;
elseif get(cpanel.uicontrols.subbuttons,'SelectedObject') == cpanel.uicontrols.sub3
    whichsub = 3;
end

% Draw panels (without data)
npanels = size(GLMP.subunit{whichsub}.spiketimes_col,2) + 1;
ypos0 = linspace(.975,.1,npanels+1);
axlims = [-.2 max(GLMP.subunit{whichsub}.stimDur)+.4];
for n = 1:numel(ypos0)-1
    if n == 1
        spanel.axes(n) = axes('parent',PSTH.spikes,'units','normalized',...
            'pos',[.025 ypos0(n+1) .95 ypos0(n)-ypos0(n+1)-.025],...
            'XTick',[],'XTickLabel',[],'XLim',axlims,'box','on');
    elseif n == numel(ypos0)-1
        spanel.axes(n) = axes('parent',PSTH.spikes,'units','normalized',...
            'pos',[.025 ypos0(n+1) .95 ypos0(n)-ypos0(n+1)-.025],...
            'YTick',[],'YTickLabel',[],...
            'XLim',axlims,'box','on');
        xlabel('Time from Stimulus Onset (ms)')
    else
        spanel.axes(n) = axes('parent',PSTH.spikes,'units','normalized',...
            'pos',[.025 ypos0(n+1) .95 ypos0(n)-ypos0(n+1)-.025],...
            'YTick',[],'YTickLabel',[],'XLim',axlims,...
            'XTick',[],'XTickLabel',[],'box','on');
    end
end

% Save variables
set(PSTH.controls,'UserData',cpanel);
set(PSTH.spikes,'UserData',spanel);
set(gcf,'UserData',PSTH);

end


function resetPSTHmulti
global GLMP
% Reset the figure to display many (4?) PSTHs

% Load variables
PSTH = get(gcf,'UserData');
cpanel = get(PSTH.controls,'UserData');
spanel = get(PSTH.spikes,'UserData');

% Turn on multi stimulus comparisson buttons
set(cpanel.uicontrols.stim1,'enable','on','UserData',[])
set(cpanel.uicontrols.stim2,'enable','on','UserData',[])
set(cpanel.uicontrols.stim3,'enable','on','UserData',[])
set(cpanel.uicontrols.stim4,'enable','on','UserData',[])

% Set the correct callback function
set(cpanel.axes.stim,'ButtonDownFcn',@PSTHStimulusSelectionMulti)

% Draw axes (without data)
npanels = 4;
whichsub = cpanel.subselect;
axlims = [-.2 max(GLMP.subunit{whichsub}.stimDur)+.4];
ypos0 = linspace(.975,.1,npanels+1);
for n = 1:npanels-1
    spanel.axes(n) = axes('parent',PSTH.spikes,'units','normalized',...
        'pos',[.025 ypos0(n+1) .95 ypos0(n)-ypos0(n+1)-.025],...
        'XTick',[],'XTickLabel',[],'box','on','XLim',axlims);
end
spanel.axes(npanels) = axes('parent',PSTH.spikes,'units','normalized',...
    'pos',[.025 ypos0(npanels+1) .95 ypos0(npanels)-ypos0(npanels+1)-.025],...
    'box','on','color',[1 .8 1],'XLim',axlims);
xlabel('Time from Stimulus Onset (ms)')
set(spanel.axes(1),'color',[1 .8 .8]);
set(spanel.axes(2),'color',[.8 1 .8]);
set(spanel.axes(3),'color',[.8 .8 1]);

% Save variables
set(PSTH.controls,'UserData',cpanel);
set(PSTH.spikes,'UserData',spanel);
set(gcf,'UserData',PSTH);

end


function PSTHsubsel(~,eventdata)
global GLMP
% Set up subunit selection

% Load variables
PSTH = get(gcf,'UserData');
cpanel = get(PSTH.controls,'UserData');
spanel = get(PSTH.spikes,'UserData');

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

singvmult = get(cpanel.uicontrols.comparebuttons.SelectedObject,'string');
if strcmp(singvmult,'Compare Stimuli')
    set(cpanel.axes.stim,'ButtonDownFcn',@PSTHStimulusSelectionMulti)
else
    set(cpanel.axes.stim,'ButtonDownFcn',@PSTHStimulusSelectionSingle)
end

% Clear histogram panels
for n = 1:numel(spanel.axes)
    axes(spanel.axes(n)); cla;
end
spanel.hist = [];

% Save variables
set(PSTH.controls,'UserData',cpanel);
set(PSTH.spikes,'UserData',spanel);
set(gcf,'UserData',PSTH);

end

function PSTHStimulusSelectionSingle(~,~)
global GLMP
% Selection of new stimulus as an '*' and erasing old stimulus selection.
% Also plots rasters. Consider making this a separate function?

h = gca;
whichpt = get(h,'CurrentPoint');
whichpt = whichpt(1,[1 2]);

% Load variables
PSTH = get(gcf,'UserData');
cpanel = get(PSTH.controls,'UserData');
%spanel = get(PSTH.spikes,'UserData');

whichsub = cpanel.subselect;
temp = [GLMP.subunit{whichsub}.uniqueLcc GLMP.subunit{whichsub}.uniqueMcc];
[~,idx] = min(sum(abs(temp - repmat(whichpt,size(temp,1),1)),2));
cla; hold on;
cpanel.axes.selstim = plot(temp(idx,1),temp(idx,2),'k*');
if whichsub == 1
    cpanel.axes.stim = polar(GLMP.subunit{1}.uniquetheta,GLMP.subunit{1}.uniquerho,'ro',...
        'parent',cpanel.axes.stimmap);
elseif whichsub == 2
    cpanel.axes.stim = polar(GLMP.subunit{2}.uniquetheta,GLMP.subunit{2}.uniquerho,'go',...
        'parent',cpanel.axes.stimmap);
elseif whichsub == 3
    cpanel.axes.stim = polar(GLMP.subunit{3}.uniquetheta,GLMP.subunit{3}.uniquerho,'bo',...
        'parent',cpanel.axes.stimmap);
end
set(cpanel.axes.stim,'ButtonDownFcn',@PSTHStimulusSelectionSingle)
set(gca,'UserData',cpanel)

% Set up PSTH Panel
PSTH = get(gcf,'UserData');
spanel = get(PSTH.spikes,'UserData');
axlims = [-.2 max(GLMP.subunit{whichsub}.stimDur)+.4];
if ~isfield(spanel,'axes')
    npanels = size(GLMP.subunit{whichsub}.spiketimes_col,2) + 1;
    ypos0 = linspace(.975,.1,npanels+1);
    for n = 1:numel(ypos0)-1
        if n == 1
            spanel.axes(n) = axes('parent',PSTH.spikes,'units','normalized',...
                'pos',[.025 ypos0(n+1) .95 ypos0(n)-ypos0(n+1)-.025],...
                'XTick',[],'XTickLabel',[],'box','on');
        elseif n == numel(ypos0)-1
            spanel.axes(n) = axes('parent',PSTH.spikes,'units','normalized',...
                'pos',[.025 ypos0(n+1) .95 ypos0(n)-ypos0(n+1)-.025],...
                'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'box','on');
        else
            spanel.axes(n) = axes('parent',PSTH.spikes,'units','normalized',...
                'pos',[.025 ypos0(n+1) .95 ypos0(n)-ypos0(n+1)-.025],...
                'YTick',[],'YTickLabel',[],'box','on');                
                % 'XTick',[],'XTickLabel',[], 
        end
        
    end
end

% Plot PSTH and individual rasters
for n = 1:numel(spanel.axes)
    axes(spanel.axes(n)); cla; hold on;
    set(spanel.axes(n),'Color',[1 1 1])
    if n == 1
        tpts = cat(1,GLMP.subunit{whichsub}.spiketimes_col{idx,:});
        if ~isempty(tpts)
            tpts(tpts<axlims(1)) = [];
            tpts(tpts>axlims(2)) = [];
            spanel.hist = hist(tpts,linspace(axlims(1),axlims(2),50));
            hist(tpts,linspace(axlims(1),axlims(2),50));
            plot([0 0],[0 max(spanel.hist)],'g-')
            stimDur = mean(GLMP.subunit{whichsub}.stimoff - GLMP.subunit{whichsub}.stimon);
            plot([stimDur stimDur],[0 max(spanel.hist)],'r-')
            plot([GLMP.countingwin(1) GLMP.countingwin(1)],[0 max(spanel.hist)],'--m')
            plot([GLMP.countingwin(2) GLMP.countingwin(2)],[0 max(spanel.hist)],'--m')
            set(spanel.axes(n),'XTick',[]);
        end
        set(spanel.axes(n),'XTickLabel',[],'xlim',axlims); hold on;
    else
        
        if n == numel(spanel.axes)
            set(spanel.axes(n),'xlim',axlims,'XtickLabelMode','auto','XTick',axlims(1):.1:axlims(2)); hold on;
            xlabel('Time From Stimulus Onset (s)')
        else
            set(spanel.axes(n),'YTickLabel',[],'XTickLabel',[],'xlim',axlims,...
                'XTick',[]); hold on;
        end
        tpts = GLMP.subunit{whichsub}.spiketimes_col{idx,n-1};
        if ~isempty(tpts)
            plot([tpts tpts]',repmat([-.5 .5],size(tpts,1),1)','k')
        end
        plot(axlims,[0 0],'k'); hold on;
        plot([0 0],[-.5 .5],'g-');
        stimDur = mean(GLMP.subunit{whichsub}.stimoff - GLMP.subunit{whichsub}.stimon);
        plot([stimDur stimDur],[-.5 .5],'r-')
        plot([GLMP.countingwin(1) GLMP.countingwin(1)],[-.5 .5],'--m')
        plot([GLMP.countingwin(2) GLMP.countingwin(2)],[-.5 .5],'--m')
        if size(tpts,2) == 0
            set(spanel.axes(n),'Color',[.5 .5 .5])
        end
    end
end

% Save variables
set(PSTH.controls,'UserData',cpanel);
set(PSTH.spikes,'UserData',spanel);
set(gcf,'UserData',PSTH);

end

function PSTHStimulusSelectionMulti(~,~)
global GLMP
% Marks current stimulus selection with an * and plots the PSTH

h = gca;
cla; hold on;
whichpt = get(h,'CurrentPoint');
whichpt = whichpt(1,[1 2]);

% Load variables
PSTH = get(gcf,'UserData');
cpanel = get(PSTH.controls,'UserData');
spanel = get(PSTH.spikes,'UserData');
spanel.hist = cell(4,1);

whichsub = cpanel.subselect;
temp = [GLMP.subunit{whichsub}.uniqueLcc GLMP.subunit{whichsub}.uniqueMcc];
[~,idx] = min(sum(abs(temp - repmat(whichpt,size(temp,1),1)),2));

% Determine which stimulus button is pressed
tempstim = get(cpanel.uicontrols.multiStimSel,'SelectedObject');
if tempstim == cpanel.uicontrols.stim1
    whichstim = 1;
    set(cpanel.uicontrols.stim1,'UserData',temp(idx,:))
elseif tempstim == cpanel.uicontrols.stim2
    whichstim = 2;
    set(cpanel.uicontrols.stim2,'UserData',temp(idx,:))
elseif tempstim == cpanel.uicontrols.stim3
    whichstim = 3;
    set(cpanel.uicontrols.stim3,'UserData',temp(idx,:))
elseif tempstim == cpanel.uicontrols.stim4
    whichstim = 4;
    set(cpanel.uicontrols.stim4,'UserData',temp(idx,:))
end

% Plot each of the stimulus selections
tempstim = get(cpanel.uicontrols.stim1,'UserData');
if ~isempty(tempstim)
    cpanel.axes.selstim(1) = plot(tempstim(1),tempstim(2),'r*');
end
tempstim = get(cpanel.uicontrols.stim2,'UserData');
if ~ isempty(tempstim)
    cpanel.axes.selstim(2) = plot(tempstim(1),tempstim(2),'g*');
end
tempstim = get(cpanel.uicontrols.stim3,'UserData');
if ~ isempty(tempstim)
    cpanel.axes.selstim(3) = plot(tempstim(1),tempstim(2),'b*');
end
tempstim = get(cpanel.uicontrols.stim4,'UserData');
if ~ isempty(tempstim)
    cpanel.axes.selstim(4) = plot(tempstim(1),tempstim(2),'m*');
end

% Plot all of the stimuli
if whichsub == 1
    cpanel.axes.stim = polar(GLMP.subunit{1}.uniquetheta,GLMP.subunit{1}.uniquerho,'ro',...
        'parent',cpanel.axes.stimmap);
elseif whichsub == 2
    cpanel.axes.stim = polar(GLMP.subunit{2}.uniquetheta,GLMP.subunit{2}.uniquerho,'go',...
        'parent',cpanel.axes.stimmap);
elseif whichsub == 3
    cpanel.axes.stim = polar(GLMP.subunit{3}.uniquetheta,GLMP.subunit{3}.uniquerho,'bo',...
        'parent',cpanel.axes.stimmap);
end

% Plot PSTH
PSTH = get(gcf,'UserData');
%spanel = get(PSTH.spikes,'UserData');
tpts = cat(1,GLMP.subunit{whichsub}.spiketimes_col{idx,:});
axes(spanel.axes(whichstim));
cla; hold on;
axlims = [-.2 max(GLMP.subunit{whichsub}.stimDur)+.4];
if ~isempty(tpts)
    tpts(tpts<axlims(1)) = [];
    tpts(tpts>axlims(2)) = [];
    spanel.hist{whichstim} = hist(tpts,linspace(axlims(1),axlims(2),50));
    hist(tpts,linspace(axlims(1),axlims(2),50));
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    maxhist = max(cat(2,spanel.hist{:}));
    plot([0 0],[0 maxhist],'g-')
    stimDur = mean(GLMP.subunit{whichsub}.stimoff - GLMP.subunit{whichsub}.stimon);
    plot([stimDur stimDur],[0 maxhist],'r-')
    set(spanel.axes(whichstim),'XTick',[]);
end
set(spanel.axes(whichstim),'XTickLabel',[],'xlim',axlims);
set(spanel.axes(end),'XtickLabelMode','auto','XTick',axlims(1):.1:axlims(2));
set(cpanel.axes.stim,'ButtonDownFcn',@PSTHStimulusSelectionMulti)

% Save variables
set(PSTH.controls,'UserData',cpanel);
set(PSTH.spikes,'UserData',spanel);
set(gcf,'UserData',PSTH);

end

function PlotAllRasters(~,~)
global GLMP

% Load figure variables
PSTH = get(gcf,'UserData');
cpanel = get(PSTH.controls,'UserData');

% Determine which subunit is selected
if get(cpanel.uicontrols.subbuttons,'SelectedObject') == cpanel.uicontrols.sub1
    whichsub = 1;
elseif get(cpanel.uicontrols.subbuttons,'SelectedObject') == cpanel.uicontrols.sub2
    whichsub = 2;
elseif get(cpanel.uicontrols.subbuttons,'SelectedObject') == cpanel.uicontrols.sub3
    whichsub = 3;
end

% Set up for rasters by trial
nrows = numel(GLMP.subunit{whichsub}.Lcc);
rowcoords = linspace(1,0,nrows*2+1);
rowcents = rowcoords(2:2:end);
rowspace = rowcoords(2)-rowcoords(1);
axlims = [-.2 max(GLMP.subunit{whichsub}.stimDur)+.4];

% Set up for psth
binsize = .01; % in seconds
bins = axlims(1):binsize:axlims(2);
PSTH = histc(cat(1,GLMP.subunit{whichsub}.normspiketimes{:}),bins);
PSTH = (PSTH./numel(GLMP.subunit{whichsub}.Lcc)) ./ binsize;

% Set up for BL psth
blfrs = GLMP.subunit{whichsub}.blfr;
blbins = 0:1:max(GLMP.subunit{whichsub}.blfr);
blfrpsth = histc(blfrs,blbins);

% Set up figure
rasterfig = figure(454); clf;
set(gcf,'units','normalized','pos',[.2 .1 .6 .8],'numbertitle','off');

% plot hist
stimpsthax = axes('parent',rasterfig,'units','normalized','pos',[.3 .85 .65 .1]);
psth = bar(bins,PSTH,'facecolor',[.5 0 .9],'edgecolor','k');
set(gca,'xlim',[min(bins) max(bins)])
title('PSTH')
xlabel('Time from Stim Onset (ms)')
ylabel('Mean Firing Rate (sp/s)')

% Plot rasters
rastersax = axes('parent',rasterfig,'units','normalized','pos',[.3 .05 .65 .725]); hold on; box on;
title('Rasters by Trial')
xlim(axlims)
ylim([0 1])
plot([GLMP.countingwin(1) GLMP.countingwin(1)],[0 1],'m-')
plot([GLMP.countingwin(2) GLMP.countingwin(2)],[0 1],'m-')
for n = 1:nrows
    tpts = GLMP.subunit{whichsub}.normspiketimes{n};
    if ~isempty(tpts)
        plot(repmat([tpts tpts],numel(tpts),1),repmat([rowcents(n)-rowspace rowcents(n)+rowspace],numel(tpts,1)),'k-')
    end
end
set(gca,'YTick',[0 .25 .5 .75 1],'YTickLabel',[nrows round(3*nrows/4) round(nrows/2) round(nrows/4) 1])
xlabel('Normalized Spike Time')
ylabel('Trial Number')

% Plot BL dist
bldistax = axes('parent',rasterfig,'units','normalized','pos',[.05 .6 .2 .35]); box on; hold on;
bar(blbins,blfrpsth,'facecolor',[.5 0 .9],'edgecolor','k');
if min(blbins) < max(blbins)
    set(gca,'xlim',[min(blbins) max(blbins)])
else
    set(gca,'xlim',[min(blbins) min(blbins)+1])
end
title('PSTH of Baseline FR')
xlabel('Baseline Firing Rate (sp/s)')
ylabel('Number of trials')

% Plot rank order
rankorderax = axes('parent',rasterfig,'units','normalized','pos',[.05 .05 .2 .45]); cla; box on; hold on;
trs = (1:numel(GLMP.subunit{whichsub}.blfr))';
[rho,pval] = corr(trs,blfrs,'type','Spearman')
plot(trs,blfrs,'ko')
ylabel('Firing Rate (sp/s)')
xlim([1 max(trs)])
if max(blfrs) > 0
    ylim([0 max(blfrs)])
else
    ylim([0 1])
end

disp('Rasters Computed')

end
