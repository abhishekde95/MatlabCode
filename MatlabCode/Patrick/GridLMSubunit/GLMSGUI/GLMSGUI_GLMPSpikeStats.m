function GLMSGUI_GLMPSpikeStats(varargin)
global GLMP

% This GUI should explore the timing and organization of spike selection.

% Set up figure
figure(7548); clf
set(gcf,'units','pixels','pos',[300 100 1000 700],'NumberTitle','off',...
    'Name',['Latency Analysis (' GLMP.datafile ')']);

% Set up panels
spikefig = get(gcf,'UserData');
spikefig.statspanel = uipanel('Pos',[.01 .36 .49 .63],'Parent',gcf,...
    'title','Statistics');
spikefig.conpanel = uipanel('pos',[.01 .01 .49 .34],'parent',gcf,...
    'title','Control Panel');
spikefig.histpanel = uipanel('pos',[.51 .01 .49 .98],'parent',gcf,...
    'title','Histograms');
spikefig.scatfig = [];

% Control Panel
% Set up color direction selector
cpanel.axes.stimmap = axes('parent',spikefig.conpanel,'units','pixels',...
    'pos',[20 10 200 200]);
cpanel.axes.stim = polar(GLMP.subunit{1}.uniquetheta,GLMP.subunit{1}.uniquerho,'ro',...
    'parent',cpanel.axes.stimmap);
set(cpanel.axes.stim,'ButtonDownFcn',@SelectColorDir)

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
cpanel.uicontrols.subunitmenu = uicontrol('parent',spikefig.conpanel,...
    'style','popupmenu','string',strs,'fontsize',12,...
    'units','normalized','pos',[.75 .8 .225 .1],...
    'Callback',@subsel);
if ~isempty(varargin) && ~ishandle(varargin{1})
    set(cpanel.uicontrols.subunitmenu,'Value',varargin{1})
end

% Set up binsize selector
cpanel.uicontrols.binsizeselector = uicontrol('parent',spikefig.conpanel,...
    'style','edit','string','.001','fontsize',12,...
    'units','normalized','pos',[.8 .5 .1 .1]);
cpanel.labels.binsizeselector = uicontrol('parent',spikefig.conpanel,...
    'style','text','string','Bin Size (in s)','fontsize',12,...
    'units','normalized','pos',[.75 .6 .2 .1]);

% Set up Gaussian filter selector
cpanel.uicontrols.gaussSizeselector = uicontrol('parent',spikefig.conpanel,...
    'style','edit','string','.01','fontsize',12,...
    'units','normalized','pos',[.8 .2 .1 .1]);
cpanel.labels.gaussSizeselector = uicontrol('parent',spikefig.conpanel,...
    'style','text','string','Gauss Filter Size (in s)','fontsize',12,...
    'units','normalized','pos',[.7 .3 .3 .1]);

% Scatterplot
cpanel.uicontrols.scatterplotbutton = uicontrol('parent',spikefig.conpanel,...
    'style','pushbutton','string','Scatterplot','fontsize',12,...
    'units','normalized','pos',[.5 .8 .15 .1],'Callback',@Scatterplot);
cpanel.degpm = 20;
cpanel.conthreshCol = .025;
cpanel.conthreshLum = .25;

% Stats Panel
% PSTH for all stimuli
spanel.axes.bigPSTH = axes('parent',spikefig.statspanel,'units','normalized',...
    'pos',[.1 .6 .375 .35]);

% Response vs Contrast
spanel.axes.response = axes('parent',spikefig.statspanel,...
    'units','normalized','pos',[.575 .6 .375 .35],'box','on');
title('Response vs. Cone Contrast')
xlabel('Cone Contrast')
ylabel('Response (mean sp/bin)')

% Latency vs Contrast
spanel.axes.latency = axes('parent',spikefig.statspanel,...
    'units','normalized','pos',[.575 .1 .375 .35],'box','on');
title('Latency vs. Cone Contrast')
xlabel('Cone Contrast')
ylabel('Latency (s)')

% Response vs Latency
spanel.axes.respvlat = axes('parent',spikefig.statspanel,...
    'units','normalized','pos',[.1 .1 .375 .35],'box','on');
title('Change in Latency and Response')
xlabel('Cone Contrast')
ylabel('Normalized to CC')

% Histogram Panel
% Set up panels for individual stimuli
naxes = 8;
yaxis = linspace(.975,.075,naxes+1);
yaxis = yaxis(2:end);
yheight = abs(mean(diff(yaxis))) * .8;
hpanel.labels.contrast = uicontrol('parent',spikefig.histpanel,...
    'units','normalized','style','text','pos',[.05 .95 .1 .025],...
    'string','Contrast','fontsize',12);
for n = 1:naxes
    hpanel.axes.contrast(n) = uicontrol('parent',spikefig.histpanel,...
        'units','normalized','style','edit','pos',[.05 yaxis(n)+.025 .1 .05]);        
    hpanel.axes.indivPSTH(n) = axes('parent',spikefig.histpanel,...
        'units','normalized','pos',[.25 yaxis(n) .3 yheight],'box','on');
    if n == naxes
        xlabel('Time (ms)')
    elseif n == 1
        title('PSTHs')
        set(gca,'xticklabel',[])
    else
        set(gca,'xticklabel',[])
    end
    hpanel.axes.smoothedPSTH(n) = axes('parent',spikefig.histpanel,...
        'units','normalized','pos',[.65 yaxis(n) .3 yheight],'box','on');
    if n == naxes
        xlabel('Time (ms)')
    elseif n == 1
        title('Smoothed PSTHs')
        set(gca,'xticklabel',[])
    else
        set(gca,'xticklabel',[])
    end
end

% Save figure variables
set(spikefig.conpanel,'userdata',cpanel)
set(spikefig.histpanel,'userdata',hpanel)
set(spikefig.statspanel,'userdata',spanel)
set(gcf,'userdata',spikefig)

BuildBaselineDist()
Scatterplot()


end

function SelectColorDir(~,~)
global GLMP

% Grab current point before anything else
h = gca;
whichpt = get(h,'CurrentPoint');
whichpt = whichpt(1,[1 2]);

% Baseline distribution
BuildBaselineDist()

% Load variables
spikefig = get(gcf,'UserData');
cpanel = get(spikefig.conpanel,'UserData');
hpanel = get(spikefig.histpanel,'UserData');
spanel = get(spikefig.statspanel,'userdata');

% Clear previous data
for n = 1:numel(hpanel.axes.smoothedPSTH)
    axes(hpanel.axes.smoothedPSTH(n)); cla;
    axes(hpanel.axes.indivPSTH(n)); cla;
end

% Grab other variables
temp = cpanel.uicontrols.subunitmenu.Value;
whichsub = str2double(cpanel.uicontrols.subunitmenu.String{temp}(end));
if isnan(whichsub)
    whichsub = 3;
end
[theta,~] = cart2pol(whichpt(1),whichpt(2));
[~,idx] = min(abs(GLMP.subunit{whichsub}.theta - theta));
idx = find(GLMP.subunit{whichsub}.uniquetheta==GLMP.subunit{whichsub}.theta(idx));
theta = GLMP.subunit{whichsub}.uniquetheta(idx(1));

% plot smimulus selection
axes(cpanel.axes.stimmap)
cla; hold on;
cpanel.single.selstim(1) = polar([theta theta],[0 max(GLMP.subunit{whichsub}.rho)],'-k');
cpanel.single.selstim(2) = polar(GLMP.subunit{whichsub}.uniquetheta(idx),GLMP.subunit{whichsub}.uniquerho(idx),'k*');
cols = [1 0 0; 0 1 0; 0 0 1; 1 0 1];
cpanel.axes.stim = polar(GLMP.subunit{whichsub}.uniquetheta,GLMP.subunit{whichsub}.uniquerho,...
    'o','parent',cpanel.axes.stimmap);
set(cpanel.axes.stim,'color',cols(whichsub,:))
set(cpanel.axes.stim,'ButtonDownFcn',@SelectColorDir)
set(gca,'UserData',cpanel)

% Populate histograms 
binsize = str2double(cpanel.uicontrols.binsizeselector.String); % in seconds
bins = -.2:binsize:.5;
[rhovals,i] = sort(GLMP.subunit{whichsub}.uniquerho(idx),'descend');
idx = idx(i);

if numel(rhovals) > 8
    rhovals = [rhovals(1:7); rhovals(end)];
end
for n = 1:numel(rhovals)
    spktimes = GLMP.subunit{whichsub}.spiketimes_cat(idx(n));
    hpanel.binnedspikes{n} = histc(spktimes{:},bins)./numel(GLMP.subunit{whichsub}.uniqueIdx{idx(n)});
    if any(hpanel.binnedspikes{n} > 1)
        keyboard
    end
    axes(hpanel.axes.indivPSTH(n)); cla; grid on; hold on;
    hpanel.psth{n} = bar(bins,hpanel.binnedspikes{n},'facecolor',[.5 0 .9],'edgecolor','k');
    xlim([min(bins) max(bins)*1.1]);
end

% Plot smoothed responses
gaussSize = str2double(cpanel.uicontrols.gaussSizeselector.String);
gaussSize = ceil(gaussSize/binsize);
gaussfilter = gausswin(gaussSize,3);
gaussfilter = gaussfilter./sum(gaussfilter); %normalizing
for n = 1:numel(rhovals)
    hpanel.filteredSpikes{n} = conv(hpanel.binnedspikes{n}',gaussfilter,'same');
    axes(hpanel.axes.smoothedPSTH(n)); cla; grid on; hold on;
    hpanel.smoothedpsth{n} = plot(bins,hpanel.filteredSpikes{n},'k-');
    nrpts = numel(GLMP.subunit{whichsub}.uniqueIdx{idx(n)});
    xlim([min(bins) max(bins)]);
    plot([min(bins) max(bins)],[hpanel.threshnspikes/sqrt(nrpts) hpanel.threshnspikes/sqrt(nrpts)],'r--')
end

% Find latency
hpanel.latency.firstbin = nan(numel(rhovals),1);
hpanel.response.firstbin.meannspikes = nan(numel(rhovals),1);
hpanel.response.firstbin.maxnspikes = nan(numel(rhovals),1);
hpanel.response.firstbin.sumnspikes = nan(numel(rhovals),1);
tbin1 = find(bins >= .04,1);
tbin2 = find(bins >= .2,1);
%lastlatbin = find(bins==0); % for monotonic conditional

% Latency by first bin to exceed threshold
for n = 1:numel(rhovals)
    nrpts = numel(GLMP.subunit{whichsub}.uniqueIdx{idx(n)});
    temp = find(hpanel.filteredSpikes{n} > hpanel.threshnspikes/sqrt(nrpts));
    %temp = temp(temp>=lastlatbin); % monotonic conditional
    temp = temp(temp>=tbin1 & temp<=tbin2); % specifying latency must be within some window
    set(hpanel.axes.contrast(n),'string',rhovals(n))
    if ~isempty(temp)
        latencybin = temp(1);
        %lastlatbin = latencybin; % monotonic conditional
        %lastlatbin = find(bins==0); % just sometime after t=0
        hpanel.latency.firstbin(n) = bins(latencybin);
        axes(hpanel.axes.smoothedPSTH(n));
        plot([bins(latencybin) bins(latencybin)],[0 max(hpanel.filteredSpikes{n})],'r-')
        
        % Find mean number of spikes for 200 ms after latency
        L = bins>=bins(latencybin) & bins < (bins(latencybin)+.2);
        hpanel.response.firstbin.meannspikes(n) = mean(hpanel.filteredSpikes{n}(L));
        hpanel.response.firstbin.maxnspikes(n) = max(hpanel.filteredSpikes{n}(L));
        hpanel.response.firstbin.sumnspikes(n) = sum(hpanel.filteredSpikes{n}(L));
    end
end
tics = get(hpanel.axes.smoothedPSTH(1),'xtick');
lims = get(hpanel.axes.smoothedPSTH(1),'xlim');
set(hpanel.axes.smoothedPSTH(end),'xtick',tics,'xlim',lims); grid on;


% Latency by specifying a temporal form
hpanel.latency.gausswin = nan(numel(rhovals),1);
hpanel.response.gausswin.sumnspikes = nan(numel(rhovals),1);
hpanel.response.gausswin.meannspikes = nan(numel(rhovals),1);
hpanel.response.gausswin.maxnspikes = nan(numel(rhovals),1);
whichbins = bins >= 0 & bins <= .2;
gaussfilter = gausswin(sum(whichbins),2); % temportal form is a gaussian
tbin1 = find(bins >=  0,1);
tbin2 = find(bins >= .3,1);
for n = 1:numel(rhovals)
    
    % latency
    gf = gaussfilter .* max(hpanel.filteredSpikes{n});
    hpanel.tempwindowspikes{n} = conv(hpanel.filteredSpikes{n}(tbin1:tbin2)',gf,'valid');
    [m,i] = max(hpanel.tempwindowspikes{n});
    i = i+tbin1;
    if m ~= 0 
        hpanel.latency.gausswin(n) = bins(i);
    end
    
    % response
    if m ~= 0
        L = false(size(bins));
        L(i:i+200) = true;
        hpanel.response.gausswin.sumnspikes(n) = sum(hpanel.filteredSpikes{n}(L));
        hpanel.response.gausswin.meannspikes(n) = mean(hpanel.filteredSpikes{n}(L));
        hpanel.response.gausswin.maxnspikes(n) = max(hpanel.filteredSpikes{n}(L));
        
        % Plot
        axes(hpanel.axes.smoothedPSTH(n));
        plot(bins(L),gf,'m-')
    end
end

% Plot latency
axes(spanel.axes.latency); cla; hold on; grid on;
plot(rhovals,hpanel.latency.firstbin,'r*--')
plot(rhovals,hpanel.latency.gausswin,'m*--')
xlabel('Cone Contrast')
ylabel('Latency (s)')

% Plot normalized cc vs latency and cc vs mean spikes
% current idea: normalize so that latency goes between min cc and max cc
% (to compare with linear cc)
latency = hpanel.latency.firstbin - min(hpanel.latency.firstbin); % reduce min response to 0
normlatency = latency./max(latency); % normalize max to 1
normlatency = 1 - normlatency; % flip
scaledlatency = normlatency * (max(rhovals)-min(rhovals)) + min(rhovals); % normalize to contrasts
% same for latency measured with Gaussian temporal window
latency =  hpanel.latency.gausswin - min(hpanel.latency.gausswin);
normlatency = latency./max(latency);
normlatency = 1 - normlatency;
scaledlatency2 = normlatency * (max(rhovals)-min(rhovals)) + min(rhovals);
% same for spikes
tempnsp = hpanel.response.firstbin.meannspikes;
nsp = tempnsp - min(tempnsp);
normnsp = nsp./max(nsp);
scalednsp = normnsp * (max(rhovals)-min(rhovals)) + min(rhovals);
% same for spikes in gaussian window
tempnsp2 = hpanel.response.gausswin.meannspikes;
nsp = tempnsp2 - min(tempnsp2);
normnsp = nsp./max(nsp);
scalednsp2 = normnsp * (max(rhovals)-min(rhovals)) + min(rhovals);

% latency vs response
axes(spanel.axes.respvlat); cla; hold on; grid on;
plot(rhovals,scaledlatency,'r*--')
plot(rhovals,scaledlatency2,'m*--')
plot(rhovals,scalednsp,'*g--')
plot(rhovals,scalednsp2,'*c--')
legend('Latency','Response','Location','SouthEast')
plot(rhovals,rhovals,'k--')

% Plot response vs contrast
axes(spanel.axes.response); cla; hold on; grid on;
plot(rhovals,tempnsp,'g*--')
plot(rhovals,tempnsp2,'c*--')
xlabel('Cone Contrast')
ylabel('Mean Smoothed Spikes/Bin')
title('Contrast vs Response')

% Save variables
set(spikefig.conpanel,'UserData',cpanel)
set(spikefig.histpanel,'UserData',hpanel)
set(spikefig.statspanel,'userdata',spanel)
set(gcf,'UserData',spikefig)

end

function BuildBaselineDist()
global GLMP

% Load variables
spikefig = get(gcf,'UserData');
cpanel = get(spikefig.conpanel,'UserData');
hpanel = get(spikefig.histpanel,'UserData');
spanel = get(spikefig.statspanel,'userdata');

% Which sub
% Grab other variables
temp = cpanel.uicontrols.subunitmenu.Value;
whichsub = str2double(cpanel.uicontrols.subunitmenu.String{temp}(end));
if isnan(whichsub)
    whichsub = 3;
end

% PSTH 
axes(spanel.axes.bigPSTH); cla; hold on;
binsize = str2double(cpanel.uicontrols.binsizeselector.String); % in seconds
psthbins = -.2:binsize:.5;
PSTH = histc(cat(1,GLMP.subunit{whichsub}.normspiketimes{:}),psthbins);
PSTH = (PSTH./numel(GLMP.subunit{whichsub}.normspiketimes));
psthpanel.psth = bar(psthbins,PSTH,'facecolor',[.5 0 .9],'edgecolor','k');
set(gca,'xlim',[min(psthbins) max(psthbins)]); grid on;
title('PSTH')
xlabel('Time from Stim Onset (ms)')
ylabel('Mean Spikes/Bin')

% Bin and smooth baseline 
binsize = str2double(cpanel.uicontrols.binsizeselector.String); % in seconds
blbins = -.2:binsize:0;
gaussSize = str2double(cpanel.uicontrols.gaussSizeselector.String);
gaussSize = ceil(gaussSize/binsize);
gaussfilter = gausswin(gaussSize,3);
gaussfilter = gaussfilter./sum(gaussfilter); %normalizing
for n = 1:numel(GLMP.normspiketimes)
    hpanel.blpsth{n} = shiftdim(histc(GLMP.normspiketimes{n},blbins))';
end
hpanel.meanblbins = mean(cat(1,hpanel.blpsth{:}),1);
hpanel.smoothblbins = conv(hpanel.meanblbins,gaussfilter,'same');

% Define threshold for signal vs baseline noise
threshpercent = 99;
hpanel.threshnspikes = prctile(hpanel.smoothblbins,threshpercent);
plot([psthbins(1) psthbins(end)],[hpanel.threshnspikes hpanel.threshnspikes],'r--')
%nrepeats = floor(numel(GLMP.subunit{whichsub}.Mcc)/numel(GLMP.subunit{whichsub}.uniqueMcc));
hpanel.threshnspikes = (hpanel.threshnspikes * sqrt(numel(GLMP.subunit{whichsub}.rho)));% / sqrt(nrepeats));

% Save variables
set(spikefig.conpanel,'UserData',cpanel)
set(spikefig.histpanel,'UserData',hpanel)
set(spikefig.statspanel,'userdata',spanel)
set(gcf,'UserData',spikefig)

end

function subsel(~,~)
global GLMP

% Load variables
spikefig = get(gcf,'UserData');
cpanel = get(spikefig.conpanel,'UserData');
hpanel = get(spikefig.histpanel,'UserData');
spanel = get(spikefig.statspanel,'userdata');

% Find selected subunit and map stimuli
temp = cpanel.uicontrols.subunitmenu.Value;
subval = str2double(cpanel.uicontrols.subunitmenu.String{temp}(end));
if isnan(subval)
    subval = 3;
end
cols = {'ro' 'go' 'bo'};
axes(cpanel.axes.stimmap); cla;
cpanel.axes.stim = polar(GLMP.subunit{subval}.uniquetheta,...
    GLMP.subunit{subval}.uniquerho,cols{subval},'parent',cpanel.axes.stimmap);
set(cpanel.axes.stim,'ButtonDownFcn',@SelectColorDir)

% Clear the hist panels
for n = 1:numel(hpanel.axes.indivPSTH)
    set(hpanel.axes.contrast(n),'string',[])
    axes(hpanel.axes.indivPSTH(n)); cla;
    axes(hpanel.axes.smoothedPSTH(n)); cla;
end
axes(spanel.axes.response); cla;
axes(spanel.axes.latency); cla;
axes(spanel.axes.respvlat); cla;
BuildBaselineDist();

% Save variables
set(spikefig.conpanel,'UserData',cpanel)
set(spikefig.histpanel,'UserData',hpanel)
set(spikefig.statspanel,'userdata',spanel)
set(gcf,'UserData',spikefig)

end

function Scatterplot(~,~)
global GLMP

% Load variables
spikefig = get(gcf,'UserData');
cpanel = get(spikefig.conpanel,'UserData');
hpanel = get(spikefig.histpanel,'UserData');

% Find selected subunit
temp = cpanel.uicontrols.subunitmenu.Value;
whichsub = str2double(cpanel.uicontrols.subunitmenu.String{temp}(end));
if isnan(whichsub)
    whichsub = 3;
end

% Set up variables
scatfig.firstsp.latency = nan(numel(GLMP.subunit{whichsub}.uniquerho),1);
scatfig.firstsp.response = nan(numel(GLMP.subunit{whichsub}.uniquerho),1);
scatfig.gauss.latency = nan(numel(GLMP.subunit{whichsub}.uniquerho),1);
scatfig.gauss.response = nan(numel(GLMP.subunit{whichsub}.uniquerho),1);
scatfig.contrast = GLMP.subunit{whichsub}.uniquerho;

% Latency + response by first spike in some time window
binsize = str2double(cpanel.uicontrols.binsizeselector.String); % in seconds
bins = -.2:binsize:.5;
gaussSize = str2double(cpanel.uicontrols.gaussSizeselector.String);
gaussSize = ceil(gaussSize/binsize);
gaussfilter = gausswin(gaussSize,3);
gaussfilter = gaussfilter./sum(gaussfilter); %normalizing
tbin1 = find(bins >= .04,1);
tbin2 = find(bins >= .2,1);
for n = 1:numel(GLMP.subunit{whichsub}.uniquerho)
    nrpts = numel(GLMP.subunit{whichsub}.uniqueIdx{n});
    spktimes = GLMP.subunit{whichsub}.spiketimes_cat{n};
    scatfig.binnedspikes{n} = histc(spktimes,bins)./numel(GLMP.subunit{whichsub}.uniqueIdx{n});
    scatfig.filteredSpikes{n} = conv(scatfig.binnedspikes{n}',gaussfilter,'same');
    temp = find(scatfig.filteredSpikes{n} > (hpanel.threshnspikes/sqrt(nrpts)));
    temp = temp(temp>=tbin1 & temp<=tbin2); % specifying latency must be within some window
    if ~isempty(temp)
        latencybin = temp(1); % the first spike determines latency
        scatfig.firstsp.latency(n) = bins(latencybin);
        L = bins >= bins(latencybin) & bins <= (bins(latencybin)+.2);
        scatfig.firstsp.response(n) = mean(scatfig.filteredSpikes{n}(L));
    end
end

% Latency and response by specifying a temporal form
rhovals = GLMP.subunit{whichsub}.uniquerho;
binsize = str2double(cpanel.uicontrols.binsizeselector.String); % in seconds
bins = -.2:binsize:.5;
whichbins = bins >= 0 & bins <= .2;
gaussfilter = gausswin(sum(whichbins),3); % temportal form is a gaussian
tbin1 = find(bins >=  0,1);
tbin2 = find(bins >= .4,1);
for n = 1:numel(rhovals)
    
    % Filter spikes
    %spktimes = GLMP.subunit{whichsub}.spiketimes_cat{n};
    %scatfig.binnedspikes{n} = histc(spktimes,bins)./numel(GLMP.subunit{whichsub}.uniqueIdx{n});
    %scatfig.filteredSpikes{n} = conv(scatfig.binnedspikes{n}',gaussfilter,'same');
    
    % latency
    gf = gaussfilter .* max(scatfig.filteredSpikes{n});
    scatfig.tempwindowspikes{n} = conv(scatfig.filteredSpikes{n}(tbin1:tbin2)',gf,'valid');
    [m,i] = max(scatfig.tempwindowspikes{n});
    i = i+tbin1;
    if m ~= 0 
        scatfig.gauss.latency(n) = bins(i);
    end
    
    % response
    if m ~= 0
        L = false(size(bins));
        L(i:i+200) = true;
        scatfig.gauss.response(n) = mean(scatfig.filteredSpikes{n}(L));
    end
end

% Plot lum and col points
%thetas = GLMP.subunit{whichsub}.uniquetheta;
%lumL = (thetas > (pi/4 - degpm/180*pi) & thetas < (pi/4 + degpm/180*pi))...
%    | (thetas > (-3*pi/4 - degpm/180*pi) & thetas < (-3*pi/4 + degpm/180*pi));
%colL = (thetas > (3*pi/4 - degpm/180*pi) & thetas < (3*pi/4 + degpm/180*pi))...
%    | (thetas > (-pi/4 - degpm/180*pi) & thetas < (-pi/4 + degpm/180*pi));
scatfig.degpm = cpanel.degpm;

% Insead of using just certain spokes, project onto col and lum axes
[colx,coly] = pol2cart(-pi/4,1);
[lumx,lumy] = pol2cart(pi/4,1);
scatfig.colprojcc =  [GLMP.subunit{whichsub}.uniqueLcc GLMP.subunit{whichsub}.uniqueMcc] * [colx coly]';
scatfig.lumprojcc =  [GLMP.subunit{whichsub}.uniqueLcc GLMP.subunit{whichsub}.uniqueMcc] * [lumx lumy]';
contcolL = abs(scatfig.colprojcc) > cpanel.conthreshCol;
contlumL = abs(scatfig.lumprojcc) > cpanel.conthreshLum;
%spthreshL = scatfig.response > prctile(cat(1,scatfig.firstsp.response,scatfig.gauss.response),20);
%colL = contcolL & spthreshL;
%lumL = contlumL & spthreshL;

% Set up figure
figure(3765); clf;
set(gcf,'units','normalized','pos',[.3 .3 .6 .6])
scatfig.axes.resplat = axes('parent',gcf,'units','normalized',...
    'pos',[.05 .55 .9 .4],'box','on'); grid on; hold on;
xlabel('Response (mean sp/bin)')
ylabel('Latency (s)')
title('Response vs Latency')
scatfig.axes.contlat = axes('parent',gcf,'units','normalized',...
    'pos',[.05 .05 .9 .4],'box','on'); grid on; hold on;
xlabel('Contrast')
ylabel('Latency (s)')
title('Contrast vs Latency')

% Plot Latency vs Response
axes(scatfig.axes.resplat); cla; hold on;
p = nan(4,1);

% First Spike: Luminance
spthreshL = scatfig.firstsp.response > prctile(scatfig.firstsp.response,20);
lumL = contlumL & spthreshL;
if sum(lumL) > 0
    x = cat(2,ones(sum(lumL),1),scatfig.firstsp.response(lumL));
    y = scatfig.firstsp.latency(lumL);
    [b,bint,~,~,stats] = regress(y,x);
    scatfig.regress.firstsp.lum.coefs = b;
    scatfig.regress.firstsp.lum.coefsint = bint;
    scatfig.regress.firstsp.lum.pval = stats(3);
    xvals = linspace(0,max(x(:,2)),50);
    yvals = b(2).*xvals + b(1);
    p(1) = plot(xvals,yvals,'--','color',[1 0 1]);
    plot(x(:,2),y,'o','color',[1 0 1])
else
    scatfig.regress.firstsp.lum.coefs = nan(2,1);
    scatfig.regress.firstsp.lum.coefsint = nan(2,2);
    scatfig.regress.firstsp.lum.pval = nan;
end

% First Spike: Chromatic
colL = contcolL & spthreshL;
if sum(colL) > 0
    x = cat(2,ones(sum(colL),1),scatfig.firstsp.response(colL));
    y = scatfig.firstsp.latency(colL);
    [b,bint,~,~,stats] = regress(y,x);
    scatfig.regress.firstsp.col.coefs = b;
    scatfig.regress.firstsp.col.coefsint = bint;
    scatfig.regress.firstsp.col.pval = stats(3);
    xvals = linspace(0,max(x(:,2)),50);
    yvals = b(2).*xvals + b(1);
    p(3) = plot(xvals,yvals,'--','color',[0 1 0]);
    plot(x(:,2),y,'o','color',[0 1 0])
else
    scatfig.regress.firstsp.col.coefs = nan(2,1);
    scatfig.regress.firstsp.col.coefsint = nan(2,2);
    scatfig.regress.firstsp.col.pval = nan;
end

% Gaussian: Luminance
spthreshL = scatfig.gauss.response > prctile(scatfig.gauss.response,20);
lumL = contlumL & spthreshL;
if sum(lumL) > 0
    x = cat(2,ones(sum(lumL),1),scatfig.gauss.response(lumL));
    y = scatfig.gauss.latency(lumL);
    [b,bint,~,~,stats] = regress(y,x);
    scatfig.regress.gauss.lum.coefs = b;
    scatfig.regress.gauss.lum.coefsint = bint;
    scatfig.regress.gauss.lum.pval = stats(3);
    xvals = linspace(0,max(x(:,2)),50);
    yvals = b(2).*xvals + b(1);
    p(2) = plot(xvals,yvals,'*-','color',[1 0 1]);
    plot(x(:,2),y,'*','color',[1 0 1])
else
    scatfig.regress.gauss.lum.coefs = nan(2,1);
    scatfig.regress.gauss.lum.coefsint = nan(2,2);
    scatfig.regress.gauss.lum.pval = nan;
end

% Gaussian: Chromatic
colL = contcolL & spthreshL;
if sum(colL) > 0
    x = cat(2,ones(sum(colL),1),scatfig.gauss.response(colL));
    y = scatfig.gauss.latency(colL);
    [b,bint,~,~,stats] = regress(y,x);
    scatfig.regress.gauss.col.coefs = b;
    scatfig.regress.gauss.col.coefsint = bint;
    scatfig.regress.gauss.col.pval = stats(3);
    xvals = linspace(0,max(x(:,2)),50);
    yvals = b(2).*xvals + b(1);
    p(4) = plot(xvals,yvals,'*-','color',[.2 1 0]);
    plot(x(:,2),y,'*','color',[0 1 0])
else
    scatfig.regress.gauss.col.coefs = nan(2,1);
    scatfig.regress.gauss.col.coefsint = nan(2,2);
    scatfig.regress.gauss.col.pval = nan;
end

% Add legend
strs = {'Lum: First Spike','Lum: Gaussian Window',...
    'Col:First Spike','Col: gaussian Window'};
legend(p(~isnan(p)),strs{~isnan(p)})

% Save figure variables
spikefig.scatfig = scatfig; % Different structure bc scatfig is a figure, not a uipanel with userdata
set(7548,'userdata',spikefig);

end


