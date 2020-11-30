function GLMSPopGUI_Gabors(varargin)
global GLMSPopData

figure(16); clf;
set(gcf,'units','normalized','pos',[.05 .1 .9 .8],'numbertitle','off','name','Gabors');

% Set up panels
gaborfig.conpanel = uipanel(gcf,'units','normalized','pos',[.01 .01 .19 .98]);
gaborfig.table = uitable('parent',gaborfig.conpanel,'units','normalized','pos',[.01 .5 .98 .49],...
    'BackgroundColor',[1 1 1]);
gaborfig.dispanel = uipanel(gcf,'units','normalized','pos',[.21 .01 .78 .98]);

% Control panel
conpanel.uicontrols.overview = uicontrol('parent',gaborfig.conpanel,'units','normalized',...
    'pos',[.01 .01 .98 .08],'style','pushbutton','string','Overview',...
    'fontsize',12,'callback',@GLMSGUI_Overview);
conpanel.selectedidx = [];
conpanel.uicontrols.reanalyze = uicontrol('parent',gaborfig.conpanel,'units','normalized',...
    'pos',[.01 .11 .98 .08],'style','pushbutton','string','Reanalyze All Data',...
    'fontsize',12,'callback',@reanalyzeall);
conpanel.uicontrols.surface = uicontrol('parent',gaborfig.conpanel,'units','normalized',...
    'pos',[.01 .21 .98 .08],'style','pushbutton','string','Surface Fits',...
    'fontsize',12,'callback',@surfacefit);
conpanel.uicontrols.showpop = uicontrol('parent',gaborfig.conpanel,'units','normalized',...
    'pos',[.01 .31 .98 .08],'style','pushbutton','string','Population Tuning Diffs',...
    'fontsize',12,'callback',@poptuningdiffs);
conpanel.uicontrols.corranal = uicontrol('parent',gaborfig.conpanel,'units','normalized',...
    'pos',[.01 .41 .98 .08],'style','pushbutton','string','Correlation',...
    'fontsize',12,'callback',@corranal);

% Set up disppanel
dispanel.axes.flash.hist = axes('parent',gaborfig.dispanel,'units','normalized',...
    'pos',[.05 .55 .25 .4],'box','on','xtick',[],'ytick',[]); axis equal square;
title('Flash PSTH')
dispanel.axes.gabor.hist = axes('parent',gaborfig.dispanel,'units','normalized',...
    'pos',[.05 .05 .25 .4],'box','on','xtick',[],'ytick',[]); axis equal square;
title('Gabor PSTH')
%dispanel.axes.flash.sta = axes('parent',gaborfig.dispanel,'units','normalized',...
%    'pos',[.375 .55 .25 .4],'box','on','xtick',[],'ytick',[]); axis equal square;
dispanel.axes.flash.sta(2) = axes('parent',gaborfig.dispanel,'units','normalized',...
    'pos',[.375 .55 .1 .15],'box','on','xtick',[],'ytick',[]); %axis equal square;
dispanel.axes.flash.sta(1) = axes('parent',gaborfig.dispanel,'units','normalized',...
    'pos',[.5 .75 .1 .15],'box','on','xtick',[],'ytick',[]); %axis equal square;
dispanel.axes.flash.sta(3) = axes('parent',gaborfig.dispanel,'units','normalized',...
    'pos',[.5 .55 .1 .15],'box','on','xtick',[],'ytick',[]); %axis equal square;
dispanel.axes.flash.sta(4) = axes('parent',gaborfig.dispanel,'units','normalized',...
    'pos',[.375 .75 .1 .15],'box','on','xtick',[],'ytick',[]); %axis equal square;
%title('WN STA')
dispanel.axes.gabor.sta = axes('parent',gaborfig.dispanel,'units','normalized',...
    'pos',[.375 .05 .25 .4],'box','on','xtick',[],'ytick',[]); axis equal square;
title('Gabor')
dispanel.axes.flash.bubbleplot = axes('parent',gaborfig.dispanel,'units','normalized',...
    'pos',[.7 .55 .25 .4],'box','on'); axis equal square;
title('Punctate Flashes')
dispanel.axes.gabor.bubbleplot = axes('parent',gaborfig.dispanel,'units','normalized',...
    'pos',[.7 .05 .25 .4],'box','on'); axis equal square;
title('Gabors')

% Save user variables
set(gaborfig.conpanel,'userdata',conpanel);
set(gaborfig.dispanel,'userdata',dispanel);
set(gcf,'userdata',gaborfig)

% Run Setup
UnpackPopulation();

if ~isempty(varargin)
    for n = 1:numel(varargin)
        eval(varargin{n})
    end
end


end

function UnpackPopulation(~,~)
global GLMSPopData datatypes

% Grab saved population data
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
%load([library 'GLMSPopData.mat'],'GLMSPopData')

% Load user data
gaborfig = get(gcf,'userdata');
conpanel = get(gaborfig.conpanel,'userdata');
dispanel = get(gaborfig.dispanel,'userdata');
datatypes = GLMSPopData(1,:);

% Some variables for culling gabors
maxdur = 100;
minnstim = 10;
minmod = 10; % in sp/s


% Call another analysis to get exclusion criteia (just easier than
% rebuilding)
GLMSPopGUI_Params
ParamsFig = get(552,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');
excludeL = ~poppanel.modL | ~poppanel.kappaL | poppanel.RSL | ~poppanel.uniqL;

% Find which datafiles have Gabor, and do initial weeding
glmps = GLMSPopData(2:end,strcmp(datatypes,'GLMP'));
gaborL = zeros(numel(glmps),1);
for n = 1:numel(glmps)
    if excludeL(n) == 0 % datafile has not already been excluded
        if any(glmps{n}.epoch == 3) % showed any Gabors
            if mean(glmps{n}.subunit{3}.stimDur) < maxdur % Weed out long duration gabors
                if numel(glmps{n}.subunit{3}.Lcc) > minnstim % Those with too few stimuli
                    if (max(glmps{n}.subunit{3}.meanfr)-mean(glmps{n}.subunit{3}.blfr)) > minmod % Must drive at least 10 sp/s
                        gaborL(n) = 1;
                    end
                end
            end
        end
    end
end
conpanel.gaboridx = find(gaborL);

% Symmetrize gabor and flash responses for correlation measurement
% Preallocate space for symmetrized data
conpanel.corranal.corrcoef = nan(size(conpanel.gaboridx));
conpanel.corranal.stim = cell(size(conpanel.gaboridx));
conpanel.corranal.flash.allresps = cell(size(conpanel.gaboridx));
conpanel.corranal.flash.maxlumresp = nan(size(conpanel.gaboridx));
conpanel.corranal.flash.maxcolresp = nan(size(conpanel.gaboridx));
conpanel.corranal.flash.mod = nan(size(conpanel.gaboridx));
conpanel.corranal.gabor.allresps = cell(size(conpanel.gaboridx));
conpanel.corranal.gabor.maxlumresp = nan(size(conpanel.gaboridx));
conpanel.corranal.gabor.maxcolresp = nan(size(conpanel.gaboridx));
conpanel.corranal.gabor.mod = nan(size(conpanel.gaboridx));

% Which datafiles have symmetrized resps that are less than minmod?
for n = 1:numel(conpanel.gaboridx)
    
    % Pull out data from current file
    idx = conpanel.gaboridx(n);
    GLMP = GLMSPopData{idx+1,strcmp(datatypes,'GLMP')};
    sub = GLMSPopData{idx+1,strcmp(datatypes,'Subunit')};
    glmp = GLMP.subunit{sub};
    gabor = GLMP.subunit{3};
    
    % Symmetrize responses across origin
    glmpsymresps = nan(size(glmp.uniqueLcc));
    for i = 1:numel(glmp.uniqueLcc)
        stim = [glmp.uniqueLcc(i) glmp.uniqueMcc(i)];
        symidx = all(softEq([glmp.Lcc glmp.Mcc],stim),2) | ...
            all(softEq([glmp.Lcc glmp.Mcc],-stim),2);
        glmpsymresps(i) = mean(glmp.fr(symidx));
    end
    
    % Symmetrize responses across origin
    gaborsymresps = nan(size(gabor.uniqueLcc));
    for i = 1:numel(gabor.uniqueLcc)
        stim = [gabor.uniqueLcc(i) gabor.uniqueMcc(i)];
        symidx = all(softEq([gabor.Lcc gabor.Mcc],stim),2) | ...
            all(softEq([gabor.Lcc gabor.Mcc],-stim),2);
        gaborsymresps(i) = mean(gabor.fr(symidx));
        gabor.fr(symidx)
    end
    
    % Find only commonly tested stimuli
    [c,ai,bi] = intersect([glmp.uniquetheta glmp.uniquerho],[gabor.uniquetheta gabor.uniquerho],'rows');
    glmpresps = glmpsymresps(ai);
    gaborresps = gaborsymresps(bi);
        
    % Catalgoue resp to highest lum and col contrast
    temp = find(c(:,1)==pi/4 | c(:,1)==-3*pi/4);
    [~,tempidx] = max(c(temp,2));
    lumidx = temp(tempidx);
    
    temp = find(c(:,1)==-pi/4 | c(:,1)==3*pi/4);
    [~,tempidx] = max(c(temp,2));
    colidx = temp(tempidx);
    
    
    % GLMP bubble plot
%     axlim = max(glmp.Lcc);
%     maxresp = max(cat(1,glmpresps,gaborresps));
%     figure(1); clf;
%     for i = 1:size(c,1)
%         mn = glmpresps(i)/maxresp*30+3; hold on;
%         h = polar(c(i,1),c(i,2),'ko');
%         set(h,'MarkerFaceColor','none','MarkerSize',mn)
%     end
%     box on; grid off;
%     axis equal square tight
%     set(gca,'CameraPosition',[0 0 275.3526],'xlim',[-axlim axlim],'ylim',[-axlim axlim])
% 
%     % Gabor bubble plot
%     figure(2); clf;
%     for i = 1:size(c,1)
%         mn = gaborresps(i)/maxresp*30+3; hold on;
%         h = polar(c(i,1),c(i,2),'ko');
%         set(h,'MarkerFaceColor','none','MarkerSize',mn)
%     end
%     box on; grid off;
%     axis equal square tight
%     set(gca,'xlim',[-axlim axlim],'ylim',[-axlim axlim])
    
    if max(gaborresps) - mean(gabor.blfr) < minmod
        gaborL(conpanel.gaboridx(n)) = 0;
    else
        conpanel.corranal.corrcoef(n) = corr(glmpresps,gaborresps);
        conpanel.corranal.stim{n} = c;
        conpanel.corranal.flash.allresps{n} = glmpresps;
        conpanel.corranal.flash.maxlumresp(n) = glmpresps(lumidx);
        conpanel.corranal.flash.maxcolresp(n) = glmpresps(colidx);
        conpanel.corranal.flash.mod(n) = max(glmpresps) - mean(glmp.blfr);
        conpanel.corranal.gabor.allresps{n} = gaborresps;
        conpanel.corranal.gabor.maxlumresp(n) = gaborresps(lumidx);
        conpanel.corranal.gabor.maxcolresp(n) = gaborresps(colidx);
        conpanel.corranal.gabor.mod(n) = max(gaborresps) - mean(gabor.blfr);
    end
    
end
conpanel.gaboridx = find(gaborL);
L = ~isnan(conpanel.corranal.corrcoef);
conpanel.corranal.corrcoef = conpanel.corranal.corrcoef(L);
conpanel.corranal.stim = conpanel.corranal.stim(L);
conpanel.corranal.flash.allresps = conpanel.corranal.flash.allresps(L);
conpanel.corranal.flash.maxlumresp = conpanel.corranal.flash.maxlumresp(L);
conpanel.corranal.flash.maxcolresp = conpanel.corranal.flash.maxcolresp(L);
conpanel.corranal.flash.mod = conpanel.corranal.flash.mod(L);
conpanel.corranal.gabor.allresps = conpanel.corranal.gabor.allresps(L);
conpanel.corranal.gabor.maxlumresp = conpanel.corranal.gabor.maxlumresp(L);
conpanel.corranal.gabor.maxcolresp = conpanel.corranal.gabor.maxcolresp(L);
conpanel.corranal.gabor.mod = conpanel.corranal.gabor.mod(L);

% Build table
filenames = GLMSPopData(conpanel.gaboridx+1,strcmp(datatypes,'Datafile'));
data(:,1) = filenames;
nds = cell(numel(conpanel.gaboridx),1);
[nds{poppanel.oneD.L}] = deal('1D');
[nds{~poppanel.oneD.L}] = deal('2D');
data(:,2) = nds(conpanel.gaboridx);
colname = {'Datafile' 'n-D'};
colformat = cell(1,size(data,2));
colformat(:) = {'char'};
set(gaborfig.table,'data',data,'columnname',colname,'rowname',conpanel.gaboridx,...
    'columnformat',colformat,'CellSelectionCallback',@cellselect);

% Unpack for easy reference
dispanel.gaborfits = cell(numel(conpanel.gaboridx),1);
dispanel.flashfits = cell(numel(conpanel.gaboridx),1);
dispanel.gaborparams = nan(numel(conpanel.gaboridx),8);
dispanel.flashparams = nan(numel(conpanel.gaboridx),8);
dispanel.oneDL = nan(numel(conpanel.gaboridx),1);

for n = 1:numel(conpanel.gaboridx)
    dispanel.flashfits{n} = GLMSPopData{conpanel.gaboridx(n)+1,strcmp(datatypes,'Surface Parameters')};
    dispanel.gaborfits{n} = GLMSPopData{conpanel.gaboridx(n)+1,strcmp(datatypes,'Gabor')};
    if dispanel.flashfits{n}.diffnormLL < .07
        dispanel.oneDL(n) = true;
        dispanel.flashparams(n,:) = dispanel.flashfits{n}.oneD.parvals;
        dispanel.gaborparams(n,:) = dispanel.gaborfits{n}.oneD.parvals;
    else
        dispanel.oneDL(n) = false;
        dispanel.flashparams(n,:) = dispanel.flashfits{n}.twoD.parvals;
        dispanel.gaborparams(n,:) = dispanel.gaborfits{n}.twoD.parvals;
    end 
end


set(gaborfig.conpanel,'userdata',conpanel);
set(gaborfig.dispanel,'userdata',dispanel);
set(16,'userdata',gaborfig)

end


%%% For cell selection

function cellselect(a,b)
global GLMSPopData DN GLMP 

tic

% Load figure and pop variables
gaborfig = get(16,'userdata');
conpanel = get(gaborfig.conpanel,'userdata');

% Set aside the index (remember to +1 for referencing GLMSPopData)
if ~isempty(b)
    if ~isempty(b.Indices)
        if conpanel.selectedidx == b.Indices(1)
            return
        else
            conpanel.selectedidx = b.Indices(1);
        end
    end
elseif isnumeric(a)
    conpanel.selectedidx = find(conpanel.gaboridx == a);
else
    return
end

% Pull out data from pop struct
datanames = GLMSPopData(1,:);
idx = conpanel.gaboridx(conpanel.selectedidx);
conpanel.popidx = idx+1;
DN = GLMSPopData{conpanel.popidx,strcmp(datanames,'DN')};
GLMP = GLMSPopData{conpanel.popidx,strcmp(datanames,'GLMP')};

% Save figure variables
set(gaborfig.conpanel,'userdata',conpanel)
set(gcf,'userdata',gaborfig)

% Populate different displays
showHists()
showDNpanel()
showGabor()
BubblePlot()

toc

end

function showHists()
global GLMP GLMSPopData

% Load user data
gaborfig = get(gcf,'userdata');
conpanel = get(gaborfig.conpanel,'userdata');
dispanel = get(gaborfig.dispanel,'userdata');
datatypes = GLMSPopData(1,:);

% Flash PSTH
sub = GLMSPopData{conpanel.popidx,strcmp(datatypes,'Subunit')};
binsize = .001; % in seconds
bins = -.2:binsize:.5;
gaussSize = .01;
gaussSize = ceil(gaussSize/binsize);
gaussfilter = gausswin(gaussSize,3);
gaussfilter = gaussfilter./sum(gaussfilter); %normalizing
PSTH = histc(cat(1,GLMP.subunit{sub}.normspiketimes{:}),bins);
PSTH = (PSTH./numel(GLMP.subunit{sub}.Lcc)) ./ binsize;
smoothPSTH = conv(PSTH,gaussfilter,'same');
nrows = numel(GLMP.subunit{sub}.Lcc);
rowcoords = linspace(max(smoothPSTH),0,nrows*2+1);
rowcents = rowcoords(2:2:end);
rowspace = rowcoords(2)-rowcoords(1);

% plot GLMP hist
axes(dispanel.axes.flash.hist); 
yyaxis right; cla; 
yyaxis left; cla; hold on; box on;
plot(bins,smoothPSTH,'linestyle','-','color',[.5 .5 .5])
plot([bins(1) bins(end)],[GLMP.blfrthresh GLMP.blfrthresh],'r--')
plot([GLMP.countingwin(1) GLMP.countingwin(1)],[0 max(smoothPSTH)],'r--')
plot([GLMP.countingwin(2) GLMP.countingwin(2)],[0 max(smoothPSTH)],'r--')
x = [];
y = [];
thetarhos = cat(2,GLMP.subunit{sub}.theta,GLMP.subunit{sub}.rho);
[~,order] = sortrows(thetarhos);
order = flipud(order);
for n = 1:numel(order)
    tpts = GLMP.subunit{sub}.normspiketimes{order(n)};
    if ~isempty(tpts)
        x = cat(1,x,repmat(tpts,1,2));
        y = cat(1,y,repmat([rowcents(n)-rowspace rowcents(n)+rowspace],numel(tpts),1));
    end
end

yyaxis left
plot(fliplr(x'),fliplr(y'),'k-')
xlim([bins(1) bins(end)])
ylim([0 max(smoothPSTH)])
set(gca,'xlim',[min(bins) max(bins)],'ytick',0:10:1000)
title('PSTH')
xlabel('Time from Stim Onset (ms)')
ylabel('Mean Firing Rate (sp/s)')

yyaxis right;
ylabel('Color Direction (rad)')
coldirs = linspace(-180,180,5);
set(gca,'ylim',[min(coldirs) max(coldirs)],'YTick',coldirs,...
    'XTick',bins(1):.1:bins(end))
h = ylabel('Angle in LM plane');
set(h,'rotation',-90)

% gabor PSTH
sub = 3;
binsize = .001; % in seconds
bins = -.2:binsize:max(GLMP.subunit{3}.stimDur)+.2;
gaussSize = .01;
gaussSize = ceil(gaussSize/binsize);
gaussfilter = gausswin(gaussSize,3);
gaussfilter = gaussfilter./sum(gaussfilter); %normalizing
PSTH = histc(cat(1,GLMP.subunit{sub}.normspiketimes{:}),bins);
PSTH = (PSTH./numel(GLMP.subunit{sub}.Lcc)) ./ binsize;
smoothPSTH = conv(PSTH,gaussfilter,'same');
nrows = numel(GLMP.subunit{sub}.Lcc);
rowcoords = linspace(max(smoothPSTH),0,nrows*2+1);
rowcents = rowcoords(2:2:end);
rowspace = rowcoords(2)-rowcoords(1);

% plot gabor hist
axes(dispanel.axes.gabor.hist);
yyaxis right; cla; 
yyaxis left; cla; hold on; box on;
plot(bins,smoothPSTH,'linestyle','-','color',[.5 .5 .5])
meandur = mean(GLMP.subunit{sub}.stimDur);
plot([bins(1) bins(end)],[GLMP.blfrthresh GLMP.blfrthresh],'r--')
plot([0 0],[0 max(smoothPSTH)],'r--')
plot([meandur meandur],[0 max(smoothPSTH)],'r--')
x = [];
y = [];
thetarhos = cat(2,GLMP.subunit{sub}.theta,GLMP.subunit{sub}.rho);
[~,order] = sortrows(thetarhos);
%uniq = flipud(uniq);
order = flipud(order);
for n = 1:numel(order)
    tpts = GLMP.subunit{sub}.normspiketimes{order(n)};
    if ~isempty(tpts)
        x = cat(1,x,repmat(tpts,1,2));
        y = cat(1,y,repmat([rowcents(n)-rowspace rowcents(n)+rowspace],numel(tpts),1));
    end
end

yyaxis left
plot(fliplr(x'),fliplr(y'),'k-')
xlim([bins(1) bins(end)])
ylim([0 max(smoothPSTH)])
set(gca,'xlim',[min(bins) max(bins)],'ytick',0:10:1000)
title('PSTH')
xlabel('Time from Stim Onset (ms)')
ylabel('Mean Firing Rate (sp/s)')

yyaxis right
ylabel('Color Direction (rad)')
coldirs = linspace(-180,180,5);
set(gca,'ylim',[min(coldirs) max(coldirs)],'YTick',coldirs)
ylabel('Angle in LM plane')
set(gca,'XTick',bins(1):.1:bins(end))


% Save figure variables
set(gaborfig.conpanel,'userdata',conpanel)
set(gaborfig.dispanel,'userdata',dispanel)
set(gcf,'userdata',gaborfig)

end

function showGabor()
global DN GLMP

% Load user data
gaborfig = get(gcf,'userdata');
%conpanel = get(gaborfig.conpanel,'userdata');
dispanel = get(gaborfig.dispanel,'userdata');

% Define spatial and temporal profile of Gabor
dvaPerStixel = DN.DVAPerStix(1);
pixperdeg = 800; % this is arbitrary to set the meshgrid
pixperstix = round(pixperdeg*dvaPerStixel);
substixX = GLMP.subunit{1}.gridX{1};
substixY = GLMP.subunit{1}.gridY{1};
subpixX = [];
subpixY = [];
for n = 1:numel(substixX)
    tempX = substixX(n)*pixperstix:substixX(n)*pixperstix+pixperstix;
    tempY = substixY(n)*pixperstix:substixY(n)*pixperstix+pixperstix;
    [temppixX, temppixY] = meshgrid(tempX,tempY);
    subpixY = cat(1,subpixY,temppixY(:));
    subpixX = cat(1,subpixX,temppixX(:));
end

coeff = pca([subpixX subpixY]);
PC_1 = coeff(:,1);
PC_end = [PC_1(2); -PC_1(1)];
theta = mod(atan2(PC_1(2),-PC_1(1)),pi); % Preferred direction. 0 means horizonal
proj = [subpixX subpixY] * PC_end;
halfperiod = range(proj); % in pixels
sf = 1/(halfperiod * 2 / pixperdeg); % in cycles/deg
sigma = .4; % what is this?
%nframesplateau = 26;
%nframesramp = 12;
phi = 0;
gamma = 1;
%tf = 3;
nsigmas = 3;
lambda = 1/sf;
%zeropadding = zeros(1,10);
%ramp = linspace(0,1,nframesramp);
%plateau = ones(1,nframesplateau);
%temporalprofile = [zeropadding ramp plateau fliplr(ramp) zeropadding];
%nframes = numel(temporalprofile);
stimsizeindeg = sigma*nsigmas; % This is half stim size (we we go out +/- nsigmas)
stimsizeinpix = round(stimsizeindeg*pixperdeg);  % full stim size in doublewide pixels
[x,y] = meshgrid(linspace(-1,1,stimsizeinpix), linspace(-1,1,stimsizeinpix));% x and y are in dva
xprime = x*cos(-theta) + y*sin(-theta);
yprime = -x*sin(-theta) + y*cos(-theta);
% deltaphase = tf * 2 * pi / DN.framerate;
% gabor = ones(stimsizeinpix,stimsizeinpix,3,nframes);
% for a = 1:length(temporalprofile)
%     phase = phi+(a-1)*deltaphase;
%     gabor(:,:,:,a) = repmat(temporalprofile(a)*exp(-(xprime.^2 + gamma.^2 .* yprime.^2)./ (2*sigma^2)).*cos(2*pi*yprime./lambda+phase),[1 1 3]);
% end
gabor = repmat(exp(-(xprime.^2 + gamma.^2 .* yprime.^2)./ (2*sigma^2)).*cos(2*pi*yprime./lambda+phi),[1 1 3]);
gabor = gabor - min(gabor(:));
gabor = gabor./max(gabor(:));
dispanel.gabormov = gabor;

% Show the center frame of movie
axes(dispanel.axes.gabor.sta); cla; hold on; axis tight equal
%dispanel.gaborimg = image(x(:),y(:),gabor(:,:,:,floor(numel(temporalprofile)/2)));
dispanel.gaborimg = image(x(:),y(:),gabor);
plot(subpixX/pixperdeg,subpixY/pixperdeg,'r','marker','square')
%set(dispanel.gaborimg,'ButtonDownFcn',@playmov)
title('Gabor')
set(gca,'xtick',linspace(-1,1,5),'xticklabel',linspace(-stimsizeindeg,stimsizeindeg,5),...
    'ytick',linspace(-1,1,5),'yticklabel',linspace(-stimsizeindeg,stimsizeindeg,5),...
    'Ydir','normal');

%     function playmov(~,~)
%         
%         % Show gabor movie
%         axes(dispanel.axes.gabor.sta); cla; hold on;
%         set(gca,'nextplot','replacechildren')
%         for i = 1:length(temporalprofile)
%             image(gabor(:,:,:,i)); axis tight
%             %plot(subpixX/pixperdeg,subpixY/pixperdeg,'ro')
%             pause(1/DN.framerate)
%         end
%         
%         % Reset to still image
%         axes(dispanel.axes.gabor.sta); cla;
%         set(dispanel.axes.gabor.sta,'XTick',[],'Ytick',[])
%         dispanel.gaborimg = image(gabor(:,:,:,floor(numel(temporalprofile)/2)));
%         set(dispanel.gaborimg,'ButtonDownFcn',@playmov)
%         title('Gabor')
%         
%     end

% Save image variables
set(gaborfig.dispanel,'userdata',dispanel)
%set(gaborfig.conpanel,'userdata',conpanel)
set(gcf,'userdata',gaborfig)


end

function showDNpanel()
global DN GLMP

% Load user data
gaborfig = get(gcf,'userdata');
%conpanel = get(gaborfig.conpanel,'userdata');
dispanel = get(gaborfig.dispanel,'userdata');

nspikes = DN.stats.pLpM.nspikes; % number of spikes
nstix = numel(GLMP.subunit{1}.gridX{1});  % number of stixels

% Pull out the correct STS info
pLpM = DN.stats.pLpM.STS(1:DN.NStixGrid(1)^2,:);
mLmM = abs(DN.stats.mLmM.STS(1:DN.NStixGrid(1)^2,:));
pLmM = abs(DN.stats.pLmM.STS(1:DN.NStixGrid(1)^2,:));
mLpM = abs(DN.stats.mLpM.STS(1:DN.NStixGrid(1)^2,:));

% Narrow down to the subunit stixels
subx = cell(2,1);
suby = subx;
subidx = suby;
for n = 1:2
    if ~isempty(GLMP.subunit{n})
        subx{n} = GLMP.subunit{n}.gridX{1} + ceil(DN.NStixGrid(1)/2);
        suby{n} = -GLMP.subunit{n}.gridY{1} + ceil(DN.NStixGrid(1)/2);
        subidx{n} = sub2ind([DN.NStixGrid(1) DN.NStixGrid(1)],subx{n},suby{n});
    end
end

% Pull out just subunit stixels on each color channel and time bin
pLpMstix = sum(pLpM(subidx{1},:),1);
mLmMstix = sum(mLmM(subidx{1},:),1);
pLmMstix = sum(pLmM(subidx{1},:),1);
%mLpMstix = sum(1,mLpM(subidx,:)); % leaving this one out bc the distribution is rank 3
stixSTS = cat(1,pLpMstix,mLmMstix,pLmMstix);
 
% This code just makes the mean vector and covariance matrix
mu = nspikes * nstix * .25;
sigma = repmat(-nspikes * nstix * .25^2,3,3);
sigma(1,1) = nspikes * nstix * .25 * .75; % This is n*p*(1-p)
sigma(2,2) = nspikes * nstix * .25 * .75;
sigma(3,3) = nspikes * nstix * .25 * .75;

% Get significance for selected stixels at each time point
ps = nan(1,size(stixSTS,2));
X = ps;
for n = 1:size(stixSTS,2)
    X(n) = (stixSTS(:,n) - mu)' * (sigma \ (stixSTS(:,n) - mu));
    dispanel.signif(n) = 1-chi2pdf(X(n),3);
end
[~,tpt] = max(dispanel.signif); % find frame with greatest deviation from mean (or least prob of being rando draw from dist)

% Pull out selected DN frame and convert STS to STA
STApLpM = pLpM(:,tpt)./nspikes;
STAmLmM = mLmM(:,tpt)./nspikes;
STApLmM = pLmM(:,tpt)./nspikes;
STAmLpM = mLpM(:,tpt)./nspikes;

% Convert from color channels to LMS
lumcc = sqrt((DN.lumCC(1).^2)/2);
colcc = sqrt((DN.colCC(1).^2)/2);
Lconemod = STApLpM * lumcc;
Lconemod = Lconemod + (STAmLmM * -lumcc); 
Lconemod = Lconemod + (STApLmM * colcc);
Lconemod = Lconemod + (STAmLpM * -colcc);
Mconemod = STApLpM * lumcc;
Mconemod = Mconemod + (STAmLmM * -lumcc);
Mconemod = Mconemod + (STApLmM * -colcc);
Mconemod = Mconemod + (STAmLpM * colcc);
Sconemod = zeros(size(Lconemod));

% convert from LMS to RGB
conemods = cat(2,Lconemod,Mconemod,Sconemod);
gunmods = conemods * inv(GLMP.M); % make sure to put M mtx in DN struct

% normalize STA
mincc = min(gunmods(:));
temp = gunmods - mincc;
normRGB = temp ./ max(temp(:));

% reshape STAs
Rmod = reshape(normRGB(:,1),[sqrt(size(gunmods,1)) sqrt(size(gunmods,1))]);
Gmod = reshape(normRGB(:,2),[sqrt(size(gunmods,1)) sqrt(size(gunmods,1))]);
Bmod = reshape(normRGB(:,3),[sqrt(size(gunmods,1)) sqrt(size(gunmods,1))]);
dispanel.gunSTAmat = cat(3,Rmod,Gmod,Bmod);

% Display STA
% axes(dispanel.axes.flash.sta); hold on; cla; axis tight
% dispanel.STA = image(dispanel.gunSTAmat);
% set(dispanel.STA,'buttondownfcn',@showsubs)
% set(dispanel.axes.flash.sta,'xtick',[],'ytick',[],'Ydir','reverse')
% title('WN STA')

stas = cat(2,STApLpM,STAmLmM,STApLmM,STAmLpM);
for n = 1:4
    normSTA = (stas(:,n) - min(stas(:)))/(max(stas(:))-min(stas(:)));
    STA = repmat(reshape(normSTA,[sqrt(size(stas,1)) sqrt(size(stas,1))]),[1 1 3]);
    cla(dispanel.axes.flash.sta(n));
    dispanel.STA(n) = image(STA,'parent',dispanel.axes.flash.sta(n));
    %axis square equal
    set(dispanel.STA(n),'buttondownfcn',@showsubs)
    set(dispanel.axes.flash.sta(n),'xtick',[],'ytick',[])%,'Ydir','reverse')
end

subon = 0;
showsubs()

    function showsubs(~,~)
        
        % clear axes
%         axes(dispanel.axes.flash.sta); cla; hold on;
%         set(dispanel.axes.flash.sta,'xtick',[],'ytick',[],'ydir','reverse')
%         dispanel.STA = image(dispanel.gunSTAmat);
%         set(dispanel.STA,'buttondownfcn',@showsubs)
%         axis tight
%         title('WN STA')
        
        % if subs were on, leave unplotted. If they were off, plot them.
        if subon == 0
            cols = [1 0 0; 0 1 0];
            for q = 1:4
                for i = 1:2
                    if ~isempty(GLMP.subunit{i})
                        hold(dispanel.axes.flash.sta(q),'on');
                        plot(subx{i},suby{i},'o','color',cols(i,:),'parent',dispanel.axes.flash.sta(q))
                    end
                end
            end
            subon = 1;
        else
%             for q = 1:4
%                 dispanel.STA(q) = image(STA,'parent',dispanel.axes.flash.sta(q)); hold on;
%             end
            subon = 0;
        end
    end

% Save image variables
set(gaborfig.dispanel,'userdata',dispanel)
set(gcf,'userdata',gaborfig)

end

function BubblePlot()
global GLMSPopData GLMP datatypes

% Load user data
gaborfig = get(gcf,'userdata');
conpanel = get(gaborfig.conpanel,'userdata');
dispanel = get(gaborfig.dispanel,'userdata');
datatypes = GLMSPopData(1,:);

% Plot Gabor data
% surfparams = GLMSPopData{conpanel.popidx,strcmp(datatypes,'Surface Parameters')};
% sub = GLMSPopData{conpanel.popidx,strcmp(datatypes,'Subunit')};
% whichnd = GLMSPopData{conpanel.popidx,strcmp(datatypes,'nD')};
% if strcmp(whichnd,'1D')
%     nd = 'oneD';
%     PD = [-max(GLMP.subunit{sub}.rho) 0; max(GLMP.subunit{sub}.rho) 0];
% else
%     nd = 'twoD';
%     PD = [-max(GLMP.subunit{sub}.rho) 0; max(GLMP.subunit{sub}.rho) 0;...
%         0 max(GLMP.subunit{sub}.rho); 0 -max(GLMP.subunit{sub}.rho)];
% end
% 
% axes(dispanel.axes.flash.bubbleplot); cla;
% params = surfparams.(nd).parvals;
% surftype = surfparams.(nd).surftype;
% meannsp = GLMP.subunit{sub}.meannspikes;
% ticks = -.6:.2:.6;
% colormap('cool')

% % Define thetas and rhos for the whole grid
% thetas = linspace(-pi,pi,361)';
% majorax = max(GLMP.subunit{sub}.uniquerho(GLMP.subunit{sub}.uniquetheta==pi/4))*1.05;
% minorax = max(GLMP.subunit{sub}.uniquerho(GLMP.subunit{sub}.uniquetheta==-pi/4))*1.5;
% nom = majorax * minorax;
% denom = (majorax*sin(thetas-pi/4)).^2 + (minorax*cos(thetas-pi/4)).^2;
% rhos = nom ./ sqrt(denom);
% scalars = 0:.05:1;
% rhosgrid = rhos * scalars;
% thetasgrid = repmat(thetas,1,numel(scalars));
% [x,y] = pol2cart(thetasgrid,rhosgrid);
% surface = ComputeNakaRushtonJPW(params,[x(:) y(:)],surftype);
% surface = reshape(surface,size(rhosgrid));
% axlim = max(x(:));

% % (1) Surface
% h = surf(x,y,surface); hold on;
% set(h,'edgecolor','none');
% alpha(.3);
% set(gca,'XTick',[],'YTick',[],'xlim',[-axlim axlim],'ylim',[-axlim axlim]); 
% 
% % (2) contours and pd
% rotmat = [cos(params(end-1)) sin(params(end-1)); -sin(params(end-1)) cos(params(end-1))];
% rotPD = PD * rotmat;
% contour3(x,y,surface,5,'linewidth',2); hold on; % masked grid
% if strcmp(nd,'oneD')
%     plot(rotPD(:,1),rotPD(:,2),'k-')
% else
%     plot(rotPD(1:2,1),rotPD(1:2,2),'k-')
%     plot(rotPD(3:4,1),rotPD(3:4,2),'k-')
% end
% set(gca,'tickdir','out','XTick',ticks,'YTick',ticks,'xlim',[-axlim axlim],'ylim',[-axlim axlim]);

% % (3) Bubbles
% polar(0,axlim,'k*'); hold on; % this is stupid, just to get axis limits right
% for i = 1:numel(GLMP.subunit{sub}.uniqueLcc)
%     mn = meannsp(i)/max(meannsp)*30+3; hold on;
%     h = polar(GLMP.subunit{sub}.uniquetheta(i),GLMP.subunit{sub}.uniquerho(i),'ko');
%     set(h,'MarkerFaceColor','none','MarkerSize',mn)
% end
% set(gca,'XTick',[],'YTick',[],'xlim',[-axlim axlim],'ylim',[-axlim axlim]); 
% 
% 
% % axis stuff
% box on; grid off;
% axis equal square tight
% set(gca,'CameraPosition',[0 0 275.3526],'xlim',[-axlim axlim],'ylim',[-axlim axlim])


% % Plot Gabor data
% gaborparams = GLMSPopData{conpanel.popidx,strcmp(datatypes,'Gabor')};
% sub = 3;
% gaborparams.diffnormLL
% if gaborparams.diffnormLL < .07
%     nd = 'oneD'
% else
%     nd = 'twoD'
% end
% 
% axes(dispanel.axes.gabor.bubbleplot); cla;
% params = gaborparams.fitcomps.(nd).parvals(2,:);
% surftype = gaborparams.(nd).surftype;
% meannsp = GLMP.subunit{sub}.meannspikes;
% ticks = -.6:.2:.6;
% colormap('cool')
% 
% 
% % Define thetas and rhos for the whole grid
% thetas = linspace(-pi,pi,361)';
% majorax = max(GLMP.subunit{sub}.uniquerho(GLMP.subunit{sub}.uniquetheta==pi/4))*1.05;
% minorax = max(GLMP.subunit{sub}.uniquerho(GLMP.subunit{sub}.uniquetheta==-pi/4))*1.5;
% nom = majorax * minorax;
% denom = (majorax*sin(thetas-pi/4)).^2 + (minorax*cos(thetas-pi/4)).^2;
% rhos = nom ./ sqrt(denom);
% scalars = 0:.05:1;
% rhosgrid = rhos * scalars;
% thetasgrid = repmat(thetas,1,numel(scalars));
% [x,y] = pol2cart(thetasgrid,rhosgrid);
% surface = ComputeNakaRushtonJPW(params,[x(:) y(:)],surftype);
% surface = reshape(surface,size(rhosgrid));
% axlim = max(x(:));
% 
% % (1) Surface
% h = surf(x,y,surface); hold on;
% set(h,'edgecolor','none');
% alpha(.3);
% set(gca,'XTick',[],'YTick',[],'xlim',[-axlim axlim],'ylim',[-axlim axlim]); 
% 
% % (2) contours and pd
% rotmat = [cos(params(end-1)) sin(params(end-1)); -sin(params(end-1)) cos(params(end-1))];
% rotPD = PD * rotmat;
% contour3(x,y,surface,5,'linewidth',2); hold on; % masked grid
% if strcmp(nd,'oneD')
%     plot(rotPD(:,1),rotPD(:,2),'k-')
% else
%     plot(rotPD(1:2,1),rotPD(1:2,2),'k-')
%     plot(rotPD(3:4,1),rotPD(3:4,2),'k-')
% end
% set(gca,'tickdir','out','XTick',ticks,'YTick',ticks,'xlim',[-axlim axlim],'ylim',[-axlim axlim]);
% 
% 
% % (3) Bubbles
% polar(0,axlim,'k*'); hold on; % this is stupid, just to get axis limits right
% for i = 1:numel(GLMP.subunit{sub}.uniqueLcc)
%     mn = meannsp(i)/max(meannsp)*30+3; hold on;
%     h = polar(GLMP.subunit{sub}.uniquetheta(i),GLMP.subunit{sub}.uniquerho(i),'ko');
%     set(h,'MarkerFaceColor','none','MarkerSize',mn)
% end
% set(gca,'XTick',[],'YTick',[],'xlim',[-axlim axlim],'ylim',[-axlim axlim]); 


% % axis stuff
% box on; grid off;
% axis equal square tight
% set(gca,'CameraPosition',[0 0 275.3526],'xlim',[-axlim axlim],'ylim',[-axlim axlim])

    stim = conpanel.corranal.stim{conpanel.selectedidx};
    flashresps = conpanel.corranal.flash.allresps{conpanel.selectedidx};
    gaborresps = conpanel.corranal.gabor.allresps{conpanel.selectedidx};
    maxresp = max(cat(1,flashresps,gaborresps));
    axlim = max(GLMP.subunit{1}.Lcc);
    
    % GLMP bubble plot
    axes(dispanel.axes.flash.bubbleplot); cla;
    for i = 1:size(stim,1)
        mn = flashresps(i)/maxresp*30+3; hold on;
        h = polar(stim(i,1),stim(i,2),'ko');
        set(h,'MarkerFaceColor','none','MarkerSize',mn)
    end
    box on; grid off;
    axis equal square tight
    set(gca,'CameraPosition',[0 0 275.3526],'xlim',[-axlim axlim],'ylim',[-axlim axlim])
    legend(['max resp = ' num2str(max(flashresps))],'location','northwest')
    
    
    % Gabor bubble plot
    axes(dispanel.axes.gabor.bubbleplot); cla;
    for i = 1:size(stim,1)
        mn = gaborresps(i)/maxresp*30+3; hold on;
        h = polar(stim(i,1),stim(i,2),'ko');
        set(h,'MarkerFaceColor','none','MarkerSize',mn)
    end
    box on; grid off;
    axis equal square tight
    set(gca,'xlim',[-axlim axlim],'ylim',[-axlim axlim])
    legend(['max resp = ' num2str(max(gaborresps))],'location','northwest')


    disp(['Corr = ' num2str(conpanel.corranal.corrcoef(conpanel.selectedidx))])

end


%%% Analysis

function surfacefit(~,~,~)
global GLMSPopData

% Load user data
gaborfig = get(gcf,'userdata');
conpanel = get(gaborfig.conpanel,'userdata');
dispanel = get(gaborfig.dispanel,'userdata');
datatypes = GLMSPopData(1,:);
sub = 3;

GLMSGUI_Surface([],[],'Analyze',sub)
    
% Get parameter values from figure variables
SurfFig = get(60,'userdata');
fitspanel = get(SurfFig.fitspanel,'userdata');

% Organize parameters into population structure and save
oneD.surftype = fitspanel.surftype;
twoD.surftype = fitspanel.surftype;
oneD.parvals = fitspanel.params.oneD;
twoD.parvals = fitspanel.params.twoD;
oneD.LL = fitspanel.LL.oneD;
twoD.LL = fitspanel.LL.twoD;
oneD.normLL = fitspanel.normLL.oneDnormLL;
twoD.normLL = fitspanel.normLL.twoDnormLL;
gaborparams = struct('oneD',oneD,'twoD',twoD);
gaborparams.diffnormLL = twoD.normLL - oneD.normLL;
gaborparams.fitcomps = fitspanel.fitcomps;

GLMSPopData{conpanel.popidx,strcmp(datatypes,'Gabor')} = gaborparams;

if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
save([library 'GLMSPopData.mat'],'GLMSPopData')


% Save image variables
set(gaborfig.dispanel,'userdata',dispanel)
set(gaborfig.conpanel,'userdata',conpanel)
set(gcf,'userdata',gaborfig)

end

function reanalyzeall(~,~,~)
global GLMSPopData GLMP datatypes


% Load user data
gaborfig = get(gcf,'userdata');
conpanel = get(gaborfig.conpanel,'userdata');

% Load, Organize, and Save Rawdata
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/nex files/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\nex files\';
end
datanames = GLMSPopData(1,:);  

if ~any(strcmp(datanames,'Gabor'))
    GLMSPopData{1,end+1} = 'Gabor';
end

for n = 1:numel(conpanel.gaboridx)
    disp(['Analyzing ' num2str(n) ' of ' num2str(numel(conpanel.gaboridx))])
    
    conpanel.popidx = conpanel.gaboridx(n)+1; %+1 for GLMSPopData indexing (for reorganizing all datafiles)
    %GLMP = GLMSPopData{conpanel.popidx,strcmp(datanames,'GLMP')};
    %set(gaborfig.conpanel,'userdata',conpanel);
    %surfacefit();
    
    datafile = GLMSPopData{conpanel.popidx,strcmp(datatypes,'Datafile')};
    rawdata = nex2stro([library char(datafile) '.nex']);
    [GLMP,~] = OrganizeRawGLMSData(rawdata);
    %GLMSPopData{idx,strcmp(datanames,'DN')} = DN;
    GLMSPopData{conpanel.popidx,strcmp(datanames,'GLMP')} = GLMP;
end

% Save and display data
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
save([library 'GLMSPopData'],'GLMSPopData')


disp([char(GLMP.datafile) ' Raw Data Reorganized.'])


end

function poptuningdiffs(~,~)

% Load figure data
gaborfig = get(gcf,'userdata');
conpanel = get(gaborfig.conpanel,'userdata');
dispanel = get(gaborfig.dispanel,'userdata');

% Load differences in tuning
dispanel.oneDL = logical(dispanel.oneDL);
oneDflashPD = dispanel.flashparams(dispanel.oneDL,end-1);
oneDgaborPD = dispanel.gaborparams(dispanel.oneDL,end-1);
twoDflashPD = dispanel.flashparams(~dispanel.oneDL,end-1);
twoDgaborPD = dispanel.gaborparams(~dispanel.oneDL,end-1);

oneDdiffs = oneDflashPD - oneDgaborPD;
twoDdiffs = twoDflashPD - twoDgaborPD;

figure(100); clf;
ax = axes('parent',gcf,'units','normalized','pos',[.1 .6 .8 .3]);
cla; hold on; box on; grid on;
temp = mod(oneDdiffs,pi);
temp(temp>pi/2) = temp(temp>pi/2)-pi;
hist(temp,'parent',ax)
xlabel('Gabor - Flash PDs')
ylabel('# of neurons')
title('1D Flash vs Gabor PDs')

ax = axes('parent',gcf,'units','normalized','pos',[.1 .1 .8 .3]);
cla; hold on; box on; grid on;
temp = mod(twoDdiffs,pi/2);
temp(temp>pi/4) = temp(temp>pi/4)-(pi/2);
hist(temp,'parent',ax)
xlabel('Flash - Gabor PAs')
ylabel('# of neurons')
title('2D Flash vs Gabor PAs')


% Save figure data
set(gaborfig.dispanel,'userdata',dispanel)
set(gaborfig.conpanel,'userdata',conpanel)
set(gcf,'userdata',gaborfig)

end

function corranal(~,~)
global GLMSPopData datatypes

% Load user data
gaborfig = get(gcf,'userdata');
conpanel = get(gaborfig.conpanel,'userdata');
dispanel = get(gaborfig.dispanel,'userdata');

% % Preallocate space
% conpanel.corranal.corrcoef = nan(size(conpanel.gaboridx));
% conpanel.corranal.stim = cell(size(conpanel.gaboridx));
% conpanel.corranal.flash.allresps = cell(size(conpanel.gaboridx));
% conpanel.corranal.flash.maxlumresp = nan(size(conpanel.gaboridx));
% conpanel.corranal.flash.maxcolresp = nan(size(conpanel.gaboridx));
% conpanel.corranal.flash.mod = nan(size(conpanel.gaboridx));
% conpanel.corranal.gabor.allresps = cell(size(conpanel.gaboridx));
% conpanel.corranal.gabor.maxlumresp = nan(size(conpanel.gaboridx));
% conpanel.corranal.gabor.maxcolresp = nan(size(conpanel.gaboridx));
% conpanel.corranal.gabor.mod = nan(size(conpanel.gaboridx));
% 
% numel(conpanel.gaboridx)
% 
% for n = 1:numel(conpanel.gaboridx)
%     
%     idx = conpanel.gaboridx(n);
%     GLMP = GLMSPopData{idx+1,strcmp(datatypes,'GLMP')};
%     sub = GLMSPopData{idx+1,strcmp(datatypes,'Subunit')};
%     glmp = GLMP.subunit{sub};
%     gabor = GLMP.subunit{3};
%     axlim = max(glmp.Lcc);
%     
%     % Symmetrize responses across origin
%     glmpsymresps = nan(size(glmp.uniqueLcc));
%     for i = 1:numel(glmp.uniqueLcc)
%         stim = [glmp.uniqueLcc(i) glmp.uniqueMcc(i)];
%         symidx = all(softEq([glmp.Lcc glmp.Mcc],stim),2) | ...
%             all(softEq([glmp.Lcc glmp.Mcc],-stim),2);
%         glmpsymresps(i) = mean(glmp.fr(symidx));
%     end
%     
%     % Symmetrize responses across origin
%     gaborsymresps = nan(size(gabor.uniqueLcc));
%     for i = 1:numel(gabor.uniqueLcc)
%         stim = [gabor.uniqueLcc(i) gabor.uniqueMcc(i)];
%         symidx = all(softEq([gabor.Lcc gabor.Mcc],stim),2) | ...
%             all(softEq([gabor.Lcc gabor.Mcc],-stim),2);
%         gaborsymresps(i) = mean(gabor.fr(symidx));
%         gabor.fr(symidx)
%     end
%     
%     
%     % Find only commonly tested stimuli
%     [c,ai,bi] = intersect([glmp.uniquetheta glmp.uniquerho],[gabor.uniquetheta gabor.uniquerho],'rows');
%     glmpresps = glmpsymresps(ai);
%     gaborresps = gaborsymresps(bi);
%     
%     maxresp = max(cat(1,glmpresps,gaborresps));
%     
%     % Catalgoue resp to highest lum and col contrast
%     temp = find(c(:,1)==pi/4 | c(:,1)==-3*pi/4);
%     [~,tempidx] = max(c(temp,2));
%     lumidx = temp(tempidx);
%     
%     temp = find(c(:,1)==-pi/4 | c(:,1)==3*pi/4);
%     [~,tempidx] = max(c(temp,2));
%     colidx = temp(tempidx);
%     
%     
%     % GLMP bubble plot
% %     figure(1); clf;
% %     for i = 1:size(c,1)
% %         mn = glmpresps(i)/maxresp*30+3; hold on;
% %         h = polar(c(i,1),c(i,2),'ko');
% %         set(h,'MarkerFaceColor','none','MarkerSize',mn)
% %     end
% %     box on; grid off;
% %     axis equal square tight
% %     set(gca,'CameraPosition',[0 0 275.3526],'xlim',[-axlim axlim],'ylim',[-axlim axlim])
% % 
% %     % Gabor bubble plot
% %     figure(2); clf;
% %     for i = 1:size(c,1)
% %         mn = gaborresps(i)/maxresp*30+3; hold on;
% %         h = polar(c(i,1),c(i,2),'ko');
% %         set(h,'MarkerFaceColor','none','MarkerSize',mn)
% %     end
% %     box on; grid off;
% %     axis equal square tight
% %     set(gca,'xlim',[-axlim axlim],'ylim',[-axlim axlim])
%     
%     conpanel.corranal.corrcoef(n) = corr(glmpresps,gaborresps);
%     conpanel.corranal.stim{n} = c;
%     conpanel.corranal.flash.allresps{n} = glmpresps;
%     conpanel.corranal.flash.maxlumresp(n) = glmpresps(lumidx);
%     conpanel.corranal.flash.maxcolresp(n) = glmpresps(colidx);
%     conpanel.corranal.flash.mod(n) = max(glmpresps) - mean(glmp.blfr);
%     conpanel.corranal.gabor.allresps{n} = gaborresps;
%     conpanel.corranal.gabor.maxlumresp(n) = gaborresps(lumidx);
%     conpanel.corranal.gabor.maxcolresp(n) = gaborresps(colidx);
%     conpanel.corranal.gabor.mod(n) = max(gaborresps) - mean(gabor.blfr)
%     if conpanel.corranal.gabor.mod(n) < 10
%         keyboard
%     end
%     
%     max(conpanel.corranal.flash.allresps{n});
%     max(conpanel.corranal.gabor.allresps{n});
% end

figure(3); clf; hold on;
histogram(conpanel.corranal.corrcoef,-1:.05:1);
set(gca,'xlim',[-1 1])

lummax = max(cat(1,conpanel.corranal.flash.maxlumresp,conpanel.corranal.gabor.maxlumresp));
colmax = max(cat(1,conpanel.corranal.flash.maxcolresp,conpanel.corranal.gabor.maxcolresp));
figure(4); clf; hold on;
scatter(conpanel.corranal.flash.maxlumresp,conpanel.corranal.flash.maxcolresp,'ko')
scatter(conpanel.corranal.gabor.maxlumresp,conpanel.corranal.gabor.maxcolresp,'ro')
xlim([0 lummax])
ylim([0 colmax])
set(gca,'xscale','log','yscale','log');

% Save figure data
set(gaborfig.dispanel,'userdata',dispanel)
set(gaborfig.conpanel,'userdata',conpanel)
set(gcf,'userdata',gaborfig)

end

