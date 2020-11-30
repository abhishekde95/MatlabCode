% This is a new script where I would like to perform some additional
% analysis on Greg's old DO data as it has way more number of spikes from
% the WhiteNoise experiments
% Author - Abhishek De, 11/17, hacked most of the code from WNPop.m, section 6
% Trying to make this a GUI: might take some time but worth it 

function [] = WNColorOppAnalysis()
clearvars; close all;
global row col subunitcmap
subunitcmap = [0.89412 0.10196 0.1098;0.21569 0.49412 0.72157;0.30196 0.68627 0.2902;0.59608 0.30588 0.63922;1 0.49804 0;1 1 0.2;0.65098 0.33725 0.15686;0.96863 0.50588 0.74902;0.6 0.6 0.6];
row = 1;
col = 1;
SetupFig();
SetUpTable();
end

function SetupFig()
gcf = figure;
set(gcf,'Name','WNColorOppAnalysis');
set(gcf,'DefaultAxesUnits','pixels')
set(gcf,'position',[10 10 1050 800]);
set(gcf,'ButtonDownFcn','drawnow');  % in case user drags figure
clf;
a = get(gcf,'UserData');
a.uicontrols.WN_check = uicontrol('style','pushbutton','Callback',@calcWN_check, 'string','WN_check','Position',[500 775 70 20]);    
set(gcf,'UserData',a);
SetUpAxes();
end

function SetUpAxes()
global maxT subunitcmap
AXESWIDTH = 50;
a = get(gcf,'UserData');
h1 = []; h6 = [];
maxT = 15;
nframes = maxT;
figpos = get(gcf,'Position');
for i = 1:nframes
    h = axes('position',[(1+1/nframes)*AXESWIDTH*(i-1)+200 715 AXESWIDTH AXESWIDTH]);
    set(gca,'XTick',[],'YTick',[],'Box','on');
    axis image;
    h1 = [h1; h];  % STA checkerboard
    
    h = axes('position',[(1+1/nframes)*AXESWIDTH*(i-1)+200 615 AXESWIDTH AXESWIDTH]);
    set(gca,'XTick',[],'YTick',[],'Box','on');
    axis image;
    h6 = [h6; h];  % subunit STA
end
h2 = axes('position',[200 400 3*AXESWIDTH 3*AXESWIDTH]);set(gca,'XTick',[],'YTick',[],'Box','on'); axis image;
h3 = axes('position',[400 400 3*AXESWIDTH 3*AXESWIDTH]);set(gca,'XTick',[],'YTick',[],'Box','on'); axis image;
h4 = axes('position',[600 400 3*AXESWIDTH 3*AXESWIDTH]);set(gca,'XTick',[],'YTick',[],'Box','on'); axis image;
h5 = axes('position',[825 400 3*AXESWIDTH 3*AXESWIDTH]);set(gca,'XTick',[],'YTick',[],'Box','on'); axis image;
h7 = axes('position',[200 100 3*AXESWIDTH 3*AXESWIDTH]);set(gca,'XTick',[],'YTick',[],'Box','on'); axis image;
h8 = axes('position',[450 100 5*AXESWIDTH 5*AXESWIDTH]);set(gca,'XTick',[],'YTick',[],'Box','on'); axis image;
a.axeshandles.STA_check = h1;
a.axeshandles.STA_singleframe = h2;
a.axeshandles.basisvec1 = h3;
a.axeshandles.basisvec2 = h4;
a.axeshandles.FRmap = h5;
a.axeshandles.FRcontour = h7;
a.axeshandles.STA_subunit = h6;
a.axeshandles.RGBs = h8;
a.axeshandles.numspikes = axes('position',[600 785 1 1]);

maxsubunits = size(subunitcmap, 1);
subunit_pos = get(a.axeshandles.STA_singleframe, 'Position');
FRmap_pos = get(a.axeshandles.FRmap, 'Position');
UD = get(a.axeshandles.STA_singleframe);
UD.mask_display = zeros(10,10);
set(a.axeshandles.STA_singleframe,'UserData',UD);
axes_w = subunit_pos(3);
slider_w = axes_w; slider_h = 15;
slsubunit_pos = [subunit_pos(1) subunit_pos(2)+axes_w+5 slider_w slider_h];
updatemask_h = 22;
updatemask_pos = [subunit_pos(1) subunit_pos(2)-updatemask_h-5 axes_w updatemask_h];
FRmap_h = updatemask_h;
calcFRmap_pos = [FRmap_pos(1) FRmap_pos(2)-FRmap_h-50 axes_w updatemask_h];
text_w = 22; text_h = text_w;
clearmask_pos = [updatemask_pos(1) updatemask_pos(2)-updatemask_h-20 updatemask_pos(3:4)];
a.uicontrols.txsubunit = uicontrol('Style','edit','ForegroundColor',subunitcmap(1,:),'FontSize',12,'FontWeight','bold','String','1','Enable','inactive','Position',[slider_w/2-text_w/2+slsubunit_pos(1) slsubunit_pos(2)+slider_h+5 text_w text_h]);
a.uicontrols.slsubunit = uicontrol('Style','slider','Callback',{@update_subunit_text,a.uicontrols.txsubunit},'Min',1,'Max',9,'SliderStep',[1/(maxsubunits-1) 1/(maxsubunits-1)],'Value',1,'Position',slsubunit_pos);
a.uicontrols.updatemask = uicontrol('style','pushbutton','string','Update Mask','Callback',@update_mask,'Position', updatemask_pos);
a.uicontrols.clearmask = uicontrol('style','pushbutton','string','Clear Mask','Position',clearmask_pos,'Callback',@clear_subunit_mask);
a.uicontrols.calcFRmap = uicontrol('style','pushbutton','string','Calc FRmap','Callback',@calcFRmap,'Position',calcFRmap_pos);

a.axeshandles.label_basisvec1 = axes('position',[450 subunit_pos(2)+axes_w+15 1 1]);
a.axeshandles.label_basisvec2 = axes('position',[650 subunit_pos(2)+axes_w+15 1 1]);
a.axeshandles.label_FRmap = axes('position',[850 subunit_pos(2)+axes_w+15 1 1]);
a.axeshandles.label_FRcontour = axes('position',[230 275 1 1]);
a.axeshandles.label_RGBs = axes('position',[520 370 1 1]);
axes(a.axeshandles.label_basisvec1);text(0,0,'Basisvec 1','FontSize',10);
axes(a.axeshandles.label_basisvec2);text(0,0,'Basisvec 2','FontSize',10);
axes(a.axeshandles.label_FRmap);text(0,0,'FR map','FontSize',10);
axes(a.axeshandles.label_FRcontour);text(0,0,'Contour plot','FontSize',10);
axes(a.axeshandles.label_RGBs);text(0,0,'RGB subunits','FontSize',10);
set(gcf,'UserData',a);
drawnow;
end

function SetUpTable(~,~)
[fnames, spikenums] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\WhiteNoise\ColorOpponent.txt');
filename = [];
for ii = 1:numel(fnames)
    tmp = fnames{ii};
    filename = [filename ;cellstr(tmp{1})];
end
a = get(gcf,'UserData');
a.uitable = uitable('Parent',gcf,'Position',[20 275 150 500],'ColumnName',{'Filename'},'Data',cellstr(filename),'CellSelectionCallBack',@StoreRC,'ColumnWidth', {98,'auto','auto','auto'});
set(gcf,'UserData',a);
end

function calcWN_check(~,~)
global row col maxT stro
a = get(gcf,'UserData');
% clearing all the figures 
cla(a.axeshandles.STA_singleframe,'reset'); set(a.axeshandles.STA_singleframe,'XTick',[],'YTick',[],'Box','on');
cla(a.axeshandles.basisvec1,'reset'); set(a.axeshandles.basisvec1,'XTick',[],'YTick',[],'Box','on');
cla(a.axeshandles.basisvec2,'reset'); set(a.axeshandles.basisvec2,'XTick',[],'YTick',[],'Box','on');
cla(a.axeshandles.FRmap,'reset'); set(a.axeshandles.FRmap,'XTick',[],'YTick',[],'Box','on');
cla(a.axeshandles.FRcontour,'reset'); set(a.axeshandles.FRcontour,'XTick',[],'YTick',[],'Box','on');
cla(a.axeshandles.RGBs,'reset');set(a.axeshandles.RGBs,'XTick',[],'YTick',[],'Box','on');
cla(a.axeshandles.numspikes,'reset');


stro = nex2stro(findfile(char(a.uitable.Data(row,col))));
framerate = stro.sum.exptParams.framerate;
nstixperside = stro.sum.exptParams.nstixperside;
ntrials = length(stro.sum.absTrialNum);
stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
stimoffidx = find(strcmp(stro.sum.trialFields(1,:),'all_off'));
nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
sigmaidxs = strmatch('sigma',stro.sum.trialFields(1,:));
L = stro.trial(:,noisetypeidx) == 1;
stro.ras(~L,:) = [];
stro.trial(~L,:) = [];
out = getWhtnsStats(stro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, stro.sum.rasterCells{1});
tmpstro = [];
STAs = out{1};
nspikes = out{3};
energy = sum(STAs.^2);
whichframe = find(energy == max(energy));
muvect = reshape(repmat([.5 .5 .5],nstixperside^2,1),nstixperside^2*3,1);
UD  = get(a.axeshandles.STA_singleframe,'UserData');
UD.STAs = STAs;
UD.whichframe = whichframe;
UD.mask_display = zeros(10,10);
normfactor = 0.5/((max(abs(STAs(:))))*1.05); % This step is necessary to constrain the values within [-0.5, 0.5]
% The plotting begins here
for i = 1:size(STAs,2) % evaluates for each frame
    STA = normfactor*(STAs(:,i))+muvect; % This makes the values fall back within a range of 0 and 1.
    STA = reshape(STA,[nstixperside nstixperside 3]); % Decomposing the STA into R, G and B plane
    axes(a.axeshandles.STA_check(i)); image(STA); set(gca,'XTick',[],'YTick',[]); axis square;
    if i == whichframe
        image(STA,'Parent', a.axeshandles.STA_singleframe, 'ButtonDownFcn', @update_subunit_mask); set(gca,'XTick',[],'YTick',[]); axis square;
    end
end
set(a.axeshandles.STA_singleframe,'UserData',UD);
cla(a.axeshandles.numspikes,'reset');
axes(a.axeshandles.numspikes);text(0,0,num2str(nspikes),'FontSize',10);
set(a.axeshandles.STA_singleframe,'UserData',UD);
set(gcf,'UserData',a);
end

function StoreRC(~,evt)
a = get(gcf,'UserData');
global row col
row = evt.Indices(1);
col = evt.Indices(2);
set(gcf,'UserData',a);
end

function update_subunit_text(h, ~, htext)
global subunitcmap
n = round(get(h, 'Value'));
set(h, 'Value', n);
set(htext, 'String', n, 'ForegroundColor', subunitcmap(n,:));
end

function update_subunit_mask(himage, ~) 
axsubunit = get(himage, 'Parent');
st = get(gcf, 'SelectionType');
UD = get(axsubunit, 'UserData');
whichpt = get(axsubunit, 'CurrentPoint');
whichpt = round(whichpt(1,1:2));
whichpt = min(max(1, whichpt), 10); % a hack here

if strcmp(st, 'normal') % left click
    UD.mask_display(whichpt(2), whichpt(1)) = current_subunit; % mask gets stored the UserData of a.axeshandles.STA_singleframe
elseif strcmp(st, 'alt') % right click
    UD.mask_display(whichpt(2), whichpt(1)) = 0;
end
set(axsubunit, 'UserData', UD);
paint_subunit_selection(axsubunit);

    function n = current_subunit
        a = get(gcf, 'UserData');
        n = get(a.uicontrols.slsubunit, 'Value');
    end
end

function paint_subunit_selection(ax)
UD = get(ax, 'UserData');
paint_rectangles(UD.mask_display, ax);
set(ax, 'UserData', UD);
end

function paint_rectangles(mask, ax)
global subunitcmap
if verLessThan('matlab', '8.4.0.150421') % R2014b - Graphics changes
    HT = {'HitTest' 'off'}; % old method
else
    HT = {'PickableParts' 'none'}; % new method
end
axes(ax);
delete(findobj(get(ax, 'children'), 'type', 'rectangle'));
if any(~isnan(mask(:)))
    maxsubunits = size(subunitcmap, 1);
    for n = 1:maxsubunits
        [I,J] = find(mask == n);
        for pos = [J-.5 I-.5 ones(length(I), 2)]'
            rectangle('Position', pos, 'FaceColor', subunitcmap(n,:), HT{:});
        end
    end
end

end

function update_mask(~,~)
a = get(gcf,'UserData');
UD = get(a.axeshandles.STA_singleframe,'UserData');
mask = UD.mask_display;
STAs = UD.STAs;
whichframe = UD.whichframe;
newmask = repmat(mask(:),[3 1]); % for R, G and B channels
basisvec1 = zeros(size(STAs)); basisvec2 = zeros(size(STAs)); 
basisvec1 = STAs.*repmat(logical(newmask==1),[1 size(STAs,2)]);
basisvec2 = STAs.*repmat(logical(newmask==2),[1 size(STAs,2)]);
muvect = reshape(repmat([.5 .5 .5],size(mask,1)^2,1),size(mask,1)^2*3,1);
UD_bv1 = get(a.axeshandles.basisvec1,'UserData');
UD_bv2 = get(a.axeshandles.basisvec2,'UserData');
UD_bv1.vec = basisvec1;
UD_bv2.vec = basisvec2;


% displaying basisvec1
normfactor1 = 0.5/((max(abs(basisvec1(:,whichframe))))*1.05); % This step is necessary to constrain the values within [-0.5, 0.5]
dispbasisvec1 = normfactor1*(basisvec1(:,whichframe))+muvect; % This makes the values fall back within a range of 0 and 1.
dispbasisvec1 = reshape(dispbasisvec1,[size(mask,1) size(mask,2) 3]);
axes(a.axeshandles.basisvec1); image(dispbasisvec1);set(gca,'XTick',[],'YTick',[]);

% displaying basisvec2
normfactor2 = 0.5/((max(abs(basisvec2(:,whichframe))))*1.05); % This step is necessary to constrain the values within [-0.5, 0.5]
dispbasisvec2 = normfactor2*(basisvec2(:,whichframe))+muvect; % This makes the values fall back within a range of 0 and 1.
dispbasisvec2 = reshape(dispbasisvec2,[size(mask,1) size(mask,2) 3]);
axes(a.axeshandles.basisvec2); image(dispbasisvec2);set(gca,'XTick',[],'YTick',[]);

R1_ind = find(newmask(1:100)==1);
G1_ind = R1_ind + 100;
B1_ind = R1_ind + 200;
R2_ind = find(newmask(1:100)==2);
G2_ind = R2_ind + 100;
B2_ind = R2_ind + 200;

axes(a.axeshandles.RGBs);plot3(basisvec1(1:100,whichframe),basisvec1(101:200,whichframe),basisvec1(201:300,whichframe),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]); hold on;
plot3(mean(basisvec1(R1_ind,whichframe)),mean(basisvec1(G1_ind,whichframe)),mean(basisvec1(B1_ind,whichframe)),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
plot3(basisvec2(1:100,whichframe),basisvec2(101:200,whichframe),basisvec2(201:300,whichframe),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); 
plot3(mean(basisvec2(R2_ind,whichframe)),mean(basisvec2(G2_ind,whichframe)),mean(basisvec2(B2_ind,whichframe)),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); 
plot3([-20 20],[0 0],[0 0],'-k','Linewidth',2); % x- axis 
plot3([0 0],[-20 20],[0 0],'-k','Linewidth',2); % y- axis
plot3([0 0],[0 0],[-20 20],'-k','Linewidth',2); % z- axis
xlabel('R'),ylabel('G'),zlabel('B'); hold off;

subunitSTAs = basisvec1 + basisvec2;
normfactor = 0.5/((max(abs(subunitSTAs(:))))*1.05); % This step is necessary to constrain the values within [-0.5, 0.5]
for i = 1:size(STAs,2) % evaluates for each frame
    subunitSTA = normfactor*(subunitSTAs(:,i))+muvect; % This makes the values fall back within a range of 0 and 1.
    subunitSTA = reshape(subunitSTA,[size(mask,1) size(mask,2) 3]); % Decomposing the STA into R, G and B plane
    axes(a.axeshandles.STA_subunit(i)); image(subunitSTA); set(gca,'XTick',[],'YTick',[]); axis square;
end
set(a.axeshandles.basisvec1,'UserData', UD_bv1);
set(a.axeshandles.basisvec2,'UserData', UD_bv2);
set(gcf,'UserData',a);
end

function calcFRmap(~,~)
% This is the part where I project all the stimuli onto the basisvectors
global stro maxT
a = get(gcf,'UserData');
spikename = getSpikenum(stro);
spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
framerate = stro.sum.exptParams.framerate;
nstixperside = stro.sum.exptParams.nstixperside;
ntrials = size(stro.trial,1);
stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
stimoffidx = find(strcmp(stro.sum.trialFields(1,:),'all_off'));
nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
sigmaidxs = strmatch('sigma',stro.sum.trialFields(1,:));
seedidx = strcmp(stro.sum.trialFields(1,:),'seed'); %get seed index from trialFields
muidxs = [find(strcmp(stro.sum.trialFields(1,:),'mu1')),... %get mu indices into vector from trialFields
         find(strcmp(stro.sum.trialFields(1,:),'mu2')),...
         find(strcmp(stro.sum.trialFields(1,:),'mu3'))];
sigmaidxs = [find(strcmp(stro.sum.trialFields(1,:),'sigma1')),... %get sigma indices into vector from trialFields
            find(strcmp(stro.sum.trialFields(1,:),'sigma2')),...
            find(strcmp(stro.sum.trialFields(1,:),'sigma3'))];
nrandnums_perchannel = nstixperside^2;
gammaTable = stro.sum.exptParams.gamma_table; %get gamma_table
gammaTable = reshape(gammaTable,length(gammaTable)/3,3); %reshapse gamma_table into three columns
invgammaTable = InvertGamma(gammaTable,1); %invert gamma_table
ngammasteps = size(invgammaTable,1); %get number of rows of gamma_table (65536)
xx = linspace(stro.sum.exptParams.gauss_locut/1000, stro.sum.exptParams.gauss_hicut/1000, ...
    ngammasteps); %dividing gauss by 1000 and making equal intervals so that there are 65536 values
yy = norminv(xx');
msperframe = 1000/stro.sum.exptParams.framerate;
L = stro.trial(:,noisetypeidx) == 1;
stro.ras(~L,:) = [];
stro.trial(~L,:) = [];
UD_bv1 = get(a.axeshandles.basisvec1,'UserData');
UD_bv2 = get(a.axeshandles.basisvec2,'UserData');
vec1 = UD_bv1.vec;
vec1_rev = flipdim(vec1,2);
vec2 = UD_bv2.vec;
vec2_rev = flipdim(vec2,2);
initargs = {[vec1_rev(:), vec2_rev(:)] , 0, sum(stro.trial(:,nframesidx)), [nrandnums_perchannel 3 maxT]};
STPROJmod('init',initargs); % initialising the STPROJmod
for k = 1:ntrials-10
    nframes = stro.trial(k,nframesidx);
    if (nframes == 0)
        continue;
    end
    seed = stro.trial(k,seedidx);
    mu = stro.trial(k,muidxs)/1000;
    sigma = stro.trial(k,sigmaidxs)/1000;
    invnormcdf = bsxfun(@plus, bsxfun(@times, yy, sigma), mu);
    randnums = getEJrandnums(3*nrandnums_perchannel*nframes, seed); % random numbers for 300 x 9 pixels
    randnums = reshape(randnums, 3*nrandnums_perchannel, nframes);
    for gun = 1:3
        idxs = (1:nrandnums_perchannel)+nrandnums_perchannel*(gun-1);
        randnums(idxs,:) = reshape(invnormcdf(randnums(idxs,:)+1,gun),[length(idxs) nframes]);
    end
    rgbs = randnums;
    t_stimon = stro.trial(k, stimonidx);
    spiketimes = (stro.ras{k,spikeidx}-t_stimon)*1000; % observing spiketimes in milliseconds
    frametimes = linspace(0, nframes*msperframe, nframes)+(msperframe/2)';
    % Deleting the spikes taking place before the first 9 frames as I need to look at the 9 preceding frames
    spiketimes(spiketimes < maxT*msperframe) = [];
    % Deleting the spikes that take place after the stimulus was
    % removed as it would imply that these spikes do not result from
    % the stimulus shown on the screen
    spiketimes(spiketimes > frametimes(end)) = [];
    n = hist(spiketimes, frametimes);
%     keyboard;
    STPROJmod(rgbs(:),n);
end
out = STPROJmod('return');
projs = out{1};
Lspike = out{2};
clear STPROJmod;
clear out;
X_label = 'Basisvec1';
Y_label = 'Basisvec2';
min_val = floor(min(projs(:))*100)/100;
max_val = ceil(max(projs(:))*100)/100;
num_bins = 15;
bin_interval = (max_val-min_val)/num_bins;
nbins1 = min_val:bin_interval:max_val;
% choose the bins such that they lie between 5 and 95 percentile
nbins = linspace(prctile(projs(:),5), prctile(projs(:),95),numel(nbins1)+2);
new_nbins = [nbins(1)-mean(diff(nbins)), nbins, mean(diff(nbins)) + nbins(end)];
[n_spike,out_spike] = hist3([projs(Lspike>0,1), projs(Lspike>0,2)],{new_nbins,new_nbins});
[n_raw,out_raw] = hist3([projs(:,1), projs(:,2)],{new_nbins,new_nbins});
n_spike = n_spike(2:end-1,2:end-1);
n_raw = n_raw(2:end-1,2:end-1);
non_lin = n_spike./n_raw;
non_lin_pa = padarray(non_lin,[2 2],'replicate'); % pad the array with the border elements
filt = fspecial('gaussian',5,1.0); % building a gaussian filter
non_lin_blurred = conv2(non_lin_pa,filt,'same'); % convolving a gaussian filter with the firing rate map
non_lin_blurred = non_lin_blurred(3:end-2,3:end-2);
xmin = min(nbins); xmax = max(nbins);
axes(a.axeshandles.FRmap);imagesc([xmin xmax],[xmin xmax],non_lin); xlabel(Y_label), ylabel(X_label); 
axes(a.axeshandles.FRcontour);contour(nbins',nbins,non_lin_blurred);xlabel(Y_label), ylabel(X_label);set(gca,'YDir','Reverse');
set(gcf,'UserData',a);

end

function clear_subunit_mask(~,~)
a = get(gcf,'UserData');
UD = get(a.axeshandles.STA_singleframe,'UserData');
UD.mask_display = zeros(10,10);
cla(a.axeshandles.basisvec1,'reset'); set(a.axeshandles.basisvec1,'XTick',[],'YTick',[],'Box','on');
cla(a.axeshandles.basisvec2,'reset'); set(a.axeshandles.basisvec2,'XTick',[],'YTick',[],'Box','on');
cla(a.axeshandles.FRmap,'reset'); set(a.axeshandles.FRmap,'XTick',[],'YTick',[],'Box','on');
cla(a.axeshandles.FRcontour,'reset'); set(a.axeshandles.FRcontour,'XTick',[],'YTick',[],'Box','on');
cla(a.axeshandles.RGBs,'reset');set(a.axeshandles.RGBs,'XTick',[],'YTick',[],'Box','on');
set(a.axeshandles.STA_singleframe,'UserData',UD);
set(gcf,'UserData',a);
paint_subunit_selection(a.axeshandles.STA_singleframe)

end