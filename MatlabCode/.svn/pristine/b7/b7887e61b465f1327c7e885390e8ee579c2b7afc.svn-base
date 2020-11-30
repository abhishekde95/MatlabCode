function [] = WhiteNoisethreshGUI()
% Pilot script for better visualisation of WhiteNoisethresh Data
%close all; 
clearvars;
global row col
row = 1;
col = 1;
SetupFig();
SetUpTable();
end

function Parse_stro(stro)
global spikename maskidx spikeidx neurothreshidx nstixperside ngammasteps seedidx nframesidx correctidx
global fponidx stimoffidx stimonidx muidxs sigmaidxs basisvecidx weightsidx fpacqidx targetspikerateidx basisvecdiridx latencyidx reversalflagidx parentverticesidx
global msperframe ntrials maxT xx yy M linepredtol stepsizescale stepsize nreversals oogscale
global gl
spikename = 'sig001a';%getSpikenum(stro);
maskidx = strcmp(stro.sum.rasterCells(1,:), 'subunit_mask');
spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
basisvecidx = strcmp(stro.sum.rasterCells(1,:),'basis_vec');
weightsidx = strcmp(stro.sum.rasterCells(1,:),'weights');
parentverticesidx = strcmp(stro.sum.rasterCells(1,:),'parentvertices');
nstixperside = stro.sum.exptParams.nstixperside;
ngammasteps = 2^16; % 65536
linepredtol = stro.sum.exptParams.linepredtol;
stepsizescale = stro.sum.exptParams.stepsizescale;
stepsize = stro.sum.exptParams.stepsize;
nreversals = stro.sum.exptParams.nreversals;
oogscale = stro.sum.exptParams.oogscale;
seedidx = strcmp(stro.sum.trialFields(1,:),'seed');
nframesidx = strcmp(stro.sum.trialFields(1,:),'num_frames');
stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
stimoffidx = strcmp(stro.sum.trialFields(1,:),'stim_off');
fponidx = strcmp(stro.sum.trialFields(1,:),'fp_on');
fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
basisvecdiridx = strcmp(stro.sum.trialFields(1,:),'weights_idx');
neurothreshidx = strcmp(stro.sum.trialFields(1,:),'neurothresh'); % when exactly the neurothresh trials started
targetspikerateidx = strcmp(stro.sum.trialFields(1,:),'targetspikerate');
correctidx = strcmp(stro.sum.trialFields(1,:),'correct');
muidxs = [find(strcmp(stro.sum.trialFields(1,:),'mu1')), ...
    find(strcmp(stro.sum.trialFields(1,:),'mu2')), ...
    find(strcmp(stro.sum.trialFields(1,:),'mu3'))];
sigmaidxs = [find(strcmp(stro.sum.trialFields(1,:),'sigma1')), ...
    find(strcmp(stro.sum.trialFields(1,:),'sigma2')), ...
    find(strcmp(stro.sum.trialFields(1,:),'sigma3'))];
latencyidx = strcmp(stro.sum.trialFields(1,:),'latency');
reversalflagidx = strcmp(stro.sum.trialFields(1,:),'reversalflag');

msperframe = 1000/stro.sum.exptParams.framerate;
ntrials = size(stro.trial,1);
maxT = 12; % this represents the temporal part in the spatiotemporal receptive field
xx = linspace(stro.sum.exptParams.gauss_locut/1000, stro.sum.exptParams.gauss_hicut/1000,ngammasteps); % xx represents the probabilities. For more info, have a look at the MATLAB 'norminv' function.
yy = norminv(xx'); % defining norminv to extract the values for which the cdf values range between gauss_locut and gauss_hicut

fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S 
mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
M = inv(M');

mask_changes = [2];
all_masks = stro.ras(:,maskidx);
Fx = @(xi) any(isnan(xi)); % function that finds 'NaN' in a cell array
inds = find(cellfun(Fx,stro.ras(:,basisvecidx))==0);
if isempty(inds)
    inds = size(stro.trial,1)-1;
end
last_wntrial =  inds(1)-1;
for k = 3:last_wntrial
    if isequal(all_masks{k}, all_masks{k-1}) %|| all(all_masks{k} == 0) && any(isnan(all_masks{k-1}))
        continue
    else
        mask_changes = [mask_changes k-1 k]; %#ok<AGROW>
    end
end
if mask_changes(end) == last_wntrial 
    mask_changes(end) = [];
else
    mask_changes = [mask_changes  last_wntrial];
end
mask_changes = reshape(mask_changes , 2, []);

gl.mask_changes = mask_changes;
gl.stro = stro;
gl.inds = inds;
end


% Function for setting initial figures
function SetupFig()
gcf = figure;
set(gcf,'Name','WhiteNoisethreshGUI');
set(gcf,'DefaultAxesUnits','pixels')
set(gcf,'position',[10 10 1400 800]);
set(gcf,'ButtonDownFcn','drawnow');  % in case user drags figure
clf;
a = get(gcf,'UserData');
a.uicontrols.WN_check = uicontrol('style','pushbutton','Callback',@calcWN_check, 'string','WN_check','Position',[500 775 70 20]);
a.uicontrols.WN_check = uicontrol('style','pushbutton','Callback',@calcWN_checksigmap, 'string','Sigmap','Position',[750 775 70 20]);
a.uicontrols.WN_subunit = uicontrol('style','pushbutton','Callback',@calcWN_subunit, 'string','WN_subunit','Position',[500 605 70 20]);
a.uicontrols.NT_iso_resp = uicontrol('style','pushbutton','Callback',@calcNT_isoresp, 'string','NT_iso_resp','Position',[450 420 70 20]);
a.uicontrols.basis_vec = uicontrol('style','pushbutton','Callback',@calcbasis_vec, 'string','basis_vec','Position',[150 420 70 20]);
a.uicontrols.execute_file = uicontrol('style','pushbutton','Callback',@execute_file, 'string','Execute_file','Position',[50 420 70 20]);
a.uicontrols.logrvstheta = uicontrol('style','pushbutton','Callback',@logrvstheta, 'string','LogRvsTheta','Position',[915 420 70 20]);
a.uicontrols.BaselineFR = uicontrol('style','pushbutton','Callback',@BaselineFR, 'string','BaselineFR','Position',[1270 480 70 20]);
a.uicontrols.NTmode = uicontrol('style','popupmenu','string',{'subunit','STAvsPC1','completelist','DO','Lum','SO','htc','isorespfiles'},'Position',[20 775 150 20],'Callback',@Loadfilenames);     
set(gcf,'UserData',a);
SetUpAxes();
end

function SetUpAxes()
AXESWIDTH = 50;
a = get(gcf,'UserData');
nframes = 12;
h1 = []; h2 = []; h3 = []; h4 = [];
figpos = get(gcf,'Position');
for i = 1:nframes
    h = axes('position',[(1+1/nframes)*AXESWIDTH*(i-1)+200 715 AXESWIDTH AXESWIDTH]);
    set(gca,'XTick',[],'YTick',[],'Box','on');
    axis image;
    h1 = [h1; h];  % STA checkerboard
    
    h = axes('position',[(1+1/nframes)*AXESWIDTH*(i-1)+200 645 AXESWIDTH AXESWIDTH]);
    set(h,'XTick',[],'YTick',[],'Box','on');
    axis image;
    h2 = [h2; h];  % PC checkerboard
    
    h = axes('position',[(1+1/nframes)*AXESWIDTH*(i-1)+200 545 AXESWIDTH AXESWIDTH]);
    set(gca,'XTick',[],'YTick',[],'Box','on');
    axis image;
    h3 = [h3; h];  % STA subunit
    
    h = axes('position',[(1+1/nframes)*AXESWIDTH*(i-1)+200 475 AXESWIDTH AXESWIDTH]);
    set(h,'XTick',[],'YTick',[],'Box','on');
    axis image;
    h4 = [h4; h];  % PC subunit
end
a.axeshandles.STA_check = h1;
a.axeshandles.PC_check = h2;
a.axeshandles.STA_subunit = h3;
a.axeshandles.PC_subunit = h4;
a.axeshandles.basisvec1 = axes('position',[40 250 3*AXESWIDTH 3*AXESWIDTH]); set(gca,'XTick',[],'YTick',[]);
a.axeshandles.basisvec2 = axes('position',[40 50 3*AXESWIDTH 3*AXESWIDTH]); set(gca,'XTick',[],'YTick',[]);
a.axeshandles.display_basisvec_combo = axes('position',[665 250 3*AXESWIDTH 3*AXESWIDTH]); set(gca,'XTick',[],'YTick',[]); axis image;
a.axeshandles.display_basisvec_subunitRGB = axes('position',[665 50 3*AXESWIDTH 3*AXESWIDTH]); set(gca,'XTick',[],'YTick',[]); axis image; 

h5 = []; h6 = []; h7 = []; h8 = []; h9 = []; h10 = []; h11 = [];

for j = 1:2 % At max for 2 different firing rates
    h = axes('position',[250 250-(200*(j-1)) 3*AXESWIDTH 3*AXESWIDTH]); set(h,'XTick',[],'YTick',[]);
    axis image;
    h5 = [h5; h];  % TFR dirs
    
    h = axes('position',[460 250-(200*(j-1)) 3*AXESWIDTH 3*AXESWIDTH]); set(h,'XTick',[],'YTick',[]);
    axis image;
    h6 = [h6; h];  % iso resp contours
 
    % For displaying gun signals for the subunits (STA and PC1)
    h = axes('position',[885 520+135*(j-1) 4.5*AXESWIDTH 2*AXESWIDTH]); set(h,'XTick',[],'YTick',[]);
    h7 = [h7; h]; % STA subunits
    
    h = axes('position',[1160 520+135*(j-1) 4.5*AXESWIDTH 2*AXESWIDTH]); set(h,'XTick',[],'YTick',[]);
    h8 = [h8; h]; % PC1 subunits

     % At max for 2 different firing rates
    h = axes('position',[870 250-(200*(j-1)) 3*AXESWIDTH 3*AXESWIDTH]); set(h,'XTick',[],'YTick',[]);
    axis image;
    h9 = [h9; h];  % Log R vs Theta

    % WN and NT Baseline FR raster plots
    h = axes('position',[1250 50+220*(j-1) 2.5*AXESWIDTH 3.5*AXESWIDTH]); set(h,'XTick',[],'YTick',[]);
    h10 = [h10; h];  % Baseline Raster plots
    
    % 2-D FR (firing rate) map for as a function of contrast for 2 TFRs  
    h = axes('position',[1060 250-(200*(j-1)) 3*AXESWIDTH 3*AXESWIDTH]); set(h,'XTick',[],'YTick',[]);
    axis image;
    h11 = [h11; h];  
end

a.axeshandles.TFR_dirs = h5;    
a.axeshandles.TFR_isorespsurface = h6;
a.axeshandles.STA_subunits_temp = h7;
a.axeshandles.PC_subunits_temp = h8;
a.axeshandles.LogRvsTheta = h9;
a.axeshandles.BaselineFR = h10;
a.axeshandles.Baselinespikerate = h11;
a.axeshandles.text_STA_subunits = axes('position',[1000 790 1 1]); text(0,0,['STA']);
a.axeshandles.text_PC_subunits = axes('position',[1200 790 1 1]); text(0,0,['PC']);
a.axeshandles.WN_Baselinerasters = axes('position',[1290 250 1 1]); text(0,0,['WN']);
a.axeshandles.NT_Baselinerasters = axes('position',[1290 460 1 1]); text(0,0,['NT']);
a.axeshandles.file_comments = axes('Position',[280 460 1 1]); 
a.axeshandles.filename = axes('position',[600 785 1 1]);
set(gcf,'UserData',a);
drawnow;
end

function SetUpTable(~,~)
global comments NTmode filename
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM WNthresh');
comments = fetch(conn,'SELECT comments FROM WNthresh');
NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
close(conn);
a = get(gcf,'UserData');
a.uitable = uitable('Parent',gcf,'Position',[20 475 150 300],'ColumnName',{'Filename'},'Data',cellstr(filename),'CellSelectionCallBack',@StoreRC,'ColumnWidth', {98,'auto','auto','auto'});
set(gcf,'UserData',a);
end

function calcWN_check(~,~)
global gl 
global maskidx spikeidx nstixperside seedidx nframesidx stimonidx muidxs sigmaidxs msperframe maxT yy
global STAscheck STCscheck nspikescheck
stro = gl.stro;
a = get(gcf,'UserData');
mask_changes = gl.mask_changes(:,1);
for mask_span = mask_changes
    STCOVmex('init', {nstixperside^2 3 maxT});
    for k = mask_span(1):mask_span(2)
        nframes = stro.trial(k,nframesidx);
        if nframes == 0, continue; end
        
        seed = stro.trial(k,seedidx);
        mu = stro.trial(k,muidxs)/1000;
        sigma = stro.trial(k,sigmaidxs)/1000;
        
        org_mask = stro.ras{k,maskidx}; % useful for subunits computation
        nrandnums_perchannel = nstixperside^2; % In this case, it is the 100 pixels which are flickering when no subunits are selected
        % assuming Gaussian gun noise only, random number generator
        % routine as a mexfile (getEJrandnums.mexw64)
        invnormcdf = bsxfun(@plus, bsxfun(@times, yy, sigma), mu);
        %         figure(3),plot(invnormcdf);
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
        spiketimes(spiketimes < maxT*msperframe) = [];
        spiketimes(spiketimes > frametimes(end)) = [];
        n = hist(spiketimes, frametimes);
        STCOVmex(rgbs(:),n);
    end
    out = STCOVmex('return');
    STS = out{1};
    STCross = out{2};
    nspikes = out{3};
    clear STCOVmex;
    clear out;
    STAs = STS/nspikes;
    STCs = zeros(size(STCross));
    for i = 1:maxT
        tmp = STS(:,i)*STS(:,i)';
        STCs(:,i) = (nspikes.*STCross(:,i)-tmp(:))/(nspikes*(nspikes-1));
    end
    STAscheck = STAs;
    STCscheck = STCs;
    nspikescheck = nspikes;
    an_mask = [];
    an_mask = zeros(nstixperside, nstixperside);   % Pixel mask.  Someday make a tool to make this non-zero.
    Lmask = logical(repmat(~an_mask(:),[3 1]));
    PCs = [];
    for i = 1:size(STCs,2)
        STC = reshape(STCs(:,i), 3*nstixperside^2, 3*nstixperside^2);
        subSTC = STC(Lmask, Lmask);
        subSTA = STAs(Lmask,i);
        % subtracting the samples from the STA to ensure the PCs are
        % orthogonal to the the STA.
        P = eye(size(subSTC)) - subSTA*inv(subSTA'*subSTA)*subSTA';
        subSTC = P*subSTC*P';
        [tmp,d] = eig(subSTC);
        v = repmat(double(Lmask),[1 size(tmp,2)]);
        v(Lmask,:) = tmp;
        [~, idxs] = sort(diag(d));
        v = v(:,idxs);
        v = v(:,end);  % Collecting the first 3 principle components
        PCs = cat(3, PCs, v);
    end
    muvect = reshape(repmat([.5 .5 .5],nstixperside^2,1),nstixperside^2*3,1); % creating a 300 x 1 array with each entry as 0.5
    
    % The plotting begins here
    for i = 1:size(STAs,2) % evaluates for each frame
        normfactor = 0.5/((max(abs(STAs(:))))*1.05); % This step is necessary to constrain the values within [-0.5, 0.5]
        STA = normfactor*(STAs(:,i))+muvect; % This makes the values fall back within a range of 0 and 1.
        STA = reshape(STA,[nstixperside nstixperside 3]); % Decomposing the STA into R, G and B plane
        axes(a.axeshandles.STA_check(i)); image(STA); set(gca,'XTick',[],'YTick',[]); axis square;
        PC_int = squeeze(PCs(:,1,i));
        PC_int = 0.5*(PC_int)/(max(abs(PC_int)));
        normfactor = 1;
        PC = normfactor*PC_int+muvect;
        PC = reshape(PC,[nstixperside nstixperside 3]);
        axes(a.axeshandles.PC_check(i)); image(PC); set(gca,'XTick',[],'YTick',[]); axis square;
    end   
end

set(gcf,'UserData',a);
end

function calcWN_checksigmap(~,~)
global gl nstixperside 
global STAscheck STCscheck nspikescheck
a = get(gcf,'UserData');
stro = gl.stro;
STAs = STAscheck;
STCs = STCscheck;
nspikes = nspikescheck;

noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
L = stro.trial(:,noisetypeidx)==1;
mu1idx = find(strcmp(stro.sum.trialFields(1,:),'mu1'));
mu2idx = find(strcmp(stro.sum.trialFields(1,:),'mu2'));
mu3idx = find(strcmp(stro.sum.trialFields(1,:),'mu3'));
sigma1idx = find(strcmp(stro.sum.trialFields(1,:),'sigma1'));
sigma2idx = find(strcmp(stro.sum.trialFields(1,:),'sigma2'));
sigma3idx = find(strcmp(stro.sum.trialFields(1,:),'sigma3'));
muvect = unique(stro.trial(L,[mu1idx mu2idx mu3idx]),'rows')/1000;
sigmavect = unique(stro.trial(L,[sigma1idx sigma2idx sigma3idx]),'rows')/1000;
sigmavect(all(any(sigmavect == 0),2),:) = [];
gausslims = [stro.sum.exptParams.gauss_locut stro.sum.exptParams.gauss_hicut]/1000;
mumat = repmat(reshape(repmat(muvect,nstixperside^2,1),[nstixperside^2*3, 1]),[1,size(STAs,2)]);
sigmamat = repmat(reshape(repmat(sigmavect,nstixperside^2,1),[nstixperside^2* 3, 1]),[1,size(STAs,2)]);
zscoremeans = (STAs-mumat)./(sigmamat/sqrt(nspikes));
% Doing the calculations for the variances
% Only considering one correction factor per dimension (assuming variances on green and blue guns are same as red gun)

NPOINTS = 65536;
x = linspace(gausslims(1),gausslims(2),NPOINTS);
Fx = norminv(x)*sigmavect(1);
sigmacorrectionfactor = std(Fx)./sigmavect(1);
for i = 1:size(STCs,2)
    STC = reshape(STCs(:,i),[sqrt(length(STCs(:,i))),sqrt(length(STCs(:,i)))]);
    STVs(:,i) = diag(STC);
end
muvar = (sigmavect(1)*sigmacorrectionfactor)^2;
sigmavar = muvar*sqrt(2/nspikes); % Check out this page for more info https://groups.google.com/forum/#!topic/sci.stat.math/dsgmWBLJoHc
zscorevars = (STVs-muvar)./sigmavar;
maxzscoremeans = max(abs(zscoremeans(:)));
maxzscorevars = max(abs(zscorevars(:)));
alpha = 0.025;
crit = norminv(1-alpha,0,1);

% Plotting
FigHandle = figure(3);
set(FigHandle, 'Position', [50, 70, 1000, 505]);
for i = 1:size(STAs,2)
    % First STAs
    zmat = reshape(zscoremeans(:,i),[nstixperside nstixperside 3]);
    for j = 1:3
        pmat = logical(abs(zmat(:,:,j))>crit);
        im = zmat(:,:,j)./(2*maxzscoremeans)+.5; % ensures that all values lie between 0 and 1
        im = repmat(im,[1 1 3]);
        sigidxs = find(pmat);
        im(sigidxs) = .5;  % red to .5 where sig.  Looks red on dark and cyan on bright.
        subplot(6,size(STAs,2),(j-1)*size(STAs,2)+i);
        image(im);
        axis image;
        if (i == 1)
            if (j == 1)
                ylabel('STA-R');
            elseif (j == 2)
                ylabel('STA-G');
            elseif (j == 3)
                ylabel('STA-B');
            end
        end
        set(gca,'XTick',[],'YTick',[]);
    end
    % Then STVs
    zmat = reshape(zscorevars(:,i),[nstixperside nstixperside 3]);
    for j = 4:6
        pmat = logical(abs(zmat(:,:,j-3))>crit);
        im = zmat(:,:,j-3)./(2*maxzscorevars)+.5; % ensures that all values lie between 0 and 1
        im = repmat(im,[1 1 3]);
        sigidxs = find(pmat);
        im(sigidxs) = .5;
        subplot(6,size(STCs,2),(j-1)*size(STCs,2)+i);
        image(im);
        axis image;
        if (i == 1)
            if (j == 4)
                ylabel('STV-R');
            elseif (j == 5)
                ylabel('STV-G');
            elseif (j == 6)
                ylabel('STV-B');
            end
        end
        set(gca,'XTick',[],'YTick',[]);
    end
end
set(gcf,'UserData',a);
end

function calcWN_subunit(~,~)
global gl maskidx maxT seedidx muidxs sigmaidxs nstixperside nframesidx yy stimonidx spikeidx msperframe
a = get(gcf,'UserData');
stro = gl.stro;
mask_span = gl.mask_changes(:,2);
        
% Just store enough space to accomodate the 9 frames of the basis vector
st_mask = stro.ras{mask_span(1),maskidx}; % subunit mask
st_mask(st_mask == 0) = Inf;
[stIdxs,~,~] = unique(st_mask); % now the Infs map to nsubunits+1
num_subunits = length(stIdxs)-any(isinf(stIdxs)); % nsubunits, like subunits A and B
STCOV_st('init', {num_subunits 3 maxT});
for k = mask_span(1):mask_span(2)
    nframes = stro.trial(k,nframesidx);
    if (nframes == 0)
        continue;
    end
    seed = stro.trial(k,seedidx);
    mu = stro.trial(k,muidxs)/1000;
    sigma = stro.trial(k,sigmaidxs)/1000;
    % org_mask tells u if u have updated the mask or not. If org_mask is non-zero it means at this particular trial
    % u have selected the subunits and need to analyse its computation
    org_mask = stro.ras{k,maskidx};
    if any(org_mask)
        org_mask(org_mask == 0) = Inf;
        [subunitIdxs,~,mask] = unique(org_mask); % now the Infs map to nsubunits+1
        nrandnums_perchannel = length(subunitIdxs)-any(isinf(subunitIdxs)); % nsubunits, like subunits A and B
        mask = [mask; mask+max(mask); mask+2*max(mask)]; %#ok<AGROW>
    else
        nrandnums_perchannel = nstixperside^2; % In this case, it is the 100 pixels which are flickering when no subunits are selected
    end
    
    % assuming Gaussian gun noise only, random number generator
    % routine as a mexfile (getEJrandnums.mexw64)
    invnormcdf = bsxfun(@plus, bsxfun(@times, yy, sigma), mu);
    randnums = getEJrandnums(3*nrandnums_perchannel*nframes, seed); % random numbers for 300 x 9 pixels
    % This is the extracted colors for subunits/pixels using the seed number
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
    STCOV_st(rgbs(:),n);
end

out = STCOV_st('return'); % returns the covariance matrix on frame by frame basis
STS = out{1};  % A (dimension) x 9(frames) matrix
STCross = out{2};  % A (dimension x frames)x (dimension x frames)  matrix
nspikes = out{3}; % Number of spikes in the given file
clear STCOV_st;
clear out;
% Coverting the STS and the STCross into STA and STC respectively
STAs = STS/nspikes;
tmp = STS(:)*STS(:)';
STCs = (nspikes.*STCross-tmp)/(nspikes*(nspikes-1));

% Obtaining the eigenvectors and their corresponding eigenvalues  
% subtracting the samples from the STA to ensure the PCs are
% orthogonal to the the STA.
P = eye(size(STCs)) - STAs(:)*inv(STAs(:)'*STAs(:))*STAs(:)'; % WHAT DOES THIS LINE MEAN
STCs = P*STCs*P';
[tmp,d] = eig(STCs);
eig_PC = sort(diag(d)); % storing all the eigenvalues
v = real(tmp);
[~, idxs] = sort(diag(d));
v = v(:,idxs);
suppresive_PC = 2;
PCs = v(:,end);  % Collecting the first principle component
% -----------------------------------------------------------------------------------------------------
% Flipping the STAs and the PCs such that the last frame appears first and the first frame appears last
tmp1 = [];
for i=1:size(PCs,2)
    tmp2 = fliplr(reshape(PCs(:,i),[nrandnums_perchannel*3 maxT]));
    tmp1 = [tmp1, tmp2(:)];
end
PCs = tmp1;
STAs = fliplr(STAs); 
% -----------------------------------------------------------------------------------------------------

muvect = repmat([.5 ;.5; .5],nstixperside^2,1); % creating a 300 x 1 array with each entry as 0.5
% The plotting begins here
PC_int = PCs(:,1,i);
PC_int = 0.5*(PC_int)/(max(abs(PC_int))*1.05);
PC_int = reshape(PC_int,[nrandnums_perchannel*3 maxT]);
for i = 1:size(STAs,2) % evaluates for each frame
    normfactor = 0.5/((max(abs(STAs(:))))*1.05); % This step is necessary to constrain the values within [-0.5, 0.5]
    STA = STAs(:,i);
    STA = expand_vector(STA,nrandnums_perchannel,mask,1);
    STA = normfactor*(STA)+muvect;  % This makes the values fall back within a range of 0 and 1.
    STA = reshape(STA,[nstixperside nstixperside 3]); % Decomposing the STA into R, G and B plane
    axes(a.axeshandles.STA_subunit(i)); image(STA); set(gca,'XTick',[],'YTick',[]); axis square;
    
    PC = expand_vector(PC_int(:,i),nrandnums_perchannel,mask,1);
    normfactor = 1;
    PC = normfactor*PC+muvect;
    PC = reshape(PC,[nstixperside nstixperside 3]);
    axes(a.axeshandles.PC_subunit(i));image(PC); set(gca,'XTick',[],'YTick',[]); axis square;
end

% Need to plot the temporal variation of subunits - STA and PC1 
subunits = nrandnums_perchannel;
temp_plot = zeros(3,maxT,subunits);
vectors = [STAs(:), PCs(:)];
tmp_axes = [a.axeshandles.STA_subunits_temp a.axeshandles.PC_subunits_temp];
for j = 1:size(vectors,2)
    vec_temp = reshape(vectors(:,j),[nrandnums_perchannel*3 maxT]);
    for i = 1:subunits
        temp_plot(1,:,i) = vec_temp(i,:); % R
        temp_plot(2,:,i) = vec_temp(subunits+i,:); % G
        temp_plot(3,:,i) = vec_temp(2*subunits+i,:); % B
        axes(tmp_axes(i,j));
        plot(temp_plot(1,:,i),'r','LineWidth',2); hold on;
        plot(temp_plot(2,:,i),'g','LineWidth',2); hold on;
        plot(temp_plot(3,:,i),'LineWidth',2);hold off;
    end
    clear vec_temp
end
set(gcf,'UserData',a);
gl.mask = mask;

end

function calcNT_isoresp(~,~)
global gl neurothreshidx targetspikerateidx correctidx basisvecdiridx weightsidx stimonidx stimoffidx latencyidx
stro = gl.stro;
inds = gl.inds;
basisvec = gl.basisvec;
bkgnd_monitor = gl.bkgnd_monitor;
a = get(gcf,'UserData');
t_offset = stro.trial(end,latencyidx)/1000;
neurothreshmode = stro.trial(:,neurothreshidx);
basisvec_dropidx = inds(end); 
neurothresh_startidx = max(find(neurothreshmode==1,1),basisvec_dropidx+1);
num_targetspikerates = unique(stro.trial(neurothresh_startidx:end,targetspikerateidx));
norms = cell(1,numel(num_targetspikerates));
completed_search_alongdir = cell(1,numel(num_targetspikerates));
for jj = 1: numel(num_targetspikerates)
    idxs = find(~isnan(stro.trial(:,correctidx)) & stro.trial(:,targetspikerateidx)==num_targetspikerates(jj));
    idxs(idxs<=neurothresh_startidx) = [];
    different_weights = unique(stro.trial(idxs,basisvecdiridx));
    tmp_norm = [];
    tmp_completed_search_alongdir = [];

    for kk = 1:numel(different_weights)
        idxs1 = find(stro.trial(:,basisvecdiridx) == different_weights(kk));
        idxs1(idxs1<neurothresh_startidx) = [];
        raster_data = stro.ras(idxs1,1);
        tmp_norm = [tmp_norm; stro.ras{idxs1(end),weightsidx}'];
        for ii = 1:size(raster_data,1)
            tmp = raster_data{ii} ;
            spikes = tmp(tmp>stro.trial(idxs1(ii),stimonidx)+t_offset & tmp < stro.trial(idxs1(ii),stimoffidx));
            spikes = spikes - stro.trial(idxs1(ii),stimonidx)-t_offset;
        end
        [flag, gamutViolation] = Check_ifcompletely_probed(stro,basisvec,bkgnd_monitor,idxs1,tmp_norm(end,:));
        tmp_completed_search_alongdir = [tmp_completed_search_alongdir; flag gamutViolation];
    end
    norms{jj} = tmp_norm;
    completed_search_alongdir{jj} = tmp_completed_search_alongdir;
end

% Plotting the end norm or contrast values for each search direction and Converting the end norms into polar coordinates
% Need to write a small function to check for Gamut Violation based on the
% reversalflagidx
basisvec1 = basisvec{1}-bkgnd_monitor;
basisvec2 = basisvec{2}-bkgnd_monitor;
oog1 = min(abs((1-bkgnd_monitor(:))./basisvec1(:)));
oog2 = min(abs((1-bkgnd_monitor(:))./basisvec2(:)));

% keyboard;
lo = -5.0; hi = 5.0;
color = ['g','k'];
for ii = 1:size(norms,2)
    tmp = norms{ii};
    completed_dir = completed_search_alongdir{ii};
    probed_dirs = logical(completed_dir(:,1)==1); % only including the directions that have been completely probed 
    oog_idx = find(completed_dir(:,1)==1 & completed_dir(:,2)==1); % probed and out of gamut
    not_oog_idx = find(completed_dir(:,1)==1 & completed_dir(:,2)==0);
%     keyboard
    fact = 0.5/sqrt(tmp(probed_dirs,1).^2 + tmp(probed_dirs,2).^2); % factor needed to extract unit vector
    axes(a.axeshandles.TFR_dirs(ii)); compass(fact.*tmp(probed_dirs,1),fact.*tmp(probed_dirs,2),color(ii)); set(gca,'LineWidth',1);
    title(num_targetspikerates(ii));
    [THETA,RHO] = cart2pol(tmp(:,1),tmp(:,2));
    ind = (1:numel(THETA))';
    r = fliplr(linspace(0,1,numel(ind)));
    b = fliplr(r);
    THETA = THETA * (180/pi);
    % Earlier points in time are blue in color and later points in time are red in color 
    axes(a.axeshandles.TFR_isorespsurface(ii));
%     keyboard
    for jj = 1: numel(ind)
        if ~isempty(find(not_oog_idx==ind(jj)))
            plot(tmp(ind(jj),1), tmp(ind(jj),2),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[r(jj) 0 b(jj)],'PickableParts','none','MarkerEdgeColor',[r(jj) 0 b(jj)]); hold on;
        end
    end
    plot(0,0,'k*');
    set(gca,'Xlim',[lo,hi],'Ylim',[lo,hi],'ButtonDownFcn',@Display_basisvec_combo); grid on;
    axis square; xlabel('Basisvec 1'); ylabel('Basisvec 2');
%     plot(tmp(ind(ismember(ind,oog_idx)),1), tmp(ind(ismember(ind,oog_idx)),2),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor','g','PickableParts','none');
    if ~isempty(oog_idx)
        for jj = 1: numel(oog_idx)
            line([0 tmp(ind(ismember(ind,oog_idx(jj))),1)], [0 tmp(ind(ismember(ind,oog_idx(jj))),2)],'Color','black');
        end
    end
    line([oog1 oog1], [lo hi],'Color','green');
    line([lo hi], [oog2 oog2],'Color','green');
    hold off;
    end
% keyboard
set(gcf,'UserData',a);

end

function calcbasis_vec(~,~)
global gl neurothreshidx targetspikerateidx basisvecidx nstixperside
neurothreshmode = gl.stro.trial(:,neurothreshidx);
basisvec_dropidx = gl.inds(end); 
neurothresh_startidx = max(find(neurothreshmode==1,1),basisvec_dropidx+1);
num_targetspikerates = unique(gl.stro.trial(neurothresh_startidx:end,targetspikerateidx));
vect = gl.stro.ras{basisvec_dropidx,basisvecidx};
basisvec_size = nstixperside*nstixperside*3;
numvect = (numel(vect)/basisvec_size)-1;
basisvec = cell(1,numvect);
a = get(gcf,'UserData'); 
for ii = 1:numvect
    tmp_vec = vect((ii-1)*basisvec_size+1:basisvec_size*ii) + vect(numvect*basisvec_size+1:basisvec_size*(numvect+1));
    basisvec{ii} = reshape(tmp_vec,[nstixperside nstixperside 3]);
end
bkgnd_monitor = reshape(vect(numvect*basisvec_size+1:basisvec_size*(numvect+1)),[nstixperside nstixperside 3]);
grayimage = 0.5*ones(size(basisvec{1}));
basisvec1tobedisplayed = (basisvec{1}-bkgnd_monitor);
basisvec2tobedisplayed = (basisvec{2}-bkgnd_monitor);
basisvec1tobedisplayed = 0.5*basisvec1tobedisplayed/(0.01+ max(abs(basisvec1tobedisplayed(:))));
basisvec2tobedisplayed = 0.5*basisvec2tobedisplayed/(0.01+ max(abs(basisvec2tobedisplayed(:))));
basisvec1tobedisplayed = basisvec1tobedisplayed + grayimage;
basisvec2tobedisplayed = basisvec2tobedisplayed + grayimage;
image(basisvec1tobedisplayed,'Parent',a.axeshandles.basisvec1); set(a.axeshandles.basisvec1,'XTick',[],'YTick',[],'Box','on');
image(basisvec2tobedisplayed,'Parent',a.axeshandles.basisvec2); set(a.axeshandles.basisvec2,'XTick',[],'YTick',[],'Box','on');
set(gcf,'UserData',a);
gl.basisvec = basisvec;
gl.bkgnd_monitor = bkgnd_monitor;
end

function logrvstheta(~,~)
global gl neurothreshidx targetspikerateidx correctidx basisvecdiridx weightsidx stimonidx stimoffidx latencyidx parentverticesidx reversalflagidx
stro = gl.stro;
inds = gl.inds;
basisvec = gl.basisvec;
bkgnd_monitor = gl.bkgnd_monitor;
a = get(gcf,'UserData');
neurothreshmode = stro.trial(:,neurothreshidx);
basisvec_dropidx = inds(end); 
neurothresh_startidx = max(find(neurothreshmode==1,1),basisvec_dropidx+1);
num_targetspikerates = unique(stro.trial(neurothresh_startidx:end,targetspikerateidx));
norms = cell(1,numel(num_targetspikerates));
completed_search_alongdir = cell(1,numel(num_targetspikerates));
t_offset = stro.trial(end,latencyidx)/1000;
for jj = 1: numel(num_targetspikerates)
    idxs = find(~isnan(stro.trial(:,correctidx)) & stro.trial(:,targetspikerateidx)==num_targetspikerates(jj));
    idxs(idxs<=neurothresh_startidx) = [];
    different_weights = unique(stro.trial(idxs,basisvecdiridx));
    tmp_norm = [];
    tmp_completed_search_alongdir = [];

    for kk = 1:numel(different_weights)
        idxs1 = find(stro.trial(:,basisvecdiridx) == different_weights(kk));
        idxs1(idxs1<neurothresh_startidx) = [];
        raster_data = stro.ras(idxs1,1);
        tmp_norm = [tmp_norm; stro.ras{idxs1(end),weightsidx}'];
        for ii = 1:size(raster_data,1)
            tmp = raster_data{ii} ;
            spikes = tmp(tmp>stro.trial(idxs1(ii),stimonidx)+t_offset & tmp < stro.trial(idxs1(ii),stimoffidx));
            spikes = spikes - stro.trial(idxs1(ii),stimonidx)-t_offset;
        end
        [flag, gamutViolation] = Check_ifcompletely_probed(stro,basisvec,bkgnd_monitor,idxs1,tmp_norm(end,:));
        tmp_completed_search_alongdir = [tmp_completed_search_alongdir; flag gamutViolation];
    end
    norms{jj} = tmp_norm;
    completed_search_alongdir{jj} = tmp_completed_search_alongdir;
end

% Plotting the end norm or contrast values for each search direction and Converting the end norms into polar coordinates
% Need to write a small function to check for Gamut Violation based on the
% reversalflagidx
color = ['g','k'];
for ii = 1:size(norms,2)
    tmp = norms{ii};
    completed_dir = completed_search_alongdir{ii};
    probed_dirs = logical(completed_dir(:,1)==1); % only including the directions that have been completely probed 
    oog_idx = find(completed_dir(:,1)==1 & completed_dir(:,2)==1); % probed and out of gamut
    not_oog_idx = find(completed_dir(:,1)==1 & completed_dir(:,2)==0);
    fact = 0.5./sqrt(tmp(probed_dirs,1).^2 + tmp(probed_dirs,2).^2); % factor needed to extract unit vector
    [THETA,RHO] = cart2pol(tmp(:,1),tmp(:,2));
    ind = (1:numel(THETA))';
    r = fliplr(linspace(0,1,numel(ind)));
    b = fliplr(r);
    THETA = THETA * (180/pi);
    % Earlier points in time are blue in color and later points in time are red in color 
    axes(a.axeshandles.LogRvsTheta(ii));
    for jj = 1: numel(ind)
        if ~isempty(find(not_oog_idx==ind(jj)))
            plot(THETA(ind(jj)), log(RHO(ind(jj))),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[r(jj) 0 b(jj)],'MarkerEdgeColor',[r(jj) 0 b(jj)]); hold on;
        elseif ~isempty(find(oog_idx==ind(jj)))
            plot(THETA(ind(jj)), log(RHO((jj))),'o','MarkerSize',5,'LineWidth',1.5,'MarkerFaceColor','g','MarkerEdgeColor',[r(jj) 0 b(jj)]); hold on;
        end
    end
    grid on; axis square; xlabel('Theta'); ylabel('Log R');
    set(gca,'Xlim',[-180,180]); hold off;
end

% Plotting the firing rate map
for ss = 1:numel(num_targetspikerates)
    tmp_n = [];
    tmp_wts = [];
    num_dur =[];
    firing_rate = [];
    for ii = 1:size(norms{ss},1)
        dir = ii; 
        idxs1 = find(stro.trial(:,basisvecdiridx) == dir);
        for jj = 1:numel(idxs1)
            tmp_wts = [tmp_wts; stro.ras{idxs1(jj),weightsidx}'];
            tmp_n = [tmp_n; norm(stro.ras{idxs1(jj),weightsidx})];
        end
        raster_data = stro.ras(idxs1,1);
        for ii = 1:size(raster_data,1)
            tmp = raster_data{ii} ;
            % FR calculation
            spikes = tmp(tmp<stro.trial(idxs1(ii),stimoffidx) & tmp>stro.trial(idxs1(ii),stimonidx)+t_offset);
            num_dur = [num_dur; (stro.trial(idxs1(ii), stimoffidx)- stro.trial(idxs1(ii),stimonidx))-t_offset];
            firing_rate = [firing_rate; numel(spikes)/num_dur(end)];
        end
    end
    st = tpaps(tmp_wts',firing_rate');
    % An extended figure just for more clarity
    figure(2), set(gcf,'Name','FR map'); subplot(1,numel(num_targetspikerates),ss);
    plot3(tmp_wts(:,1),tmp_wts(:,2),firing_rate,'o','LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerSize',4); hold on;
    fnplt(st);hold off;
end
set(gcf,'UserData',a);
end

function BaselineFR(~,~)
global gl neurothreshidx stimonidx fpacqidx
stro = gl.stro;
inds = gl.inds;
basisvec = gl.basisvec;
a = get(gcf,'UserData');
neurothreshmode = stro.trial(:,neurothreshidx);

% Estimating the baseline firing rate for non-Neurothresh trials
% The number of spikes between the fpon and stimon codes
for i = 1:2
    idxs = find(neurothreshmode==(i-1));
    raster_data = stro.ras(idxs,1);
    num_spikes =[];
    num_dur = [];
    spikerate = [];
    add_time = 1;
    for ii = 1:size(raster_data,1)
        tmp = raster_data{ii} ;
        spikes_disp = tmp(tmp<stro.trial(idxs(ii),stimonidx)+ add_time & tmp>stro.trial(idxs(ii),fpacqidx));
        spikes = tmp(tmp<stro.trial(idxs(ii),stimonidx) & tmp>stro.trial(idxs(ii),fpacqidx));
        axes(a.axeshandles.BaselineFR(i));plot(spikes_disp-stro.trial(idxs(ii),fpacqidx),(ii-1)*ones(1,length(spikes_disp)),'k.'); hold on;
        num_spikes = [num_spikes; numel(spikes)];
        num_dur = [num_dur; (stro.trial(idxs(ii),stimonidx)- stro.trial(idxs(ii),fpacqidx))];
        firing_rate = numel(spikes)/(stro.trial(idxs(ii),stimonidx)- stro.trial(idxs(ii),fpacqidx)); % add time not added as I need am only counting spikes between fpacqidx and stimonidx to obtain the firing rate
        spikerate = [spikerate; firing_rate];
    end
    spikerate = num_spikes./num_dur; set(gca,'Xlim',[0-0.05 max(num_dur)+ add_time])
    line([0 0],[0 ii-1],'Color',[1 0 0])
    line([min(num_dur) min(num_dur)],[0 ii-1],'Color',[1 0 0])
    hold off;
    % plotting the baseline FR for each trial
    axes(a.axeshandles.Baselinespikerate(3-i)),plot(1:1:size(raster_data,1),spikerate,'Linewidth',2); xlabel ('Trials'), ylabel('Spikerate'); hold off;
end

end

function StoreRC(~,evt)
a = get(gcf,'UserData');
global row col comments
row = evt.Indices(1);
col = evt.Indices(2);
cla(a.axeshandles.file_comments,'reset');
axes(a.axeshandles.file_comments), text(0,0,comments(row));
set(gcf,'UserData',a);

end

function execute_file(~,~)
global row col maxT
a = get(gcf,'UserData');
stro = nex2stro(findfile(char(a.uitable.Data(row,col))));
Parse_stro(stro);
calcWN_check();
calcbasis_vec();
for ii = 1:maxT
    cla(a.axeshandles.STA_subunit(ii),'reset');set(a.axeshandles.STA_subunit(ii),'XTick',[],'YTick',[],'Box','on');
    cla(a.axeshandles.PC_subunit(ii),'reset'); set(a.axeshandles.PC_subunit(ii),'XTick',[],'YTick',[],'Box','on');
end
for jj = 1:2
    cla(a.axeshandles.TFR_dirs(jj),'reset'); set(a.axeshandles.TFR_dirs(jj),'XTick',[],'YTick',[],'Box','on');    
    cla(a.axeshandles.TFR_isorespsurface(jj),'reset'); set(a.axeshandles.TFR_isorespsurface(jj),'XTick',[],'YTick',[],'Box','on');
    cla(a.axeshandles.STA_subunits_temp(jj),'reset'); set(a.axeshandles.STA_subunits_temp(jj),'XTick',[],'YTick',[],'Box','on');
    cla(a.axeshandles.PC_subunits_temp(jj),'reset'); set(a.axeshandles.PC_subunits_temp(jj),'XTick',[],'YTick',[],'Box','on');
    cla(a.axeshandles.LogRvsTheta(jj),'reset'); set(a.axeshandles.LogRvsTheta(jj),'XTick',[],'YTick',[],'Box','on');
    cla(a.axeshandles.BaselineFR(jj),'reset'); set(a.axeshandles.BaselineFR(jj),'XTick',[],'YTick',[],'Box','on');
    cla(a.axeshandles.Baselinespikerate(jj),'reset'); set(a.axeshandles.Baselinespikerate(jj),'XTick',[],'YTick',[],'Box','on');
end
cla(a.axeshandles.filename,'reset');
cla(a.axeshandles.display_basisvec_combo,'reset'); set(a.axeshandles.display_basisvec_combo,'XTick',[],'YTick',[],'Box','on');
cla(a.axeshandles.display_basisvec_subunitRGB,'reset'); set(a.axeshandles.display_basisvec_combo,'XTick',[],'YTick',[],'Box','on');
axes(a.axeshandles.filename);text(0,0,[char(a.uitable.Data(row,col))],'FontSize',10);

set(gcf,'UserData',a);
end

function Display_basisvec_combo(~,evt)
% keyboard;
global gl nstixperside M
% pts = get(gca,'CurrentPoint');
x = evt.IntersectionPoint(1); y = evt.IntersectionPoint(2);
a = get(gcf,'UserData');
basisvec = gl.basisvec;
bkgnd_monitor = gl.bkgnd_monitor;
mask = gl.mask;
vec = (basisvec{1}-bkgnd_monitor)*x + (basisvec{2}-bkgnd_monitor)*y + bkgnd_monitor;
whichbutton = find(strcmp(get(gcf,'SelectionType'),{'normal','alt'}));
if any(vec(:)>1)
    disp('Image out of Gamut');
else
    axes(a.axeshandles.display_basisvec_combo); image(vec); set(gca,'XTick',[],'YTick',[]);
    mask = reshape(mask,[nstixperside nstixperside 3]);
    [s1x,s1y] = find(mask(:,:,1)==1,1);
    [s2x,s2y] = find(mask(:,:,1)==2,1);
    S1RGB = squeeze(vec(s1x,s1y,:) - bkgnd_monitor(s1x,s1y,:));
    S2RGB = squeeze(vec(s2x,s2y,:) - bkgnd_monitor(s2x,s2y,:));
    if whichbutton == 1
        axes(a.axeshandles.display_basisvec_subunitRGB); bar([S1RGB S2RGB]); % blue 1st vector, yellow 2nd vector
        set(gca,'XTick',[1 2 3],'XTickLabel',{'R','G','B'});
    elseif whichbutton == 2
        % convert the RGB to LMS
        S1LMS = M*S1RGB;
        S2LMS = M*S2RGB;
        axes(a.axeshandles.display_basisvec_subunitRGB); bar([S1LMS S2LMS]); 
        set(gca,'XTick',[1 2 3],'XTickLabel',{'L','M','S'});  
    end
end
set(gcf,'UserData',a);
end

function Loadfilenames(~,evt)
global NTmode comments filename
a = get(gcf,'UserData');
NTmodestring = get(a.uicontrols.NTmode,'String');
NTmodeval = get(a.uicontrols.NTmode,'Value');
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
comments = fetch(conn,'SELECT comments FROM WNthresh');
close(conn);
load filename_c.mat
load filename_l.mat
load S1LMS.mat
load S2LMS.mat
load newOCidx.mat
load newLMidx.mat
load newLUMidx.mat
load newSOidx.mat
load newhardtoclassifyidx.mat
isorespfiles = [filename_c; filename_l];
% [OCidx, LMidx, LUMidx, SOidx, hardtoclassifyidx] = classifycells(S1LMS,S2LMS);
OCidx = newOCidx;
LMidx = newLMidx;
LUMidx = newLUMidx;
SOidx = newSOidx;
hardtoclassifyidx = newhardtoclassifyidx;
if NTmodeval < 3
    new_idxs = find(strcmp(NTmode,NTmodestring(NTmodeval)));
    comments = comments(new_idxs);
    set(a.uitable,'Data',cellstr(filename(new_idxs)));
elseif NTmodeval == 3
    set(a.uitable,'Data',cellstr(filename));
elseif NTmodeval == 4
    set(a.uitable,'Data',cellstr(isorespfiles([OCidx LMidx])));
elseif NTmodeval == 5
    set(a.uitable,'Data',cellstr(isorespfiles([LUMidx])));
elseif NTmodeval == 6
    set(a.uitable,'Data',cellstr(isorespfiles([SOidx])));
elseif NTmodeval == 7
    set(a.uitable,'Data',cellstr(isorespfiles([hardtoclassifyidx])));
elseif NTmodeval == 8
    set(a.uitable,'Data',cellstr(isorespfiles([[OCidx LMidx LUMidx SOidx hardtoclassifyidx]])));
end
set(gcf,'UserData',a);
end
