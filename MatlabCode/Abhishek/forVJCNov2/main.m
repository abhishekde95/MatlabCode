% Figures for VJC on Nov 2nd
% author - Abhishek De, 10/17

%% starting with figure 1 
close all; clearvars;
stro = nex2stro(findfile('M020317005.nex')); % example luminance cell
global spikename maskidx spikeidx neurothreshidx nstixperside ngammasteps seedidx nframesidx correctidx
global fponidx stimoffidx stimonidx muidxs sigmaidxs basisvecidx weightsidx fpacqidx targetspikerateidx basisvecdiridx latencyidx reversalflagidx parentverticesidx
global msperframe ntrials maxT xx yy M linepredtol stepsizescale stepsize nreversals oogscale
spikename = getSpikenum(stro);
maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
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
maxT = 15; % this represents the temporal part in the spatiotemporal receptive field
xx = linspace(stro.sum.exptParams.gauss_locut/1000, stro.sum.exptParams.gauss_hicut/1000,ngammasteps); % xx represents the probabilities. For more info, have a look at the MATLAB 'norminv' function.
yy = norminv(xx'); % defining norminv to extract the values for which the cdf values range between gauss_locut and gauss_hicut
t_offset = stro.trial(end,latencyidx)/1000;

% Obtaining the M matrix, code extracted from Greg, fitting a cubic spline
% using the command 'spline'. 'SplineRaw' only availabe through
% psychtoolbox which I currently don't have now.
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

if isfield(stro.sum.exptParams,'nrepframes')
    if ~isnan(stro.sum.exptParams.nrepframes)
        nvblsperstimupdate = stro.sum.exptParams.nrepframes;
    else
        nvblsperstimupdate = 1;
    end
else
    nvblsperstimupdate = 1;
end
for mask_span = mask_changes(:,1)
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
%         frametimes = linspace(0, nframes*msperframe, nframes)+(msperframe/2)';
        frametimes = linspace(0, nframes*msperframe*nvblsperstimupdate, nframes)+(msperframe*nvblsperstimupdate/2)';
        spiketimes(spiketimes < maxT*msperframe) = [];
        spiketimes(spiketimes > frametimes(end)) = [];
        n = hist(spiketimes, frametimes);
        STCOVmex(rgbs(:),n);
    end
    out = STCOVmex('return');
    STS = out{1};
    nspikes = out{3};
    clear STCOVmex;
    clear out;
    STAs = STS/nspikes;
    an_mask = [];
    an_mask = zeros(nstixperside, nstixperside);   % Pixel mask.  Someday make a tool to make this non-zero.
    Lmask = logical(repmat(~an_mask(:),[3 1]));
    muvect = reshape(repmat([.5 .5 .5],nstixperside^2,1),nstixperside^2*3,1); % creating a 300 x 1 array with each entry as 0.5
    
    % The plotting begins here
    normfactor = 0.5/((max(abs(STAs(:))))*1.05); % This step is necessary to constrain the values within [-0.5, 0.5]
    STA = normfactor*(STAs(:,round(t_offset*1000/msperframe)))+muvect; % This makes the values fall back within a range of 0 and 1.
    STA = reshape(STA,[nstixperside nstixperside 3]); % Decomposing the STA into R, G and B plane
        
end

% Determining when Neurothresh mode was active, plotting the basis vector, working correctly
neurothreshmode = stro.trial(:,neurothreshidx);
basisvec_dropidx = inds(end); 
neurothresh_startidx = max(find(neurothreshmode==1,1),basisvec_dropidx+1); % There is possibility of pressing the Neurothresh button multiple times
num_targetspikerates = unique(stro.trial(neurothresh_startidx:end,targetspikerateidx));
vect = stro.ras{basisvec_dropidx,basisvecidx};
basisvec_size = nstixperside*nstixperside*3;
numvect = (numel(vect)/basisvec_size)-1;
basisvec = cell(1,numvect);
% Actual basis vec
for ii = 1:numvect
    tmp_vec = vect((ii-1)*basisvec_size+1:basisvec_size*ii) + vect(numvect*basisvec_size+1:basisvec_size*(numvect+1));
    basisvec{ii} = reshape(tmp_vec,[nstixperside nstixperside 3]);
end
bkgnd_monitor = reshape(vect(numvect*basisvec_size+1:basisvec_size*(numvect+1)),[nstixperside nstixperside 3]);

org_mask = stro.ras{mask_changes(1,2),maskidx};
org_mask(org_mask==Inf) = 0;
org_mask = reshape(org_mask,[10 10]);
ind1 = org_mask == 1;
ind2 = org_mask == 2;
outline1 = bwperim(ind1);
[a1,b1] = find(outline1 == 1);
a1 = a1 + 0.5; b1 = b1 + 0.5;
outline2 = bwperim(ind2);
[a2,b2] = find(outline2 == 1);
a2 = a2 + 0.5; b2 = b2 + 0.5;

subunitSTA = basisvec{1} + basisvec{2}-2*bkgnd_monitor;
subunitSTA = (0.5*subunitSTA./(max(subunitSTA(:))+0.01)) + 0.5;
figure(2),subplot(244), image(subunitSTA); set(gca,'XTick',[],'YTick',[]); axis square;
figure(2),subplot(248), image(STA); set(gca,'XTick',[],'YTick',[]); axis square;
dirs = [1 0; sqrt(1.5) 0.5; sqrt(0.5) sqrt(0.5); 0.5 sqrt(1.5); 0 1];
for ii = 1: size(dirs,1)
    newvec = dirs(ii,1)*(basisvec{1}-bkgnd_monitor) + dirs(ii,2)*(basisvec{2}-bkgnd_monitor);
    newvec = (0.5*newvec./(max(abs(newvec(:)))+0.01)) + 0.5;
    figure(2),subplot(1,5,ii),image(newvec); set(gca,'XTick',[],'YTick',[]); axis square;
end

for jj = 1:4
    newvec =(basisvec{1}-bkgnd_monitor);
    newvec = ((0.5*newvec./(max(abs(newvec(:)))+0.01))/(jj)) + 0.5;
    figure(3),subplot(2,4,jj),image(newvec); set(gca,'XTick',[],'YTick',[]); axis square;
    newvec = (basisvec{2}-bkgnd_monitor);
    newvec = ((0.5*newvec./(max(abs(newvec(:)))+0.01))/(jj)) + 0.5;
    figure(3),subplot(2,4,4+jj),image(newvec); set(gca,'XTick',[],'YTick',[]); axis square;
end
figure(2),set(gcf,'PaperPositionMode','auto'); 
figure(3),set(gcf,'PaperPositionMode','auto');
%%
weight_direction_mat = [];
for i=neurothresh_startidx:size(stro.trial,1)
    weight_direction_mat = [weight_direction_mat; stro.trial(i,basisvecdiridx) stro.ras{i,weightsidx}'/norm(stro.ras{i,weightsidx}')];
end
[~,idx] = sort(weight_direction_mat(:,1)); % sorting it according to the increasing order of weight direction indexes
weight_direction_mat1 = weight_direction_mat(idx,:);

% plotting the rasterplots for neurothresh trials
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
        % flag = 0, incompletely probed 
        % flag = 1, completely probed 
        % gamutViolation = 1, out of gamut point
        tmp_completed_search_alongdir = [tmp_completed_search_alongdir; flag gamutViolation];
    end
    norms{jj} = tmp_norm;
    completed_search_alongdir{jj} = tmp_completed_search_alongdir;
end

% Refer NTpreprocess written by Greg, 
% Plotting the end norm or contrast values for each search direction and Converting the end norms into polar coordinates
% Need to write a small function to check for Gamut Violation based on the
% reversalflagidx


color = ['g-';'k-';'r-';'b-';'m-'];
lo = -1.0; hi = 1.0;
for ii = 1:size(norms,2)
    tmp = norms{ii};
    completed_dir = completed_search_alongdir{ii};
    probed_dirs = logical(completed_dir(:,1)==1); % only including the directions that have been completely probed 
    oog_idx = find(completed_dir(:,1)==1 & completed_dir(:,2)==1); % probed and out of gamut
    not_oog_idx = find(completed_dir(:,1)==1 & completed_dir(:,2)==0);
    fact = 0.5/sqrt(tmp(probed_dirs,1).^2 + tmp(probed_dirs,2).^2); % factor needed to extract unit vector
    [THETA,RHO] = cart2pol(tmp(:,1),tmp(:,2));
    ind = (1:numel(THETA))';
    r = fliplr(linspace(0,1,numel(ind)));
    b = fliplr(r);
    THETA = THETA * (180/pi);
    % Earlier points in time are blue in color and later points in time are red in color 
    for jj = 1: numel(ind)
        if ~isempty(find(not_oog_idx==ind(jj)))
            figure(4+ii),plot(tmp(ind(jj),1), tmp(ind(jj),2),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[r(jj) 0 b(jj)],'PickableParts','none','MarkerEdgeColor',[r(jj) 0 b(jj)]); hold on;
        end
    end
    plot(0,0,'k*');
    set(gca,'Xlim',[lo,hi],'Ylim',[lo,hi]);grid on;
    axis square; xlabel('Basisvec 1'); ylabel('Basisvec 2');
    if ~isempty(oog_idx)
        for jj = 1: numel(oog_idx)
            line([0 tmp(ind(ismember(ind,oog_idx(jj))),1)], [0 tmp(ind(ismember(ind,oog_idx(jj))),2)],'Color','black');
        end
    end
    hold off;
end
figure(5), axis square;
figure(6), axis square;


