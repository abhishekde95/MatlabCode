% A new scipt for checking the effective contrast onto the RF using WN
% checkerboard, WN subunit and isopresp stimuli.
% Author - Abhishek De, 5/23

% Section 1: Compute the STA of an example DO cell

% Section 2: Stimulus comparison between the WN checkerboard, WN subunit and the gamut edge 

% Section 3: PLotting the firing rate map of an example DO cell along with isoresponse contours from a GLM fit 

%% Section 1: Compute the STA of an example DO cell

close all; clearvars;
plot_counter = 1;

global nstixperside maxT nframesidx
% Loading all the files
try 
    % Using the JDBC connection
    conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
    filename = fetch(conn,'SELECT filename FROM WNthresh');
    NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
    spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
    close(conn);
    filename = filename(strcmp(string(NTmode),"subunit"));
    NTmode = NTmode(strcmp(string(NTmode),"subunit"));
    spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));

catch
    csv_filename = '/Users/abhishekde/Desktop/MatlabCode/Abhishek/CSV_PHPmyadmin_files/WNthresh.csv';
    [filename, NTmode, spikeIdx] = get_WNthreshdata_from_csvfile(csv_filename, 'subunit');
    spikeidx_NT = str2num(cell2mat(spikeIdx));
end


% Include Gun noise data for statistical tests in order to estimate the RF size
channels = 3;
NPOINTS = 65536;
resize_fact2 = 1;
CHI2CRIT = 0.95; % For flaging a stixel as significant (adding gun and cone noise z-scores, squared)
maxT = 22;
crit = chi2inv(CHI2CRIT,300); % 3 color channels

% 24- example BY DO cell
% 31 - example RG DO cell
WN = nex2stro(findfile(char(filename(24))));
framerate = WN.sum.exptParams.framerate;
nstixperside = WN.sum.exptParams.nstixperside;
ntrials = length(WN.sum.absTrialNum);
seedidx = strcmp(WN.sum.trialFields(1,:),'seed');
stimonidx = find(strcmp(WN.sum.trialFields(1,:),'stim_on'));
stimoffidx = find(strcmp(WN.sum.trialFields(1,:),'all_off'));
nframesidx = find(strcmp(WN.sum.trialFields(1,:),'num_frames'));
noisetypeidx = find(strcmp(WN.sum.trialFields(1,:),'noise_type'));
sigmaidxs = strmatch('sigma',WN.sum.trialFields(1,:));
hepidx = find(strcmp(WN.sum.rasterCells(1,:),'AD11'));
vepidx = find(strcmp(WN.sum.rasterCells(1,:),'AD12'));
maskidx = strcmp(WN.sum.rasterCells(1,:), 'subunit_mask');
anlgStartTimeidx = find(strcmp(WN.sum.rasterCells(1,:),'anlgStartTime'));
mu1idx = find(strcmp(WN.sum.trialFields(1,:),'mu1'));
mu2idx = find(strcmp(WN.sum.trialFields(1,:),'mu2'));
mu3idx = find(strcmp(WN.sum.trialFields(1,:),'mu3'));
sigma1idx = find(strcmp(WN.sum.trialFields(1,:),'sigma1'));
sigma2idx = find(strcmp(WN.sum.trialFields(1,:),'sigma2'));
sigma3idx = find(strcmp(WN.sum.trialFields(1,:),'sigma3'));
maskidx = strcmp(WN.sum.rasterCells(1,:),'subunit_mask');
basisvecidx = strcmp(WN.sum.rasterCells(1,:),'basis_vec');
latencyidx = strcmp(WN.sum.trialFields(1,:),'latency');
neurothreshidx = strcmp(WN.sum.trialFields(1,:),'neurothresh');
muidxs = [find(strcmp(WN.sum.trialFields(1,:),'mu1')), ...
    find(strcmp(WN.sum.trialFields(1,:),'mu2')), ...
    find(strcmp(WN.sum.trialFields(1,:),'mu3'))];
sigmaidxs = [find(strcmp(WN.sum.trialFields(1,:),'sigma1')), ...
    find(strcmp(WN.sum.trialFields(1,:),'sigma2')), ...
    find(strcmp(WN.sum.trialFields(1,:),'sigma3'))];

L = WN.trial(:,noisetypeidx)==1;
gammaTable = WN.sum.exptParams.gamma_table;
gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
invgamma = InvertGamma(gammaTable, 0);
ngammasteps = 2^16; % 65536; %get number of rows of gamma_table (65536)
sigmavect = unique(WN.trial(L,[sigma1idx sigma2idx sigma3idx]),'rows')/1000;
sigmavect(all(any(sigmavect == 0),2),:) = [];
gausslims = [WN.sum.exptParams.gauss_locut WN.sum.exptParams.gauss_hicut]/1000;
x = linspace(gausslims(1),gausslims(2),NPOINTS);
Fx = norminv(x)*sigmavect(1);
sigmacorrectionfactor = std(Fx)./sigmavect(1);
muvar = (sigmavect(1)*sigmacorrectionfactor)^2;
xx = linspace(WN.sum.exptParams.gauss_locut/1000, WN.sum.exptParams.gauss_hicut/1000,ngammasteps); % xx represents the probabilities. For more info, have a look at the MATLAB 'norminv' function.
yy = norminv(xx'); % defining norminv to extract the values for which the cdf values range between gauss_locut and gauss_hicut

% Getting the background rgb/lms
% Calculating the M matrix
fundamentals = WN.sum.exptParams.fundamentals;
mon_spd = WN.sum.exptParams.mon_spd;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

ridx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_r'));
gidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_g'));
bidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_b'));
bkgndRGB = [mode(WN.trial(:,ridx)), mode(WN.trial(:,gidx)), mode(WN.trial(:,bidx))];
bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
bkgndlms = M*bkgndrgb;
Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
Mrgbtocc = inv(Mrgbtocc');

mask_changes = [2 size(WN.trial,1)];
if any(basisvecidx)
    mask_changes = [2];
    all_masks = WN.ras(:,maskidx);
    Fx = @(xi) any(isnan(xi)); % function that finds 'NaN' in a cell array
    inds = find(cellfun(Fx,WN.ras(:,basisvecidx))==0);
    if isempty(inds)
        inds = size(WN.trial,1)-1;
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
    
    
    % Making a new WN subunit structure from the WN structure
    idxs_subunit = zeros(size(WN.trial,1),1);
    idxs_subunit(mask_changes(2,1):mask_changes(2,2)) = 1;
    idxs_subunit = logical(idxs_subunit);
    
    WN_subunit = WN;
    WN_subunit.ras(~idxs_subunit,:) = [];
    WN_subunit.trial(~idxs_subunit,:) = [];
    
    % Same as above but for WN checkerboard stimulus
    idxs_check = zeros(size(WN.trial,1),1);
    idxs_check(mask_changes(1,1):mask_changes(1,2)) = 1;
    idxs_check = logical(idxs_check);
    
    WN_check = WN;
    WN_check.ras(~idxs_check,:) = [];
    WN_check.trial(~idxs_check,:) = [];
end
spikename = 'sig001a';
spikeidx = strcmp(WN.sum.rasterCells(1,:),spikename);
msperframe = 1000/WN.sum.exptParams.framerate;

% Calculating STA from the subunit WN
out_subunit = getWhtnsStats_AD(WN_subunit,maxT,'STCOVmex',{2,3,maxT},2,spikename);
STS_gun = out_subunit{1}; nspikes_gun = out_subunit{3}; clear out_subunit;
STAs_subunit = STS_gun/nspikes_gun;

%% Section 2: Stimulus comparison between the WN checkerboard, WN subunit and the gamut edge
nrandnums_perchannel = 2;
initargs = cell(1,2);

% Creating initargs for the WN check noise 
mask = WN.ras{mask_changes(1,2),maskidx};
mask1 = mask;
mask1(mask==0) = nrandnums_perchannel+1;
mask2 = [mask1; mask1+max(mask1); mask1+2*max(mask1)];
reformed_vec_forcheckstim = expand_vector(STAs_subunit,nrandnums_perchannel,mask2,maxT);

subunit1_mask = repmat(logical(mask(1:100)==1),[3*maxT 1]);
subunit2_mask = repmat(logical(mask(1:100)==2),[3*maxT 1]);
basis = [reformed_vec_forcheckstim(:).*subunit1_mask reformed_vec_forcheckstim(:).*subunit2_mask];
initargs{1} = {basis , 0, sum(WN_check.trial(:,nframesidx)), [nstixperside.^2 3 maxT]};

% Creating initargs for the WN subunit noise 
initargs{2} = define_basis_vec(nrandnums_perchannel,WN_subunit,STAs_subunit,[],0);

for ii = 1:2
    STPROJmod('init',initargs{ii}); % initialising the STPROJmod
    
    for k = mask_changes(1,ii):mask_changes(2,ii)
        nframes = WN.trial(k,nframesidx);
        if (nframes == 0)
            continue;
        end
        seed = WN.trial(k,seedidx);
        mu = WN.trial(k,muidxs)/1000;
        sigma = WN.trial(k,sigmaidxs)/1000;
        
        % org_mask tells u if u have updated the mask or not. If org_mask is non-zero it means at this particular trial
        % u have selected the subunits and need to analyse its computation
        org_mask = WN.ras{k,maskidx};
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
        t_stimon = WN.trial(k, stimonidx);
        spiketimes = (WN.ras{k,spikeidx}-t_stimon)*1000; % observing spiketimes in milliseconds
        frametimes = linspace(0, nframes*msperframe, nframes)+(msperframe/2)';
        spiketimes(spiketimes < maxT*msperframe) = [];
        spiketimes(spiketimes > frametimes(end)) = [];
        n = hist(spiketimes, frametimes);
        STPROJmod(rgbs(:),n);
        
    end
    out = STPROJmod('return');
    projs = out{1};
    Lspike = out{2};
    clear STPROJmod out;
    
    % Rotate the data by X degrees; X=0 to 180 
    [X_mod,Y_mod] = rotate_data(projs(:,1),projs(:,2),0);
    projs = [X_mod Y_mod];
    
    % Calculating some stuff that will be useful for plotting 
    min_val = floor(min(projs(:))*100)/100;
    max_val = ceil(max(projs(:))*100)/100;
    num_bins = 15;
    bin_interval = (max_val-min_val)/num_bins;
    nbins1 = min_val:bin_interval:max_val;
    nbins = linspace(prctile(projs(:),5), prctile(projs(:),95),numel(nbins1)+2);
    xmin = min(nbins); xmax = max(nbins);
    new_nbins = [nbins(1)-mean(diff(nbins)), nbins, mean(diff(nbins)) + nbins(end)];
    [n_spike,out_spike] = hist3([projs(Lspike>0,1), projs(Lspike>0,2)],{new_nbins,new_nbins});
    [n_raw,out_raw] = hist3([projs(:,1), projs(:,2)],{new_nbins,new_nbins});
    n_spike = n_spike(2:end-1,2:end-1);
    n_raw = n_raw(2:end-1,2:end-1);
    n_raw(n_raw==0) = 1; % You can either enter 'NaN'or '1' to avoid division by zero
    non_lin = n_spike./n_raw;
    
    if ii == 1
        projs_check = projs;
        Lspike_check = Lspike;
        non_lin_check = non_lin;
        xmin_check = xmin;
        xmax_check = xmax;
        
    else 
        projs_subunit = projs;
        Lspike_subunit = Lspike;
        non_lin_subunit = non_lin;
        xmin_subunit = xmin;
        xmax_subunit = xmax;
    end

end

% Finding out the maximum projection using 
STA_subunit1 = STAs_subunit; STA_subunit2 = STAs_subunit;
STA_subunit1(2:2:end,:) = deal(0);
STA_subunit2(1:2:end,:) = deal(0);

% Creating bkgnrgb matrix 
k = repmat(bkgndrgb',[2 1]);
bkgnd = repmat(k(:),[1 maxT]);

% Creating radial direction with each norm =1 
N = 200;
theta = linspace(0,2*pi,N);
[x,y] = pol2cart(theta,ones(size(theta)));

fact = zeros(N,2);
for ii = 1:N
    image = x(ii)*STA_subunit1 + y(ii)*STA_subunit2 + bkgnd;
    tmp = min([min(abs((1-bkgnd(:))./(image(:)-bkgnd(:)))) min(abs((0-bkgnd(:))./(image(:)-bkgnd(:))))]);
    new_im = tmp*(image-bkgnd);
    proj_val = [dot(new_im(:),STA_subunit1(:)) dot(new_im(:),STA_subunit2(:))];
    proj_val = [proj_val(1)/norm(STA_subunit1(:)) proj_val(2)/norm(STA_subunit2(:))];
    fact(ii,:) = proj_val;
end

% Rotating the data by 45 degrees
[X_mod,Y_mod] = rotate_data(fact(:,1),fact(:,2),0);
fact = [X_mod Y_mod];

% Plotting 
figure(plot_counter); 
subplot(221);  hold on;
plot(projs_subunit(:,1),projs_subunit(:,2),'r.'); plot(projs_check(:,1),projs_check(:,2),'k.');
axis square; set(gca,'Tickdir','out','Xlim',[-1.3 1.3],'XTick',-1.2:0.6:1.2,'Ylim',[-1.3 1.3],'YTick',-1.2:0.6:1.2); xlabel('Projection magnitude S1'); ylabel('Projection magnitude S2'); title('Stimulus comparison');

subplot(224); imagesc([xmin_subunit xmax_subunit],[xmin_subunit xmax_subunit],non_lin_subunit); hold on; 
axis square; set(gca,'Tickdir','out','Xlim',[-1.3 1.3],'XTick',-1.2:0.6:1.2,'Ylim',[-1.3 1.3],'YTick',-1.2:0.6:1.2); xlabel('Projection magnitude S1'); ylabel('Projection magnitude S2'); title('WN subunit: Firing rate map'); axis xy; colormap('jet');

subplot(221); hold on; plot(fact(:,1),fact(:,2),'-g','Linewidth',2); plot([-1.3 1.3],[0 0],'k'); plot([0 0],[-1.3 1.3],'k');
legend('Subunit','Check','Gamut edge'); hold off;
subplot(224); hold on; plot(fact(:,1),fact(:,2),'-g','Linewidth',2); plot([-1.3 1.3],[0 0],'k'); plot([0 0],[-1.3 1.3],'k'); hold off;


% Visualizing the density of the WN checkerboard and subunit stimulus projection distributions 
bins = -1.3:0.02:1.3;

[WNcheckhist,~] = hist3([projs_check(:,1), projs_check(:,2)],{bins,bins});
subplot(222); imagesc([bins(1) bins(end)],[bins(1) bins(end)],WNcheckhist); hold on; 
hold on; plot(fact(:,1),fact(:,2),'-g','Linewidth',2); plot([-1.3 1.3],[0 0],'k'); plot([0 0],[-1.3 1.3],'k');
axis square; set(gca,'Tickdir','out','Xlim',[-1.3 1.3],'XTick',-1.2:0.6:1.2,'Ylim',[-1.3 1.3],'YTick',-1.2:0.6:1.2); xlabel('Projection magnitude S1'); ylabel('Projection magnitude S2'); title('WN check stim distribution'); axis xy; colormap('gray');


[WNsubunithist,~] = hist3([projs_subunit(:,1), projs_subunit(:,2)],{bins,bins});
subplot(223); imagesc([bins(1) bins(end)],[bins(1) bins(end)],WNsubunithist); hold on; 
hold on; plot(fact(:,1),fact(:,2),'-g','Linewidth',2); plot([-1.3 1.3],[0 0],'k'); plot([0 0],[-1.3 1.3],'k');
axis square; set(gca,'Tickdir','out','Xlim',[-1.3 1.3],'XTick',-1.2:0.6:1.2,'Ylim',[-1.3 1.3],'YTick',-1.2:0.6:1.2); xlabel('Projection magnitude S1'); ylabel('Projection magnitude S2'); title('WN subunit stim distribution'); axis xy; colormap('gray');

plot_counter = plot_counter + 1;


%%  Section 3: PLotting the firing rate map of an example DO cell along with isoresponse contours from a GLM fit 

% Fitting a GLM to the actual RGB triplet data
ndiv = 15;
modelbins = linspace(xmin_subunit,xmax_subunit,ndiv);
[modelbinsX, modelbinsY] = meshgrid(modelbins);
mdllin =  fitglm(projs_subunit,logical(Lspike_subunit),'linear','Distribution','binomial','Link','logit');
predlin = predict(mdllin,[modelbinsX(:) modelbinsY(:)]); % perdiction from GLM

figure(plot_counter);
subplot(221); imagesc([xmin_subunit xmax_subunit],[xmin_subunit xmax_subunit],non_lin_subunit); hold on; 
axis square; set(gca,'Tickdir','out','Xlim',[-0.4 0.4],'XTick',-0.4:0.2:0.4,'Ylim',[-0.4 0.4],'YTick',-0.4:0.2:0.4); xlabel('Projection magnitude S1-S2'); 
ylabel('Projection magnitude S1+S2'); title('WN subunit: Firing rate map'); axis xy; colormap('gray');

% NOTE: For some reason, I am not able to ablign the orientation of the firing
% rate map and the GLM based firing rate map fits and contours
GLMresult = reshape(predlin,[ndiv ndiv]);

subplot(222); imagesc([xmin_subunit xmax_subunit],[xmin_subunit xmax_subunit],GLMresult');hold on; 
axis square; set(gca,'Tickdir','out','Xlim',[-0.4 0.4],'XTick',-0.4:0.2:0.4,'Ylim',[-0.4 0.4],'YTick',-0.4:0.2:0.4); xlabel('Projection magnitude S1-S2'); 
ylabel('Projection magnitude S1+S2'); title('Firing rate map: GLM'); axis xy; colormap('gray');

subplot(223); contour(modelbinsX,modelbinsY,GLMresult',5,'Linewidth',2);hold on; 
axis square; set(gca,'Tickdir','out','Xlim',[-0.4 0.4],'XTick',-0.4:0.2:0.4,'Ylim',[-0.4 0.4],'YTick',-0.4:0.2:0.4); xlabel('Projection magnitude S1-S2'); 
ylabel('Projection magnitude S1+S2'); title('Contour map: GLM'); axis xy; colormap('gray');

plot_counter = plot_counter + 1;

