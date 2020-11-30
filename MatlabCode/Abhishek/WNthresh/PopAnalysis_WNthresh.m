% PopAnalysis for WNthresh Data
% Author - Abhishek De, 06/17
% Obsolete code, Abhishek - 3/18
close all; clearvars;
plot_counter = 1;
mode = 3; % 1 - all cells in the database, 2 - all color cells in Neurothresh, 3 - all luminance cells in Neurothresh
if mode == 1
    conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
    filename1 = fetch(conn,'SELECT filename FROM WNthresh');
    filename2 = fetch(conn,'SELECT filename FROM WNSubunit');
    NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
    STAmode = fetch(conn,'SELECT mode FROM WNSubunit');
    close(conn);
    subunit = find(strcmp(NTmode,'subunit'));
    tmp_STAmode = find(strcmp(STAmode, 'STA'));filename = ['M122416004.nex';'M122116001.nex';'M122016003.nex';'M121616004.nex';'M121216001.nex'; ...
                'M120916003.nex';'M123116004.nex';'M010217001.nex';'P072816001.nex';'P080216003.nex'; ...
                'P080316001.nex';'M010317003.nex';'M013017001.nex';'M020517001.nex';'M021317001.nex'; ...
                'M050617003.nex';'M061517001.nex';'M080817002.nex';'M081417003.nex';'M090117002.nex'; ...
                'P090117002.nex';'P090217001.nex';'M090617002.nex';'M090717002.nex';'M090817001.nex'; ...
                'M091117002.nex'];
          
    filename = cellstr(filename);
elseif mode == 3
    filename = ['M122616001.nex';'M121916001.nex';'M120216001.nex';'M122916001.nex';'M011217002.nex'; ...
                'M011417003.nex';'M013117001.nex';'M020317005.nex';'M021617003.nex';'M022717003.nex'; ...
                'M030117001.nex';'M050517001.nex';'M050617004.nex';'M050717001.nex';'M061217001.nex'; ...
                'M061317005.nex';'M061417001.nex';'M070117003.nex';'M073117002.nex';'M080217002.nex'; ...
                'M082517004.nex';'P082917001.nex'];
    filename = cellstr(filename);
    filename = [cellstr(filename1(subunit,:)); cellstr(filename2(tmp_STAmode,:))];
elseif mode == 2
    
end
    
% Need to convert the RGB values to LMS values and calculate following
% things :-
% 1) Check if there is any time difference between the latency of individual subunits; basically see if you can group/classify cells this way 
% 2) Plot the RF locations 
% 3) Calculate the correlation of color between the adjacent subunits
t_offset = [];
time_diff = [];
time_diff_zero_idx = [];
time_diff_nonzero_idx = [];
RF_loc = [];
subunit_chrom_corr = [];
Pangu = [];
Maui = [];
global nstixperside 
RGBS1 = cell(1,numel(filename));
RGBS2 = cell(1,numel(filename));
for ii = 1:numel(filename)
    ii
    fileofinterest = char(filename(ii,:));
    if strcmp(fileofinterest(1),'M') == 1
        Maui = [Maui; ii];
    else
        Pangu = [Pangu; ii];
    end
    stro = nex2stro(findfile(fileofinterest));
    spikename = 'sig001a';%getSpikenum(stro);
    maskidx = strcmp(stro.sum.rasterCells(1,:), 'subunit_mask');
    spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
    nstixperside = stro.sum.exptParams.nstixperside;
    ngammasteps = 2^16; % 65536
    seedidx = strcmp(stro.sum.trialFields(1,:),'seed');
    nframesidx = strcmp(stro.sum.trialFields(1,:),'num_frames');
    stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
    basisvecidx = strcmp(stro.sum.rasterCells(1,:),'basis_vec');
    %     latencyidx = strcmp(stro.sum.trialFields(1,:),'latency');
    muidxs = [find(strcmp(stro.sum.trialFields(1,:),'mu1')), ...
        find(strcmp(stro.sum.trialFields(1,:),'mu2')), ...
        find(strcmp(stro.sum.trialFields(1,:),'mu3'))];
    sigmaidxs = [find(strcmp(stro.sum.trialFields(1,:),'sigma1')), ...
        find(strcmp(stro.sum.trialFields(1,:),'sigma2')), ...
        find(strcmp(stro.sum.trialFields(1,:),'sigma3'))];
    msperframe = 1000/stro.sum.exptParams.framerate;
    ntrials = size(stro.trial,1);
    maxT = 15; % this represents the temporal part in the spatiotemporal receptive field
    xx = linspace(stro.sum.exptParams.gauss_locut/1000, stro.sum.exptParams.gauss_hicut/1000,ngammasteps); % xx represents the probabilities. For more info, have a look at the MATLAB 'norminv' function.
    yy = norminv(xx'); % defining norminv to extract the values for which the cdf values range between gauss_locut and gauss_hicut
    
    mask_changes = [2];
    all_masks = stro.ras(:,maskidx);
    if numel(unique(basisvecidx))>1
        Fx = @(xi) any(isnan(xi)); % function that finds 'NaN' in a cell array
        inds = find(cellfun(Fx,stro.ras(:,basisvecidx))==0);
        if isempty(inds)
            inds = size(stro.trial,1)-1;
        end
        last_wntrial =  inds(1)-1;
    else
        last_wntrial = ntrials;
    end
    for k = 3:last_wntrial
        if isequal(all_masks{k}, all_masks{k-1}) %|| all(all_masks{k} == 0) && any(isnan(all_masks{k-1}))
            continue
        else
            mask_changes = [mask_changes k-1 k];
        end
    end
    if mask_changes(end) == last_wntrial
        mask_changes(end) = [];
    else
        mask_changes = [mask_changes  last_wntrial];
    end
    mask_changes = reshape(mask_changes , 2, []);
    trial_span = mask_changes(:,2);
    % keyboard
    st_mask = stro.ras{trial_span(1),maskidx}; % subunit mask
    st_mask(st_mask == 0) = Inf;
    [stIdxs,~,~] = unique(st_mask); % now the Infs map to nsubunits+1
    num_subunits = length(stIdxs)-any(isinf(stIdxs)); % nsubunits, like subunits A and B
    if ~num_subunits
        num_subunits = nstixperside^2;
    end
    STCOV_st('init', {num_subunits 3 maxT});
    for k = trial_span(1):trial_span(2)
        nframes = stro.trial(k,nframesidx);
        if (nframes == 0)
            continue;
        end
        seed = stro.trial(k,seedidx);
        mu = stro.trial(k,muidxs)/1000;
        sigma = stro.trial(k,sigmaidxs)/1000;
        org_mask = stro.ras{k,maskidx};
        if any(org_mask)
            org_mask(org_mask == 0) = Inf;
            [subunitIdxs,~,mask] = unique(org_mask); % now the Infs map to nsubunits+1
            nrandnums_perchannel = length(subunitIdxs)-any(isinf(subunitIdxs)); % nsubunits, like subunits A and B
            mask = [mask; mask+max(mask); mask+2*max(mask)]; %#ok<AGROW>
        else
            nrandnums_perchannel = nstixperside^2; % In this case, it is the 100 pixels which are flickering when no subunits are selected
        end
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
        spiketimes(spiketimes < maxT*msperframe) = [];
        spiketimes(spiketimes > frametimes(end)) = [];
        n = hist(spiketimes, frametimes);
        STCOV_st(rgbs(:),n);
    end
    
    out = STCOV_st('return'); % returns the covariance matrix on frame by frame basis
    STS = out{1};  % A (dimension) x 9(frames) matrix
    STCross = out{2};  % A (dimension x frames)x (dimension x frames)  matrix
    nspikes = out{3}; % Number of spikes in the given file
    clear STCOV_st out
    % Coverting the STS and the STCross into STA and STC respectively
    STAs = STS/nspikes;
    tmp = STS(:)*STS(:)';
    
    % Calculating the latency of the peak energy
    latencyidx = strcmp(stro.sum.trialFields(1,:),'latency');
        
    % Flipping the STAs such that the last frame appears first and the first frame appears last
    STAs = fliplr(STAs);
    basisvec1 = STAs(1:2:5,:); basisvec2 = STAs(2:2:6,:);
    basisvec1_energy = sum(basisvec1.^2); basisvec2_energy = sum(basisvec2.^2);
    [~,ind1] = max(basisvec1_energy);
    [~,ind2] = max(basisvec2_energy);
    time_diff = [time_diff; (ind1-ind2)*msperframe];
%     if isnan(stro.trial(end,latencyidx))
    [~,ind3] = max(sum((basisvec1+basisvec2).^2));
    t_offset = [t_offset; ind3*msperframe];
%     else
%         t_offset = [t_offset; stro.trial(end,latencyidx)];
%     end
    
    frame_to_select = round(t_offset(end)/msperframe);
    vec = expand_vector(STAs(:,frame_to_select),nrandnums_perchannel,mask,1);
    vec = reshape(vec,[nstixperside nstixperside 3]);
    im = (0.5*vec/(0.0001 + max(abs(vec(:))))) + 0.5;
    
    
    if abs(time_diff(end))>0
        time_diff_nonzero_idx = [time_diff_nonzero_idx; ii];
        figure(plot_counter),subplot(8,10,numel(time_diff_nonzero_idx)), image(im); set(gca,'XTick',[],'YTick',[]);drawnow;
    else
        time_diff_zero_idx = [time_diff_zero_idx; ii];
        figure(plot_counter+1),subplot(8,10,numel(time_diff_zero_idx)), image(im); set(gca,'XTick',[],'YTick',[]);drawnow;
    end
    
    % Acquiring the receptive field location
    RF_loc = [RF_loc; stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y]; 
    
    % Calculating the correlation of the gun weights between the adjacent
    % subunits
    tmp_RGBS1 = STAs(1:2:5,frame_to_select);
    tmp_RGBS2 = STAs(2:2:6,frame_to_select);
    subunit_chrom_corr = [subunit_chrom_corr; corr(tmp_RGBS1,tmp_RGBS2)];
    RGBS1{ii} = tmp_RGBS1;
    RGBS2{ii} = tmp_RGBS2;
end

plot_counter = plot_counter + 2;
% Histogram of time to reach the peak
figure(plot_counter), hist(t_offset); xlabel('ms'), ylabel('Frequency'), title('Time to peak');  
plot_counter = plot_counter + 1;

% Histogram of the latency differences between the V1 subunits
figure(plot_counter), hist(time_diff); xlabel('Time difference'), ylabel('Frequency'); title('Latency Differences between subunits');
plot_counter = plot_counter + 1;

% Plot of the RF location
figure(plot_counter), plot(RF_loc(Pangu,1), RF_loc(Pangu,2),'o','LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerSize',5); hold on;
plot(RF_loc(Maui,1), RF_loc(Maui,2),'o','LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerSize',5); 
xlabel('X'), ylabel('Y'); title('RF locations'); plot_counter = plot_counter + 1; hold off;

% Plot of the gun weights correlation between the adjacent subunits
figure(plot_counter), hist(subunit_chrom_corr); xlabel('Correlation values'); ylabel('Frequency'); title('Gun Weight correlation between subunits');
plot_counter = plot_counter + 1;

%% Trying something new with the CIE chromaticity plots
load T_stiles10.mat
bkgnd_r = strcmp(stro.sum.trialFields(1,:),'bkgnd_r');
bkgnd_g = strcmp(stro.sum.trialFields(1,:),'bkgnd_g');
bkgnd_b = strcmp(stro.sum.trialFields(1,:),'bkgnd_b');
bkgnd_col = [stro.trial(1,bkgnd_r); stro.trial(1,bkgnd_g); stro.trial(1,bkgnd_b)]/255;
bkgnd_xyY = XYZToxyY(SRGBPrimaryToXYZ(bkgnd_col));

xyYs = [];
for ii = 1:numel(filename)
    xyYs = [xyYs XYZToxyY(SRGBPrimaryToXYZ(RGBS1{ii}+bkgnd_col)) XYZToxyY(SRGBPrimaryToXYZ(RGBS1{ii}+bkgnd_col))];
end 

stiles10xyY = XYZToxyY(SRGBPrimaryToXYZ(T_stiles10));
figure(plot_counter),scatter(stiles10xyY(1,:),stiles10xyY(2,:),'b','filled'); hold on;
scatter(bkgnd_xyY(1),bkgnd_xyY(2),'k','filled'); 
scatter(xyYs(1,:),xyYs(2,:),'g','filled'); xlabel('x'); ylabel('y'); title('Chromaticity plot'); hold off; plot_counter = plot_counter + 1;