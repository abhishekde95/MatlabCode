% This script is mostly to see how well the chromatic signature of the STAs
% obtained from WhiteNoise check and subunit stimuli match
% Author  - Abhishek 06/17
close all; clearvars;
plot_counter = 1;
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename1 = fetch(conn,'SELECT filename FROM WNthresh');
filename2 = fetch(conn,'SELECT filename FROM WNSubunit');
NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
STAmode = fetch(conn,'SELECT mode FROM WNSubunit');
close(conn);
subunit = find(strcmp(NTmode,'subunit'));
tmp_STAmode = find(strcmp(STAmode, 'STA'));
filename = [cellstr(filename1(subunit,:)); cellstr(filename2(tmp_STAmode,:))];
STA_subunit_corr = [];
STAs = cell(numel(filename),2);

for mm = 1:numel(filename)
    disp(mm);
    stro = nex2stro(findfile(char(filename(mm,:))));
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
    
    for jj = 1:2
        trial_span = mask_changes(:,jj);
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
        nspikes = out{3}; % Number of spikes in the given file
        clear STCOV_st out
        % Coverting the STS and the STCross into STA and STC respectively
        tmp_STAs = STS/nspikes;
        tmp_STAs = fliplr(tmp_STAs);
        % Flipping the STAs such that the last frame appears first and the first frame appears last
        
        STAs{mm,jj} = tmp_STAs;
        
    end
    
    ind1 = find(org_mask==1);
    ind2 = find(org_mask==2);
    tmp_STAs = STAs{mm,1};
    tmp_STAs = [sum(tmp_STAs(ind1,:)); sum(tmp_STAs(ind2,:)); sum(tmp_STAs(ind1+100,:));sum(tmp_STAs(ind2+100,:)); sum(tmp_STAs(ind1+200,:)); sum(tmp_STAs(ind2+200,:))];
    STAs{mm,1} = tmp_STAs;
    vec1 = STAs{mm,1};
    vec2 = STAs{mm,2};
    STA_subunit_corr = [STA_subunit_corr; corr(vec1(:),vec2(:))];
end

figure(1),hist(STA_subunit_corr), xlabel('Correlation'), ylabel('Frequency'), title('Check & Subunit spatiochromotemporal correlation');
%% works if you wanna see the data for 1 cell

select_cell_idx = 97;
figure(2), subplot(211); plot(STAs{select_cell_idx,1}','Linewidth',2);
subplot(212); plot(STAs{select_cell_idx,2}','Linewidth',2);

% Performing a head to head comparison of the the two STAs
muvect = repmat([.5 ;.5; .5],nstixperside^2,1);
FigHandle = figure(3);
set(FigHandle, 'Position', [50, 70, 1449, 705]);
for jj = 1:2
    tmp_STA = STAs{select_cell_idx,jj};
    for i = 1:size(STAs{select_cell_idx,1},2) % evaluates for each frame
        normfactor = 0.5/((max(abs(tmp_STA(:)))) + 0.0001); % This step is necessary to constrain the values within [-0.5, 0.5]
        STA = tmp_STA(:,i);
        STA = expand_vector(STA,2,mask,1);
        STA = normfactor*(STA)+muvect;  % This makes the values fall back within a range of 0 and 1.
        STA = reshape(STA,[nstixperside nstixperside 3]); % Decomposing the STA into R, G and B plane
        subplot(2,size(tmp_STA,2),(jj-1)*size(STAs{1},2)+i);
        image(STA); % for viewing the image
        set(gca,'XTick',[],'YTick',[]); axis square;
        if (i == 1)
            ylabel('STA');
        end
    end
end

