% Writing the new script for understanding how the perception of an edge
% changes with changes in the daylight spectra.
% Author - Abhishek De, 1/17
% Daylight spectra - sky.asc

% Code begins here
% Loading the illumination spectra
clearvars; close all;
load sky.asc
wave_il = 390:4:1070; % specific to this file
illuminant = [];
hi_wave_ind = find([390:5:1070]==780);
for ii = 1:size(sky,1)
    tmp = sky(ii,:);
    splinefit = spline(wave_il, tmp, [390:5:1070]);
    illuminant = [illuminant; splinefit(1:hi_wave_ind)];
end
% The illuminant has the illumination spectras ranging from 390 to 780 in steps of 5nm.

% Loading the reflectance spectra of two surfaces and place them adjacent to each other
load munsell380_800_1.mat % rows - wavelength, columns - different reflectances
L = numel(380:1:780); % for the munsell reflectance spectras - ftp://ftp.cs.joensuu.fi/pub/color/spectra/mspec/README.txt
ind = 1:5:L;
ind = ind(3:end);
rand_idxs = randi(size(munsell,2),[2 1]); % choose 2 random numbers from 1269 possible reflectance surfaces
reflectance_spectra1 = munsell(ind,rand_idxs(1));
reflectance_spectra1 = reflectance_spectra1';
reflectance_spectra2 = munsell(ind,rand_idxs(2));
% reflectance_spectra2 = reflectance_spectra2';
reflectance_spectra2 = reflectance_spectra1;

% Loading the cone action spectra and monitor spectral distributions
load fundamentals.mat
load mon_spd;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);  % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
fundamentals = fundamentals(3:end,:);
mon_spd = mon_spd(:,3:end);
M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone excitations

hi = size(illuminant,1);
% Declaring the image plane
L_cone_exc_ratios = zeros(1,hi);
M_cone_exc_ratios = zeros(1,hi);
S_cone_exc_ratios = zeros(1,hi);
RGBsurf1 = [];
RGBsurf2 = [];

wave = 390:5:780;
for ii = 1:1:hi
    % Calculating the net light energy entering the eye of an edge
    net_spectra1 = illuminant(ii,:).*reflectance_spectra1;
    net_spectra2 = illuminant(ii,:).*reflectance_spectra2;
    
    % Calculating the L,M,S excitations/quantal catches from the incident light
    L_exc1 = net_spectra1 * fundamentals(:,1);
    M_exc1 = net_spectra1 * fundamentals(:,2);
    S_exc1 = net_spectra1 * fundamentals(:,3);
    L_exc2 = net_spectra2 * fundamentals(:,1);
    M_exc2 = net_spectra2 * fundamentals(:,2);
    S_exc2 = net_spectra2 * fundamentals(:,3);
    
    % Calculating the RGB intensities from the L,M,S excitations after
    % normalizng the responses, 0 - don't normalize, 1 - do normalize
    [L_norm1,M_norm1,S_norm1,L_norm2,M_norm2,S_norm2] = normalize_resp(L_exc1,M_exc1,S_exc1,L_exc2,M_exc2,S_exc2,0);
    RGB1 = inv(M)*([L_norm1; M_norm1; S_norm1]);
    RGB_norm1 = (RGB1)/norm(RGB1);
    RGB2 = inv(M)*([L_norm2; M_norm2; S_norm2]);
    RGB_norm2 = (RGB2)/norm(RGB2);
    L_cone_exc_ratios(ii) = L_norm1/L_norm2;
    M_cone_exc_ratios(ii) = M_norm1/M_norm2;
    S_cone_exc_ratios(ii) = S_norm1/S_norm2;
    RGBsurf1 = [RGBsurf1; RGB_norm1'];
    RGBsurf2 = [RGBsurf2; RGB_norm2'];   
end


%% Try projecting the incoming light onto the STAs of the cells from your WNthresh project and see if the predictions match
% the isoresponse contour acquired from your experiment
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM WNthresh');
comments = fetch(conn,'SELECT comments FROM WNthresh');
NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
close(conn);
subunitexp_indices = find(strcmp(NTmode,'subunit')==1);
global spikename maskidx spikeidx neurothreshidx nstixperside ngammasteps seedidx nframesidx correctidx
global fponidx stimoffidx stimonidx muidxs sigmaidxs basisvecidx weightsidx fpacqidx targetspikerateidx basisvecdiridx latencyidx reversalflagidx parentverticesidx
global msperframe ntrials maxT xx yy linepredtol stepsizescale stepsize nreversals oogscale

fig_idx = 3;
fact = 0;
num_neurons = 5;
rel_filenames = filename(subunitexp_indices);
% rel_filenames = rel_filenames(1:15);
num_figs = numel(rel_filenames)/num_neurons;
color = ['g-';'k-';'r-';'b-';'m-'];
lo = -5.0; hi = 5.0;
for bb = fig_idx:fig_idx+ceil(num_figs)
    figure(bb);
    fact = fact + 1 ;
    end_num = fact*num_neurons;
    if end_num > numel(rel_filenames)
        end_num = numel(rel_filenames);
    end
    for aa = (fact-1)*num_neurons+1:end_num
        stro = nex2stro(findfile(char(rel_filenames(aa))));
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
        mask_span = mask_changes(:,2);
        
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
        basis_vec = STAs;
        tmp_basis = basis_vec;
        frametimes = [1:maxT]*msperframe;
        basisvec_energy =  sum(tmp_basis.^2);
        max_energy_idx = find(basisvec_energy == max(basisvec_energy(1,1:10)));
        desired_im = tmp_basis(:,max_energy_idx);
        
        RF = zeros(nstixperside*nstixperside,3);
        S1 = find(org_mask==1);
        S2 = find(org_mask==2);
        RF(S1,1) = desired_im(1); RF(S1,2) = desired_im(3); RF(S1,3) = desired_im(5);
        RF(S2,1) = desired_im(2); RF(S2,2) = desired_im(4); RF(S2,3) = desired_im(6);
        RGBS1  = zeros(3,1); RGBS2  = zeros(3,1);
        RGBS1(1) = desired_im(1); RGBS1(2) = desired_im(3); RGBS1(3) = desired_im(5);
        RGBS2(1) = desired_im(2); RGBS2(3) = desired_im(4); RGBS2(3) = desired_im(6);
        Subunit_act1 = RGBsurf1 * RGBS1;
        Subunit_act2 = RGBsurf2 * RGBS2;
        Subunit_act3 = RGBsurf1 * RGBS2;
        Subunit_act4 = RGBsurf2 * RGBS1;
        rev_aa = aa - (fact-1)*num_neurons;
        subplot(num_neurons,5,5*rev_aa-4); image(reshape((0.5*RF/(max(abs(RF(:)))+0.01)+0.5),[nstixperside nstixperside 3])); set(gca,'XTick',[],'YTick',[]); axis square;
        subplot(num_neurons,5,5*rev_aa-3), plot(Subunit_act1, Subunit_act2,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1]); axis square;
        %set(gca,'XTick',[],'YTick',[],'gridlinestyle','-');
        xlabel('S1'); ylabel('S2'); grid on;
        subplot(num_neurons,5,5*rev_aa-2), plot(Subunit_act3, Subunit_act4,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1]); axis square;
        set(gca,'XTick',[],'YTick',[],'gridlinestyle','-'); xlabel('S1'); ylabel('S2'); grid on;
        
        neurothreshmode = stro.trial(:,neurothreshidx);
        basisvec_dropidx = inds(end);
        neurothresh_startidx = max(find(neurothreshmode==1,1),basisvec_dropidx+1);
        num_targetspikerates = unique(stro.trial(neurothresh_startidx:end,targetspikerateidx));
        % plotting the rasterplots for neurothresh trials
        norms = cell(1,1);
        completed_search_alongdir = cell(1,1);
        vect = stro.ras{basisvec_dropidx,basisvecidx};
        t_offset = stro.trial(end,latencyidx)/1000;
        basisvec_size = nstixperside*nstixperside*3;
        numvect = (numel(vect)/basisvec_size)-1;
        basisvec = cell(1,numvect);
        % Actual basis vec
        for ii = 1:numvect
            tmp_vec = vect((ii-1)*basisvec_size+1:basisvec_size*ii) + vect(numvect*basisvec_size+1:basisvec_size*(numvect+1));
            basisvec{ii} = reshape(tmp_vec,[nstixperside nstixperside 3]); 
        end
        bkgnd_monitor = reshape(vect(numvect*basisvec_size+1:basisvec_size*(numvect+1)),[nstixperside nstixperside 3]);
        for jj = 1:1
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
        color = ['g','k'];
        for ii = 1:size(norms,2)
            tmp = norms{ii};
            completed_dir = completed_search_alongdir{ii};
            tmp = tmp(completed_dir(:,1)==1,:);
            oog_idx = find(completed_dir(:,2)==1);
            [THETA,RHO] = cart2pol(tmp(:,1),tmp(:,2));
            ind = (1:numel(THETA))';
            r = fliplr(linspace(0,1,numel(ind)));
            b = fliplr(r);
            THETA = THETA * (180/pi);
            % Earlier points in time are blue in color and later points in time are red in color
            subplot(num_neurons,5,5*rev_aa-1); 
            for jj = 1: numel(ind)
                plot(tmp(ind(jj),1), tmp(ind(jj),2),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[r(jj) 0 b(jj)]); hold on;
            end
            grid on;
            set(gca,'Xlim',[lo,hi],'Ylim',[lo,hi]); colormap(gca,winter); 
            axis square; plot(tmp(ind(ismember(ind,oog_idx)),1), tmp(ind(ismember(ind,oog_idx)),2),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor','g'); set(gca,'XTick',[],'YTick',[]);hold off;
        end
        subplot(num_neurons,5,5*rev_aa); bar([RGBS1,RGBS2]); axis square; set(gca,'YTick',[]);
        drawnow;
    end
end