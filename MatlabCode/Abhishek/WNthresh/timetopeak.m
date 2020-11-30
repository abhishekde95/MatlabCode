% Author - Abhishek De, 08/16
% Latency to peak energy
clear all;
close all;
echo off;

%***********************************
% Data from NUT
%***********************************
file_info_wn1 = {'Abhi_V1_Nut\N021915003.nex', 'STA';...
                 'Abhi_V1_Nut\N022715002.nex', 'STA';...
                 'Abhi_V1_Nut\N040415002.nex', 'STA';...
                 'Abhi_V1_Nut\N050415001.nex', 'STA';...
                 'Abhi_V1_Nut\N041115001.nex', 'STA';...
                 'Abhi_V1_Nut\N050215001-header-merge.nex', 'STA';...
                 'Abhi_V1_Nut\N050215002.nex', 'STA';...
                 'Abhi_V1_Nut\N050415002.nex', 'STA';... % has blue STA, PC1
                 'Abhi_V1_Nut\N060515002.nex', 'PC1';...
                 'Abhi_V1_Nut\N060615002.nex', 'STA';... % blue yellow STA
                 'Abhi_V1_Nut\N060615005.nex', 'STA';... % blue STA, luminance PC1
                 'Abhi_V1_Nut\N060815001.nex', 'STA';... % blue STA
                 'Abhi_V1_Nut\N060815006.nex', 'STA';... % red-green DO 
                 'Abhi_V1_Nut\N061615001.nex', 'STA';... % red-green DO
                 'Abhi_V1_Nut\N061615003.nex', 'PC1';...
                 'Abhi_V1_Nut\N061715002.nex', 'STA';... % blue STA
                 'Abhi_V1_Nut\N061715004.nex', 'STA';...
                 'Abhi_V1_Nut\N062315002.nex', 'STA';... %  yellowish STA
                 'Abhi_V1_Nut\N062915002.nex', 'PC1';... 
                 'Abhi_V1_Nut\N063015002.nex', 'STA';... % blue-yellow STA, PC1, PC2
                 'Abhi_V1_Nut\N063015003.nex', 'PC1';...
                 'Abhi_V1_Nut\N070115001.nex', 'PC1';... % blue spot in STA, PC1
                 'Abhi_V1_Nut\N070715005.nex', 'STA';...
                 'Abhi_V1_Nut\N071015002.nex', 'PC1';... % blue signature STA, PC1
                 'Abhi_V1_Nut\N071015003.nex', 'PC1';... % STA, PC1, PC2
                 'Abhi_V1_Nut\N071015004.nex', 'PC1';...
                 'Abhi_V1_Nut\N071615002.nex', 'STA';...
                 'Abhi_V1_Nut\N071715001.nex', 'PC1';...
                 'Abhi_V1_Nut\N071815001.nex', 'STA';...
                 'Abhi_V1_Nut\N072215003.nex', 'PC1';...
                 'Abhi_V1_Nut\N072615002.nex', 'STA';...
                 'Abhi_V1_Nut\N072815001.nex', 'STA';... % orange-yellow smudge
                 'Abhi_V1_Nut\N073015003.nex', 'STA';... % red-green DO 
                 'Abhi_V1_Nut\N073015004.nex', 'STA';...
                 'Abhi_V1_Nut\N080715002.nex', 'STA';... % could be a blue-yellow DO cell
                 'Abhi_V1_Nut\N080715003.nex', 'PC1';...
                 'Abhi_V1_Nut\N081015001.nex', 'STA';...
                 'Abhi_V1_Nut\N081115002.nex', 'STA';...
                 'Abhi_V1_Nut\N081215002.nex', 'PC1';...
                 'Abhi_V1_Nut\N081515001.nex', 'STA';...
                 'Abhi_V1_Nut\N081515003.nex', 'STA';...
                 'Abhi_V1_Nut\N081715001.nex', 'STA';...
                 'Abhi_V1_Nut\N081715002.nex', 'STA';... % also has PC1
                 'Abhi_V1_Nut\N082015002.nex','STA'};
    
%***********************************
% Data from PANGU
%***********************************
% Neurothresh files
file_info_nt = {'Abhi_V1_Pangu_WhiteNoiseThresh\P050116001.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P052016002.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P052316005.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P052616001.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P052816003.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P052916001.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P072716002.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P072816001.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P072816004.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P080116002.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P080216003.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P080316001.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P080516003.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P080716001.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P080816001.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P080816003.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P081216004.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P081516002.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P081616002.nex', 'STA'}; % blue-yellow STA, luminance PC1

% Neurothresh and WhiteNoise files
file_info_wn2 = {
    'Abhi_V1_Pangu_WhiteNoiseThresh\P032916001.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P033116001.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P042816002.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P051916001.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P051616005.nex', 'PC1';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P052316004.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P052816001.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P072816002.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P080216001.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P080216002.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P080616001.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P081016004.nex', 'STA';...
    'Abhi_V1_Pangu_WhiteNoiseThresh\P081116001.nex', 'STA'};

file_info = cat(1,file_info_wn1,file_info_nt, file_info_wn2);
online_latency = zeros(size(file_info,1),1);
offline_latency = zeros(size(file_info,1),1);
filepath = 'C:\Users\Abhishek\Documents\MATLAB\';
basisvec_energy = [];
normalized_basisvec_energy = [];

for ii = 1:size(file_info,1)
    stro = nex2stro(strcat(filepath,file_info{ii,1}));
    spikename = getSpikenum(stro);
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
            mask_changes = [mask_changes k-1 k]; %#ok<AGROW>
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
    STCs = (nspikes.*STCross-tmp)/(nspikes*(nspikes-1));
    P = eye(size(STCs)) - STAs(:)*inv(STAs(:)'*STAs(:))*STAs(:)'; % WHAT DOES THIS LINE MEAN
    STCs = P*STCs*P';
    [tmp,d] = eig(STCs);
    v = real(tmp);
    [~, idxs] = sort(diag(d));
    v = v(:,idxs);
    PCs = v(:,end); 
    
    % Flipping the STAs and the PCs such that the last frame appears first and the first frame appears last
    PCs = fliplr(reshape(PCs,[nrandnums_perchannel*3 maxT]));
    STAs = fliplr(STAs);
    
    if strcmp(file_info{ii,2},'STA')
        tmp_basis = STAs;
    else
        tmp_basis = PCs;
    end
    frametimes = [1:maxT]*msperframe;
    basisvec_energy =  [basisvec_energy; sum(tmp_basis.^2)];
    normalized_basisvec_energy = [normalized_basisvec_energy; basisvec_energy(ii,:)/max(basisvec_energy(ii,:))];
    offline_latency(ii) = frametimes(basisvec_energy(ii,:) == max(basisvec_energy(ii,1:10)));
%     online_latency(ii) = stro.trial(end,latencyidx);
end

mean_offline_latency = mean(offline_latency);
std_offline_latency = std(offline_latency);

% Plotting the results
figure(1), subplot(131);hist(offline_latency), xlabel('Time to reach peak'), ylabel('Frequency'); title('Offline Latency');
subplot(132); plot(frametimes,normalized_basisvec_energy','Linewidth',2), xlabel('Time in ms'), ylabel('Energy'), title ('Normalized energy');
subplot(133), errorbar(frametimes,mean(normalized_basisvec_energy,1),std(normalized_basisvec_energy,1)/size(file_info,1),'Linewidth',2);
xlabel('Time in ms'); ylabel('Energy'), title('Mean normalized energy with SE'); 
%subplot(122); hist(online_latency), xlabel('Time to reach peak'), ylabel('Frequency'); title('Online Latency');
