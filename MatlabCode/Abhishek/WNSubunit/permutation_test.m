function [plot_counter,prctileval,z_scoreval] = permutation_test(stro,mask_changes,plot_counter,eig_val,plotmode,whichframe)
% Perform the permutation test
% Made some changes - Abhishek, 1/18
global maskidx spikeidx nstixperside seedidx nframesidx stimonidx muidxs sigmaidxs msperframe maxT yy
vals = []; % To store the eigenvalues for all the iterations of the permutation test
trial_span = mask_changes(:,2);
num_iters = 1000;
for iter = 1:num_iters
    for mask_span = trial_span
        % Just store enough space to accomodate the 9 frames of the basis vector
        st_mask = stro.ras{trial_span(1),maskidx}; % subunit mask
        st_mask(st_mask == 0) = Inf;
        [stIdxs,~,~] = unique(st_mask); % now the Infs map to nsubunits+1
        num_stunits = length(stIdxs)-any(isinf(stIdxs)); % nsubunits, like subunits A and B
        STCOVmex('init', {num_stunits 3 maxT});
    end
    for k = mask_span(1):mask_span(2)
        nframes = stro.trial(k,nframesidx);
        if nframes == 0, continue; end
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
        % Deleting the spikes taking place before the first 9 frames as I need to look at the 9 preceding frames
        spiketimes(spiketimes < maxT*msperframe) = [];
        spiketimes(spiketimes > frametimes(end)) = [];
        n = hist(spiketimes, frametimes);
        STCOVmex(rgbs(:),circshift(n',randi(numel(n))));
    end
    out = STCOVmex('return');
    STS = out{1};  % A 6 x 9(frames) matrix
    STCross = out{2};  % A 36 x 9(frames) matrix
    nspikes = out{3}; % Number of spikes in the given file
    clear STCOVmex;
    clear out;
    STAs = STS/nspikes;
    tmp = STS(:,whichframe)*STS(:,whichframe)';
    STCross = reshape(STCross(:,whichframe),[num_stunits*3 num_stunits*3]);
    STCs = (nspikes.*STCross-tmp)/(nspikes*(nspikes-1));
    P = eye(size(STCs)) - STAs(:,whichframe)*inv(STAs(:,whichframe)'*STAs(:,whichframe))*STAs(:,whichframe)'; % Subtracting such that the STA is orthogonal to PC
    STCs = P*STCs*P';
    [~,d] = eig(STCs);
    eig_PC = sort(diag(d));
    vals = [vals; eig_PC(end)]; % storing all the eigenvalues
    
end
if plotmode 
    figure(plot_counter); hist(vals,30); hold on; plot(eig_val,0,'kv','MarkerFaceColor','g');
    plot([prctile(vals,95) prctile(vals,95)],[0 50],'-k','Linewidth',4); xlabel('eig values'); ylabel('Frequency');
    title('Permutation test'); hold off;
    plot_counter = plot_counter + 1;
end
prctileval = (find(sort(vals)>eig_val,1)/num_iters)*100;
z_scoreval = (eig_val - mean(vals))/std(vals);
if isempty(prctileval)
    prctileval = 100;
end
clear vals
end

