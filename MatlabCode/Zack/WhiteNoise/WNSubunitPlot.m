function WNSubunitPlot(stro)
spikename = getSpikenum(stro);
maskidx = strcmp(stro.sum.rasterCells, 'subunit_mask');
spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
nstixperside = stro.sum.exptParams.nstixperside;
ngammasteps = 2^16;
seedidx = strcmp(stro.sum.trialFields(1,:),'seed');
nframesidx = strcmp(stro.sum.trialFields(1,:),'num_frames');
stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
muidxs = [find(strcmp(stro.sum.trialFields(1,:),'mu1')), ...
    find(strcmp(stro.sum.trialFields(1,:),'mu2')), ...
    find(strcmp(stro.sum.trialFields(1,:),'mu3'))];
sigmaidxs = [find(strcmp(stro.sum.trialFields(1,:),'sigma1')), ...
    find(strcmp(stro.sum.trialFields(1,:),'sigma2')), ...
    find(strcmp(stro.sum.trialFields(1,:),'sigma3'))];

msperframe = 1000/stro.sum.exptParams.framerate;
ntrials = size(stro.trial,1);

maxT = 9;
xx = linspace(stro.sum.exptParams.gauss_locut/1000, stro.sum.exptParams.gauss_hicut/1000, ...
    ngammasteps);
yy = norminv(xx');

% determine when the mask changed
mask_changes = 1;
all_masks = stro.ras(:,maskidx);
for k = 2:ntrials
    if isequal(all_masks{k}, all_masks{k-1})
        continue
    else
        mask_changes = [mask_changes k-1 k]; %#ok<AGROW>
    end
end
% each column of `mask_changes` holds the start and end trial index where each
% mask was used
mask_changes = reshape([mask_changes ntrials], 2, []);

% get rid of transient mask changes from old files
single_trials = diff(mask_changes, 2) == 0;
mask_changes(:,single_trials) = [];

hfig = figure();
set(hfig, 'Position', [500 0 850 size(mask_changes,2)*100]);

plot_counter = 0;
for mask_span = mask_changes
    STCOVmex('init', {nstixperside^2 3 maxT});
    plot_counter = plot_counter + 1;
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
            nrandnums_perchannel = length(subunitIdxs)-any(isinf(subunitIdxs)); % nsubunits
            mask = [mask; mask+max(mask); mask+2*max(mask)]; %#ok<AGROW>
        else
            nrandnums_perchannel = nstixperside^2;
        end
        
        % assuming Gaussian gun noise only
        invnormcdf = bsxfun(@plus, bsxfun(@times, yy, sigma), mu);
        randnums = getEJrandnums(3*nrandnums_perchannel*nframes, seed);
        randnums = reshape(randnums, 3*nrandnums_perchannel, nframes);
        for gun = 1:3
            idxs = (1:nrandnums_perchannel)+nrandnums_perchannel*(gun-1);
            randnums(idxs,:) = reshape(invnormcdf(randnums(idxs,:)+1,gun),[length(idxs) nframes]);
        end
        if any(org_mask) % subunit regime
            if any(isinf(subunitIdxs))
                % intersperse mu's RGB components (the background) at the end of
                % each run of Rs, Gs, and Bs
                randnums = [
                    randnums(1:nrandnums_perchannel,:);
                    mu(1)+zeros(1,nframes);
                    randnums(nrandnums_perchannel+1:2*nrandnums_perchannel,:);
                    mu(2)+zeros(1,nframes);
                    randnums(2*nrandnums_perchannel+1:3*nrandnums_perchannel,:);
                    mu(3)+zeros(1,nframes);
                    ];
            end
            rgbs = randnums(mask,:);
        else
            rgbs = randnums;
        end
        
        t_stimon = stro.trial(k, stimonidx);
        spiketimes = (stro.ras{k,spikeidx}-t_stimon)*1000;
        frametimes = linspace(0, nframes*msperframe, nframes)+(msperframe/2)';
        spiketimes(spiketimes < maxT*msperframe) = [];
        spiketimes(spiketimes > frametimes(end)) = [];
        n = hist(spiketimes, frametimes);
        STCOVmex(rgbs(:),n);
    end
    
    out = STCOVmex('return');
    clear STCOVmex;
    
    STAs = out{1};
    clear out;
    STAs = STAs./(2*max(abs(STAs(:))))+.5;
    STAs = reshape(STAs,[nstixperside nstixperside 3 maxT]);
    
    htime = subplot(size(mask_changes,2), maxT+1, plot_counter*(maxT+1), 'Parent', hfig);
    set(htime,'TickLength',[0 0],'YTick',[],'XTickLabel',[],'XLim',[.5 maxT+.5]); axis(htime, 'square');
    for j = 1:size(STAs,4)
        subplot(size(mask_changes,2), maxT+1, (plot_counter-1)*(maxT+1)+j, 'Parent', hfig);
        image(STAs(:,:,:,j), 'ButtonDownFcn', {@paint_time_course,hfig,htime,maxT,STAs,nstixperside});
        set(gca,'XTick',[],'YTick',[]); axis square;
    end
    drawnow();
end

function paint_time_course(h, ~, hfig, htime, maxT, STAs, nstixperside)
st = get(hfig, 'SelectionType');
whichpt = get(get(h, 'Parent'), 'CurrentPoint');
if strcmp(st, 'normal')
    whichpt = round(whichpt(1,1:2));
    whichpt = min(max(1, whichpt), nstixperside);
    hlines = plot(htime, squeeze(STAs(whichpt(2),whichpt(1),:,:))', 'linewidth', 2);
    set(hlines, {'Color'}, {'r' 'g' 'b'}');
    set(htime,'TickLength',[0 0],'YTick',[],'XTickLabel',[],'XLim',[.5 maxT+.5]);
    axis(htime, 'square');
end
