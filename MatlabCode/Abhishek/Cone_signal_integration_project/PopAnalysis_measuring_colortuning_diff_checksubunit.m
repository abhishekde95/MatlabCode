% A script for doing some preliminary analysis for cone signal integration project
% Author - Abhishek De, 9/20,
% Some of the code can be found in
% WNthresh/PopAnalysis_WNthresh_anglecalc.m

% Notes: This will be a more rigorous analysis than I have previously attempted
% 1) First calculate the weighted STA from checkerboad and subunit whitenoise stimuli.
% 2) Store the RGB values from individual subfields for checkerboard and subunit whitenoise stimuli
% 3) Do the SVD on the weighted STAs. Will be used for comparing the color tuning between the two STAs

close all; clearvars;
plot_counter = 1;

% Loading the files 
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM WNthresh');
NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
close(conn);
filename = filename(strcmp(string(NTmode),"subunit"));
NTmode = NTmode(strcmp(string(NTmode),"subunit"));
spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));

RGBsubunits = cell(1,numel(filename));
RGBcheck = cell(1,numel(filename));
latencySTAcheck = [];
latencySTAsubunit = [];
latencydiffWNchecksubunit = [];
conewtscheck_svd = [];
conewtssubunit_svd = [];
anglediffWNchecksubunit = []; % Stores the differences in color tuning - angle between 2 6-D RGB vectors
anglediffWNchecksubunit_SVD_LMS = []; % Stores the diffeneces in color tuning - angle between 2 3-D LMS vectors 

for ii = 1:numel(filename)
    ind = ii;
    fileofinterest = char(filename(ind,:));
    disp(ii);

    stro = nex2stro(findfile(fileofinterest));
    spikename = 'sig001a'; %getSpikenum(stro);
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
    gammaTable = stro.sum.exptParams.gamma_table;
    gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
    
    % Calculating the M matrix for files
    fundamentals = stro.sum.exptParams.fundamentals;
    mon_spd = stro.sum.exptParams.mon_spd;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    
    % Getting the background rgb/lms
    ridx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_r'));
    gidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_g'));
    bidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_b'));
    bkgndRGB = [mode(stro.trial(:,ridx)), mode(stro.trial(:,gidx)), mode(stro.trial(:,bidx))];
    bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
    bkgndlms = M*bkgndrgb;
    
    Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
    Mrgbtocc = inv(Mrgbtocc');
    
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
    for jj = 1:2 % 1st column is for checkerboard noise, 2nd column is for the subunit noise
        trial_span = mask_changes(:,jj);
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
        nspikes = out{3};
        
        if jj == 1
            STAcheck = STS/nspikes;
            STAcheck = fliplr(STAcheck);
            [~,whichframe_check] = max(sum(STAcheck.^2,1));
            id = max(sum(STAcheck.^2,1) == max(sum(STAcheck.^2,1)));
    
            if whichframe_check~=1
                id(whichframe_check-1)= 1;
            end
            if whichframe_check <=maxT-1
                id(whichframe_check+1)=1;
            end
            
            STAweights = sqrt(sum(STAcheck(:,id).^2));
            STAweights = STAweights./sum(STAweights);
            weighted_STAcheck = STAcheck(:,id)*STAweights';
            
            
        else
            STAsubunit = STS/nspikes; % Storing the STA subunit
            STAsubunit = fliplr(STAsubunit); 
            [~,whichframe_subunit] = max(sum(STAsubunit.^2,1));
            id = max(sum(STAsubunit.^2,1) == max(sum(STAsubunit.^2,1)));
            
            if whichframe_subunit~=1
                id(whichframe_subunit-1)= 1;
            end
            if whichframe_subunit <=maxT-1
                id(whichframe_subunit+1)=1;
            end
            
            STAweights = sqrt(sum(STAsubunit(:,id).^2));
            STAweights = STAweights./sum(STAweights);
            weighted_STAsubunit = STAsubunit(:,id)*STAweights';
            
            % SVD on chekerboard STA 
            [u,~,~] = svd(reshape(weighted_STAsubunit,[2 3])');
            conewtssubunit_svd = [conewtssubunit_svd Mrgbtocc * u(:,1)];
            
        end
        clear STCOV_st out
    end
    
    % Calculating stuff latency and the color tuning of individual subunits
    RGBsubunits{ii} = [weighted_STAsubunit(1:2:5); weighted_STAsubunit(2:2:6)];
    latencySTAsubunit = [latencySTAsubunit; whichframe_subunit*msperframe];
    
    % Calculating latency STA check
    latencySTAcheck = [latencySTAcheck; whichframe_check*msperframe];
    
    % Latency differences between the Checkerboard and Subunit WN stim 
    latencydiffWNchecksubunit = [latencydiffWNchecksubunit; (whichframe_subunit-whichframe_check)*msperframe];
    subunit1R = mean(weighted_STAcheck(mask==1)); subunit2R = mean(weighted_STAcheck(mask==2));
    subunit1G = mean(weighted_STAcheck(mask==4)); subunit2G = mean(weighted_STAcheck(mask==5));
    subunit1B = mean(weighted_STAcheck(mask==7)); subunit2B = mean(weighted_STAcheck(mask==8));
    
    % SVD on chekerboard STA 
    [u,~,~] = svd([subunit1R subunit2R; subunit1G subunit2G; subunit1B subunit2B]);
    conewtscheck_svd = [conewtscheck_svd Mrgbtocc * u(:,1)];
    
    % Differences between color tuning between the checkerboard and subunit WN stim - 2 6-D RGB vectors 
    RGBcheck{ii} = [subunit1R; subunit1G; subunit1B; subunit2R; subunit2G; subunit2B];
    tmp = min([180*acos(dot(RGBcheck{ii},RGBsubunits{ii})/(norm(RGBcheck{ii})*norm(RGBsubunits{ii})))/pi; 180*acos(dot(RGBcheck{ii},-1*RGBsubunits{ii})/(norm(RGBcheck{ii})*norm(-1*RGBsubunits{ii})))/pi]);
    anglediffWNchecksubunit = [anglediffWNchecksubunit; tmp ];
    
    % Differences between color tuning between the checkerboard and subunit WN stim - 2 3-D LMS vectors
    tmp = min([180*acos(dot(conewtscheck_svd(:,end),conewtssubunit_svd(:,end))/(norm(conewtscheck_svd(:,end))*norm(conewtssubunit_svd(:,end))))/pi; 180*acos(dot(conewtscheck_svd(:,end),-1*conewtssubunit_svd(:,end))/(norm(conewtscheck_svd(:,end))*norm(-1*conewtssubunit_svd(:,end))))/pi]);
    anglediffWNchecksubunit_SVD_LMS = [anglediffWNchecksubunit_SVD_LMS; tmp];
end


savevariables = 1;
if savevariables 
    save RGBsubunits RGBsubunits
    save RGBcheck RGBcheck 
    save latencySTAcheck latencySTAcheck
    save latencySTAsubunit latencySTAsubunit
    save latencydiffWNchecksubunit latencydiffWNchecksubunit 
    save conewtscheck_svd conewtscheck_svd 
    save conewtssubunit_svd conewtssubunit_svd 
    save anglediffWNchecksubunit anglediffWNchecksubunit
    save anglediffWNchecksubunit_SVD_LMS anglediffWNchecksubunit_SVD_LMS 
end

%% Analysis of color tuning 

% Plotting the SVD derived cone wts for checkerboard and subunit whitenoise

LMS_check = repmat(sign(conewtscheck_svd(2,:)),[3 1]).*(conewtscheck_svd./repmat(sum(abs(conewtscheck_svd),1),[3 1]));
LMS_subunit = repmat(sign(conewtssubunit_svd(2,:)),[3 1]).*(conewtssubunit_svd./repmat(sum(abs(conewtssubunit_svd),1),[3 1]));

figure(plot_counter);
subplot(221); plot(LMS_check(1,:), LMS_check(2,:),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');  
title ('Checkerboard SVD'); xlabel('L'), ylabel('M'); hold off;
subplot(222); plot(LMS_subunit(1,:), LMS_subunit(2,:),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');  
title ('Subunit SVD'); xlabel('L'), ylabel('M'); hold off;

% Checking the difference between the color tuning between subunit and checkerboard noise using 2 different methods 
subplot(223); plot(anglediffWNchecksubunit, anglediffWNchecksubunit_SVD_LMS, 'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot([0 90],[00 90],'k'); hold off; axis square; set(gca,'Tickdir','out','Xlim',[0 90],'Ylim',[0 90],'XTick',0:30:90,'YTick',0:30:90); xlabel('6-D RGB method'); ylabel('3-D SVD LMS method');
title('Changes in color tuning');

% Latency measurements using checkerboard and subunit whitenoise STA
subplot(224); plot(latencySTAcheck, latencySTAsubunit, 'o', 'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot([0 200],[0 200],'k'); axis square; set(gca,'Tickdir','out', 'Xlim',[0 200], 'Ylim',[0 200],'XTick',0:50:200, 'YTick',0:50:200); 
ylabel('Latency subunit STA'); xlabel('Latency check STA'); title('Latency Analysis'); hold off;

plot_counter = plot_counter + 1;

[r,p] = corr(anglediffWNchecksubunit,anglediffWNchecksubunit_SVD_LMS,'type','Spearman');

