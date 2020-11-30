% This is an new script which is very similar to th PopAnalysis_WNthresh_comparebetween two methodst
% Author - Abhishek De, 1/18
close all; clearvars;
plot_counter = 1;
load filename_c.mat % Color cells
load filename_l.mat % Luminance cells
filename = [filename_c; filename_l];
plot_counter = 1;
STAsub1 = cell(numel(filename),1);
STAsub2 = cell(numel(filename),1);
STCsub1 = cell(numel(filename),1);
STCsub2 = cell(numel(filename),1);
STAcheck1 = cell(numel(filename),1);
STAcheck2 = cell(numel(filename),1);
STCcheck1 = cell(numel(filename),1);
STCcheck2 = cell(numel(filename),1);
totspikes = zeros(numel(filename),2);
azimuthcheck1 = cell(numel(filename),1);
elevationcheck1 = cell(numel(filename),1);
azimuthcheck2 = cell(numel(filename),1);
elevationcheck2 = cell(numel(filename),1);
azimuthsub1 = cell(numel(filename),1);
elevationsub1 = cell(numel(filename),1);
azimuthsub2 = cell(numel(filename),1);
elevationsub2 = cell(numel(filename),1);
numcells = numel(filename);
numsubplots = ceil(sqrt(numcells));
plotmode = 1; % 0 - don't plot, 1 - plot
mahalanobisdist = zeros(numel(filename),2);
bhattacharyadist =zeros(numel(filename),2); % measures the similarity between 2 probability distributions
Fstatistic = zeros(numel(filename),2);
Pval = zeros(numel(filename),2);
threshold = 0.05;
confidence_interval = 0.99;
whichoctants = cell(numel(filename),2);

for ii = 1:numcells
    disp(ii);
    stro = nex2stro(findfile(char(filename(ii,:))));
    global spikename maskidx spikeidx nstixperside ngammasteps seedidx nframesidx
    global fponidx stimoffidx stimonidx muidxs sigmaidxs basisvecidx fpacqidx latencyidx basisvecdiridx neurothreshidx targetspikerateidx
    global msperframe ntrials maxT xx yy M
    spikename = 'sig001a';%getSpikenum(stro);
    maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
    spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
    nstixperside = stro.sum.exptParams.nstixperside;
    ngammasteps = 2^16; % 65536
    seedidx = strcmp(stro.sum.trialFields(1,:),'seed');
    nframesidx = strcmp(stro.sum.trialFields(1,:),'num_frames');
    stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
    stimoffidx = strcmp(stro.sum.trialFields(1,:),'stim_off');
    fponidx = strcmp(stro.sum.trialFields(1,:),'fp_on');
    fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
    muidxs = [find(strcmp(stro.sum.trialFields(1,:),'mu1')), ...
        find(strcmp(stro.sum.trialFields(1,:),'mu2')), ...
        find(strcmp(stro.sum.trialFields(1,:),'mu3'))];
    sigmaidxs = [find(strcmp(stro.sum.trialFields(1,:),'sigma1')), ...
        find(strcmp(stro.sum.trialFields(1,:),'sigma2')), ...
        find(strcmp(stro.sum.trialFields(1,:),'sigma3'))];
    latencyidx = strcmp(stro.sum.trialFields(1,:),'latency');
    basisvecidx = strcmp(stro.sum.rasterCells(1,:),'basis_vec');
    msperframe = 1000/stro.sum.exptParams.framerate;
    basisvecdiridx = strcmp(stro.sum.trialFields(1,:),'weights_idx');
    neurothreshidx = strcmp(stro.sum.trialFields(1,:),'neurothresh'); % when exactly the neurothresh trials started
    targetspikerateidx = strcmp(stro.sum.trialFields(1,:),'targetspikerate');
    ntrials = size(stro.trial,1);
    maxT = 15; % this represents the temporal part in the spatiotemporal receptive field
    xx = linspace(stro.sum.exptParams.gauss_locut/1000, stro.sum.exptParams.gauss_hicut/1000,ngammasteps); % xx represents the probabilities. For more info, have a look at the MATLAB 'norminv' function.
    yy = norminv(xx'); % defining norminv to extract the values for which the cdf values range between gauss_locut and gauss_hicut
    
    % Obtaining the M matrix, code extracted from Greg, fitting a cubic spline
    % using the command 'spline'. 'SplineRaw' only availabe through
    % psychtoolbox which I currently don't have now.
    if ii == 1
        fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
        fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
        mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
        mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
        mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
        M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
        M = inv(M');
    end
    
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
    sequence = [2 1]; % first subunit stimuli and then the checkerboard stimuli
    
    for jj = 1:numel(sequence)
        masktrials = mask_changes(:,sequence(jj));
        st_mask = stro.ras{masktrials(1),maskidx}; % subunit mask
        st_mask(st_mask == 0) = Inf;
        [stIdxs,~,~] = unique(st_mask); % now the Infs map to nsubunits+1
        num_subunits = length(stIdxs)-any(isinf(stIdxs)); % nsubunits, like subunits A and B
        if jj == 1
            STCOVmex('init', {num_subunits 3 maxT});
        else
            STCOVmex('init', {numel(newidxstoSTCOVmex)/3 3 maxT});
        end
        
        for k = masktrials(1):masktrials(2)
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
            if jj == 1
                rgbs = randnums;
            else
                rgbs = randnums(newidxstoSTCOVmex,:);
            end
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
            STCOVmex(rgbs(:),n);
        end
        out = STCOVmex('return'); % returns the covariance matrix on frame by frame basis
        STS = out{1};  % A (dimension) x 9(frames) matrix
        STCross = out{2};  % A (dimension x frames)x (dimension x frames)  matrix
        nspikes = out{3}; % Number of spikes in the given file
        totspikes(ii,jj) = nspikes;
        clear STCOVmex out;
        STAs = STS./nspikes;
        if jj == 1
            % WhiteNoise subunit as a stimuli
            
            [~,whichframe] = max(sum(STAs.^2,1)); % finding the frame with the maximum energy
            STStmp = STS(:,whichframe)*STS(:,whichframe)';
            STCrosstmp = reshape(STCross(:,whichframe),size(STStmp));
            STCs = (nspikes*STCrosstmp-STStmp)/(nspikes*(nspikes-1));
            STAsub1{ii} = STAs(1:2:6,whichframe);
            STAsub2{ii} = STAs(2:2:6,whichframe);
            STCsub1{ii} = [STCs(1,1) STCs(1,3) STCs(1,5); STCs(3,1) STCs(3,3) STCs(3,5); STCs(5,1) STCs(5,3) STCs(5,5)];
            STCsub2{ii} = [STCs(2,2) STCs(2,4) STCs(2,6); STCs(4,2) STCs(4,4) STCs(4,6); STCs(6,2) STCs(6,4) STCs(6,6)];
            newmask = repmat(org_mask,[3 1]);
            pixelsinsubunit1 = find(newmask==1);
            pixelsinsubunit2 = find(newmask==2);
            newidxstoSTCOVmex = [pixelsinsubunit1; pixelsinsubunit2];
        else
            % WhiteNoise checkerboard as a stimuli
            idx1 = reshape(1:1:numel(pixelsinsubunit1),[numel(pixelsinsubunit1)/3 3])';
            idx2 = reshape(numel(pixelsinsubunit1)+1:1:numel(pixelsinsubunit1)+numel(pixelsinsubunit2),[numel(pixelsinsubunit2)/3 3])';
            STAcheck1{ii} = [mean(STAs(idx1(1,:),whichframe)); mean(STAs(idx1(2,:),whichframe)); mean(STAs(idx1(3,:),whichframe))];
            STAcheck2{ii} = [mean(STAs(idx2(1,:),whichframe)); mean(STAs(idx2(2,:),whichframe)); mean(STAs(idx2(3,:),whichframe))];
            STStmp = STS(:,whichframe)*STS(:,whichframe)';
            STCrosstmp = reshape(STCross(:,whichframe),size(STStmp));
            STCs = (nspikes*STCrosstmp-STStmp)/(nspikes*(nspikes-1));
            STC1 = STCs(1:numel(pixelsinsubunit1),1:numel(pixelsinsubunit1));
%             A1 = repmat(eye(3),[1 numel(pixelsinsubunit1)/3])./(numel(pixelsinsubunit1)/3);
            A1 = zeros(3,numel(pixelsinsubunit1)); A1(1,1:numel(pixelsinsubunit1)/3) = deal(3/numel(pixelsinsubunit1));
            A1(2,(numel(pixelsinsubunit1)/3)+1:2*numel(pixelsinsubunit1)/3) = deal(3/numel(pixelsinsubunit1));
            A1(3,(2*numel(pixelsinsubunit1)/3)+1:numel(pixelsinsubunit1)) = deal(3/numel(pixelsinsubunit1));
            STCcheck1{ii} = A1*STC1*A1';
            STC2 = STCs(numel(pixelsinsubunit1)+1:end,numel(pixelsinsubunit1)+1:end);
%             A2 = repmat(eye(3),[1 numel(pixelsinsubunit2)/3])./(numel(pixelsinsubunit2)/3);
            A2 = zeros(3,numel(pixelsinsubunit2)); A2(1,1:numel(pixelsinsubunit2)/3) = deal(3/numel(pixelsinsubunit2));
            A2(2,(numel(pixelsinsubunit2)/3)+1:2*numel(pixelsinsubunit2)/3) = deal(3/numel(pixelsinsubunit2));
            A2(3,(2*numel(pixelsinsubunit2)/3)+1:numel(pixelsinsubunit2)) = deal(3/numel(pixelsinsubunit2));
            STCcheck2{ii} = A2*STC2*A2';
        end     
    end
    mfcheck1 = 1/norm(STAcheck1{ii});% multiplication factor 1 : check
    mfcheck2 = 1/norm(STAcheck2{ii});% multiplication factor 2 : check
    mfsub1 = 1/norm(STAsub1{ii});% multiplication factor 1 : subunit
    mfsub2 = 1/norm(STAsub2{ii});% multiplication factor 2 : subunit
    
    % Drawing random numbers from a nultivariate distribution
    A_new1 = mfcheck1 * eye(3);
    A_new2 = mfcheck2 * eye(3);
    A_new3 = mfsub1 * eye(3);
    A_new4 = mfsub2 * eye(3);
    STAcheck1new = A_new1*STAcheck1{ii}; STCcheck1new = A_new1*STCcheck1{ii}*A_new1';
    STAcheck2new = A_new2*STAcheck2{ii}; STCcheck2new = A_new2*STCcheck2{ii}*A_new2';
    STAsub1new = A_new3*STAsub1{ii}; STCsub1new = A_new3*STCsub1{ii}*A_new3';
    STAsub2new = A_new4*STAsub2{ii}; STCsub2new = A_new4*STCsub2{ii}*A_new4';
    
    n1 = totspikes(ii,1); n2 = totspikes(ii,2); p = 3; n = n1 + n2;
    % Drawing an ellipsoid by calculating 95% confidence interval around the mean
    pooledcovariance1 = (STCsub1new/n1) + (STCcheck1new/n2);
    pooledcovariance2 = (STCsub2new/n1) + (STCcheck2new/n2);
    pts1 = error_ellipse('C',pooledcovariance1,'mu',STAcheck1new-STAsub1new,'conf',confidence_interval);
    pts2 = error_ellipse('C',pooledcovariance2,'mu',STAcheck2new-STAsub2new,'conf',confidence_interval);
    whichoctants{ii,1} = findwhichoctant(pts1);
    whichoctants{ii,2} = findwhichoctant(pts2);
    
    mahalanobisdist(ii,1) = sqrt((STAcheck1new-STAsub1new)'*inv(pooledcovariance1)*(STAcheck1new-STAsub1new));
    mahalanobisdist(ii,2) = sqrt((STAcheck2new-STAsub2new)'*inv(pooledcovariance2)*(STAcheck2new-STAsub2new));
    bhattacharyadist(ii,1) = (((STAcheck1new-STAsub1new)'*inv(pooledcovariance1/2)*(STAcheck1new-STAsub1new))/8) + 0.5*log(det(pooledcovariance1/2)/sqrt(det(STCsub1new/n1)*det(STCcheck1new/n2)));
    bhattacharyadist(ii,2) = (((STAcheck2new-STAsub2new)'*inv(pooledcovariance2/2)*(STAcheck2new-STAsub2new))/8) + 0.5*log(det(pooledcovariance2/2)/sqrt(det(STCsub2new/n1)*det(STCcheck2new/n2)));

    if plotmode    
        % Subunit 1
        figure(plot_counter),subplot(numsubplots,numsubplots,ii),plot3(pts1(:,1),pts1(:,2),pts1(:,3),'o','MarkerSize',1,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on;        if Pval(ii,1)< 0.05
        line([-1 1],[0 0],[0 0],'Linewidth',2); line([0 0],[-1 1],[0 0],'Linewidth',2); line([0 0],[0 0],[-1 1],'Linewidth',2);
        if numel(whichoctants{ii,1}) < 8
            set(gca,'XColor','red','YColor','red','ZColor','red');
        end
        grid on; hold off;
        % Subunit 2
        figure(plot_counter + 1),subplot(numsubplots,numsubplots,ii),plot3(pts2(:,1),pts2(:,2),pts2(:,3),'o','MarkerSize',1,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on;
        line([-1 1],[0 0],[0 0],'Linewidth',2); line([0 0],[-1 1],[0 0],'Linewidth',2); line([0 0],[0 0],[-1 1],'Linewidth',2);
        if numel(whichoctants{ii,2}) < 8
            set(gca,'XColor','red','YColor','red','ZColor','red');
        end
        grid on; hold off;     
        end
    end
end
plot_counter = plot_counter + 2;
% saving variables
savevariables = 1;
if savevariables == 1
    save mahalanobisdist mahalanobisdist
    save bhattacharyadist bhattacharyadist
end

%% Using pretty corr to see any correlation among variables 
N = numel(filename_c);
load logresidualratio.mat
load anglebwvectorsRGB.mat
load absS.mat
load vals.mat
load z_scores.mat % from permutation test of eigenvalues
load mahalanobisdist.mat
load bhattacharyadist.mat
% for Color cells
prettycorr([vals(1:N),z_scores(1:N),logresidualratio(1:N),anglebwvectorsRGB(1:N),absS(1:N)',mahalanobisdist(1:N,1),mahalanobisdist(1:N,2),bhattacharyadist(1:N,1),bhattacharyadist(1:N,2)],{'eig1prctile','zscores','LogResRat','AngleRGB','absS','mahadist1','mahadist2','bhatdist1','bhatdist2'});
set(gcf,'Name','Color Cells');
% for Luminance cells 
prettycorr([vals(N+1:end),z_scores(N+1:end),logresidualratio(N+1:end),anglebwvectorsRGB(N+1:end),absS(N+1:end)',mahalanobisdist(N+1:end,1),mahalanobisdist(N+1:end,2),bhattacharyadist(N+1:end,1),bhattacharyadist(N+1:end,2)],{'eig1prctile','zscores','LogResRat','AngleRGB','absS','mahadist1','mahadist2','bhatdist1','bhatdist2'});
set(gcf,'Name','Luminance Cells');
