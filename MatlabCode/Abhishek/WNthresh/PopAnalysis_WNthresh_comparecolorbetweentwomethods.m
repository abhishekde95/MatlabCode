% This is new analysis I am writing after my discussion with Greg about comparing the RGB of the STA subunits obtained using both checkerboard and subunit stimuli.
% Author - Abhishek De, 1/18
% Use this method of classifying cells
close all; clearvars;

plot_counter = 1;
load('T_vos1978_Y');
load RHO_all.mat
load THETA_all.mat
load oog_idx_all.mat
load not_oog_idx_all.mat
load significantsubunits.mat

conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM WNthresh');
NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
close(conn);
filename = filename(strcmp(string(NTmode),"subunit"));
NTmode = NTmode(strcmp(string(NTmode),"subunit"));
spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));

significantsubunits = significantsubunits;
numcells = numel(filename);
numsubplots = ceil(sqrt(numcells));
plotmode = 1; % 0 - don't plot, 1 - plot
Vlambda = T_vos1978_Y';
plot_counter = 1;
STAsub1 = cell(numcells,1);
STAsub2 = cell(numcells,1);
STCsub1 = cell(numcells,1);
STCsub2 = cell(numcells,1);
STAcheck1 = cell(numcells,1);
STAcheck2 = cell(numcells,1);
STCcheck1 = cell(numcells,1);
STCcheck2 = cell(numcells,1);
totspikes = zeros(numcells,2);
azimuthcheck1 = cell(numcells,1);
elevationcheck1 = cell(numcells,1);
azimuthcheck2 = cell(numcells,1);
elevationcheck2 = cell(numcells,1);
azimuthsub1 = cell(numcells,1);
elevationsub1 = cell(numcells,1);
azimuthsub2 = cell(numcells,1);
elevationsub2 = cell(numcells,1);
alpha = 0.001;
confidence_interval = 1-alpha;
azimuthbins = -180:5:180;
elevationbins = -90:5:90;
rvalchecksub = zeros(numcells,1);
pvalchecksub = zeros(numcells,1);
rvalchecksub1 = zeros(numcells,1);
pvalchecksub1 = zeros(numcells,1);
whichoctants = cell(numcells,1);
doesitbelongtoluminanceoctant = zeros(numcells,1);
distance1 = zeros(numcells,1); distance2 = zeros(numcells,1);
rad1 = zeros(numcells,1); rad2 = zeros(numcells,1);
rad1LMS = zeros(numcells,1); rad2LMS = zeros(numcells,1);
WwSTAsub1dist = zeros(numcells,1); WwSTAsub2dist = zeros(numcells,1); 
distSTA12LMS = zeros(numcells,1); distSTA21LMS = zeros(numcells,1); 

% Storing the identity of cells
luminancecells = zeros(numcells,1);
DOcells = zeros(numcells,1);
SOcells = zeros(numcells,1);
hardtoclassifycells = zeros(numcells,1);
STAsub1LMS_R = cell(numcells,1);
STAsub2LMS_R = cell(numcells,1);
for ii = 1:numcells
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
    latencyidx = strcmp(stro.sum.trialFields(1,:),'latency');
    ntrials = size(stro.trial,1);
    maxT = 15; % this represents the temporal part in the spatiotemporal receptive field
    xx = linspace(stro.sum.exptParams.gauss_locut/1000, stro.sum.exptParams.gauss_hicut/1000,ngammasteps); % xx represents the probabilities. For more info, have a look at the MATLAB 'norminv' function.
    yy = norminv(xx'); % defining norminv to extract the values for which the cdf values range between gauss_locut and gauss_hicut
    
    % Calculating M matrix and projection of Vlambda in LMS and RGB space
    fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
    mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
    M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
    M = inv(M');
    LMSweighting = fundamentals'*Vlambda;
    LMSweighting = LMSweighting/norm(LMSweighting);
    RGBVlambda = inv(M)*LMSweighting; RGBVlambda = RGBVlambda/(8*norm(RGBVlambda));
    
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
            %             [~,whichframe] = max(sum(STAs.^2,1)); % finding the frame with the maximum energy
            t_offset = stro.trial(end,latencyidx);
            whichframe = round(t_offset/msperframe);
            STStmp = STS(:,whichframe)*STS(:,whichframe)';
            STCrosstmp = reshape(STCross(:,whichframe),size(STStmp));
            STCs = (nspikes*STCrosstmp-STStmp)/(nspikes*(nspikes-1));
            STAsub1{ii} = STAs(1:2:6,whichframe);
            STAsub2{ii} = STAs(2:2:6,whichframe);
            STCsub1{ii} = [STCs(1,1) STCs(1,3) STCs(1,5); STCs(3,1) STCs(3,3) STCs(3,5); STCs(5,1) STCs(5,3) STCs(5,5)]/nspikes; % standard error variance
            STCsub2{ii} = [STCs(2,2) STCs(2,4) STCs(2,6); STCs(4,2) STCs(4,4) STCs(4,6); STCs(6,2) STCs(6,4) STCs(6,6)]/nspikes;
            newmask = repmat(org_mask,[3 1]);
            pixelsinsubunit1 = find(newmask==1);
            pixelsinsubunit2 = find(newmask==2);
            newidxstoSTCOVmex = [pixelsinsubunit1; pixelsinsubunit2];
            
        else
            
            % WhiteNoise checkerboard as a stimuli
            mf1 = norm(STAsub1{ii})/norm(STAcheck1{ii});% multiplication factor 1
            mf2 = norm(STAsub1{ii})/norm(STAcheck1{ii}); % multiplication factor 2
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
            STCcheck1{ii} = A1*STC1*A1'/nspikes;
            STC2 = STCs(numel(pixelsinsubunit1)+1:end,numel(pixelsinsubunit1)+1:end);
            %             A2 = repmat(eye(3),[1 numel(pixelsinsubunit2)/3])./(numel(pixelsinsubunit2)/3);
            A2 = zeros(3,numel(pixelsinsubunit2)); A2(1,1:numel(pixelsinsubunit2)/3) = deal(3/numel(pixelsinsubunit2));
            A2(2,(numel(pixelsinsubunit2)/3)+1:2*numel(pixelsinsubunit2)/3) = deal(3/numel(pixelsinsubunit2));
            A2(3,(2*numel(pixelsinsubunit2)/3)+1:numel(pixelsinsubunit2)) = deal(3/numel(pixelsinsubunit2));
            STCcheck2{ii} = A2*STC2*A2'/nspikes;
        end
    end
    
    % Drawing an ellipsoid by calculating 95% confidence interval around the mean
    ptscheck1 = error_ellipse('C',STCcheck1{ii},'mu',STAcheck1{ii},'conf',confidence_interval);
    ptscheck2 = error_ellipse('C',STCcheck2{ii},'mu',STAcheck2{ii},'conf',confidence_interval);
    ptssub1 = error_ellipse('C',STCsub1{ii},'mu',STAsub1{ii},'conf',confidence_interval);
    ptssub2 = error_ellipse('C',STCsub2{ii},'mu',STAsub2{ii},'conf',confidence_interval);
    ptsLMSsub1 = M*ptssub1';
    ptsLMSsub2 = M*ptssub2';
    STAsub1LMS = M*STAsub1{ii};
    STAsub2LMS = M*STAsub2{ii};
    
    % Using the whitened matrix to figure out if the luminance function is within the 95% confidence interval
    [distance1(ii), rad1(ii), WwSTAsub1dist(ii)] = calc_whetherptisinside_confinterval(STAsub1{ii},ptssub1,RGBVlambda,STCsub1{ii},nspikes,alpha);
    [distance2(ii), rad2(ii), WwSTAsub2dist(ii)] = calc_whetherptisinside_confinterval(STAsub2{ii},ptssub2,RGBVlambda,STCsub2{ii},nspikes,alpha);
    [distSTA12LMS(ii), rad1LMS(ii), ~] = calc_whetherptisinside_confinterval(STAsub1LMS,ptsLMSsub1',STAsub2LMS,M*STCsub1{ii}*M',nspikes,alpha);
    [distSTA21LMS(ii), rad2LMS(ii), ~] = calc_whetherptisinside_confinterval(STAsub2LMS,ptsLMSsub2',STAsub1LMS,M*STCsub2{ii}*M',nspikes,alpha);

    [azimuthcheck1{ii},elevationcheck1{ii},~] = cart2sph(ptscheck1(:,1),ptscheck1(:,2),ptscheck1(:,3));
    [azimuthcheck2{ii},elevationcheck2{ii},~] = cart2sph(ptscheck2(:,1),ptscheck2(:,2),ptscheck2(:,3));
    [azimuthsub1{ii},elevationsub1{ii},~] = cart2sph(ptssub1(:,1),ptssub1(:,2),ptssub1(:,3));
    [azimuthsub2{ii},elevationsub2{ii},~] = cart2sph(ptssub2(:,1),ptssub2(:,2),ptssub2(:,3));
    
    %method 1
    [histcheck1,~] = hist3([azimuthcheck1{ii}*180/pi elevationcheck1{ii}*180/pi],{azimuthbins', elevationbins'});
    [histsub1,~] = hist3([azimuthsub1{ii}*180/pi elevationsub1{ii}*180/pi],{azimuthbins', elevationbins'});
    [histcheck2,~] = hist3([azimuthcheck2{ii}*180/pi elevationcheck2{ii}*180/pi],{azimuthbins', elevationbins'});
    [histsub2,~] = hist3([azimuthsub2{ii}*180/pi elevationsub2{ii}*180/pi],{azimuthbins', elevationbins'});
    [r,p] = corr([[histcheck1(:);histcheck2(:)] [histsub1(:);histsub1(:)]],'type','Pearson');
    rvalchecksub(ii) = r(1,2);
    pvalchecksub(ii) = p(1,2);
    
    %method 2
    [count1,~,~,~] = histcn([azimuthcheck1{ii} elevationcheck1{ii} azimuthcheck2{ii} elevationcheck2{ii}]*(180/pi),azimuthbins,elevationbins,azimuthbins,elevationbins);
    [count2,~,~,~] = histcn([azimuthsub1{ii} elevationsub1{ii} azimuthsub2{ii} elevationsub2{ii}]*(180/pi),azimuthbins,elevationbins,azimuthbins,elevationbins);
    [r,p] = corr([count1(:) count2(:)],'type','Pearson');
    rvalchecksub1(ii) = r(1,2);
    pvalchecksub1(ii) = max(p(1,2),eps);
    
    % Obtaining the basis vectors / STA
    neurothreshmode = stro.trial(:,neurothreshidx);
    basisvec_dropidx = inds(end);
    neurothresh_startidx = max(find(neurothreshmode==1,1),basisvec_dropidx+1); % There is possibility of pressing the Neurothresh button multiple times
    num_targetspikerates = unique(stro.trial(neurothresh_startidx:end,targetspikerateidx));
    t_offset = stro.trial(end,latencyidx)/1000;
    vect = stro.ras{basisvec_dropidx,basisvecidx};
    basisvec_size = nstixperside*nstixperside*3;
    numvect = (numel(vect)/basisvec_size)-1;
    basisvec = cell(1,numvect);
    for jj = 1:numvect
        tmp_vec = vect((jj-1)*basisvec_size+1:basisvec_size*jj) + vect(numvect*basisvec_size+1:basisvec_size*(numvect+1));
        basisvec{jj} = reshape(tmp_vec,[nstixperside nstixperside 3]);
    end
    bkgnd_monitor = reshape(vect(numvect*basisvec_size+1:basisvec_size*(numvect+1)),[nstixperside nstixperside 3]);
    basisvec1 = basisvec{1}-bkgnd_monitor;
    basisvec2 = basisvec{2}-bkgnd_monitor;
    tmp = basisvec1 + basisvec2;
    im = (0.5 * tmp./(max(abs(tmp(:))) + 0.01)) + 0.5;
        
    clear count1 count2
    % Converting the RGB values to LMS values and then checking for where
    % if the difference of signals between the subunit is cone-opponent
    ptssubcombinergb = error_ellipse('C',STCsub1{ii}+STCsub2{ii},'mu',STAsub1{ii}-STAsub2{ii},'conf',confidence_interval);
    ptssubcombinelms = M*ptssubcombinergb';
    whichoctants{ii} = findwhichoctant(ptssubcombinelms');
    if any(whichoctants{ii}==0) || any(whichoctants{ii}==7)
        doesitbelongtoluminanceoctant(ii) = 1;
    end
    
    % Finding the appropriate rotation matrix such that the cone weights of
    % the first subunit are alinged to the diagnonal in the first quadrant
    diagonalvector = norm(STAsub1LMS)*[1/sqrt(3); 1/sqrt(3); 1/sqrt(3)];
    R=fcn_RotationFromTwoVectors(STAsub1LMS,diagonalvector);
    STAsub1LMS_rotated = R*STAsub1LMS; STAsub1LMS_rotated = STAsub1LMS_rotated/norm(STAsub1LMS_rotated);
    STAsub2LMS_rotated = R*STAsub2LMS; STAsub2LMS_rotated = STAsub2LMS_rotated/norm(STAsub2LMS_rotated);
    
    STAsub1LMS_R{ii} = STAsub1LMS_rotated;
    STAsub2LMS_R{ii} = STAsub2LMS_rotated;
    
    if plotmode
        % plotting for subunits RGB
        figure(plot_counter); subplot(numsubplots,numsubplots,ii); plot3(ptssub1(:,1),ptssub1(:,2),ptssub1(:,3),'o','MarkerSize',2,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]); hold on;
        plot3(ptssub2(:,1),ptssub2(:,2),ptssub2(:,3),'o','MarkerSize',2,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]); hold on;
        plot3([RGBVlambda(1) -1*RGBVlambda(1)],[RGBVlambda(2) -1*RGBVlambda(2)],[RGBVlambda(3) -1*RGBVlambda(3)],'k','Linewidth',2);
        line([-0.03 0.03],[0 0],[0 0],'Linewidth',2); line([0 0],[-0.03 0.03],[0 0],'Linewidth',2); line([0 0],[0 0],[-0.03 0.03],'Linewidth',2);
        set(gca,'Xlim',[-0.03 0.03],'Ylim',[-0.03 0.03],'Zlim',[-0.03 0.03]); % title(char(filename(idx(ii))));
        axis equal; axis xy; grid on; xlabel('R'),ylabel('G'),zlabel('B'); hold off;
        
        % plotting for subunit 1,2:LMS
        figure(plot_counter+1);subplot(numsubplots,numsubplots,ii); plot3(ptsLMSsub1(1,:),ptsLMSsub1(2,:),ptsLMSsub1(3,:),'o','MarkerSize',2,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]); hold on;
        plot3(ptsLMSsub2(1,:),ptsLMSsub2(2,:),ptsLMSsub2(3,:),'o','MarkerSize',2,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]); hold on;
        plot3([LMSweighting(1) -1*LMSweighting(1)],[LMSweighting(2) -1*LMSweighting(2)],[LMSweighting(3) -1*LMSweighting(3)],'k','Linewidth',2);
        line([-0.5 0.5],[0 0],[0 0],'Linewidth',2); line([0 0],[-0.5 0.5],[0 0],'Linewidth',2); line([0 0],[0 0],[-0.5 0.5],'Linewidth',2);
        axis equal; axis xy; grid on; xlabel('L'),ylabel('M'),zlabel('S'); hold off;
            
        figure(plot_counter+2);  subplot(numsubplots,numsubplots,ii); image(im); set(gca,'XTick',[],'YTick',[]); 
        lum = rad1(ii)>distance1(ii) & rad2(ii)>distance2(ii) & significantsubunits(ii) == 2 & WwSTAsub2dist(ii) > rad2(ii) & WwSTAsub1dist(ii) > rad1(ii);
        if lum
            % Condition for Luminance cell
            set(gca,'XColor',[0 1 0],'YColor',[0 1 0]);
            luminancecells(ii) = 1;
        elseif ~lum & findwhichoctant(STAsub2LMS_rotated')==0
            % Condition for DO cell 
            set(gca,'XColor',[1 0 0],'YColor',[1 0 0]);
            DOcells(ii) = 1;
        elseif ~lum & findwhichoctant(STAsub2LMS_rotated')==7 % Same octant as the first subunit
            % Condition for SO cell
            SOcells(ii) = 1;
        else 
            hardtoclassifycells(ii) = 1;
        end
        axis square; hold off;
        
        % Rotating the Subunits in the LMS space
        figure(plot_counter+3);subplot(numsubplots,numsubplots,ii); plot3([0; STAsub1LMS_rotated(1,:)],[0;STAsub1LMS_rotated(2,:)],[0; STAsub1LMS_rotated(3,:)],'r','LineWidth',3,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]); hold on;
        plot3([0; STAsub2LMS_rotated(1,:)],[0;STAsub2LMS_rotated(2,:)],[0;STAsub2LMS_rotated(3,:)],'g','LineWidth',3,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]); hold on;
        line([-0.5 0.5],[0 0],[0 0],'Linewidth',2); line([0 0],[-0.5 0.5],[0 0],'Linewidth',2); line([0 0],[0 0],[-0.5 0.5],'Linewidth',2);
        axis equal; axis xy; grid on; xlabel('RL'),ylabel('RM'),zlabel('RS'); hold off;
        
        % Visualizing the angles between the 2 vectors for Lum, DO and hardtoclassify cells
        figure(plot_counter+4);
        if lum
            plot3([0; STAsub1LMS_rotated(1,:)],[0;STAsub1LMS_rotated(2,:)],[0; STAsub1LMS_rotated(3,:)],'Color',[0 0 0],'LineWidth',3); hold on;
            plot3([0; STAsub2LMS_rotated(1,:)],[0;STAsub2LMS_rotated(2,:)],[0;STAsub2LMS_rotated(3,:)],'Color',[0 0 0],'LineWidth',3); hold on;
        elseif ~lum & findwhichoctant(STAsub2LMS_rotated')==0
            % Condition for DO cell
            plot3([0; STAsub1LMS_rotated(1,:)],[0;STAsub1LMS_rotated(2,:)],[0; STAsub1LMS_rotated(3,:)],'Color',[1 0 0],'LineWidth',3); hold on;
            plot3([0; STAsub2LMS_rotated(1,:)],[0;STAsub2LMS_rotated(2,:)],[0;STAsub2LMS_rotated(3,:)],'Color',[1 0 0],'LineWidth',3); hold on;
       
        elseif ~lum & findwhichoctant(STAsub2LMS_rotated')==7 % Same octant as the first subunit
            % Condition for SO cell
            plot3([0; STAsub1LMS_rotated(1,:)],[0;STAsub1LMS_rotated(2,:)],[0; STAsub1LMS_rotated(3,:)],'Color',[0 1 0],'LineWidth',3); hold on;
            plot3([0; STAsub2LMS_rotated(1,:)],[0;STAsub2LMS_rotated(2,:)],[0;STAsub2LMS_rotated(3,:)],'Color',[0 1 0],'LineWidth',3); hold on; 
        else
            plot3([0; STAsub1LMS_rotated(1,:)],[0;STAsub1LMS_rotated(2,:)],[0; STAsub1LMS_rotated(3,:)],'Color',[0 1 1],'LineWidth',3); hold on;
            plot3([0; STAsub2LMS_rotated(1,:)],[0;STAsub2LMS_rotated(2,:)],[0;STAsub2LMS_rotated(3,:)],'Color',[0 1 1],'LineWidth',3); hold on;  
        end
        
        
        
    end
end
if plotmode
    figure(plot_counter); set(gcf,'Name','Subunit 1 2:RGB');
    figure(plot_counter+1); set(gcf,'Name','Subunit 1 2:LMS');
    figure(plot_counter+2); set(gcf,'Name','STA');
    figure(plot_counter+3); set(gcf,'Name','Rotated Subunit 1 2:LMS');
    figure(plot_counter+4); set(gcf,'Name','Rotated subunits of different cells'); hold on; line([-1 1],[0 0],[0 0],'Linewidth',3); line([0 0],[-1 1],[0 0],'Linewidth',3); line([0 0],[0 0],[-1 1],'Linewidth',3);
    axis equal; axis xy; grid on; xlabel('RL'),ylabel('RM'),zlabel('RS'); hold off;
    plot_counter = plot_counter + 4;
end

% Storing variables
newLUMidx = find(luminancecells);
newDOidx = find(DOcells);
newhardtoclassifyidx = find(hardtoclassifycells);
newSOidx = find(SOcells);

% Saving the variables
savevariables = 1;
if savevariables == 1
    save STAsub1 STAsub1
    save STAsub2 STAsub2
    save STCsub1 STCsub1
    save STCsub2 STCsub2
    save STAcheck1 STAcheck1
    save STAcheck2 STAcheck2
    save STCcheck1 STCcheck1
    save STCcheck2 STCcheck2
    save totspikes totspikes
    save azimuthcheck1 azimuthcheck1
    save azimuthcheck2 azimuthcheck2
    save azimuthsub1 azimuthsub1
    save azimuthsub2 azimuthsub2
    save elevationcheck1 elevationcheck1
    save elevationcheck2 elevationcheck2
    save elevationsub1 elevationsub1
    save elevationsub2 elevationsub2
    save rvalchecksub rvalchecksub
    save pvalchecksub pvalchecksub
    save rvalchecksub1 rvalchecksub1
    save pvalchecksub1 pvalchecksub1
    save whichoctants whichoctants
    save doesitbelongtoluminanceoctant doesitbelongtoluminanceoctant
    save newDOidx newDOidx
    save newLUMidx newLUMidx
    save newSOidx newSOidx
    save newhardtoclassifyidx newhardtoclassifyidx 
end
%% Using pretty corr to see any correlation amog variables
N = numel(filename_c);
load logresidualratio.mat
load anglebwvectorsRGB.mat
load absS.mat
load vals.mat
load z_scores.mat
load rvalchecksub.mat
load pvalchecksub.mat
load rvalchecksub1.mat
load pvalchecksub1.mat
% for Color cells
prettycorr([vals(1:N),z_scores(1:N),logresidualratio(1:N),anglebwvectorsRGB(1:N),absS(1:N)',rvalchecksub(1:N),log10(pvalchecksub(1:N)),rvalchecksub1(1:N),log10(pvalchecksub1(1:N))],{'eig1prctile','zscores','LogResRat','AngleRGB','absS','rvalchecksub','logpvalchecksub','rvalchecksub1','logpvalchecksub1'});
set(gcf,'Name','Color Cells');
% for Luminance cells
prettycorr([vals(N+1:end),z_scores(N+1:end),logresidualratio(N+1:end),anglebwvectorsRGB(N+1:end),absS(N+1:end)',rvalchecksub(N+1:end),log10(pvalchecksub(N+1:end)),rvalchecksub1(N+1:end),log10(pvalchecksub1(N+1:end))],{'eig1prctile','zscores','LogResRat','AngleRGB','absS','rvalchecksub','logpvalchecksub','rvalchecksub1','logpvalchecksub1'});
set(gcf,'Name','Luminance Cells');

