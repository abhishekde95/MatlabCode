% This exercise is based on suggestions of the committee to come up with a well-principled way of selecting the subunits
% Author - Abhishek De, 6/18, 
% modified 3/20, I picked the script from PopAnalysis_Crossvalidation.m in MatlabCode1/Abhishek/Spatialstructureanalysis_GregWN

close all; clearvars;
load fundamentals.mat 
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
resize_fact = 1;
% Include Gun noise data for statistical tests in order to estimate the RF size
channels = 3;
NPOINTS = 65536;
resize_fact2 = 1;
maxT = 9;

% Loading the Neurothresh files
% Loading all the files 
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM WNthresh');
NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
close(conn);
filename = filename(strcmp(string(NTmode),"subunit"));
NTmode = NTmode(strcmp(string(NTmode),"subunit"));
spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));


% Merging all the files in the list
Input_List = filename;
spikeIdx = spikeidx_NT;

numcells = numel(Input_List);
SSE = cell(numcells,3); 
Deviation = cell(numcells,3);
Peraccuracy = cell(numcells,3); 
Winning_model = zeros(numcells,1);

npartitions = 5; % 5-Fold cross validation
for ii = 1:numcells
    disp(ii);

    WN = nex2stro(findfile(char(Input_List(ii))));
    framerate = WN.sum.exptParams.framerate;
    nstixperside = WN.sum.exptParams.nstixperside;
    ntrials = length(WN.sum.absTrialNum);
    stimonidx = find(strcmp(WN.sum.trialFields(1,:),'stim_on'));
    stimoffidx = find(strcmp(WN.sum.trialFields(1,:),'all_off'));
    nframesidx = find(strcmp(WN.sum.trialFields(1,:),'num_frames'));
    noisetypeidx = find(strcmp(WN.sum.trialFields(1,:),'noise_type'));
    sigmaidxs = strmatch('sigma',WN.sum.trialFields(1,:));
    hepidx = find(strcmp(WN.sum.rasterCells(1,:),'AD11'));
    vepidx = find(strcmp(WN.sum.rasterCells(1,:),'AD12'));
    maskidx = strcmp(WN.sum.rasterCells(1,:), 'subunit_mask');
    anlgStartTimeidx = find(strcmp(WN.sum.rasterCells(1,:),'anlgStartTime'));
    L = WN.trial(:,noisetypeidx)==1;
    mu1idx = find(strcmp(WN.sum.trialFields(1,:),'mu1'));
    mu2idx = find(strcmp(WN.sum.trialFields(1,:),'mu2'));
    mu3idx = find(strcmp(WN.sum.trialFields(1,:),'mu3'));
    sigma1idx = find(strcmp(WN.sum.trialFields(1,:),'sigma1'));
    sigma2idx = find(strcmp(WN.sum.trialFields(1,:),'sigma2'));
    sigma3idx = find(strcmp(WN.sum.trialFields(1,:),'sigma3'));
    maskidx = strcmp(WN.sum.rasterCells(1,:),'subunit_mask');
    basisvecidx = strcmp(WN.sum.rasterCells(1,:),'basis_vec');
    gammaTable = WN.sum.exptParams.gamma_table;
    gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
    gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
    invgamma = InvertGamma(gammaTable, 0);
    sigmavect = unique(WN.trial(L,[sigma1idx sigma2idx sigma3idx]),'rows')/1000;
    sigmavect(all(any(sigmavect == 0),2),:) = [];
    gausslims = [WN.sum.exptParams.gauss_locut WN.sum.exptParams.gauss_hicut]/1000;
    x = linspace(gausslims(1),gausslims(2),NPOINTS);
    Fx = norminv(x)*sigmavect(1);
    sigmacorrectionfactor = std(Fx)./sigmavect(1);
    muvar = (sigmavect(1)*sigmacorrectionfactor)^2;
    
    % Getting the background rgb/lms
    ridx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_r'));
    gidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_g'));
    bidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_b'));
    bkgndRGB = [mode(WN.trial(:,ridx)), mode(WN.trial(:,gidx)), mode(WN.trial(:,bidx))];
    bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
    bkgndlms = M*bkgndrgb;
   
    WN.ras(~L ,:) = []; % modiftying the WN structure
    WN.trial(~L,:) = []; % modiftying the WN structure
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
        mask_changes  = mask_changes(:,1);  
        
        idxs = zeros(size(WN.trial,1),1);
        idxs(mask_changes(2,1)+1:end) = 1;
        idxs = logical(idxs);
        WN.ras(idxs,:) = []; % modiftying the WN structure 
        WN.trial(idxs,:) = []; % modiftying the WN structure
    end    
    spikename = 'sig001a'; 
    
    % Create partitions for cross validations
    c = cvpartition(mask_changes(1):mask_changes(2),'KFold',npartitions);
    SSE_Crescent = []; SSE_Gabor = []; SSE_DOG = [];
    Peraccuracy_Crescent = []; Peraccuracy_Gabor = []; Peraccuracy_DOG = [];
    Deviation_Crescent = []; Deviation_Gabor = []; Deviation_DOG = [];
    for jj = 1:c.NumTestSets
        % Calculating STA and STC for frames which triggered spikes
        stro1 = WN; stro1.trial(~c.training(jj),:) = []; stro1.ras(~c.training(jj),:) = []; % training data set
        stro2 = WN; stro2.trial(~c.test(jj),:) = []; stro2.ras(~c.test(jj),:) = []; % testing data set
        
        out_train = getWhtnsStats(stro1,maxT,'STAmex',{nstixperside^2,3,maxT},spikename);
        STS_train = out_train{1}; nspikes_train = out_train{3}; clear out_train;
        STAs_train = STS_train/nspikes_train;
        
        % Code for Statistical testing begins here
        s_train = std(STAs_train(:,1));
        STAs_train_z = STAs_train./s_train;
        grandz = zeros([nstixperside nstixperside]);
        maxzs = [];
        for i = 1:maxT
            tmp_train = reshape(STAs_train_z(:,i),[nstixperside nstixperside 3]); % This is the only place in the code where I use data from gun noise
            grandz = grandz+sum(tmp_train.^2,3);
            maxzs = [maxzs; sum(sum(tmp_train(:,:,1).^2)) sum(sum(tmp_train(:,:,2).^2)) sum(sum(tmp_train(:,:,3).^2))];
        end
        peakframe1 = max(maxzs(:,1)) == maxzs(:,1);
        peakframe2 = max(maxzs(:,2)) == maxzs(:,2);
        peakframe3 = max(maxzs(:,3)) == maxzs(:,3);
        peakframe = max(sum(maxzs,2)) == sum(maxzs,2);
        id = find(peakframe==1);
        latency = find(peakframe)*1000/stro1.sum.exptParams.framerate;
        if id~=1
            peakframe(id-1)= 1;
        end
        if id <=maxT-1
            peakframe(id+1)=1;
        end
        STAweights = sqrt(sum(STAs_train(:,peakframe).^2));
        STAweights = STAweights./sum(STAweights);
        tmpSTA = STAs_train(:,peakframe)*STAweights';
        tmpSTA = reshape(tmpSTA, [nstixperside^2 3]);
        [u1,~,v1] = svd(tmpSTA');
        m = reshape(v1(:,1),[nstixperside nstixperside]);
        
        [out2,fittedGabor,Rsquare_Gabor] = gaborfit_AD(m); % fitting Gabor
        [out3,fittedDOG,Rsquare_DOG] = DOGfit(m);
        
        [tmp_out1,fittedCrescent1,Rsquare_Crescent1] = Crescentfit_AD(m,out3); % Fitting the crescent shaped model, takes input from the DOG model
        [tmp_out2,fittedCrescent2,Rsquare_Crescent2] = Crescentfit_AD(m,[]);
        if tmp_out1.fval < tmp_out2.fval
            fittedCrescent = fittedCrescent1;
            out1 = tmp_out1;
            Rsquare_Crescent = Rsquare_Crescent1;
        else
            fittedCrescent = fittedCrescent2;
            out1 = tmp_out2;
            Rsquare_Crescent = Rsquare_Crescent2;
        end
        
        % Trying to project the spatiotempotal STA (derived from training
        % set) onto the test set % Based on Rich Pang's suggestion
        temporal_profile = sqrt(sum(STAs_train.^2));
        tC = fittedCrescent(:)*u1(:,1)'; tC = flipdim(tC(:)*temporal_profile,2); initargsC = {tC(:), 0, sum(stro2.trial(:,nframesidx)), [nstixperside^2 3 maxT]};
        tG = fittedGabor(:)*u1(:,1)'; tG = flipdim(tG(:)*temporal_profile,2); initargsG = {tG(:), 0, sum(stro2.trial(:,nframesidx)), [nstixperside^2 3 maxT]};
        tD = fittedDOG(:)*u1(:,1)'; tD = flipdim(tD(:)*temporal_profile,2); initargsD = {tD(:), 0, sum(stro2.trial(:,nframesidx)), [nstixperside^2 3 maxT]};
        outproj_C = getWhtnsStats(stro2,maxT,'STPROJmod',initargsC,spikename); % Crescent
        projC = outproj_C{1}; spikeC = outproj_C{2}; clear outproj_C; 
        rocC = roc(projC(spikeC==0), projC(spikeC==1));
        if rocC<0.5
            rocC = 1-rocC;
        end
        outproj_G = getWhtnsStats(stro2,maxT,'STPROJmod',initargsG,spikename); % Gabor
        projG = outproj_G{1}; spikeG = outproj_G{2}; clear outproj_G;
        rocG = roc(projG(spikeG==0), projG(spikeG==1));
        if rocG<0.5
            rocG = 1-rocG;
        end
        outproj_D = getWhtnsStats(stro2,maxT,'STPROJmod',initargsD,spikename); % DoG
        projD = outproj_D{1}; spikeD = outproj_D{2}; clear outproj_D;
        rocD = roc(projD(spikeD==0), projD(spikeD==1));
        if rocD<0.5
            rocD = 1-rocD;
        end
        Peraccuracy_Crescent = [Peraccuracy_Crescent; rocC]; 
        Peraccuracy_Gabor = [Peraccuracy_Gabor; rocG];
        Peraccuracy_DOG = [Peraccuracy_DOG; rocD];
        
        % Now cross validating 
        out_test = getWhtnsStats(stro2,maxT,'STAmex',{nstixperside^2,3,maxT},spikename);
        STS_test = out_test{1}; nspikes_test = out_test{3}; clear out_test;
        STAs_test = STS_test/nspikes_test;
        
        % Code for Statistical testing begins here
        s_test = std(STAs_test(:,1));
        STAs_test_z = STAs_test./s_test;
        grandz = zeros([nstixperside nstixperside]);
        maxzs = [];
        for i = 1:maxT
            tmp_test = reshape(STAs_test_z(:,i),[nstixperside nstixperside 3]); % This is the only place in the code where I use data from gun noise
            grandz = grandz+sum(tmp_test.^2,3);
            maxzs = [maxzs; sum(sum(tmp_test(:,:,1).^2)) sum(sum(tmp_test(:,:,2).^2)) sum(sum(tmp_test(:,:,3).^2))];
        end
        peakframe1 = max(maxzs(:,1)) == maxzs(:,1);
        peakframe2 = max(maxzs(:,2)) == maxzs(:,2);
        peakframe3 = max(maxzs(:,3)) == maxzs(:,3);
        peakframe = max(sum(maxzs,2)) == sum(maxzs,2);
        id = find(peakframe==1);
        latency = find(peakframe)*1000/stro2.sum.exptParams.framerate;
        if id~=1
            peakframe(id-1)= 1;
        end
        if id <=maxT-1
            peakframe(id+1)=1;
        end
        STAweights = sqrt(sum(STAs_test(:,peakframe).^2));
        STAweights = STAweights./sum(STAweights);
        tmpSTA = STAs_test(:,peakframe)*STAweights';
        testSTA = reshape(tmpSTA,[nstixperside.^2 3]);
        projSTA = testSTA*u1(:,1);
        projSTA = reshape(projSTA,[nstixperside nstixperside]);
        projSTA = projSTA/norm(projSTA(:));
        
        % Storing the SSE(s) of the model fits
        fittedCrescent = fittedCrescent/norm(fittedCrescent(:));
        fittedGabor = fittedGabor/norm(fittedGabor(:));
        fittedDOG = fittedDOG/norm(fittedDOG(:));
        SSE_Crescent = [SSE_Crescent; min([sum((fittedCrescent(:)-projSTA(:)).^2) sum((fittedCrescent(:)+projSTA(:)).^2)])];
        SSE_Gabor = [SSE_Gabor; min([sum((fittedGabor(:)-projSTA(:)).^2) sum((fittedGabor(:)+projSTA(:)).^2)])];
        SSE_DOG = [SSE_DOG; min([sum((fittedDOG(:)-projSTA(:)).^2) sum((fittedDOG(:)+projSTA(:)).^2)])];      
        Deviation_Crescent = [Deviation_Crescent; acos(abs(corr(fittedCrescent(:),projSTA(:))))*180/pi];
        Deviation_Gabor = [Deviation_Gabor; acos(abs(corr(fittedGabor(:),projSTA(:))))*180/pi];
        Deviation_DOG = [Deviation_DOG; acos(abs(corr(fittedDOG(:),projSTA(:))))*180/pi];
        

    end
    SSE{ii,1} = SSE_Crescent; SSE{ii,2} = SSE_Gabor; SSE{ii,3} = SSE_DOG;
    Peraccuracy{ii,1} = Peraccuracy_Crescent; Peraccuracy{ii,2} = Peraccuracy_Gabor; Peraccuracy{ii,3} = Peraccuracy_DOG;
    Deviation{ii,1} = Deviation_Crescent; Deviation{ii,2} = Deviation_Gabor; Deviation{ii,3} = Deviation_DOG;
    [~,idx] = min([mean(SSE_Crescent) mean(SSE_Gabor) mean(SSE_DOG)]);
    Winning_model(ii) = idx;
end

saveresults = 0;

if saveresults
   save IsoresponseRF_SSE SSE
   save IsoresponseRF_Peraccuracy Peraccuracy
   save IsoresponseRF_Deviation Deviation
end

