% Population of doing waveform analyses
% Author - Abhishek De, 10/19

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
CHI2CRIT = 0.95; % For flaging a stixel as significant (adding gun and cone noise z-scores, squared)
maxT = 9;
crit = chi2inv(CHI2CRIT,300); % 3 color channels
spikename_options = ['sig001a'; 'sig001b'];
[filename_Lum, spikeIdx_Lum] = fnamesFromTxt2('Lum.txt');
[filename_ColorOpponent, spikeIdx_ColorOpponent] = fnamesFromTxt2('ColorOpponent.txt');

% Loading the Neurothresh files
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
tmp = fetch(conn,'SELECT filename FROM WNthresh'); tmp_filename = tmp.filename;
tmp = fetch(conn,'SELECT NTmode FROM WNthresh'); NTmode = tmp.NTmode;
tmp = fetch(conn,'SELECT spikeidx FROM WNthresh'); spikeidx_NT = tmp.spikeidx;
close(conn);
tmp_filename = tmp_filename(strcmp('subunit',NTmode));
spikeidx_NT = spikeidx_NT(strcmp('subunit',NTmode));
filename_NT = cell(size(tmp_filename));
for kk = 1:size(tmp_filename)
    filename_NT(kk) = {tmp_filename(kk)};
end

% Loading the just the subunit files 
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
tmp = fetch(conn,'SELECT filename FROM WNSubunit'); tmp_filename = tmp.filename;
tmp = fetch(conn,'SELECT mode FROM WNSubunit'); subunit_mode = tmp.mode;
tmp = fetch(conn,'SELECT spikeidx FROM WNSubunit'); spikeidx_subunit = tmp.spikeidx;
close(conn);
tmp_filename = tmp_filename(strcmp('STA',subunit_mode));
spikeidx_subunit = spikeidx_subunit(strcmp('STA',subunit_mode));
filename_subunit = cell(size(tmp_filename));
for kk = 1:size(tmp_filename)
    filename_subunit(kk) = {tmp_filename(kk)};
end

% Merging all the files in the list
Input_List = [filename_Lum; filename_ColorOpponent; filename_NT; filename_subunit];
spikeIdx = [spikeIdx_Lum; spikeIdx_ColorOpponent; spikeidx_NT; spikeidx_subunit];
plotresults = 0;
computeNLI = 1;
numcells = numel(Input_List);
files_not_working = [];
count = 1;
Singleopponent_waveform = [];
timedur = 0.2;
for ii = 1:numcells
    flag = 0;
    disp(ii);
    filename = char(Input_List{ii}{1}); % acquiring the filename (1st column) from the List
    
    WN = {};
    for jj = 1:size(Input_List{ii},2)
        try 
            tmpstro = nex2stro(findfile(char(Input_List{ii}(jj))));
        catch 
            files_not_working = [files_not_working; Input_List{ii}];
            flag = 1;
            break;
        end
        if (isempty(WN))
            WN = tmpstro;
        else
            WN = strocat(WN, tmpstro);
        end
        if ~any(strcmp(WN.sum.rasterCells,'sig001a_wf'))
            files_not_working = [files_not_working; Input_List{ii}];
            flag = 1;
        end
    end
    if flag 
        continue;
    end
    Output_List{count,1} = filename;
    framerate = WN.sum.exptParams.framerate;
    nstixperside = WN.sum.exptParams.nstixperside;
    ntrials = length(WN.sum.absTrialNum);
    fponidx = find(strcmp(WN.sum.trialFields(1,:),'fp_on'));
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
    spikeidx = spikeIdx(ii);
    spikename = spikename_options(spikeidx,:); 
    spikeidx = strcmp(WN.sum.rasterCells(1,:),spikename);
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
    % storing backgroung LMS excitation info
    Output_List{count,16} = bkgndlms;

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
    
    % Calculating baselinefiring rate from GetSacData
    baselineFR = [];
    try
        for jj = mask_changes(1):mask_changes(end) 
            stimontime = WN.trial(jj,stimonidx);
            numspikes = sum(WN.ras{jj,spikeidx}>stimontime-timedur & WN.ras{jj,spikeidx}<stimontime);
            baselineFR = [baselineFR; numspikes/timedur];
        end
    catch
        baselineFR = NaN;
    end
    Output_List{count,18} = mean(baselineFR); 
    
    % Calculating STA and STC for frames which triggered spikes
    out_gun = getWhtnsStats(WN,maxT,'STCOVmex',{nstixperside^2,3,maxT},spikename);
    STS_gun = out_gun{1}; STCross_gun = out_gun{2}; nspikes_gun = out_gun{3}; clear out_gun;
    Output_List{count,15} = nspikes_gun;
    STAs_gun = STS_gun/nspikes_gun;
            
    % Code for Statistical testing begins here 
    s_gun = std(STAs_gun(:,1));
    STAs_gun_z = STAs_gun./s_gun;
    
    % Spatial map
    grandz = zeros([nstixperside nstixperside]);
    maxzs = [];
    for i = 1:maxT
        tmp_gun = reshape(STAs_gun_z(:,i),[nstixperside nstixperside 3]); % This is the only place in the code where I use data from gun noise
        grandz = grandz+sum(tmp_gun.^2,3);
        maxzs = [maxzs; sum(sum(tmp_gun(:,:,1).^2)) sum(sum(tmp_gun(:,:,2).^2)) sum(sum(tmp_gun(:,:,3).^2))];
    end
    peakframe1 = max(maxzs(:,1)) == maxzs(:,1);
    peakframe2 = max(maxzs(:,2)) == maxzs(:,2);
    peakframe3 = max(maxzs(:,3)) == maxzs(:,3);
    peakframe = max(sum(maxzs,2)) == sum(maxzs,2);
    id = find(peakframe==1);
    latency = find(peakframe)*1000/WN.sum.exptParams.framerate;
    if id~=1
        peakframe(id-1)= 1;
    end
    if id <=maxT-1
        peakframe(id+1)=1;
    end
    STAweights = sqrt(sum(STAs_gun(:,peakframe).^2));
    STAweights = STAweights./sum(STAweights);
    tmpSTA = STAs_gun(:,peakframe)*STAweights';
    Output_List{count,2} = reshape(tmpSTA, [nstixperside^2 3]);
    
    
    STAgunmat_pf = Output_List{count,2};
    [u1,s,v1] = svd(STAgunmat_pf');
    m = reshape(v1(:,1),[nstixperside nstixperside]);
    Output_List{count,4} = sign(m).*(abs(m).^(1.0));
    Output_List{count,5} = u1(:,1); % storing R,G,B values
    eigvals = [s(1,1) s(2,2) s(3,3)];
    Output_List{count,17} = eigvals./sum(eigvals);
    
    % Checking if single-opponent or not
    ft = fft2(Output_List{count,4});
    if max(abs(ft(:)))== abs(ft(1))
        Singleopponent_waveform = [Singleopponent_waveform; 1];
    else
        Singleopponent_waveform = [Singleopponent_waveform; 0];
    end
   
    % storing latency
    Output_List{count,6} = latency;
    
    % Storing the zscore based on spatiotemporal energy
    zscoremeans = STAs_gun/std(STAs_gun(:,1));
    Output_List{count,7} = sum(zscoremeans.^2,1);
    
    % Storing the eccentricity
    Output_List{count,8} = [WN.sum.exptParams.rf_x WN.sum.exptParams.rf_y];
    
    % Storing other essential elements
    tmp = STS_gun(:,id)*STS_gun(:,id)';
    STCs = (nspikes_gun.*reshape(STCross_gun(:,id),[300 300])-tmp)/(nspikes_gun*(nspikes_gun-1));
    P = eye(size(STCs)) - STAs_gun(:,id)*inv(STAs_gun(:,id)'*STAs_gun(:,id))*STAs_gun(:,id)';
    subSTCs = P*STCs*P';
    [v,d] = eig(subSTCs);
    [~, idxs] = sort(diag(d));
    Output_List{count,9} = diag(d);
    Output_List{count,10} = real(reshape(v(:,idxs(end)),[nstixperside nstixperside 3])); % excitatory filter
    Output_List{count,11} = real(reshape(v(:,idxs(1)),[nstixperside nstixperside 3])); % suppressive filter
    
    % Storing variance in the noise, using for an analysis suggested by Ringach et al., 2002
    [un,~,vn] = svd(reshape(STAs_gun(:,1),[nstixperside^2 3])');
    noise_image = sign(m(:)'*vn(:,1))*reshape(vn(:,1),[nstixperside nstixperside]);
    Output_List{count,14} = noise_image;
    
    % Calculating NLI (Horwitz et al; 2007)
    % Calculating maximally informative dimension
    if computeNLI   
%         Calculating STA and STC for all frames
        out_all = getWhtnsStats(WN,0,'STCOVmex',{nstixperside^2,3,1},spikename,[],1);
        STS_all = out_all{1}; STCross_all = out_all{2}; nspikes_all = out_all{3}; clear out_all;
        STAs_all = STS_all/nspikes_all;
        tmp_all = STS_all*STS_all';
        STCs_all = (nspikes_all.*reshape(STCross_all,[300 300])-tmp_all)/(nspikes_all*(nspikes_all-1));
        Output_List{count,3} = getRFfromSTA(tmpSTA,1,0.75); % arbitrary threshold of 70%
        fact = Output_List{count,3};
        fact = diag(repmat(fact(:),[3 1]));
        [b, ~, ~] = compiSTAC(fact*STAs_gun(:,id),fact*STCs*fact',fact*STAs_all,fact*STCs_all*fact', 1); % Running it through J Pillow's compute iSTAC code
        b = sign(dot(b,tmpSTA))*b;
        filts = b; whichframe = id;
        
        initargs = {filts, whichframe, sum(WN.trial(:,nframesidx)), [nstixperside^2 3 1]};     % projecting onto the maximally informative dimension
        out = getWhtnsStats(WN,whichframe,'STPROJmod',initargs,spikename);
        proj = out{1}; Lspike = out{2}; clear out;
        lowerbound = prctile(proj,5);
        upperbound = prctile(proj,95);
        L = logical(proj < lowerbound | proj > upperbound);
        proj(L) = []; Lspike(L) = [];
        bins = linspace(min(proj), max(proj),15)';
        [Na,~] = hist(proj,bins);
        [Ns,~] = hist(proj(Lspike > 0),bins);
        fr = Ns./Na;
        [B_lin, BINT_lin, R, RINT, STATS_lin] = regress(fr', [ones(length(bins),1), bins]);
        [B_quad, BINT_quad, R, RINT, STATS_quad] = regress(fr', [ones(length(bins),1), bins.^2]);
        [B_both, BINT_both, R, RINT, STATS_both] = regress(fr', [ones(length(bins),1), bins, bins.^2]);
        NLI = (STATS_quad(1)-STATS_lin(1))/STATS_both(1);
        Output_List{count,12} = filts;
        Output_List{count,13} = NLI;
    end
    count = count + 1;
end

savevariables = 0;
if savevariables 
    save Output_List_waveform Output_List
    save Singleopponent_waveform Singleopponent_waveform
end


