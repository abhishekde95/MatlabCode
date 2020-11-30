% Computing the NLI for the subunit white noise (based on Greg's suggestion)
% Author - Abhishek De, 5/20

close all; clearvars;
plot_counter = 1;

% Calculating the M matrix
load fundamentals.mat 
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% Include Gun noise data for statistical tests in order to estimate the RF size
channels = 3;
NPOINTS = 65536;
resize_fact2 = 1;
CHI2CRIT = 0.95; % For flaging a stixel as significant (adding gun and cone noise z-scores, squared)
maxT = 15;
crit = chi2inv(CHI2CRIT,300); % 3 color channels
count = 1;

% Loading the Neurothresh files
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
spikename_options = ['sig001a'; 'sig001b'];

WNsubunitNLI = [];
FV_simple = [];
FV_complex = [];
% computing the NLI 
for ii = 1:numcells 
    disp(ii);
    WN = nex2stro(findfile(char(Input_List{ii})));
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
    mu1idx = find(strcmp(WN.sum.trialFields(1,:),'mu1'));
    mu2idx = find(strcmp(WN.sum.trialFields(1,:),'mu2'));
    mu3idx = find(strcmp(WN.sum.trialFields(1,:),'mu3'));
    sigma1idx = find(strcmp(WN.sum.trialFields(1,:),'sigma1'));
    sigma2idx = find(strcmp(WN.sum.trialFields(1,:),'sigma2'));
    sigma3idx = find(strcmp(WN.sum.trialFields(1,:),'sigma3'));
    maskidx = strcmp(WN.sum.rasterCells(1,:),'subunit_mask');
    basisvecidx = strcmp(WN.sum.rasterCells(1,:),'basis_vec');
    latencyidx = strcmp(WN.sum.trialFields(1,:),'latency');
    neurothreshidx = strcmp(WN.sum.trialFields(1,:),'neurothresh');
    
    L = WN.trial(:,noisetypeidx)==1;
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
    Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
    Mrgbtocc = inv(Mrgbtocc');
   
  
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
        
        % Obtaining the basis vector from the isoresponse data
        neurothreshmode = WN.trial(:,neurothreshidx);
        basisvec_dropidx = inds(end);
        vect = WN.ras{basisvec_dropidx,basisvecidx};
        basisvec_size = nstixperside*nstixperside*3;
        numvect = (numel(vect)/basisvec_size)-1;
        tmpbasisvec = sum(reshape(vect(1:basisvec_size*numvect),[basisvec_size numvect]),2);
        normfactor = 0.5/((max(abs(tmpbasisvec(:))))*1.05);
        tmpbasisvec = normfactor*(tmpbasisvec)+0.5;
        tmpbasisvec = reshape(tmpbasisvec,[nstixperside nstixperside 3]);
        basisvec{ii} = tmpbasisvec; 
        
        % Modifying the STRO/WN structure after I have pulled out the
        % basis vector
        idxs = zeros(size(WN.trial,1),1);
        idxs(mask_changes(2,1):mask_changes(2,2)) = 1;
        idxs = logical(idxs);
        WN.ras(~idxs,:) = []; % modiftying the WN structure 
        WN.trial(~idxs,:) = []; % modiftying the WN structure
    end
    spikename = 'sig001a';
    
    % Calculating STA and STC for frames which triggered spikes
    out_gun = getWhtnsStats_AD(WN,maxT,'STCOVmex',{2,3,maxT},2,spikename);
    STS_gun = out_gun{1}; STCross_gun = out_gun{2}; nspikes_gun = out_gun{3}; clear out_gun;
    STAs_gun = STS_gun/nspikes_gun;
            
    % Code for Statistical testing begins here 
    s_gun = std(STAs_gun(:,1));
    STAs_gun_z = STAs_gun./s_gun;
    
    % Finding peak frame
    maxzs = sum(STAs_gun_z.^2,1);
    peakframe = maxzs == max(maxzs);
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
    
    % Calculating STC or spike-trigerred ensemble
    tmp = STS_gun(:,id)*STS_gun(:,id)';
    STCs = (nspikes_gun.*reshape(STCross_gun(:,id),[6 6])-tmp)/(nspikes_gun*(nspikes_gun-1));
    
    % Raw ensemble stats 
    out_all = getWhtnsStats_AD(WN,0,'STCOVmex',{2,3,1},2,spikename,[],1);
    STS_all = out_all{1}; STCross_all = out_all{2}; nspikes_all = out_all{3}; clear out_all;
    STAs_all = STS_all/nspikes_all;
    tmp_all = STS_all*STS_all';
    STCs_all = (nspikes_all.*reshape(STCross_all,[6 6])-tmp_all)/(nspikes_all*(nspikes_all-1));
    
    % Need to make some changes here 
    fact = diag(repmat([1;1],[3 1]));
    [b, ~, ~] = compiSTAC(fact*STAs_gun(:,id),fact*STCs*fact',fact*STAs_all,fact*STCs_all*fact', 1); % Running it through J Pillow's compute iSTAC code
    b = sign(dot(b,tmpSTA))*b;
    filts = b; whichframe = id;
    
    initargs = {filts, whichframe, sum(WN.trial(:,nframesidx)), [2 3 1]};     % projecting onto the maximally informative dimension
    out = getWhtnsStats_AD(WN,whichframe,'STPROJmod',initargs,2,spikename);
    proj = out{1}; Lspike = out{2}; clear out;
    lowerbound = prctile(proj,5);
    upperbound = prctile(proj,95);
    L = logical(proj < lowerbound | proj > upperbound);
    proj(L) = []; Lspike(L) = [];
    bins = linspace(min(proj), max(proj),15)';
    [Na,~] = hist(proj,bins);
    [Ns,~] = hist(proj(Lspike > 0),bins);
    fr = Ns./Na; % Feature function/ contrast-response function
    [B_lin, BINT_lin, R, RINT, STATS_lin] = regress(fr', [ones(length(bins),1), bins]);
    [B_quad, BINT_quad, R, RINT, STATS_quad] = regress(fr', [ones(length(bins),1), bins.^2]);
    [B_both, BINT_both, R, RINT, STATS_both] = regress(fr', [ones(length(bins),1), bins, bins.^2]);
    NLI = (STATS_quad(1)-STATS_lin(1))/STATS_both(1); % NLI 
    
    WNsubunitNLI = [WNsubunitNLI; NLI];
    
    fr = (fr-min(fr));
    fr = fr/max(abs(fr));
    
    if fr(1)>fr(end)
        fr = fliplr(fr);
    end
    
    if NLI<0
        FV_simple = [FV_simple; fr];
    else
        FV_complex = [FV_complex; fr];
    end
    
end

savevariables = 1;
if savevariables 
    save WNsubunitNLI WNsubunitNLI
end


%% Figure: Plotting the feature vectors of simple and complex cells based on the NLI criterion
load vals.mat

figure(plot_counter);
subplot(221); errorbar(mean(FV_simple,1),std(FV_simple,1)/sqrt(size(FV_simple,1)),'-ko'); axis square; xlabel('proj'); ylabel('Normalized firing rate'); 
set(gca,'Tickdir','out','Ylim',[0 1],'YTick',0:0.25:1.0); title('NLI<0'); text(3,0.75,strcat('N=',num2str(size(FV_simple,1))));
subplot(222); errorbar(mean(FV_complex,1),std(FV_complex,1)/sqrt(size(FV_complex,1)),'-ko'); axis square; xlabel('proj'); ylabel('Normalized firing rate'); 
set(gca,'Tickdir','out','Ylim',[0 1],'YTick',0:0.25:1.0); title('NLI>0'); text(3,0.75,strcat('N=',num2str(size(FV_complex,1))));
subplot(223); plot(WNsubunitNLI,vals,'o','MarkerSize',4,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); axis square; 
set(gca,'Tickdir','out','Xlim',[-1 1],'Ylim',[1 100],'YScale','log'); xlabel('NLI'); ylabel('PC1 confidence val'); 
plot_counter = plot_counter + 1;

