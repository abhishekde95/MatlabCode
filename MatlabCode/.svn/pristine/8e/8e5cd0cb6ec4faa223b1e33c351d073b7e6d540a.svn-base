% This is a variant of the code I wrote before and am trying it out for Greg's old WN data
% Comparing the Gabor and DOG fits for the WN data
% Author - Abhishek De, 1/19
close all; clearvars;
plot_counter = 1;
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
tmp_filename = fetch(conn,'SELECT filename FROM WNthresh');
NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
close(conn);
tmp_filename = tmp_filename(strcmp('subunit',NTmode));
spikeidx_NT = spikeidx_NT(strcmp('subunit',NTmode));
filename_NT = cell(size(tmp_filename));
for kk = 1:size(tmp_filename)
    filename_NT(kk) = {tmp_filename(kk)};
end

% Loading the just the subunit files 
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
tmp_filename = fetch(conn,'SELECT filename FROM WNSubunit');
subunit_mode = fetch(conn,'SELECT mode FROM WNSubunit');
spikeidx_subunit = cell2mat(fetch(conn,'SELECT spikeidx FROM WNSubunit'));
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
files_not_working_idxs = [];
count = 1;
CHI2CRIT = .95;
nrows = 10;
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
            files_not_working_idxs = [files_not_working_idxs; ii];
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
            files_not_working_idxs = [files_not_working_idxs; ii];
            flag = 1;
        end
    end
    if flag 
        continue;
    end
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
        
    spikeidx = spikeIdx(ii);
    spikename = spikename_options(spikeidx,:); 
    
    % Calculating STA and STC for frames which triggered spikes
    out_gun = getWhtnsStats(WN,maxT,'STCOVmex',{nstixperside^2,3,maxT},spikename);
    STS_gun = out_gun{1}; STCross_gun = out_gun{2}; nspikes_gun = out_gun{3}; clear out_gun;
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
    normfactor = 0.5/(max(abs(tmpSTA(:)))+0.01);
    im = normfactor*tmpSTA + 0.5;
    im = reshape(im,[nstixperside nstixperside channels]);
    
    % Trying out higher-order singular value decomposition
    X = tensor(reshape(STAs_gun,[nstixperside^2 channels maxT]));
    T = hosvd(X,0.01);
       
    Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
    Mrgbtocc = inv(Mrgbtocc');
    conewts_svd = Mrgbtocc * T.u{2}(:,1);
    conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
    conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);
    
    figure(plot_counter),subplot(nrows,4,4*count-3);imagesc(reshape(T.u{1}(:,1),[nstixperside nstixperside])); colormap('gray'); axis square; set(gca,'XTick',[],'YTick',[]); 
    subplot(nrows,4,4*count-2);plot(T.u{3}(:,1),'k','Linewidth',2); axis square; set(gca,'Tickdir','out'); 
    subplot(nrows,4,4*count-1); bar(conewts_svd','FaceColor','k'); set(gca,'Ylim',[-1 1],'Xlim',[0 4],'XTick',[1 2 3],'XTicklabel',{'L','M','S'},'Tickdir','out'); axis square; 
    subplot(nrows,4,4*count); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
     
    count = count + 1;
    if count == nrows+1
        count = 1;
        plot_counter = plot_counter + 1;
    end
    
end


