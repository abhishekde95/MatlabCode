% Classifying cells based on SVD 
% Similar to what I implemented for the spatial-structure project 
% Author - Abhishek De, 
close all; clearvars;
plot_counter = 1;

% Loading the z scores
load z_scores.mat
load vals.mat

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
try
    % Directly acquiring the filenames from the SQL database
    conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
    filename = fetch(conn,'SELECT filename FROM WNthresh');
    NTmode = fetch(conn,'SELECT NTmode FROM WNthresh');
    spikeidx_NT = cell2mat(fetch(conn,'SELECT spikeidx FROM WNthresh'));
    close(conn);
    
    filename = filename(strcmp(string(NTmode),"subunit"));
    NTmode = NTmode(strcmp(string(NTmode),"subunit"));
    spikeidx_NT = spikeidx_NT(strcmp(string(NTmode),"subunit"));
catch 
    C = readtable('WNthresh.csv');
    NTmode = C.Var5;
    filename = C.Var2(strcmp(NTmode,"subunit"));
    spikeidx_NT = C.Var3(strcmp(NTmode,"subunit"));
    
end
    

% Merging all the files in the list
Input_List = filename;
spikeIdx = spikeidx_NT;
numcells = numel(Input_List);
spikename_options = ['sig001a'; 'sig001b'];
variance_accounted = [];
conewts_svd = [];
numsubplots = ceil(sqrt(numcells));
singular_values = [];
relativestrength_of_subunits = [];
basisvec = cell(1,numcells);
SVDSTA = cell(1,numcells);
angulardifference_RGB = [];
angulardifference_LMS = [];
S1RGB_svd = [];
S2RGB_svd = [];
S1LMS_svd = [];
S2LMS_svd = [];
RGBspatialweightingfunction_svd = [];

for ii = 1:numcells
    disp(ii);
       
    WN = nex2stro(findfile(char(Input_List{ii}),'/Users/abhishekde/Google Drive/Data_Physiology'));
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
    tmpSTA = reshape(STAs_gun(:,peakframe)*STAweights',[2 3]);
    [u1,s,v1] = svd(tmpSTA);
    RGB_svd = v1(:,1);
    relativestrength_of_subunits = [relativestrength_of_subunits u1(:,1)];
    
    % Template for STA
    template = all_masks{mask_changes(1,2)};
    R = template; G = template; B = template;
    R(R==1) = tmpSTA(1,1); R(R==2) = tmpSTA(2,1);
    G(G==1) = tmpSTA(1,2); G(G==2) = tmpSTA(2,2);
    B(B==1) = tmpSTA(1,3); B(B==2) = tmpSTA(2,3);
    tmpSVDSTA = [R; G; B];
    normfactor = 0.5/((max(abs(tmpSVDSTA(:))))*1.05);
    tmpSVDSTA = normfactor*(tmpSVDSTA)+0.5;
    tmpSVDSTA = reshape(tmpSVDSTA,[nstixperside nstixperside 3]);
    SVDSTA{ii} = tmpSVDSTA; 
    
    % Converting into LMS cone weights 
    tmp_conewts_svd = Mrgbtocc * RGB_svd;
    tmp_conewts_svd = tmp_conewts_svd./repmat(sum(abs(tmp_conewts_svd),1),[3 1]);
    tmp_conewts_svd = tmp_conewts_svd .* repmat(sign(tmp_conewts_svd(2,:)),[3 1]);
    conewts_svd = [conewts_svd tmp_conewts_svd];
    variance_accounted = [variance_accounted; s(1,1)^2/(s(1,1)^2+s(2,2)^2)];
    singular_values = [singular_values; s(1,1) s(2,2)];
    
    % Calculating the angular difference 
    RGB1 = tmpSTA(1,:); RGB2 = tmpSTA(2,:);
    LMS1 = Mrgbtocc*RGB1'; LMS2 = Mrgbtocc*RGB2'; 
    angulardifference_RGB = [angulardifference_RGB; acos(dot(RGB1,RGB2)/(norm(RGB1)*norm(RGB2)))*180/pi];
    angulardifference_LMS = [angulardifference_LMS; acos(dot(LMS1,LMS2)/(norm(LMS1)*norm(LMS2)))*180/pi];
    S1RGB_svd = [S1RGB_svd RGB1'];
    S2RGB_svd = [S2RGB_svd RGB2'];
    S1LMS_svd = [S1LMS_svd LMS1];
    S2LMS_svd = [S2LMS_svd LMS2];
    RGBspatialweightingfunction_svd = [RGBspatialweightingfunction_svd; u1(:,1)'];
    

end
plot_counter = plot_counter + 2;

% Saving the SVD derived cone weights
savevariables = 1;
if savevariables 
    save conewts_svd conewts_svd
    save variance_accounted variance_accounted
    save singular_values singular_values
    save angulardifferences_RGB angulardifference_RGB
    save angulardifferences_LMS angulardifference_LMS
    save S1RGB_svd S1RGB_svd 
    save S2RGB_svd S2RGB_svd
    save S1LMS_svd S1LMS_svd 
    save S2LMS_svd S2LMS_svd
    save RGBspatialweightingfunction_svd RGBspatialweightingfunction_svd
end

%% Classifying cells the same way as for spatial structure paper 
thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];

LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts Sconedominated_conewts];
hardtoclassifyidx = [Other_conewts];
hardtoclassifyidx = [hardtoclassifyidx LUMidx(vals(LUMidx)>=95) DOidx(vals(DOidx)>=95)];
LUMidx = LUMidx(vals(LUMidx)<95);
DOidx = DOidx(vals(DOidx)<95);

% Loading the old cell classification
load newLUMidx.mat
load newDOidx.mat
load newSOidx.mat
load newhardtoclassifyidx.mat
oldLUMidx = newLUMidx';
oldDOidx = newDOidx';
oldhardtoclassifyidx = [newSOidx' newhardtoclassifyidx'];

% Loading the LMS cone weights of each subunit 
load S1LMS.mat
load S2LMS.mat

% Checking the correlation with non-linearity indices 
% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat

% Load the integration within the subunit data
load AUROClinS1_CV.mat
load AUROCquadS1_CV.mat
load AUROClinS2_CV.mat
load AUROCquadS2_CV.mat

% Loading the baseline FR stats
load baselineFRstats.mat

% For storing median of differences/ratios
RSSEisoresp_medianofratios = [];
Withinsubunits_medianofdifferences = [];
r = []; p = [];

for ii = 1:numel(RSSE_linearmodel)   
    % computation for calculating median of differences/ratios
    RSSEisoresp_medianofratios = [RSSEisoresp_medianofratios; median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii})];
    Withinsubunits_medianofdifferences = [Withinsubunits_medianofdifferences; median([median(AUROCquad1{ii}-AUROClin1{ii}) median(AUROCquad2{ii}-AUROClin2{ii})])];
    
    
    % Computing the Spearman's r for the baseline FR
    [tmp_r, tmp_p] = corr((1:numel(baselineFRstats{ii}))',baselineFRstats{ii},'type','Pearson');
    r = [r; tmp_r];
    p = [p; tmp_p];
end


%*****************************************************%
% Plotting the cone weights
figure(plot_counter); set(gcf,'Name','Cone wts') 
subplot(221); plot(conewts_svd(1,:),conewts_svd(2,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');  
xlabel('L'),ylabel('M'); title('Cone weights - ALL cells'); 
subplot(222); plot(conewts_svd(1,hardtoclassifyidx),conewts_svd(2,hardtoclassifyidx),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,LUMidx),conewts_svd(2,LUMidx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,DOidx),conewts_svd(2,DOidx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');  
xlabel('L'), ylabel('M'); title('Cone weights - LUM, DO, HTC cells');
subplot(223); plot(conewts_svd(1,oldhardtoclassifyidx),conewts_svd(2,oldhardtoclassifyidx),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,oldLUMidx),conewts_svd(2,oldLUMidx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,oldDOidx),conewts_svd(2,oldDOidx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');  
xlabel('L'), ylabel('M'); title('Old cone wt classification - LUM, DO, HTC cells');
subplot(224); plot(S1LMS(1,oldhardtoclassifyidx),S1LMS(2,oldhardtoclassifyidx),'o','MarkerSize',3,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(S2LMS(1,oldhardtoclassifyidx),S2LMS(2,oldhardtoclassifyidx),'o','MarkerSize',3,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
plot(S1LMS(1,oldLUMidx),S1LMS(2,oldLUMidx),'o','MarkerSize',4,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
plot(S2LMS(1,oldLUMidx),S2LMS(2,oldLUMidx),'o','MarkerSize',4,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(S1LMS(1,oldDOidx),S1LMS(2,oldDOidx),'o','MarkerSize',4,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(S2LMS(1,oldDOidx),S2LMS(2,oldDOidx),'o','MarkerSize',4,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[-1 1],'XTick',-1:0.5:1,'YTick',-1:0.5:1,'Tickdir','out'); 
plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 -1],[0 0],'k'); plot([1 0],[0 -1],'k'); plot([0 -1],[-1 0],'k');  xlabel('L'), ylabel('M');
title('Subunit cone wts & OLD classification '); 
plot_counter = plot_counter + 1;


% Plotting the results for SVD based cone weight classification including the PC1 z-scores 
figure(plot_counter);
subplot(321); histogram(RSSEisoresp_medianofratios(hardtoclassifyidx),logspace(-1,3,31),'displaystyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(median(RSSEisoresp_medianofratios(hardtoclassifyidx)),20,'v','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); 
histogram(RSSEisoresp_medianofratios(LUMidx),logspace(-1,3,31),'displaystyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); 
plot(median(RSSEisoresp_medianofratios(LUMidx)),20,'v','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
histogram(RSSEisoresp_medianofratios(DOidx),logspace(-1,3,31),'displaystyle','stairs','EdgeColor',[1 0 0],'Linewidth',2);
plot(median(RSSEisoresp_medianofratios(DOidx)),17,'v','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 1000],'XTick',[0.1 1 10 100 1000]); ylabel('Count'); title('Isoresponse'); xlabel('median CV ratio of errors'); axis square; hold off;
subplot(322); histogram(100*(Withinsubunits_medianofdifferences(hardtoclassifyidx)),linspace(-2,8,31),'displaystyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(median(100*(Withinsubunits_medianofdifferences(hardtoclassifyidx))),20,'v','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
histogram(100*(Withinsubunits_medianofdifferences(LUMidx)),linspace(-2,8,31),'displaystyle','stairs','EdgeColor',[0 0 0],'Linewidth',2);
plot(median(100*(Withinsubunits_medianofdifferences(LUMidx))),20,'v','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
histogram(100*(Withinsubunits_medianofdifferences(DOidx)),linspace(-2,8,31),'displaystyle','stairs','EdgeColor',[1 0 0],'Linewidth',2);
plot(median(100*(Withinsubunits_medianofdifferences(DOidx))),17,'v','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-2 10]); xlabel('median CV GQM-GLM AUROC'); ylabel('Count'); title('Within subunits'); axis square; hold off;
subplot(323); plot(z_scores(hardtoclassifyidx),RSSEisoresp_medianofratios(hardtoclassifyidx),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(z_scores(LUMidx),RSSEisoresp_medianofratios(LUMidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(z_scores(DOidx),RSSEisoresp_medianofratios(DOidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out'); xlabel('Z scores'); ylabel('(Lin/Quad) errors');axis square; hold off;
subplot(324); plot(vals(hardtoclassifyidx),RSSEisoresp_medianofratios(hardtoclassifyidx),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(vals(LUMidx),RSSEisoresp_medianofratios(LUMidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(vals(DOidx),RSSEisoresp_medianofratios(DOidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out'); xlabel('Z-score Percentiles'); ylabel('(Lin/Quad) errors');axis square; hold off;
subplot(325); plot(abs(r(hardtoclassifyidx)),RSSEisoresp_medianofratios(hardtoclassifyidx),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(abs(r(LUMidx)),RSSEisoresp_medianofratios(LUMidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
plot(abs(r(DOidx)),RSSEisoresp_medianofratios(DOidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[0 0.6],'Ylim',[0.1 1000],'YScale','log'); xlabel('Spearman-r baseline FR'); ylabel('median CV ratio of errors');
subplot(326); plot(abs(r(hardtoclassifyidx)),100*(Withinsubunits_medianofdifferences(hardtoclassifyidx)),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(abs(r(LUMidx)),100*(Withinsubunits_medianofdifferences(LUMidx)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
plot(abs(r(DOidx)),100*(Withinsubunits_medianofdifferences(DOidx)),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[0 0.6],'Ylim',[-2 8]); xlabel('Spearman-r baseline FR'); ylabel('median CV GQM-GLM AUROC');
plot_counter = plot_counter + 1;

% Need to plot the STA/basis vectors of LUM, DO and HTC cells 
sequence_svd = [LUMidx DOidx hardtoclassifyidx];
sequence_subunitswts = [oldLUMidx oldDOidx oldhardtoclassifyidx];
figure(plot_counter); set(gcf,'Name','SVD');
figure(plot_counter+1); set(gcf,'Name','SVD STA'); 
figure(plot_counter+2); set(gcf,'Name','Subunit wts'); 
for ii = 1:numcells
    
    % Plotting the SVD based classified cells 
    vec1 = basisvec{sequence_svd(ii)};
    figure(plot_counter); subplot(numsubplots,numsubplots,ii); image(vec1); set(gca,'XTick',[],'YTick',[]); axis square; 
    if ismember(sequence_svd(ii),LUMidx)
        set(gca,'XColor','g','YColor','g');
    elseif ismember(sequence_svd(ii),DOidx)
        set(gca,'XColor','r','YColor','r');
    else
    end
    
    % Plotting the SVD STA 
    vec2 = SVDSTA{sequence_svd(ii)};
    figure(plot_counter+1); subplot(numsubplots,numsubplots,ii); image(vec2); set(gca,'XTick',[],'YTick',[]); axis square; 
    if ismember(sequence_svd(ii),LUMidx)
        set(gca,'XColor','g','YColor','g');
    elseif ismember(sequence_svd(ii),DOidx)
        set(gca,'XColor','r','YColor','r');
    else
    end
    
    % Plotting the subunit cone wt based classified cells
    vec3 = basisvec{sequence_subunitswts(ii)};
    figure(plot_counter+2); subplot(numsubplots,numsubplots,ii); image(vec3); set(gca,'XTick',[],'YTick',[]); axis square;
    if ismember(sequence_subunitswts(ii),newLUMidx)
        set(gca,'XColor','g','YColor','g');
    elseif ismember(sequence_subunitswts(ii),newDOidx)
        set(gca,'XColor','r','YColor','r');
    else
    end
    
end
plot_counter = plot_counter + 3;


%% Classifying cells using sone wts and NLI: making it parallel to the spatial structure analysis  
thresh = 0.8;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];
Sconesensitive = conewts_svd(:,Sconedominated_conewts);
Sconedominated_conewts(sign(Sconesensitive(1,:))==1 & sign(Sconesensitive(3,:))==1) = [];


load WNsubunitNLI.mat
LUMidx = LumIds_conewts;
DOidx = [ColorOpponentIds_conewts Sconedominated_conewts];
hardtoclassifyidx = [Other_conewts];
hardtoclassifyidx = [hardtoclassifyidx LUMidx(WNsubunitNLI(LUMidx)>=0) DOidx(WNsubunitNLI(DOidx)>=0)];
LUMidx = LUMidx(WNsubunitNLI(LUMidx)<0);
DOidx = DOidx(WNsubunitNLI(DOidx)<0);

% Loading the old cell classification
load newLUMidx.mat
load newDOidx.mat
load newSOidx.mat
load newhardtoclassifyidx.mat
oldLUMidx = newLUMidx';
oldDOidx = newDOidx';
oldhardtoclassifyidx = [newSOidx' newhardtoclassifyidx'];

% Loading the LMS cone weights of each subunit 
load S1LMS.mat
load S2LMS.mat

% Checking the correlation with non-linearity indices 
% Load the isoresponse data
load RSSE_linearmodel_CV.mat % Robust regression
load RSSE_quadmodel_CV.mat

% Load the integration within the subunit data
load AUROClinS1_CV.mat
load AUROCquadS1_CV.mat
load AUROClinS2_CV.mat
load AUROCquadS2_CV.mat

% Loading the baseline FR stats
load baselineFRstats.mat

% For storing median of differences/ratios
RSSEisoresp_medianofratios = [];
Withinsubunits_medianofdifferences = [];
r = []; p = [];

for ii = 1:numel(RSSE_linearmodel)   
    % computation for calculating median of differences/ratios
    RSSEisoresp_medianofratios = [RSSEisoresp_medianofratios; median(RSSE_linearmodel{ii}./RSSE_quadmodel{ii})];
    Withinsubunits_medianofdifferences = [Withinsubunits_medianofdifferences; median([median(AUROCquad1{ii}-AUROClin1{ii}) median(AUROCquad2{ii}-AUROClin2{ii})])];
    
    
    % Computing the Spearman's r for the baseline FR
    [tmp_r, tmp_p] = corr((1:numel(baselineFRstats{ii}))',baselineFRstats{ii},'type','Pearson');
    r = [r; tmp_r];
    p = [p; tmp_p];
end


%*****************************************************%
% Plotting the cone weights
figure(plot_counter); set(gcf,'Name','Cone wts') 
subplot(221); plot(conewts_svd(1,:),conewts_svd(2,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');  
xlabel('L'),ylabel('M'); title('Cone weights - ALL cells'); 
subplot(222); plot(conewts_svd(1,hardtoclassifyidx),conewts_svd(2,hardtoclassifyidx),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,LUMidx),conewts_svd(2,LUMidx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,DOidx),conewts_svd(2,DOidx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');  
xlabel('L'), ylabel('M'); title('Cone weights - LUM, DO, HTC cells');
subplot(223); plot(conewts_svd(1,oldhardtoclassifyidx),conewts_svd(2,oldhardtoclassifyidx),'o','MarkerSize',3,'LineWidth',0.5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,oldLUMidx),conewts_svd(2,oldLUMidx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
plot(conewts_svd(1,oldDOidx),conewts_svd(2,oldDOidx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[0 1],'XTick',-1:0.5:1,'YTick',0:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');  
xlabel('L'), ylabel('M'); title('Old cone wt classification - LUM, DO, HTC cells');
subplot(224); plot(S1LMS(1,oldhardtoclassifyidx),S1LMS(2,oldhardtoclassifyidx),'o','MarkerSize',3,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(S2LMS(1,oldhardtoclassifyidx),S2LMS(2,oldhardtoclassifyidx),'o','MarkerSize',3,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
plot(S1LMS(1,oldLUMidx),S1LMS(2,oldLUMidx),'o','MarkerSize',4,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
plot(S2LMS(1,oldLUMidx),S2LMS(2,oldLUMidx),'o','MarkerSize',4,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(S1LMS(1,oldDOidx),S1LMS(2,oldDOidx),'o','MarkerSize',4,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(S2LMS(1,oldDOidx),S2LMS(2,oldDOidx),'o','MarkerSize',4,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[-1 1],'XTick',-1:0.5:1,'YTick',-1:0.5:1,'Tickdir','out'); 
plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 -1],[0 0],'k'); plot([1 0],[0 -1],'k'); plot([0 -1],[-1 0],'k');  xlabel('L'), ylabel('M');
title('Subunit cone wts & OLD classification '); 
plot_counter = plot_counter + 1;


% Plotting the results for SVD based cone weight classification including the PC1 z-scores 
figure(plot_counter);
subplot(321); histogram(RSSEisoresp_medianofratios(hardtoclassifyidx),logspace(-1,3,31),'displaystyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(median(RSSEisoresp_medianofratios(hardtoclassifyidx)),20,'v','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); 
histogram(RSSEisoresp_medianofratios(LUMidx),logspace(-1,3,31),'displaystyle','stairs','EdgeColor',[0 0 0],'Linewidth',2); 
plot(median(RSSEisoresp_medianofratios(LUMidx)),20,'v','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
histogram(RSSEisoresp_medianofratios(DOidx),logspace(-1,3,31),'displaystyle','stairs','EdgeColor',[1 0 0],'Linewidth',2);
plot(median(RSSEisoresp_medianofratios(DOidx)),17,'v','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','XScale','log','Xlim',[0.1 1000],'XTick',[0.1 1 10 100 1000]); ylabel('Count'); title('Isoresponse'); xlabel('median CV ratio of errors'); axis square; hold off;
subplot(322); histogram(100*(Withinsubunits_medianofdifferences(hardtoclassifyidx)),linspace(-2,8,31),'displaystyle','stairs','EdgeColor',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(median(100*(Withinsubunits_medianofdifferences(hardtoclassifyidx))),20,'v','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
histogram(100*(Withinsubunits_medianofdifferences(LUMidx)),linspace(-2,8,31),'displaystyle','stairs','EdgeColor',[0 0 0],'Linewidth',2);
plot(median(100*(Withinsubunits_medianofdifferences(LUMidx))),20,'v','MarkerSize',8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
histogram(100*(Withinsubunits_medianofdifferences(DOidx)),linspace(-2,8,31),'displaystyle','stairs','EdgeColor',[1 0 0],'Linewidth',2);
plot(median(100*(Withinsubunits_medianofdifferences(DOidx))),17,'v','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-2 10]); xlabel('median CV GQM-GLM AUROC'); ylabel('Count'); title('Within subunits'); axis square; hold off;
subplot(323); plot(z_scores(hardtoclassifyidx),RSSEisoresp_medianofratios(hardtoclassifyidx),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(z_scores(LUMidx),RSSEisoresp_medianofratios(LUMidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(z_scores(DOidx),RSSEisoresp_medianofratios(DOidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out'); xlabel('Z scores'); ylabel('(Lin/Quad) errors');axis square; hold off;
subplot(324); plot(vals(hardtoclassifyidx),RSSEisoresp_medianofratios(hardtoclassifyidx),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(vals(LUMidx),RSSEisoresp_medianofratios(LUMidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
plot(vals(DOidx),RSSEisoresp_medianofratios(DOidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out'); xlabel('Z-score Percentiles'); ylabel('(Lin/Quad) errors');axis square; hold off;
subplot(325); plot(abs(r(hardtoclassifyidx)),RSSEisoresp_medianofratios(hardtoclassifyidx),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(abs(r(LUMidx)),RSSEisoresp_medianofratios(LUMidx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
plot(abs(r(DOidx)),RSSEisoresp_medianofratios(DOidx),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[0 0.6],'Ylim',[0.1 1000],'YScale','log'); xlabel('Spearman-r baseline FR'); ylabel('median CV ratio of errors');
subplot(326); plot(abs(r(hardtoclassifyidx)),100*(Withinsubunits_medianofdifferences(hardtoclassifyidx)),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]); hold on;
plot(abs(r(LUMidx)),100*(Withinsubunits_medianofdifferences(LUMidx)),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); 
plot(abs(r(DOidx)),100*(Withinsubunits_medianofdifferences(DOidx)),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[0 0.6],'Ylim',[-2 8]); xlabel('Spearman-r baseline FR'); ylabel('median CV GQM-GLM AUROC');
plot_counter = plot_counter + 1;

% Need to plot the STA/basis vectors of LUM, DO and HTC cells 
sequence_svd = [LUMidx DOidx hardtoclassifyidx];
sequence_subunitswts = [oldLUMidx oldDOidx oldhardtoclassifyidx];
figure(plot_counter); set(gcf,'Name','SVD');
figure(plot_counter+1); set(gcf,'Name','SVD STA'); 
figure(plot_counter+2); set(gcf,'Name','Subunit wts'); 
for ii = 1:numcells
    
    % Plotting the SVD based classified cells 
    vec1 = basisvec{sequence_svd(ii)};
    figure(plot_counter); subplot(numsubplots,numsubplots,ii); image(vec1); set(gca,'XTick',[],'YTick',[]); axis square; 
    if ismember(sequence_svd(ii),LUMidx)
        set(gca,'XColor','g','YColor','g');
    elseif ismember(sequence_svd(ii),DOidx)
        set(gca,'XColor','r','YColor','r');
    else
    end
    
    % Plotting the SVD STA 
    vec2 = SVDSTA{sequence_svd(ii)};
    figure(plot_counter+1); subplot(numsubplots,numsubplots,ii); image(vec2); set(gca,'XTick',[],'YTick',[]); axis square; 
    if ismember(sequence_svd(ii),LUMidx)
        set(gca,'XColor','g','YColor','g');
    elseif ismember(sequence_svd(ii),DOidx)
        set(gca,'XColor','r','YColor','r');
    else
    end
    
    % Plotting the subunit cone wt based classified cells
    vec3 = basisvec{sequence_subunitswts(ii)};
    figure(plot_counter+2); subplot(numsubplots,numsubplots,ii); image(vec3); set(gca,'XTick',[],'YTick',[]); axis square;
    if ismember(sequence_subunitswts(ii),newLUMidx)
        set(gca,'XColor','g','YColor','g');
    elseif ismember(sequence_subunitswts(ii),newDOidx)
        set(gca,'XColor','r','YColor','r');
    else
    end
    
end
plot_counter = plot_counter + 3;

