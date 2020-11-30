% Creating figures for Greg's Janelia Color Vision talk
% Author - Abhishek De, 3/19
%% Section 1: Gabor and DOG fits to the spatial RFs
close all; clearvars;
load Output_ListWN.mat
load Singleopponent.mat
load Crescenterror.mat
load Gaborerror.mat
load DOGerror.mat
load modelfits.mat
plot_counter = 1;
N = 100;
numcrescentparams = 8;
numgaborparams = 8;
numDOGparams = 6;
crit = chi2inv(0.9999,300);
Z = cell2mat(Output_List(:,7));
Zmax = max(Z(:,2:7),[],2);
Z_cellsofinterest = Zmax>crit;
Output_List(~Z_cellsofinterest,:) = [];
NLI = cell2mat(Output_List(:,13));
simplecells = NLI<0;
Singleopponent = logical(Singleopponent);
ind = find(~Singleopponent & simplecells);

CrescentBIC = N*log(Crescenterror(~Singleopponent & simplecells,:)/N) + numcrescentparams*log(N);
GaborBIC = N*log(Gaborerror(~Singleopponent & simplecells,:)/N) + numgaborparams*log(N);
DOGBIC = N*log(DOGerror(~Singleopponent & simplecells,:)/N) + numDOGparams*log(N);
BICs = [CrescentBIC GaborBIC DOGBIC];
[minval,I] = min(BICs,[],2);

% calculating the M matrix
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
M = inv(M');

RGB_svd = cell2mat(Output_List(~Singleopponent & simplecells,5)');
conewts_svd = M * RGB_svd;
conewts_svd = conewts_svd./repmat(sum(abs(conewts_svd),1),[3 1]);
conewts_svd = conewts_svd .* repmat(sign(conewts_svd(2,:)),[3 1]);

thresh = 0.7;
LumIds_conewts = find(conewts_svd(1,:) + conewts_svd(2,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==2 & conewts_svd(1,:)>0.1 & conewts_svd(2,:)>0.1);
ColorOpponentIds_conewts = find(conewts_svd(2,:) - conewts_svd(1,:) >thresh & sum(sign(conewts_svd(1:2,:)),1)==0 & sqrt((conewts_svd(2,:)-0.5).^2 + (conewts_svd(1,:)+0.5).^2)<0.3);
Sconedominated_conewts = find(abs(conewts_svd(3,:))>1-thresh);
Other_conewts = 1:size(conewts_svd,2); Other_conewts([LumIds_conewts ColorOpponentIds_conewts Sconedominated_conewts]) = [];

% Next, I am plotting the SVD derived luminance, color-opponent and other cells
ind = find(~Singleopponent & simplecells);
Lumind = ind(LumIds_conewts); Lumsubplot = ceil(sqrt(numel(Lumind)));
COind = ind(ColorOpponentIds_conewts); Sconeind = ind(Sconedominated_conewts);
DOind = [COind; Sconeind]; DOsubplot = ceil(sqrt(numel(DOind)));

% Figure 1: Need to handpick 16 Lum & 16 DO cells
% plotting STAs of lumninance simple cells
figure(plot_counter); set(gcf,'Name','L+M');
count = 1;
indices = reshape(1:64,[8 8]);
STA_ind = indices(1:4,1:4); STA_ind = STA_ind(:);
RF_ind = indices(1:4,5:end); RF_ind = RF_ind(:);
Gaborfit_ind = indices(5:end,1:4); Gaborfit_ind = Gaborfit_ind(:);
DOGfit_ind = indices(5:end,5:end); DOGfit_ind = DOGfit_ind(:);
for ii = [6 8 9 11 13 19 26 28 31 32 33 44 45 47 60 66]
    idx = Lumind(ii);
    tmp_vec_gun = Output_List{idx,2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*tmp_vec_gun + 0.5;
    im = reshape(im,[10 10 3]);
    figure(plot_counter); subplot(8,8,STA_ind(count)); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
    im = Output_List{idx,4};
    figure(plot_counter); subplot(8,8,RF_ind(count)); image(255*(im./(2*max(abs(im(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255))
    im = modelfits{idx,2};
    figure(plot_counter); subplot(8,8,Gaborfit_ind(count));image(255*(im./(2*max(abs(im(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255))
    im = modelfits{idx,3};
    figure(plot_counter); subplot(8,8,DOGfit_ind(count)); image(255*(im./(2*max(abs(im(:))))+.5));  set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255))
    count = count + 1;
end
plot_counter = plot_counter + 1;

% plotting STAs of DO cells
figure(plot_counter); set(gcf,'Name','DO');
count = 1;
for ii = [5 11 16 26 28 30 39 60 86 100 107 115 117 119 121 127]
    idx = DOind(ii);
    tmp_vec_gun = Output_List{idx,2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*tmp_vec_gun + 0.5;
    im = reshape(im,[10 10 3]);
    figure(plot_counter); subplot(8,8,STA_ind(count)); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
    im = Output_List{idx,4};
    figure(plot_counter); subplot(8,8,RF_ind(count)); image(255*(im./(2*max(abs(im(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255))
    im = modelfits{idx,2};
    figure(plot_counter); subplot(8,8,Gaborfit_ind(count));image(255*(im./(2*max(abs(im(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255))
    im = modelfits{idx,3};
    figure(plot_counter); subplot(8,8,DOGfit_ind(count)); image(255*(im./(2*max(abs(im(:))))+.5));  set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255))
    count = count + 1;
end
plot_counter = plot_counter + 1;

% Figure 2: plotting the BIC (Gabor) vs BIC (DOG). Simple cells in one color, DO cells(L-M and S-dominated) in another
fnames = Output_List(~Singleopponent & simplecells,1);
figure(plot_counter); set(gcf,'Name','Scatter plot DoG BIC vs Gabor BIC: simple & DO cells');
subplot(131); hold on;
for ii = 1:length(fnames)
    if any(ismember(LumIds_conewts,ii))
        h(ii) = plot(DOGBIC(ii), GaborBIC(ii),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', char(fnames{ii,:}),''')']);
    elseif any(ismember(ColorOpponentIds_conewts,ii))
        h(ii) = plot(DOGBIC(ii), GaborBIC(ii),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', char(fnames{ii,:}),''')']);
    elseif any(ismember(Sconedominated_conewts,ii))
        h(ii) = plot(DOGBIC(ii), GaborBIC(ii),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', char(fnames{ii,:}),''')']);
    end
end
plot([-1200 -200],[-1200 -200],'k');
xlabel('DoG BIC'); ylabel('Gabor BIC'); set(gca,'TickDir','out','Ylim',[-1200 -200],'Xlim',[-1200 -200],'XTick',-1200:200:-200,'YTick',-1200:200:-200); axis square; hold off;
subplot(132); hold on;
for ii = 1:length(fnames)
    if any(ismember(LumIds_conewts,ii))
        h(ii) = plot(CrescentBIC(ii), GaborBIC(ii),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', char(fnames{ii,:}),''')']);
    elseif any(ismember(ColorOpponentIds_conewts,ii))
        h(ii) = plot(CrescentBIC(ii), GaborBIC(ii),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', char(fnames{ii,:}),''')']);
    elseif any(ismember(Sconedominated_conewts,ii))
        h(ii) = plot(CrescentBIC(ii), GaborBIC(ii),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', char(fnames{ii,:}),''')']);
    end
end
plot([-1200 -200],[-1200 -200],'k');
xlabel('Crescent BIC'); ylabel('Gabor BIC'); set(gca,'TickDir','out','Ylim',[-1200 -200],'Xlim',[-1200 -200],'XTick',-1200:200:-200,'YTick',-1200:200:-200); axis square; hold off;
fnames = Output_List(~Singleopponent & simplecells,1);
figure(plot_counter); set(gcf,'Name','Scatter plot DoG BIC vs Gabor BIC: simple & DO cells');
subplot(131); hold on;
for ii = 1:length(fnames)
    if any(ismember(LumIds_conewts,ii))
        h(ii) = plot(DOGBIC(ii), GaborBIC(ii),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', char(fnames{ii,:}),''')']);
    elseif any(ismember(ColorOpponentIds_conewts,ii))
        h(ii) = plot(DOGBIC(ii), GaborBIC(ii),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', char(fnames{ii,:}),''')']);
    elseif any(ismember(Sconedominated_conewts,ii))
        h(ii) = plot(DOGBIC(ii), GaborBIC(ii),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', char(fnames{ii,:}),''')']);
    end
end
plot([-1200 -200],[-1200 -200],'k');
xlabel('DoG BIC'); ylabel('Gabor BIC'); set(gca,'TickDir','out','Ylim',[-1200 -200],'Xlim',[-1200 -200],'XTick',-1200:200:-200,'YTick',-1200:200:-200); axis square; hold off;
subplot(133); hold on;
for ii = 1:length(fnames)
    if any(ismember(LumIds_conewts,ii))
        h(ii) = plot(DOGBIC(ii), CrescentBIC(ii),'o','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', char(fnames{ii,:}),''')']);
    elseif any(ismember(ColorOpponentIds_conewts,ii))
        h(ii) = plot(DOGBIC(ii), CrescentBIC(ii),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', char(fnames{ii,:}),''')']);
    elseif any(ismember(Sconedominated_conewts,ii))
        h(ii) = plot(DOGBIC(ii), CrescentBIC(ii),'o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]);
        set(h(ii),'ButtonDownFcn',['disp(''', char(fnames{ii,:}),''')']);
    end
end
plot([-1200 -200],[-1200 -200],'k');
xlabel('DOG BIC'); ylabel('Crescent BIC'); set(gca,'TickDir','out','Ylim',[-1200 -200],'Xlim',[-1200 -200],'XTick',-1200:200:-200,'YTick',-1200:200:-200); axis square; hold off;
plot_counter = plot_counter + 1;

plot_counter = plot_counter + 1;


%% Section 2: Generating isoresponse contour examples: 1 DO cell & 1 luminance simple cell
load filename_c.mat
load filename_l.mat
filename = [filename_c; filename_l];
load newLMidx.mat
load newOCidx.mat
load newLUMidx.mat
DOidx = [newOCidx newLMidx];
DO_filenames = filename([DOidx]);
lum_filenames = filename([newLUMidx]);
% good luminance cells : 2,5,11,12,14,15,17,18,19,20,22
if ~exist('plot_counter')
    plot_counter = 1;
end
fnames = ['M090817001.nex'];%; lum_filenames{22}];
for zz = 1:size(fnames,1)
    stro = nex2stro(findfile(fnames(zz,:)));
    global spikename maskidx spikeidx neurothreshidx nstixperside ngammasteps seedidx nframesidx correctidx
    global fponidx stimoffidx stimonidx muidxs sigmaidxs basisvecidx weightsidx fpacqidx targetspikerateidx basisvecdiridx latencyidx reversalflagidx parentverticesidx
    global msperframe ntrials maxT xx yy M linepredtol stepsizescale stepsize nreversals oogscale
    spikename = 'sig001a';
    maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
    spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
    basisvecidx = strcmp(stro.sum.rasterCells(1,:),'basis_vec');
    weightsidx = strcmp(stro.sum.rasterCells(1,:),'weights');
    parentverticesidx = strcmp(stro.sum.rasterCells(1,:),'parentvertices');
    nstixperside = stro.sum.exptParams.nstixperside;
    ngammasteps = 2^16; % 65536
    linepredtol = stro.sum.exptParams.linepredtol;
    stepsizescale = stro.sum.exptParams.stepsizescale;
    stepsize = stro.sum.exptParams.stepsize;
    nreversals = stro.sum.exptParams.nreversals;
    oogscale = stro.sum.exptParams.oogscale;
    seedidx = strcmp(stro.sum.trialFields(1,:),'seed');
    nframesidx = strcmp(stro.sum.trialFields(1,:),'num_frames');
    stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
    stimoffidx = strcmp(stro.sum.trialFields(1,:),'stim_off');
    fponidx = strcmp(stro.sum.trialFields(1,:),'fp_on');
    fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
    basisvecdiridx = strcmp(stro.sum.trialFields(1,:),'weights_idx');
    neurothreshidx = strcmp(stro.sum.trialFields(1,:),'neurothresh'); % when exactly the neurothresh trials started
    targetspikerateidx = strcmp(stro.sum.trialFields(1,:),'targetspikerate');
    correctidx = strcmp(stro.sum.trialFields(1,:),'correct');
    muidxs = [find(strcmp(stro.sum.trialFields(1,:),'mu1')), ...
        find(strcmp(stro.sum.trialFields(1,:),'mu2')), ...
        find(strcmp(stro.sum.trialFields(1,:),'mu3'))];
    sigmaidxs = [find(strcmp(stro.sum.trialFields(1,:),'sigma1')), ...
        find(strcmp(stro.sum.trialFields(1,:),'sigma2')), ...
        find(strcmp(stro.sum.trialFields(1,:),'sigma3'))];
    latencyidx = strcmp(stro.sum.trialFields(1,:),'latency');
    reversalflagidx = strcmp(stro.sum.trialFields(1,:),'reversalflag');
    msperframe = 1000/stro.sum.exptParams.framerate;
    ntrials = size(stro.trial,1);
    maxT = 10; % this represents the temporal part in the spatiotemporal receptive field
    xx = linspace(stro.sum.exptParams.gauss_locut/1000, stro.sum.exptParams.gauss_hicut/1000,ngammasteps); % xx represents the probabilities. For more info, have a look at the MATLAB 'norminv' function.
    yy = norminv(xx'); % defining norminv to extract the values for which the cdf values range between gauss_locut and gauss_hicut
    t_offset = stro.trial(end,latencyidx)/1000;
    
    % Obtaining the M matrix, code extracted from Greg, fitting a cubic spline
    % using the command 'spline'. 'SplineRaw' only availabe through
    % psychtoolbox which I currently don't have now.
    fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
    mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
    M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
    M = inv(M');
    
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
    
    if isfield(stro.sum.exptParams,'nrepframes')
        if ~isnan(stro.sum.exptParams.nrepframes)
            nvblsperstimupdate = stro.sum.exptParams.nrepframes;
        else
            nvblsperstimupdate = 1;
        end
    else
        nvblsperstimupdate = 1;
    end
    
    mask_changes = reshape(mask_changes , 2, []);
    for ii = 1:2
        idxs = zeros(size(stro.trial,1),1);
        idxs(mask_changes(1,ii):mask_changes(2,ii)) = 1;
        idxs = logical(idxs);
        WN = stro;
        WN.ras(~idxs,:) = []; WN.trial(~idxs,:) = [];
        if ii == 1
            nrandnumsperchannel = nstixperside^2;
        else
            mask = stro.ras{mask_changes(1,ii),maskidx}; % subunit mask
            mask(mask == 0) = Inf;
            [stIdxs,~,~] = unique(mask); % now the Infs map to nsubunits+1
            num_subunits = length(stIdxs)-any(isinf(stIdxs)); % nsubunits, like subunits A and B
            mask(isinf(mask)) = num_subunits + 1;
            nrandnumsperchannel = num_subunits;
        end
        out_gun = getWhtnsStats(WN,maxT,'STCOVmex',{nrandnumsperchannel,3,maxT},spikename);
        STS_gun = out_gun{1}; STCross_gun = out_gun{2}; nspikes_gun = out_gun{3}; clear out_gun;
        STAs_gun = STS_gun/nspikes_gun;
        energy = sum(STAs_gun.^2,1);
        peakframe = energy == max(energy);
        id = find(peakframe==1);
        latency = find(peakframe)*1000/stro.sum.exptParams.framerate;
        if id~=1
            peakframe(id-1)= 1;
        end
        if id <=maxT-1
            peakframe(id+1)=1;
        end
        STAweights = sqrt(sum(STAs_gun(:,peakframe).^2));
        STAweights = STAweights./sum(STAweights);
        tmpSTA = STAs_gun(:,peakframe)*STAweights'; % weighted combination of peak and its adjacent frames
        tmpSTA2 = STAs_gun(:,id); % just the peak frame
        if ii == 1
            [u1,~,v1] = svd(reshape(tmpSTA,[nstixperside^2 3])');
            SpatialRF = reshape(v1(:,1),[nstixperside nstixperside]);
        elseif ii == 2
            G_mask = [mask; mask+max(mask); mask+2*max(mask)];
            tmpSTA = expand_vector(tmpSTA,num_subunits,G_mask,1);
            tmpSTA2 = expand_vector(tmpSTA2,num_subunits,G_mask,1);
        end
        normfactor = 0.5/(max(abs(tmpSTA(:)))+0.01);
        tmpSTA = normfactor*tmpSTA + 0.5;
        tmpSTA = reshape(tmpSTA,[nstixperside nstixperside 3]);
        normfactor = 0.5/(max(abs(tmpSTA2(:)))+0.01);
        tmpSTA2 = normfactor*tmpSTA2 + 0.5;
        tmpSTA2 = reshape(tmpSTA2,[nstixperside nstixperside 3]);
        figure(plot_counter),subplot(3,3,3*ii-2); image(tmpSTA); set(gca,'XTick',[],'YTick',[]); axis square;
        subplot(3,3,3*ii-1); plot(energy,'k','Linewidth',2); set(gca,'tickdir','out'); axis square;
        subplot(3,3,3*ii); image(tmpSTA2); set(gca,'XTick',[],'YTick',[]); axis square;
    end
    
    % Determining when Neurothresh mode was active, plotting the basis vector, working correctly
    neurothreshmode = stro.trial(:,neurothreshidx);
    basisvec_dropidx = inds(end);
    neurothresh_startidx = max(find(neurothreshmode==1,1),basisvec_dropidx+1); % There is possibility of pressing the Neurothresh button multiple times
    num_targetspikerates = unique(stro.trial(neurothresh_startidx:end,targetspikerateidx));
    vect = stro.ras{basisvec_dropidx,basisvecidx};
    basisvec_size = nstixperside*nstixperside*3;
    numvect = (numel(vect)/basisvec_size)-1;
    basisvec = cell(1,numvect);
    % Actual basis vec
    for ii = 1:numvect
        tmp_vec = vect((ii-1)*basisvec_size+1:basisvec_size*ii) + vect(numvect*basisvec_size+1:basisvec_size*(numvect+1));
        basisvec{ii} = reshape(tmp_vec,[nstixperside nstixperside 3]);
    end
    bkgnd_monitor = reshape(vect(numvect*basisvec_size+1:basisvec_size*(numvect+1)),[nstixperside nstixperside 3]);
    tmp1 = unique(basisvec{1}-bkgnd_monitor,'stable');
    tmp2 = unique(basisvec{2}-bkgnd_monitor,'stable');
    tmpLMS1 = M*tmp1(tmp1~=0);
    S1LMS = tmpLMS1/sum(abs(tmpLMS1));
    tmpLMS2 = M*tmp2(tmp2~=0);
    S2LMS = tmpLMS2/sum(abs(tmpLMS2));
    figure(plot_counter); subplot(3,3,7);
    h = bar([S1LMS(1) S2LMS(1); S1LMS(2) S2LMS(2); S1LMS(3) S2LMS(3)]); set(h(1),'FaceColor',[0 0.5 1]); set(h(2),'FaceColor',[0.5 0.5 0.5]);
    set(gca,'XTick',[1 2 3],'XTickLabel',{'L','M','S'},'TickDir','out'); axis square;
    subplot(3,3,8); image(255*(SpatialRF./(2*max(abs(SpatialRF(:))))+.5)); set(gca,'Xlim',[0 10],'Ylim',[0 10],'XTick',[],'YTick',[]); axis image; colormap(gray(255))
    plot_counter = plot_counter + 1;
    
    % plotting the rasterplots for neurothresh trials
    tmp_wts = [];
    spikecounts = [];
    idxs = find(~isnan(stro.trial(:,correctidx)) & stro.trial(:,targetspikerateidx)==num_targetspikerates(1));
    idxs(idxs<=neurothresh_startidx) = [];
    different_weights = unique(stro.trial(idxs,basisvecdiridx));
    tmp_norm = [];
    tmp_completed_search_alongdir = [];
    staircasetermination_idxs = [];
    for kk = 1:numel(different_weights)
        idxs1 = find(stro.trial(:,basisvecdiridx) == different_weights(kk));
        idxs1(idxs1<neurothresh_startidx) = [];
        raster_data = stro.ras(idxs1,1);
        tmp_norm = [tmp_norm; stro.ras{idxs1(end),weightsidx}'];
        tmp_wts = [tmp_wts;cell2mat(stro.ras(idxs1,weightsidx)')'];
        staircasetermination_idxs = [staircasetermination_idxs; idxs1(end)];
        for ii = 1:size(raster_data,1)
            tmp = raster_data{ii} ;
            spikes = tmp(tmp>stro.trial(idxs1(ii),stimonidx)+t_offset & tmp < stro.trial(idxs1(ii),stimoffidx));
            spikes = spikes - stro.trial(idxs1(ii),stimonidx)-t_offset;
            spikecounts = [spikecounts; numel(spikes)];
        end
        [flag, gamutViolation] = Check_ifcompletely_probed(stro,basisvec,bkgnd_monitor,idxs1,tmp_norm(end,:));
        % flag = 0, incompletely probed
        % flag = 1, completely probed
        % gamutViolation = 1, out of gamut point
        tmp_completed_search_alongdir = [tmp_completed_search_alongdir; flag gamutViolation];
    end
    
    tmp = tmp_norm;
    completed_dir = tmp_completed_search_alongdir;
    probed_dirs = logical(completed_dir(:,1)==1); % only including the directions that have been completely probed
    oog_idx = find(completed_dir(:,1)==1 & completed_dir(:,2)==1); % probed and out of gamut
    not_oog_idx = find(completed_dir(:,1)==1 & completed_dir(:,2)==0);
    fact = 0.5./sqrt(tmp(probed_dirs,1).^2 + tmp(probed_dirs,2).^2); % factor needed to extract unit vector
    [THETA,RHO] = cart2pol(tmp(:,1),tmp(:,2));
    if any(THETA>(3*pi/4))
        allthetas = linspace(-pi,pi,100);
        newtheta = linspace(-pi,pi,101);
    else
        allthetas = linspace(-pi/4,3*pi/4,100);
        newtheta = linspace(-pi/4,3*pi/4,101);
    end
    outofgamut = zeros(size(THETA));
    outofgamut(oog_idx) = 1;
    outofgamut = logical(outofgamut);
    [x_orig, y_orig] = pol2cart(THETA,RHO);
    
    % Fitting the linear model
    initguess1 = [x_orig(not_oog_idx) y_orig(not_oog_idx)]\ones(numel(x_orig(not_oog_idx)),1);
    [final_model1] = linefit_AD2(RHO, THETA,not_oog_idx,outofgamut,initguess1');
    rho1 = 1./(final_model1*[cos(allthetas); sin(allthetas)]);
    LOOGtmp1= rho1<0;
    [x_lin,y_lin] = pol2cart(allthetas(~LOOGtmp1),rho1(~LOOGtmp1));
    initguess3 = [0 0 0 final_model1];
    [final_model3] = quadfit_AD2(RHO, THETA, not_oog_idx,outofgamut,initguess3);
    [final_model3] = conicsectionfit(x_orig(not_oog_idx),y_orig(not_oog_idx),initguess3); % passing on the values to conicsection fit which fits in x-y plane
    [final_model3] = quadfit_AD2(RHO, THETA, not_oog_idx,outofgamut,final_model3); % repassing the values to quadfit which fits in r-theta plane
    [x_quad,y_quad,rho3] = calc_xyvalues(allthetas, final_model3);
    L = rho3>0 & rho3==real(rho3);
    [x_quad2,y_quad2] = pol2cart(newtheta(L),rho3(L)');
    
    % Plotting the figures
    figure(plot_counter), subplot(3,3,5); plot(x_lin,y_lin,'k','Linewidth',2); hold on;
    plot(x_orig(not_oog_idx), y_orig(not_oog_idx),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'PickableParts','none','MarkerEdgeColor',[1 1 1]);
    if any(outofgamut)
        plot(upsample(x_orig(outofgamut),2),upsample(y_orig(outofgamut),2),'color',[0.5 0.5 0.5]);
    end
    set(gca,'XLim',[-2 2],'YLim',[-2 2],'Tickdir','out'); drawnow; axis square; hold off;
    
    figure(plot_counter+1),
    for ii = 1:numel(not_oog_idx)
        h(ii) = plot(x_orig(not_oog_idx(ii)), y_orig(not_oog_idx(ii)),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]); hold on;
        set(h(ii),'ButtonDownFcn',{@dispimage,x_orig(not_oog_idx(ii)),y_orig(not_oog_idx(ii)),basisvec,bkgnd_monitor});%(x_orig(not_oog_idx(ii)),y_orig(not_oog_idx(ii)),basisvec,bkgnd_monitor)});
    end
    plot(x_lin,y_lin,'k','Linewidth',2); hold on;
    
    if any(outofgamut)
        plot(upsample(x_orig(outofgamut),2),upsample(y_orig(outofgamut),2),'color',[0.5 0.5 0.5]);
    end
    set(gca,'XLim',[-2 2],'YLim',[-2 2],'Tickdir','out'); drawnow; axis square; hold off;
    
    % Now, I am going to plot the
    wts  = [-1 1;0 1;1 1;-1 0;1 0;-1 -1;0 -1;1 -1];
    subplotidxs = [1;2;3;4;6;7;8;9];
    for ii = 1:numel(subplotidxs)
        vec = wts(ii,1)*(basisvec{1}-bkgnd_monitor) + wts(ii,2)*(basisvec{2}-bkgnd_monitor);
        normfactor = 0.5/(max(abs(vec(:)))+0.01);
        vec = normfactor*vec + 0.5;
        figure(plot_counter),subplot(3,3,subplotidxs(ii)); image(vec); set(gca,'XTick',[],'YTick',[]); axis square;
    end
    plot_counter = plot_counter + 2;
end

%% Doing some population analyses
load RSSE_linearmodel.mat
load RSSE_quadmodel.mat
Q_rob = log(RSSE_linearmodel./RSSE_quadmodel);
bins = 0:0.5:10;
Q_rob(Q_rob>max(bins)) = max(bins);
figure(plot_counter); set(gcf,'Name','Double-opponent vs Simple');
subplot(311); histogram(Q_rob(newLUMidx),bins,'Normalization','probability','FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]); hold on;
histogram(Q_rob(DOidx),bins,'Normalization','probability','FaceColor',[1 0 0],'EdgeColor',[1 1 1]); 
plot(median(Q_rob(DOidx)),0,'kv','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
plot(median(Q_rob(newLUMidx)),0,'kv','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Ylim',[0 0.6],'YTick',[0:0.3:0.6]); xlabel('Q'); legend('Lum','DO'); ylabel('Proportion of cells'); 
subplot(312); histogram(Q_rob(newLUMidx),bins,'Normalization','probability','FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]); hold on;
plot(median(Q_rob(newLUMidx)),0,'kv','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Ylim',[0 0.6],'YTick',[0:0.3:0.6]); xlabel('Q'); legend('Lum'); ylabel('Proportion of cells'); 
subplot(313); histogram(Q_rob(DOidx),bins,'Normalization','probability','FaceColor',[1 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(Q_rob(DOidx)),0,'kv','MarkerSize',8,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Ylim',[0 0.6],'YTick',[0:0.3:0.6]); xlabel('Q'); legend('DO'); ylabel('Proportion of cells'); 
set(gcf,'renderer','painters');
% subplot(122); plot(log(RSSE_linearmodel(DOidx)),log(RSSE_quadmodel(DOidx)),'o','MarkerSize',6','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 1 1]); hold on;
% plot(log(RSSE_linearmodel(newLUMidx)),log(RSSE_quadmodel(newLUMidx)),'o','MarkerSize',6','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); 
% axis square; set(gca,'Tickdir','out','Xlim',[-15 5],'Ylim',[-15 5]); legend('DO','Lum'); xlabel('Linear SSE'); ylabel('Quadratic SSE'); 
plot_counter = plot_counter + 1;