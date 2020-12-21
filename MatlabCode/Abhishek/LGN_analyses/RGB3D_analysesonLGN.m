% A derivative of PopCellAnalysis_RGB3D_crossval.m
% Testing GLM and GQM analyses on LGN 
% Author - Abhishek De, 3/20

close all; clearvars;
plot_counter = 1;
[Input_List,unitidx,~,neuron_id] = fnamesFromTxt('WhiteNoiseLGN_forIS');

nspikes = [];
AUROClin = cell(1,numel(Input_List)); 
AUROCquad = cell(1,numel(Input_List));

files_not_working = [];
files_not_working_idxs = [];
spikeoptions = [{'sig001a'}; {'sig001b'}];
numsubplots = ceil(sqrt(numel(Input_List)));
plotSTAs = 1;
conewtsLGN_svd = [];
numstixelsperpixels = [];
numpixelsperside = [];

for ii = 1:numel(Input_List)
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
    end
    if flag 
        continue;
    end
    stro = WN;    
    global maskidx spikeidx nstixperside ngammasteps bkgnd_ridx bkgnd_gidx bkgnd_bidx seedidx
    global correctidx nframesidx stimonidx fponidx stimoffidx fpacqidx muidxs spikename noisetypeidx
    global sigmaidxs latencyidx basisvecidx msperframe basisvecdiridx neurothreshidx maxT xx yy gammaTable
    spikename = char(spikeoptions(unitidx(ii)));
    maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
    spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
    nstixperside = stro.sum.exptParams.nstixperside;
    ngammasteps = 2^16; % 65536
    bkgnd_ridx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_r'));
    bkgnd_gidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_g'));
    bkgnd_bidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_b'));
    seedidx = strcmp(stro.sum.trialFields(1,:),'seed');
    correctidx = strcmp(stro.sum.trialFields(1,:),'correct');
    nframesidx = strcmp(stro.sum.trialFields(1,:),'num_frames');
    noisetypeidx = strcmp(stro.sum.trialFields(1,:),'noise_type');
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
    msperframe = 1000/stro.sum.exptParams.framerate;
    gammaTable = stro.sum.exptParams.gamma_table;
    gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
    
    ntrials = size(stro.trial,1);
    maxT = 15; % this represents the temporal part in the spatiotemporal receptive field
    xx = linspace(stro.sum.exptParams.gauss_locut/1000, stro.sum.exptParams.gauss_hicut/1000,ngammasteps); % xx represents the probabilities. For more info, have a look at the MATLAB 'norminv' function.
    yy = norminv(xx'); % defining norminv to extract the values for which the cdf values range between gauss_locut and gauss_hicut
    if isfield(stro.sum.exptParams,'nrepframes')
        if isnan(stro.sum.exptParams.nrepframes)
            nvblsperstimupdate = 1;
        else
            nvblsperstimupdate = stro.sum.exptParams.nrepframes;
        end
    else
        nvblsperstimupdate = 1;
    end
    
    % Storing some stixelsperpixels and nstixperside
    numpixelsperside = [numpixelsperside; stro.sum.exptParams.nstixperside];
    numstixelsperpixels = [numstixelsperpixels; stro.sum.exptParams.npixperstix];
    
    % Reconstructing the M matrix
    fundamentals = stro.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = stro.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    
    % Getting the background rgb/lms
    ridx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_r'));
    gidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_g'));
    bidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_b'));
    bkgndRGB = [mode(WN.trial(:,ridx)), mode(WN.trial(:,gidx)), mode(WN.trial(:,bidx))];
    bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
    bkgndlms = M*bkgndrgb;
    Mrgbtocc = diag(1./bkgndlms)*M; % M can be considered to be in cone excitation differences
    Mrgbtocc = inv(Mrgbtocc');
    
    % Extracting the peakframe 
    out_gun = getWhtnsStats(stro,maxT,'STAmex', {nstixperside^2, 3, maxT}, spikename);
    STAs = out_gun{1}/out_gun{3}; 
    nspikes = [nspikes; out_gun{3}]; 
    clear out_gun;
    [~,peakframe] = find(sum(STAs.^2,1) == max(sum(STAs.^2,1))); 
    tmp = sum(reshape(STAs(:,peakframe).^2,[nstixperside nstixperside 3]),3);
    [~,peakpixelloc] = max(tmp(:));
    peakpixelloc = [peakpixelloc; peakpixelloc+nstixperside.^2; peakpixelloc+2*nstixperside.^2];
    
    % Calculating the SVD derived cone weights 
    peakenergy = sum(STAs.^2,1) == max(sum(STAs.^2,1));
    id = find(peakenergy==1); 
    if id~=1
        peakenergy(id-1)= 1;
    end
    if id <=maxT-1
        peakenergy(id+1)=1;
    end
    
    STAweights = sqrt(sum(STAs(:,peakenergy).^2));
    STAweights = STAweights./sum(STAweights);
    tmpSTA = reshape(STAs(:,peakenergy)*STAweights',[nstixperside.^2 3]);
    %[u1,s,v1] = svd(tmpSTA);
%     RGB_svd = v1(:,1);
    RGB_svd = tmpSTA(peakpixelloc(1),:)';
    tmp_conewts_svd = Mrgbtocc * RGB_svd;
    tmp_conewts_svd = tmp_conewts_svd./repmat(sum(abs(tmp_conewts_svd),1),[3 1]);
%     tmp_conewts_svd = tmp_conewts_svd .* repmat(sign(tmp_conewts_svd(2,:)),[3 1]);
    conewtsLGN_svd = [conewtsLGN_svd tmp_conewts_svd];

    % Acquiring all the R,G,B stimuli at the peak pixel location
    cum_rgbs = [];
    cum_n = [];
    for i = 1:ntrials %get values and insert into given column
        nframes = stro.trial(i,nframesidx);
        if nframes == 0, continue; end
        if isnan(nframes), continue; end
        seed = stro.trial(i,seedidx);
        noisetype = stro.trial(i,noisetypeidx);
        mu = stro.trial(i,muidxs)/1000;
        sigma = stro.trial(i,sigmaidxs)/1000;
        
        (noisetype == 1)  % Gaussian gun
        invnormcdf = bsxfun(@plus, bsxfun(@times, yy, sigma), mu);
        randnums = getEJrandnums(3*nstixperside^2*nframes, seed);
        randnums = reshape(randnums, [nstixperside^2*3, nframes]);
        for j = 1:3
            idxs = (1:nstixperside^2)+nstixperside^2*(j-1);
            randnums(idxs,:) = reshape(invnormcdf(randnums(idxs,:)+1,j),[length(idxs),nframes]); % order in drawing on each frame frame is: red (all pixels), green (all pixels), blue (all pixels).
        end
        rgbs = randnums(peakpixelloc,:);
        t_stimon = stro.trial(i, stimonidx);
        spiketimes = (stro.ras{i,spikeidx}-t_stimon)*1000;  % converting to ms
        frametimes = linspace(0, nframes*msperframe*nvblsperstimupdate, nframes)+(msperframe*nvblsperstimupdate/2)';
        spiketimes(spiketimes < maxT*msperframe) = [];
        spiketimes(spiketimes > frametimes(end)) = [];
        n = hist(spiketimes, frametimes);
        % Stroring all the rgbs and spikes
        cum_rgbs = [cum_rgbs rgbs];
        cum_n = [cum_n n];  
    end
    
    basisvecST = STAs(peakpixelloc,:); % temporal phosphor dynamics of basisvec 1
    stindices = find(cum_n>0);
    stindices = stindices - peakframe+1;
    rawrgb_onbasisvec = cum_rgbs;
    response = -1*ones(1,size(rawrgb_onbasisvec,2)); response(stindices) = 1;
    
    % Need to implement the cross-validation here: 2 cross vals for 2 subunits
    response(response==-1) = 0; 
    response = logical(response);
    C = cvpartition(response,'KFold',10);
    tmp_AUROClin = []; tmp_AUROCquad = [];

    % Comparing GLM and GQM fits
    for jj = 1:C.NumTestSets
        mdl1lin =  fitglm(rawrgb_onbasisvec(:,C.training(jj))',response(C.training(jj))','linear','Distribution','binomial','Link','logit');
        mdl1quad =  fitglm(rawrgb_onbasisvec(:,C.training(jj))',response(C.training(jj))','quadratic','Distribution','binomial','Link','logit');
        
        % Performing some additional analyses on the model fits (GLM and GQM)
        predlin = predict(mdl1lin,rawrgb_onbasisvec(:,C.test(jj))'); % perdiction from GLM subunit 1
        predquad = predict(mdl1quad,rawrgb_onbasisvec(:,C.test(jj))'); % perdiction from GQM subunit 1
        
        % Quantifying accuracy using AUROC 
        tmp_AUROClin = [tmp_AUROClin; rocN(predlin(response(C.test(jj))),predlin(~response(C.test(jj))))];
        tmp_AUROCquad = [tmp_AUROCquad; rocN(predquad(response(C.test(jj))),predquad(~response(C.test(jj))))];
    end
    AUROClin{ii} = tmp_AUROClin;
    AUROCquad{ii} = tmp_AUROCquad;   
    
    % Plotting the STAs
    if plotSTAs
        im = reshape(STAs(:,peakframe),[nstixperside nstixperside 3]);
        im = 0.5*im/(max(abs(im(:)))+eps) + 0.5;
        figure(plot_counter); subplot(numsubplots,numsubplots,ii); image(im); axis square; set(gca,'XTick',[],'YTick',[]);
    end
    
end
plot_counter = plot_counter + 1;


%% For storing median of differences/ratios
% Changing the metric to conesignal NLI
% AD, 8/20
ConesignalNLI_LGN = [];

for ii = 1:numel(AUROClin) 
    Error_lin = (1-AUROClin{ii});
    Error_quad = (1-AUROCquad{ii});
    ConesignalNLI_LGN = [ConesignalNLI_LGN; log10(median(Error_lin./Error_quad))];   
end

% Plotting the cone weights along with the NLI
figure(plot_counter); 
subplot(121); plot(conewtsLGN_svd(1,:),conewtsLGN_svd(2,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]); hold on;
axis equal; set(gca,'Xlim',[-1 1],'Ylim',[-1 1],'XTick',-1:0.5:1,'YTick',-1:0.5:1,'Tickdir','out'); plot([-1 0],[0 1],'k'); plot([0 1],[1 0],'k'); plot([-1 1],[0 0],'k');
plot([-1 0],[0 -1],'k'); plot([0 1],[-1 0],'k'); plot([0 0],[-1 1],'k');
xlabel('L'),ylabel('M'); title('Cone weights - LGN'); 
subplot(122); histogram(ConesignalNLI_LGN,-0.04:0.005:0.1,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); hold on;
plot(median(ConesignalNLI_LGN),60,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
set(gca,'Tickdir','out','Xlim',[-0.03 0.06],'XTick',-0.03:0.01:0.06,'Ylim',[0 60],'YTick',[0 15 30 45 60]); xlabel('Cone signal NLI'); ylabel('Count'); title('LGN'); axis square; hold off;
plot_counter = plot_counter + 1;

save ConesignalNLI_LGN ConesignalNLI_LGN
save AUROClin_LGN AUROClin
save AUROCquad_LGN AUROCquad
save conewtsLGN_svd conewtsLGN_svd


[r,p] = corr(abs(conewtsLGN_svd(3,:))', ConesignalNLI_LGN, 'type', 'Spearman');