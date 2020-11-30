% This project is about looking at if the DO spatial structure can be explained by a gabor or a DOG
% Author - Abhishek De, 7/18
close all; clearvars;
load filename_c.mat
load filename_l.mat
load filenamegun.mat
load S1LMS.mat
load S2LMS.mat
load fundamentals.mat 
load mon_spd.mat
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filenamesub = fetch(conn,'SELECT filename FROM WNSubunit');
comments = fetch(conn,'SELECT comments FROM WNSubunit');
mode1 = fetch(conn,'SELECT mode FROM WNSubunit');
close(conn);
filenamesub = filenamesub(strcmp(mode1,'STA'));
Input_List = [filename_c;filename_l];
% Input_List = [filenamesub];
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = reshape(mon_spd,[length(mon_spd)/3,3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;
resize_fact = 1;
% Include Gun noise data for statistical tests in order to estimate the RF size
channels = 3;
Gaborerror = [];
DOGerror = [];
kk = 1;
num_rows = 10; % Number of cells in a figure
C = 5;
resize_fact2 = 1;
numspikes = [];
for ii = 1:10% numel(Input_List)
    global nstixperside
    filename = Input_List{ii}; % acquiring the filename from the List
    Output_List{ii,1} = filename;
    WN=nex2stro(findfile(filename));
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
    eyestart_t = [WN.ras{:,anlgStartTimeidx}]';
    eyesampperiod = 1/WN.sum.analog.storeRates{1};
    gammaTable = WN.sum.exptParams.gamma_table;
    gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
    gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
    invgamma = InvertGamma(gammaTable, 0);
    basisvecidx = strcmp(WN.sum.rasterCells(1,:),'basis_vec');
    
    % Getting the background rgb/lms
    ridx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_r'));
    gidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_g'));
    bidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_b'));
    bkgndRGB = [mode(WN.trial(:,ridx)), mode(WN.trial(:,gidx)), mode(WN.trial(:,bidx))];
    bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
    bkgndlms = M*bkgndrgb;
    mask_changes = [2];
    all_masks = WN.ras(:,maskidx);
    Fx = @(xi) any(isnan(xi)); % function that finds 'NaN' in a cell array
    basisvec_dropidx = [];
    if find(basisvecidx)
        inds = find(cellfun(Fx,WN.ras(:,basisvecidx))==0);  
    else
        inds = [];
    end
    if isempty(inds)
        inds = size(WN.trial,1)-1;
    else
         basisvec_dropidx = inds(end);
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
    spikename = getSpikenum(WN,'first'); % Getting the spikes present in the first channel
    spikeidx = find(strcmp(WN.sum.rasterCells(1,:),spikename));
    maxT = 15;
    tmpstro_gun = WN; % creating a temporary stro that will be used for analysing stimulus presented in gun space
    tmpstro_gun.ras(mask_changes(2,1)+1:end,:) = [];
    tmpstro_gun.trial(mask_changes(2,1)+1:end,:) = [];
    out_gun = getWhtnsStats(tmpstro_gun,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
    STAs_gun = out_gun{1}; STCs_gun = out_gun{2}; nspikes_gun = out_gun{3};
    numspikes = [numspikes; nspikes_gun];
    
    tmpstro_gunsub = WN; % creating a temporary stro that will be used for analysing subunit stimulus presented in gun space
    tmpstro_gunsub.ras(1:mask_changes(2,1),:) = [];
    tmpstro_gunsub.trial(1:mask_changes(2,1),:) = [];
    st_mask = WN.ras{mask_changes(1,2),maskidx}; % subunit mask                    
    st_mask(st_mask == 0) = Inf;
    [stIdxs,~,mask] = unique(st_mask); % now the Infs map to nsubunits+1
    num_subunits = length(stIdxs)-any(isinf(stIdxs)); % nsubunits, like subunits A and B          
    mask = [mask; mask+max(mask); mask+2*max(mask)]; %#ok<AGROW>
    out_gunsub = getWhtnsStats(tmpstro_gunsub,maxT,'STCOVmex', {num_subunits, 3, maxT}, spikename);
    STAs_gunsub = out_gunsub{1}; STCs_gunsub = out_gunsub{2}; nspikes_gunsub = out_gunsub{3};
    STAs_gunsub_STimg = expand_vector(STAs_gunsub,num_subunits,mask,maxT);
    % Code for Statistical testing begins here 
    
    s_gun = std(STAs_gun(:,1));
    STAs_gun_z = STAs_gun./s_gun;
    
    % Spatial map
    grandz = zeros([nstixperside nstixperside]);
    maxzs = [];
    for i = 1:maxT
        tmp_gun = reshape(STAs_gun_z(:,i),[nstixperside nstixperside 3]); % This is the only place in the code where I use data from gun noise
        maxzs = [maxzs; sum(sum(tmp_gun(:,:,1).^2)) sum(sum(tmp_gun(:,:,2).^2)) sum(sum(tmp_gun(:,:,3).^2))];
    end
    peakframe1 = max(maxzs(:,1)) == maxzs(:,1);
    peakframe2 = max(maxzs(:,2)) == maxzs(:,2);
    peakframe3 = max(maxzs(:,3)) == maxzs(:,3);
    peakframe = max(sum(maxzs,2)) == sum(maxzs,2);
    peakframeST = max(sum(STAs_gunsub.^2,1)) == sum(STAs_gunsub.^2,1);
    id = find(peakframe==1);
    if id~=1
        peakframe(id-1)= 1;
    end
    if id <=maxT-1
        peakframe(id+1)=1;
    end
       
    % Now get largest contiguous block
    tmpSTA = reshape(STAs_gun, [nstixperside^2 3 maxT]);
    Output_List{ii,2} = squeeze(sum(tmpSTA(:,:,peakframe),3));
    tmpSTA = permute(tmpSTA, [2 1 3]);
    STAgunmat_pf = Output_List{ii,2};
    [u1,~,v1] = svd(STAgunmat_pf');
    m = reshape(v1(:,1),[nstixperside nstixperside]);
    %     m = m/(max(abs(m(:)))+0.001); 
      
    Output_List{ii,3} = sign(m).*(abs(m).^(1.0));
    [out1,Output_List{ii,4}] = gaborfit_AD(imresize(Output_List{ii,3},resize_fact)); % fitting gabor
    Gaborerror = [Gaborerror; out1.fval];
    [out2,Output_List{ii,5}] = DOGfit(imresize(Output_List{ii,3},resize_fact)); % fitting DOG
    DOGerror = [DOGerror; out2.fval];
    Output_List{ii,6} = u1(:,1);
    
    
    tmp = reshape(squeeze(sum(STAs_gunsub_STimg(:,peakframeST),2)),[nstixperside nstixperside 3]);
    tmp = (0.5*tmp/(max(tmp(:))+0.01))+0.5;
    Output_List{ii,7} = tmp;
    tmp1= inv(M')*STAs_gunsub(1:num_subunits:end,peakframeST); 
    tmp2 = inv(M')*STAs_gunsub(2:num_subunits:end,peakframeST);
    if size(tmp1,2)==15
%         keyboard;
    end
    Output_List{ii,8} = tmp1/norm(tmp1); % LMS weights, subunit 1 
    Output_List{ii,9} = tmp2/norm(tmp2); % LMS weights, subunit 2
%     keyboard;
    figure(kk);
    
    % gunnoise based STA
    tmp_vec_gun = Output_List{ii,2};
    normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
    im = normfactor*Output_List{ii,2} + 0.5;
    im = reshape(im,[nstixperside nstixperside 3]);
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+1); image(im); set(gca,'XTick',[],'YTick',[]); axis square;
    if mod(ii,num_rows) ==1 
        title('Gun');
    end
    
    % Peak Frame
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+2); imagesc(imresize(Output_List{ii,3},resize_fact2)); set(gca,'XTick',[],'YTick',[]); axis square;
    if mod(ii,num_rows) ==1 
        title('RF');
    end
   
    % Gabor Fit R 
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+3); imagesc(Output_List{ii,4}); set(gca,'XTick',[],'YTick',[]); axis square;
    if mod(ii,num_rows) ==1 
        title('Gabor');
    end
    
    % DOG fit 
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+4); imagesc(Output_List{ii,5}); set(gca,'XTick',[],'YTick',[]); axis square;
    if mod(ii,num_rows) ==1 
        title('DOG');
    end
    
    % with subunits
    subplot(num_rows,C,((ii-(kk-1)*num_rows)-1)*C+5); image(Output_List{ii,7}); set(gca,'XTick',[],'YTick',[]); axis square;
    if mod(ii,num_rows) ==1 
        title('sub');
    end
    
    % Subunits 
    
    if mod(ii,num_rows) == 0
        kk = kk + 1;
    end
       
end


%% Need to come up with a summary statistics for which is a better fit, This statistic probably won't work since the models are not nested.
% Finding a way to separate color cells from luminance cells
S1LMS = cell2mat(Output_List(:,8)');
S2LMS = cell2mat(Output_List(:,9)');
idx = find(sum(sign(S1LMS).*sign(S2LMS),1)==-3);
OCidx = idx(find(sign(S1LMS(2,idx)).*sign(S1LMS(3,idx))==1 & sign(S1LMS(1,idx)).*sign(S1LMS(3,idx))==-1));
LMidx = idx(find(sign(S1LMS(1,idx)).*sign(S1LMS(3,idx))==1 & sign(S1LMS(2,idx)).*sign(S1LMS(3,idx))==-1));
LUMidx = idx(find(sign(S1LMS(1,idx)).*sign(S1LMS(2,idx))==1));
hardtoclassifyidx = find(sum(sign(S1LMS).*sign(S2LMS),1)>-3);
DOidx = [OCidx LMidx];
N = nstixperside.^2;
numgaborparams = 8;
numDOGparams = 6;
bins = -150:10:150;
GaborBIC = N*log(Gaborerror/N) + numgaborparams*log(N);
DOGBIC = N*log(DOGerror/N) + numDOGparams*log(N);
deltaBIC = DOGBIC-GaborBIC;
plot_counter = kk+1;
figure(plot_counter); set(gcf,'Name','model comparison');
subplot(221); plot(GaborBIC(DOidx),DOGBIC(DOidx),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[1 1 1]); hold on; 
plot(GaborBIC(LUMidx),DOGBIC(LUMidx),'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]); axis square;
set(gca,'Xlim',[-600 -100],'Ylim',[-600 -100],'TickDir','out'); line([-600 -100],[-600 -100]);xlabel('Gabor BIC'); ylabel('DOG BIC'); hold off;
subplot(222); histogram(deltaBIC(DOidx),bins); hold on; histogram(deltaBIC(LUMidx),bins); set(gca,'TickDir','out'); axis square; xlabel('delta BIC'); ylabel('count'); legend('DO','Lum'); hold off;
subplot(223); stem(deltaBIC); hold on;line([0 numel(Input_List)+1],[6 6]); set(gca,'TickDir','out'); axis square; ylabel('delta BIC'); hold off;
plot_counter = plot_counter+1;

%%



