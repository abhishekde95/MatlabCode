% In this script, I want to compare Gratings and WN data
% NLI vs. F1/F0
% Author: Abhishek De, 11/19

close all; clearvars;
[filename,spikeIdx] = fnamesFromTxt2('Gratings&WN.txt');
plot_counter = 1;
numcells = numel(filename);
Gratings_ID = 150; % pulled the info from the files
WN_ID = 100; % pulled the info from the files
channels = 3;
NPOINTS = 65536;
CHI2CRIT = 0.95; % For flaging a stixel as significant (adding gun and cone noise z-scores, squared)
spikename_options = ['sig001a'; 'sig001b'];
maxT = 9;
F1overF0 = [];
NLI = [];
totspikes = [];
ngammasteps = 2^16; % 65536
numpixels = [];
global maskidx spikeidx nstixperside seedidx nframesidx stimonidx muidxs sigmaidxs msperframe maxT yy
for ii = 1:numcells
    disp(ii);
    flag = 0;
    stro = []; ID = [];
    for jj = 1:size(filename{ii},2)
        try
            stro = [stro; nex2stro(findfile(char(filename{ii}(jj))))];
            ID = [ID; stro(end).sum.paradigmID];
            if stro(end).sum.paradigmID == WN_ID
                noisetypeidx = find(strcmp(stro(end).sum.trialFields(1,:),'noise_type'));
                L = stro(end).trial(:,noisetypeidx)==1;
                if ~any(L)
                    flag = 1;
                end
            end
        catch
            flag = 1;
            break;
        end
    end
    if flag
        continue;
    end
    
    WNidx = find(ismember(ID,WN_ID));
    Gratingsidx = find(ismember(ID,Gratings_ID));
    
    % First: WN Analysis
    if ~isempty(WNidx)
        WN = {};
        for kk = 1:numel(WNidx)
            if kk ==1
                WN = stro(WNidx(kk));
            else
                WN = strocat(WN, stro(WNidx(kk)));
            end
        end
        framerate = WN.sum.exptParams.framerate;
        nstixperside = WN.sum.exptParams.nstixperside;
        ntrials = length(WN.sum.absTrialNum);
        stimonidx = find(strcmp(WN.sum.trialFields(1,:),'stim_on'));
        stimoffidx = find(strcmp(WN.sum.trialFields(1,:),'all_off'));
        nframesidx = find(strcmp(WN.sum.trialFields(1,:),'num_frames'));
        noisetypeidx = find(strcmp(WN.sum.trialFields(1,:),'noise_type'));
        seedidx = strcmp(WN.sum.trialFields(1,:),'seed');
        sigmaidxs = strmatch('sigma',WN.sum.trialFields(1,:));
        hepidx = find(strcmp(WN.sum.rasterCells(1,:),'AD11'));
        vepidx = find(strcmp(WN.sum.rasterCells(1,:),'AD12'));
        maskidx = strcmp(WN.sum.rasterCells(1,:), 'subunit_mask');
        anlgStartTimeidx = find(strcmp(WN.sum.rasterCells(1,:),'anlgStartTime'));
        L = WN.trial(:,noisetypeidx)==1;
        maskidx = strcmp(WN.sum.rasterCells(1,:),'subunit_mask');
        basisvecidx = strcmp(WN.sum.rasterCells(1,:),'basis_vec');
        muidxs = [find(strcmp(WN.sum.trialFields(1,:),'mu1')), ...
            find(strcmp(WN.sum.trialFields(1,:),'mu2')), ...
            find(strcmp(WN.sum.trialFields(1,:),'mu3'))];
        msperframe = 1000/WN.sum.exptParams.framerate;
        spikeidx = spikeIdx(ii);
        spikename = spikename_options(spikeidx,:);
        xx = linspace(WN.sum.exptParams.gauss_locut/1000, WN.sum.exptParams.gauss_hicut/1000,ngammasteps); % xx represents the probabilities. For more info, have a look at the MATLAB 'norminv' function.
        yy = norminv(xx');
        
       
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
                
        % Calculating STA and STC for frames which triggered spikes
        out_gun = getWhtnsStats(WN,maxT,'STCOVmex',{nstixperside^2,3,maxT},spikename);
        STS_gun = out_gun{1}; STCross_gun = out_gun{2}; nspikes_gun = out_gun{3}; clear out_gun;
        STAs_gun = STS_gun/nspikes_gun;
        totspikes = [totspikes; nspikes_gun];
        
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
        tmp = STS_gun(:,id)*STS_gun(:,id)';
        STCs = (nspikes_gun.*reshape(STCross_gun(:,id),[300 300])-tmp)/(nspikes_gun*(nspikes_gun-1));
        
        
%         P = eye(size(STCs)) - STAs_gun(:,id)*inv(STAs_gun(:,id)'*STAs_gun(:,id))*STAs_gun(:,id)'; % Subtracting such that the STA is orthogonal to PC
%         STCs = P*STCs*P';
%         [tmp,d] = eig(STCs);
%         [eig_PC,ind] = sort(diag(d)); % storing all the eigenvalues
%         
%         [~,tmp_percentileval,tmpzscore] = permutation_test(WN,[1 1; mask_changes]',plot_counter,eig_PC(end),0,id);
%         keyboard;
        % Calculating NLI (Horwitz et al; 2007)
        % Calculating maximally informative dimension
        try
            
            out_all = getWhtnsStats(WN,maxT-id,'STCOVmex',{nstixperside^2,3,1},spikename,[],1);
            STS_all = out_all{1}; STCross_all = out_all{2}; nspikes_all = out_all{3}; clear out_all;
            STAs_all = STS_all/nspikes_all;
            tmp_all = STS_all*STS_all';
            STCs_all = (nspikes_all.*reshape(STCross_all,[300 300])-tmp_all)/(nspikes_all*(nspikes_all-1));
            fact = getRFfromSTA(tmpSTA,1,0.0); % arbitrary threshold of 70%
%             fact = getRFfromSTA(tmpSTA,0);
            tmpnumpixels = sum(fact(:));
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
            tmpNLI = (STATS_quad(1)-STATS_lin(1))/STATS_both(1);
            
            tmp_vec_gun = filts;
            normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
            im = normfactor*tmp_vec_gun + 0.5;
            im = reshape(im,[10 10 3]);
            figure(plot_counter); subplot(8,8,ii); image(im); axis square; set(gca,'XTick',[],'YTick',[]);
            if tmpNLI>0
                set(gca,'XColor','r','YColor','r');
            end
           
            tmp_vec_gun = tmpSTA;
            normfactor = 0.5/(max(abs(tmp_vec_gun(:)))+0.01);
            im = normfactor*tmp_vec_gun + 0.5;
            im = reshape(im,[10 10 3]);
            figure(plot_counter+1); subplot(8,8,ii); image(im); axis square; set(gca,'XTick',[],'YTick',[]);
            if tmpNLI>0
                set(gca,'XColor','r','YColor','r');
            end
        catch
            tmpNLI = nan;
            tmpnumpixels = nan;
        end
        NLI = [NLI; tmpNLI]; 
        numpixels = [numpixels; tmpnumpixels];
    end
    
    % Second: Gratings Analysis
    if ~isempty(Gratingsidx)
        GT = {};
        for kk = 1:numel(Gratingsidx)
            if kk == 1
                GT = stro(Gratingsidx(kk));
            else
                GT = strocat(GT, stro(Gratingsidx(kk)));
            end
        end
        orients = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'orient'));
        sfs = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'sf'));
        diams = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'diam'));
        protocols = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'protocol'));
        framerate = GT.sum.exptParams.framerate;
        nframes = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'nframes'));
        stimon_t = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'stim_on'));
        stimoff_t= GT.trial(:,strcmp(GT.sum.trialFields(1,:),'stim_off')); 
        spikeidx = strcmp(GT.sum.rasterCells(1,:),spikename_options(spikeIdx(ii),:));
        spikerates = [];
        baselines = [];  baseline_t = 0.5;
        for i = 1:size(GT.trial,1)
            spiketimes = GT.ras{i,spikeidx};
            nspikes = sum(spiketimes > stimon_t(i) & spiketimes < stimoff_t(i));
            spikerates = [spikerates; nspikes./(stimoff_t(i)-stimon_t(i))];
            nspikes = sum(spiketimes > stimon_t(i)-baseline_t & spiketimes < stimon_t(i));
            baselines = [baselines; nspikes./baseline_t];
        end
        binwidth = .0025;
        % This smoothing filter may be unnecessary
        %filterwidth = .000012; % s
        %smoothingfilter = normpdf(linspace(-3,3,6*filterwidth./binwidth),0,1);
        %smoothingfilter = smoothingfilter./sum(smoothingfilter);
        
        Lprotocol = protocols == 4;
        start_t = stimon_t(Lprotocol);
        end_t = stimoff_t(Lprotocol);
        minstimdur = min(end_t-start_t);
        if (any(protocols == 4))
            Lcc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'lcont'));
            Mcc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'mcont'));
            Scc = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'scont'));
            
            colordirections = [Lcc(Lprotocol) Mcc(Lprotocol) Scc(Lprotocol)];
            uniquecolordirs = unique(colordirections,'rows');
            bins = [0:binwidth:minstimdur];
            data = [];
            for i = 1:size(uniquecolordirs,1)
                L = Lprotocol &...
                    Lcc == uniquecolordirs(i,1) & ...
                    Mcc == uniquecolordirs(i,2) &...
                    Scc == uniquecolordirs(i,3);
                start_t = stimon_t(L);
                end_t = stimoff_t(L);
                spiketimes = GT.ras(L,1);
                totalspikes = [];
                for j = 1:sum(L)
                    totalspikes = [totalspikes; spiketimes{j}-start_t(j)];
                end
                totalspikes(totalspikes > minstimdur) = [];
                psth = histc(totalspikes',bins)/binwidth/sum(L);
                
                %   psth = conv(psth,smoothingfilter,'full');
                %   psth(1:length(smoothingfilter)/2-1) = [];
                %   psth(end:-1:end-length(smoothingfilter)/2) = [];
                
                k = unique(GT.trial(L,strcmp(GT.sum.trialFields(1,:), 'tf')));
                if (length(k) > 1)
                    error('Too many temporal frequencies!');
                end
                basis0 = ones(length(psth),1);
                basis1 = exp(-2*pi*sqrt(-1)*k*[0:length(psth)-1]/length(psth))';
                F0 = abs((psth-mean(baselines))*basis0);
                F1 = 2*abs((psth-mean(baselines))*basis1);
                % Not sure what the logic for the the "2*" in the line above
                % but it makes a half-wave rectified sinewave have a modulation
                % ratio of pi/2, which is what it's supposed to have.
                %   F0 = abs(sum(psth));
                %   F1 = abs(sum((psth)*basis1));
                data = [data; F1 F0];
            end
            
            L = data(:,2)/length(psth) > 10;
            if (any(L))
                unweighted = mean(data(L,1)./data(L,2));
                weighted = sum(data(L,2).*(data(L,1)./data(L,2)))/sum(data(L,2));
                
            else
                weighted = nan;
            end
            F1overF0 = [F1overF0; weighted];
        end
    end
end

plot_counter  = plot_counter + 2;
 
%%
idx = ~(isnan(NLI) | isnan(F1overF0));
[r,p] = corr(NLI(idx),F1overF0(idx),'rows','complete','type','Spearman');
bins = logspace(2,5,21);
figure(plot_counter);
subplot(221); plot(NLI(idx),F1overF0(idx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[-1 1],'XTick',-1:1:1,'Ylim',[0 3],'YTick',0:1:3); xlabel('NLI (WN)'); ylabel('F1/F0 (Gratings)'); 
title(strcat('r=',num2str(r,2),',p=',num2str(p,3),',n=36')); 
subplot(222); histogram(totspikes,bins,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]); axis square; 
set(gca,'Tickdir','out','Xlim',[100 100000],'XScale','log'); xlabel('Total spikes in WN'); ylabel('count'); 
subplot(223); plot(totspikes(idx),F1overF0(idx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[100 100000],'XScale','log','Ylim',[0 3],'YTick',0:1:3); xlabel('Total Spikes in WN'); ylabel('F1/F0 (Gratings)');
subplot(224); plot(totspikes(idx),NLI(idx),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
axis square; set(gca,'Tickdir','out','Xlim',[100 100000],'XScale','log','Ylim',[-1 1],'YTick',-1:1:1); xlabel('Total Spikes in WN'); ylabel('NLI (WN)');
plot_counter = plot_counter + 1;

