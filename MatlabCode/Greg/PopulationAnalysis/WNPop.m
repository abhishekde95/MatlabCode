% Population analysis of white noise data
% Contents
% 1) Comparing cone weight estimates from gun and cone noise.
%
% 2) Trying automatically go through a whole directory of nex files
% and pull out the white noise files.  This might want to go into
% CellScreen.m someday.
%
% 3) Looking at PSTH synced on fixational saccades
%
% 4) Looking at static nonlinearities (based on STAs) for groups of cells.
%
% 5) Looking at STAs omitting spikes, either around the time of 
% a fixational eye movement or randomly distributed.
%
% 6) Plots STAs and Gabor function fits for a group of cells.
%
% 6.1) Plotting STAs in space-space and space-time.  One way to look for
% direction selective cells.
%
% 7) Just like section 6 but also calculates saccade-triggered PSTHs.  This
% is for asking whether there is a systematic relationship between the
% modulations in activity induced by microsaccades and the SF or
% orientation tuning of the cell (Sedna, at least, tends to make eye
% movements along a particular axis).
%
% 7.1) Does the STA latency correlate with the latency of the dip or the
% peak in the saccade-triggered PSTH?  (Must run code in section 7 above
% before running this code).
%
% 7.2) Does anything about the saccade-triggered PSTH correlate with the
% cone weights of the cell?  (Specifically, do L-M neurons show less
% suppression?)
%
% 8) For each cell, calculating the saccade-triggered PSTH for
% microsaccades roughly aligned with and orthogonal to the cell's preferred
% orientation (obtained from a Gabor fit).
%
% 9) Comparing stimulus tuning measured with white noise and gratings.
% This has been superceded by CompareSTAandGrating.m.
%
% 10) STAs and PC1 for a group of cells.
%
% 10.1) STAs (across time) for a group of cells.
%
% 10.2) Comparing cone weights at the single hottest stixel vs. cone
% weights computed from all of the stixels in the RF.
%
% 11) Looking at the effects of microsaccades on the firing rate of the
% cell as a function of microsaccade amnplitude.  Just looking at
% microsaccades above and below the median amplitude.
%
% 11.1) Microsaccade effects on firing rate as a function of saccade
% amplitude.
%
% 12) Comparing neuronal visual response latency between STA and replays
%
% 13) Looking at all the (temporalchromatic STAs) for a group of LGN
% neurons (M,P,K).
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 1
% Estimates of cone weights from gun and cone noise.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datapath = nexfilepath;
[filenames,spikecds] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\WhiteNoise\ConeVGun2.txt');
maxT = 8;
CHI2CRIT = .95; % For flaging a stixel as significant (adding gun and cone noise z-scores, squared)
data = [];
for i = 1:size(filenames,1)
    stro = {};
    for j = 1:length(filenames{i})
        stro{j} = nex2stro(char(findfile(filenames{i}(j))));
    end
    stro = strocat(stro);
 
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    sigmaidxs = strmatch('sigma',stro.sum.trialFields(1,:));
    spikename = ['sig001',char(double('a'+spikecds(i)-1))];
    nstixperside = stro.sum.exptParams.nstixperside;
    
    % Reconstructing the M matrix and gamma table
    fundamentals = stro.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = stro.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    gammaTable = stro.sum.exptParams.gamma_table;
    gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
    
    % Getting the background rgb/lms
    ridx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_r'));
    gidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_g'));
    bidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_b'));
    bkgndRGB = [mode(stro.trial(:,ridx)), mode(stro.trial(:,gidx)), mode(stro.trial(:,bidx))];
    bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
    bkgndlms = M*bkgndrgb;
    
    gunsigmas = unique(stro.trial(stro.trial(:,noisetypeidx) == 1, sigmaidxs),'rows')/1000;  % normalized gun intensity units
    conesigmas = unique(stro.trial(stro.trial(:,noisetypeidx) == 2, sigmaidxs),'rows')/1000;  % cone excitation units
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Getting STAs
    for noisetype = 1:2
        tmpstro = stro;
        L = stro.trial(:,noisetypeidx) == noisetype;
        tmpstro.ras(~L,:) = [];
        tmpstro.trial(~L,:) = [];
        out = getWhtnsStats(tmpstro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
        if (noisetype == 1)
            STAs_gun = out{1};
            STCs_gun = out{2};
            nspikes_gun = out{3};
        elseif (noisetype == 2)
            STAs_cone = out{1};
            STCs_cone = out{2};
            nspikes_cone = out{3};
        end
    end
    clear tmpstro;
    
    s_gun = std(STAs_gun(:,1));
    STAs_gun_z = STAs_gun./s_gun;
    s_cone = std(STAs_cone(:,1));
    STAs_cone_z = STAs_cone./s_cone;
    
    % Spatial map
    grandz = zeros([nstixperside nstixperside]);
    maxzs = [];
    for i = 1:maxT
        tmp_gun = reshape(STAs_gun_z(:,i),[nstixperside nstixperside 3]);
        tmp_cone = reshape(STAs_cone_z(:,i),[nstixperside nstixperside 3]);
        grandz = grandz+sum(tmp_gun.^2,3)+sum(tmp_cone.^2,3);
        maxzs(i) = sum(sum(sum(tmp_gun.^2,3)+sum(tmp_cone.^2,3)));
    end
    crit = chi2inv(CHI2CRIT,6*maxT); % 6 = 3 color channels * 2 stimulus sets
    L = grandz > crit;
    % Now get largest contiguous block
    [tmpi,tmpj] = ind2sub(size(L),find(L));
    ij = [tmpi,tmpj];
    T = clusterdata(ij,'linkage','single','distance','euclidean','criterion','distance','cutoff',sqrt(2));
    
    clusternmembers = [];
    for k =1:max(T)
        clusternmembers(k) = sum(T == k);
    end
    dominantcluster = find(clusternmembers == max(clusternmembers));
    
    clustermat = zeros(nstixperside, nstixperside);
    clustermat(sub2ind(size(clustermat),ij(T==dominantcluster,1),ij(T==dominantcluster,2))) = 1;
    
    % Then get convex hull
    dominantclusteridxs = ij(T==dominantcluster,:);    
    K=convhull(dominantclusteridxs(:,1), dominantclusteridxs(:,2));
    tmp = dominantclusteridxs(K,:);
    [x,y] = meshgrid(1:nstixperside,1:nstixperside);
    inRF = reshape(inpolygon(x(:),y(:),dominantclusteridxs(K,2),dominantclusteridxs(K,1)),[nstixperside, nstixperside]);
    
    % Comparing estimates of cone weights
    whichpix = find(inRF);
    tmpSTA = reshape(STAs_gun, [nstixperside^2 3 maxT]);
    tmpSTA = permute(tmpSTA, [2 1 3]);
    tmpSTA = tmpSTA(:,whichpix,:);
    STAgunmat = reshape(tmpSTA,[3 length(whichpix)*maxT]);
    tmpSTA = reshape(STAs_cone, [nstixperside^2 3 maxT]);
    tmpSTA = permute(tmpSTA, [2 1 3]);
    tmpSTA = tmpSTA(:,whichpix,:);
    STAconemat = reshape(tmpSTA,[3 length(whichpix)*maxT]);
    conesigmas = stro.trial(stro.trial(:,noisetypeidx) == 2, sigmaidxs)/1000;  % cone excitation units
    
    % First, the gun noise STA
    [u,~,v] = svd(STAgunmat);
    if (sum(v(:,1)) < 0)
        u = -u;
    end
    coneweights = inv(M')*u(:,1);
    %coneweights = coneweights.*bkgndlms;   % Cone specific adaptation (on the weights?)
    %coneweights = coneweights.*(M*[1 1 1]');
    coneweights_gun = coneweights./sum(abs(coneweights));
    
    % Cone weights from cone noise expressed in contrast units
    %[u,s,v] = svd(STAconemat(:,whichpix)./repmat(lmscontrasts',[1 length(whichpix)]));%
    [u,s,v] = svd(STAconemat);
    if (coneweights_gun'*u(:,1)< 0)
        u = -u;
    end
    u(:,1) = u(:,1)./mode(abs(conesigmas))';  % Think about this line
    coneweights_cone = u(:,1)./sum(abs(u(:,1)));
    figure; axes; hold on;
    bar([coneweights_cone coneweights_gun]);
    legend({'cone','gun'});
    set(gca,'XTick',[1 2 3],'XTickLabel',{'L','M','S'});
    data = [data; coneweights_cone' coneweights_gun'];
end

% 3D cone weight plot
figure; axes; hold on;
for i = 1:size(data,1)
   plot3(data(i,1),data(i,2),data(i,3),'ko','MarkerFaceColor','black');
   plot3(data(i,4),data(i,5),data(i,6),'ro','MarkerFaceColor','red');
   plot3(data(i,[1,4]),data(i,[2,5]),data(i,[3,6]),'k-')
end
axis equal vis3d;

% Normalized cone weight plot
figure; axes; hold on;
plot(data(:,[1,4])',data(:,[2,5])','k-');
for S = 0:1
    if S == 0
        L = data(:,3) > 0;
        plot(data(L,1),data(L,2),'ks','MarkerFaceColor','black');
         L = data(:,6) > 0;
        plot(data(L,4),data(L,5),'ro','MarkerFaceColor','red');
    else
        L = data(:,3) < 0;
        plot(data(L,1),data(L,2),'ks','MarkerFaceColor','white');
        L = data(:,6) < 0;
        plot(data(L,4),data(L,5),'ro','MarkerFaceColor','white');
    end
end
plot([-1 0 1 0 -1],[0 -1 0 1 0],'k-');
axis square;
xlabel('Normalized L-cone weight');
ylabel('Normalized M-cone weight');

% % Differences in cone weights
% Ldiff = data(:,1)-data(:,4);
% Mdiff = data(:,2)-data(:,5);
% Sdiff = data(:,3)-data(:,6);
% 
% figure;
% subplot(2,2,1); hist(Ldiff);
% subplot(2,2,2); hist(Mdiff);
% subplot(2,2,3); hist(Sdiff);

figure; axes; hold on;
plot(data(:,1),data(:,4),'ro','MarkerFaceColor','red');
plot(data(:,2),data(:,5),'go','MarkerFaceColor','green');
plot(data(:,3),data(:,6),'bo','MarkerFaceColor','blue');
plot([-1 1],[-1 1],'k-');
set(gca,'Xlim',[-1 1],'Ylim',[-1 1]);
axis square;
xlabel('Normalized cone weights (cone noise)');
ylabel('Normalized cone weights (gun noise)');



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 2
% Starting at some point in the directory tree, going through it
% and any child directories (one level down only), and pulling out
% all the nex files.  For Nex files that are paradigm 100 (white noise)
% and have 200 or more spikes (in gun noise trials), we calculate the STA and store
% the "hottest" frame".
%
% Assumes a 10 x 10 grid Arg!
%%%%%%%%%%%%%%%%%%%%%%%%%%
PARADIGMID = 100;
PLOTPC1 = false; 
% datapath = 'N:\NexFiles\Greg';
% eval(['cd (''',datapath,''')']);
% % First making a gigantic dirstruct
% dirstruct = dir;
% for i = 1:size(dirstruct,1)
%     if (dirstruct(i).isdir && ~strncmp(dirstruct(i).name, '.',1))
%         eval(['cd (''',dirstruct(i).name,''')']);
%         dirstruct = [dirstruct; dir];
%         eval(['cd ..']);
%     end
% end

conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
fnamesfromDB = fetch(conn, 'SELECT * FROM WhiteNoise_old');
fnamesfromDB = fnamesfromDB(:,2);
close(conn)

STAdata = {};
PC1data = {};
filenames = {};
maxT = 9;
opts.disp = 0;
opts.isreal = 1;
opts.issym = 1;
counter = 1;
for i = 1:size(fnamesfromDB,1)
    if (length(fnamesfromDB) < 4)
        continue;
    end
    if (all(fnamesfromDB{i}(end-3:end) == '.nex'))
        if (getparadigmID(findfile(fnamesfromDB{i}))~=PARADIGMID)
            continue;
        end
        try
        stro = nex2stro(findfile(fnamesfromDB{i}));
        catch
            disp(['Problem with ',dirstruct(i).name]);
            continue
        end
        if (~isfield(stro,'sum'))  % If there's a error
            continue
        end

        % This bit below should be modular
        for spikename = {'sig001a','sig001b'}
            if (~strcmp(spikename,stro.sum.rasterCells))
                continue;
            end
            noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
            nstixperside = stro.sum.exptParams.nstixperside;
            tmpstro = stro;
            L = stro.trial(:,noisetypeidx) == 1;  % Only looking at gun noise for now
            tmpstro.ras(~L,:) = [];
            tmpstro.trial(~L,:) = [];
            out = getWhtnsStats(tmpstro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, char(spikename));
            STS = out{1};
            STCross = out{2};
            n = out{3};
            if (n > 200)
                energy = sum(out{1}.^2);
                whichframe = find(energy == max(energy));
                STA = STS(:,whichframe)./n;
                STAdata{counter} = reshape(STA,[nstixperside nstixperside 3]);
                
                % GDLH
                tmp = STA*STA';
                STCframe = (n.*STCross(:,whichframe)-tmp(:))/(n*(n-1));

                STC = reshape(STCframe,3*nstixperside^2, 3*nstixperside^2);
                P = eye(size(STC))-STA*inv(STA'*STA)*STA';
                [v,d] = eigs(P*STC*P',1,'LM',opts);
                PC1data{counter} = reshape(v,[nstixperside nstixperside 3]);;
                slashpos = strfind(char(stro.sum.fileName),filesep);
                filenames{counter} = [stro.sum.fileName(slashpos(end)+1:length(stro.sum.fileName)-4),'   ',num2str(find(strcmp(spikename, {'sig001a','sig001b'})))];
                counter=counter+1
            end
        end
    end
end
% Getting filenames
for i = 1:length(filenames)
    disp(char(filenames{i}));
end

% Plotting
%ncells = find(isnan(STAdata(1,:)),1)-1;
ncells = length(STAdata);
NCELLSPERFIG = 64;
for j = 1:ceil(ncells/NCELLSPERFIG)  % loop over figures
    for k = 1:2 % STA, PC1
        figure;
        for i = 1:NCELLSPERFIG
            idx = (j-1)*NCELLSPERFIG+i;
            if (idx > ncells);
                break
            end
            
            subplot(ceil(sqrt(NCELLSPERFIG)),ceil(sqrt(NCELLSPERFIG)),i)
            
            if (k == 1)
                normfactor = .5./max(abs(STAdata{idx}(:)));
                im = normfactor*STAdata{idx};
            else
                if PLOTPC1
                    normfactor = .5./max(abs(PC1data{idx}(:)));
                    im = normfactor*PC1data{idx};
                else
                    normfactor = .5./max(abs(STAdata{idx}(:)));
                    im = normfactor*STAdata{idx};
                    im = cat(3,im,zeros(size(im,1),size(im,2)));
                    im = reshape(im,size(im,1),4*size(im,2));
                    im = [im(:,1:size(im,2)/2); im(:,size(im,2)/2+1:end)];
                end
            end
            if ndims(im) == 3
                h = image(im+.5);
            else
                h = image(64*(im+.5));
                colormap(gray(64));
            end
            if ~PLOTPC1 & k == 2
                hold on;
                plot([0,size(im,1)],[size(im,1)/2 size(im,1)/2]+.5,'k-');
                plot([size(im,1)/2 size(im,1)/2]+.5,[0,size(im,1)],'k-');
            end
            set(h,'ButtonDownFcn',['disp(''',char(filenames{idx}),''')']);
            axis square; set(gca,'XTick',[],'YTick',[]);
        end
    end
end


%%
% Section 3

% Looking at PSTHs synced to fixational saccades.
filelist = 'LGN.txt';
[filenames,spikenums] = fnamesFromTxt2(fullfile(nexfilepath, 'nexfilelists','Greg','WhiteNoise',filelist));

BINWIDTH = .01;
timebins = [-.25:BINWIDTH:.5];
secstoskip = 0;
data = nan*ones(length(filenames), length(timebins));
for i = 1:size(filenames,1)
    stro = {};
    for f = 1:size(filenames{i},1)
        tmpstro = nex2stro(findfile(char(filenames{i}(f))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
    
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
    nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
    
    if (isempty(stro.sum.analog.sigid))
        continue;
    end
    sacstats = getSacData(stro);
    close;
    L = logical(stro.trial(:,noisetypeidx)~= 0);

    PSTH = zeros(1,length(timebins));
    for j = find(L')
        stimon_t = stro.trial(j,stimonidx);
        numframes = stro.trial(j,nframesidx);
        stimoff_t = stimon_t+numframes/stro.sum.exptParams.framerate;
        st = sacstats.starttimes{j};
        Lsac = (st > stimon_t+secstoskip) & (st < stimoff_t-timebins(end)); %Don't include msacs at end of trial
        if any(Lsac)
            Lsac(sacstats.amplitudes{j}(Lsac) > 2) = [];
            spiketimes = [];
            for k = find(Lsac')
                tmp = stro.ras{j,spikenums(i)}-sacstats.starttimes{j}(k);
                spiketimes = [spiketimes; tmp((tmp > timebins(1)-BINWIDTH/2) & (tmp < timebins(end)+BINWIDTH/2))];
            end
            [n,x] = hist(spiketimes, timebins);
            PSTH = PSTH + n./(BINWIDTH*sum(Lsac));
        end
    end
   % close all;
    data(i,:) = PSTH./sum(L);
end

figure;
subplot(2,1,1);
peakfr = max(data,[],2);
imagesc(data./repmat(peakfr,1,length(timebins)));
set(gca,'XTick',find(rem(timebins,.1) == 0));
set(gca,'XTickLabel',timebins(get(gca,'Xtick')));
title([filelist,' ',date]);
subplot(2,1,2);
%errorbar(timebins,nanmean(data),nanstd(data)./sqrt(sum(~isnan(data(:,1)))),'k-')
plot(timebins, nanmean(data),'k-');
ylabel('sp/sec');
xlabel('time (s)');

%%
% Section 4
% Looking at static nonlinearities for a bunch of cells.  This is to see
% whether there is a systematic difference in the shape of this
% nonlinearity with color direction or preferred spatial frequency (we
% might expect a 1/f spectrum to lead to similar levels of activation
% across neurons).
datapath = 'N:\NexFiles\Kali';
[filenames, spikenums] = fnamesFromTxt2('N:\NexFiles\Kali\Lum.txt');
niters = 10;
CVprop = .7;
nbins = 7;
mns = []; stds = [];
for i = 1:size(filenames,1)
    stro = {};
    for f = 1:size(filenames{i},1)
        tmpstro = nex2stro(findfile(char(filenames{i}(f))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
    nstixperside = stro.sum.exptParams.nstixperside;
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
    Lgaussnoise = stro.trial(:,noisetypeidx) == 1;
    binlim = 2*mode(stro.trial(Lgaussnoise,9))/1000;
    ntrials = size(stro.trial,1);
    ntrialsSTA = round(ntrials*CVprop);
    frs = [];
    for iter = 1:niters
        trialidxs = randperm(ntrials);
        tmpstro = stro;
        tmpstro.trial = stro.trial(trialidxs(1:ntrialsSTA),:);
        tmpstro.ras = stro.ras(trialidxs(1:ntrialsSTA),:);
        out = getWhtnsStats(tmpstro, maxT, 'STCOVmex', {nstixperside^2, 3, maxT});
        STAs = out{1};
        
        energy = sum(STAs.^2);
        whichframe = find(energy == max(energy));
        
        basis = STAs(:,whichframe);
        basis = basis./sqrt(sum(basis.^2));

        tmpstro = stro;
        tmpstro.trial = stro.trial(trialidxs(ntrialsSTA+1:end),:);
        tmpstro.ras = stro.ras(trialidxs(ntrialsSTA+1:end),:);
        nframestotal = sum(stro.trial(:,nframesidx));
        initargs = {basis, whichframe-1, nframestotal, [nstixperside^2 3 1]};
        out = getWhtnsStats(tmpstro,whichframe,'STPROJmod',initargs);
        projs = out{1};
        Lspike = out{2};
        
        x = linspace(-binlim, binlim, nbins+2);
        Na = hist(projs(:,1), x);
        Ns = zeros(size(Na));
        for j = 1:max(Lspike)
            Ns = Ns + hist(projs(Lspike >=j,1), x);
        end
        frs = [frs; Ns(2:end-1)./Na(2:end-1)]
    end
    mns = [mns; nanmean(frs)];
    stds = [stds; nanstd(frs)]; 
end

data = [];
normfactor = max(range(mns'))
for i = 1:size(filenames,1)
    subplot(ceil(sqrt(size(filenames,1))),ceil(sqrt(size(filenames,1))),i);
    plot(x(2:end-1),mns(i,:),'b.-');
    set(gca,'YLim',[0 max(mns(i,:))*1.1]);
    set(gca,'XTick',[],'YTick',[0 max(mns(i,:))*1.1]);
    set(gca,'Color',repmat(range(mns(i,:))/normfactor,1 ,3));
    data = [data; range(mns(i,:))];
end

%%
% Section 5
% Looking at STAs from little pieces of trials, either right before 
% a fixational eye movement or right afterward.  This will answer the
% question whether the suppression/rebound observed following a fixational 
% saccade actually changes the signal to noise ratio.

[filenames, spikenums] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\Lum.txt');
offset_t = 0;
delta_t = .2;
maxT = 9;
data = [];
for i = 1:size(filenames,1)
    stro = {};
    for f = 1:size(filenames{i},1)
        tmpstro = nex2stro(findfile(char(filenames{i}(f))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
    if (isempty(stro.sum.analog.sigid))
        continue;
    end
    
    nstixperside = stro.sum.exptParams.nstixperside;
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
    stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));

    Lnoise = stro.trial(:,noisetypeidx) == 1;
    stro.trial(~Lnoise,:) = [];
    stro.ras(~Lnoise,:) = [];
    stimon_t = stro.trial(:,stimonidx);
    numframes = stro.trial(:,nframesidx);
    stimoff_t = stimon_t+numframes/stro.sum.exptParams.framerate;
    sacstats = getSacData(stro);
    close;
    spikes = stro.ras(:,spikenums(i));
    spikename = stro.sum.rasterCells{spikenums(i)};
    
    % First just getting the STA and nspikes (STAs0)
    out = getWhtnsStats(stro, maxT, 'STCOVmex', {nstixperside^2, 3, maxT},spikename);
    STAs0 = out{1}./out{3};
    nspikes0 = out{3};
    energy = sum(STAs0.^2);
    whichframe = find(energy == max(energy));
    
    %peakframeSTA = STAs0(:,whichframe);
    %peakframeenergy = squeeze(sum(reshape(peakframeSTA,[10 10 3]).^2,3));
    %centerpix = find(peakframeenergy == max(abs(peakframeenergy(:))));
 
    % Now omitting spikes around the time of saccades (STAs1)
    oldspikes = spikes;
    for j = 1:size(stro.trial,1)
        for t = sacstats.starttimes{j}'
            L = spikes{j} > t+offset_t & spikes{j} < t+offset_t+delta_t;
            spikes{j}(L) = [];
        end
    end
    stro.ras(:,spikenums(i)) = spikes;
    out = getWhtnsStats(stro, maxT, 'STCOVmex', {nstixperside^2, 3, maxT},spikename);
    STAs1 = out{1}./out{3};
    nspikes1 = out{3};
  
    % Now omitting that same number of spikes from random times in the data
    % (STAs2)
    spikes = oldspikes;
    nspikes2omit = nspikes0-nspikes1;
    sizes = [];
    for j = 1:length(oldspikes)
        sizes = [sizes; size(spikes{j},1)];
    end
    cumsumsizes = cumsum(sizes);
    while nspikes2omit > 0
        candidate = unidrnd(sum(sizes));
        whichtrial = find(cumsumsizes > candidate,1,'first');
        if isempty(whichtrial)
            whichtrial = length(cumsumsizes);
        end
        spikeidx = cumsumsizes(whichtrial)-candidate+1;
        if (length(spikes{whichtrial}) < spikeidx)
            continue
        end
        tmp = spikes{whichtrial}(spikeidx);
        if (tmp > stimon_t(whichtrial)+maxT/stro.sum.exptParams.framerate & tmp < stimoff_t(whichtrial))
            spikes{whichtrial}(spikeidx) = [];
            nspikes2omit = nspikes2omit - 1;
        end
    end
    stro.ras(:,spikenums(i)) = spikes;
    out = getWhtnsStats(stro, maxT, 'STCOVmex', {nstixperside^2, 3, maxT},spikename);
    STAs2 = out{1}./out{3};
    nspikes2 = out{3};
    
    % Now only *keeping* spikes around the time of saccades (STAs3)
    spikes = oldspikes;
    for j = 1:size(stro.trial,1)
        L = logical(zeros(length(spikes{j}),1));
        for t = sacstats.starttimes{j}'
            L(spikes{j} > t+offset_t & spikes{j} < t+offset_t+delta_t) = 1;
        end
        spikes{j} = spikes{j}(L);
    end
    stro.ras(:,spikenums(i)) = spikes;
    out = getWhtnsStats(stro, maxT, 'STCOVmex', {nstixperside^2, 3, maxT},spikename);
    STAs3 = out{1}./out{3};
    nspikes3 = out{3};
    
    lim = max(abs([STAs1(:);STAs3(:)]));
    figure;
    subplot(3,1,1);
%    plot(STAs1(centerpix+[0 100 200],:)','-')
%    plot(STAs1(:)); set(gca,'YLim',[-lim lim]);
    plot(STAs1(:,whichframe)); set(gca,'YLim',[-lim lim]);
    subplot(3,1,2)
%    plot(STAs3(centerpix+[0 100 200],:)','-')
%    plot(STAs3(:));  set(gca,'YLim',[-lim lim]);
   plot(STAs3(:,whichframe)); set(gca,'YLim',[-lim lim]);

    subplot(3,1,3)
%    plot(STAs3(centerpix+[0 100 200],:)'-STAs1(centerpix+[0 100 200],:)','-')
%    plot(STAs1(:)-STAs3(:));  set(gca,'YLim',[-lim lim]);
    plot(STAs1(:,whichframe)-STAs3(:,whichframe));  set(gca,'YLim',[-lim lim]);
    
    SNR0 = norm(STAs0(:,whichframe))./norm(STAs0(:,1));
    SNR1 = norm(STAs1(:,whichframe))./norm(STAs1(:,1));
    SNR2 = norm(STAs2(:,whichframe))./norm(STAs2(:,1));
    data = [data; SNR0 SNR1 SNR2]
end
% Key: 
% STAs0 uses all the data
% STAs1 omits spikes near the time of a saccade
% STAs2 omits the same number of spikes as STAs1 but randomly in time
% STAs3 only uses spikes near the time of the saccade



figure;
subplot(2,1,1);
hist(data(:,2)-data(:,1),20)
[h,p] = ttest(data(:,2)-data(:,1)) % mean is -: SNR is better with all the spikes
subplot(2,1,2);
hist(data(:,2)-data(:,3),20) % mean is +: We're better off losing spikes near saccades
[h,p] = ttest(data(:,2)-data(:,3))


%%
% Section 6
%
% Plots STAs and gabor fits for a bunch of cells. 

ONEFIGPERCELL = 0;
%[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\WhiteNoise\LGN.txt');
[fnames, spikeIdx] = fnamesFromTxt('WhiteNoiseLGN_forIS')
nExpts = size(fnames,1);
sigmas = nan(1,nExpts);
lambdas = nan(1,nExpts);
gammas = nan(1,nExpts);
thetas = nan(1,nExpts);
success = nan(1, nExpts);
if (~ONEFIGPERCELL)
    h1 = figure;
    h2 = figure;
end
for a = 1:nExpts;
    stro = {};
    for i = 1:size(fnames{a},2)
        tmpstro = nex2stro(findfile(char(fnames{a}(i))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
    maxT = 8;
    framerate = stro.sum.exptParams.framerate;
    nstixperside = stro.sum.exptParams.nstixperside;
    ntrials = length(stro.sum.absTrialNum);
    stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoffidx = find(strcmp(stro.sum.trialFields(1,:),'all_off'));
    nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    sigmaidxs = strmatch('sigma',stro.sum.trialFields(1,:));

    mondist = 100; % cm
    screenwidth = 36; %cm
    screenwidthinpix = 1024; % singlewide pixels
    pixpercm = screenwidthinpix/screenwidth;
    cmperdeg = screenwidth/(2*atan2(screenwidth/2, mondist)*180/pi);
    pixperdeg = pixpercm*cmperdeg;
    stixperdeg = pixperdeg/(2*stro.sum.exptParams.npixperstix);

    L = stro.trial(:,noisetypeidx) == 1;
    stro.ras(~L,:) = [];
    stro.trial(~L,:) = [];

    out = getWhtnsStats(stro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, stro.sum.rasterCells{spikeIdx(a)});
    tmpstro = [];
    STAs = out{1};
  
    energy = sum(STAs.^2);
    whichframe = find(energy == max(energy));

    % fitting the gabor
    im = reshape(STAs(:,whichframe),nstixperside*nstixperside,3);
    [u,s,v] = svd(im');
    im = reshape(v(:,1),[nstixperside nstixperside]);
    struct = gaborfit(im);

    sigmas(a) = struct.sigma./stixperdeg;
    lambdas(a) = struct.lambda./stixperdeg;
    gammas(a) = struct.gamma;
    thetas(a) = struct.theta;
    success(a) = struct.exitflag;

    % Plotting STAs
    maxSTA = max(STAs(:,whichframe));
    minSTA = min(STAs(:,whichframe));

    potentialnormfactors = [(.5-eps)./(maxSTA-.5); (-.5+eps)./(minSTA-.5)];
    % 'eps' in above line is a kludge that is required for avoiding
    % out of bounds errors.
    potentialnormfactors(potentialnormfactors < 0) = []; % if min > mu or max < mu
    normfactor = min(potentialnormfactors);
    muvect = reshape(repmat([.5 .5 .5],nstixperside^2,1),nstixperside^2*3,1);

    % Plotting
    if (ONEFIGPERCELL)
        figure;
        subplot(2,1,1);
    else
        figure(h1);
        subplot(ceil(sqrt(nExpts)),ceil(sqrt(nExpts)),a);
    end
    STA = normfactor*(STAs(:,whichframe)-muvect)+muvect;
    STA = reshape(STA,[nstixperside nstixperside 3]);
    image(STA);
    set(gca,'XTick',[],'YTick',[]); axis square;

   % if (sum(v(:,1)) < 0)
   %     gaborrgb = -u(:,1);
   % else
   %     gaborrgb = u(:,1);
   % end
    gaborrgb = u(:,1);
    gaborimage = DrawGaborEdge([.5 .5 .5], gaborrgb/2, [0 0 0], struct.theta, struct.lambda, struct.sigma, struct.gamma, struct.phi, struct.xoffset, struct.yoffset, 0, 0, .99, 10);
    if (ONEFIGPERCELL)
        subplot(2,1,1);
    else
        figure(h2);
        subplot(ceil(sqrt(nExpts)),ceil(sqrt(nExpts)),a);
    end
    image(gaborimage);
    set(gca,'XTick',[],'YTick',[]); axis square;
    if (ONEFIGPERCELL)
        title(sprintf('file: %s \n success: %d', fnames(a,:), success(a)));
    end
    drawnow;
end

%%
% Section 6.1
%
% Plots STAs and Gabor fits in space-space and space-time

[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\Lum.txt');
nExpts = size(fnames,1);
sigma_s = nan(1,nExpts);
lambda_s = nan(1,nExpts);
gamma_s = nan(1,nExpts);
theta_s = nan(1,nExpts);
sigma_t = nan(1,nExpts);
lambda_t = nan(1,nExpts);
gamma_t = nan(1,nExpts);
theta_t = nan(1,nExpts);

success = nan(1, nExpts);
h1 = figure;
h2 = figure;
h3 = figure;
h4 = figure;

for a = 1:nExpts;
    stro = {};
    for i = 1:size(fnames{a},2)
        tmpstro = nex2stro(findfile(char(fnames{a}(i))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
    maxT = 8;
    framerate = stro.sum.exptParams.framerate;
    nstixperside = stro.sum.exptParams.nstixperside;
    ntrials = length(stro.sum.absTrialNum);
    stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoffidx = find(strcmp(stro.sum.trialFields(1,:),'all_off'));
    nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    sigmaidxs = strmatch('sigma',stro.sum.trialFields(1,:));

    mondist = 100; % cm
    screenwidth = 36; %cm
    screenwidthinpix = 1024; % singlewide pixels
    pixpercm = screenwidthinpix/screenwidth;
    cmperdeg = screenwidth/(2*atan2(screenwidth/2, mondist)*180/pi);
    pixperdeg = pixpercm*cmperdeg;
    stixperdeg = pixperdeg/(2*stro.sum.exptParams.npixperstix);

    L = stro.trial(:,noisetypeidx) == 1;
    stro.ras(~L,:) = [];
    stro.trial(~L,:) = [];

    out = getWhtnsStats(stro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, stro.sum.rasterCells{spikeIdx(a)});
    tmpstro = [];
    STAs = out{1};
    energy = sum(STAs.^2);
    whichframe = find(energy == max(energy));
 
    % fitting the gabor
    im = reshape(STAs(:,whichframe),nstixperside*nstixperside,3);
    [u,s,v] = svd(im');
    im = reshape(v(:,1),[nstixperside nstixperside]);
    struct = gaborfit(im);

    sigma_s(a) = struct.sigma./stixperdeg;
    lambda_s(a) = struct.lambda./stixperdeg;
    gamma_s(a) = struct.gamma;
    theta_s(a) = struct.theta;
    success_s(a) = struct.exitflag;

    STAs = reshape(STAs,[nstixperside nstixperside 3 maxT]);
    STA  = squeeze(STAs(:,:,:,whichframe));
    energy = sum(STA.^2,3);
    peakpixelidx = find(energy == max(energy(:)));
    [i,j] = ind2sub(size(energy),peakpixelidx);
    ang = mod(thetas(a),2*pi);
    if ((ang > pi/4) & (ang < 3*pi/4)) | ((ang > 5*pi/4) & (ang < 7*pi/4))
        % Use horizontal section
        STA_t = squeeze(STAs(i,:,:,:));
    else
        % Use vertical section
        STA_t = squeeze(STAs(:,j,:,:));
    end
    STA_t = permute(STA_t,[1 3 2]);
  
    % Plotting STAs
    maxSTA = max(STAs(:));
    minSTA = min(STAs(:));

    potentialnormfactors = [(.5-eps)./(maxSTA-.5); (-.5+eps)./(minSTA-.5)];
    % 'eps' in above line is a kludge that is required for avoiding
    % out of bounds errors.
    potentialnormfactors(potentialnormfactors < 0) = []; % if min > mu or max < mu
    normfactor = min(potentialnormfactors);

    % Plotting space-space STA
    figure(h1);
    subplot(ceil(sqrt(nExpts)),ceil(sqrt(nExpts)),a);
    STA = normfactor*(STAs(:,:,:,whichframe)-.5)+.5;
    image(STA);
    set(gca,'XTick',[],'YTick',[]); axis square;

    % Plotting space-time STA
    figure(h2);
    subplot(ceil(sqrt(nExpts)),ceil(sqrt(nExpts)),a);
    STA = normfactor*(STA_t-.5)+.5;
    image(STA);
    set(gca,'XTick',[],'YTick',[]); axis square;
    
    % Plotting the Gabor
    figure(h3);
    gaborrgb = u(:,1);
    gaborimage = DrawGaborEdge([.5 .5 .5], gaborrgb/2, [0 0 0], struct.theta, struct.lambda, struct.sigma, struct.gamma, struct.phi, struct.xoffset, struct.yoffset, 0, 0, .99, 10);
    subplot(ceil(sqrt(nExpts)),ceil(sqrt(nExpts)),a);
    image(gaborimage);
    set(gca,'XTick',[],'YTick',[]); axis square;
    if (~success_s(a))
        set(gca,'Xcolor',[1 1 0],'Ycolor',[1 1 0]);
    end
    drawnow;
end

%%
% Section 7
%
% Looking for relationships between saccade-related activity modulation and
% receptive field parameters.

[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\LvsM.txt');
nExpts = size(fnames, 1);
h1 = figure;
h2 = figure;
data = [];
BINWIDTH = .01;
timebins = [-.25:BINWIDTH:.5];
for a = 1:nExpts;
    stro = {};
    for i = 1:size(fnames{a},2)
        tmpstro = nex2stro(findfile(char(fnames{a}(i))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
    if (isempty(stro.sum.analog.sigid))
        continue;
    end
    maxT = 8;
    framerate = stro.sum.exptParams.framerate;
    nstixperside = stro.sum.exptParams.nstixperside;
    ntrials = length(stro.sum.absTrialNum);
    stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoffidx = find(strcmp(stro.sum.trialFields(1,:),'all_off'));
    nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    sigmaidxs = strmatch('sigma',stro.sum.trialFields(1,:));

    mondist = 100; % cm
    screenwidth = 36; %cm
    screenwidthinpix = 1024; % singlewide pixels
    pixpercm = screenwidthinpix/screenwidth;
    cmperdeg = screenwidth/(2*atan2(screenwidth/2, mondist)*180/pi);
    pixperdeg = pixpercm*cmperdeg;
    stixperdeg = pixperdeg/(2*stro.sum.exptParams.npixperstix);

    L = stro.trial(:,noisetypeidx) == 1;
    stro.ras(~L,:) = [];
    stro.trial(~L,:) = [];

    out = getWhtnsStats(stro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, stro.sum.rasterCells{spikeIdx(a)});
    tmpstro = [];
    STAs = out{1};
  
    energy = sum(STAs.^2);
    whichframe = find(energy == max(energy));

    % fitting the gabor
    im = reshape(STAs(:,whichframe),nstixperside*nstixperside,3);
    [u,s,v] = svd(im');
    im = reshape(v(:,1),[nstixperside nstixperside]);
    struct = gaborfit(im);

    data(a).sigma = struct.sigma./stixperdeg;
    data(a).lambda = struct.lambda./stixperdeg;
    data(a).gamma = struct.gamma;
    data(a).theta = struct.theta;
    data(a).success = struct.exitflag;
    data(a).name = stro.sum.fileName;
    data(a).whichframe = whichframe;

    sacstats = getSacData(stro);
    close;

    PSTH = zeros(1,length(timebins));
    for j = 1:size(stro.trial,1)
        stimon_t = stro.trial(j,stimonidx);
        numframes = stro.trial(j,nframesidx);
        stimoff_t = stimon_t+numframes/stro.sum.exptParams.framerate;
        st = sacstats.starttimes{j};
        Lsac = (st > stimon_t) & (st < stimoff_t-.1); %.1 to avoid a saccade that leaves the window
        if any(Lsac)
            Lsac(sacstats.amplitudes{j}(Lsac) > 2) = [];
            spiketimes = [];
            for k = find(Lsac')
                tmp = stro.ras{j,spikeIdx(a)}-sacstats.starttimes{j}(k);
                spiketimes = [spiketimes; tmp((tmp > timebins(1)-BINWIDTH/2) & (tmp < timebins(end)+BINWIDTH/2))];
            end
            [n,x] = hist(spiketimes, timebins);
            PSTH = PSTH + n./(BINWIDTH*sum(Lsac));
        end
    end
    data(a).PSTH = PSTH./sum(L);
    
    % Plotting STAs
    maxSTA = max(STAs(:,whichframe));
    minSTA = min(STAs(:,whichframe));

    potentialnormfactors = [(.5-eps)./(maxSTA-.5); (-.5+eps)./(minSTA-.5)];
    potentialnormfactors(potentialnormfactors < 0) = []; % if min > mu or max < mu
    normfactor = min(potentialnormfactors);
    muvect = reshape(repmat([.5 .5 .5],nstixperside^2,1),nstixperside^2*3,1);

    % Plotting
    figure(h1);
    subplot(ceil(sqrt(nExpts)),ceil(sqrt(nExpts)),a);
    STA = normfactor*(STAs(:,whichframe)-muvect)+muvect;
    STA = reshape(STA,[nstixperside nstixperside 3]);
    image(STA);
    set(gca,'XTick',[],'YTick',[]); axis square;
    
    gaborrgb = u(:,1);
    gaborimage = DrawGaborEdge([.5 .5 .5], gaborrgb/2, [0 0 0], struct.theta, struct.lambda, struct.sigma, struct.gamma, struct.phi, struct.xoffset, struct.yoffset, 0, 0, .99, 10);
    figure(h2);
    subplot(ceil(sqrt(nExpts)),ceil(sqrt(nExpts)),a);
    image(gaborimage);
    set(gca,'XTick',[],'YTick',[]); axis square;
    drawnow;
end

RFparams = []; PSTHs = [];
for i = 1:length(data)
    PSTHs(i,:) = (data(i).PSTH-mean(data(i).PSTH))./norm(data(i).PSTH);
    RFparams(i,:) = [data(i).sigma data(i).lambda data(i).gamma data(i).theta];
end

sacmodindex = sum(PSTHs(:,timebins>.1 & timebins<.2),2)-sum(PSTHs(:,timebins>0 & timebins<.1),2);
template = mean(PSTHs);
template = template./norm(template);
projs = PSTHs*template';
[rho, pval] = corr([sacmodindex, projs, RFparams],'type','spearman')

figure
L = RFparams(:,2) < 100;
plot(sacmodindex(L), RFparams(L,2),'k.');
% Cells with large lambdas (low preferred spatial frequencies) tend to show
% a lot of modulation.

% Below, for plotting PSTHs and names of individual cells.
% i = 77; figure(2); plot(timebins,PSTHs(i,:)); data(i).name
%%
% Section 7.1 (Continued from above)
% Does neural response latency correlate with postsaccadic modulation?
latency = [data(:).whichframe];
smoothedpsths = conv2(PSTHs,[1 1 1 1]/4);
window = [30:50];
sacpsthlat = [];
for i = 1:size(PSTHs,1)
   minidx = find(smoothedpsths(i,window) == min(smoothedpsths(i,window)));
   maxidx = find(smoothedpsths(i,window) == max(smoothedpsths(i,window)));
   sacpsthlat(i,:) = [minidx maxidx];
end
L = sacpsthlat(:,1) < sacpsthlat(:,2);

figure; subplot(2,2,1); hold on;
plot(sacpsthlat(~L,1), sacpsthlat(~L,2),'k.')
plot(sacpsthlat(L,1), sacpsthlat(L,2),'r.')
xlabel('Dip Latency'); ylabel('Peak latency');
[r,p] = corr([latency(L)', sacpsthlat(L,:)]);

subplot(2,2,2); hold on;
plot(latency(L), sacpsthlat(L,1),'r.');
xlabel('STA Latency'); ylabel('Dip latency'); title(['r = ',num2str(r(1,2)),' p = ',num2str(p(1,2))]);
lsline;

subplot(2,2,3); hold on;
plot(latency(L), sacpsthlat(L,2),'r.');
xlabel('STA Latency'); ylabel('Peak latency'); title(['r = ',num2str(r(1,3)),' p = ',num2str(p(1,3))]);
lsline;


%%
% Section 7.2
% Does postsaccadic modulation vary with S-cone weight?

filenamelist = {'N:\NexFiles\nexfilelists\Greg\LvsM.txt',...
    'N:\NexFiles\nexfilelists\Greg\Lum.txt'};
fnames = []; spikeIdx = []; whichlist = [];
for i = 1:length(filenamelist)
    [tmpfnames tmpspikeIdx] = fnamesFromTxt2(char(filenamelist{i}));
    fnames = [fnames; tmpfnames];
    spikeIdx = [spikeIdx; tmpspikeIdx];
    whichlist = [whichlist; i*ones(length(tmpfnames),1)];
end

nExpts = size(fnames, 1);
data = [];
BINWIDTH = .01;
timebins = [-.25:BINWIDTH:.5];
for a = 1:nExpts;
    stro = {};
    for i = 1:size(fnames{a},2)
        tmpstro = nex2stro(findfile(char(fnames{a}(i))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
%     if (isempty(stro.sum.analog.sigid))
%         continue;
%     end
    maxT = 8;
    framerate = stro.sum.exptParams.framerate;
    nstixperside = stro.sum.exptParams.nstixperside;
    ntrials = length(stro.sum.absTrialNum);
    stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoffidx = find(strcmp(stro.sum.trialFields(1,:),'all_off'));
    nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    sigmaidxs = strmatch('sigma',stro.sum.trialFields(1,:));
    
    % Reconstructing the M matrix
    fundamentals = stro.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = stro.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    
    L = stro.trial(:,noisetypeidx) == 1;
    stro.ras(~L,:) = [];
    stro.trial(~L,:) = [];

    out = getWhtnsStats(stro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, stro.sum.rasterCells{spikeIdx(a)});
    tmpstro = [];
    STAs = out{1};
  
    energy = sum(STAs.^2);
    whichframe = find(energy == max(energy));
    im = reshape(STAs(:,whichframe),nstixperside*nstixperside,3);
    [u,s,v] = svd(im');
    data(a).coneweights = inv(M')*u(:,1);
    
    sacstats = getSacData(stro);
    close;

    PSTH = zeros(1,length(timebins));
    for j = 1:size(stro.trial,1)
        stimon_t = stro.trial(j,stimonidx);
        numframes = stro.trial(j,nframesidx);
        stimoff_t = stimon_t+numframes/stro.sum.exptParams.framerate;
        st = sacstats.starttimes{j};
        Lsac = (st > stimon_t) & (st < stimoff_t-.1); %.1 to avoid a saccade that leaves the window
        if any(Lsac)
            Lsac(sacstats.amplitudes{j}(Lsac) > 2) = [];
            spiketimes = [];
            for k = find(Lsac')
                tmp = stro.ras{j,spikeIdx(a)}-sacstats.starttimes{j}(k);
                spiketimes = [spiketimes; tmp((tmp > timebins(1)-BINWIDTH/2) & (tmp < timebins(end)+BINWIDTH/2))];
            end
            [n,x] = hist(spiketimes, timebins);
            PSTH = PSTH + n./(BINWIDTH*sum(Lsac));
        end
    end
    data(a).PSTH = PSTH./sum(L);
end
coneweights = [data.coneweights]./repmat(sum(abs([data.coneweights])),[3,1]);

PSTHs = [];
for i = 1:length(data)
    PSTHs(i,:) = (data(i).PSTH-mean(data(i).PSTH))./norm(data(i).PSTH);
end

sacmodindex = sum(PSTHs(:,timebins>0 & timebins<-.12),2)-sum(PSTHs(:,timebins>0 & timebins<.2),2);
template = mean(PSTHs);
template = template./norm(template);
projs = PSTHs*template';
[rho, pval] = corr([sacmodindex, projs, abs(coneweights(3,:))'],'type','spearman')  % S
[rho, pval] = corr([sacmodindex, projs, abs(coneweights(1,:)-coneweights(2,:))'],'type','spearman')  % L-M
% There's basically no correlation between these indices of saccade
% modulation and how much L-M there is in the STA.

figure; axes; hold on;
imagesc(PSTHs);
plot([1 size(PSTHs,2)],[find(whichlist == 1, 1, 'last' ) find(whichlist == 1, 1, 'last' )]+.5,'y-');
colormap(gray);
axis tight;
set(gca,'Ytick',[]);
y = [timebins(1):.25:timebins(end)];
x = (y-timebins(1))/BINWIDTH+1;
set(gca,'XTick',x);
set(gca,'XTickLabel',y);

[h,p] = ttest2(sacmodindex(whichlist == 1), sacmodindex(whichlist == 2))
[h,p] = ttest2(projs(whichlist == 1), projs(whichlist == 2))

%%
% Section 8
%
% Looking to see whether there is a difference in the saccade-triggered
% PSTH when the saccade is aligned with (or orthogonal to) the direction of
% a fixational saccade.
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\Lum.txt');
% inline: saccade is aligned with the preferred orientation
% orthogonal1: saccade is orthogonal preferred orientation in one direction
% orthogonal2: saccade is orthogonal preferred orientation in opposite
% direction
nExpts = size(fnames, 1);
inlinedata = [];
orthodata1 = [];
orthodata2 = [];
BINWIDTH = .01;
timebins = [-.25:BINWIDTH:.5];
for a = 1:nExpts;
    stro = {};
    for i = 1:size(fnames{a},2)
        tmpstro = nex2stro(findfile(char(fnames{a}(i))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
    if (isempty(stro.sum.analog.sigid))
        continue;
    end
    maxT = 8;
    framerate = stro.sum.exptParams.framerate;
    nstixperside = stro.sum.exptParams.nstixperside;
    ntrials = length(stro.sum.absTrialNum);
    stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoffidx = find(strcmp(stro.sum.trialFields(1,:),'all_off'));
    nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    sigmaidxs = strmatch('sigma',stro.sum.trialFields(1,:));

    L = stro.trial(:,noisetypeidx) == 1;
    stro.ras(~L,:) = [];
    stro.trial(~L,:) = [];

    out = getWhtnsStats(stro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, stro.sum.rasterCells{spikeIdx(a)});
    tmpstro = [];
    STAs = out{1};
  
    energy = sum(STAs.^2);
    whichframe = find(energy == max(energy));

    % fitting the gabor
    im = reshape(STAs(:,whichframe),nstixperside*nstixperside,3);
    [u,s,v] = svd(im');
    im = reshape(v(:,1),[nstixperside nstixperside]);
    gaborstruct = gaborfit(im);
    if(gaborstruct.lambda > 20)
        disp('Omitting very low SF cell');
        continue;
    end
    gaborstruct.lambda
    
    sacstats = getSacData(stro);
    close;
    
    PSTH = zeros(3,length(timebins));
    saccadecounters = zeros(3,1);
    for j = 1:size(stro.trial,1)
        stimon_t = stro.trial(j,stimonidx);
        numframes = stro.trial(j,nframesidx);
        stimoff_t = stimon_t+numframes/stro.sum.exptParams.framerate;
        st = sacstats.starttimes{j};
        dir = sacstats.directions{j};
        for k = 1:3
            angulardiff = mod(dir-gaborstruct.theta,2*pi);
            if (k == 1)  % inline
                Ldir = ((angulardiff < pi/4) | (angulardiff > 7*pi/4)) | ((angulardiff > 3*pi/4 & angulardiff < 5*pi/4));
            elseif (k == 2) % orthogonal 1
                Ldir = (angulardiff > pi/4 & angulardiff < 3*pi/4);
            elseif (k == 3) % orthogonal 2 
                Ldir = (angulardiff > 5* pi/4 & angulardiff < 7*pi/4);
            end
            Lsac = Ldir & (st > stimon_t) & (st < stimoff_t-.1); %.1 to avoid a saccade that leaves the window
            if any(Lsac)
                spiketimes = [];
                for m = find(Lsac')
                    tmp = stro.ras{j,spikeIdx(i)}-sacstats.starttimes{j}(m);
                    spiketimes = [spiketimes; tmp((tmp > timebins(1)-BINWIDTH/2) & (tmp < timebins(end)+BINWIDTH/2))];
                end
                [n,x] = hist(spiketimes, timebins);
                PSTH(k,:) = PSTH(k,:)+n;
                saccadecounters(k) = saccadecounters(k)+sum(Lsac);
            end
        end
    end
    inlinedata(a,:)= PSTH(1,:)./(BINWIDTH*saccadecounters(1));
    orthodata1(a,:)= PSTH(2,:)./(BINWIDTH*saccadecounters(2));
    orthodata2(a,:)= PSTH(3,:)./(BINWIDTH*saccadecounters(3));
end
peakfr = max([inlinedata orthodata1 orthodata2],[],2);
figure; subplot(6,1,1);
imagesc(inlinedata./repmat(peakfr,1,length(timebins)));
set(gca,'XTick',find(rem(timebins,.1) == 0));
set(gca,'XTickLabel',timebins(get(gca,'Xtick')));
subplot(6,1,2);
plot(timebins, nanmean(inlinedata./repmat(peakfr,1,length(timebins))));

subplot(6,1,3);
imagesc(orthodata1./repmat(peakfr,1,length(timebins)));
set(gca,'XTick',find(rem(timebins,.1) == 0));
set(gca,'XTickLabel',timebins(get(gca,'Xtick')));
subplot(6,1,4);
plot(timebins, nanmean(orthodata1./repmat(peakfr,1,length(timebins))));

subplot(6,1,5);
imagesc(orthodata1./repmat(peakfr,1,length(timebins)));
set(gca,'XTick',find(rem(timebins,.1) == 0));
set(gca,'XTickLabel',timebins(get(gca,'Xtick')));
subplot(6,1,6);
plot(timebins, nanmean(orthodata2./repmat(peakfr,1,length(timebins))));

% Paired analysis
figure
imagesc(orthodata1-orthodata2)
imagesc(inlinedata./repmat(peakfr,1,length(timebins))-orthodata2./repmat(peakfr,1,length(timebins)))
imagesc(inlinedata./repmat(peakfr,1,length(timebins))-orthodata1./repmat(peakfr,1,length(timebins))-orthodata2./repmat(peakfr,1,length(timebins)))
%
% Looking at projections on a template
L = ~all(inlinedata == 0,2)
[u,s,v] = svd([inlinedata(L,:); orthodata1(L,:); orthodata2(L,:)])
u = reshape(u(:,1),[size(u(:,1),1)/3,3])
hist(u,15);
% They are all about the same.
figure; axes; hold on;
plot(u(:,1),mean([u(:,2),u(:,3)],2),'k.');
plot([0 .3],[0 .3],'k-');
%%
% Section 9
% Comparing stimulus tuning measured with gratings and white noise.

data = [];
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\Gratings&WN.txt');
MINRESP = 10;  % At least 10 sp/sec to preferred grating.

for a = 1:size(fnames,1)
    WN = {}; GT = {};
    for i = 1:size(fnames{a},2)
        tmpstro = nex2stro(findfile(char(fnames{a}(i))));
        if tmpstro.sum.paradigmID == 150
            GT = tmpstro;
        elseif (isempty(WN))
            WN = tmpstro;
        else
            WN = strocat(WN, tmpstro);
        end
    end
 
    maxT = 8;
    nstixperside = WN.sum.exptParams.nstixperside;
    noisetypeidx = find(strcmp(WN.sum.trialFields(1,:),'noise_type'));
    mondist = 100; % cm
    screenwidth = 36; %cm
    screenwidthinpix = 1024; % singlewide pixels
    pixpercm = screenwidthinpix/screenwidth;
    cmperdeg = screenwidth/(2*atan2(screenwidth/2, mondist)*180/pi);
    pixperdeg = pixpercm*cmperdeg;
    stixperdeg = pixperdeg/(2*WN.sum.exptParams.npixperstix);
    % Reconstructing the M matrix and gamma table
    fundamentals = WN.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = WN.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    
    Lgun = WN.trial(:,noisetypeidx) == 1;
    Lcone = WN.trial(:,noisetypeidx) == 2;
 
    if (sum(Lgun) > sum(Lcone))
        noisetype = 1;
        WN.ras(~Lgun,:) = [];
        WN.trial(~Lgun,:) = [];
        
        out = getWhtnsStats(WN,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, WN.sum.rasterCells{spikeIdx(a)});
        STAs = out{1};
        
        energy = sum(STAs.^2);
        whichframe = find(energy == max(energy));
        
        % fitting the gabor
        im = reshape(STAs(:,whichframe),nstixperside*nstixperside,3);
        [u,s,v] = svd(im');
        im = reshape(v(:,1),[nstixperside nstixperside]);
        gaborstruct = gaborfit(im);
        gaborstruct.sf = stixperdeg./gaborstruct.lambda;
        gaborstruct.rgb = u(:,1);
        gaborstruct.lms = inv(M')*u(:,1);
        gaborstruct.lms = gaborstruct.lms./sum(abs(gaborstruct.lms));
    else
        noisetype = 2;
        WN.ras(~Lcone,:) = [];
        WN.trial(~Lcone,:) = [];
        
        conesigmas = [WN.trial(1,strcmp(WN.sum.trialFields(1,:),'sigma1'));...
            WN.trial(1,strcmp(WN.sum.trialFields(1,:),'sigma2'));...
            WN.trial(1,strcmp(WN.sum.trialFields(1,:),'sigma3'))];
        out = getWhtnsStats(WN,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, WN.sum.rasterCells{spikeIdx(a)});
        STAs = out{1};
        
        energy = sum(STAs.^2);
        whichframe = find(energy == max(energy));
        
        % fitting the gabor
        im = reshape(STAs(:,whichframe),nstixperside*nstixperside,3);
        [u,s,v] = svd(im');
        im = reshape(v(:,1),[nstixperside nstixperside]);
        gaborstruct = gaborfit(im);
        gaborstruct.sf = stixperdeg./gaborstruct.lambda;
        u(:,1) = u(:,1)./conesigmas;
        gaborstruct.lms = u(:,1)./sum(abs(u(:,1)));
        gaborstruct.rgb = inv(M')*u(:,1);
        gaborstruct.rgb = gaborstruct.lms./sum(abs(gaborstruct.lms));
    end
    
    gratingstruct = getGratingTuning(GT, spikeIdx(a));
    if (gratingstruct.prefcolor * gaborstruct.lms < 0)
        gratingstruct.prefcolor = -gratingstruct.prefcolor;
    end
    if (max(gratingstruct.colresp(:,1))-gratingstruct.baselines(1) < MINRESP)
        data(a,:) = nan(1,12);
    else
        data(a,:) = [gaborstruct.sigma gaborstruct.theta gratingstruct.preforient...
                  gaborstruct.sf gratingstruct.prefSF gaborstruct.lms'...
                  gratingstruct.prefcolor noisetype]
    end
end

% First orientation tuning
figure; axes; hold on;
plot([0 pi],[0 pi],'k:');
data(:,3) = mod(data(:,3),pi);
L = data(:,3) > pi/2;
data(L,3) = data(L,3)-pi;
for i =1:size(data,1)
    if (data(i,4) > .7)  % No orientation if lowpass
        h = plot(mod(data(i,2),pi),mod(data(i,3),pi),'k.','MarkerSize',20);
        set(h,'ButtonDownFcn',['disp(''',num2str(i),': ',char(fnames{i}(1)),''')']);
        set(h,'ButtonDownFcn',['disp(''',num2str(i),': ',char(fnames{i}(1)),''')']);
    end
end
axis square;
set(gca,'XLim',[0 pi],'YLim',[0 pi]);
xlabel('STA orientation (rad)'); ylabel('Grating orientation (rad)');

% Spatial frequency tuning
figure; axes; hold on;
plot([0 4],[0 4],'k:');
for i = 1:size(data,1)
    h = plot(data(i,4),data(i,5),'k.','MarkerSize',20);
    set(h,'ButtonDownFcn',['disp(''',num2str(i),': ',char(fnames{i}(1)),''')']);
end
xlabel('STA spatial freq (dva)'); ylabel('Grating spatial freq (dva)');
axis square;
set(gca,'XLim',[0 4],'YLim',[0 4]);

figure; axes; hold on;
set(gca,'XLim',[-1 1],'Ylim',[-1 1]);
plot([1 0 -1 0 1],[0 1 0 -1 0],'k:');
axis square;

for i = 1:size(data,1)
    plot([data(i,6) data(i,9)], [data(i,7) data(i,10)],'k-');
    h1 = plot(data(i,6),data(i,7),'k.','MarkerSize',10);  % white noise
    h2 = plot(data(i,9),data(i,10),'m.','MarkerSize',20); % gratings
    set(h1,'ButtonDownFcn',['disp(''',num2str(i),': ',char(fnames{i}(1)),''')']);
    set(h2,'ButtonDownFcn',['disp(''',num2str(i),': ',char(fnames{i}(1)),''')']);
end
legend([h1;h2],{'STA','Gratings'});

% Now let's see how many cells have cone weights that actually switch sign.
Lnan = ~any(isnan(data(:,[6 7 8 9 10 11])),2);
STA_cw = data(:,[6 7 8]);
Grat_cw = data(:,[9 10 11]);
L = any(sign(STA_cw(:,[1 2])) ~= sign(Grat_cw(:,[1 2])),2)
Ldiff = L&Lnan;
Lsame = ~L&Lnan;

%%
% Section 10
%
% STAs and PC1s for a bunch of cells (and spatial maps from SVD) .
s = load('T_cones_smj10');
fundamentals = s.T_cones_smj10;
[fnames, spikeIdx] = fnamesFromTxt2(fullfile(nexfilepath,'nexfilelists','Greg','WhiteNoise','ColorOpponent.txt'));
%[fnames, spikeIdx] = fnamesFromTxt2(fullfile(nexfilepath,'nexfilelists','Greg','WhiteNoise','Lum.txt'));
%[fnames, spikeIdx] = fnamesFromTxt('WhiteNoiseLGN_forIS','cellClass',{'K'},'subjID',{'A'});

nExpts = size(fnames, 1);
h1 = figure;
h2 = figure;
h3 = figure;

data = [];
for a = 1:nExpts
    fnames{a}
    stro = {};
    for i = 1:size(fnames{a},2)
        tmpstro = nex2stro(findfile(char(fnames{a}(i))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
    maxT = 8;
    framerate = stro.sum.exptParams.framerate;
    nstixperside = stro.sum.exptParams.nstixperside;
    edgepixelmat = zeros(nstixperside,nstixperside);
    edgepixelmat(:,1) = 1; edgepixelmat(:,end) = 1;
    edgepixelmat(1,:) = 1; edgepixelmat(end,:) = 1;
    edgepixels = edgepixelmat(:);

    ntrials = length(stro.sum.absTrialNum);
    stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoffidx = find(strcmp(stro.sum.trialFields(1,:),'all_off'));
    nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    sigmaidxs = strmatch('sigma',stro.sum.trialFields(1,:));
    bkgndRGB = stro.trial(end,strncmp(stro.sum.trialFields(1,:),'bkgnd',5));
    gammaTable = stro.sum.exptParams.gamma_table;
    gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
    gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
    invgamma = InvertGamma(gammaTable, 0);
    bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];

    monspd = reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3);
    if size(monspd,1) == 101
       monspd = SplineSpd([380:4:780]',monspd,[380:5:780]');
    end
    if size(monspd,1) ~= 81
        error('monspd is the wrong size');
    end
    
    % Hacky way of making the rgb -> cc matrix
    M = fundamentals*monspd; % If rgb is represented in deltas from bkgnd then Mrgbtocc*rgb is in cone contrast.
    bkgndlms = M*bkgndrgb;
    M = diag(1./bkgndlms)*M;

    L = stro.trial(:,noisetypeidx) == 1;
    stro.ras(~L,:) = [];
    stro.trial(~L,:) = [];

    out = getWhtnsStats(stro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, stro.sum.rasterCells{spikeIdx(a)});
    STAs = out{1};
    STCs = out{2};
    nspikes = out{3};
    
    energy = sum(STAs.^2);
    whichframe = find(energy == max(energy));
    
    energy = sum(reshape(STAs(:,whichframe).^2,nstixperside^2,3),2);
    whichpixel = find(energy == max(energy));
    
    STArgb = STAs(whichpixel+[0, nstixperside^2, 2*nstixperside^2],:);
    energy = sum(STArgb.^2);
    whichframe = find(energy == max(energy));
    whichframes = whichframe+[-2:2];
    whichframes(whichframes < 1 | whichframes > maxT) = [];
    STArgbtrunc = STArgb(:,whichframes);
    [u,~,~] = svd(STArgbtrunc);
    signal = mean(STArgb(:,whichframe).^2);
    noise = mean(STAs(logical(edgepixels), whichframe).^2);
    
    if u(:,1)'*STArgb(:,whichframe) < 0
        u(:,1) = -u(:,1);
    end
    STAlms = inv(M')*u(:,1);
    
    STC = reshape(STCs(:,whichframe),sqrt(size(STCs,1)),sqrt(size(STCs,1)));
    [PC1,d] = eigs(STC,1);
    
    % Normalizing images
    template = reshape([1:nstixperside^2],nstixperside,nstixperside);
    edgepixels = [template(:,1); template(1,[2:end-1])'; template(end,[2:end-1])'; template(:,end)];
    edgepixelidxs = [edgepixels; edgepixels+nstixperside^2; edgepixels+2*(nstixperside^2)];
    PCelements = PC1(edgepixelidxs,:,:);
    PCsds = std(PCelements);    % One std calculated per PC
    PC1 = PC1.*repmat(std(STAs(:,1))./PCsds,[size(PC1,1),1,1]);

    rowidxs = reshape([1:3*nstixperside^2],[nstixperside^2 3]);
    maxes = []; mins = [];
    imagevectors = [STAs(:,whichframe), PC1];
    for i = 1:3
        maxes = [maxes; max(max(imagevectors(rowidxs(:,i),:)))];
        mins = [mins; min(min(imagevectors(rowidxs(:,i),:)))];
    end
    potentialnormfactors = [(1-[.5; .5; .5]-eps)./(maxes-[.5; .5; .5]); (-[.5; .5; .5]+eps)./(mins-[.5; .5; .5])];
    % 'eps' in above line is a kludge that is required for avoiding
    % out of bounds errors.
    potentialnormfactors(potentialnormfactors < 0) = []; % if min > mu or max < mu
    normfactor = min(potentialnormfactors);    
    muvect = reshape(repmat([.5 .5 .5],nstixperside^2,1),nstixperside^2*3,1);

    % Plotting
    figure(h1);
    subplot(ceil(sqrt(nExpts)),ceil(sqrt(nExpts)),a);
    STA = normfactor*(STAs(:,whichframe)-muvect)+muvect;
    STA = reshape(STA,[nstixperside nstixperside 3]);
    image(STA);
    set(gca,'XTick',[],'YTick',[]); axis square;
    %title(num2str(signal/noise));

    figure(h2);
    subplot(ceil(sqrt(nExpts)),ceil(sqrt(nExpts)),a);
    v = normfactor*(PC1-muvect)+muvect;
    v = reshape(v,[nstixperside nstixperside 3]);
    image(v);
    set(gca,'XTick',[],'YTick',[]); axis square;
    
    figure(h3);
    [u,~,~] = svd(reshape(STAs(:,whichframe),nstixperside^2,3));
    subplot(ceil(sqrt(nExpts)),ceil(sqrt(nExpts)),a);
    image(reshape(255*(u(:,1)./(2.1*max(abs(u(:,1))))+.5),nstixperside,nstixperside));
    colormap(gray(255))
    set(gca,'XTick',[],'YTick',[]); axis square;
    drawnow;
   
    data = [data; STAlms', signal/noise, stro.sum.exptParams.rf_x/10 stro.sum.exptParams.rf_y/10];
end
coneweights = data(:,[1 2 3]);
normalized_coneweights = coneweights./repmat(sum(abs(coneweights),2),1,3);
rflocations = data(:,[5 6]);
LgoodSTA = data(:,4) > 10 ;% & sqrt(sum(rflocations.^2,2)) > 5; % arbitrary threshold

figure; subplot(2,1,1); hold on;
for i = 1:size(rflocations,1)
    plot(rflocations(i,1),rflocations(i,2),'.');
   % text(rflocations(i,1),rflocations(i,2),num2str(i));
end
axis equal
axis square
L_posS = normalized_coneweights(:,3) > 0;
subplot(2,1,2); hold on;
plot(normalized_coneweights(LgoodSTA&L_posS,1),normalized_coneweights(LgoodSTA&L_posS,2),'ko','MarkerFaceColor','black')
plot(normalized_coneweights(LgoodSTA&~L_posS,1),normalized_coneweights(LgoodSTA&~L_posS,2),'ko','MarkerFaceColor','white')

axis square;
plot([-1 0 1 0 -1],[0 -1 0 1 0],'k-')

% Any change when we moved to ProPixx?
Lpropixx = zeros(nExpts,1);
for a = 1:nExpts
    filename = fnames{a}{1};
    if str2num(filename([6 7])) >= 16 % 2016
        Lpropixx(a) = 1;
    end
end

figure; axes; hold on;
plot(normalized_coneweights(LgoodSTA&Lpropixx,1),normalized_coneweights(LgoodSTA&Lpropixx,2),'ks','MarkerFaceColor','black');
plot(normalized_coneweights(LgoodSTA&~Lpropixx,1),normalized_coneweights(LgoodSTA&~Lpropixx,2),'ms','MarkerFaceColor','magenta');
axis square;
plot([-1 0 1 0 -1],[0 -1 0 1 0],'k-')


%%
% Section 10.1
%
% STAs (over time) for a bunch of cells. 
maxT = 7;
ncellsperfig = 5;
[fnames, spikeIdx] = fnamesFromTxt('WhiteNoiseLGN_forIS','cellClass',{'K'});

nExpts = size(fnames, 1);
counter = 1;
for a = 1:nExpts
    if (mod(a,ncellsperfig) == 1)
        figure;
        counter = 1;
    else
        counter = counter+1;
    end
    stro = {};
    for i = 1:size(fnames{a},2)
        tmpstro = nex2stro(findfile(char(fnames{a}(i)),'~/Desktop/LGN data'));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
    framerate = stro.sum.exptParams.framerate;
    nstixperside = stro.sum.exptParams.nstixperside;
    ntrials = length(stro.sum.absTrialNum);
    stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoffidx = find(strcmp(stro.sum.trialFields(1,:),'all_off'));
    nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    sigmaidxs = strmatch('sigma',stro.sum.trialFields(1,:));

    L = stro.trial(:,noisetypeidx) == 1;
    stro.ras(~L,:) = [];
    stro.trial(~L,:) = [];

    out = getWhtnsStats(stro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, stro.sum.rasterCells{spikeIdx(a)});
    STAs = out{1};
    normSTAs = .5*STAs./max(abs(STAs(:)))+.5;
    for whichframe = 1:maxT
        % Plotting
        subplot(ncellsperfig,maxT,(counter-1)*maxT+whichframe);
        STA = normSTAs(:,whichframe);
        nstixperside = sqrt(length(STA)/3);
        STA = reshape(STA,[nstixperside nstixperside 3]);
        image(STA);
        set(gca,'XTick',[],'YTick',[]); axis square;
        drawnow;
    end
end

%%
% Section 10.2
% Looking at cone weights derived from the single hottest pixel or from all
% of the pixels in the cRF.
CHI2CRIT = .8; % For finding stixels in RF
s = load('T_cones_smj10');
fundamentals = s.T_cones_smj10;
[fnames, spikeIdx] = fnamesFromTxt2(fullfile(nexfilepath,'nexfilelists','Greg','WhiteNoise','LGN_Utu.txt'));
nExpts = size(fnames, 1);
data = [];
for a = 1:nExpts
    stro = {};
    for i = 1:size(fnames{a},2)
        tmpstro = nex2stro(findfile(char(fnames{a}(i))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
    maxT = 8;
    framerate = stro.sum.exptParams.framerate;
    nstixperside = stro.sum.exptParams.nstixperside;
    edgepixelmat = zeros(nstixperside,nstixperside);
    edgepixelmat(:,1) = 1; edgepixelmat(:,end) = 1;
    edgepixelmat(1,:) = 1; edgepixelmat(end,:) = 1;
    edgepixels = edgepixelmat(:);

    ntrials = length(stro.sum.absTrialNum);
    stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoffidx = find(strcmp(stro.sum.trialFields(1,:),'all_off'));
    nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    sigmaidxs = strmatch('sigma',stro.sum.trialFields(1,:));
    bkgndRGB = stro.trial(end,strncmp(stro.sum.trialFields(1,:),'bkgnd',5));
    gammaTable = stro.sum.exptParams.gamma_table;
    gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
    gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
    invgamma = InvertGamma(gammaTable, 0);
    bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];

    monspd = reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3);
    if size(monspd,1) == 101
       monspd = SplineSpd([380:4:780]',monspd,[380:5:780]');
    end
    if size(monspd,1) ~= 81
        error('monspd is the wrong size');
    end
    
    % Hacky way of making the rgb -> cc matrix
    M = fundamentals*monspd; % If rgb is represented in deltas from bkgnd then Mrgbtocc*rgb is in cone contrast.
    bkgndlms = M*bkgndrgb;
    M = diag(1./bkgndlms)*M;

    L = stro.trial(:,noisetypeidx) == 1;
    stro.ras(~L,:) = [];
    stro.trial(~L,:) = [];

    out = getWhtnsStats(stro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, stro.sum.rasterCells{spikeIdx(a)});
    STAs = out{1};
    nspikes = out{3};

    % Spatial map
    STAs_z = STAs./std(STAs(:,end)); % Using last frame of STA (STS) as a noise estimate. Fine for LGN.
    grandz = zeros([nstixperside nstixperside]);
    for framecounter = 1:maxT
        tmp = reshape(STAs_z(:,framecounter),[nstixperside nstixperside 3]);
        grandz = grandz+sum(tmp.^2,3);
    end
    crit = chi2inv(CHI2CRIT,3*maxT); % 3 color channels
    L = grandz > crit;
    % Now get the block that contains the single hottest stixel
    whichstixelidx = find(grandz(:) == max(grandz(:)));
    [whichstixel_i,whichstixel_j] = ind2sub(size(grandz),whichstixelidx);
    [tmpi,tmpj] = ind2sub(size(L),find(L));
    ij = [tmpi,tmpj];
    T = clusterdata(ij,'linkage','single','distance','euclidean','criterion','distance','cutoff',sqrt(2));
    dominantcluster = T(ij(:,1) == whichstixel_i & ij(:,2) == whichstixel_j);
    clustermat = zeros(nstixperside, nstixperside);
    clustermat(sub2ind(size(clustermat),ij(T==dominantcluster,1),ij(T==dominantcluster,2))) = 1;
    
    % Then get convex hull
    dominantclusteridxs = ij(T==dominantcluster,:); 
    if size(dominantclusteridxs,1) <= 2 % Only one stixel in RF
        whichstixels = sub2ind(size(clustermat),dominantclusteridxs(:,1),dominantclusteridxs(:,2));
    else
        K=convhull(dominantclusteridxs(:,1), dominantclusteridxs(:,2));
        tmp = dominantclusteridxs(K,:);
        [x,y] = meshgrid(1:nstixperside,1:nstixperside);
        inRF = reshape(inpolygon(x(:),y(:),dominantclusteridxs(K,2),dominantclusteridxs(K,1)),[nstixperside, nstixperside]);
        whichframes = find(max(abs(STAs_z)) == max(abs(STAs_z(:))))+[-1 0 1];
        whichstixels = find(inRF);
    end
    tmpSTA = reshape(STAs, [nstixperside^2 3 maxT]);
    tmpSTA = permute(tmpSTA, [2 1 3]);
    
    [u1,~,v1] = svd(squeeze(tmpSTA(:,whichstixelidx,whichframes)));  % single hottest stixel
    %[u2,~,v2] = svd(squeeze(reshape(tmpSTA(:,whichstixels,whichframes),[3 length(whichstixels)*length(whichframes)]))); % whole RF
    [u2,~,v2] = svd(squeeze(mean(tmpSTA(:,whichstixels,whichframes),2))); % whole RF

    if max(abs(v1(:,1))) ~= max(v1(:,1))
        u1 = -u1;
    end
    if u1(:,1)'*u2(:,1) < 0
        u2 = -u2;
    end
    
    single_stixel_cone_weights = inv(M')*u1(:,1);
    whole_RF_cone_weights = inv(M')*u2(:,1);
    data = [data; single_stixel_cone_weights'./sum(abs(single_stixel_cone_weights)) whole_RF_cone_weights'./sum(abs(whole_RF_cone_weights)) max(grandz(:))];
end

L = data(:,end) > 80; % arbitrary threshold on STA SNR
figure; subplot(2,1,1); hold on;
plot(data(L,1),data(L,2),'b.');
plot(data(L,4),data(L,5),'bo');
for i = find(L)'
   plot([data(i,1),data(i,4)],[data(i,2),data(i,5)],'k-')
end
plot([-1 0 1 0 -1],[0 -1 0 1 0],'k-');
axis equal;
axis square;


conepurityidx = [];
for i = 1:size(data,1)
    conepurityidx(i,1) = data(i,[1 2])./norm(data(i,[1 2]))*[1; 1]/sqrt(2);
    conepurityidx(i,2) = data(i,[4 5])./norm(data(i,[4 5]))*[1; 1]/sqrt(2);
end
Lopp = sign(data(:,1)) ~= sign(data(:,2));

subplot(2,1,2); hold on;
plot(conepurityidx(L&Lopp,1),conepurityidx(L&Lopp,2),'ko','MarkerFaceColor','black','Markersize',4)
plot([-1 1],[-1 1],'k-');
xlabel('cone purity single stix');
ylabel('cone purity multiple stix');
axis square;

%%
% Section 11
% Looking for a relationship between the amplitudes of fixational
% saccades and the modulation of firing rate.  For now just comparing
% microsaccades above and below the median amplitude.
datapath = 'N:\NexFiles\Kali';
filelist = 'N:\NexFiles\nexfilelists\Greg\LvsM.txt';
[filenames1, spikenums1] = fnamesFromTxt2(filelist);
filelist = 'N:\NexFiles\nexfilelists\Greg\Lum.txt';
[filenames2, spikenums2] = fnamesFromTxt2(filelist);
filenames = cat(1,filenames1,filenames2);
spikenums = cat(1,spikenums1,spikenums2);

BINWIDTH = .01;
timebins = [-.25:BINWIDTH:.5];
data = nan*ones(length(filenames), length(timebins),2);
for i = 1:size(filenames,1)
    stro = {};
    for f = 1:size(filenames{i},1)
        tmpstro = nex2stro(findfile(char(filenames{i}(f))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
    
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
    nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
    
    if (isempty(stro.sum.analog.sigid))
        continue;
    end
    sacstats = getSacData(stro);
    close;
    amps = [];
    L = logical(stro.trial(:,noisetypeidx)~= 0);
    for j = find(L')
        stimon_t = stro.trial(j,stimonidx);
        numframes = stro.trial(j,nframesidx);
        stimoff_t = stimon_t+numframes/stro.sum.exptParams.framerate;
        st = sacstats.starttimes{j};
        Lsac = (st > stimon_t) & (st < stimoff_t-timebins(end)); %Don't include msacs at end of trial
        amps = [amps; sacstats.amplitudes{j}(Lsac)]; 
    end
    medianamp = median(amps);
    medianamp = 0.5; % Reviewer 1's suggested threshold

    PSTH_small = zeros(1,length(timebins));
    PSTH_big = zeros(1,length(timebins));
    n_small = 0;
    n_big = 0;
    for j = find(L')
        stimon_t = stro.trial(j,stimonidx);
        numframes = stro.trial(j,nframesidx);
        stimoff_t = stimon_t+numframes/stro.sum.exptParams.framerate;
        st = sacstats.starttimes{j};
        Lsac = (st > stimon_t) & (st < stimoff_t-timebins(end)); %Don't include msacs at end of trial
        if any(Lsac)
            if (any( Lsac(sacstats.amplitudes{j}(Lsac) > 2)))
                keyboard
            end
            Lsac(sacstats.amplitudes{j}(Lsac) > 2) = [];
            for k = find(Lsac')
                tmp = stro.ras{j,spikenums(i)}-sacstats.starttimes{j}(k);
                spiketimes = tmp((tmp > timebins(1)-BINWIDTH/2) & (tmp < timebins(end)+BINWIDTH/2));
                [n,x] = hist(spiketimes, timebins);
                if (sacstats.amplitudes{j}(k) >= medianamp)
                    PSTH_big = PSTH_big + n./BINWIDTH;
                    n_big = n_big + 1;
                else
                    PSTH_small = PSTH_small + n./BINWIDTH;
                    n_small = n_small + 1;
                end
            end
        end
    end
    data(i,:,1) = PSTH_small./n_small;
    data(i,:,2) = PSTH_big./n_big;
end
baseline = squeeze(mean(data(:,timebins < 0,:),2));
dip = squeeze(mean(data(:,timebins >= .05 & timebins <= .1,:),2));
rebound = squeeze(mean(data(:,timebins >= .12 & timebins <= .3,:),2));

[h,p] = ttest(baseline(:,1)-baseline(:,2))
[h,p] = ttest(dip(:,1)-dip(:,2))
[h,p] = ttest(rebound(:,1)-rebound(:,2))
[h,p] = ttest(dip(:,1)./baseline(:,1)-dip(:,2)./baseline(:,2))
hist(rebound(:,1)-rebound(:,2),20)

% Rebound is delayed for larger amplitude saccades?
plot(mean(data(:,:,1)-data(:,:,2)))

%%
% Section 11.1
% Microsaccade firing rates suppression as a function of saccade amplitude
% Regressions of firing rate as a function of saccade amplitude

datapath = 'N:\NexFiles\Kali';
filelist = 'N:\NexFiles\nexfilelists\Greg\LvsM.txt';
[filenames1, spikenums1] = fnamesFromTxt2(filelist);
filelist = 'N:\NexFiles\nexfilelists\Greg\Lum.txt';
[filenames2, spikenums2] = fnamesFromTxt2(filelist);
filenames = cat(1,filenames1,filenames2);
spikenums = cat(1,spikenums1,spikenums2);

timewin = [0.05 0.1];

data = [];
for i = 1:size(filenames,1)
    stro = {};
    for f = 1:size(filenames{i},1)
        tmpstro = nex2stro(findfile(char(filenames{i}(f))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
    
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
    nframesidx = find(strcmp(stro.sum.trialFields(1,:),'num_frames'));
    
    if (isempty(stro.sum.analog.sigid))
        continue;
    end
    sacstats = getSacData(stro);
    close;
    tmpdata = [];
    L = logical(stro.trial(:,noisetypeidx)~= 0);
    for j = find(L')
        stimon_t = stro.trial(j,stimonidx);
        numframes = stro.trial(j,nframesidx);
        stimoff_t = stimon_t+numframes/stro.sum.exptParams.framerate;
        st = sacstats.starttimes{j};
        Lsac = (st > stimon_t) & (st < stimoff_t-.2); %Don't include msacs at end of trial
        if (any(Lsac))
            for k = find(Lsac)'
                tmp = stro.ras{j,spikenums(i)}-sacstats.starttimes{j}(k);
                spikecount = sum(tmp > timewin(1) & tmp < timewin(2));
                tmpdata = [tmpdata; sacstats.amplitudes{j}(k) spikecount];
            end
        end
    end
    [b,bint] = regress(tmpdata(:,2),[tmpdata(:,1) ones(size(tmpdata,1),1)]);
    [r,p] = corr(tmpdata);
    data = [data; b(1) sign(bint(1,1))== sign(bint(1,2)) r(1,2) p(1,2)];
end
[h,p] = ttest(data(:,1))  % Regression slopes
[h,p] = ttest(data(:,3))  % Correlations

figure; axes; hold on;
L = data(:,4) < 0.05
[n1,x] = hist(data(:,3),linspace(-.3,.2,20));
[n2,x] = hist(data(L,3),linspace(-.3,.2,20));
bar(x,n1,'FaceColor','black');
bar(x,n2,'FaceColor','red');
xlabel('r(sac amp, fr)');
ylabel('count');


figure; axes; hold on;
L = data(:,2) == 1;
[n1,x] = hist(data(:,1),linspace(-3,3,20));
[n2,x] = hist(data(L,1),linspace(-3,3,20));
bar(x,n1,'FaceColor','black');
bar(x,n2,'FaceColor','red');
xlabel('Regression slope');
ylabel('count');

%%
% Section 12
% Comparing neuronal visual response latency between STA and replays
[filenames, spikenums] = fnamesFromTxt2([nexfilepath,'/nexfilelists/Greg/WhiteNoise/WNreplay.txt']);
maxT = 9; % How many frames back to look
VERBOSE = 1; % Make one figure per cell for diagnostics
%MAXISI = 0; % Need at least two PSTH bins exceeding the threshold within MAXISI bins (Do I still need this?)

data = nan*ones(length(filenames), 3);
for i = 1:size(filenames,1)
    stro = {};
    for f = 1:size(filenames{i},1)
        tmpstro = nex2stro(findfile(char(filenames{i}(f))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
    framerate = stro.sum.exptParams.framerate;
    nstixperside = stro.sum.exptParams.nstixperside;
    ntrials = length(stro.sum.absTrialNum);
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    stimonidx = find(strcmp(stro.sum.trialFields(1,:),'stim_on'));
    spikename = stro.sum.rasterCells{spikenums(i)};
    spikeidx = find(strcmp(stro.sum.rasterCells(1,:),spikename));
    Lconenoise = stro.trial(:,noisetypeidx) == 2;
    stro.ras(Lconenoise,:) = [];
    stro.trial(Lconenoise,:) = [];
    out = getWhtnsStats(stro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
    STAs = out{1};
    STCs = out{2};
    nspikes = out{3};
    
    % Below copied and pasted from WNAnalyis
    frametimes = [1:maxT]*1/framerate;
    tmpSTA = reshape(STAs, [nstixperside^2 3 maxT]);
    STAenergy = squeeze(sum(tmpSTA.^2, 2));
    % Picking pixels whose energy (R^2+G^2+B^2) exceeds the first-frame
    % baseline by >3 stds.
    baseline_energy = STAenergy(1:nstixperside^2,1);
    energy_threshold = mean(baseline_energy)+3*std(baseline_energy);
    Lwhichpix = any(STAenergy>energy_threshold,2);
    STAtemporal_waveform = sum(STAenergy(Lwhichpix,:));
    
    % Debugging plots
    if (VERBOSE)
        figure; subplot(2,2,1); hold on;
        image(reshape(Lwhichpix,nstixperside, nstixperside));
        colormap(gray(2)); axis image;
        subplot(2,2,2); plot(STAtemporal_waveform,'ko-');
        rastrialcounter = 0;
    end
    peak_t = frametimes(STAtemporal_waveform == max(STAtemporal_waveform));
    centroid_t = (STAtemporal_waveform/sum(STAtemporal_waveform))*frametimes';
    
    % OK, now getting a measure of latency from the replay trials.
    imageidx = find(strcmp(stro.sum.rasterCells(1,:),'synth_image'));
    allimages = stro.ras(:,imageidx);
    Lsynthimg = [];
    for j = 1:length(allimages)
        Lsynthimg = [Lsynthimg; ~isnan(allimages{j}(1))];
    end
    
    ntimebins = 500;
    timestep = 5; % ms
    timebins_baseline = [-ntimebins:timestep:-1]/1000; % ms
    timebins_driven = [1:timestep:ntimebins]/1000; % ms
    PSTH_baseline = zeros(1,length(timebins_baseline));
    PSTH_driven = zeros(1,length(timebins_driven));
    for j = find(Lsynthimg)'
        stimon_t = stro.trial(j,stimonidx);
        spikes = stro.ras{j,spikeidx}-stimon_t; % Spike times are wrt stimulus onset
        PSTH_baseline = PSTH_baseline + hist(spikes(spikes < 0 & spikes > timebins_baseline(1)),timebins_baseline);
        PSTH_driven = PSTH_driven + hist(spikes(spikes>0 & spikes<timebins_driven(end)),timebins_driven);
        if (VERBOSE)
            rastrialcounter = rastrialcounter+1;
            subplot(2,2,3); hold on; plot(spikes,rastrialcounter*ones(length(spikes),1),'k.'); % Debugging
        end
    end
    
    % PSTHs are sums across trials. Threshold is mean (of spike counts) across time.
    threshold = max(1,mean(PSTH_baseline)+3*std(PSTH_baseline)); % At least one spike
    abovethreshold = find(PSTH_driven>threshold);
    isi = diff(abovethreshold); % interevent interval
    %replay_lat =
    %timebins_driven(abovethreshold(find(isi<=MAXISI/timestep,1))); % Do I need this?
    replay_lat = timebins_driven(abovethreshold(1));
    if (VERBOSE)
        subplot(2,2,3); % Setting the y axis limits correctly on the previous plot.
        set(gca,'ylim',[1 rastrialcounter])
        title(stro.sum.fileName(end-13:end));
        subplot(2,2,4); hold on;
        plot(timebins_driven,PSTH_driven);
        plot([replay_lat replay_lat],[0 max(PSTH_driven)],'k-')
    end
    data(i,:) = [replay_lat peak_t centroid_t];
end
prettycorr(data,{'Replay','Peak','Centroid'}); 

figure; subplot(2,1,1); hold on;
plot(data(:,1),data(:,3),'ko','MarkerFaceColor','black');
plot([0 0.1],[0 0.1]);
xlabel('Replay latency (s)');
ylabel('STA centroid latency (s)');
axis square;

subplot(2,1,2); hold on;
uniquepairs = unique(data(:,[1 2]),'rows');
for i = 1:size(uniquepairs,1)
    count = sum(data(:,1) == uniquepairs(i,1) & data(:,2) == uniquepairs(i,2));
    plot(uniquepairs(i,1),uniquepairs(i,2),'k.','markersize',count*5);
end
plot([0 0.1],[0 0.1]);
xlabel('Replay latency (s)');
ylabel('STA peak latency (s)');
axis square;

%%
% Section 13
% Getting STAs for all of the LGN cells of a certain type
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
whitenoisedata = fetch(conn, 'SELECT fileID, rfX, rfY, neuron, spikeCode FROM WhiteNoiseLGN_forIS WHERE quality = ''1''');
isosampdata = fetch(conn, 'SELECT fileID, rfX, rfY, neuron, cellClass FROM IsoSamp_LGN WHERE quality = ''1''');
close(conn);

% Only looking at the white noise files for cells that are quality 1 in IsoSamp
% Getting filenames. One row per neuron. Potentially multiple files per
% neuron.
uniqueWNneurons = unique([whitenoisedata{:,4}]);
WNfilenames = {};
for neuronid = uniqueWNneurons
    LWN =  [whitenoisedata{:,4}] == neuronid;
    Lisosamp = [isosampdata{:,4}] == neuronid;
    if any(Lisosamp)
        tmp = [];
        for i = find(LWN)       
        	tmp = [tmp, whitenoisedata(i,1)];
        end
        WNfilenames{length(WNfilenames)+1} = tmp;
        
        % Sanity checks
        tmp=[];
        for i = find(Lisosamp)
            tmp=[tmp;isosampdata{i,2:4}];
        end
        for i = find(LWN)
            tmp=[tmp;whitenoisedata{i,2:4}];
        end
        if any(std(tmp))
            error('White noise and IsoSamp entries into database disagree');
        end
    end
end

% Now doing the actual data
data=[];
filenamestems = {};
Ms = [];
for i = 1:length(WNfilenames)
    stro=[];
    for j = 1:length(WNfilenames{i})
        stro = strocat(stro,nex2stro(findfile(WNfilenames{i}{j})));
    end
    spikename = whitenoisedata{strcmp(whitenoisedata(:,1),WNfilenames{i}(1)),5};

    nstixperside = stro.sum.exptParams.nstixperside;
    ntrials = length(stro.sum.absTrialNum);
    spikeidx = find(strcmp(stro.sum.rasterCells(1,:),spikename));
    
    % Reconstructing the M matrix and gamma table
    fundamentals = stro.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = stro.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    Ms(:,:,i) = fundamentals'*mon_spd;

    maxT = 6;
    out = getWhtnsStats(stro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
    STA = out{1}/out{3};
    STA = reshape(STA,[nstixperside.^2  3 maxT]);
 
    noise = sum(STA(:,:,1).^2,2)*maxT;
    energy = sum(sum(STA.^2,2),3);
    whichpix_tmp = energy > max(noise) | energy == max(energy);
    
    if sum(whichpix_tmp) > 1
        [tmpi,tmpj] = ind2sub([nstixperside nstixperside],find(whichpix_tmp));
        ij = [tmpi,tmpj];
        T = clusterdata(ij,'linkage','single','distance','euclidean','criterion','distance','cutoff',sqrt(2));
        
        clusternmembers = [];
        for k =1:max(T)
            clusternmembers(k) = sum(T == k);
        end
        dominantcluster = find(clusternmembers == max(clusternmembers));
        
        clustermat = zeros(nstixperside, nstixperside);
        clustermat(sub2ind(size(clustermat),ij(T==dominantcluster,1),ij(T==dominantcluster,2))) = 1;
        
        if sum(clustermat(:)) > 2 & length(unique(ij(T==dominantcluster,1))) > 1 & length(unique(ij(T==dominantcluster,2))) > 1
            % Then get convex hull
            dominantclusteridxs = ij(T==dominantcluster,:);
            K=convhull(dominantclusteridxs(:,1), dominantclusteridxs(:,2));
            tmp = dominantclusteridxs(K,:);
            [x,y] = meshgrid(1:nstixperside,1:nstixperside);
            inRF = reshape(inpolygon(x(:),y(:),dominantclusteridxs(K,2),dominantclusteridxs(K,1)),[nstixperside, nstixperside]);
        else
            inRF = clustermat;
        end
    else
        inRF = reshape(whichpix_tmp,nstixperside,nstixperside);
    end
    reshape(inRF,nstixperside,nstixperside)
    
%    whichpix_tmp = conv2(reshape(whichpix,nstixperside,nstixperside),ones(3),'same');
    temporalSTA = squeeze(STA(logical(inRF(:)),:,:));
    if ndims(temporalSTA) > 2
        temporalSTA = mean(temporalSTA,1);
    end
    %temporalSTA = squeeze(STA(whichpix,:,:));
    
    data(:,:,i) = squeeze(temporalSTA);
    filenames{i} =  WNfilenames{i}{1};
end

% --------------
% Setting up clusters by hand and seeing where they fall in the t-SNE plot
% 1) 1 = P, 2 = M, 3 = K
% 2) 0 = OFF, 1 = ON
% 3) 0 = L, 1 = M
% --------------
celltypemat = nan(size(data,3),3);
for j = 1:size(data,3)
    sta = squeeze(data(:,:,j));
    [u,s,v] = svd(sta);
    L = strcmp(filenames{j},whitenoisedata(:,1));
    neuronidx = whitenoisedata{L,4};
    celltype = isosampdata{[isosampdata{:,4}] == neuronidx,5};
    switch(celltype)
        case {'P'}
            celltypemat(j,1) = 1;
        case {'M'}
            celltypemat(j,1) = 2;
        case {'K'}
            celltypemat(j,1) = 3;
        otherwise
            celltypemat(j,1) = 4;
    end
    rgb = u(:,1);
    temporal = v(:,1);
    if max(abs(rgb)) == -min(rgb)
        rgb = -rgb;
        temporal = -temporal;
    end
    if max(abs(temporal)) == max(temporal)
        celltypemat(j,2) = 1; % ON
    else
        celltypemat(j,2) = 0; % OFF
    end
    M = Ms(:,:,j);
    lm_ratio = sta(1,2)/sta(2,2); % Just using second frame for cone weights
    if abs(lm_ratio-M(2,1)/M(2,2)) < abs(lm_ratio-M(1,1)/M(1,2))
        celltypemat(j,3) = 1; % M
    else
        celltypemat(j,3) = 0; % L
    end
end

clusterids = {'Parvo L-cone ON','Parvo M-cone ON','Parvo L-cone OFF','Parvo M-cone OFF','Magno ON','Magno OFF','K'};
idx = zeros(size(celltypemat,1),1);
nclusters = length(clusterids);
for i = 1:nclusters
    done = 0;
    L = zeros(size(celltypemat,1),1);
    if strfind(clusterids{i},'Parvo')
        L = L | celltypemat(:,1) == 1;
    elseif strfind(clusterids{i},'Magno')
        L = L | celltypemat(:,1) == 2;
    elseif strfind(clusterids{i},'K')
        L = L | celltypemat(:,1) == 3;
        done = 1;
    end
    if ~done
        if strfind(clusterids{i},'ON')
            L = L & celltypemat(:,2) == 1;
        else
            L = L & celltypemat(:,2) == 0;
        end
        if strfind(clusterids{i},'Parvo')
            if strfind(clusterids{i},'L-cone')
                L = L & celltypemat(:,3) == 0;
            elseif strfind(clusterids{i},'M-cone')
                L = L & celltypemat(:,3) == 1;
            end
        end
    end
    idx(logical(L)) = i;
end

tab = reshape(permute(data,[2 1 3]),size(data,1)*size(data,2),size(data,3)); % guns, time, cell
y = tsne(tab','Algorithm','exact','Distance','euclidean','LearnRate',1,'Perplexity',5,'Exaggeration',1);

figure; axes; hold on;
colors = [0 0 0; 1 0 0; 0 .5 0; 0 0 1; 1 0 1; 0 1 1; 1 .5 0];
for i = 1:nclusters
    for j = find(idx == i)'
        h = plot(y(j,1),y(j,2),'ko','MarkerFaceColor',colors(i,:));
        set(h,'ButtonDownFcn',['disp([''',clusterids{i},': ',filenames{j},'''])'])
    end
end

olddata = data;
% Reducing the dimensionality of each STA via SVD
newdata = zeros(size(data,1)+size(data,2),size(data,3));
for i = 1:size(data,3)
    sta = squeeze(data(:,:,i));
    [u,s,v] = svd(sta);
    if sta(:,2)'*u(:,1) < 0
        u = -u;
        v = -v;
    end
    newdata(:,i) = [u(:,1);v(:,1)];
end
y = tsne(newdata','Algorithm','exact','Distance','euclidean','LearnRate',1,'Perplexity',5,'Exaggeration',1);

figure; axes; hold on;
for i = 1:nclusters
    for j = find(idx == i)'
        h = plot(y(j,1),y(j,2),'ko','MarkerFaceColor',colors(i,:));
        set(h,'ButtonDownFcn',['disp([''',clusterids{i},': ',filenames{j},'''])'])
    end
end

% STAs for each cluster
figure;
framestep = stro.sum.exptParams.nrepframes/stro.sum.exptParams.framerate*1000; % ms/frame
t = [0:size(data,2)-1]*framestep;
for i = 1:nclusters
    subplot(ceil(sqrt(nclusters)),ceil(sqrt(nclusters)),i); hold on;
    set(gca,'ColorOrder',[1 0 0; 0 1 0; 0 0 1],'Ycolor',colors(i,:),'Xcolor',colors(i,:));
    for j = find(idx == i)'
        plot(t,squeeze(data(:,:,j))');
    end
    set(gca,'Xlim',[t(1) t(end)]);
end
xlabel('Time (ms)');
ylabel('STA (A.U.)');

